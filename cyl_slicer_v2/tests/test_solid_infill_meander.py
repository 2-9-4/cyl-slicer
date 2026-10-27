from pathlib import Path
import math
import re
import unittest
from unittest.mock import patch

from shapely.geometry import box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import (
    _orient_path_for_current_position,
    emit_gcode,
    parse_extrusion_segments,
)
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath
from cyl_slicer_v2.pipeline import slice_to_gcode


class SolidInfillMeanderTests(unittest.TestCase):
    def test_parallel_solid_passes_meander_in_both_directions(self):
        for radial in (False, True):
            for limit in (1024, 0):
                with self.subTest(radial=radial, ordering_limit=limit):
                    paths = [Toolpath("SOLID_INFILL", 30.0,
                        [(0, row), (1, row)] if radial else
                        [(row / 30, 0), (row / 30, 30)])
                        for row in (0.0, 0.4, 0.8, 1.2)]
                    region = box(0, 0, 1, 1.2) if radial else box(0, 0, 1.2 / 30, 30)
                    layer = LayerToolpaths(0, 30.0, {"SOLID_INFILL": paths}, travel_region=region)
                    config = SlicerConfig(axial_min=0, extrusion_per_mm=0.1)
                    with patch("cyl_slicer_v2.gcode.PATH_ORDER_OPTIMIZE_LIMIT", limit):
                        code = emit_gcode([layer], config)
                    segments = parse_extrusion_segments(code)
                    self.assertEqual(len(segments), 4)
                    for index, (_, _, x0, y0, x1, y1) in enumerate(segments):
                        delta = math.radians(y1 - y0) * 30 if radial else x1 - x0
                        self.assertAlmostEqual(delta, 30 * (-1) ** index, delta=.001)
                        if index:
                            previous = segments[index - 1]
                            travel = math.hypot(x0 - previous[4], math.radians(y0 - previous[5]) * 30)
                            self.assertAlmostEqual(travel, .4, delta=.001)
                    self.assertNotIn("; z hop", code)
                    self.assertEqual(code.count("; retract"), 1)  # Layer end only.
                    self.assertNotIn("; unretract", code)
                    moves = [line for line in code.splitlines() if line.startswith("G1 X")]
                    self.assertAlmostEqual(float(re.search(r"\bE([0-9.]+)", moves[-1])[1]), 12)

    def test_reversal_preserves_segment_widths_and_deposited_filament(self):
        points = [(math.radians(deg), 0) for deg in (350, 355, 370)]
        path = Toolpath("SOLID_INFILL", 30, points, segment_widths=[.25, .4])
        oriented = _orient_path_for_current_position(path, GcodeState(theta_accum=points[-1][0]))
        self.assertEqual(oriented.points, list(reversed(points)))
        self.assertEqual(oriented.segment_widths, [.4, .25])
        self.assertEqual(path.points, points)
        self.assertEqual(path.segment_widths, [.25, .4])
        config = SlicerConfig(retract_mm=0)
        codes = [emit_gcode([LayerToolpaths(0, 30, {"SOLID_INFILL": [p]})], config)
                 for p in (path, oriented)]
        deposited = [float(re.findall(r"\bE([0-9.]+)", code)[-1]) for code in codes]
        self.assertAlmostEqual(*deposited, places=5)
        reversed_segments = parse_extrusion_segments(codes[1])
        self.assertAlmostEqual(sum(s[5] - s[3] for s in reversed_segments), -20, places=3)

    def test_rotary_meander_preserves_seam_and_long_winding(self):
        for start, end in ((math.radians(350), math.radians(370)), (0, 4.5), (0, 2 * math.pi)):
            with self.subTest(start=start, end=end):
                paths = [Toolpath("SOLID_INFILL", 30, [(start, x), (end, x)]) for x in (0, .4)]
                code = emit_gcode([LayerToolpaths(0, 30, {"SOLID_INFILL": paths})],
                                  SlicerConfig(retract_mm=0, extrusion_per_mm=.1))
                segments = parse_extrusion_segments(code)
                self.assertEqual(len(segments), 2)
                deltas = [s[5] - s[3] for s in segments]
                # A complete turn has coincident endpoints, so either sign is
                # equally near. Shorter passes must reverse on the second row.
                if end - start < 2 * math.pi:
                    self.assertAlmostEqual(deltas[0], -deltas[1], places=3)
                for delta in deltas:
                    self.assertAlmostEqual(abs(delta), math.degrees(end - start), places=3)
                deposited = float(re.findall(r"\bE([0-9.]+)", code)[-1])
                self.assertAlmostEqual(deposited, 2 * (end - start) * 30 * .1, places=5)

    def test_rotary_support_passes_meander_without_changing_extrusion(self):
        for role in ("SUPPORT", "SUPPORT_INTERFACE"):
            for start, end in ((math.radians(350), math.radians(370)), (0, 4.5), (0, 2 * math.pi)):
                with self.subTest(role=role, start=start, end=end):
                    paths = [Toolpath(role, 30, [(start, x), (end, x)]) for x in (0, .4)]
                    region = box(start, 0, end, .4)
                    code = emit_gcode([LayerToolpaths(0, 30, {role: paths}, travel_region=region)],
                                      SlicerConfig(extrusion_per_mm=.1))
                    segments = parse_extrusion_segments(code)
                    self.assertEqual(len(segments), 2)
                    deltas = [s[5] - s[3] for s in segments]
                    if end - start < 2 * math.pi:
                        self.assertAlmostEqual(deltas[0], -deltas[1], places=3)
                    for delta in deltas:
                        self.assertAlmostEqual(abs(delta), math.degrees(end - start), places=3)
                    self.assertEqual(code.count("; retract"), 1)  # Layer end only.
                    e = float(re.findall(r"\bE([0-9.]+)", code)[-1])
                    self.assertAlmostEqual(e + 1.4, 2 * (end - start) * 30 * .1, places=5)

    def test_rotary_support_reversal_keeps_segment_widths(self):
        for role in ("SUPPORT", "SUPPORT_INTERFACE"):
            path = Toolpath(role, 30, [(0, 0), (.1, 0), (.3, 0)], segment_widths=[.25, .4])
            oriented = _orient_path_for_current_position(path, GcodeState(theta_accum=.3))
            self.assertEqual(oriented.points, list(reversed(path.points)))
            self.assertEqual(oriented.segment_widths, [.4, .25])
            config = SlicerConfig(retract_mm=0)
            codes = [emit_gcode([LayerToolpaths(0, 30, {role: [p]})], config) for p in (path, oriented)]
            amounts = [float(re.findall(r"\bE([0-9.]+)", code)[-1]) for code in codes]
            self.assertAlmostEqual(*amounts, places=5)

    def test_default_retraction_coupon_layer_one_meanders(self):
        root = Path(__file__).resolve().parents[2]
        code = slice_to_gcode(load_mesh(root / "calibration/retraction_posts.stl"),
                              SlicerConfig(), max_layers=7)
        segments = [s for s in parse_extrusion_segments(code)
                    if s[0] == 1 and s[1] == "SOLID_INFILL"]
        self.assertGreater(len(segments), 2)
        for previous, current in zip(segments, segments[1:]):
            self.assertLess((previous[5] - previous[3]) * (current[5] - current[3]), 0)
            self.assertAlmostEqual(previous[5], current[3], places=3)
        shell = code.split(";LAYER:1\n")[1].split(";LAYER:2\n")[0]
        passes = [part.split(";TYPE:")[0] for part in shell.split(";TYPE:SOLID_INFILL")[1:]]
        self.assertTrue(all("; retract" not in part for part in passes[:-1]))
        self.assertTrue(all("; unretract" not in part for part in passes[1:]))

    def test_default_coupon_top_surfaces_meander_in_both_axes(self):
        root = Path(__file__).resolve().parents[2]
        mesh = load_mesh(root / "calibration/retraction_posts.stl")
        config = SlicerConfig()
        codes = [slice_to_gcode(mesh, config, max_layers=17) for config in
                 (config, SlicerConfig(monotonic_top_skin=True))]
        for index in (10, 11):
            with self.subTest(layer=index):
                fills = [[s for s in parse_extrusion_segments(code)
                          if s[0] == index and s[1] in {"SOLID_INFILL", "TOP_SKIN"}]
                         for code in codes]
                self.assertGreater(len(fills[0]), 4)
                # The first four selected stripes lie in one uninterrupted
                # band between the post footprints.
                for previous, current in zip(fills[0][:3], fills[0][1:4]):
                    if index == 10:
                        self.assertLess((previous[4] - previous[2]) * (current[4] - current[2]), 0)
                        self.assertAlmostEqual(previous[4], current[2], places=3)
                    else:
                        self.assertLess((previous[5] - previous[3]) * (current[5] - current[3]), 0)
                        self.assertAlmostEqual(previous[5], current[3], places=3)
                radius = config.bed_radius + (index + .5) * config.layer_height
                lengths = [sum(math.hypot(s[4] - s[2], math.radians(s[5] - s[3]) * radius) for s in fill)
                           for fill in fills]
                self.assertAlmostEqual(*lengths, delta=.01)
                self.assertTrue(all(s[1] == "SOLID_INFILL" for s in fills[0]))
                self.assertTrue(all(s[1] == "TOP_SKIN" for s in fills[1]))


if __name__ == "__main__":
    unittest.main()
