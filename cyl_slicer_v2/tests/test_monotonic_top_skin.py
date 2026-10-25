from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.gcode import _ordered_layer_paths, emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.models import GcodeState, LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


def slab(index=0, *, top=None, bottom=None):
    occupied = box(0, 0, 1, 10)
    blank = GeometryCollection()
    return LayerSlice(index, 30.1 + .2 * index, occupied,
                      top if top is not None else blank,
                      bottom if bottom is not None else blank, blank, blank)


def config(**kwargs):
    return SlicerConfig(wall_count=1, layer_height=.2, nozzle_width=.4,
                        min_path_length=0, bed_radius=29, **kwargs)


def total_length(paths):
    return sum(math.hypot((a[0]-b[0])*path.radius, a[1]-b[1])
        for path in paths for a, b in zip(path.points, path.points[1:]))


class MonotonicTopSkinTests(unittest.TestCase):
    def test_disabled_keeps_solid_infill_role(self):
        layer = build_layer_toolpaths(slab(top=box(0,0,1,10)),
                                     config(monotonic_top_skin=False))
        self.assertTrue(layer.paths_by_role.get("SOLID_INFILL"))
        self.assertFalse(layer.paths_by_role.get("TOP_SKIN"))

    def test_top_only_has_oneway_stripes_and_same_deposited_geometry(self):
        original = build_layer_toolpaths(slab(top=box(0,0,1,10)),
                                         config(monotonic_top_skin=False))
        mono = build_layer_toolpaths(slab(top=box(0,0,1,10)),
                                     config(monotonic_top_skin=True))
        old_lines = original.paths_by_role["SOLID_INFILL"]
        top_lines = mono.paths_by_role["TOP_SKIN"]
        self.assertFalse(mono.paths_by_role.get("SOLID_INFILL"))
        self.assertAlmostEqual(total_length(old_lines), total_length(top_lines), places=7)
        seq = _ordered_layer_paths(mono, GcodeState(), config(monotonic_top_skin=True))
        skin = [p for p in seq if p.role == "TOP_SKIN"]
        self.assertTrue(skin)
        stripes = [round(p.points[0][0], 9) for p in skin]
        self.assertEqual(stripes, sorted(stripes))
        for path in skin:
            self.assertGreater(path.points[-1][1], path.points[0][1],
                               "all exposed stripes must print in the same axial direction")

    def test_alternate_top_skin_direction_is_consistent_on_rotary_passes(self):
        c = config(monotonic_top_skin=True, solid_infill_pattern="alternate")
        for index, expect_rotary in ((0, False), (1, True)):
            lp = build_layer_toolpaths(slab(index, top=box(0,0,1,10)), c)
            seq = [p for p in _ordered_layer_paths(lp, GcodeState(), c) if p.role == "TOP_SKIN"]
            self.assertTrue(seq)
            coords = [p.points[0][1] if expect_rotary else p.points[0][0] for p in seq]
            self.assertEqual(coords, sorted(coords))
            for p in seq:
                if expect_rotary:
                    self.assertGreater(p.points[-1][0], p.points[0][0])
                else:
                    self.assertGreater(p.points[-1][1], p.points[0][1])

    def test_top_bottom_overlap_does_not_double_extrude(self):
        region = box(0,0,1,10)
        only = build_layer_toolpaths(slab(top=region), config(monotonic_top_skin=True))
        both = build_layer_toolpaths(slab(top=region, bottom=region),
                                     config(monotonic_top_skin=True))
        self.assertFalse(both.paths_by_role.get("SOLID_INFILL"))
        self.assertAlmostEqual(
            total_length(only.paths_by_role.get("TOP_SKIN", [])),
            total_length(both.paths_by_role.get("TOP_SKIN", [])), places=7)

    def test_bottom_solid_stays_independent_bidirectional_role(self):
        lp = build_layer_toolpaths(slab(bottom=box(0,0,1,10)),
                                   config(monotonic_top_skin=True))
        self.assertFalse(lp.paths_by_role.get("TOP_SKIN"))
        self.assertTrue(lp.paths_by_role.get("SOLID_INFILL"))

    def test_gcode_preserves_total_extruded_distance(self):
        region = box(0,0,1,10)
        traditional = build_layer_toolpaths(slab(top=region), config())
        monotonic = build_layer_toolpaths(slab(top=region), config(monotonic_top_skin=True))
        a = emit_gcode([traditional], config(retract_mm=0))
        b = emit_gcode([monotonic], config(monotonic_top_skin=True, retract_mm=0))
        self.assertIn(";MONOTONIC_TOP_SKIN:true", b)
        self.assertIn(";TYPE:TOP_SKIN", b)
        self.assertNotIn(";TYPE:TOP_SKIN", a)
        def physical_length(gcode):
            return sum(math.hypot(x1-x0, math.radians(y1-y0)*30.1)
                for _, role, x0, y0, x1, y1 in parse_extrusion_segments(gcode)
                if role in {"SOLID_INFILL", "TOP_SKIN"})
        self.assertAlmostEqual(physical_length(a), physical_length(b), delta=.01)

    def test_rejects_nonboolean_flag(self):
        with self.assertRaisesRegex(ValueError, "monotonic_top_skin"):
            config(monotonic_top_skin="yes").validate()


if __name__ == "__main__":
    unittest.main()
