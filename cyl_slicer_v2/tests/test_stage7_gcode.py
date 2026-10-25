from __future__ import annotations

from pathlib import Path
import math
import re
import unittest
from unittest.mock import patch

from shapely.geometry import GeometryCollection, LineString, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.calibration import RetractionSweep
from cyl_slicer_v2.coordinates import cylindrical_distance
from cyl_slicer_v2.gcode import _orient_path_for_current_position, _ordered_layer_paths, emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.intersection_backend import sample_intersection_layer
from cyl_slicer_v2.layers import layer_radii, slice_mesh
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import GcodeState, LayerSlice, LayerToolpaths, Toolpath
from cyl_slicer_v2.pipeline import slice_to_gcode
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, build_toolpaths


class GcodeTests(unittest.TestCase):
    def test_closed_perimeter_entry_projects_onto_nearest_edge_and_preserves_length(self):
        path = Toolpath("INTERNAL_PERIMETER", 10.0,
                        [(0.04, 0.4), (0.2, 0.4), (0.2, 9.6), (0.04, 9.6), (0.04, 0.4)])
        oriented = _orient_path_for_current_position(path, GcodeState(theta_accum=0.25, axial=5.0))
        self.assertEqual(oriented.points[0], oriented.points[-1])
        self.assertAlmostEqual(oriented.points[0][0], 0.2)
        self.assertAlmostEqual(oriented.points[0][1], 5.0)
        length = lambda p: sum(cylindrical_distance(a, b, p.radius) for a, b in zip(p.points, p.points[1:]))
        self.assertAlmostEqual(length(oriented), length(path), places=10)
        signed_area = lambda p: sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(p.points, p.points[1:]))
        self.assertAlmostEqual(signed_area(oriented), signed_area(path), places=10)
        self.assertEqual(path.points[0], (0.04, 0.4), "Do not mutate input geometry")

    def test_adjacent_closed_perimeters_do_not_retract_for_distant_original_starts(self):
        layer = LayerToolpaths(0, 10.0, {
            "EXTERNAL_PERIMETER": [Toolpath("EXTERNAL_PERIMETER", 10.0,
                [(0, 0), (0.24, 0), (0.24, 10), (0, 10), (0, 0)])],
            "INTERNAL_PERIMETER": [Toolpath("INTERNAL_PERIMETER", 10.0,
                [(0.2, 9.6), (0.04, 9.6), (0.04, 0.4), (0.2, 0.4), (0.2, 9.6)])],
        }, travel_region=box(0, 0, 0.24, 10))
        code = emit_gcode([layer], SlicerConfig())
        self.assertEqual(code.count("; retract"), 1)
        self.assertNotIn("; unretract", code)

    def test_next_island_is_selected_from_current_position(self):
        a = Toolpath("EXTERNAL_PERIMETER", 30, [(0, 0), (0.2, 0)])
        b = Toolpath("EXTERNAL_PERIMETER", 30, [(1, 20), (1.2, 20)])
        layer = LayerToolpaths(0, 30, {"EXTERNAL_PERIMETER": [a, b]}, island_paths=[
            {"EXTERNAL_PERIMETER": [a]}, {"EXTERNAL_PERIMETER": [b]},
        ])
        ordered = _ordered_layer_paths(layer, GcodeState(theta_accum=1, axial=20))
        self.assertIs(ordered[0], b)
        self.assertIs(ordered[1], a)

    def test_coupon_finishes_each_island_and_retracts_at_calibration_boundary(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(layer_backend="intersection", layer_height=0.2, bed_radius=30,
                              wall_count=3, feed_mm_per_min=4500, travel_mm_per_min=9000, hotend_temp=195,
                              retraction_sweep=RetractionSweep("retract-mm", (0, 1, 2, 3, 4)))
        code = slice_to_gcode(load_mesh(root / "calibration/retraction_posts.stl"), config, max_layers=30)
        shell23 = code.split(";LAYER:23\n")[1].split(";LAYER:24\n")[0]
        shell24 = code.split(";LAYER:24\n")[1].split(";LAYER:25\n")[0]
        self.assertIn("value=0 unit=mm height=4.8000", shell23)
        self.assertNotIn("; retract", shell23)
        self.assertIn("value=1 unit=mm height=5.0000", shell24)
        self.assertEqual(shell24.count("; retract"), 4, "Three island transfers and the layer end")
        islands = shell24.split(";TYPE:EXTERNAL_PERIMETER")[1:]
        self.assertEqual(len(islands), 4)
        for island in islands:
            self.assertEqual(island.count(";TYPE:INTERNAL_PERIMETER"), 2)
            self.assertIn(";TYPE:INFILL_", island)
            self.assertEqual(island.count("; retract"), 1)
        for line in shell24.splitlines():
            if line.startswith(("G0 ", "G1 X")):
                self.assertRegex(line, r"\bF\d+(?:\.\d+)?",
                                 "Travel/print commands carry a positive Marlin feedrate")
                feed = float(line.rsplit("F", 1)[1].split()[0])
                self.assertGreater(feed, 0)

    def test_short_gap_between_islands_still_retracts(self):
        layer = LayerToolpaths(0, 30, {"INFILL_RADIAL": [
            Toolpath("INFILL_RADIAL", 30, [(0, 0), (1, 0)]),
            Toolpath("INFILL_RADIAL", 30, [(0, 0.4), (1, 0.4)]),
        ]}, travel_region=box(0, -0.05, 1, 0.05).union(box(0, 0.35, 1, 0.45)))
        code = emit_gcode([layer], SlicerConfig())
        self.assertEqual(code.count("; retract"), 2)
        self.assertEqual(code.count("; unretract"), 1)

    def test_part_cooling_starts_at_nominal_radial_height_once(self):
        layers = [LayerToolpaths(i-2, radius) for i, radius in enumerate(
            (30.0, 31.7, 31.75, 32.0))]
        code = emit_gcode(layers, SlicerConfig(bed_radius=30.0, z_hop_height=5.0))
        self.assertLess(code.index("M106 S0 ; part cooling off for adhesion"), code.index("M109"))
        self.assertEqual(code.count("M106 S255 ; part cooling on"), 1)
        self.assertLess(code.index(";LAYER:0"), code.index("M106 S255"))
        self.assertLess(code.index("M106 S255"), code.index(";LAYER:1"))
        self.assertLess(code.index("M400 ; wait for end lift"), code.index("M106 S0 ; fan off"))

    def test_part_cooling_can_be_disabled_or_set_to_partial_from_first_shell(self):
        layers = [LayerToolpaths(0, 30.2)]
        disabled = emit_gcode(layers, SlicerConfig(bed_radius=30, fan_speed=0, fan_start_height=0))
        self.assertNotIn("; part cooling on", disabled)
        partial = emit_gcode(layers, SlicerConfig(bed_radius=30, fan_speed=128, fan_start_height=0))
        self.assertIn("M106 S128 ; part cooling on", partial)

    def test_part_cooling_rejects_invalid_settings(self):
        for kwargs in ({"fan_speed": -1}, {"fan_speed": 256}, {"fan_speed": 0.5},
                       {"fan_start_height": -1}, {"fan_start_height": float("nan")}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                SlicerConfig(**kwargs).validate()

    def test_end_lift_adds_ten_mm_to_final_radial_height_before_shutdown(self):
        for final_radius in (32.0, 46.9, 75.0):
            with self.subTest(final_radius=final_radius):
                config = SlicerConfig(bed_radius=30.0)
                layers = [LayerToolpaths(0, 30.2), LayerToolpaths(1, final_radius)]
                code = emit_gcode(layers, config)
                expected = f"G0 Z{final_radius - 30.0 + config.layer_height / 2 + 10.0:.4f} F600.0 ; end lift"
                self.assertIn(expected, code)
                self.assertLess(code.index("G90 ; absolute positioning for end lift"),
                                code.index(expected))
                wait = code.index("M400 ; wait for end lift")
                self.assertLess(code.index(expected), wait)
                self.assertLess(wait, code.index("M104 S0"))
                self.assertLess(wait, code.index("M84"))

    def test_empty_slice_does_not_move_to_an_assumed_park_height(self):
        code = emit_gcode([], SlicerConfig(bed_radius=30.0))
        self.assertNotIn("G0 Z", code)
        self.assertIn("M104 S0", code)
        self.assertIn("M84", code)

    def test_gcode_contains_layers_and_types(self):
        gcode = slice_to_gcode(make_fixture("tube"), SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0), max_layers=2)
        self.assertIn(";LAYER:0", gcode)
        self.assertIn(";TYPE:EXTERNAL_PERIMETER", gcode)
        self.assertIn(";TYPE:SOLID_INFILL", gcode)

    def test_extrusion_is_monotonic_between_layer_resets(self):
        gcode = slice_to_gcode(make_fixture("solid_box"), SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0), max_layers=2)
        values = []
        for line in gcode.splitlines():
            if line.startswith("G92 ") and "E0" in line:
                self.assertEqual(values, sorted(values))
                values = []
                continue
            if line.startswith("G1 ") and (" X" in line or " Y" in line):
                match = re.search(r"\bE(-?\d+(?:\.\d+)?)", line)
                if match:
                    values.append(float(match.group(1)))
        self.assertEqual(values, sorted(values))

    def test_retractions_are_emitted_between_paths(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 10.0)]),
                    Toolpath("SUPPORT", 10.0, [(0.1, 0.0), (0.1, 10.0)]),
                ]
            },
        )

        gcode = emit_gcode([layer], SlicerConfig(max_axial_segment=20.0, retract_mm=1.0))

        self.assertIn("; retract", gcode)
        self.assertIn("; unretract", gcode)
        self.assertRegex(gcode, r"(?m)^G1 E-?[0-9.]+ F[0-9.]+ ; retract$")

    def test_long_travel_moves_use_z_hop(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 1.0)]),
                    Toolpath("SUPPORT", 10.0, [(0.0, 10.0), (0.0, 11.0)]),
                ]
            },
        )

        gcode = emit_gcode(
            [layer],
            SlicerConfig(z_hop_height=0.5, z_hop_min_travel=2.0, max_axial_segment=20.0, retract_mm=0.0),
        )

        self.assertIn("G0 Z10.7500 F3000.0 ; z hop", gcode)
        self.assertIn("G0 Z10.2500 F3000.0 ; z unhop", gcode)
        self.assertLess(gcode.index("; z hop"), gcode.index("G0 X25.000 Y0.000"))

    def test_short_travel_moves_do_not_z_hop(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 1.0)]),
                    Toolpath("SUPPORT", 10.0, [(0.01, 1.0), (0.01, 2.0)]),
                ]
            },
        )

        gcode = emit_gcode(
            [layer],
            SlicerConfig(z_hop_height=0.5, z_hop_min_travel=5.0, max_axial_segment=20.0, retract_mm=0.0),
        )

        self.assertNotIn("; z hop", gcode)
        self.assertNotIn("; z unhop", gcode)

    def test_perimeter_crossing_travel_moves_use_z_hop_below_long_travel_threshold(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "EXTERNAL_PERIMETER": [
                    Toolpath(
                        "EXTERNAL_PERIMETER",
                        10.0,
                        [(0.4, -1.0), (0.6, -1.0), (0.6, 1.0), (0.4, 1.0), (0.4, -1.0)],
                    )
                ],
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.7, 0.0), (0.7, 0.5)]),
                ],
            },
        )

        gcode = emit_gcode(
            [layer],
            SlicerConfig(
                z_hop_height=0.5,
                z_hop_min_travel=100.0,
                max_axial_segment=20.0,
                retract_mm=0.0,
                min_path_length=0.0,
            ),
        )

        self.assertEqual(gcode.count("; z hop"), 1)
        self.assertIn("G0 Z10.7500 F3000.0 ; z hop\nG0 X16.000 Y40.107", gcode)

    def test_short_same_role_travels_do_not_retract(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 1.0)]),
                    Toolpath("SUPPORT", 10.0, [(0.01, 1.0), (0.01, 2.0)]),
                ]
            },
        )

        gcode = emit_gcode([layer], SlicerConfig(nozzle_width=0.4, max_axial_segment=20.0, retract_mm=1.0))

        self.assertEqual(gcode.count("; retract"), 1)
        self.assertEqual(gcode.count("; unretract"), 0)

    def test_adjacent_infill_travel_along_perimeter_does_not_hop(self):
        # The travel passes a tessellation vertex on the inner wall. Exercise
        # both a continuous polyline and walls supplied as separate fragments.
        for fragmented in (False, True):
            with self.subTest(fragmented=fragmented):
                wall = [(0.0, 0.0), (0.02, 0.0), (0.04, 0.0),
                        (0.04, 10.0), (0.02, 10.0), (0.0, 10.0), (0.0, 0.0)]
                walls = [[a, b] for a, b in zip(wall, wall[1:])] if fragmented else [wall]
                layer = LayerToolpaths(0, 10.0, {
                    "INTERNAL_PERIMETER": [Toolpath("INTERNAL_PERIMETER", 10.0, pts) for pts in walls],
                    "SOLID_INFILL": [
                        Toolpath("SOLID_INFILL", 10.0, [(0.0, 0.0), (0.0, 10.0)]),
                        Toolpath("SOLID_INFILL", 10.0, [(0.04, 0.0), (0.04, 10.0)]),
                    ],
                })
                code = emit_gcode([layer], SlicerConfig(z_hop_min_travel=100.0))
                second_infill = code.split(";TYPE:SOLID_INFILL")[2]
                self.assertNotIn("; z hop", second_infill)
                self.assertNotIn("; unretract", second_infill)

    def test_straight_moves_are_single_blocks_with_correct_extrusion(self):
        for points in (
            [(0.0, 0.0), (0.0, 0.02), (0.0, 0.5), (0.0, 0.52), (0.0, 30.0)],
            [(0.0, 0.0), (0.02, 0.0), (0.5, 0.0), (0.52, 0.0), (1.0, 0.0)],
            [(0.0, 0.0), (0.02, 0.6), (0.5, 15.0), (0.52, 15.6), (1.0, 30.0)],
        ):
            with self.subTest(points=points):
                layer = LayerToolpaths(0, 30.0, {
                    "EXTERNAL_PERIMETER": [Toolpath("EXTERNAL_PERIMETER", 30.0, points)],
                })
                code = emit_gcode([layer], SlicerConfig(axial_min=0.0, retract_mm=0.0, extrusion_per_mm=0.1))
                segments = parse_extrusion_segments(code)
                self.assertEqual(len(segments), 1)
                self.assertAlmostEqual(segments[0][4], points[-1][1], places=3)
                self.assertAlmostEqual(segments[0][5], math.degrees(points[-1][0]), places=3)
                expected_e = math.hypot(points[-1][0] * 30.0, points[-1][1]) * 0.1
                self.assertAlmostEqual(float(re.findall(r"\bE([0-9.]+)", code)[-1]), expected_e, places=5)

    def test_coalescing_preserves_corners_reversals_and_full_turns(self):
        paths = [
            Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.1, 0.0), (0.2, 0.0), (0.2, 1.0), (0.1, 1.0), (0.0, 1.0)]),
            Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 1.0), (0.0, 0.0)]),
            Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (math.pi, 0.0), (2 * math.pi, 0.0)]),
        ]
        for path, count, length in zip(paths, (3, 2, 1), (5.0, 2.0, 20 * math.pi)):
            with self.subTest(points=path.points):
                code = emit_gcode([LayerToolpaths(0, 10.0, {"SUPPORT": [path]})],
                                  SlicerConfig(axial_min=0.0, retract_mm=0.0, extrusion_per_mm=0.1))
                self.assertEqual(len(parse_extrusion_segments(code)), count)
                self.assertAlmostEqual(float(re.findall(r"\bE([0-9.]+)", code)[-1]), length * 0.1, places=5)

    def test_retraction_coupon_first_layer_has_no_infill_hops_or_fragmentation(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(layer_backend="intersection", layer_height=0.2, bed_radius=30.0, wall_count=3)
        code = slice_to_gcode(load_mesh(root / "calibration/retraction_posts.stl"), config, max_layers=1)
        infill = code.split(";TYPE:SOLID_INFILL")[1:]
        self.assertEqual(len(infill), 76)
        # A nearest-point perimeter seam may leave a long travel to the first
        # skin pass. Transfers between subsequent adjacent passes must not hop.
        self.assertTrue(all("; z hop" not in path for path in infill[1:]))
        self.assertTrue(all(len(parse_extrusion_segments(";TYPE:SOLID_INFILL\n" + path)) == 1 for path in infill))
        moves = [line for line in code.splitlines() if line.startswith("G1 X")]
        self.assertAlmostEqual(float(re.search(r"\bE([0-9.]+)", moves[-1])[1]), 94.77276, places=5)

    def test_layer_changes_reset_extrusion_and_rotary_axes(self):
        gcode = slice_to_gcode(make_fixture("solid_box"), SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0), max_layers=2)
        layer_count = len(re.findall(r"(?m)^;LAYER:", gcode))
        reset_count = len(re.findall(r"(?m)^G92 E0 Y0 ; reset E/Y accumulators for layer$", gcode))
        self.assertEqual(reset_count, layer_count)

    def test_parse_extrusion_segments_applies_g92_rotary_offsets(self):
        gcode = "\n".join(
            [
                ";LAYER:0",
                ";TYPE:EXTERNAL_PERIMETER",
                "G0 X0.000 Y350.000",
                "G1 X0.000 Y355.000 E1.00000",
                "G1 E0.00000 ; retract",
                "G1 E1.00000 ; unretract",
                "G92 E0 Y0",
                "G0 X0.000 Y5.000",
                "G1 X0.000 Y10.000 E1.00000",
                "",
            ]
        )

        segments = parse_extrusion_segments(gcode)

        self.assertEqual(
            segments,
            [
                (0, "EXTERNAL_PERIMETER", 0.0, 350.0, 0.0, 355.0),
                (0, "EXTERNAL_PERIMETER", 0.0, 360.0, 0.0, 365.0),
            ],
        )

    def test_forward_only_roles_do_not_move_backward_in_y(self):
        gcode = slice_to_gcode(make_fixture("helix_like_ribbon"), SlicerConfig(theta_step_deg=20, axial_step=2.0, layer_height=1.0), max_layers=2)
        segments = parse_extrusion_segments(gcode)
        self.assertTrue(segments)
        self.assertTrue(
            all(
                end_y >= start_y
                for _, role, _, start_y, _, end_y in segments
                if role
                not in {
                    "EXTERNAL_PERIMETER",
                    "INTERNAL_PERIMETER",
                    "SUPPORT",
                    "SUPPORT_INTERFACE",
                    "INFILL_RADIAL",
                }
            )
        )

    def test_constant_theta_paths_are_oriented_from_current_position(self):
        config = SlicerConfig(max_axial_segment=20.0)
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 10.0)]),
                    Toolpath("SUPPORT", 10.0, [(0.1, 0.0), (0.1, 10.0)]),
                ]
            },
        )

        segments = [segment for segment in parse_extrusion_segments(emit_gcode([layer], config)) if segment[1] == "SUPPORT"]

        self.assertEqual(len(segments), 2)
        self.assertLess(segments[0][2], segments[0][4])
        self.assertGreater(segments[1][2], segments[1][4])

    def test_radial_infill_snakes_with_short_travels_and_correct_extrusion(self):
        paths = [Toolpath("INFILL_RADIAL", 30.0, [(0.0, x), (1.0, x)]) for x in (0.0, 0.4, 0.8, 1.2)]
        layer = LayerToolpaths(0, 30.0, {"INFILL_RADIAL": paths})
        # Also exercise the fallback used when there are too many paths to sort.
        for limit in (1024, 0):
            with self.subTest(ordering_limit=limit), patch("cyl_slicer_v2.gcode.PATH_ORDER_OPTIMIZE_LIMIT", limit):
                code = emit_gcode([layer], SlicerConfig(axial_min=0.0, extrusion_per_mm=0.1))
                segments = parse_extrusion_segments(code)
                self.assertEqual(len(segments), 4)
                for index, (_, _, x0, y0, x1, y1) in enumerate(segments):
                    self.assertAlmostEqual(x0, index * 0.4)
                    self.assertEqual(x0, x1)
                    self.assertAlmostEqual(y1 - y0, math.degrees(1.0) * (-1) ** index, places=3)
                    if index:
                        self.assertAlmostEqual(y0, segments[index - 1][5])
                self.assertNotIn("; z hop", code)
                self.assertEqual(code.count("; retract"), 1)
                self.assertNotIn("; unretract", code)
                moves = [line for line in code.splitlines() if line.startswith("G1 X")]
                self.assertAlmostEqual(float(re.search(r"\bE([0-9.]+)", moves[-1])[1]), 12.0)

    def test_reversed_radial_infill_preserves_seam_crossings_and_winding(self):
        for start, end in ((math.radians(350), math.radians(370)), (0.0, 4.5)):
            with self.subTest(start=start, end=end):
                layer = LayerToolpaths(0, 30.0, {"INFILL_RADIAL": [
                    Toolpath("INFILL_RADIAL", 30.0, [(start, 0.0), (end, 0.0)]),
                    Toolpath("INFILL_RADIAL", 30.0, [(start, 0.4), (end, 0.4)]),
                ]})
                code = emit_gcode([layer], SlicerConfig(extrusion_per_mm=0.1, retract_mm=0.0))
                segments = parse_extrusion_segments(code)
                self.assertEqual(len(segments), 2)
                self.assertAlmostEqual(segments[0][5] - segments[0][3], math.degrees(end - start), places=3)
                self.assertAlmostEqual(segments[1][5] - segments[1][3], -math.degrees(end - start), places=3)
                deposited = float(re.findall(r"\bE([0-9.]+)", code)[-1])
                self.assertAlmostEqual(deposited, 2 * (end - start) * 30 * 0.1, places=5)

        reverse_turn = Toolpath("INFILL_RADIAL", 30.0, [(2 * math.pi, 0.0), (0.0, 0.0)])
        code = emit_gcode([LayerToolpaths(0, 30.0, {"INFILL_RADIAL": [reverse_turn]})],
                          SlicerConfig(retract_mm=0.0, extrusion_per_mm=0.1))
        segment, = parse_extrusion_segments(code)
        self.assertAlmostEqual(segment[5] - segment[3], -360.0)
        self.assertAlmostEqual(float(re.findall(r"\bE([0-9.]+)", code)[-1]), 6 * math.pi, places=5)

    def test_retraction_coupon_layer_three_radial_infill_snakes(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(layer_backend="intersection", layer_height=0.2, bed_radius=30.0, wall_count=3)
        # Leave enough later shells that layer 3 is not capped by max_layers.
        code = slice_to_gcode(load_mesh(root / "calibration/retraction_posts.stl"), config, max_layers=7)
        segments = [s for s in parse_extrusion_segments(code) if s[0] == 3 and s[1] == "INFILL_RADIAL"]
        self.assertEqual(len(segments), 17)
        for previous, current in zip(segments, segments[1:]):
            self.assertLess((previous[5] - previous[3]) * (current[5] - current[3]), 0.0)
            self.assertAlmostEqual(current[3], previous[5])
            self.assertAlmostEqual(abs(current[2] - previous[4]), 2.0)
        shell = code.split(";LAYER:3\n")[1].split(";LAYER:4\n")[0]
        radial = [part.split(";TYPE:")[0] for part in shell.split(";TYPE:INFILL_RADIAL")[1:]]
        self.assertTrue(all("; retract" not in part for part in radial[:-1]))
        self.assertTrue(all("; unretract" not in part for part in radial[1:]))

    def test_neighboring_infill_passes_in_material_do_not_retract(self):
        for role in ("INFILL_RADIAL", "INFILL_AXIAL", "SUPPORT"):
            for limit in (1024, 0):
                with self.subTest(role=role, ordering_limit=limit), patch("cyl_slicer_v2.gcode.PATH_ORDER_OPTIMIZE_LIMIT", limit):
                    spacing = 1.5 if role == "SUPPORT" else 2.0
                    if role == "INFILL_RADIAL":
                        paths = [Toolpath(role, 30.0, [(0, x), (1, x)]) for x in (0, spacing, 2 * spacing)]
                        region = box(0, 0, 1, 2 * spacing)
                    else:
                        paths = [Toolpath(role, 30.0, [(t, 0), (t, 6)]) for t in (0, spacing / 30, 2 * spacing / 30)]
                        region = box(0, 0, 2 * spacing / 30, 6)
                    layer = LayerToolpaths(0, 30.0, {role: paths}, travel_region=region)
                    code = emit_gcode([layer], SlicerConfig())
                    self.assertEqual(code.count("; retract"), 1, "Only retract at layer end")
                    self.assertNotIn("; unretract", code)

    def test_short_infill_travels_across_gaps_or_holes_still_retract(self):
        paths = [Toolpath("INFILL_RADIAL", 30.0, [(0, x), (1, x)]) for x in (0, 2)]
        regions = [
            box(-0.1, -0.1, 1.1, 0.1).union(box(-0.1, 1.9, 1.1, 2.1)),
            box(-0.1, -0.1, 1.1, 2.1).difference(box(0.9, 0.4, 1.05, 1.6)),
        ]
        for region in regions:
            with self.subTest(region=region.geom_type):
                layer = LayerToolpaths(0, 30.0, {"INFILL_RADIAL": paths}, travel_region=region)
                code = emit_gcode([layer], SlicerConfig())
                self.assertEqual(code.count("; retract"), 2)
                self.assertEqual(code.count("; unretract"), 1)

    def test_long_infill_travel_inside_material_still_retracts(self):
        paths = [Toolpath("INFILL_RADIAL", 30.0, [(0, x), (1, x)]) for x in (0, 6)]
        layer = LayerToolpaths(0, 30.0, {"INFILL_RADIAL": paths}, travel_region=box(0, 0, 1, 6))
        code = emit_gcode([layer], SlicerConfig())
        self.assertEqual(code.count("; retract"), 2)
        self.assertEqual(code.count("; unretract"), 1)

    def test_infill_role_change_does_not_force_retraction(self):
        layer = LayerToolpaths(0, 30.0, {
            "INFILL_RADIAL": [Toolpath("INFILL_RADIAL", 30.0, [(0, 0), (1, 0)])],
            "INFILL_AXIAL": [Toolpath("INFILL_AXIAL", 30.0, [(1, 0), (1, 5)])],
        }, travel_region=box(0, 0, 1, 5))
        code = emit_gcode([layer], SlicerConfig())
        self.assertEqual(code.count("; retract"), 1)
        self.assertNotIn("; unretract", code)

    def test_neighboring_infill_across_periodic_seam_does_not_retract(self):
        turn = 2 * math.pi
        layer = LayerToolpaths(0, 30.0, {"INFILL_AXIAL": [
            Toolpath("INFILL_AXIAL", 30.0, [(turn - 0.02, 0), (turn - 0.02, 5)]),
            Toolpath("INFILL_AXIAL", 30.0, [(0.02, 0), (0.02, 5)]),
        ]}, travel_region=box(turn - 0.04, 0, turn, 5).union(box(0, 0, 0.04, 5)))
        code = emit_gcode([layer], SlicerConfig(infill_spacing=1.2))
        self.assertEqual(code.count("; retract"), 1)
        self.assertNotIn("; unretract", code)

    def test_paths_within_role_are_ordered_from_current_position(self):
        config = SlicerConfig(axial_min=0.0, max_axial_segment=20.0)
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "EXTERNAL_PERIMETER": [
                    Toolpath("EXTERNAL_PERIMETER", 10.0, [(0.0, 10.0), (0.1, 10.0)]),
                    Toolpath("EXTERNAL_PERIMETER", 10.0, [(0.0, 0.0), (0.1, 0.0)]),
                ]
            },
        )

        segments = [
            segment
            for segment in parse_extrusion_segments(emit_gcode([layer], config))
            if segment[1] == "EXTERNAL_PERIMETER"
        ]

        self.assertTrue(segments)
        self.assertEqual(segments[0][2], 0.0)

    def test_short_perimeter_segments_are_emitted(self):
        config = SlicerConfig(axial_min=0.0, min_path_length=0.6, max_axial_segment=20.0, retract_mm=0.0)
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "EXTERNAL_PERIMETER": [
                    Toolpath("EXTERNAL_PERIMETER", 10.0, [(0.0, 0.0), (0.01, 0.0)]),
                ],
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.1, 0.0), (0.11, 0.0)]),
                ],
            },
        )

        segments = parse_extrusion_segments(emit_gcode([layer], config))

        self.assertTrue(any(segment[1] == "EXTERNAL_PERIMETER" for segment in segments))
        self.assertFalse(any(segment[1] == "SUPPORT" for segment in segments))

    def test_support_meander_preserves_reverse_radial_leg(self):
        config = SlicerConfig(max_axial_segment=20.0)
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={
                "SUPPORT": [
                    Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
                ]
            },
        )

        segments = [segment for segment in parse_extrusion_segments(emit_gcode([layer], config)) if segment[1] == "SUPPORT"]

        self.assertTrue(any(y1 < y0 for _, _, _, y0, _, y1 in segments))
        self.assertTrue(any(abs(y1 - y0) < 1e-6 and abs(x1 - x0) > 0.5 for _, _, x0, y0, x1, y1 in segments))

    def test_emitted_x_values_are_shifted_by_axial_min(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={"SUPPORT": [Toolpath("SUPPORT", 10.0, [(0.0, 0.0), (0.0, 10.0)])]},
        )

        default_gcode = emit_gcode([layer], SlicerConfig(max_axial_segment=20.0))
        shifted_gcode = emit_gcode([layer], SlicerConfig(axial_min=2.5, max_axial_segment=20.0))

        self.assertIn("G0 X15.000 Y0.000", default_gcode)
        self.assertIn("G1 X25.000 Y0.000", default_gcode)
        self.assertIn("G0 X2.500 Y0.000", shifted_gcode)
        self.assertIn("G1 X12.500 Y0.000", shifted_gcode)

    def test_axial_min_sets_minimum_emitted_x_for_negative_model_coordinates(self):
        layer = LayerToolpaths(
            layer_index=0,
            radius=10.0,
            paths_by_role={"SUPPORT": [Toolpath("SUPPORT", 10.0, [(0.0, -16.0), (0.0, 4.0)])]},
        )

        gcode = emit_gcode([layer], SlicerConfig(axial_min=30.0, max_axial_segment=30.0))

        x_values = [float(value) for value in re.findall(r"\bX(-?\d+(?:\.\d+)?)", gcode)]
        self.assertEqual(min(x_values), 30.0)
        self.assertEqual(max(x_values), 50.0)

    def test_emitted_segments_follow_planned_geometry(self):
        for fixture in ("helix_like_ribbon", "tube"):
            with self.subTest(fixture=fixture):
                config = SlicerConfig(theta_step_deg=20, axial_step=2.0, layer_height=1.0)
                layers = build_toolpaths(slice_mesh(make_fixture(fixture), config, max_layers=2), config)
                offset = config.axial_min - min(p[1] for layer in layers for paths in layer.paths_by_role.values()
                                                for path in paths for p in path.points)
                sources = {
                    (layer.layer_index, role): (layer.radius, unary_union([
                        LineString([(t * layer.radius, x + offset) for t, x in path.points]) for path in paths
                    ]).buffer(0.002))
                    for layer in layers for role, paths in layer.paths_by_role.items()
                }
                for layer, role, x0, y0, x1, y1 in parse_extrusion_segments(emit_gcode(layers, config)):
                    radius, source = sources[layer, role]
                    t0, t1 = math.radians(y0), math.radians(y1)
                    period = 2 * math.pi * radius
                    shift = round((source.centroid.x - (t0 + t1) * radius / 2) / period) * period
                    candidates = [LineString([(t0 * radius + shift + extra, x0),
                                               (t1 * radius + shift + extra, x1)])
                                  for extra in (-period, 0.0, period)]
                    self.assertTrue(any(source.covers(segment) for segment in candidates),
                                    (fixture, layer, role, x0, y0, x1, y1))

    def test_viewer_format_comments_and_axes_are_present(self):
        gcode = slice_to_gcode(make_fixture("tube"), SlicerConfig(theta_step_deg=15, axial_step=2.0, layer_height=1.0), max_layers=1)
        self.assertIn(";Z:", gcode)
        self.assertRegex(gcode, r"(?m)^G0 Z", msg="expected radial Z move")
        self.assertRegex(gcode, r"G1 X[-0-9.]+ Y[-0-9.]+ E", msg="expected X/Y/E extrusion move")

    def test_cubeandring_layers_64_to_66_have_no_full_turn_perimeter_segments(self):
        root = Path(__file__).resolve().parents[2]
        gcode = slice_to_gcode(load_mesh(root / "cubeandring.stl"), SlicerConfig(), max_layers=67)
        bad = []
        for layer, role, x0, y0, x1, y1 in parse_extrusion_segments(gcode):
            if layer in {64, 65, 66} and role in {"EXTERNAL_PERIMETER", "INTERNAL_PERIMETER"}:
                if abs(y1 - y0) >= 300.0:
                    bad.append((layer, role, x0, y0, x1, y1))
        self.assertEqual(bad, [])

    def test_single_rotated_cube_final_layer_does_not_spin_for_tiny_triangle(self):
        root = Path(__file__).resolve().parents[2]
        config = SlicerConfig(layer_backend="intersection", layer_height=0.25, bed_radius=30.0, slice_workers=1, toolpath_workers=1)
        gcode = slice_to_gcode(load_mesh(root / "single_rotated_cube.stl"), config)
        y_values = []
        in_final_layer = False
        for line in gcode.splitlines():
            if line.startswith(";LAYER:"):
                in_final_layer = line == ";LAYER:134"
                continue
            if not in_final_layer or not line.startswith(("G0 ", "G1 ")):
                continue
            match = re.search(r"\bY(-?\d+(?:\.\d+)?)", line)
            if match:
                y_values.append(float(match.group(1)))
        if y_values:
            self.assertLess(max(y_values) - min(y_values), 5.0)

    def test_rotatecube_debug_gcode_does_not_emit_periodic_seam_wall_segments(self):
        mesh_path = Path(__file__).resolve().parents[2].parent / "rotatecube2-debugging.stl"
        if not mesh_path.exists():
            self.skipTest(f"{mesh_path} not present")
        config = SlicerConfig(
            layer_backend="intersection",
            layer_height=0.25,
            bed_radius=30.0,
            slice_workers=1,
            toolpath_workers=1,
            intersection_refine_factor=0,
        )
        mesh = load_mesh(mesh_path)
        radius = layer_radii(mesh, config)[17]
        occupied, window, full = sample_intersection_layer(mesh, radius, config)
        # Only layer 17's contour matters here; avoid slicing earlier layers or
        # generating solid infill, support and gap rafts for the whole model.
        empty = GeometryCollection()
        layer = LayerSlice(
            index=17, radius=radius, occupied=occupied,
            top_surfaces=empty, bottom_surfaces=empty, needs_support=empty, support=empty,
            theta_window=window, is_full_wrap=full,
        )
        self.assertEqual(window, (0.0, 2.0 * math.pi))
        for theta in window:
            seam = LineString([(theta, occupied.bounds[1]), (theta, occupied.bounds[3])])
            self.assertGreater(occupied.intersection(seam).length, 0.0)
        gcode = emit_gcode([build_layer_toolpaths(layer, config)], config)
        perimeters = [s for s in parse_extrusion_segments(gcode) if s[0] == 17 and s[1] == "EXTERNAL_PERIMETER"]
        self.assertTrue(perimeters, "Expected external perimeter extrusion in the target layer")

        seam_segments = []
        for _, _, x0, y0, x1, y1 in perimeters:
            a = math.radians(y0 % 360.0)
            b = math.radians(y1 % 360.0)
            if min(a, 2.0 * math.pi - a) < 1e-4 and min(b, 2.0 * math.pi - b) < 1e-4:
                seam_segments.append((x0, y0, x1, y1))

        self.assertEqual(seam_segments, [])


if __name__ == "__main__":
    unittest.main()
