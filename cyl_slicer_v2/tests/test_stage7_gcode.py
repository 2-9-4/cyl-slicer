from __future__ import annotations

from pathlib import Path
import math
import re
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.gcode import emit_gcode, parse_extrusion_segments
from cyl_slicer_v2.mesh_io import load_mesh
from cyl_slicer_v2.models import LayerToolpaths, Toolpath
from cyl_slicer_v2.pipeline import slice_to_gcode


class GcodeTests(unittest.TestCase):
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

        self.assertIn("G0 Z10.5000 F3000.0 ; z hop", gcode)
        self.assertIn("G0 Z10.0000 F3000.0 ; z unhop", gcode)
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
            SlicerConfig(z_hop_height=0.5, z_hop_min_travel=100.0, max_axial_segment=20.0, retract_mm=0.0),
        )

        self.assertEqual(gcode.count("; z hop"), 1)
        self.assertIn("G0 Z10.5000 F3000.0 ; z hop\nG0 X16.000 Y40.107", gcode)

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

    def test_model_extrusion_segments_do_not_move_backward_in_y(self):
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

    def test_no_thread_jump_like_segments(self):
        for fixture in ("helix_like_ribbon", "tube"):
            with self.subTest(fixture=fixture):
                gcode = slice_to_gcode(make_fixture(fixture), SlicerConfig(theta_step_deg=20, axial_step=2.0, layer_height=1.0), max_layers=2)
                flags = []
                for _, role, x0, y0, x1, y1 in parse_extrusion_segments(gcode):
                    if abs(x1 - x0) > 2.0 and abs(y1 - y0) < 30.0:
                        flags.append(role)
                self.assertEqual(flags, [])

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
        gcode = slice_to_gcode(load_mesh(mesh_path), config, max_layers=18)

        seam_segments = []
        for layer, role, x0, y0, x1, y1 in parse_extrusion_segments(gcode):
            if layer != 17 or role != "EXTERNAL_PERIMETER":
                continue
            a = math.radians(y0 % 360.0)
            b = math.radians(y1 % 360.0)
            if min(a, 2.0 * math.pi - a) < 1e-4 and min(b, 2.0 * math.pi - b) < 1e-4:
                seam_segments.append((x0, y0, x1, y1))

        self.assertEqual(seam_segments, [])


if __name__ == "__main__":
    unittest.main()
