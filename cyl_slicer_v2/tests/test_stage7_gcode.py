from __future__ import annotations

from pathlib import Path
import math
import re
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.fixtures import make_fixture
from cyl_slicer_v2.gcode import parse_extrusion_segments
from cyl_slicer_v2.mesh_io import load_mesh
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
            if line.startswith("G1 "):
                match = re.search(r"\bE(-?\d+(?:\.\d+)?)", line)
                if match:
                    values.append(float(match.group(1)))
        self.assertEqual(values, sorted(values))

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

    def test_no_extrusion_segment_moves_backward_in_y(self):
        gcode = slice_to_gcode(make_fixture("helix_like_ribbon"), SlicerConfig(theta_step_deg=20, axial_step=2.0, layer_height=1.0), max_layers=2)
        segments = parse_extrusion_segments(gcode)
        self.assertTrue(segments)
        self.assertTrue(all(end_y >= start_y for _, _, _, start_y, _, end_y in segments))

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
        self.assertTrue(y_values)
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
