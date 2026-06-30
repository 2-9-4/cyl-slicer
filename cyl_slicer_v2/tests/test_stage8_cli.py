from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from cyl_slicer_v2.cli import build_parser


class CliTests(unittest.TestCase):
    def run_cli(self, *args: str) -> str:
        root = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "out.gcode"
            cmd = [
                sys.executable,
                "-m",
                "cyl_slicer_v2.cli",
                *args,
                "--out",
                str(out),
                "--max-layers",
                "1",
                "--theta-step-deg",
                "20",
                "--axial-step",
                "2.0",
            ]
            subprocess.run(cmd, cwd=root, check=True, capture_output=True, text=True)
            return out.read_text()

    def test_cli_slices_generated_fixtures(self):
        for name in ("tube", "offset_box_pair", "helix_like_ribbon"):
            with self.subTest(name=name):
                gcode = self.run_cli("--fixture", name)
                self.assertIn(";start of slice", gcode)
                self.assertIn(";end of slice", gcode)

    def test_cli_smoke_slices_existing_stls(self):
        root = Path(__file__).resolve().parents[2]
        for rel in ("tube.stl", "cubeandring.stl", "helix.stl"):
            with self.subTest(stl=rel):
                gcode = self.run_cli(str(root / rel))
                self.assertIn(";start of slice", gcode)
                self.assertIn(";end of slice", gcode)

    def test_cli_accepts_build_cylinder_radius(self):
        gcode = self.run_cli("--fixture", "tube", "--build-cylinder-radius", "8.0")
        self.assertIn(";LAYER:-", gcode)
        self.assertIn(";TYPE:SUPPORT", gcode)

    def test_cli_accepts_top_and_bottom_layers(self):
        args = build_parser().parse_args(
            [
                "--fixture",
                "tube",
                "--top-layers",
                "4",
                "--bottom-layers",
                "2",
                "--filament-diameter",
                "2.85",
                "--raft-layers",
                "3",
                "--raft-grid-spacing",
                "2.5",
                "--raft-axial-gap",
                "1.2",
                "--support-gap-layers",
                "2",
                "--support-interface-layers",
                "3",
                "--support-interface-spacing",
                "0.8",
                "--support-meander-max-angle-deg",
                "12.5",
                "--support-cross-spacing",
                "6.0",
                "--min-path-length",
                "0.7",
                "--retract-mm",
                "0.9",
                "--z-hop-height",
                "0.6",
                "--z-hop-min-travel",
                "3.5",
                "--no-z-hop-perimeter-crossing",
                "--axial-min",
                "12.5",
            ]
        )

        self.assertEqual(args.top_layers, 4)
        self.assertEqual(args.bottom_layers, 2)
        self.assertEqual(args.filament_diameter, 2.85)
        self.assertEqual(args.raft_layers, 3)
        self.assertEqual(args.raft_grid_spacing, 2.5)
        self.assertEqual(args.raft_axial_gap, 1.2)
        self.assertEqual(args.support_gap_layers, 2)
        self.assertEqual(args.support_interface_layers, 3)
        self.assertEqual(args.support_interface_spacing, 0.8)
        self.assertEqual(args.support_meander_max_angle_deg, 12.5)
        self.assertEqual(args.support_cross_spacing, 6.0)
        self.assertEqual(args.min_path_length, 0.7)
        self.assertEqual(args.retract_mm, 0.9)
        self.assertEqual(args.z_hop_height, 0.6)
        self.assertEqual(args.z_hop_min_travel, 3.5)
        self.assertFalse(args.z_hop_on_perimeter_crossing)
        self.assertEqual(args.axial_min, 12.5)


if __name__ == "__main__":
    unittest.main()
