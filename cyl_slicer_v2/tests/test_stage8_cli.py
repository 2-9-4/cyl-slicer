from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


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


if __name__ == "__main__":
    unittest.main()
