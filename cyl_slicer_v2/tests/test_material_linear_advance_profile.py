from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from cyl_slicer_v2.cli import _selected_material_advance_k, build_parser, main


class MaterialScopedLinearAdvanceTests(unittest.TestCase):
    def test_only_selected_material_uses_a_user_supplied_k(self):
        parser = build_parser()
        args = parser.parse_args([
            "--material-profile", "petg",
            "--linear-advance-material-k", "pla=.04",
            "--linear-advance-material-k", "petg=.12"])
        self.assertAlmostEqual(_selected_material_advance_k(args, parser), .12)
        args = parser.parse_args(["--material-profile","abs",
            "--linear-advance-material-k","petg=.12"])
        self.assertIsNone(_selected_material_advance_k(args, parser))

    def test_rejects_bad_or_duplicated_calibration_records(self):
        parser = build_parser()
        for entry in ("petg", "unknown=.12", "pla=nan", "pla=inf", "pla=-1",
                      "pla=11", "pla=abc"):
            with self.subTest(entry=entry), self.assertRaises(SystemExit):
                _selected_material_advance_k(parser.parse_args([
                    "--linear-advance-material-k",entry]), parser)
        args = parser.parse_args(["--linear-advance-material-k","pla=.1",
                                   "--linear-advance-material-k","pla=.2"])
        with self.assertRaises(SystemExit):
            _selected_material_advance_k(args,parser)

    def test_selected_k_needs_explicit_firmware(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)/"sample.gcode"
            with self.assertRaises(SystemExit):
                main(["--fixture","tube","--max-layers","1","--out",str(output),
                      "--material-profile","pla","--linear-advance-material-k","pla=.13"])
            self.assertFalse(output.exists())

    def test_k_is_scoped_to_matching_material_and_changes_no_motion(self):
        with tempfile.TemporaryDirectory() as tmp:
            def make(name, extra):
                filename = Path(tmp)/(name+".gcode")
                self.assertEqual(main(["--fixture", "tube", "--max-layers", "1",
                                       "--out",str(filename), *extra]), 0)
                return filename.read_text()
            base = ["--material-profile","petg"]
            original = make("base", base)
            selected = make("chosen", [*base, "--marlin-linear-advance",
                       "--linear-advance-material-k","pla=.05",
                       "--linear-advance-material-k","petg=.16"])
            self.assertIn(";LINEAR_ADVANCE:Marlin M900 K0.16",selected)
            self.assertIn("M900 K0.16",selected)
            self.assertIn("M900 K0 ; clear job-scoped linear advance",selected)
            self.assertNotIn("\nM900 K",original)
            active_moves = [line for line in selected.splitlines()
                            if line.startswith(("G0 ","G1 ","G92 "))]
            default_moves = [line for line in original.splitlines()
                            if line.startswith(("G0 ","G1 ","G92 "))]
            self.assertEqual(active_moves, default_moves)
            unmatched = make("none", ["--material-profile","abs",
                   "--marlin-linear-advance",
                   "--linear-advance-material-k","petg=.16"])
            self.assertNotIn("\nM900 K",unmatched)
            explicit = make("override", [*base, "--marlin-linear-advance",
                   "--linear-advance-material-k","petg=.16",
                   "--linear-advance-k","0.02"])
            self.assertIn("M900 K0.02",explicit)
            self.assertNotIn("M900 K0.16",explicit)


if __name__=="__main__":
    unittest.main()
