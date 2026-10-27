from __future__ import annotations

import math
import unittest

from cyl_slicer_v2.cli import CLI_PROFILE_OVERRIDES, _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode


QUALITY_EXPECTED = {
    "layer_backend": "intersection",
    "bed_radius": 30,
    "nozzle_diameter": 0.4,
    "filament_diameter": 1.75,
    "layer_height": 0.16,
    "nozzle_width": 0.42,
    "bead_model": "rounded",
    "wall_count": 4,
    "wall_order": "outer-first",
    "top_layers": 5,
    "bottom_layers": 5,
    "infill_spacing": 1.5,
    "sparse_infill_pattern": "alternate",
    "solid_infill_pattern": "alternate",
    "monotonic_top_skin": False,
    "min_path_length": 0.1,
    "intersection_curve_error_mm": 0.01,
    "support_min_island_area": 0.25,
    "support_min_region_area": 1.0,
    "support_bearing_min_area": 1.0,
    "support_merge_gap": 2.0,
    "raft_layers": 4,
    "raft_grid_spacing": None,
    "raft_axial_gap": 1.0,
    "raft_merge_gap_mm": 2.0,
    "support_interface_layers": 3,
    "support_interface_spacing": 0.4,
    "support_gap_layers": 1,
    "feed_mm_per_min": 2700,
    "travel_mm_per_min": 9000,
    "cooling_min_layer_time_s": 8,
    "cooling_min_print_speed_mm_per_min": 300,
    "fan_speed": 255,
    "fan_start_height": 1,
    "fan_ramp_height_mm": 1,
    "hotend_temp": 195,
    "material_profile": "pla",
    "retract_speed_mm_per_min": 1200,  # 20 mm/s CLI
    "retract_mm": 1.4,
}


class QualityDefaultTests(unittest.TestCase):
    def test_python_defaults_are_exact_requested_quality_profile(self):
        config = SlicerConfig()
        config.validate()
        for field, expected in QUALITY_EXPECTED.items():
            with self.subTest(field=field):
                self.assertEqual(getattr(config, field), expected)

    def test_cli_quality_defaults_agree_with_python(self):
        args = build_parser().parse_args([])
        self.assertEqual(args.profile, "quality")
        for field, expected in QUALITY_EXPECTED.items():
            with self.subTest(field=field):
                resolved = (args.filament_diameter if field == "filament_diameter"
                            else _profile_value(args, field))
                self.assertEqual(resolved, expected)

    def test_legacy_cli_is_complete_pinned_old_recipe(self):
        args = build_parser().parse_args(["--profile", "legacy-cli"])
        for field, expected in CLI_PROFILE_OVERRIDES.items():
            with self.subTest(field=field):
                self.assertEqual(_profile_value(args, field), expected)
        self.assertFalse(_profile_value(args, "monotonic_top_skin"))
        self.assertEqual(_profile_value(args, "layer_backend"), "raster")

    def test_library_alias_matches_quality(self):
        args = build_parser().parse_args(["--profile", "library"])
        for field, expected in QUALITY_EXPECTED.items():
            if field == "filament_diameter":
                continue
            with self.subTest(field=field):
                self.assertEqual(_profile_value(args, field), expected)

    def test_explicit_overrides_work_in_both_profiles_and_allow_false(self):
        options = ["--layer-backend", "raster", "--layer-height", "0.2",
                   "--line-width", "0.5", "--nozzle-diameter", "0.6",
                   "--build-cylinder-radius", "42", "--print-feed", "1234",
                   "--travel-feed", "4321", "--retract-speed", "12",
                   "--retract-mm", "0", "--no-monotonic-top-skin",
                   "--fan-speed", "0", "--cooling-min-layer-time-s", "0",
                   "--sparse-infill-pattern", "both", "--raft-merge-gap-mm", "0",
                   "--support-merge-gap", "0"]
        want = dict(layer_backend="raster", layer_height=.2, nozzle_width=.5,
                    nozzle_diameter=.6, bed_radius=42, feed_mm_per_min=1234,
                    travel_mm_per_min=4321, retract_speed_mm_per_min=720,
                    retract_mm=0, monotonic_top_skin=False, fan_speed=0,
                    cooling_min_layer_time_s=0, sparse_infill_pattern="both",
                    raft_merge_gap_mm=0, support_merge_gap=0)
        for profile in ("quality", "legacy-cli", "library"):
            args = build_parser().parse_args(["--profile", profile, *options])
            for field, value in want.items():
                with self.subTest(profile=profile, field=field):
                    self.assertEqual(_profile_value(args, field), value)

    def test_default_no_unsupported_pressure_advance_and_geometry_is_valid(self):
        c = SlicerConfig()
        self.assertIsNone(c.linear_advance_k)
        self.assertIsNone(c.linear_advance_firmware)
        self.assertFalse(c.intersection_strict_recovery)
        self.assertEqual(c.seam_policy, "nearest")
        self.assertGreater(c.bead_area, 0)
        self.assertGreater(c.dense_line_spacing, 0)
        c.validate()
        code = emit_gcode([], c)
        self.assertIn(";LAYER_HEIGHT:0.16", code)
        self.assertIn(";BEAD_MODEL:rounded", code)
        self.assertIn(";LINE_WIDTH:0.42", code)
        self.assertIn(";MATERIAL_PROFILE:pla", code)
        self.assertIn(";NOZZLE_DIAMETER:0.4", code)
        self.assertNotIn("\nM900 K", code)

    def test_monotonic_top_skin_remains_an_explicit_override(self):
        for profile in ("quality", "legacy-cli", "library"):
            for flag, expected in (("--monotonic-top-skin", True),
                                   ("--no-monotonic-top-skin", False)):
                with self.subTest(profile=profile, flag=flag):
                    args = build_parser().parse_args(["--profile", profile, flag])
                    self.assertEqual(_profile_value(args, "monotonic_top_skin"), expected)


if __name__ == "__main__":
    unittest.main()
