from __future__ import annotations

import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _effective_material_fan_pwm, _layer_fan_pwm, emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


class MaterialFanPolicyTests(unittest.TestCase):
    def test_quality_and_legacy_defaults_remain_manual(self):
        self.assertEqual(SlicerConfig().material_fan_policy, "manual")
        for profile in ("quality", "legacy-cli", "library"):
            args = build_parser().parse_args(["--profile", profile])
            self.assertEqual(_profile_value(args, "material_fan_policy"), "manual")

    def test_preset_ceilings_never_increase_user_fan(self):
        caps = {"pla":255, "unprofiled":255, "petg":180, "tpu":180, "abs":0}
        for material, cap in caps.items():
            for user in (0, 80, 200, 255):
                c = SlicerConfig(material_profile=material, material_fan_policy="preset",
                                 fan_speed=user)
                self.assertEqual(_effective_material_fan_pwm(c), min(cap, user))

    def test_manual_policy_ignores_material(self):
        for material in ("pla", "petg", "abs", "tpu"):
            c = SlicerConfig(material_profile=material, fan_speed=217)
            self.assertEqual(_effective_material_fan_pwm(c), 217)
            self.assertEqual(_layer_fan_pwm(4, c), 217)

    def test_ramp_uses_capped_pwm_and_start_height(self):
        c = SlicerConfig(material_profile="petg", material_fan_policy="preset",
                         fan_speed=255, fan_start_height=1, fan_ramp_height_mm=2)
        self.assertEqual(_layer_fan_pwm(.5, c), 0)
        self.assertEqual(_layer_fan_pwm(1.5, c), 45)
        self.assertEqual(_layer_fan_pwm(3.5, c), 180)

    def test_bridge_override_restores_material_capped_fan(self):
        c = SlicerConfig(material_profile="petg", material_fan_policy="preset",
                         fan_speed=255, bridge_mode="print", bridge_fan_pwm=230,
                         bed_radius=30, retract_mm=0)
        p = Toolpath("BRIDGE", 32, [(0,0),(0,20)])
        g = emit_gcode([LayerToolpaths(0,32,{"BRIDGE":[p]})], c)
        self.assertIn("M106 S180 ; part cooling on", g)
        self.assertIn("M106 S230 ; bridge fan override", g)
        self.assertIn("M106 S180 ; restore layer fan after bridge", g)

    def test_cli_opt_in_and_validation(self):
        args = build_parser().parse_args(["--material-profile", "petg",
                                          "--material-fan-policy", "preset",
                                          "--fan-speed", "125"])
        self.assertEqual(_profile_value(args, "material_fan_policy"), "preset")
        self.assertEqual(args.material_profile, "petg")
        with self.assertRaisesRegex(ValueError, "material_fan_policy"):
            SlicerConfig(material_fan_policy="magic").validate()


if __name__ == "__main__":
    unittest.main()
