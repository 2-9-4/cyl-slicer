from __future__ import annotations

from dataclasses import replace
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _role_fan_override, emit_gcode
from cyl_slicer_v2.models import LayerToolpaths, Toolpath


def fixture():
    r = 32
    paths = {}
    for i, role in enumerate(("SUPPORT", "SUPPORT_INTERFACE", "SOLID_INFILL", "BRIDGE")):
        paths[role] = [Toolpath(role, r, [(i*.2, 0), (i*.2, 10)])]
    return LayerToolpaths(0, r, paths)


def command_lines(code):
    return [line for line in code.splitlines()
            if line.startswith(("G0 ", "G1 ", "G92 ", "M201 ", "M204 "))]


class SupportRoleFanTests(unittest.TestCase):
    def config(self, **kwargs):
        return SlicerConfig(bed_radius=30, retract_mm=0,
                            cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0,
                            fan_start_height=0, fan_ramp_height_mm=0,
                            **kwargs)

    def test_default_preserves_all_motion_and_m106_events(self):
        base = self.config(bridge_mode="off")
        modified = replace(base, support_fan_pwm=None, support_interface_fan_pwm=None)
        self.assertEqual(emit_gcode([fixture()], base), emit_gcode([fixture()], modified))
        self.assertEqual(_role_fan_override("SUPPORT", base)[0], None)

    def test_support_body_and_interface_have_independent_pwm_and_restore(self):
        settings = self.config(fan_speed=200,
                               support_fan_pwm=0, support_interface_fan_pwm=255)
        code = emit_gcode([fixture()], settings)
        self.assertIn(";SUPPORT_FAN_PWM:0", code)
        self.assertIn(";SUPPORT_INTERFACE_FAN_PWM:255", code)
        self.assertIn("M106 S0 ; support fan override", code)
        self.assertIn("M106 S255 ; support_interface fan override", code)
        self.assertIn("M106 S200 ; restore layer fan after support", code)
        self.assertIn("M106 S200 ; restore layer fan after support_interface", code)
        self.assertEqual(command_lines(code),
                         command_lines(emit_gcode([fixture()], replace(
                             settings, support_fan_pwm=None,
                             support_interface_fan_pwm=None))))

    def test_matching_pwm_does_not_repeat_m106(self):
        cfg = self.config(fan_speed=100, support_fan_pwm=100,
                          support_interface_fan_pwm=100)
        code = emit_gcode([fixture()], cfg)
        self.assertNotIn("fan override", code)

    def test_material_preset_caps_support_not_explicit_bridge(self):
        c = self.config(material_profile="petg", material_fan_policy="preset",
                        support_fan_pwm=255, support_interface_fan_pwm=240,
                        bridge_mode="print", bridge_fan_pwm=230)
        self.assertEqual(_role_fan_override("SUPPORT", c)[0], 180)
        self.assertEqual(_role_fan_override("SUPPORT_INTERFACE", c)[0], 180)
        self.assertEqual(_role_fan_override("BRIDGE", c)[0], 230)
        code = emit_gcode([fixture()], c)
        self.assertIn("M106 S230 ; bridge fan override", code)
        abs_cfg = replace(c, material_profile="abs")
        self.assertEqual(_role_fan_override("SUPPORT", abs_cfg)[0], 0)

    def test_validation_and_cli_profiles(self):
        for name in ("support_fan_pwm", "support_interface_fan_pwm"):
            for bad in (-1, 256, True, .5, "255"):
                with self.subTest(name=name, bad=bad), self.assertRaisesRegex(ValueError, name):
                    self.config(**{name: bad}).validate()
        for profile in ("quality","legacy-cli","library"):
            args = build_parser().parse_args(["--profile", profile])
            self.assertIsNone(_profile_value(args, "support_fan_pwm"))
            self.assertIsNone(_profile_value(args, "support_interface_fan_pwm"))
        args = build_parser().parse_args(["--support-fan-pwm", "0",
                                          "--support-interface-fan-pwm", "185"])
        self.assertEqual(_profile_value(args, "support_fan_pwm"), 0)
        self.assertEqual(_profile_value(args, "support_interface_fan_pwm"), 185)


if __name__ == "__main__":
    unittest.main()
