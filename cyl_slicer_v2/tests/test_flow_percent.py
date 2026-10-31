import unittest
from dataclasses import replace

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _cooling_feed_for_target, emit_toolpath
from cyl_slicer_v2.models import GcodeState, Toolpath


class FlowPercentTests(unittest.TestCase):
    def test_cli_defaults_overrides_and_validation(self):
        for profile in ('quality', 'legacy-cli'):
            args = build_parser().parse_args(['--profile', profile])
            self.assertEqual(_profile_value(args, 'flow_percent'), 100)
            args = build_parser().parse_args(['--profile', profile, '--flow-percent', '98'])
            self.assertEqual(_profile_value(args, 'flow_percent'), 98)
        for value in (0, -1, float('nan'), float('inf')):
            with self.assertRaisesRegex(ValueError, 'flow_percent'):
                SlicerConfig(flow_percent=value).validate()

    def test_deposition_scales_for_roles_and_layer_thickness(self):
        base = SlicerConfig(flow_percent=98, extrusion_per_mm=.1,
                            support_flow_ratio=.8, bridge_flow_ratio=.9,
                            extrusion_feed_ramp_mm=3, retract_mm=1.4)
        for height in (.16, .2):
            config = base.for_deposited_thickness(height)
            for role, ratio in (('EXTERNAL_PERIMETER', 1), ('SUPPORT', .8), ('BRIDGE', .9)):
                with self.subTest(height=height, role=role):
                    state = GcodeState(radial=30.3)
                    path = Toolpath(role, 30.1, [(0, 0), (0, 10)])
                    lines = emit_toolpath(path, state, config, 15, retract_after=True)
                    self.assertAlmostEqual(state.extrusion + config.retract_mm,
                                           10 * config.extrusion_per_mm * .98 * ratio, places=9)
                    self.assertTrue(any('; retract' in line for line in lines))

    def test_flow_limit_and_cooling_use_scaled_deposition(self):
        base = SlicerConfig(extrusion_per_mm=.1, flow_percent=100,
                            max_volumetric_flow_mm3_s=1, feed_mm_per_min=3000)
        low = replace(base, flow_percent=50)
        _, base_seconds = _cooling_feed_for_target({('INFILL', 1): 100}, base, 0)
        _, low_seconds = _cooling_feed_for_target({('INFILL', 1): 100}, low, 0)
        self.assertGreater(base_seconds, low_seconds)
        self.assertAlmostEqual(base_seconds / low_seconds, 2)
