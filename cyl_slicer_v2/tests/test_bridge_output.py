from __future__ import annotations

import math
import unittest

from shapely.geometry import GeometryCollection, box
from shapely.ops import unary_union

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _role_print_settings, emit_gcode
from cyl_slicer_v2.layers import _assign_printable_bridges
from cyl_slicer_v2.models import LayerSlice, LayerToolpaths, Toolpath
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


EMPTY = GeometryCollection()


def samples():
    whole = box(0,0,1,10)
    underneath = unary_union([box(0,0,1,4),box(0,6,1,10)])
    gap = box(0,4,1,6)
    lower = LayerSlice(0,30.08,underneath,EMPTY,EMPTY,EMPTY,EMPTY)
    upper = LayerSlice(1,30.24,whole,EMPTY,gap,gap,EMPTY,
                       support_interface_target=gap)
    return lower,upper


class BridgeOutputTests(unittest.TestCase):
    def test_bridge_mode_is_explicit_and_default_off(self):
        c = SlicerConfig()
        self.assertEqual(c.bridge_mode,"off")
        lower,upper = samples()
        _assign_printable_bridges([lower,upper],c)
        self.assertIsNone(upper.bridge_regions)
        self.assertGreater(upper.needs_support.area,0)

    def test_support_only_removed_under_confirmed_two_bearing_gap(self):
        config = SlicerConfig(bridge_mode="print",wall_count=2)
        lower,upper = samples()
        original = upper.needs_support.area
        _assign_printable_bridges([lower,upper],config)
        self.assertIsNotNone(upper.bridge_regions)
        self.assertGreater(upper.bridge_regions["axial"].area,0)
        self.assertTrue(upper.bridge_regions["radial"].is_empty)
        self.assertLess(upper.needs_support.area,original)
        self.assertLess(upper.support_interface_target.area,original)
        # All printed bridges are still a subset of the actual model slab.
        self.assertTrue(upper.occupied.covers(upper.bridge_regions["axial"]))

    def test_bridge_infill_oriented_to_two_anchors_and_not_double_printed(self):
        config = SlicerConfig(bridge_mode="print",wall_count=2)
        lower,upper = samples()
        _assign_printable_bridges([lower,upper],config)
        roles = build_layer_toolpaths(upper,config).paths_by_role
        bridge = roles.get("BRIDGE",[])
        self.assertTrue(bridge)
        for p in bridge:
            # Axial bridge strands span the support gap along axial X.
            self.assertAlmostEqual(p.points[0][0],p.points[-1][0],places=7)
            self.assertNotAlmostEqual(p.points[0][1],p.points[-1][1],places=7)
        # The bridge-specific bottom surface is excluded from SOLID_INFILL.
        # Normal top skin or other model walls may still coexist, but must not
        # print the same inner corridor centerline.
        for solid in roles.get("SOLID_INFILL",[]) + roles.get("TOP_SKIN",[]):
            self.assertFalse(
                any(abs(a[0]-b[0]) < 1e-8 and
                    abs(a[1]-b[1]) < 1e-8
                    for bridge_path in bridge
                    for a in solid.points for b in bridge_path.points)
            )

    def test_physical_bridge_feed_flow_and_temporary_fan_override(self):
        c = SlicerConfig(bridge_mode="print",bridge_fan_pwm=220,
            bridge_print_feed_mm_per_min=900, bridge_flow_ratio=.9,
            retract_mm=0,bed_radius=30)
        self.assertEqual(_role_print_settings("BRIDGE",c),(900,.9))
        bridge = Toolpath("BRIDGE",30.24,[(.4,4.2),(.4,5.8)])
        lp = LayerToolpaths(1,30.24,{"BRIDGE":[bridge]})
        code = emit_gcode([lp],c)
        self.assertIn(";TYPE:BRIDGE",code)
        self.assertIn("M106 S220 ; bridge fan override",code)
        self.assertIn("M106 S0 ; restore layer fan after bridge",code)
        extrusions=[line for line in code.splitlines() if line.startswith("G1 X") and " E" in line]
        self.assertEqual(len(extrusions),1)
        self.assertIn("F900.0",extrusions[0])
        filament=float(extrusions[0].split(" E")[1].split()[0])
        self.assertAlmostEqual(filament,1.6*c.extrusion_per_mm*.9,delta=1e-4)

    def test_long_span_still_supported(self):
        c=SlicerConfig(bridge_mode="print",bridge_max_span_mm=.5)
        lower,upper=samples()
        _assign_printable_bridges([lower,upper],c)
        self.assertIsNone(upper.bridge_regions)
        self.assertGreater(upper.needs_support.area,0)

    def test_invalid_bridge_parameters(self):
        for kwargs in ({"bridge_mode":"all"},
            {"bridge_max_span_mm":0},{"bridge_anchor_overlap_mm":-1},
            {"bridge_print_feed_mm_per_min":0},
            {"bridge_flow_ratio":math.nan},
            {"bridge_fan_pwm":256}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                SlicerConfig(**kwargs).validate()


if __name__=="__main__":
    unittest.main()
