from __future__ import annotations

from dataclasses import replace
import math
import unittest

from shapely.geometry import GeometryCollection, Point, box
from shapely.ops import unary_union

from cyl_slicer_v2.cli import build_parser, _profile_value
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _role_print_settings, emit_gcode
from cyl_slicer_v2.layers import add_bottom_support_interfaces, propagate_support
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


EMPTY = GeometryCollection()
TURN = 2*math.pi


def stack():
    model = box(.2, 0, .65, 10)
    on_model = box(.25, 1, .55, 9)
    detached = box(1.2, 1, 1.5, 9)
    support = unary_union([on_model, detached])
    layers = [LayerSlice(i, 30.08 + i*.16,
                         model if i < 2 else EMPTY, EMPTY, EMPTY, EMPTY,
                         support if i >= 2 else EMPTY,
                         theta_window=(0, TURN))
              for i in range(6)]
    return layers, on_model, detached


class BottomSupportInterfaceTests(unittest.TestCase):
    def cfg(self, **kwargs):
        return replace(SlicerConfig(support_bottom_interface_layers=2,
                   support_bottom_interface_spacing=.7,
                   min_path_length=0, support_meander_max_angle_deg=360,
                   generate_gap_support=False), **kwargs)

    def test_default_zero_leaves_all_existing_support_untouched(self):
        layers, on_model, detached = stack()
        add_bottom_support_interfaces(layers, SlicerConfig())
        self.assertTrue(all(s.support_bottom_interface is None for s in layers))
        self.assertEqual(SlicerConfig().support_bottom_interface_layers,0)
        self.assertIsNone(SlicerConfig().support_bottom_interface_spacing)

    def test_model_top_only_seeds_matching_support_and_exact_depth(self):
        layers, on_model, detached = stack()
        add_bottom_support_interfaces(layers,self.cfg())
        for index in (2,3):
            interface=layers[index].support_bottom_interface
            self.assertFalse(interface.is_empty)
            self.assertAlmostEqual(interface.area,on_model.area,delta=1e-9)
            self.assertTrue(interface.covers(Point(.4,4)))
            self.assertTrue(interface.disjoint(detached))
        for index in (0,1,4,5):
            self.assertIsNone(layers[index].support_bottom_interface)
        self.assertTrue(all(l.support.equals(
            unary_union([on_model,detached])) for l in layers[2:]))

    def test_bottom_gap_uses_real_radial_faces_and_does_not_cut_other_support(self):
        layers, on_model, detached = stack()
        cfg = self.cfg(support_bottom_contact_gap_mm=.24)
        add_bottom_support_interfaces(layers,cfg)
        # Model upper deposited face = 30.32; support lower faces
        # 30.32, 30.48, 30.64, 30.80; first legal base is slab 4.
        for i in (2,3):
            self.assertTrue(layers[i].support.equals(detached))
            self.assertTrue(layers[i].support_bottom_interface is None)
        for i in (4,5):
            self.assertTrue(layers[i].support.covers(on_model))
            self.assertTrue(layers[i].support_bottom_interface.equals(on_model))
            inner = layers[i].radius - cfg.layer_height/2
            self.assertGreaterEqual(inner-(layers[1].radius+cfg.layer_height/2),.24-1e-9)
        for layer in layers:
            self.assertTrue(layer.occupied.equals(
                box(.2,0,.65,10) if layer.index<2 else EMPTY))

    def test_exact_physical_boundaries_and_independent_gap_without_interfaces(self):
        for gap, cut_count in ((0,0),(.16,1),(.160001,2),(.32,2),(.320001,3)):
            with self.subTest(gap=gap):
                layers,on_model,detached=stack()
                add_bottom_support_interfaces(layers,self.cfg(
                    support_bottom_contact_gap_mm=gap,
                    support_bottom_interface_layers=0))
                for i in range(2,6):
                    expect=detached if i<2+cut_count else unary_union([on_model,detached])
                    self.assertTrue(layers[i].support.equals(expect))
                    self.assertIsNone(layers[i].support_bottom_interface)
        # Too large a gap may eliminate the entire contacting support branch.
        layers,on_model,detached=stack()
        add_bottom_support_interfaces(layers,self.cfg(
            support_bottom_contact_gap_mm=1.0))
        self.assertTrue(all(layer.support.equals(detached) for layer in layers[2:]))
        self.assertTrue(all(layer.support_bottom_interface is None for layer in layers))

    def test_variable_slab_heights_do_not_round_up_to_nominal_count(self):
        heights=(.16,.12,.08,.2,.1,.1)
        edge=30.0
        layers,on_model,detached=stack()
        for index,height in enumerate(heights):
            layers[index].deposited_thickness=height
            layers[index].radius=edge+height/2
            edge+=height
        # Model ends at 30.28. Support begins at the faces:
        # 30.28, 30.36, 30.56, 30.66. For .12-mm gap only
        # slab 2 must be carved; slab 3 starts 0.08 mm above
        # and must also be carved. Slab 4 starts 0.28 mm above.
        add_bottom_support_interfaces(layers,self.cfg(support_bottom_contact_gap_mm=.12))
        self.assertTrue(layers[2].support.equals(detached))
        self.assertTrue(layers[3].support.equals(detached))
        self.assertTrue(layers[4].support_bottom_interface.equals(on_model))
        self.assertTrue(layers[5].support_bottom_interface.equals(on_model))

    def test_gap_carves_preexisting_top_interface_only_on_contacting_footprint(self):
        layers,on_model,detached=stack()
        for i in (2,3):
            layers[i].support_interface=unary_union([on_model,detached])
        add_bottom_support_interfaces(layers,self.cfg(support_bottom_contact_gap_mm=.2))
        for i in (2,3):
            self.assertTrue(layers[i].support_interface.equals(detached))
        self.assertTrue(layers[4].support_bottom_interface.equals(on_model))

    def test_no_jumping_gaps_and_shifted_periodic_contacts(self):
        layers, _, _ = stack()
        layers[2].support = box(1.2,1,1.5,9)
        add_bottom_support_interfaces(layers,self.cfg())
        self.assertTrue(all(l.support_bottom_interface is None for l in layers))
        shifted, on_model, _ = stack()
        shifted[1].occupied = box(.2 + TURN,0,.65+TURN,10)
        add_bottom_support_interfaces(shifted,self.cfg())
        self.assertTrue(shifted[2].support_bottom_interface.equals(on_model))
        self.assertTrue(shifted[3].support_bottom_interface.equals(on_model))

    def test_top_bottom_interfaces_do_not_double_deposit_or_bridge_wedges(self):
        layers, on_model, detached = stack()
        add_bottom_support_interfaces(layers,self.cfg())
        layers[2].support_interface = box(.3,2,.45,7)
        paths = build_layer_toolpaths(layers[2],self.cfg())
        roles=paths.paths_by_role
        self.assertTrue(roles["SUPPORT_INTERFACE"])
        self.assertTrue(roles["SUPPORT_BOTTOM_INTERFACE"])
        self.assertTrue(roles["SUPPORT"])
        top=layers[2].support_interface
        for p in roles["SUPPORT_BOTTOM_INTERFACE"]:
            for point in p.points:
                self.assertFalse(top.contains(Point(point)))
                self.assertTrue(on_model.buffer(1e-7).covers(Point(point)))
        self.assertEqual({p.role for p in roles["SUPPORT_BOTTOM_INTERFACE"]},
                         {"SUPPORT_BOTTOM_INTERFACE"})
        self.assertEqual(_role_print_settings("SUPPORT_BOTTOM_INTERFACE",self.cfg(
            support_interface_print_feed_mm_per_min=700,
            support_interface_flow_ratio=.8)), (700,.8))

    def test_exposed_top_surface_with_no_model_contact_does_not_seed(self):
        layers,on_model,detached=stack()
        layers[1].occupied = box(2.1,0,2.5,10)
        add_bottom_support_interfaces(layers,self.cfg())
        self.assertTrue(all(l.support_bottom_interface is None for l in layers))

    def test_cli_and_validation(self):
        for name in ("quality","legacy-cli","library"):
            parser=build_parser()
            args=parser.parse_args(["--profile",name])
            self.assertEqual(_profile_value(args,"support_bottom_interface_layers"),0)
            self.assertIsNone(_profile_value(args,"support_bottom_interface_spacing"))
            self.assertEqual(_profile_value(args,"support_bottom_contact_gap_mm"),0.0)
            args=parser.parse_args(["--profile",name,
                "--support-bottom-interface-layers","3",
                "--support-bottom-interface-spacing","0.64",
                "--support-bottom-contact-gap-mm","0.24"])
            self.assertEqual(_profile_value(args,"support_bottom_interface_layers"),3)
            self.assertEqual(_profile_value(args,"support_bottom_interface_spacing"),.64)
            self.assertEqual(_profile_value(args,"support_bottom_contact_gap_mm"),.24)
        for changes in ({"support_bottom_interface_layers":-1},
                        {"support_bottom_interface_layers":.5},
                        {"support_bottom_interface_spacing":0},
                        {"support_bottom_interface_spacing":math.inf},
                        {"support_bottom_interface_spacing":math.nan},
                        {"support_bottom_contact_gap_mm":-1},
                        {"support_bottom_contact_gap_mm":math.inf},
                        {"support_bottom_contact_gap_mm":math.nan}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                self.cfg(**changes).validate()


if __name__=="__main__":
    unittest.main()
