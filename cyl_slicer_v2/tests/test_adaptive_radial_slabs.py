from __future__ import annotations

import math
import unittest

import numpy as np
import trimesh
from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.cli import build_parser, _profile_value
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.layers import layer_slabs
from cyl_slicer_v2.models import LayerSlice, LayerToolpaths, Toolpath
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, gap_support_slices


def slope_mesh():
    # One radial-plane normal active ONLY near R=30.08; one near
    # R=30.32 has a purely axial normal and must not cause refinement.
    vertices=np.array([
        [30.08,0,0], [30.08,.01,0], [30.08,0,1],
        [30.36,0,0], [30.42,0,0], [30.36,.1,0],
    ])
    return trimesh.Trimesh(vertices=vertices,
                           faces=np.array([[0,1,2],[3,4,5]]),
                           process=False)


def cfg(**kwargs):
    return SlicerConfig(bed_radius=30,min_radius=30,max_radius=30.48,
                        radial_layer_mode="adaptive",
                        adaptive_cusp_error_mm=.08,
                        adaptive_min_layer_height_mm=.04,
                        **kwargs)


class AdaptiveRadialSlabTests(unittest.TestCase):
    def test_only_active_facet_radial_normal_refines_layers(self):
        mesh=slope_mesh()
        c=cfg()
        adaptive=layer_slabs(mesh,c)
        self.assertEqual([round(h,6) for _,h in adaptive],
                         [.08,.08,.16,.16])
        self.assertEqual(len(layer_slabs(mesh,
            SlicerConfig(bed_radius=30,min_radius=30,max_radius=30.48,
                         radial_layer_mode="complete-terminal"))),3)
        previous=30.0
        for radius,height in adaptive:
            self.assertAlmostEqual(radius-height*.5,previous,places=9)
            previous=radius+height*.5
            self.assertGreaterEqual(height,c.adaptive_min_layer_height_mm-1e-9)
            self.assertLessEqual(height,c.layer_height+1e-9)
        self.assertAlmostEqual(previous,30.48,places=9)

    def test_conservative_cusp_limit_and_winding_flip(self):
        mesh=slope_mesh()
        for scale in (.08,.06):
            from dataclasses import replace
            c=replace(cfg(),adaptive_cusp_error_mm=scale)
            slabs=layer_slabs(mesh,c)
            for radius,height in slabs:
                if radius-height*.5 < 30.16:
                    self.assertLessEqual(height,scale+1e-9)
            inverted=mesh.copy()
            inverted.invert()
            self.assertEqual(slabs,layer_slabs(inverted,c))

    def test_minimum_thickness_fails_instead_of_violating_requested_cusp(self):
        mesh=slope_mesh()
        from dataclasses import replace
        c=replace(cfg(),adaptive_cusp_error_mm=.04)
        with self.assertRaisesRegex(ValueError,"Adaptive cusp target"):
            layer_slabs(mesh,c)

    def test_exact_terminal_face_rebalance_is_preserved_after_refinement(self):
        mesh=slope_mesh()
        from dataclasses import replace
        c=replace(cfg(),max_radius=30.35)
        slabs=layer_slabs(mesh,c)
        self.assertAlmostEqual(sum(h for _,h in slabs),.35,places=9)
        edge=30
        for radius,h in slabs:
            self.assertAlmostEqual(radius-h/2,edge,places=9)
            edge=radius+h/2
        self.assertAlmostEqual(edge,30.35,places=9)

    def test_explicit_first_layer_remains_whole_and_calibrated(self):
        from dataclasses import replace
        c=replace(cfg(),first_layer_height=.20)
        slabs=layer_slabs(slope_mesh(),c)
        self.assertAlmostEqual(slabs[0][0],30.1)
        self.assertAlmostEqual(slabs[0][1],.20)
        self.assertAlmostEqual(sum(h for _,h in slabs),.48,places=9)

    def test_adaptive_thickness_propagates_to_toolpath_and_gcode_E(self):
        c=cfg(wall_count=1,generate_gap_support=False)
        radius,height=layer_slabs(slope_mesh(),c)[0]
        empty=GeometryCollection()
        region=box(.1,0,.6,10)
        layer=LayerSlice(0,radius,region,region,empty,empty,empty,
                         deposited_thickness=height)
        paths=build_layer_toolpaths(layer,c)
        self.assertAlmostEqual(paths.nozzle_radius,radius+height/2,places=9)
        self.assertAlmostEqual(paths.deposited_thickness,height,places=9)
        single=LayerToolpaths(0,radius,{"SOLID_INFILL":[
            Toolpath("SOLID_INFILL",radius,[(.2,0),(.2,10)])
        ]},nozzle_radius=radius+height/2,deposited_thickness=height)
        code=emit_gcode([single],c)
        self.assertIn(";DEPOSITED_THICKNESS:0.0800",code)
        moves=[line for line in code.splitlines() if line.startswith("G1 X") and " E" in line]
        self.assertTrue(moves)
        e=float(moves[-1].split(" E",1)[1].split()[0])
        self.assertAlmostEqual(e,10*c.for_deposited_thickness(height).extrusion_per_mm,delta=1e-4)

    def test_support_foundation_uses_exact_faces_in_adaptive_mode(self):
        c=cfg(support_gap_layers=0,raft_layers=0)
        empty=GeometryCollection()
        first=LayerSlice(0,30.3975,box(.1,0,.5,10),empty,empty,empty,empty,
                         deposited_thickness=.095)
        support=gap_support_slices([first],c)
        self.assertTrue(support)
        self.assertAlmostEqual(support[-1].radius+support[-1].deposited_thickness/2,
                               first.radius-first.deposited_thickness/2,places=9)

    def test_cli_profiles_and_bad_parameters(self):
        from dataclasses import replace
        for profile in ("quality","legacy-cli","library"):
            parser=build_parser()
            args=parser.parse_args(["--profile",profile])
            self.assertEqual(_profile_value(args,"radial_layer_mode"),"fixed")
            self.assertEqual(_profile_value(args,"adaptive_cusp_error_mm"),.08)
            self.assertEqual(_profile_value(args,"adaptive_min_layer_height_mm"),.04)
            explicit=parser.parse_args(["--profile",profile,"--radial-layer-mode","adaptive",
                "--adaptive-cusp-error-mm",".06","--adaptive-min-layer-height-mm",".03"])
            self.assertEqual(_profile_value(explicit,"radial_layer_mode"),"adaptive")
            self.assertEqual(_profile_value(explicit,"adaptive_cusp_error_mm"),.06)
            self.assertEqual(_profile_value(explicit,"adaptive_min_layer_height_mm"),.03)
        for options in ({"adaptive_cusp_error_mm":0},
                        {"adaptive_cusp_error_mm":math.nan},
                        {"adaptive_min_layer_height_mm":-1},
                        {"adaptive_min_layer_height_mm":math.inf},
                        {"adaptive_min_layer_height_mm":.12}):
            with self.subTest(options=options),self.assertRaises(ValueError):
                replace(cfg(),**options).validate()


if __name__=="__main__":
    unittest.main()
