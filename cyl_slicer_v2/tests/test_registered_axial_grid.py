from __future__ import annotations

from collections import Counter
from dataclasses import replace
import math
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.models import LayerSlice
from cyl_slicer_v2.toolpaths import (_registered_axial_grid, line_grid_lines,
                                      build_layer_toolpaths, build_toolpaths)


TAU = 2 * math.pi
EMPTY = GeometryCollection()


def thetas(geom, radius, spacing, config):
    return [round(line.coords[0][0], 10)
            for line in line_grid_lines(geom, radius, spacing,
                                        direction="axial", config=config)]


def slab(index, radius, model=None, support=None, interface=None):
    return LayerSlice(index, radius, model if model is not None else EMPTY,
                      EMPTY, EMPTY, EMPTY, support if support is not None else EMPTY,
                      support_interface=interface)


class RegisteredAxialGridTests(unittest.TestCase):
    def cfg(self, **kwargs):
        return replace(SlicerConfig(axial_grid_mode="fixed-outer",
                            axial_grid_outer_radius_mm=60,
                            generate_gap_support=False,
                            wall_count=1, min_path_length=0,
                            support_merge_gap=0,
                            cooling_min_layer_time_s=0,
                            cooling_min_island_print_time_s=0), **kwargs)

    def test_registered_full_wrap_has_exact_unique_lines_and_outer_pitch(self):
        geom = box(0, 0, TAU, 10)
        cfg = self.cfg()
        for spacing in (0.4, 1.5, 8):
            n = math.ceil(TAU*60/spacing - 1e-10)
            at_inner = thetas(geom, 30, spacing, cfg)
            at_outer = thetas(geom, 60, spacing, cfg)
            self.assertEqual(at_inner, at_outer)
            self.assertEqual(len(at_outer), n)  # Not n+1 doubled at 2pi.
            self.assertAlmostEqual(at_outer[0], 0)
            self.assertLessEqual(TAU*60/n, spacing + 1e-8)
            self.assertAlmostEqual(at_outer[-1], TAU*(n-1)/n, places=9)

    def test_no_island_rephasing_or_seam_doubling(self):
        cfg = self.cfg(axial_grid_phase_deg=17)
        for radius in (30, 35, 59):
            a = thetas(box(.02, 0, .4, 10), radius, 1.5, cfg)
            b = thetas(box(.02, 30, .4, 40), radius, 1.5, cfg)
            c = thetas(box(TAU+.02, 0, TAU+.4, 10), radius, 1.5, cfg)
            self.assertEqual(a, b)
            self.assertEqual(len(a), len(c))
            for theta, shifted_theta in zip(a, c):
                self.assertAlmostEqual(shifted_theta - theta, TAU, delta=2e-9)
        self.assertEqual(len(thetas(box(0,0,TAU,5),50,4,self.cfg())),
                         math.ceil(TAU*60/4-1e-10))

    def test_nested_tracks_are_registered_subsets_and_never_exceed_gap(self):
        cfg = self.cfg(axial_grid_mode="nested", axial_grid_inner_radius_mm=18)
        geom = box(0, 0, TAU, 10)
        previous = set()
        for r in (18, 22, 30, 39, 50, 60):
            live = set(thetas(geom, r, 2, cfg))
            if previous:
                self.assertTrue(previous.issubset(live))
            step, phase, count, stride = _registered_axial_grid(r, 2, cfg)
            self.assertTrue(all(round(i*step,10) in live
                                for i in range(count) if i%stride==0))
            self.assertLessEqual(stride*step*r, 2+1e-7)
            previous = live
        self.assertGreater(len(thetas(geom,60,2,cfg)),
                           len(thetas(geom,18,2,cfg)))

    def test_nested_periodicity_even_when_outer_count_is_odd(self):
        cfg = self.cfg(axial_grid_mode="nested", axial_grid_inner_radius_mm=18)
        # Initial N=189 is odd. A periodic stride-2 lattice must pad it;
        # otherwise its last-to-first gap is only half the normal pitch.
        self.assertEqual(math.ceil(TAU * 60 / 2 - 1e-10) % 2, 1)
        step, _, count, stride = _registered_axial_grid(18, 2, cfg)
        self.assertGreater(stride, 1)
        self.assertEqual(count % stride, 0)
        for lo in (.03, 1.23, 5.95):
            baseline = thetas(box(lo, 0, lo + .21, 3), 18, 2, cfg)
            for turns in (-3, -1, 1, 7):
                shifted = thetas(box(lo + turns*TAU, 0,
                                    lo + turns*TAU + .21, 3), 18, 2, cfg)
                self.assertEqual(len(baseline), len(shifted))
                for a, b in zip(baseline, shifted):
                    self.assertAlmostEqual(b - a, turns*TAU, delta=2e-9)

    def test_nested_seam_gap_matches_every_other_gap_at_all_active_strides(self):
        cfg = self.cfg(axial_grid_mode="nested", axial_grid_inner_radius_mm=7)
        geom = box(0, 0, TAU, 6)
        prior = set()
        for radius in (7, 10, 18, 30, 45, 60):
            step, _, count, stride = _registered_axial_grid(radius, 2, cfg)
            self.assertEqual(count % stride, 0)
            live = thetas(geom, radius, 2, cfg)
            self.assertEqual(len(live), count // stride)
            gaps = [(b - a) * radius for a, b in zip(live, live[1:])]
            gaps.append((live[0] + TAU - live[-1]) * radius)
            self.assertLessEqual(max(gaps), 2 + 1e-7)
            self.assertLess(max(gaps) - min(gaps), 1e-7)
            self.assertTrue(prior.issubset(set(live)))
            prior = set(live)

    def test_nested_inner_reference_is_global_including_parallel_workers(self):
        cfg = replace(self.cfg(), axial_grid_mode="nested",
                      axial_grid_inner_radius_mm=None,
                      axial_grid_outer_radius_mm=None,
                      infill_spacing=2, sparse_infill_pattern="both")
        geom = box(.02, 0, .9, 12)
        layers = [slab(i, r, model=geom) for i, r in enumerate((16, 24, 38, 59))]
        serial = build_toolpaths(layers, replace(cfg, toolpath_workers=1))
        parallel = build_toolpaths(layers, replace(cfg, toolpath_workers=2))
        for one, two in zip(serial, parallel):
            self.assertEqual([p.points for p in one.paths_by_role.get("INFILL_AXIAL", [])],
                             [p.points for p in two.paths_by_role.get("INFILL_AXIAL", [])])
        # A standalone nested call cannot silently assume an inner reference
        # larger than the geometry, or the seam may have extra-close tracks.
        with self.assertRaisesRegex(ValueError, "inner_radius_mm"):
            build_layer_toolpaths(layers[0], replace(cfg, axial_grid_outer_radius_mm=60))
        safe = replace(cfg, axial_grid_outer_radius_mm=60,
                       axial_grid_inner_radius_mm=16 - cfg.layer_height * .5)
        self.assertTrue(build_layer_toolpaths(layers[0], safe).paths_by_role["INFILL_AXIAL"])

    def test_roles_have_their_own_spacing_for_infill_and_support(self):
        model = box(.05,0,.5,20)
        support = box(.05,0,.5,20)
        cfg = self.cfg()
        cfg = replace(cfg, infill_spacing=3, sparse_infill_pattern="both",
                      support_spacing=2, support_cross_spacing=5,
                      support_interface_spacing=1.0)
        a = build_layer_toolpaths(slab(0,35,model=model), cfg)
        b = build_layer_toolpaths(slab(1,55,model=model), cfg)
        a_angles = {round(p.points[0][0], 8) for p in a.paths_by_role.get("INFILL_AXIAL", [])}
        b_angles = {round(p.points[0][0], 8) for p in b.paths_by_role.get("INFILL_AXIAL", [])}
        self.assertTrue(a_angles)
        self.assertEqual(a_angles,b_angles)
        dense_a = build_layer_toolpaths(slab(0,35,model=model), replace(cfg, solid_infill_pattern="axial"))
        self.assertTrue(dense_a.paths_by_role.get("INFILL_AXIAL"))
        c = build_layer_toolpaths(slab(0,35,support=support), cfg)
        d = build_layer_toolpaths(slab(1,55,support=support), cfg)
        c_angles = {round(p.points[0][0],8) for p in c.paths_by_role.get("SUPPORT", [])
                    if abs(p.points[-1][0] - p.points[0][0]) < 1e-7}
        d_angles = {round(p.points[0][0],8) for p in d.paths_by_role.get("SUPPORT", [])
                    if abs(p.points[-1][0] - p.points[0][0]) < 1e-7}
        self.assertEqual(c_angles,d_angles)
        self.assertTrue(c_angles)

    def test_job_reference_is_global_and_parallel_matches_serial(self):
        cfg = replace(self.cfg(), axial_grid_outer_radius_mm=None,
                      sparse_infill_pattern="both", infill_spacing=2,
                      support_spacing=3, toolpath_workers=1)
        geom = box(.2,0,.5,12)
        layers = [slab(i,r,model=geom) for i,r in enumerate((30,36,42,48))]
        serial = build_toolpaths(layers,cfg)
        parallel = build_toolpaths(layers,replace(cfg,toolpath_workers=2))
        for a,b in zip(serial,parallel):
            for role in ("EXTERNAL_PERIMETER", "INFILL_AXIAL", "INFILL_RADIAL"):
                self.assertEqual([p.points for p in a.paths_by_role.get(role,[])],
                                 [p.points for p in b.paths_by_role.get(role,[])])
        positions = [{round(p.points[0][0],8) for p in layer.paths_by_role["INFILL_AXIAL"]}
                     for layer in serial]
        self.assertEqual(positions[0],positions[-1])
        self.assertFalse(serial[0].paths_by_role["INFILL_AXIAL"] == [])
        # Per-layer mode never requires a reference radius.
        self.assertTrue(build_toolpaths(layers,replace(cfg,axial_grid_mode="per-layer")))
        with self.assertRaisesRegex(ValueError,"outer radius"):
            build_layer_toolpaths(layers[0],cfg)

    def test_fixed_outer_refuses_excessive_small_radius_density(self):
        cfg=self.cfg(axial_grid_max_density_ratio=4)
        geom=box(0,0,2*math.pi,10)
        with self.assertRaisesRegex(ValueError,"Fixed-outer axial grid would be"):
            thetas(geom,8,1.5,cfg)
        # The nested grid keeps parent's exact registration but backs off
        # the number of live tracks when radius falls by a large factor.
        nested=self.cfg(axial_grid_mode="nested", axial_grid_inner_radius_mm=8)
        self.assertTrue(thetas(geom,8,1.5,nested))
        self.assertGreater(len(thetas(geom,8,1.5,nested)),1)
        # An explicitly calibrated density budget can permit fixed-outer,
        # but the caller must choose it rather than silently overfilling.
        override=self.cfg(axial_grid_max_density_ratio=10)
        actual=thetas(geom,8,1.5,override)
        self.assertGreater(len(actual),len(thetas(geom,8,1.5,nested)))
        self.assertEqual(thetas(geom,30,1.5,cfg),thetas(geom,60,1.5,cfg))
        # Per-layer mode remains unchanged regardless of a fixed-outer guard.
        unregistered=self.cfg(axial_grid_mode="per-layer")
        self.assertTrue(thetas(geom,8,1.5,unregistered))

    def test_dense_role_uses_its_own_density_budget(self):
        cfg=self.cfg()
        # Same radius: a dense 0.4-mm role can overfill while a sparse
        # 2-mm role is also dense by the same *ratio*: the guard is
        # dimensionless and applies to both, not just infill.
        with self.assertRaisesRegex(ValueError,"Fixed-outer axial grid"):
            _registered_axial_grid(10,0.4,cfg)
        with self.assertRaisesRegex(ValueError,"Fixed-outer axial grid"):
            _registered_axial_grid(10,2.0,cfg)
        self.assertTrue(_registered_axial_grid(20,2.0,cfg))

    def test_validation_and_profiles(self):
        for profile in ("quality","legacy-cli","library"):
            args = build_parser().parse_args(["--profile",profile])
            self.assertEqual(_profile_value(args,"axial_grid_mode"),"per-layer")
            self.assertIsNone(_profile_value(args,"axial_grid_outer_radius_mm"))
            self.assertEqual(_profile_value(args,"axial_grid_max_density_ratio"),4)
            self.assertIsNone(_profile_value(args,"axial_grid_inner_radius_mm"))
            selected = build_parser().parse_args(["--profile",profile,
                        "--axial-grid-mode","nested",
                        "--axial-grid-outer-radius-mm","60",
                        "--axial-grid-inner-radius-mm","18",
                        "--axial-grid-phase-deg","25",
                        "--axial-grid-max-density-ratio","8"])
            for key,value in (("axial_grid_mode","nested"),
                              ("axial_grid_outer_radius_mm",60),
                              ("axial_grid_inner_radius_mm",18),
                              ("axial_grid_phase_deg",25),
                              ("axial_grid_max_density_ratio",8)):
                self.assertEqual(_profile_value(selected,key),value)
        for arg in ({"axial_grid_mode":"magic"},
                    {"axial_grid_outer_radius_mm":0},
                    {"axial_grid_outer_radius_mm":float("inf")},
                    {"axial_grid_inner_radius_mm":0},
                    {"axial_grid_inner_radius_mm":float("nan")},
                    {"axial_grid_inner_radius_mm":61},
                    {"axial_grid_phase_deg":360},
                    {"axial_grid_phase_deg":float("nan")},
                    {"axial_grid_max_density_ratio":0},
                    {"axial_grid_max_density_ratio":float("inf")},
                    {"axial_grid_max_density_ratio":float("nan")}):
            with self.subTest(arg=arg),self.assertRaises(ValueError):
                replace(self.cfg(),**arg).validate()


if __name__=="__main__":
    unittest.main()
