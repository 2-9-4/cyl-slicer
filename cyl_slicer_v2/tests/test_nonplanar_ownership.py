from __future__ import annotations

from dataclasses import replace
import math
import unittest

from shapely.geometry import GeometryCollection, box

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.models import LayerSlice, Toolpath
from cyl_slicer_v2.nonplanar_ownership import audit_conformal_top_skin_ownership
from cyl_slicer_v2.toolpaths import build_layer_toolpaths


EMPTY=GeometryCollection()


def fixture():
    cfg=SlicerConfig(bed_radius=30,layer_height=.16,nozzle_width=.42,
                     wall_count=1,monotonic_top_skin=True,
                     solid_infill_pattern="axial",gap_fill_enabled=False,
                     generate_gap_support=False,raft_layers=0)
    region=box(.2,1,1.2,10)
    previous=LayerSlice(0,30.08,region,EMPTY,EMPTY,EMPTY,EMPTY,
                         theta_window=(.2,1.2))
    target=LayerSlice(1,30.24,region,region,EMPTY,EMPTY,EMPTY,
                      theta_window=(.2,1.2))
    following=LayerSlice(2,30.40,EMPTY,EMPTY,EMPTY,EMPTY,EMPTY,
                         theta_window=(.2,1.2))
    layers=[previous,target,following]
    toolpaths=[build_layer_toolpaths(layer,cfg) for layer in layers]
    top=toolpaths[1].paths_by_role.get("TOP_SKIN",[])
    if not top:
        raise AssertionError("test fixture did not generate true monotonic top skin")
    old=top[len(top)//2]
    candidate=Toolpath("TOP_SKIN",target.radius,list(old.points),
        point_radii=[target.radius+.08+.002*(x-5)
                     for _,x in old.points],
        segment_areas_mm2=[cfg.bead_area]*(len(old.points)-1))
    return cfg,layers,toolpaths,candidate


def audit(cfg,layers,generated,candidate,**kwargs):
    return audit_conformal_top_skin_ownership(
        layers,generated,candidate,cfg,target_layer_index=1,**kwargs)


class ConformalTopSkinOwnershipTests(unittest.TestCase):
    def test_whole_legacy_skin_may_be_reserved_without_modifying_either_alias(self):
        cfg,layers,generated,candidate=fixture()
        paths_before=list(generated[1].paths_by_role["TOP_SKIN"])
        islands_before=[(role,tuple(id(p) for p in paths))
                        for group in generated[1].island_paths
                        for role,paths in group.items()]
        report=audit(cfg,layers,generated,candidate)
        self.assertEqual(report["status"],
            "exclusive_legacy_TOP_SKIN_replacement_candidate",report)
        self.assertFalse(report["nonplanar_extrusion_enabled"])
        self.assertFalse(report["full_head_collision_certified"])
        self.assertFalse(report["conventional_toolpaths_mutated"])
        self.assertEqual(len(report["reserved_legacy_paths"]),1)
        reserved=report["reserved_legacy_paths"][0]
        self.assertEqual(reserved["layer_index"],1)
        self.assertEqual(reserved["role"],"TOP_SKIN")
        self.assertIs(generated[1].paths_by_role["TOP_SKIN"][
            reserved["role_path_index"]],paths_before[reserved["role_path_index"]])
        self.assertEqual(list(generated[1].paths_by_role["TOP_SKIN"]),paths_before)
        self.assertEqual(islands_before,[(role,tuple(id(p) for p in paths))
                           for group in generated[1].island_paths
                           for role,paths in group.items()])
        self.assertEqual(report["backing_layers_checked"],1)
        self.assertGreater(report["original_volume_mm3"],0)
        self.assertLess(report["relative_volume_difference"],.05)

    def test_partial_original_skin_cannot_be_removed_as_whole_path(self):
        cfg,layers,generated,candidate=fixture()
        shortened=list(candidate.points)
        a,b=shortened[0],shortened[-1]
        shortened[-1]=(b[0],.5*(a[1]+b[1]))
        changed=replace(candidate,points=shortened)
        result=audit(cfg,layers,generated,changed)
        self.assertEqual(result["status"],"rejected")
        self.assertEqual(result["reason"],"partial_legacy_TOP_SKIN_overlap")
        self.assertEqual(result["reserved_legacy_paths"],[])

    def test_extra_conformal_volume_cannot_be_reserved(self):
        cfg,layers,generated,candidate=fixture()
        changed=replace(candidate,
                    segment_areas_mm2=[a*2 for a in candidate.segment_areas_mm2])
        result=audit(cfg,layers,generated,changed)
        self.assertEqual(result["reason"],
                         "conformal_volume_does_not_replace_legacy_volume_within_budget")
        self.assertGreater(result["conformal_volume_mm3"],
                           result["original_volume_mm3"]*1.8)

    def test_exposed_top_must_not_be_covered_by_a_higher_slab(self):
        cfg,layers,generated,candidate=fixture()
        layers[2]=replace(layers[2],occupied=layers[1].occupied)
        result=audit(cfg,layers,generated,candidate)
        self.assertEqual(result["reason"],
                         "candidate_bead_corridor_is_not_fully_exposed_top_skin")

    def test_buried_layer_must_have_solid_backing(self):
        cfg,layers,generated,candidate=fixture()
        layers[0]=replace(layers[0],occupied=EMPTY)
        result=audit(cfg,layers,generated,candidate)
        self.assertEqual(result["reason"],
                         "required_solid_backing_layer_is_missing")

    def test_solid_backing_must_cover_entire_bead_corridor(self):
        cfg,layers,generated,candidate=fixture()
        layers[0]=replace(layers[0],
                  occupied=box(.2,1,.7,10))
        result=audit(cfg,layers,generated,candidate)
        self.assertEqual(result["reason"],"candidate_lacks_full_solid_backing")

    def test_protected_bridge_or_wall_path_must_never_be_removed(self):
        cfg,layers,generated,candidate=fixture()
        bridge=Toolpath("BRIDGE",layers[1].radius,list(candidate.points))
        generated[1].paths_by_role["BRIDGE"]=[bridge]
        result=audit(cfg,layers,generated,candidate)
        self.assertEqual(result["reason"],
                         "conformal_candidate_conflicts_with_protected_BRIDGE")
        self.assertEqual(generated[1].paths_by_role["BRIDGE"],[bridge])

    def test_one_full_theta_shift_is_the_same_physical_ribbon(self):
        cfg,layers,generated,candidate=fixture()
        shifted=replace(candidate,points=[
            (theta+2*math.pi,axial) for theta,axial in candidate.points])
        result=audit(cfg,layers,generated,shifted)
        self.assertEqual(result["status"],
                         "exclusive_legacy_TOP_SKIN_replacement_candidate",result)

    def test_extrusion_is_not_allowed_without_distinct_top_skin_role(self):
        cfg,layers,generated,candidate=fixture()
        result=audit(replace(cfg,monotonic_top_skin=False),
                     layers,generated,candidate)
        self.assertEqual(result["reason"],
                         "requires_separate_monotonic_TOP_SKIN_role")
        self.assertFalse(result["nonplanar_extrusion_enabled"])

    def test_no_unknown_outer_skin_from_truncated_layer_stack(self):
        cfg,layers,generated,candidate=fixture()
        result=audit(cfg,layers[:2],generated[:2],candidate)
        self.assertEqual(result["reason"],
                         "outer_layer_occupancy_is_unknown_for_truncated_job")

    def test_radii_must_lie_near_one_deposited_outer_face(self):
        cfg,layers,generated,candidate=fixture()
        changed=replace(candidate,
                point_radii=[r+.5 for r in candidate.point_radii])
        result=audit(cfg,layers,generated,changed)
        self.assertEqual(result["reason"],
                         "candidate_radius_outside_single_deposited_top_slab")

    def test_bottom_solid_footprint_cannot_be_reassigned(self):
        cfg,layers,generated,candidate=fixture()
        layers[1]=replace(layers[1],bottom_surfaces=layers[1].occupied)
        result=audit(cfg,layers,generated,candidate)
        self.assertEqual(result["reason"],
                         "candidate_overlaps_undifferentiated_bottom_solid")

    def test_invalid_settings_fail_closed(self):
        cfg,layers,generated,candidate=fixture()
        for opts in (
            {"match_tolerance_mm":-1},
            {"protected_role_clearance_mm":float("nan")},
            {"radial_face_tolerance_mm":-1},
            {"max_volume_difference_fraction":1},
            {"minimum_solid_backing_layers":0},
            {"target_layer_index":1.5},
        ):
            with self.subTest(opts=opts),self.assertRaises(ValueError):
                params={"target_layer_index":1}
                params.update(opts)
                audit_conformal_top_skin_ownership(
                    layers,generated,candidate,cfg,**params)


if __name__=="__main__":
    unittest.main()
