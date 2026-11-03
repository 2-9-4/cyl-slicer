from __future__ import annotations

from dataclasses import replace
import unittest

from cyl_slicer_v2.gcode import emit_gcode
from cyl_slicer_v2.nonplanar_ownership import stage_conformal_top_skin_replacement
from cyl_slicer_v2.tests.test_nonplanar_ownership import fixture


def count_identity(paths, needle):
    return sum(item is needle for group in paths for role_paths in group.values()
               for item in role_paths)


class StagedConformalOwnershipTests(unittest.TestCase):
    def test_atomic_replacement_keeps_one_role_path_and_one_island_alias(self):
        cfg,layers,generated,candidate=fixture()
        old_lp=generated[1]
        before_roles={k:list(v) for k,v in old_lp.paths_by_role.items()}
        before_island_aliases=[
            {role:tuple(id(path) for path in paths)
             for role,paths in group.items()}
            for group in old_lp.island_paths
        ]
        transfer=stage_conformal_top_skin_replacement(
            layers,generated,candidate,cfg,target_layer_index=1)
        report=transfer.ownership_report
        staged=transfer.layer_toolpaths
        self.assertTrue(report["ownership_transfer_staged"])
        self.assertTrue(report["production_printing_blocked"])
        self.assertFalse(report["nonplanar_extrusion_enabled"])
        self.assertEqual(report["transferred_path_count"],1)
        self.assertEqual(len(staged),3)
        self.assertIs(staged[0],generated[0])
        self.assertIsNot(staged[1],generated[1])
        self.assertIs(staged[2],generated[2])
        index=report["reserved_legacy_paths"][0]["role_path_index"]
        old=before_roles["TOP_SKIN"][index]
        self.assertIs(staged[1].paths_by_role["TOP_SKIN"][index],candidate)
        self.assertEqual(len(before_roles["TOP_SKIN"]),
                         len(staged[1].paths_by_role["TOP_SKIN"]))
        self.assertEqual(count_identity(staged[1].island_paths,candidate),1)
        self.assertEqual(count_identity(staged[1].island_paths,old),0)
        self.assertEqual(count_identity(old_lp.island_paths,old),1)
        self.assertFalse(any(item is old for paths in
                             staged[1].paths_by_role.values() for item in paths))
        self.assertIs(old_lp.paths_by_role["TOP_SKIN"][index],old)
        self.assertEqual(before_island_aliases,[
            {role:tuple(id(path) for path in paths)
             for role,paths in group.items()}
            for group in old_lp.island_paths])
        # Other roles remain byte-for-byte unchanged by identity.
        for role,paths in before_roles.items():
            if role!="TOP_SKIN":
                self.assertEqual(paths,staged[1].paths_by_role[role])
                self.assertTrue(all(a is b for a,b in
                         zip(paths,staged[1].paths_by_role[role])))

    def test_staged_path_cannot_ever_be_emitted_by_old_gcode_pipeline(self):
        cfg,layers,generated,candidate=fixture()
        stage=stage_conformal_top_skin_replacement(
            layers,generated,candidate,cfg,target_layer_index=1)
        with self.assertRaisesRegex(ValueError,"Conformal XYZ"):
            emit_gcode(list(stage.layer_toolpaths),cfg)
        self.assertFalse(stage.ownership_report["full_head_collision_certified"])

    def test_partial_and_doubled_volume_are_rejected_without_mutation(self):
        cfg,layers,generated,candidate=fixture()
        for changed in (
            replace(candidate,points=[
               candidate.points[0],
               (candidate.points[-1][0],
                (candidate.points[0][1]+candidate.points[-1][1])*.5)]),
            replace(candidate,segment_areas_mm2=[
                area*2 for area in candidate.segment_areas_mm2]),
        ):
            before=tuple(tuple(id(p) for p in paths)
                         for paths in generated[1].paths_by_role.values())
            with self.subTest(candidate=changed),self.assertRaises(ValueError):
                stage_conformal_top_skin_replacement(
                    layers,generated,changed,cfg,target_layer_index=1)
            self.assertEqual(before,tuple(tuple(id(p) for p in paths)
                         for paths in generated[1].paths_by_role.values()))

    def test_missing_island_alias_is_hard_error(self):
        cfg,layers,generated,candidate=fixture()
        old=generated[1].paths_by_role["TOP_SKIN"][
             len(generated[1].paths_by_role["TOP_SKIN"])//2]
        for group in generated[1].island_paths:
            for role,paths in group.items():
                group[role]=[path for path in paths if path is not old]
        with self.assertRaisesRegex(ValueError,"island alias must occur exactly once"):
            stage_conformal_top_skin_replacement(
                layers,generated,candidate,cfg,target_layer_index=1)

    def test_duplicate_island_alias_is_hard_error(self):
        cfg,layers,generated,candidate=fixture()
        old=generated[1].paths_by_role["TOP_SKIN"][
             len(generated[1].paths_by_role["TOP_SKIN"])//2]
        next(group for group in generated[1].island_paths
             if any(path is old for paths in group.values() for path in paths))[
                 "TOP_SKIN"].append(old)
        with self.assertRaisesRegex(ValueError,"island alias must occur exactly once"):
            stage_conformal_top_skin_replacement(
                layers,generated,candidate,cfg,target_layer_index=1)

    def test_candidate_already_registered_elsewhere_cannot_be_added_twice(self):
        cfg,layers,generated,candidate=fixture()
        generated[0].paths_by_role["BRIDGE"]=[candidate]
        with self.assertRaisesRegex(ValueError,"already registered"):
            stage_conformal_top_skin_replacement(
                layers,generated,candidate,cfg,target_layer_index=1)

    def test_explicitly_empty_island_groups_have_valid_paths_by_role_fallback(self):
        cfg,layers,generated,candidate=fixture()
        generated[1]=replace(generated[1],island_paths=[])
        stage=stage_conformal_top_skin_replacement(
            layers,generated,candidate,cfg,target_layer_index=1)
        self.assertEqual(stage.layer_toolpaths[1].island_paths,[])
        self.assertEqual(sum(item is candidate
            for paths in stage.layer_toolpaths[1].paths_by_role.values()
            for item in paths),1)
        with self.assertRaises(ValueError):
            emit_gcode(list(stage.layer_toolpaths),cfg)


if __name__=="__main__":
    unittest.main()
