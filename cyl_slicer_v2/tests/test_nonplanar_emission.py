from __future__ import annotations

from dataclasses import replace
import math
import re
import unittest

from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import emit_gcode, emit_toolpath
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath
from cyl_slicer_v2.nonplanar_emission import preview_conformal_toolpath_xyz
from cyl_slicer_v2.nonplanar_motion import evaluate_conformal_move


def path(*, points=None,radii=None,areas=None,role="TOP_SKIN"):
    if points is None:
        points=[(.1,2),(.8,5)]
    if radii is None:
        radii=[30.2,31.8]
    if areas is None:
        areas=[.062]
    return Toolpath(role,30.2,points,point_radii=radii,
                    segment_areas_mm2=areas)


def cfg(**kw):
    return SlicerConfig(raft_layers=0,generate_gap_support=False,
                        cooling_min_layer_time_s=0,
                        cooling_min_island_print_time_s=0,**kw)


def moves(preview):
    rows=[]
    for line in preview.commands:
        if line.startswith("G1 "):
            fields={m.group(1):float(m.group(2))
                    for m in re.finditer(r"\b([XYZEF])(-?\d+(?:\.\d+)?)",line)}
            rows.append(fields)
    return rows


class ConformalCandidateEmissionTests(unittest.TestCase):
    def test_opt_in_model_fields_preserve_legacy_path_positional_api(self):
        legacy=Toolpath("TOP_SKIN",30.2,[(0,2),(.2,2)],False)
        self.assertIsNone(legacy.point_radii)
        self.assertIsNone(legacy.segment_areas_mm2)
        self.assertIsNone(legacy.segment_widths)
        candidate=path()
        self.assertEqual(candidate.points,[(.1,2),(.8,5)])
        self.assertEqual(candidate.point_radii,[30.2,31.8])

    def test_true_xyz_degrees_bed_z_and_physical_e_are_emitted_in_preview(self):
        candidate=path()
        report=preview_conformal_toolpath_xyz(candidate,cfg(),
                                               max_local_e_error_fraction=.002)
        self.assertEqual(report.status,"uncertified_motion_preview_only")
        self.assertFalse(report.nonplanar_extrusion_enabled)
        self.assertGreater(report.segment_count,1)
        self.assertEqual(report.source_segment_count,1)
        self.assertTrue(report.commands[0].startswith(";"))
        self.assertFalse(any(line.startswith(("G28","M104","M109","M190"))
                             for line in report.commands))
        lines=moves(report)
        self.assertEqual(len(lines),report.segment_count)
        self.assertEqual(set(lines[0]),{"X","Y","Z","E","F"})
        self.assertAlmostEqual(lines[-1]["X"],5,places=3)
        self.assertAlmostEqual(lines[-1]["Y"],math.degrees(.7),places=3)
        self.assertAlmostEqual(lines[-1]["Z"],1.8,places=4)
        self.assertTrue(all(line["F"]>0 for line in lines))
        self.assertTrue(all(b["E"]>a["E"]
                            for a,b in zip(lines,lines[1:])))
        exact=evaluate_conformal_move(
            (.1,2,30.2),(.8,5,31.8),
            physical_print_feed_mm_per_min=cfg().feed_mm_per_min,
            bead_area_mm2=.062,filament_diameter_mm=cfg().filament_diameter)
        self.assertAlmostEqual(report.exact_physical_length_mm,
                               exact.length_mm,delta=1e-9)
        self.assertAlmostEqual(report.emitted_filament_length_mm,
                               exact.filament_mm,delta=.00001)
        self.assertAlmostEqual(lines[-1]["E"],
                               report.emitted_filament_length_mm,places=5)
        self.assertLessEqual(report.maximum_local_e_error_fraction,.0020000001)

    def test_per_segment_areas_and_existing_calibration_multipliers(self):
        config=cfg(flow_percent=115)
        calibrated=replace(config,extrusion_per_mm=config.extrusion_per_mm*1.25)
        candidate=path(points=[(.1,1),(.25,3),(.4,4)],
                       radii=[30.2,30.5,30.7],
                       areas=[.055,.088])
        result=preview_conformal_toolpath_xyz(candidate,calibrated)
        values=[evaluate_conformal_move(
                 (p[0],p[1],r0),(q[0],q[1],r1),
                 physical_print_feed_mm_per_min=calibrated.feed_mm_per_min,
                 bead_area_mm2=area,filament_diameter_mm=calibrated.filament_diameter)
                for p,q,r0,r1,area in zip(
                   candidate.points,candidate.points[1:],
                   candidate.point_radii,candidate.point_radii[1:],
                   candidate.segment_areas_mm2)]
        expected=sum(v.filament_mm for v in values)*1.25*1.15
        self.assertAlmostEqual(result.ideal_filament_length_mm,expected,
                               delta=1e-9)
        self.assertAlmostEqual(result.emitted_filament_length_mm,
                               expected,delta=.00001)
        self.assertEqual(candidate.segment_areas_mm2,[.055,.088])

    def test_rounding_flow_cap_uses_true_native_xyz_not_only_xy(self):
        candidate=path(points=[(0,2),(.06,2.4)],
                       radii=[30.2,31.2],areas=[.07])
        config=cfg(max_volumetric_flow_mm3_s=.05)
        report=preview_conformal_toolpath_xyz(candidate,config,
                                               max_local_e_error_fraction=.001)
        rows=moves(report)
        filament_area=math.pi*(config.filament_diameter/2)**2
        prev=(2.0,0.0,.2,0.0)
        for row in rows:
            now=(row["X"],row["Y"],row["Z"],row["E"])
            native=math.dist(prev[:3],now[:3])
            delta_e=now[3]-prev[3]
            rate=delta_e*filament_area*row["F"]/(60*native)
            self.assertLessEqual(rate,.05+1e-9)
            self.assertGreater(now[2],0)
            prev=now

    def test_global_periodic_theta_offset_produces_same_local_machine_commands(self):
        original=path()
        base=preview_conformal_toolpath_xyz(original,cfg())
        for turn in (-3,-1,1,8):
            copy=path(points=[(t+turn*2*math.pi,x)
                              for t,x in original.points])
            shifted=preview_conformal_toolpath_xyz(copy,cfg())
            self.assertEqual(moves(base),moves(shifted))
            self.assertAlmostEqual(base.exact_physical_length_mm,
                                   shifted.exact_physical_length_mm,delta=1e-10)

    def test_one_point_is_insufficient_and_unsafe_shapes_fail(self):
        cases=[
            path(points=[(.1,2)],radii=[30.2],areas=[]),
            path(points=[(.1,2),(.2,3)],radii=[30.2],areas=[.062]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,31.2],areas=[]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,float("inf")],areas=[.06]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,29.9],areas=[.06]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,31.2],areas=[0]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,31.2],areas=[float("nan")]),
            path(points=[(.1,2),(.1,2)],radii=[30.2,30.2],areas=[.06]),
            path(points=[(.1,2),(.2,3)],radii=[30.2,31.2],areas=[.06],role="SOLID_INFILL"),
        ]
        for candidate in cases:
            with self.subTest(candidate=candidate),self.assertRaises(ValueError):
                preview_conformal_toolpath_xyz(candidate,cfg())
        with self.assertRaisesRegex(ValueError,"closed"):
            preview_conformal_toolpath_xyz(replace(path(),closed=True),cfg())
        with self.assertRaisesRegex(ValueError,"segment_widths"):
            preview_conformal_toolpath_xyz(replace(path(),segment_widths=[.3]),cfg())
        with self.assertRaises(ValueError):
            preview_conformal_toolpath_xyz(path(),cfg(),
                                           max_segments=1,
                                           max_local_e_error_fraction=.00001)
        with self.assertRaisesRegex(ValueError,"native-Y"):
            preview_conformal_toolpath_xyz(
                path(points=[(0,2),(25*math.pi,5)],
                     radii=[30.2,31.8],areas=[.06]),cfg())

    def test_production_emitter_rejects_experimental_xyz_before_any_gcode(self):
        candidate=path()
        with self.assertRaisesRegex(ValueError,"Conformal XYZ candidate"):
            emit_gcode([LayerToolpaths(0,30.2,{"TOP_SKIN":[candidate]})],cfg())
        with self.assertRaisesRegex(ValueError,"Conformal XYZ candidates"):
            emit_toolpath(candidate,GcodeState(),cfg(),0)
        # A conventional path is not blocked by the optional fields.
        ordinary=Toolpath("TOP_SKIN",30.2,[(.1,2),(.2,5)])
        code=emit_gcode([LayerToolpaths(0,30.2,{"TOP_SKIN":[ordinary]})],cfg())
        self.assertIn(";TYPE:TOP_SKIN",code)
        self.assertNotIn("CONFORMAL_XYZ",code)

    def test_sub_resolution_motion_cannot_become_e_only_at_constant_xyz(self):
        with self.assertRaises(ValueError):
            preview_conformal_toolpath_xyz(
                path(points=[(0,2),(1e-9,2)],radii=[30.2,30.2],
                     areas=[.07]),cfg())


if __name__ == "__main__":
    unittest.main()
