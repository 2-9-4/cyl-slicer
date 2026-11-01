from __future__ import annotations

import math
import unittest

from cyl_slicer_v2.tests.legacy_config import LegacySlicerConfig as SlicerConfig
from cyl_slicer_v2.coordinates import cylindrical_distance
from cyl_slicer_v2.gcode import (
    _orient_path_for_current_position,
    _ordered_layer_paths,
    _seam_angle_excluded,
    emit_gcode,
    parse_extrusion_segments,
)
from cyl_slicer_v2.models import GcodeState, LayerToolpaths, Toolpath


def rect(theta_min=0.15, theta_max=1.5, axial_min=0, axial_max=10, radius=30):
    return Toolpath("EXTERNAL_PERIMETER", radius, [
        (theta_min, axial_min), (theta_max, axial_min),
        (theta_max, axial_max), (theta_min, axial_max),
        (theta_min, axial_min),
    ], closed=True)


def length(path):
    return sum(cylindrical_distance(a, b, path.radius)
               for a, b in zip(path.points, path.points[1:]))


def oriented_area(path):
    return sum(a[0] * b[1] - b[0] * a[1]
               for a, b in zip(path.points, path.points[1:]))


class SeamPolicyTests(unittest.TestCase):
    def test_default_nearest_still_projects_onto_edge(self):
        original = rect(.1, .4)
        state = GcodeState(theta_accum=.41, axial=5)
        old = _orient_path_for_current_position(original, state)
        new = _orient_path_for_current_position(
            original, state, SlicerConfig(seam_policy="nearest")
        )
        self.assertEqual(old.points, new.points)
        self.assertAlmostEqual(new.points[0][1], 5)

    def test_aligned_seam_projects_to_target_angle_and_preserves_winding(self):
        # The desired 25-degree ray crosses two horizontal polygon edges.
        original = rect(theta_min=.1, theta_max=.9)
        c = SlicerConfig(seam_policy="aligned", seam_angle_deg=25)
        for current in (-3.0, .4, 4.0, 10.0):
            with self.subTest(current=current):
                p = _orient_path_for_current_position(original,
                    GcodeState(theta_accum=current, axial=5), c)
                self.assertAlmostEqual(math.degrees(p.points[0][0]), 25, places=6)
                self.assertAlmostEqual(p.points[0][0], p.points[-1][0])
                self.assertAlmostEqual(p.points[0][1], p.points[-1][1])
                self.assertAlmostEqual(length(p), length(original), places=9)
                self.assertAlmostEqual(oriented_area(p), oriented_area(original), places=9)
        self.assertEqual(original.points[0], (.1, 0))

    def test_corner_prefers_strong_turn_at_target_phase(self):
        original = rect(theta_min=.2, theta_max=.8)
        # Corner seams should avoid splitting a straight edge.
        p = _orient_path_for_current_position(original,
            GcodeState(theta_accum=.5, axial=5),
            SlicerConfig(seam_policy="corner", seam_angle_deg=math.degrees(.8)))
        self.assertAlmostEqual(p.points[0][0], .8)
        self.assertIn(p.points[0], original.points[:-1])
        self.assertAlmostEqual(length(p), length(original), places=9)

    def test_excluded_phase_including_wraparound(self):
        original = rect(theta_min=math.radians(350), theta_max=math.radians(370))
        c = SlicerConfig(seam_policy="aligned", seam_angle_deg=0,
                         seam_exclude_deg=(355, 5))
        p = _orient_path_for_current_position(original,
            GcodeState(theta_accum=math.radians(0), axial=3), c)
        self.assertFalse(_seam_angle_excluded(p.points[0][0], (355, 5)))
        self.assertAlmostEqual(length(p), length(original), places=8)
        self.assertAlmostEqual(oriented_area(p), oriented_area(original), places=8)
        self.assertAlmostEqual(math.degrees(p.points[0][0]) % 360, 5, places=5)

    def test_nearest_with_exclusion_avoids_forbidden_sector(self):
        original = rect(theta_min=.1, theta_max=.9)
        forbidden = (math.degrees(.35), math.degrees(.55))
        p = _orient_path_for_current_position(original,
            GcodeState(theta_accum=.45, axial=0),
            SlicerConfig(seam_exclude_deg=forbidden))
        self.assertFalse(_seam_angle_excluded(p.points[0][0], forbidden))
        self.assertAlmostEqual(length(p), length(original), places=9)

    def test_fully_excluded_closed_loop_aborts_before_print(self):
        # A loop located only in a tiny angular sector.
        original = rect(theta_min=.1, theta_max=.2)
        with self.assertRaisesRegex(ValueError, "No eligible perimeter seam"):
            _orient_path_for_current_position(original, GcodeState(),
                SlicerConfig(seam_exclude_deg=(0, 20)))

    def test_open_perimeters_are_not_rotated(self):
        path = Toolpath("EXTERNAL_PERIMETER", 30,
                        [(math.radians(350), 0), (math.radians(370), 0)])
        c = SlicerConfig(seam_policy="aligned", seam_angle_deg=180,
                         seam_exclude_deg=(0, 90))
        self.assertIs(_orient_path_for_current_position(path, GcodeState(), c), path)

    def test_gcode_length_identical_and_starts_fixed_across_layers(self):
        original = rect(theta_min=.1, theta_max=.9)
        layers = [LayerToolpaths(i, 30.5+i*.2, {
            "EXTERNAL_PERIMETER": [Toolpath(
                "EXTERNAL_PERIMETER", 30.5+i*.2, list(original.points), closed=True)]
        }) for i in range(3)]
        usual = emit_gcode(layers, SlicerConfig(bed_radius=30, retract_mm=0))
        aligned = emit_gcode(layers, SlicerConfig(bed_radius=30, retract_mm=0,
            seam_policy="aligned", seam_angle_deg=25))
        a = parse_extrusion_segments(usual)
        b = parse_extrusion_segments(aligned)
        # Each loop's deposited surface length is invariant even where the
        # edge is split into an additional move; don't compare raw G-code counts.
        for layer_index in range(3):
            for segments in (a, b):
                items = [item for item in segments if item[0] == layer_index]
                r = layers[layer_index].radius
                total = sum(math.hypot(x1 - x0, math.radians(y1-y0)*r)
                            for _, _, x0, y0, x1, y1 in items)
                self.assertAlmostEqual(total, length(layers[layer_index].paths_by_role['EXTERNAL_PERIMETER'][0]), delta=.002)
        self.assertIn(";SEAM_POLICY:aligned", aligned)
        # The seam stays at 25 degrees in the continuous Y frame on each
        # layer, with no layer-local rotary reset.
        for snippet in aligned.split(";LAYER:")[1:]:
            travels = [line for line in snippet.splitlines() if line.startswith("G0 X")]
            self.assertTrue(travels)
            self.assertFalse(any(line.startswith("G92.9 Y") for line in snippet.splitlines()))
            y = float(travels[0].split("Y", 1)[1].split()[0])
            self.assertAlmostEqual(y % 360, 25, places=2)


    def test_axial_target_picks_exact_position_on_constant_theta_edge(self):
        original = rect(theta_min=.1, theta_max=.9)
        cfg = SlicerConfig(seam_policy="aligned",
                           seam_angle_deg=math.degrees(.1),
                           seam_axial_target_mm=7.25)
        for turn in (-3, 0, 6):
            shifted = Toolpath(original.role, original.radius,
                [(t + turn*2*math.pi, x) for t,x in original.points],
                closed=True)
            result = _orient_path_for_current_position(
                shifted, GcodeState(theta_accum=0.8, axial=0), cfg)
            self.assertAlmostEqual(result.points[0][0] % (2*math.pi), .1, places=9)
            self.assertAlmostEqual(result.points[0][1], 7.25, places=9)
            self.assertAlmostEqual(length(result), length(shifted), places=8)

    def test_axial_exclusion_uses_printable_upper_edge_not_nearby_vertex(self):
        original = rect(theta_min=.1, theta_max=.9)
        cfg = SlicerConfig(seam_policy="aligned",
                           seam_angle_deg=math.degrees(.1),
                           seam_axial_target_mm=6,
                           seam_exclude_axial_mm=(0, 8))
        p = _orient_path_for_current_position(original, GcodeState(), cfg)
        self.assertAlmostEqual(p.points[0][0], .1, places=8)
        self.assertAlmostEqual(p.points[0][1], 8, places=8)
        self.assertAlmostEqual(length(p), length(original), places=8)
        with self.assertRaisesRegex(ValueError,"No eligible perimeter seam"):
            _orient_path_for_current_position(
                original, GcodeState(),
                SlicerConfig(seam_exclude_axial_mm=(-1, 11)))

    def test_variable_width_seam_cut_and_vertex_preserve_all_extrusion(self):
        from cyl_slicer_v2.gcode import _segment_width_ratio, emit_toolpath
        original = rect(theta_min=.1, theta_max=.9)
        original.segment_widths = [.22, .32, .42, .30]
        cfg = SlicerConfig(retract_mm=0, seam_policy="aligned",
                           seam_angle_deg=math.degrees(.1),
                           seam_axial_target_mm=6.25,
                           cooling_min_layer_time_s=0)
        expected = sum(
            cylindrical_distance(a,b,original.radius) *
            cfg.extrusion_per_mm * _segment_width_ratio(w,cfg)
            for a,b,w in zip(original.points,original.points[1:],original.segment_widths)
        )
        for config,position in (
            (cfg,(.1,6.25)),
            (SlicerConfig(retract_mm=0),(.9,0.0)),
            (SlicerConfig(retract_mm=0,seam_policy="corner",
                          seam_angle_deg=math.degrees(.9)),(.9,0.0)),
        ):
            state=GcodeState(theta_accum=position[0],axial=position[1],
                             radial=original.radius+.08)
            path = _orient_path_for_current_position(original,state,config)
            self.assertIsNotNone(path.segment_widths)
            self.assertEqual(len(path.segment_widths),len(path.points)-1)
            actual = sum(
                cylindrical_distance(a,b,original.radius) *
                cfg.extrusion_per_mm * _segment_width_ratio(w,cfg)
                for a,b,w in zip(path.points,path.points[1:],path.segment_widths)
            )
            self.assertAlmostEqual(actual,expected,places=8)
            self.assertAlmostEqual(length(path),length(original),places=8)
            emitted = emit_toolpath(path,
                 GcodeState(theta_accum=path.points[0][0],
                            axial=path.points[0][1],radial=original.radius+.08),
                 config, 0, retract_after=False)
            import re
            moves=[line for line in emitted if line.startswith("G1 X") and " E" in line]
            self.assertTrue(moves)
            final_e=float(re.search(r"\\bE(-?\\d+(?:\\.\\d+)?)",moves[-1]).group(1))
            self.assertAlmostEqual(final_e,expected,delta=2e-4)

    def test_axial_options_have_profile_parity_and_bad_values_rejected(self):
        from cyl_slicer_v2.cli import build_parser, _profile_value
        for profile in ("quality","library","legacy-cli"):
            args = build_parser().parse_args(["--profile",profile])
            self.assertIsNone(_profile_value(args,"seam_axial_target_mm"))
            self.assertIsNone(_profile_value(args,"seam_exclude_axial_mm"))
            explicit = build_parser().parse_args(["--profile",profile,
                "--seam-axial-target-mm","-5.5",
                "--seam-exclude-axial-mm","2","7"])
            self.assertEqual(_profile_value(explicit,"seam_axial_target_mm"),-5.5)
            self.assertEqual(tuple(explicit.seam_exclude_axial_mm),(2,7))
        for fields in (
            {"seam_axial_target_mm": math.nan},
            {"seam_axial_target_mm": math.inf},
            {"seam_exclude_axial_mm": (4,4)},
            {"seam_exclude_axial_mm": (9,1)},
            {"seam_exclude_axial_mm": (0,math.inf)},
            {"seam_exclude_axial_mm": (math.nan,2)},
        ):
            with self.subTest(fields=fields),self.assertRaises(ValueError):
                SlicerConfig(**fields).validate()

    def test_invalid_seam_inputs(self):
        cases = (
            {"seam_policy": "scatter"},
            {"seam_angle_deg": 360},
            {"seam_angle_deg": math.nan},
            {"seam_exclude_deg": (10, 10)},
            {"seam_exclude_deg": (-1, 30)},
            {"seam_exclude_deg": (350, math.inf)},
        )
        for fields in cases:
            with self.subTest(fields=fields), self.assertRaises(ValueError):
                SlicerConfig(**fields).validate()


if __name__ == "__main__":
    unittest.main()
