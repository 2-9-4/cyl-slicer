from dataclasses import replace
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from shapely import affinity
from shapely.geometry import GeometryCollection, LineString, box

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _should_retract_after_path, emit_gcode
from cyl_slicer_v2.layers import propagate_support, segment_support_wedges
from cyl_slicer_v2.models import GcodeState, LayerSlice, Toolpath
from cyl_slicer_v2.periodic import periodic_copies
from cyl_slicer_v2.toolpaths import build_layer_toolpaths, gap_support_slices


EMPTY = GeometryCollection()
TURN = 2 * math.pi
RING = box(0, 0, TURN, 8)


class SupportWedgeTests(unittest.TestCase):
    def cfg(self, **overrides):
        return replace(SlicerConfig(support_wedge_count=12,
                                   support_min_island_area=0,
                                   support_min_region_area=0,
                                   support_bearing_min_area=0), **overrides)

    def assert_slots_empty(self, shape, config, radius):
        """Check the same radial opening at every boundary, including seam."""
        phase = math.radians(config.support_wedge_phase_deg)
        for i in range(config.support_wedge_count + 1):
            boundary = (phase + i * TURN / config.support_wedge_count) % TURN
            half_gap = config.support_wedge_gap_mm / (2 * config.bed_radius)
            for shift in (-TURN, 0, TURN):
                slot = box(boundary + shift - half_gap + 1e-8, -1,
                           boundary + shift + half_gap - 1e-8, 9)
                self.assertLess(shape.intersection(slot).area * radius, 1e-7)

    def test_fixed_boundaries_periodic_windows_and_physical_gap(self):
        for phase in (0, 7, -353, 367):
            c = self.cfg(support_wedge_phase_deg=phase)
            for radius in (30.08, 45, 60):
                with self.subTest(phase=phase, radius=radius):
                    result = segment_support_wedges(RING, radius, c)
                    self.assert_slots_empty(result, c, radius)
                    self.assertEqual(len(periodic_copies(result).geoms),
                                     c.support_wedge_count * 3 + (1 if phase else 0))
                    # Geometry is radially aligned, so gap arcs widen with R.
                    self.assertTrue(result.equals(segment_support_wedges(RING, 30.08, c)))
                    expected_removed = c.support_wedge_count * c.support_wedge_gap_mm / c.bed_radius * 8
                    self.assertAlmostEqual(RING.area - result.area, expected_removed)
            shifted = affinity.translate(RING, xoff=-TURN)
            self.assertTrue(affinity.translate(segment_support_wedges(shifted, 45, c), xoff=TURN)
                            .equals_exact(segment_support_wedges(RING, 45, c), 1e-12))

    def test_support_propagation_keeps_model_contact_gap_and_unsegmented_carry(self):
        def stack():
            layers = [LayerSlice(i, 30.08 + .16*i, EMPTY, EMPTY, EMPTY, EMPTY, EMPTY)
                      for i in range(7)]
            layers[-1].occupied = RING
            layers[-1].needs_support = RING
            layers[-1].support_interface_target = RING
            return layers
        baseline, actual = stack(), stack()
        c = self.cfg(support_merge_gap=10, support_contact_gap_mm=.16)
        propagate_support(baseline, replace(c, support_wedge_count=0))
        propagate_support(actual, c)
        self.assertTrue(actual[5].support.is_empty)
        for before, after in zip(baseline, actual):
            self.assertTrue(before.occupied.equals(after.occupied))
            self.assertTrue(before.needs_support.equals(after.needs_support))
            self.assertTrue(segment_support_wedges(before.support, after.radius, c).equals(after.support))
            self.assertTrue(segment_support_wedges(before.support_interface, after.radius, c)
                            .equals(after.support_interface))
            self.assert_slots_empty(after.support, c, after.radius)
            self.assert_slots_empty(after.support_interface, c, after.radius)

    def test_raft_foundation_and_partial_slabs_keep_slots_after_merging(self):
        for mode in ("fixed", "complete-terminal"):
            c = self.cfg(radial_layer_mode=mode, raft_layers=2, raft_axial_gap=0,
                         raft_merge_gap_mm=10, support_merge_gap=10)
            first = LayerSlice(0, 31.08, RING, EMPTY, EMPTY, EMPTY, EMPTY,
                               is_full_wrap=True, theta_window=(0, TURN), deposited_thickness=.1)
            foundation = gap_support_slices([first], c)
            self.assertGreater(len(foundation), c.raft_layers)
            for layer in foundation:
                with self.subTest(mode=mode, index=layer.index):
                    self.assert_slots_empty(layer.support, c, layer.radius)
                    self.assertEqual(len(layer.support.geoms), c.support_wedge_count)
                    self.assertLessEqual(layer.radius + (layer.deposited_thickness or c.layer_height)/2,
                                         first.radius - .05 - c.layer_height + 1e-9)
            self.assertTrue(foundation[0].is_raft)
            self.assertFalse(foundation[-1].is_raft)

    def test_all_support_roles_keep_full_beads_outside_slots_and_separate_islands(self):
        c = self.cfg(support_meander_max_angle_deg=360)
        for kind in ("sparse", "interface-axial", "interface-radial", "raft-axial", "raft-radial"):
            for radius in (30.08, 45):
                layer = LayerSlice(0, radius, EMPTY, EMPTY, EMPTY, EMPTY, RING,
                                   theta_window=(0, TURN), is_full_wrap=True)
                if kind.startswith("interface"):
                    layer.support_interface = RING
                    layer.support_interface_direction = kind.split("-")[1]
                if kind.startswith("raft"):
                    layer.is_raft = True
                    layer.raft_direction = kind.split("-")[1]
                paths = build_layer_toolpaths(layer, c)
                self.assertEqual(len(paths.island_paths), c.support_wedge_count)
                self.assert_slots_empty(paths.travel_region, c, radius)
                roles = paths.paths_by_role
                expected = "SUPPORT_INTERFACE" if kind.startswith("interface") else "SUPPORT"
                self.assertTrue(roles[expected])
                for path in roles[expected]:
                    # Buffer the complete physical polyline, including turns,
                    # to catch both endpoint and cross-tie bead leakage.
                    bead = LineString([(t*radius, z) for t, z in path.points]).buffer(c.nozzle_width/2)
                    bead = affinity.scale(bead, xfact=1/radius, yfact=1, origin=(0, 0))
                    self.assert_slots_empty(bead, c, radius)
                if kind == "sparse":
                    self.assertTrue(any(abs(b[0]-a[0]) < 1e-8 and abs(b[1]-a[1]) > 1
                                        for p in roles[expected] for a, b in zip(p.points, p.points[1:])))

    def test_short_transfers_across_slot_retract_including_rotary_seam(self):
        c = self.cfg()
        radius = 40
        footprint = periodic_copies(segment_support_wedges(RING, radius, c))
        for boundary in (0, TURN/12):
            dt = c.support_wedge_gap_mm/(2*c.bed_radius) + c.nozzle_width/(2*radius)
            a, b = (boundary-dt, 4), (boundary+dt, 4)
            current = Toolpath("SUPPORT", radius, [(a[0], 3), a])
            following = Toolpath("SUPPORT_INTERFACE", radius, [b, (b[0], 5)])
            state = GcodeState(theta_accum=a[0], axial=a[1])
            self.assertTrue(_should_retract_after_path(current, following, state, c, footprint))

    def test_wedge_crossing_coordinate_seam_remains_one_support_island(self):
        c = self.cfg(support_wedge_phase_deg=7)
        layer = LayerSlice(0, 40, EMPTY, EMPTY, EMPTY, EMPTY, RING,
                           theta_window=(0, TURN), is_full_wrap=True)
        result = build_layer_toolpaths(layer, c)
        self.assertEqual(len(result.island_paths), c.support_wedge_count)
        self.assert_slots_empty(result.travel_region, c, layer.radius)

    def test_disabled_mode_preserves_geometry_and_gcode(self):
        c = self.cfg(support_wedge_count=0)
        self.assertIs(segment_support_wedges(RING, 40, c), RING)
        layer = LayerSlice(0, 40, box(.2, 0, .4, 8), EMPTY, EMPTY, EMPTY, RING)
        other = replace(c, support_wedge_gap_mm=2, support_wedge_phase_deg=23)
        self.assertEqual(emit_gcode([build_layer_toolpaths(layer, c)], c),
                         emit_gcode([build_layer_toolpaths(layer, other)], other))
        on = self.cfg()
        a, b = build_layer_toolpaths(layer, c), build_layer_toolpaths(layer, on)
        for role in ("EXTERNAL_PERIMETER", "INTERNAL_PERIMETER", "GAP_FILL", "SOLID_INFILL"):
            self.assertEqual(a.paths_by_role.get(role), b.paths_by_role.get(role))

    def test_cli_defaults_and_overrides_for_all_profiles(self):
        for profile in ("quality", "legacy-cli", "library"):
            args = build_parser().parse_args(["--profile", profile])
            for field in ("support_wedge_count", "support_wedge_gap_mm", "support_wedge_phase_deg"):
                self.assertEqual(_profile_value(args, field), getattr(SlicerConfig(), field))
            args = build_parser().parse_args(["--profile", profile, "--support-wedge-count", "8",
                                             "--support-wedge-gap-mm", "1.2",
                                             "--support-wedge-phase-deg", "-10"])
            self.assertEqual(_profile_value(args, "support_wedge_count"), 8)
            self.assertEqual(_profile_value(args, "support_wedge_gap_mm"), 1.2)
            self.assertEqual(_profile_value(args, "support_wedge_phase_deg"), -10)

    def test_cli_wires_settings_into_gcode(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "wedges.gcode"
            subprocess.run([sys.executable, "-m", "cyl_slicer_v2.cli", "--fixture", "tube",
                            "--build-cylinder-radius", "9", "--max-layers", "1",
                            "--support-wedge-count", "8", "--support-wedge-gap-mm", "1.2",
                            "--support-wedge-phase-deg", "5", "--out", str(output)],
                           check=True, capture_output=True, text=True)
            code = output.read_text()
            self.assertIn(";SUPPORT_WEDGE_COUNT:8", code)
            self.assertIn(";SUPPORT_WEDGE_GAP_MM:1.2", code)
            self.assertIn(";SUPPORT_WEDGE_PHASE_DEG:5", code)
            self.assertIn(";TYPE:SUPPORT", code)

    def test_invalid_settings_rejected(self):
        for value in (-1, 1, 361, 12.5, True):
            with self.subTest(count=value), self.assertRaisesRegex(ValueError, "support_wedge_count"):
                self.cfg(support_wedge_count=value).validate()
        for value in (-1, 0, math.nan, math.inf):
            with self.subTest(gap=value), self.assertRaisesRegex(ValueError, "support_wedge_gap_mm"):
                self.cfg(support_wedge_gap_mm=value).validate()
        for value in (math.nan, math.inf):
            with self.subTest(phase=value), self.assertRaisesRegex(ValueError, "support_wedge_phase_deg"):
                self.cfg(support_wedge_phase_deg=value).validate()
        with self.assertRaisesRegex(ValueError, "one bead width"):
            self.cfg(support_wedge_count=360).validate()
        for radius in (0, -1, math.nan, math.inf):
            with self.subTest(bed_radius=radius), self.assertRaisesRegex(ValueError, "bed_radius"):
                self.cfg(bed_radius=radius).validate()


if __name__ == "__main__":
    unittest.main()
