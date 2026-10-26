from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _ramped_extrusion_steps, emit_toolpath
from cyl_slicer_v2.models import GcodeState, Toolpath


def run(path: Toolpath, cfg: SlicerConfig):
    state = GcodeState(theta_accum=path.points[0][0],
                       axial=path.points[0][1], radial=path.radius+.08)
    lines = emit_toolpath(path, state, cfg, 0, retract_after=False)
    moves = [line for line in lines if line.startswith("G1 X") and " E" in line]
    return moves, state


def feed(line: str):
    return float(re.search(r" F([0-9.]+)", line).group(1))


class SeamFinishRampTests(unittest.TestCase):
    @staticmethod
    def loop(*, closed=True, role="EXTERNAL_PERIMETER"):
        return Toolpath(role, 30,
                        [(0, 0), (.1, 0), (.1, 8), (0, 8), (0, 0)],
                        closed=closed)

    def test_default_identical_to_original_unsplit_loop(self):
        path = self.loop()
        cfg = SlicerConfig(retract_mm=0, seam_end_ramp_mm=0,
                           extrusion_feed_ramp_mm=0)
        moves, state = run(path, cfg)
        self.assertEqual(len(moves), 4)
        self.assertAlmostEqual(state.extrusion, sum(
            math.hypot((b[0]-a[0])*path.radius, b[1]-a[1])
            for a, b in zip(path.points, path.points[1:])) * cfg.extrusion_per_mm)

    def test_closed_wall_finishes_slowly_without_any_volume_loss(self):
        path = self.loop()
        normal = SlicerConfig(feed_mm_per_min=1200, retract_mm=0)
        tuned = SlicerConfig(feed_mm_per_min=1200, retract_mm=0,
                             seam_end_ramp_mm=4, seam_end_ramp_fraction=.5)
        old, old_state = run(path, normal)
        new, new_state = run(path, tuned)
        self.assertGreater(len(new), len(old))
        self.assertLess(feed(new[-1]), 1200)
        self.assertGreaterEqual(feed(new[-1]), 600)
        self.assertEqual(feed(new[0]), feed(old[0]))
        self.assertAlmostEqual(new_state.extrusion, old_state.extrusion, places=9)
        self.assertEqual(new_state.theta_accum, old_state.theta_accum)
        self.assertEqual(new_state.axial, old_state.axial)

    def test_only_closed_perimeters_receive_finish_tuning(self):
        cfg = SlicerConfig(seam_end_ramp_mm=5, seam_end_ramp_fraction=.4)
        for path in (self.loop(closed=False),
                     self.loop(role="SUPPORT"),
                     self.loop(role="BRIDGE")):
            with self.subTest(role=path.role, closed=path.closed):
                moves, _ = run(path, cfg)
                self.assertEqual(len(moves), 4)

    def test_start_and_finish_ramps_split_short_loop_without_duplicating_e(self):
        path = Toolpath("EXTERNAL_PERIMETER", 30,
                        [(0, 0), (0, 2), (.05, 2), (.05, 0), (0, 0)],
                        closed=True)
        base = SlicerConfig(extrusion_feed_ramp_mm=0, seam_end_ramp_mm=0)
        cfg = SlicerConfig(extrusion_feed_ramp_mm=8,
                           extrusion_feed_ramp_start_fraction=.5,
                           seam_end_ramp_mm=8, seam_end_ramp_fraction=.6)
        old, old_state = run(path, base)
        new, new_state = run(path, cfg)
        self.assertGreater(len(new), len(old))
        self.assertAlmostEqual(new_state.extrusion, old_state.extrusion, places=9)
        self.assertTrue(all(feed(line) > 0 for line in new))
        steps = _ramped_extrusion_steps(path.points, path.radius, cfg, seam_end=True)
        self.assertEqual(steps[-1][0], path.points[-1])
        self.assertEqual(len(steps), len(new))

    def test_variable_width_loop_preserves_segment_areas(self):
        path = Toolpath("INTERNAL_PERIMETER", 30,
                        [(0, 0), (0, 4), (.2, 4), (.2, 0), (0, 0)],
                        closed=True, segment_widths=[.23,.34,.42,.29])
        base = SlicerConfig(retract_mm=0)
        cfg = SlicerConfig(retract_mm=0, seam_end_ramp_mm=6,
                           seam_end_ramp_fraction=.6)
        before, bs = run(path, base)
        after, as_ = run(path, cfg)
        self.assertGreater(len(after), len(before))
        self.assertAlmostEqual(as_.extrusion, bs.extrusion, places=9)

    def test_cli_defaults_and_invalid_values(self):
        defaults = SlicerConfig()
        self.assertEqual(defaults.seam_end_ramp_mm, 0)
        args = build_parser().parse_args(["--seam-end-ramp-mm", "3",
                                          "--seam-end-ramp-fraction", ".45"])
        self.assertEqual(_profile_value(args, "seam_end_ramp_mm"), 3)
        self.assertEqual(_profile_value(args, "seam_end_ramp_fraction"), .45)
        for kwargs in ({"seam_end_ramp_mm": -1},
                       {"seam_end_ramp_mm": math.inf},
                       {"seam_end_ramp_fraction": 0},
                       {"seam_end_ramp_fraction": math.nan},
                       {"seam_end_ramp_fraction": 1.1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                SlicerConfig(**kwargs).validate()


if __name__ == "__main__":
    unittest.main()
