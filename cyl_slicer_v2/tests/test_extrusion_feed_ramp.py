from __future__ import annotations

import math
import re
import unittest

from cyl_slicer_v2.cli import _profile_value, build_parser
from cyl_slicer_v2.config import SlicerConfig
from cyl_slicer_v2.gcode import _ramped_extrusion_steps, _surface_feed_for_marlin, emit_toolpath
from cyl_slicer_v2.models import GcodeState, Toolpath


def extruded(lines):
    return [line for line in lines if line.startswith("G1 X") and " E" in line]


def feed(line):
    return float(re.search(r" F([\d.]+)",line).group(1))


def emitted_e(line):
    return float(re.search(r" E([\d.]+)",line).group(1))


def run(path, cfg):
    state = GcodeState(theta_accum=path.points[0][0],
                       axial=path.points[0][1], radial=path.radius+.08)
    code = emit_toolpath(path,state,cfg,0,retract_after=False)
    return extruded(code),state


class ExtrusionFeedRampTests(unittest.TestCase):
    def test_disabled_preserves_exact_moves_and_e(self):
        c=SlicerConfig(extrusion_feed_ramp_mm=0,retract_mm=0)
        path=Toolpath("EXTERNAL_PERIMETER",30,[(0,0),(0,12)])
        moves,state=run(path,c)
        self.assertEqual(len(moves),1)
        self.assertEqual(feed(moves[0]), c.feed_mm_per_min)
        self.assertAlmostEqual(state.extrusion,12*c.extrusion_per_mm,places=9)
        self.assertNotIn("M900", "\n".join(moves))

    def test_four_ramp_stages_no_short_blocks_and_unchanged_deposited_volume(self):
        original=Toolpath("SOLID_INFILL",30,[(0,0),(0,10)])
        base=SlicerConfig(feed_mm_per_min=1200,retract_mm=0)
        ramp=SlicerConfig(feed_mm_per_min=1200,retract_mm=0,
                         extrusion_feed_ramp_mm=4,
                         extrusion_feed_ramp_start_fraction=.5)
        old,old_state=run(original,base)
        commands,new_state=run(original,ramp)
        self.assertEqual(len(old),1)
        self.assertEqual(len(commands),5)
        speeds=[feed(x) for x in commands]
        self.assertEqual(speeds[-1],1200)
        self.assertEqual(speeds,sorted(speeds))
        self.assertAlmostEqual(speeds[0],675,places=2)
        self.assertAlmostEqual(speeds[1],825,places=2)
        self.assertAlmostEqual(speeds[2],975,places=2)
        self.assertAlmostEqual(speeds[3],1125,places=2)
        self.assertTrue(all(.8<= b-a for a,b in
          zip([0]+[float(re.search(r"X([\d.]+)",x).group(1)) for x in commands[:-1]],
              [float(re.search(r"X([\d.]+)",x).group(1)) for x in commands])))
        self.assertAlmostEqual(new_state.extrusion,old_state.extrusion,places=9)
        self.assertAlmostEqual(emitted_e(commands[-1]),emitted_e(old[-1]),places=5)

    def test_mixed_axis_and_rotary_feed_changes_preserve_physical_target_ratio(self):
        radius=60.
        for path in (
            Toolpath("SOLID_INFILL",radius,[(0,0),(.2,0)]),
            Toolpath("SOLID_INFILL",radius,[(0,0),(.12,9)]),
        ):
            with self.subTest(path=path):
                c=SlicerConfig(feed_mm_per_min=1800,retract_mm=0,
                    extrusion_feed_ramp_mm=4,
                    extrusion_feed_ramp_start_fraction=.5)
                steps=_ramped_extrusion_steps(path.points,radius,c)
                lines,state=run(path,c)
                self.assertGreater(len(lines),1)
                self.assertEqual(len(steps),len(lines))
                last=path.points[0]
                calculated=0.
                feeds=[]
                for ((pos,_,fraction),line) in zip(steps,lines):
                    dt=pos[0]-last[0]
                    dx=pos[1]-last[1]
                    expected=_surface_feed_for_marlin(dx,dt,radius,1800*fraction)
                    self.assertAlmostEqual(feed(line),expected,delta=.08)
                    calculated += math.hypot(dx,dt*radius)
                    feeds.append(fraction)
                    last=pos
                self.assertAlmostEqual(calculated,12 if path.points[-1][1]==0 else math.hypot(7.2,9),places=8)
                self.assertAlmostEqual(state.extrusion,calculated*c.extrusion_per_mm,places=9)
                self.assertEqual(feeds,sorted(feeds))

    def test_variable_width_bead_areas_and_total_e_are_invariant(self):
        path=Toolpath("GAP_FILL",30,[(0,0),(0,2),(0,4),(0,6)],
                      segment_widths=[.24,.37,.30])
        c0=SlicerConfig(retract_mm=0,extrusion_feed_ramp_mm=0)
        c1=SlicerConfig(retract_mm=0,extrusion_feed_ramp_mm=4)
        old,old_state=run(path,c0)
        new,new_state=run(path,c1)
        self.assertGreater(len(new),len(old))
        self.assertAlmostEqual(new_state.extrusion,old_state.extrusion,places=9)
        self.assertAlmostEqual(emitted_e(new[-1]),emitted_e(old[-1]),places=5)
        self.assertEqual(len(old),3)

    def test_role_speed_cap_is_never_bypassed(self):
        path=Toolpath("SUPPORT_INTERFACE",30,[(0,0),(0,10)])
        c=SlicerConfig(feed_mm_per_min=2700,
            support_interface_print_feed_mm_per_min=600,
            extrusion_feed_ramp_mm=3,extrusion_feed_ramp_start_fraction=.5,
            retract_mm=0)
        commands,_=run(path,c)
        speeds=[feed(line) for line in commands]
        self.assertGreater(len(speeds),1)
        self.assertTrue(all(0 < x <=600 for x in speeds))
        self.assertEqual(speeds[-1],600)

    def test_cli_profile_parity_and_validation(self):
        assert SlicerConfig().extrusion_feed_ramp_mm==0
        for profile in ("quality","library","legacy-cli"):
            args=build_parser().parse_args(["--profile",profile])
            self.assertEqual(_profile_value(args,"extrusion_feed_ramp_mm"),0)
            self.assertEqual(_profile_value(args,"extrusion_feed_ramp_start_fraction"),.65)
        args=build_parser().parse_args(["--extrusion-feed-ramp-mm","4",
                                         "--extrusion-feed-ramp-start-fraction",".45"])
        self.assertEqual(_profile_value(args,"extrusion_feed_ramp_mm"),4)
        self.assertEqual(_profile_value(args,"extrusion_feed_ramp_start_fraction"),.45)
        for fields in ({"extrusion_feed_ramp_mm":-1},
                       {"extrusion_feed_ramp_mm":math.inf},
                       {"extrusion_feed_ramp_start_fraction":0},
                       {"extrusion_feed_ramp_start_fraction":1.2}):
            with self.subTest(fields=fields),self.assertRaises(ValueError):
                SlicerConfig(**fields).validate()


if __name__=="__main__":
    unittest.main()
