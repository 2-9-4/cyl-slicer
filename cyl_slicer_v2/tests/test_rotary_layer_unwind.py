from __future__ import annotations

import re
import unittest

from cyl_slicer_v2.tests.test_marlin_rotary_precision import MarlinRotaryModel
from cyl_slicer_v2.unwind_gcode import add_layer_unwinds


PROGRAM = """;start of slice
G90
M82
M400
G92.9 Y0
G92 E0 Y0
;LAYER:0
G0 Z0.2000 F9000
G0 X1 Y-180.000 F3600
G1 X2 Y-190.000 E1.00000 F1200
G1 E-0.40000 F1200 ; retract
;LAYER:1
G0 Z0.4000 F9000
G92 E0 ; reset extrusion accumulator for layer
G0 X1 Y-180.000 F3600
G1 E1.40000 F1200 ; unretract
G1 X2 Y-170.000 E2.00000 F1200
G1 E0.60000 F1200 ; retract
;end of slice
"""


class RotaryLayerUnwindTests(unittest.TestCase):
    def test_physical_return_preserves_print_moves_and_extrusion(self):
        altered = add_layer_unwinds(PROGRAM)
        self.assertIn(
            "G0 Z1.0000 F9000.0 ; lift before physical rotary unwind\n"
            "M400 ; complete lift before rotary unwind\n"
            "G0 Y0.000 F3600.0 ; physically return to Y0\n"
            "M400 ; complete physical rotary unwind\n;LAYER:1",
            altered,
        )
        self.assertEqual(altered.count(";PHYSICAL_LAYER_UNWIND:1"), 1)
        self.assertEqual(altered.count("G92.9 Y"), 1)
        self.assertEqual(altered.count("G92 E0 Y0"), 1)
        old_moves = [line for line in PROGRAM.splitlines() if line.startswith(("G0 ", "G1 "))]
        new_moves = [line for line in altered.splitlines() if line.startswith(("G0 ", "G1 "))]
        self.assertEqual([line for line in new_moves if "; physically return" not in line
                          and "; lift before physical" not in line], old_moves)

        def printed_steps(code: str) -> list[int]:
            model = MarlinRotaryModel()
            result = []
            for line in code.splitlines():
                model.feed(line)
                if line.startswith("G1 X") and re.search(r"\bE[-\d]", line):
                    result.append(model.physical_steps)
            return result

        self.assertEqual(printed_steps(PROGRAM), printed_steps(altered))

    def test_rejects_unretracted_or_rebased_input(self):
        with self.assertRaisesRegex(ValueError, "retracted filament"):
            add_layer_unwinds(PROGRAM.replace("G1 E-0.40000 F1200 ; retract\n", ""))
        with self.assertRaisesRegex(ValueError, "no Y coordinate resets"):
            add_layer_unwinds(PROGRAM.replace(
                "G92 E0 ; reset extrusion accumulator for layer",
                "G92 E0 Y0 ; obsolete rotary reset",
            ))


if __name__ == "__main__":
    unittest.main()
