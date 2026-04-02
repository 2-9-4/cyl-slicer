#!/usr/bin/env python3

import re
from decimal import Decimal, getcontext

getcontext().prec = 50

# Patterns
Y_PATTERN = re.compile(r'(?P<prefix>\bY)(?P<value>[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
LAYER_PATTERN = re.compile(r';LAYER:-?\d+')

def process_gcode(input_path, output_path):
    abs_mode = True  # G90 = absolute, G91 = relative
    current_y = Decimal(0)

    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            stripped = line.strip()

            # Detect positioning mode
            if 'G90' in stripped:
                abs_mode = True
            elif 'G91' in stripped:
                abs_mode = False

            # Reset accumulator on layer marker
            if LAYER_PATTERN.search(line):
                current_y = Decimal(0)

            def replace_y(match):
                nonlocal current_y

                value_str = match.group("value")
                prefix = match.group("prefix")

                try:
                    val = Decimal(value_str)

                    if abs_mode:
                        # Absolute: target position
                        current_y = val
                    else:
                        # Relative: increment
                        current_y += val

                    wrapped = current_y % Decimal(360)

                    # Emit wrapped value (preserve numeric fidelity)
                    return f"{prefix}{format(wrapped, 'f')}"

                except Exception:
                    return match.group(0)

            new_line = Y_PATTERN.sub(replace_y, line)
            outfile.write(new_line)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Wrap cumulative Y axis modulo 360 with layer resets")
    parser.add_argument("input", help="Input G-code file")
    parser.add_argument("output", help="Output G-code file")

    args = parser.parse_args()
    process_gcode(args.input, args.output)