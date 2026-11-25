# SPDX-License-Identifier: MIT
"""Compare a fresh run of code/check_dominance.py with results/check_dominance.txt under the declared policy.

usage: compare_output.py NEW RECORDED

Policy (declared before tagging; the same in README.md, ENVIRONMENT.txt and REPLAY_RECEIPT.md):
- Part A (exact arithmetic): every line identical.
- Part B (floating point): every line identical once its decimal numbers are masked, so all PASS/FAIL words, labels and
  the final line "All checks passed." match exactly; each printed decimal agrees within 1e-5 relative or one unit in its
  last printed digit, whichever is larger (the script rounds to 6 decimals, 3 for one ratio).
- The two files have the same number of lines and the last line is "All checks passed.".
Exits 1 on any difference.
"""
import re
import sys

NUM = re.compile(r'-?\d+\.\d+')


def main():
    new = open(sys.argv[1], encoding='utf-8').read().splitlines()
    old = open(sys.argv[2], encoding='utf-8').read().splitlines()
    bad = []
    if len(new) != len(old):
        bad.append(f'line count {len(new)} != {len(old)}')
    if not new or new[-1] != 'All checks passed.':
        bad.append('last line is not "All checks passed."')
    part_b = False
    exact = toler = 0
    for i, (a, b) in enumerate(zip(new, old), 1):
        if b.startswith('Part B'):
            part_b = True
        if a == b:
            exact += 1
            continue
        if not part_b:
            bad.append(f'Part A line {i} differs: {a!r} != {b!r}')
            continue
        if NUM.sub('#', a) != NUM.sub('#', b):
            bad.append(f'Part B line {i} differs outside numbers: {a!r} != {b!r}')
            continue
        for x, y in zip(NUM.findall(a), NUM.findall(b)):
            ulp = 10.0 ** -len(y.split('.')[1])
            tol = max(1e-5 * abs(float(y)), ulp) * (1 + 1e-9)
            if abs(float(x) - float(y)) > tol:
                bad.append(f'Part B line {i}: {x} vs {y} beyond tolerance {tol:g}')
        toler += 1
    if bad:
        print('DIFFERENCES:')
        for m in bad:
            print(' ', m)
        sys.exit(1)
    print(f'OK: {exact} lines identical, {toler} Part B lines within the declared tolerance; last line "All checks passed."')


if __name__ == '__main__':
    main()
