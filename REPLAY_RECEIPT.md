# Replay receipt (release 0.1.0-candidate, 25 November 2026)

The publisher replayed the approved candidate from fresh copies before adding the DOI line. Only `manuscript.tex`
(the DOI line under the author), the rebuilt `manuscript.pdf` and the release records changed afterwards; the script and
its recorded output are byte-identical to the replayed copies (checked by SHA-256 in the frozen-archive gate).

Interpreters were printed inside the run (see [ENVIRONMENT.txt](ENVIRONMENT.txt)): `uv run --python 3.12` resolved
to CPython 3.14.7 through the host's python3 shim, so the script was also run under CPython 3.12.11.

## Runs (fresh copies, cwd = package root, macOS arm64, PYTHONDONTWRITEBYTECODE=1)

| Command | Exit | Time | Output |
|---|---|---|---|
| `python3 -I -B -O code/check_dominance.py` (CPython 3.14.7) | 0 | about 4 s | last line `All checks passed.`; stdout byte-identical to results/check_dominance.txt |
| `python3 -I -B code/check_dominance.py` (CPython 3.14.7) | 0 | about 4 s | byte-identical to results/check_dominance.txt |
| `python3 -I -B -O code/check_dominance.py` (CPython 3.12.11) | 0 | about 4 s | byte-identical to results/check_dominance.txt |

Every run wrote nothing to stderr.

## Negative control (fresh copy, `python3 -I -B -O`)

| Mutation | Result |
|---|---|
| Part A2: the standard first moment of J_a + J_(n-a) expected as 1/n instead of 1/(n-1) | `FAIL: J+J first moment` (repeated for each tested pair), `27 check(s) FAILED`, exit 1 |

## Comparison policy (declared before tagging; also in README.md, ENVIRONMENT.txt and the CI workflow)

- Part A (exact arithmetic): every line identical.
- Part B (floating point): every line identical outside its printed decimals, so every PASS/FAIL line, every label and
  the final line `All checks passed.` match exactly; each printed decimal within 1e-5 relative or one unit in its last
  printed digit, whichever is larger. The script rounds printed values to 6 decimals (3 for one ratio).

The comparison script is `.github/ci/compare_output.py`; it was tested against the release output (identical), against
the negative-control output (rejected) and against a copy with one decimal moved by ten units in the last digit
(rejected).
