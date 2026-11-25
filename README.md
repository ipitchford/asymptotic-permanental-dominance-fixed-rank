# Asymptotic permanental dominance at fixed Gram rank

*Every non-trivial symmetry block averages below the permanent once the matrix is large; the standard block is the
worst, at (d-1) log n / n*

Version 0.1.0-candidate. DOI [10.5281/zenodo.24536861](https://doi.org/10.5281/zenodo.24536861) (all versions:
[10.5281/zenodo.24536860](https://doi.org/10.5281/zenodo.24536860)). Release page:
https://evidencepress.org/releases/asymptotic-permanental-dominance-fixed-rank/

**Status: unrefereed research candidate.** Author: Anonymous. Text CC0-1.0; code MIT (see [LICENSE.md](LICENSE.md)).

## What this is about

**The problem.** The permanent of a square matrix is the determinant without minus signs: add up the products of the
entries along every permutation. Take n vectors and the table A of their inner products (a positive semidefinite
"Gram" matrix). The permutations of the n vectors break into symmetry types, one for each way of drawing a Young
diagram with n boxes. Each type comes with a matrix block B_λ(A), and the all-symmetric type gives exactly the
permanent. In 1966 Lieb asked whether every other block is, on average, no bigger than the permanent; for the ordinary
blocks (the immanants) this says that the average eigenvalue of B_λ(A) is at most per A. In 2026 Lavi asked for the
same with mean squares of the eigenvalues. Write m_1 and m_2 for these two averages divided by per A and (per A)^2; both
conjectures say m_1 <= 1 and m_2 <= 1.

**Why anyone cares.** A stronger guess, "permanent on top" (the permanent beats every single eigenvalue of every
block), is false, and single eigenvalues can exceed the permanent by a factor of about log n. The averaged inequalities
are what may survive. Lavi's version would also imply Chollet's permanent inequality. In optics, Gram matrices of
rank d describe photons with d internal modes, and the blocks measure symmetry sectors of the photons.

**What was known.**
- Exact checks for small sizes: Lieb's version (ordinary immanants) up to 15 vectors (Pate through 13; Li for the
  remaining shapes at 14 and 15), Lavi's up to 5 (Lavi).
- Infinite families of shapes: Pate (1999, Theorem 2) proved Lieb's version for all shapes (p, q^u, r^v, 2^s, 1^t)
  with u <= 2 and v <= 1. Since blocks with more than d rows vanish at rank d, this already gives Lieb's version for
  every n when the vectors span at most three dimensions. Pate's Theorem 4 gives eventual dominance for the shapes
  (n+s, n^k); Li proves it for every fixed tail once the first row is long.
- Real two-dimensional vectors: both inequalities hold for every n, by the rank-two permanent-on-top theorem.

**What is new.** Fix the dimension d of the vectors (the Gram rank).
1. **Fixed rank** (Theorem 1.1). As n grows, both normalised moments tend to zero, uniformly over every non-trivial
   symmetry type and every input of rank at most d. So both inequalities hold at each fixed rank for all sufficiently
   large n. The worst type is eventually the standard one (one box moved to the second row); for complex vectors its
   moments are asymptotic to (d-1) log n / n and (d-1)(log n)^2 / n. Over the reals with d >= 3 the orders are the same
   with constants between 1/2 (respectively 1/4) and d-1; at real rank two both are of order 1/n.
2. **Growing rank** (Theorem 1.2). Both moments still tend to zero when d log n = o(sqrt n); the first moment does
   whenever d log n = o(n).
3. **The standard type at every n and d** (Theorem 1.3). Explicit bounds m_1 <= (6 sqrt(e n d) - 1)/(n - 1) and
   m_2 <= (6 sqrt(e n d)(log n + 5) - 1)/(n - 1), with a real example giving (d-1)/(n-1). So the standard first moment
   tends to zero exactly when d = o(n), and Lavi's inequality for the standard type holds whenever
   6 sqrt(e n d)(log n + 5) <= n (for d = 2, from n = 48,836 on).

**Limits.** "Sufficiently large n" is not computed, so no new case of a fixed size is proved. Lieb's version is new
only from dimension four, and there only for four-row shapes with fourth row at least 3 and second row longer than the
third; everything else at dimension at most four is due to Pate. Only the symmetric-group (immanant) form of Lieb's
conjecture is treated, not the subgroup form. Real-vector constants are known only up to a factor. Lavi's identity
gives an eventual Chollet inequality for pairs of rank at most d, but direct estimates (Ma) are much stronger, so it is
not counted as a contribution.

**Idea of the proof.** Each block is bounded by a simpler "ordered-partition" matrix whose trace is estimated by
averaging over random Gaussian directions. A new strict inequality for logarithms of sums of exponential random
variables (Lemma 3.1) shows that, whenever two rows of the Young diagram are long, the dimension of the symmetry type
grows exponentially faster than this bound. Diagrams with one dominant row are handled by counting. The standard type
is pinned down with the sharp Bapat-Sunder growth from the earlier Evidence Press release
[10.5281/zenodo.24509057](https://doi.org/10.5281/zenodo.24509057) and a self-contained sharp first-trace constant
(Proposition 4.3).

**Credit.** Rank at most three, and every rank-four shape except four rows with fourth row at least 3 and second row
longer than the third: Pate (1999), Theorem 2. Eventual dominance for the shapes (n+s, n^k): Pate, Theorem 4. Orders
up to 15: Pate and Li. Lavi's conjecture through order 5 and the identities used for m_2: Lavi. Real rank two: the
liqsweep preprint, with a complete proof in the Bapat-Sunder release. The eigenvalue bound of Lemma 2.4 is derived from
Anari, Gurvits, Oveis Gharan and Saberi (Theorem 1 and eq. (2)). The d - 1 trace constant and the growth bounds are this
paper's.

## Contents

| File | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (17 pages; pdflatex three times, no figures, no BibTeX). Section 1 is a plain introduction with a status table and a roadmap. |
| [code/check_dominance.py](code/check_dominance.py) | Finite checks of the ingredients of the proofs (Python standard library only). Part A, exact arithmetic: positivity, vanishing and trace identities of all blocks for n <= 6, the Young and Lieb trace bounds, the Plancherel identity, the real rank-two bound Tr F_A <= 4 per A, the block example J_a + J_(n-a) with both standard moments 1/(n-1), and the repeated orthonormal basis with Tr C_1 = d per A. Part B, floating point with a stated tolerance: the scalar inequality of Lemma 3.1, the trace bounds of Section 4, the chain of inequalities of Section 8 and the test-vector bound of Section 6 on small frames. It is a check of the ingredients, not part of any proof. |
| [results/check_dominance.txt](results/check_dominance.txt) | Its recorded output. |
| [refs.bib](refs.bib) | Bibliography (12 entries, from Crossref, DataCite and arXiv records). |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits. |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | What the internal review found and what was done. |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-art map, sources read and limitation triage. |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment. |
| [.github/workflows/replay.yml](.github/workflows/replay.yml), [.github/ci/compare_output.py](.github/ci/compare_output.py) | CI replay on Linux and the comparison policy. |
| [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Citation metadata, licences, SHA-256 sums. |

## Replay

From the package root (Python 3.12 or later, standard library only):

```sh
python3 -I -B -O code/check_dominance.py > out.txt    # also once without -O
python3 -I .github/ci/compare_output.py out.txt results/check_dominance.txt
```

The script never uses `assert`; it prints `FAIL: <label>` and exits 1 if any check fails, also under `python -O`. It
takes a few seconds. The last line of a successful run is `All checks passed.`

**Comparison policy** (declared before tagging; the same in the CI workflow, ENVIRONMENT.txt and REPLAY_RECEIPT.md).
Part A is exact and must match line by line. Part B is floating point: every PASS/FAIL line, every label and the final
line must match exactly, and each printed decimal must agree within 1e-5 relative or one unit in its last printed
digit, whichever is larger (values are printed to 6 decimals, one ratio to 3). On the release host (macOS) the output
was byte-identical to `results/check_dominance.txt` with and without `-O`, under CPython 3.14.7 and 3.12.11.

## Assurance

Written proofs; the script checks ingredients on finite examples. Internal AI-assisted review only: a combined
five-role Claude review panel, a proof audit of every proof by the lead, one revision round, a verification re-review
with a proof audit of the new and rewritten proofs, and a clean replay with a negative control. No Astra/Codex review
took place. This is not peer review, independent reproduction or formal verification. See
[REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).

## Cite

Anonymous (2026). *Asymptotic permanental dominance at fixed Gram rank*. Version 0.1.0-candidate. Evidence Press.
https://doi.org/10.5281/zenodo.24536861
