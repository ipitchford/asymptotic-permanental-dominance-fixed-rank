# Review summary (public, 0.1.0-candidate)

This records the internal review of the candidate: who looked, what they raised, and what was done. It is
AI-assisted internal review, not peer review or independent reproduction.

## What ran

1. **Panel (round 1).** Five reviewer roles, editor, methodology (proof audit), domain (literature, attribution and
   novelty), perspective (readability for non-specialists, optics and significance) and devil's advocate, run as one
   combined Claude pass, with a novelty search, a check of the cited sources and numerical recomputation of the key
   inequalities.
2. **Lead proof audit** of every proof in the source manuscript, before the panel reported.
3. **One revision round** with extensions: a self-contained sharp first-trace constant d - 1 (Proposition 4.3), which
   removed a dependence on an unpublished companion manuscript; an explicit every-(n, d) bound on the standard quadratic
   moment and an explicit range where Lavi's inequality holds for the standard representation (Theorem 1.3); a new
   context-first introduction with a status table; a standard-library check script with explicit failure exits.
4. **Verification re-review and proof audit** of all new and rewritten proofs: two must-fix items (a wrong number in
   the d = 2 example, 7.3e4 rather than 3.9e4, whose conclusion still held; the missing recorded output file) and the
   should-fix items applied.
5. **Clean replay** by the publisher from fresh copies, with and without `python -O` and under two interpreters, plus
   a negative control.

No second-model (Astra/Codex) adversarial review took place for this release.

## Sources

Lavi's preprint (Conjecture 1.1, Lemma 2.2, Theorem 2.1, eqs. (2.5) and (2.6), Theorem 4.1), Li (Theorems 1.3 and 1.8,
Section 2.2), Pate (1999; Theorems 1 to 4), Anari-Gurvits-Oveis Gharan-Saberi (Theorem 1, eq. (2), Definition 1), Ma
(Lemma 3.1), Djokovic (1969), Wanless (2022) and the earlier Evidence Press Bapat-Sunder release were read and checked
against the statements the paper uses. Shchesnovich (2016) and Bapat-Sunder (1986) are cited for attribution only.

## Main objections and dispositions

| Objection | Disposition |
|---|---|
| Correctness | No error found by the panel, the lead's audit or the verification re-review. |
| The sharp constant d - 1 relied on an unpublished companion manuscript | New Proposition 4.3 with a full proof (overlap-conditioned Gaussian average); the companion is no longer cited. |
| An eigenvalue bound was attributed to a theorem of Anari et al. that does not state it | Lemma 2.4 now derives it from their Theorem 1 and eq. (2) via the Schur power of a dominating diagonal matrix. |
| The threshold n_0(d) is not computed, yet the text read as if new finite cases were settled | Abstract and Theorem 1.1 say the threshold is not computed and no finite case is claimed; open problem 1. |
| Prior art sat in a short closing section; no plain introduction | New Section 1: problem, why it matters, what was known (status table), results, method and roadmap; every section opens with its purpose. |
| Rank-three coverage by Pate | Checked against Pate (1999), Theorem 2: every shape with at most three rows is covered, and at rank four every shape except four rows with fourth row at least 3 and second row longer than the third; the ordinary-immanant part is claimed as new only for those shapes. |
| Theorem 1.1(3) claimed a bound for all n over the reals at rank two where only the standard shape is controlled for all n | Restated: the standard supremum lies in [1/(n-1), 3/(n-1)] for all n; the overall supremum lies in that range for large n. |
| Status table overstated novelty at real rank two; Chollet remark scope | Table says "new except real rank 2"; the Chollet remark is scoped to pairs of rank at most d, credits Ma's stronger direct bound, and is not counted as a contribution. |
| Notation overloaded; one forward reference | Integrated bound moved into Lemma 4.1; clashing symbols renamed. |
| Checks could pass vacuously under optimisation | The script never uses `assert`; it exits 1 on any failed check, also under `python -O` (negative control). |

## Limits that remain

The open problems of Section 9 of the paper: an effective threshold n_0(d), the real-rank constants (known up to a
factor), the gap between d = o(n) and d log n = o(n) (first moment) or d log n = o(sqrt n) (second moment) for growing
rank, and the subgroup form of Lieb's conjecture. The novelty search is recorded in RESEARCH_GATES.json.
