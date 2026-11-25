# AI index: Asymptotic permanental dominance at fixed Gram rank (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with a status table and a roadmap). DOI 10.5281/zenodo.24536861 (concept 10.5281/zenodo.24536860).
Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

K = R or C. Gram matrix A = (v_i^* v_j), v_i in K^d; Gram rank = rank A; correlation matrix = unit diagonal (positive
diagonal rescaling leaves both moments unchanged). For lambda |- n: B_lambda(A) = sum_sigma a_sigma rho_lambda(sigma^-1)
(positive semidefinite), Imm_lambda(A) = Tr B_lambda(A), f^lambda = dim rho_lambda. m_1(lambda, A) = Tr B_lambda /
(f^lambda per A); m_2(lambda, A) = Tr B_lambda^2 / (f^lambda (per A)^2). M_j^K(n, d) = sup of m_j over lambda != (n)
and correlation matrices of order n over K of rank at most d. Lieb's conjecture for ordinary immanants: M_1 <= 1. Lavi's
quadratic conjecture (hexagon:2611.00429v1, Conjecture 1.1): M_2 <= 1. Standard shape (n-1, 1). Logarithms natural.

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A (main) | Fixed d >= 2: M_1^K(n,d) -> 0 and M_2^K(n,d) -> 0; shapes with n - lambda_1 >= 2 give o(1/n) uniformly; threshold n_0(d) not computed | Theorem 1.1(1); Sections 2-5 | written proof (compactness, ineffective); new for the quadratic inequality at complex d >= 2 and real d >= 3; for ordinary immanants new only at rank >= 4, and at rank 4 only for four-row shapes with lambda_4 >= 3 and lambda_2 > lambda_3 (all else at rank <= 4: Pate 1999, Theorem 2) |
| B | For large n the supremum is attained at the standard shape | Theorem 1.1(2) | written proof |
| C (main) | Over C: M_1 ~ (d-1) log n / n, M_2 ~ (d-1)(log n)^2 / n. Over R, d >= 3: orders (log n)/n and (log n)^2/n with ratio constants in [1/2, d-1] and [1/4, d-1]. Over R, d = 2: standard sup in [1/(n-1), 3/(n-1)] for all n | Theorem 1.1(3); Sections 4 and 6 | written proof; d - 1 trace constant new (Proposition 4.3); uses the Bapat-Sunder growth of the EP release 10.5281/zenodo.24509057 |
| D | Growing rank: d log n = o(sqrt n) gives M_2 -> 0 (hence M_1 -> 0); d log n = o(n) gives M_1 -> 0 | Theorem 1.2; Section 7 | written proof; new |
| E | Standard shape, every n >= 2, 1 <= d <= n: m_1 <= (6 sqrt(e n d) - 1)/(n-1), m_2 <= (6 sqrt(e n d)(log n + 5) - 1)/(n-1); a real example gives (d-1)/(n-1); Lavi's inequality for the standard shape whenever 6 sqrt(e n d)(log n + 5) <= n (d = 2: from n = 48,836) | Theorem 1.3; Section 8 | written proof; new |
| F | Strict merge inequality for expected logarithms of sums of exponential variables | Lemma 3.1 | written proof; new |
| G | Trace bound for ordered-partition matrices (Gaussian averaging, refined 1/D_{n,d} form) and the sharp first-trace inequality | Lemma 4.1, Corollary 4.2, Proposition 4.3 | written proof |
| H | Eigenvalue bound lambda_max(B_tau) <= e^{(1+gamma)k} per A | Lemma 2.4 | derived from Anari-Gurvits-Oveis Gharan-Saberi, Theorem 1 and eq. (2) |
| I | Eventual Chollet inequality for pairs of rank at most d | Section 9 remark | via Lavi eqs. (2.5), (2.6); weaker than Ma's Lemma 3.1; not counted as a contribution |

## Evidence

| Check | Script | Output | Scope |
|---|---|---|---|
| Part A1: all irreducible blocks for n <= 6, d = 2, 3, real and complex: positivity, real traces, Tr B^2 <= (Tr B)^2, Schur-Weyl vanishing, trivial block, Young trace bound, Lieb splitting, standard-block identities, Plancherel identity, dimension sum | [code/check_dominance.py](code/check_dominance.py) | [results/check_dominance.txt](results/check_dominance.txt) | exact (Gaussian integers, fractions) |
| Part A2: real rank two Tr F_A <= 4 per A (largest found 2.420276); J_a + J_(n-a) gives both standard moments 1/(n-1), 3 <= n <= 8 | same | same | exact |
| Part A3: repeated orthonormal basis gives Tr C_1 = d per A, 1 <= d <= 4, d <= n <= 7 | same | same | exact |
| Part B1-B5: Lemma 3.1 (gap 0.039721 at 0+), Lemma 4.1 two-colour form, Proposition 4.3, Theorem 1.3 chain, Section 6 frame bounds | same | same | floating point, tolerance stated |

The script checks ingredients on finite examples; it is not part of any proof. Numbers quoted on the release page:
(d-1) log n / n and (d-1)(log n)^2 / n at complex rank d; real rank two of order 1/n; d log n = o(sqrt n) and
d log n = o(n); 6 sqrt(e n d); Lieb's version known up to 15 (Pate, Li), Lavi's up to 5 (Lavi). All proved or cited.

## Limits

Threshold n_0(d) not computed: no new finite case. Only irreducible characters of S_n (not the subgroup form of
Lieb's conjecture). Real-rank constants known up to a factor. Growing-rank ranges not sharp. Open problems: Section 9.

## Replay and records

Replay commands and cross-platform policy (Part A exact; Part B decimals within 1e-5 relative or one unit in the last
printed digit): [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml) with [.github/ci/compare_output.py](.github/ci/compare_output.py).
Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json).
Bibliography: [refs.bib](refs.bib). Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes:
[MANIFEST.sha256](MANIFEST.sha256).
