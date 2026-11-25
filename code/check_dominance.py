"""Finite checks for "Asymptotic permanental dominance at fixed Gram rank".

The theorems in the paper are proved analytically. This script checks, on
explicit finite examples, the identities and inequalities that the proofs
use. Part A is exact (integer / Gaussian-integer / rational arithmetic).
Part B uses floating point with a declared tolerance and prints only
rounded values or PASS/FAIL, so its output is the same on every platform.

Run:  python3 check_dominance.py   (standard library only)
The script exits with status 1 if any check fails; it does not use assert.
"""

import cmath
import math
import random
import sys
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial

FAILURES = []


def check(cond, label):
    if not cond:
        FAILURES.append(label)
        print("FAIL:", label)


# ---------------------------------------------------------------- Gaussian integers
def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def gconj(x):
    return (x[0], -x[1])


ZERO, ONE = (0, 0), (1, 0)


def gram(vectors):
    """a_ij = v_i^* v_j for Gaussian-integer vectors."""
    n = len(vectors)
    A = [[ZERO] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = ZERO
            for p, q in zip(vectors[i], vectors[j]):
                s = gadd(s, gmul(gconj(p), q))
            A[i][j] = s
    return A


def per(A, rows, cols):
    """Permanent of the submatrix (Gaussian integers), Laplace expansion with memo."""
    rows, cols = tuple(rows), tuple(cols)
    memo = {}

    def rec(r, used):
        if r == len(rows):
            return ONE
        key = (r, used)
        if key in memo:
            return memo[key]
        s = ZERO
        for c in range(len(cols)):
            if not used >> c & 1:
                s = gadd(s, gmul(A[rows[r]][cols[c]], rec(r + 1, used | 1 << c)))
        memo[key] = s
        return s

    return rec(0, 0)


# ---------------------------------------------------------------- characters
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def mn_character(lam, mu):
    """Murnaghan-Nakayama rule, beta-set (abacus) form."""
    if not mu:
        return 1 if sum(lam) == 0 else 0
    r, rest = mu[0], mu[1:]
    L = len(lam)
    beta = [lam[i] + (L - 1 - i) for i in range(L)]
    bset = set(beta)
    total = 0
    for b in beta:
        if b - r >= 0 and (b - r) not in bset:
            height = sum(1 for c in beta if b - r < c < b)
            newbeta = sorted([c for c in beta if c != b] + [b - r], reverse=True)
            Ln = len(newbeta)
            newlam = tuple(x for x in (newbeta[i] - (Ln - 1 - i) for i in range(Ln)) if x > 0)
            total += (-1) ** height * mn_character(newlam, rest)
    return total


def cycle_type(p):
    seen, out = [False] * len(p), []
    for i in range(len(p)):
        if not seen[i]:
            j, length = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                length += 1
            out.append(length)
    return tuple(sorted(out, reverse=True))


def compose(p, q):  # (p q)(i) = p(q(i))
    return tuple(p[q[i]] for i in range(len(q)))


def ordered_partition_trace(A, sizes):
    """Tr C_alpha = sum over ordered set partitions with the given sizes of prod per A[I_l]."""
    n = len(A)

    def rec(remaining, k):
        if k == len(sizes):
            return ONE
        s = ZERO
        for I in combinations(remaining, sizes[k]):
            rest = tuple(x for x in remaining if x not in I)
            s = gadd(s, gmul(per(A, I, I), rec(rest, k + 1)))
        return s

    return rec(tuple(range(n)), 0)


def exact_block_checks():
    print("Part A1. All irreducible blocks of small Gram frames (exact).")
    rng = random.Random(20261124)
    worst = {}
    for n in range(2, 7):
        for d in (2, 3):
            for field in ("real", "complex"):
                vecs = []
                for _ in range(n):
                    if field == "real":
                        v = [(rng.randint(-3, 3), 0) for _ in range(d)]
                    else:
                        v = [(rng.randint(-3, 3), rng.randint(-3, 3)) for _ in range(d)]
                    if all(x == ZERO for x in v):
                        v[0] = (1, 0)
                    vecs.append(v)
                A = gram(vecs)
                perms = list(permutations(range(n)))
                a = {}
                for s in perms:
                    t = ONE
                    for i in range(n):
                        t = gmul(t, A[i][s[i]])
                    a[s] = t
                P = per(A, range(n), range(n))
                # convolution c(g) = sum_{s t = g} a_s a_t
                conv = {}
                for s in perms:
                    for t in perms:
                        g = compose(s, t)
                        conv[g] = gadd(conv.get(g, ZERO), gmul(a[s], a[t]))
                ct = {s: cycle_type(s) for s in perms}
                plancherel, dimsum_imm = 0, ZERO
                for lam in partitions(n):
                    chi = {mu: mn_character(lam, mu) for mu in set(ct.values())}
                    f = chi[(1,) * n]
                    imm, q2 = ZERO, ZERO
                    for s in perms:
                        c = chi[ct[s]]
                        imm = gadd(imm, (c * a[s][0], c * a[s][1]))
                        q2 = gadd(q2, (c * conv[s][0], c * conv[s][1]))
                    check(imm[1] == 0 and q2[1] == 0, f"real traces n={n} {lam}")
                    imm, q2 = imm[0], q2[0]
                    check(imm >= 0 and q2 >= 0, f"PSD traces n={n} {lam}")
                    check(q2 <= imm * imm, f"Tr B^2 <= (Tr B)^2 n={n} {lam}")
                    if len(lam) > d:
                        check(imm == 0 and q2 == 0, f"Schur-Weyl vanishing n={n} {lam} d={d}")
                    if lam == (n,):
                        check(imm == P[0] and q2 == P[0] ** 2, f"trivial block n={n}")
                    trC = ordered_partition_trace(A, lam)
                    check(trC[1] == 0 and imm <= trC[0], f"Young trace bound n={n} {lam}")
                    k = n - lam[0]
                    if k >= 1 and len(lam) >= 2:
                        trCk = ordered_partition_trace(A, (n - k, k))[0]
                        check(trC[0] <= (len(lam) - 1) ** k * trCk, f"Lieb splitting n={n} {lam}")
                    if lam == (n - 1, 1):
                        F = [[gmul(A[i][j], per(A, [r for r in range(n) if r != i], [c for c in range(n) if c != j]))
                              for j in range(n)] for i in range(n)]
                        trF = sum(F[i][i][0] for i in range(n))
                        trF2 = sum(gmul(F[i][j], F[j][i])[0] for i in range(n) for j in range(n))
                        check(imm == trF - P[0], f"standard trace identity n={n}")
                        check(q2 == trF2 - P[0] ** 2, f"standard square identity n={n}")
                    if P[0] > 0 and lam != (n,):
                        m1 = Fraction(imm, f * P[0])
                        m2 = Fraction(q2, f * P[0] ** 2)
                        key = (n, field)
                        worst[key] = max(worst.get(key, (0, 0)), (m2, m1))
                    plancherel += f * q2
                    dimsum_imm = gadd(dimsum_imm, (f * imm, 0))
                perAA = sum(a[s][0] ** 2 + a[s][1] ** 2 for s in perms)
                check(plancherel == factorial(n) * perAA, f"Plancherel n={n}")
                diag = ONE
                for i in range(n):
                    diag = gmul(diag, A[i][i])
                check(dimsum_imm[0] == factorial(n) * diag[0], f"dimension sum n={n}")
    for (n, field), (m2, m1) in sorted(worst.items()):
        print(f"  n={n} {field:7s} largest nontrivial m2 = {float(m2):.6f}, with m1 = {float(m1):.6f}")
    print("  checked: real traces, positivity, Tr B^2 <= (Tr B)^2, Schur-Weyl vanishing,")
    print("  trivial block, Young trace bound, Lieb splitting, standard-block identities,")
    print("  Plancherel identity n! per(A o conj A) = sum f Tr B^2, dimension sum.")


def ryser(M):
    n = len(M)
    if n == 0:
        return Fraction(1)
    total = Fraction(0)
    for mask in range(1, 1 << n):
        cols = [j for j in range(n) if mask >> j & 1]
        prod = Fraction(1)
        for i in range(n):
            prod *= sum(M[i][j] for j in cols)
        total += (-1) ** len(cols) * prod
    return (-1) ** n * total


def real_rank_two_trace():
    print("Part A2. Real rank two: Tr F_A <= 4 per A (exact rationals).")
    rng = random.Random(7)
    worst = Fraction(0)
    for n in range(2, 10):
        for _ in range(6):
            vs = [(rng.randint(-5, 5), rng.randint(-5, 5)) for _ in range(n)]
            vs = [v if v != (0, 0) else (1, 0) for v in vs]
            A = [[Fraction(p[0] * q[0] + p[1] * q[1]) for q in vs] for p in vs]
            P = ryser(A)
            trF = sum(A[i][i] * ryser([[A[r][c] for c in range(n) if c != i] for r in range(n) if r != i])
                      for i in range(n))
            ratio = trF / P
            check(ratio <= 4, f"real rank-two trace n={n}")
            worst = max(worst, ratio)
    print(f"  largest Tr F_A / per A found: {float(worst):.6f} (bound 4)")
    for j in range(1, 30):
        check(Fraction(j) / (Fraction(j) - Fraction(1, 2)) <= 2, "gamma ratio")
    # J_a + J_{n-a}: standard first and second moments equal 1/(n-1)
    for n in range(3, 9):
        for a in range(1, n):
            vs = [(1, 0)] * a + [(0, 1)] * (n - a)
            A = [[Fraction(p[0] * q[0] + p[1] * q[1]) for q in vs] for p in vs]
            P = ryser(A)
            F = [[A[i][j] * ryser([[A[r][c] for c in range(n) if c != j] for r in range(n) if r != i])
                  for j in range(n)] for i in range(n)]
            trF = sum(F[i][i] for i in range(n))
            trF2 = sum(F[i][j] * F[j][i] for i in range(n) for j in range(n))
            check((trF - P) / ((n - 1) * P) == Fraction(1, n - 1), "J+J first moment")
            check((trF2 - P * P) / ((n - 1) * P * P) == Fraction(1, n - 1), "J+J second moment")
    print("  J_a (+) J_(n-a): both standard moments equal 1/(n-1) for 3 <= n <= 8 (exact).")


def threshold_lower():
    print("Part A3. Repeated orthonormal basis: Tr C_1 / per A = d (exact).")
    for d in range(1, 5):
        for n in range(d, 8):
            m = [n // d + (1 if j < n % d else 0) for j in range(d)]
            vs = []
            for j in range(d):
                vs += [tuple(1 if t == j else 0 for t in range(d))] * m[j]
            A = [[Fraction(sum(p[t] * q[t] for t in range(d))) for q in vs] for p in vs]
            P = ryser(A)
            trF = sum(ryser([[A[r][c] for c in range(n) if c != i] for r in range(n) if r != i]) for i in range(n))
            check(trF / P == d, f"repeated basis d={d} n={n}")
    print("  equality Tr C_1 = d per A holds for 1 <= d <= 4, d <= n <= 7.")


# ---------------------------------------------------------------- Part B (floating point)
TOL = 1e-9


def cper(M):
    n = len(M)
    if n == 0:
        return 1.0
    total = 0j
    for mask in range(1, 1 << n):
        cols = [j for j in range(n) if mask >> j & 1]
        prod = 1 + 0j
        for i in range(n):
            prod *= sum(M[i][j] for j in cols)
        total += (-1) ** len(cols) * prod
    return ((-1) ** n * total).real


def random_frame(rng, n, d):
    vs = []
    for _ in range(n):
        v = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(d)]
        nv = math.sqrt(sum(abs(x) ** 2 for x in v))
        vs.append([x / nv for x in v])
    return [[sum(p[t].conjugate() * q[t] for t in range(d)) for q in vs] for p in vs]


def F2(t):
    return 1.0 if abs(t - 1) < 1e-12 else t * math.log(1 / t) / (1 - t)


def D(s, d):
    return comb(s + d - 1, d - 1)


def scalar_lemma():
    print("Part B1. Strict merge inequality (tolerance 1e-9).")
    ok = True
    for i in range(1, 4001):
        x = i * 0.01
        lhs = x / math.tanh(x)
        rhs = 1.5 * math.log(2 * math.cosh(2 * x / 3))
        z = 2 * x / 3
        ok &= rhs - lhs > -TOL and 2 * math.sinh(z) + math.tanh(z) - 3 * z > 0
        # derivative identity for g(z) = 2 sinh z + tanh z - 3z
        g1 = 2 * math.cosh(z) + 1 / math.cosh(z) ** 2 - 3
        g2 = (math.cosh(z) - 1) ** 2 * (2 * math.cosh(z) + 1) / math.cosh(z) ** 2
        ok &= abs(g1 - g2) <= 1e-7 * max(1.0, g1)
        if x <= 10:
            # 2 sinh^2 x (coth x - x/sinh^2 x - tanh z) = 2 sinh z + tanh z - 3z
            lhs2 = 2 * math.sinh(x) ** 2 * (1 / math.tanh(x) - x / math.sinh(x) ** 2 - math.tanh(z))
            rhs2 = 2 * math.sinh(z) + math.tanh(z) - 3 * z
            ok &= abs(lhs2 - rhs2) <= 1e-8 * max(1.0, math.sinh(x) ** 2)
    check(ok, "merge reduction on (0, 40]")
    print(f"  x coth x < (3/2) log(2 cosh(2x/3)) on (0,40]: {'PASS' if ok else 'FAIL'};"
          f" gap at 0+ = {1.5 * math.log(2) - 1:.6f}")
    # Laplace-integral sign for c >= 0
    rng = random.Random(3)
    worst = -1.0
    for _ in range(40):
        a, b, c = rng.uniform(0.05, 3), rng.uniform(0.05, 3), rng.choice([0.0, rng.uniform(0, 3)])
        C = (a ** (2 / 3) + b ** (2 / 3)) ** 1.5
        # integral of e^{-cu}(a+b-C+abu)/((1+Cu)(1+au)(1+bu)) on (0, inf), substitution u = s/(1-s)
        N, tot = 20000, 0.0
        for k in range(N):
            s = (k + 0.5) / N
            u = s / (1 - s)
            tot += math.exp(-c * u) * (a + b - C + a * b * u) / ((1 + C * u) * (1 + a * u) * (1 + b * u)) / (1 - s) ** 2
        worst = max(worst, tot / N)
    check(worst < 0, "merge integral negative")
    print(f"  E log(c+aX+bY) - E log(c+CZ) < 0 at 40 random (a,b,c): {'PASS' if worst < 0 else 'FAIL'}")


def coloured_trace_bound():
    print("Part B2. Two-colour trace bound (integrated form), random complex frames.")
    rng = random.Random(11)
    ok, worst = True, 0.0
    for n in range(3, 8):
        for d in (2, 3):
            A = random_frame(rng, n, d)
            P = cper(A)
            for k in range(1, n):
                trCk = 0.0
                for I in combinations(range(n), k):
                    J = [j for j in range(n) if j not in I]
                    trCk += cper([[A[r][c] for c in I] for r in I]) * cper([[A[r][c] for c in J] for r in J])
                best = min(D(k, d) * D(n - k, d) / (D(n, d) * comb(n, k)) * t ** (-k) * math.exp(n * F2(t))
                           for t in [10 ** (e / 20) for e in range(-80, 41)])
                r = trCk / P / best
                ok &= r <= 1 + TOL
                worst = max(worst, r)
    check(ok, "two-colour trace bound")
    print(f"  Tr C_k / per A <= bound for 3 <= n <= 7, d = 2,3, all k: {'PASS' if ok else 'FAIL'}"
          f" (largest ratio to bound {worst:.3f})")


def first_trace_constant():
    print("Part B3. Sharp first-trace inequality Tr C_1/per A <= d/(n+d-1) (E e^{nL}-1)/t.")
    rng = random.Random(5)
    ok = True
    for n in range(3, 8):
        for d in (2, 3):
            A = random_frame(rng, n, d)
            P = cper(A)
            tr1 = sum(cper([[A[r][c] for c in range(n) if c != i] for r in range(n) if r != i]) for i in range(n))

            def bound(t):
                # alpha = |z^* y|^2 ~ Beta(1, d-1); midpoint rule in alpha
                N, tot = 4000, 0.0
                for k in range(N):
                    al = (k + 0.5) / N
                    w = (d - 1) * (1 - al) ** (d - 2)
                    disc = math.sqrt((1 - t) ** 2 + 4 * t * al)
                    lp, lm = (1 + t + disc) / 2, (1 + t - disc) / 2
                    L = math.log(lp) + (F2(lm / lp) if lm > 0 else 0.0)
                    tot += w * math.exp(n * L)
                return d / (n + d - 1) * (tot / N - 1) / t

            b = min(bound(10 ** (e / 4)) for e in range(-12, 5))
            ok &= tr1 / P <= b * (1 + 1e-6)
    check(ok, "first-trace inequality")
    print(f"  holds on random complex frames, 3 <= n <= 7, d = 2,3: {'PASS' if ok else 'FAIL'}")


def threshold_upper():
    print("Part B4. Tr C_1/per A <= 3 log(per(A+I)/per A) <= 6 sqrt(e n d).")
    rng = random.Random(9)
    ok = True
    for n in range(2, 9):
        for d in (1, 2, 3, 4):
            A = random_frame(rng, n, d)
            P = cper(A)
            tr1 = sum(cper([[A[r][c] for c in range(n) if c != i] for r in range(n) if r != i]) for i in range(n))
            AI = [[A[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
            lp = math.log(cper(AI) / P)
            ok &= tr1 / P <= 3 * lp * (1 + TOL) and 3 * lp <= 6 * math.sqrt(math.e * n * d)
    check(ok, "standard threshold upper chain")
    print(f"  holds on random complex frames, 2 <= n <= 8, 1 <= d <= 4: {'PASS' if ok else 'FAIL'}")


# ---------------------------------------------------------------- Part B5: the shell in h directions
def polymul(p, q):
    out = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = tuple(x + y for x, y in zip(m1, m2))
            out[m] = out.get(m, 0) + c1 * c2
    return out


def fock_norm2(p):
    tot = 0.0
    for m, c in p.items():
        w = 1
        for e in m:
            w *= factorial(e)
        tot += abs(c) ** 2 * w
    return tot


def shell_directions():
    print("Part B5. Concentric-circle frame in h = d-1 directions: Rayleigh bound on the test space.")
    ok = True
    for h, ks, a in ((1, (3, 6), 0.75), (2, (3,), 0.75), (2, (3, 4), 0.8)):
        d = h + 1
        cs = [math.sqrt(a * k) * cmath.exp(2j * math.pi * j / k) for k in ks for j in range(k)]
        S, L = len(cs), len(ks)
        M = 3 * h * S
        n = M + h * S
        # linear forms as coefficient vectors in (u, w_1..w_h)
        forms = [[1.0 + 0j] + [0j] * h for _ in range(M)]
        group = []
        for b in range(h):
            for c in cs:
                f = [c / math.sqrt(M)] + [0j] * h
                f[1 + b] = 1.0 + 0j
                forms.append(f)
                group.append((b, c))

        def lin(f):
            return {tuple(1 if t == p else 0 for t in range(d)): f[p] for p in range(d) if f[p] != 0}

        prefix = [{(0,) * d: 1.0 + 0j}]
        for f in forms:
            prefix.append(polymul(prefix[-1], lin(f)))
        suffix = [{(0,) * d: 1.0 + 0j}]
        for f in reversed(forms):
            suffix.append(polymul(suffix[-1], lin(f)))
        suffix = suffix[::-1]
        perA = fock_norm2(prefix[-1])
        # D and ||p||^{2h}
        prodc = 1
        for c in cs:
            prodc *= c ** h
        KM = factorial(n) / M ** (h * S)
        Dval = perA / (abs(prodc) ** 2 * KM)
        pnorm2 = fock_norm2({(e,): v for e, v in enumerate(_poly_p(cs))})
        ok &= Dval <= pnorm2 ** h * (1 + 1e-9)
        # quadratic form x^* F x = sum_p ||Q_p||^2 on the test space x_{bj} = l_b / conj(c_j)
        for trial in range(3):
            l = [complex(math.cos(trial + b), math.sin(2 * trial + b)) for b in range(h)]
            x = [0j] * M + [l[b] / c.conjugate() for (b, c) in group]
            num = 0.0
            for p in range(d):
                Q = {}
                for i, f in enumerate(forms):
                    if x[i] == 0 or f[p] == 0:
                        continue
                    term = polymul(prefix[i], suffix[i + 1])
                    for m, v in term.items():
                        Q[m] = Q.get(m, 0) + x[i] * f[p] * v
                num += fock_norm2(Q)
            xx = sum(abs(v) ** 2 for v in x)
            rayleigh = num / (perA * xx)
            predicted = (M / (n * Dval)) * (L / a)
            ok &= rayleigh >= predicted * (1 - 1e-9)
            orth = sum(x)  # tests lie in the complement of the all-ones vector
            ok &= abs(orth) < 1e-9
        print(f"  h={h} circles={ks} a={a} n={n}: D={Dval:.6f} <= ||p||^(2h)={pnorm2 ** h:.6f}; "
              f"Rayleigh >= (M/(nD))(L/a)={predicted:.6f}")
    check(ok, "shell in h directions")
    print(f"  {'PASS' if ok else 'FAIL'}")


def _poly_p(cs):
    coeffs = [1.0 + 0j]
    for c in cs:
        new = [0j] * (len(coeffs) + 1)
        for e, v in enumerate(coeffs):
            new[e] += v
            new[e + 1] += v / c
        coeffs = new
    return coeffs


if __name__ == "__main__":
    exact_block_checks()
    real_rank_two_trace()
    threshold_lower()
    scalar_lemma()
    coloured_trace_bound()
    first_trace_constant()
    threshold_upper()
    shell_directions()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED")
        sys.exit(1)
    print("All checks passed.")
