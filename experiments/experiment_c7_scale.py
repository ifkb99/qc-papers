"""TODO 12e -- push C7's CIRCUIT-level density series past 24 qubits.

C7 ("generic r ⟹ Θ(2ⁿ); the small-n results are not pre-asymptotic") rests on
six circuit-level points at q = 15..24, density 0.473 → 0.498, slope 1.008
bits/qubit. `PAPER_A.md` §11.1 concedes the circuit series is modest. It stopped
at 24 because q = 27 was ~40 min and q = 30 hours on CPU; with the GPU backend
(and the int32 replay, which is what makes q = 30 fit in 20 GiB at all) they are
minutes.

Two series, each varying **only N** at fixed n_exp -- the same controlled design
C7 used, extended, not a new one:

  A  n_exp = 2 : the six logged points plus N = 77 (q=27) and N = 143 (q=30)
  B  n_exp = 3 : an independent series at a different fixed n_exp, filling in
                 q = 22, 25, 28. Restricted to beta > 1 instances, because a
                 beta = 1 instance with n_exp > v2(r) has a LOCKED support
                 (C21) and would contribute a point off the trend for a reason
                 that has nothing to do with n. Series A keeps C7's original
                 membership, including its two beta = 1 rows (N = 5 and 15,
                 both r = 4): at n_exp = 2 <= v2(r) they have not locked yet,
                 so they are on the trend legitimately.

Everything is measured on the EXACT integer path (support test `!= 0`, not
|c| > 1e-12), so the re-measured rows also answer a second question: does the
thresholding artifact noted in NOTES §GF touch these numbers?

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P0  Re-measuring the six logged rows exactly reproduces them
      (15493 / 15539 / 127936 / 1037322 / 8347241 / 8346759). Derived: at
      q <= 30 the coefficients are integers/2^q with smallest magnitude
      2^-30 = 9.3e-10, three orders above the 1e-12 threshold, so the
      threshold cannot have been binding. If this FAILS, the logged numbers
      are wrong and everything downstream needs re-checking.

  P1  The new points stay pinned near 1/2: density at q = 27, 28 and 30 lies
      in (0.4975, 0.5). The upper end is DERIVED -- C30 puts the support in
      the hyperplane <z, w> = 0, so density < 1/2 exactly, at every width.
      The lower end is EXTRAPOLATION from the monotone approach seen at
      q <= 24, and is flagged as such.

  P2  The mechanism, not just the number: w = b_msb XOR anc is a linear
      structure at every width including the new ones, i.e. the count of
      support elements with <z, w> = 1 is exactly 0. This is what forces P1's
      upper bound, so it must hold wherever P1 does.

  P3  Fitted slope of log2|support| against q over the extended series is
      within 0.01 of 1.000, i.e. Θ(2^q) over a Hilbert space 64x larger than
      the one C7 settled on.

  C1  MUST FAIL -- density is not automatically ~1/2. A beta = 1 instance
      (N=5, a=4, r=2) grown in n_exp instead of N must show density
      COLLAPSING like 2^-k, because C21 locks its support. If this also sat
      at 1/2 the instrument would be measuring the register size, not the
      circuit.

  C2  MUST FAIL -- the instrument can report Θ(1). A ripple-carry adder's low
      output bit is affine (b0 XOR a0), so its Walsh sparsity is exactly 1 at
      every width and its density is 2^-q, not 1/2.

OUTCOME (2026-08-08): 6/6. All six logged rows reproduce bit-for-bit, so the
thresholding artifact does not reach circuit level. q = 30 (N = 143) gives
|support| = 536,271,623 at density 0.49944, in 851 s. Slopes 1.0054 (A),
1.0083 (B), 1.0062 pooled, against C7's logged 1.008 over q <= 24. The C30
hyperplane holds in 13/13 circuits. Both controls fail as required.

One thing this could NOT do: q = 26, 28 are unreachable in series A, because
q = 3n + 6 there, so the modulus alone only reaches q = 24, 27, 30. Series B
exists to fill in q = 22, 25, 28 -- that is why there are two series and not
one longer one.

Run:  LAB_GPU=1 uv run python -m experiments.experiment_c7_scale
"""
from __future__ import annotations

import math
import time

import numpy as np

from lab import Experiment, build_modexp, order, v2_split, stats, support
from circuits import ripple_adder

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P0", "the six logged C7 rows reproduce exactly on the exact path")
exp.predict("P1", "density at q = 27, 28, 30 lies in (0.4975, 0.5)")
exp.predict("P2", "w = b_msb XOR anc is a linear structure at every width "
                  "(odd-count exactly 0)")
exp.predict("P3", "slope of log2|S| vs q is within 0.01 of 1.000")
exp.must_fail("C1", "beta=1 grown in n_exp: density must COLLAPSE, not stay ~1/2")
exp.must_fail("C2", "affine observable: sparsity exactly 1, density 2^-q")

LOGGED = {(5, 2): 15493, (7, 3): 15539, (15, 7): 127936,
          (21, 2): 1037322, (33, 5): 8347241, (35, 3): 8346759}


def row(N, a, n_exp):
    me = build_modexp(N=N, a=a, n_exp=n_exp)
    qc = me.build()
    w = (1 << me.b[me.m - 1]) | (1 << me.anc)      # the C30 linear structure
    t0 = time.time()
    st = stats(qc, me.x[0], masks=(w,))
    r = order(a, N)
    al, beta = v2_split(r)
    print(f"    {N:5d} {a:3d} {me.n:3d} {r:4d} {beta:5d} {qc.n:4d} "
          f"{st['count']:12,d} {st['density']:9.5f} {st['odd'][w]:9d} "
          f"{time.time() - t0:8.1f}", flush=True)
    return qc.n, st, w


def header():
    print(f"    {'N':>5} {'a':>3} {'n':>3} {'r':>4} {'beta':>5} {'q':>4} "
          f"{'|support|':>12} {'density':>9} {'<z,w>=1':>9} {'secs':>8}")


SERIES_A = [(5, 2), (7, 3), (15, 7), (21, 2), (33, 5), (35, 3),
            (77, 3), (143, 5)]
SERIES_B = [(7, 3), (21, 2), (33, 5), (35, 3), (77, 3)]

# ---------------------------------------------------------------------------
exp.section("Series A -- n_exp = 2, vary N (C7's own design, extended)")
header()
A = {}
for N, a in SERIES_A:
    assert math.gcd(a, N) == 1
    A[(N, a)] = row(N, a, 2)

exp.section("Series B -- n_exp = 3, vary N (independent, beta > 1 only)")
header()
B = {}
for N, a in SERIES_B:
    r = order(a, N)
    assert v2_split(r)[1] > 1, f"N={N} a={a} is beta=1; excluded by design"
    B[(N, a)] = row(N, a, 3)

# ---------------------------------------------------------------------------
exp.section("P0  do the six logged rows reproduce exactly?")
mism = [(k, A[k][1]["count"], v) for k, v in LOGGED.items()
        if A[k][1]["count"] != v]
for k, got, want in mism:
    exp.log(f"  MISMATCH {k}: exact {got:,} vs logged {want:,}")
exp.check("P0", not mism,
          f"{len(LOGGED) - len(mism)}/{len(LOGGED)} logged rows reproduced "
          f"bit-for-bit on the exact path")

exp.section("P1  density at the new widths")
new = [("A", (77, 3), A[(77, 3)]), ("A", (143, 5), A[(143, 5)]),
       ("B", (77, 3), B[(77, 3)])]
for tag, k, (q, st, _) in new:
    exp.log(f"series {tag} {k}: q={q} density={st['density']:.5f}")
exp.check("P1", all(0.4975 < st["density"] < 0.5 for _, _, (_, st, _) in new),
          "all three new points in (0.4975, 0.5); the < 1/2 half is C30, "
          "the > 0.4975 half was extrapolation")

exp.section("P2  the C30 hyperplane at every width")
bad = [(k, st["odd"][w]) for k, (q, st, w) in list(A.items()) + list(B.items())
       if st["odd"][w] != 0]
exp.log(f"widths checked: {len(A) + len(B)}; "
        f"nonzero <z,w>=1 counts: {bad if bad else 'none'}")
exp.check("P2", not bad,
          f"w = b_msb XOR anc is a linear structure in {len(A) + len(B)}/"
          f"{len(A) + len(B)} circuits, q = 15..30")

exp.section("P3  growth exponent over the extended series")
for tag, S in (("A", A), ("B", B)):
    qs = np.array([v[0] for v in S.values()], float)
    sp = np.array([v[1]["count"] for v in S.values()], float)
    slope = float(np.polyfit(qs, np.log2(sp), 1)[0])
    exp.log(f"series {tag}: q = {sorted(int(x) for x in qs)}, "
            f"slope {slope:.4f} bits/qubit")
    if tag == "A":
        slope_a = slope
    else:
        slope_b = slope
allq = np.array([v[0] for v in list(A.values()) + list(B.values())], float)
allsp = np.array([v[1]["count"] for v in list(A.values()) + list(B.values())],
                 float)
slope_all = float(np.polyfit(allq, np.log2(allsp), 1)[0])
exp.log(f"pooled: slope {slope_all:.4f} bits/qubit "
        f"(C7 logged 1.008 over q <= 24)")
exp.check("P3", all(abs(s - 1.0) < 0.01 for s in (slope_a, slope_b, slope_all)),
          f"A {slope_a:.4f}, B {slope_b:.4f}, pooled {slope_all:.4f}")

# ---------------------------------------------------------------------------
exp.section("C1  control -- beta = 1 grown in n_exp must COLLAPSE in density")
header()
dens_c1 = []
for ne in (2, 4, 6):
    q, st, w = row(5, 4, ne)
    dens_c1.append(st["density"])
ratios = [dens_c1[i + 1] / dens_c1[i] for i in range(len(dens_c1) - 1)]
exp.log(f"densities {[f'{d:.5f}' for d in dens_c1]}, "
        f"ratios {[f'{r:.4f}' for r in ratios]} (2^-2 = 0.25 expected)")
exp.fail_check("C1", all(d < 0.4975 for d in dens_c1[1:]),
               f"density falls to {dens_c1[-1]:.5f}: the ~1/2 reading is a "
               f"property of the circuit, not of the register size")

exp.section("C2  control -- an affine observable must read Θ(1), not Θ(2^q)")
sp_c2 = []
for nb in (3, 5, 7, 9):
    qc, lay = ripple_adder(nb)
    zs = support(qc, lay["b"][0], exact=True)
    sp_c2.append(zs.size)
    exp.log(f"{nb}-bit adder, q={qc.n}: |support|={zs.size}, "
            f"density={zs.size / (1 << qc.n):.2e}")
exp.fail_check("C2", all(s == 1 for s in sp_c2),
               f"sparsity {sp_c2} -- constant 1 while q grows, so the "
               f"instrument does report Θ(1) when Θ(1) is the truth")

exp.finish()
