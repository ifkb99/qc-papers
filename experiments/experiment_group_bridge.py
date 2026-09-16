"""Is C99's differential count special to XOR, and does a group-matched basis shrink the actual Shor observable?

STATUS (after run 1, 2026-09-14): exits nonzero ON PURPOSE. PA4 is REFUTED:
the exhaustive maxima are 16 / 55 / 48, not 16 / 56 / 56. The telescoping
obstruction and scope are in claims/C100.md and notes/GB-group-bridge.md.
All other checks pass. Do not re-scope PA4 to make this file pass.

Context. C99 counts Pauli (GF(2)^n Weyl) terms of U^dag X^a Z^b U through the
XOR difference table. Nothing in the identity needs XOR: for any finite
abelian group G with Weyl operators W(a,chi)|y> = chi(y)|y+a>,

    U^dag W(a,chi) U |y> = chi(pi(y)) |sigma(y)>,  sigma = pi^-1 T_a pi,
    coefficient on W(c,psi) = |G|^-1 sum_{y in D_c} chi(pi(y)) conj(psi(y)),
    D_c = {y : sigma(y) - y = c},  |D_c| = DDT^G_{pi^-1}(a, c),         (I_G)

so T_G(a,chi) = sum_c |supp F_G[1_{D_c} chi o pi]|. Uncertainty on G and
Cauchy-Schwarz give sum_c |G|/|D_c| <= T_G <= R_a |G|. The XOR-only
ingredient in C99 is that sigma is an involution (every element of GF(2)^n
has order 2): classes are unions of c-pairs, so R_a <= |G|/2 and T <= |G|^2/4.
In a group with elements of order > 2 the pairing disappears. For
T_(a) of order o, a 2-point class {y1,y2} has Fourier support |G| or
|G| - |G|/o(y2-y1), and singleton classes have full support |G|.

The post-inverse-QFT measurement in order finding is a combination of
Z/2^t clock shifts on the exponent register, so G = Z/2^t x GF(2)^rest is
the group matched to the real observable (TODO14's original question).

PREDICTIONS, WRITTEN BEFORE MEASURING.
(The helper `lab.differential.group_terms` was smoke-tested against its dense
reference on one random 8-element permutation before this header was
written; no prediction below was measured.)

Part A -- groups of order 8 and 16.
  PA1  (I_G) equals explicit Weyl-matrix traces for every one of the |G|^2
       labels, for G in {GF(2)^4, Z/16, Z/4xZ/4, Z/2xZ/8} and permutations
       {random (seed 20260914), +5 mod 16, a GF(2)-linear map, GF(16) inversion}.
  PA2  Group-matched collapse: +5 mod 16 has T = 1 for every label in Z/16
       but T > 1 for some label in GF(2)^4; the GF(2)-linear map has T = 1
       for every label in GF(2)^4 but T > 1 for some label in Z/16.
  PA3  Bounds sum_c |G|/|D_c| <= T <= R_a |G| hold for every a != 0 label,
       every group and permutation of PA1.
  PA4  Exhaustive over all 40320 permutations of each order-8 group, the
       maximum of T over a != 0 and all characters is 16 for GF(2)^3 and 56
       for both Z/8 and Z/2xZ/4 (R_a = 7 with full-support double class).

  CA1  The XOR pairing cap R_a <= |G|/2 must be violated in Z/8 by some
       permutation and direction (exhaustive search).
  CA2  The XOR pair rule "every 2-point class has support exactly |G|/2"
       must fail for some class in Z/16 (PA1 fixtures).

Part B -- the actual observable, ToffoliModExp(N=7, a=3, n_exp=3), 16 qubits,
full dirty space, exponent = top three qubits as one Z/8 digit. For each
shift a in 1..7, M_a = U^dag (T_a (x) I) U is ONE operator; count its terms
in the Pauli basis (radix 2^16) and in the mixed basis (radix 2^13 x Z/8).
  PB1  Every mixed-basis X-part c has exponent digit exactly a (the exponent
       register is only a control and is restored).
  PB2  For a = 4 the clock shift equals X on the top exponent qubit, so its
       Pauli count equals C99's formula count for that Pauli label.
  PB3  OPEN, no commitment. H-same: T_pauli ~= T_mixed for every a.
       H-mixed: T_mixed < T_pauli for odd a (carries in e+a cost Pauli terms).
       H-pauli: T_pauli < T_mixed. Report every a.

Run:  uv run python -m experiments.experiment_group_bridge   (from research/)
CPU; Part B parallelised over a with a process pool.
"""
from __future__ import annotations

import itertools
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from lab import Experiment
from lab.differential import (dense_group_terms, group_terms, group_classes,
                              g_add, g_sub, g_char, power_map, formula_terms,
                              permutation_matrix_terms)
from toffoli_arith import ToffoliModExp
import walsh

GROUPS16 = {"GF2^4": (2, 2, 2, 2), "Z16": (16,), "Z4xZ4": (4, 4), "Z2xZ8": (2, 8)}
GROUPS8 = {"GF2^3": (2, 2, 2), "Z8": (8,), "Z2xZ4": (2, 4)}


def gf2_linear16(rng):
    while True:
        A = rng.integers(0, 2, size=(4, 4))
        if round(abs(np.linalg.det(A))) % 2 == 1:
            break
    ys = np.arange(16)
    bits = (ys[:, None] >> np.arange(4)) & 1
    return (((bits @ A.T) % 2) * (1 << np.arange(4))).sum(axis=1).astype(np.int64)


# ---------------------------------------------------------------- exhaustive order 8
def all_perms8():
    return np.array(list(itertools.permutations(range(8))), dtype=np.int64)


def exhaustive_T(P, radix, tol=1e-9):
    """T(pi, a, k) for all permutations rows of P, all a != 0, all k; plus R_a."""
    n = P.shape[1]
    ys = np.arange(n)
    inv = np.argsort(P, axis=1)
    rows = np.arange(len(P))[:, None]
    shape = (len(P),) + tuple(reversed(radix))
    axes = tuple(range(1, len(shape)))
    T = np.zeros((len(P), n, n), dtype=np.int64)
    R = np.zeros((len(P), n), dtype=np.int64)
    for a in range(1, n):
        sigma = inv[rows, g_add(P, a, radix)]
        cvec = g_sub(sigma, ys[None, :], radix)
        R[:, a] = np.array([len(np.unique(r)) for r in cvec])
        for k in range(n):
            s = g_char(k, P, radix)
            tot = np.zeros(len(P), dtype=np.int64)
            for cval in range(1, n):
                g = np.where(cvec == cval, s, 0).reshape(shape)
                F = np.fft.fftn(g, axes=axes)
                tot += (np.abs(F).reshape(len(P), -1) > tol).sum(axis=1)
            T[:, a, k] = tot
    return T, R


# ---------------------------------------------------------------- Part B worker
def part_b_worker(args):
    perm, a, low, t = args
    D = len(perm)
    ys = np.arange(D, dtype=np.int64)
    inv = np.argsort(perm)
    radix_mix = (2,) * low + (1 << t,)
    shift = a << low                                   # exponent digit = a
    sigma = inv[g_add(perm, shift, radix_mix)]
    cvec = g_sub(sigma, ys, radix_mix)
    exp_digits = np.unique(cvec >> low)
    t0 = time.time()
    T_mix, R_mix = permutation_matrix_terms(sigma, radix_mix)
    T_pau, R_pau = permutation_matrix_terms(sigma, (2,) * (low + t))
    return a, exp_digits.tolist(), T_mix, R_mix, T_pau, R_pau, time.time() - t0


if __name__ == "__main__":
    exp = Experiment("experiment_group_bridge", doc=__doc__)
    exp.predict("PA1", "(I_G) == dense Weyl traces, all labels, 4 groups x 4 permutations")
    exp.predict("PA2", "+5 mod 16: T=1 in Z16, >1 somewhere in GF2^4; GF(2)-linear: reverse")
    exp.predict("PA3", "sum_c |G|/|D_c| <= T <= R_a|G| everywhere")
    exp.predict("PA4", "exhaustive order-8 maxima: GF2^3 16, Z8 56, Z2xZ4 56")
    exp.must_fail("CA1", "XOR pairing cap R_a <= |G|/2 violated in Z8")
    exp.must_fail("CA2", "2-point classes all |G|/2 fails in Z16")
    exp.predict("PB1", "every mixed X-part has exponent digit a")
    exp.predict("PB2", "a=4: Pauli count == C99 formula count for X on top exponent qubit")
    exp.predict("PB3", "OPEN: report T_pauli vs T_mixed for a=1..7")

    rng = np.random.default_rng(20260914)
    perms16 = {"rand16": rng.permutation(16).astype(np.int64),
               "plus5": ((np.arange(16) + 5) % 16).astype(np.int64),
               "gf2lin": gf2_linear16(rng),
               "inv16": power_map(4, 14)}

    exp.section("PA1/PA2/PA3/CA2  order-16 groups")
    bad, bnd_ok, pair_vals = 0, True, set()
    Tall = {}
    for gname, radix in GROUPS16.items():
        for pname, perm in perms16.items():
            dense, _ = dense_group_terms(perm, radix)
            for a in range(16):
                for k in range(16):
                    T, classes = group_terms(perm, a, k, radix)
                    Tall[(gname, pname, a, k)] = T
                    bad += T != dense[(a, k)]
                    if a:
                        lo = sum(16 / sz for (_, sz, _) in classes)
                        if not (lo - 1e-9 <= T <= len(classes) * 16):
                            bnd_ok = False
                        if gname == "Z16":
                            pair_vals |= {cnt for (_, sz, cnt) in classes if sz == 2}
    exp.check("PA1", bad == 0, f"{bad} mismatching labels of {4*4*256}")
    z_plus = {Tall[("Z16", "plus5", a, k)] for a in range(16) for k in range(16)}
    g_plus = max(Tall[("GF2^4", "plus5", a, k)] for a in range(16) for k in range(16))
    g_lin = {Tall[("GF2^4", "gf2lin", a, k)] for a in range(16) for k in range(16)}
    z_lin = max(Tall[("Z16", "gf2lin", a, k)] for a in range(16) for k in range(16))
    exp.check("PA2", z_plus == {1} and g_plus > 1 and g_lin == {1} and z_lin > 1,
              f"+5: Z16 {sorted(z_plus)}, GF2^4 max {g_plus}; linear: GF2^4 {sorted(g_lin)}, Z16 max {z_lin}")
    exp.check("PA3", bnd_ok, "all a != 0 labels")
    exp.fail_check("CA2", pair_vals != {8}, f"2-point class supports in Z16: {sorted(pair_vals)}")

    exp.section("PA4/CA1  exhaustive over all 40320 permutations of each order-8 group")
    P = all_perms8()
    # validate the vectorised counter against group_terms on the first 40 permutations
    vbad = 0
    maxima, rmax = {}, {}
    for gname, radix in GROUPS8.items():
        t0 = time.time()
        T, R = exhaustive_T(P, radix)
        for i in range(40):
            for a in range(1, 8):
                for k in range(8):
                    vbad += T[i, a, k] != group_terms(P[i], a, k, radix)[0]
        maxima[gname] = int(T[:, 1:, :].max())
        rmax[gname] = int(R[:, 1:].max())
        exp.log(f"{gname}: max T {maxima[gname]}, max R_a {rmax[gname]} ({time.time()-t0:.0f}s)")
    exp.log(f"vectorised counter vs group_terms on 40 perms x 56 labels x 3 groups: {vbad} mismatches")
    exp.check("PA4", vbad == 0 and maxima == {"GF2^3": 16, "Z8": 56, "Z2xZ4": 56}, f"maxima {maxima}")
    exp.fail_check("CA1", rmax["Z8"] > 4, f"max R_a in Z8 = {rmax['Z8']} vs XOR cap 4")

    exp.section("PB  ToffoliModExp(7, 3, n_exp=3): one clock-shift operator, two bases")
    me = ToffoliModExp(7, 3, n_exp=3)
    qc = me.build()
    perm = walsh.classical_permutation(qc).astype(np.int64)
    low, t = me.anc + 1, me.n_exp
    assert me.exp == list(range(low, low + t)) and qc.n == low + t
    exp.log(f"{qc.n} qubits; exponent qubits {me.exp}")
    with ProcessPoolExecutor(max_workers=7) as pool:
        res = sorted(pool.map(part_b_worker, [(perm, a, low, t) for a in range(1, 8)]))
    pb1 = all(dig == [a] for (a, dig, *_rest) in res)
    exp.check("PB1", pb1, "exponent digits of X-parts: " + str({a: dig for (a, dig, *_r) in res}))
    t4 = next(r for r in res if r[0] == 4)
    f4, _ = formula_terms(perm, 1 << me.exp[-1], 0)
    exp.check("PB2", t4[4] == len(f4), f"Pauli count {t4[4]} vs C99 formula {len(f4)}")
    rows = []
    for (a, _dig, Tm, Rm, Tp, Rp, dt) in res:
        rows.append(f"a={a}: mixed T={Tm} (R={Rm}); Pauli T={Tp} (R={Rp}); ratio {Tp/Tm:.3f} [{dt:.0f}s]")
        exp.log(rows[-1])
    exp.log(f"mixed-basis total for all seven shifts (disjoint X-parts by PB1): {sum(r[2] for r in res)}")
    exp.check("PB3", len(res) == 7, "recorded; see rows above for which hypothesis holds")
    exp.finish()
