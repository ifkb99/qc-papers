"""Is the recurring linear/nonlinear GF(2) pattern one theorem? Mostly yes.

Three results in this project have the same shape -- linear/affine ingredients
are free, nonlinear ones are expensive:

  SS L2  (C30/C32)  a linear structure caps density at 1/2; a linear wrap
                    leaves it intact, a Toffoli (degree-2) destroys it
  SS RS  (C34)      one nonlinear monomial demotes it to a CONDITIONAL
                    structure, capping density at 3/4
  SS WD  (C36)      a block gated on an affine function of its controls keeps
                    the C15 invariance; gated on OR it does not

LITERATURE (checked at SOURCE before deriving anything -- Carlet's book
"Boolean Functions for Cryptography and Coding Theory", open access at
math.univ-paris13.fr/~carlet/book-fcts-Bool-vect-crypt-codes.pdf).

  * The m = 0 case is textbook. **Proposition 29**: "Let f be any n-variable
    Boolean function. The derivative D_e f equals the null function (resp.
    function 1) if and only if the support supp(W_f) is included in
    {0_n, e}^perp (resp. in its complement)." That IS C30 (linear structure ->
    hyperplane -> density <= 1/2) and C33 (affine structure -> the opposite
    coset). Both must CITE Proposition 29, not claim it. The functions whose
    Walsh support is an affine subspace are Carlet's PARTIALLY BENT functions.
  * Closest adjacent notion: **covering sequences** (Carlet-Tarannikov 2002),
    book Prop. 60 -- f admits lambda as a covering sequence iff the
    Fourier-Hadamard transform of lambda is constant on supp(W_f). A covering
    sequence is a SINGLE GLOBAL lambda; the "partial" variant allows two
    levels. Ours has a DIFFERENT STRUCTURE VECTOR PER CELL, which is neither.
    The DCC paper body is paywalled and was NOT read; the book (same author,
    citing it on pp. 205/206/319 and reproducing Def. 47 / Prop. 60) and the
    paper's own abstract both indicate a different scope, so the risk is low
    but nonzero. Logged in NOTES.md SS GF.
  * Also adjacent, worth citing: **Maiorana-McFarland** (book SS 5.1.1) --
    functions whose restrictions to each coset are AFFINE. Same flavour,
    much stronger hypothesis, different conclusion.

Honest altitude for what follows: it is an easy COROLLARY of Proposition 29
plus the standard decomposition of the Walsh transform over a coset partition.
Both ingredients are textbook; we did not find the combination stated, and it
may well be folklore. The contribution is not the corollary, it is noticing
that this project's density caps are all instances of it -- and that our own
conjectured ladder was the wrong shape.

THE DERIVATION (done before measuring; the point of the experiment is to
falsify it, and one prediction below is designed to break the SS RS ladder).

Let C be a set of "condition" coordinates, |C| = m, splitting F_2^n into cells
H_u = {y : y_C = u}. Suppose on each cell f has a linear structure w_u with
(w_u)_C = 0, i.e. f(y ^ w_u) = f(y) for y in H_u. Writing z = (z_C, z'),

    c_z = 2^-n * SUM_u (-1)^(u . z_C) * A_u(z'),
    A_u(z') = SUM_{y in H_u} (-1)^f(y) (-1)^(y' . z')

and substituting y -> y ^ w_u inside A_u gives A_u = (-1)^(w_u . z') A_u, so

    A_u(z') = 0  whenever  w_u . z' = 1.

Hence c_z = 0 whenever w_u . z' = 1 for EVERY u. So the support avoids

    E = {z : w_u . z = 1 for all u},

a coset of codimension d = dim span{w_u} -- PROVIDED that system is
consistent. Therefore:

    density <= 1 - 2^-d     if the system {w_u . z = 1} is consistent
    no cap at all           if it is not.

Two consequences that matter, and they contradict how SS RS guessed this:

  * the cap is set by the DIMENSION OF THE SPAN of the per-cell structures,
    not by the number of conditioning levels. The "1 - 2^-(k+1) ladder" is
    only the special case where each new level contributes one new independent
    vector. Two cells sharing the SAME structure give d = 1 and a 1/2 cap, not
    3/4.
  * consistency is a PARITY condition. Any linear dependency among the w_u
    with ODD support (e.g. w_0 ^ w_1 ^ w_2 = 0) makes the system unsolvable,
    E empty, and the cap vanishes entirely -- at the same number of cells and
    the same d as a case that does have a cap.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Support avoids E exactly (0 violations), whenever E is nonempty.
  P2  Density <= 1 - 2^-d, and close to it for random plants. Sweeping
      d = 1, 2, 3 gives caps 1/2, 3/4, 7/8.
  P3  d is the span dimension, not the conditioning depth: 2 cells sharing
      one structure => d = 1 => 1/2 cap, NOT 3/4. This is the prediction that
      breaks the SS RS ladder.
  P4  The real modexp instances are cases of this. v4's two slice structures
      are msb^anc (on t0=0) and msb (on t0=1); they span d = 2, so E =
      {z : z_msb ^ z_anc = 1 and z_msb = 1} = {z_msb=1, z_anc=0} -- EXACTLY
      the quadrant C34 found empty, and cap 3/4 exactly as measured. C30 is
      the m = 0, d = 1 case.
  P5  AFFINE INVARIANCE (the scope claim, tested rather than asserted). The
      derivation used coordinate cells {y : y_C = u} and structures with zero
      components on C. Everything in sight is affine-covariant, so composing f
      with an invertible linear map L should carry cells to cosets of an
      arbitrary subspace and structures to L^-1 w_u, with the cap unchanged and
      E transformed by L^-T. If that holds, the law covers arbitrary subspace
      partitions, not just coordinate ones.
  C1  MUST-FAIL: with an ODD dependency among the w_u the cap must VANISH --
      density goes to ~1 even though d is unchanged. If a cap survives there,
      the parity half of the derivation is wrong.
  C2  MUST-FAIL: an unplanted random function must show no cap and no
      structure at all, or the plant is not doing anything.

Run:  uv run python -m experiments.experiment_gf2law     (from research/)
"""
from __future__ import annotations

import numpy as np

from lab import (Experiment, support, density, structures, rank_kernel,
                 build_modexp, verify_correctness, slice_fn, quadrant_counts)
import walsh

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "support avoids E exactly whenever E is nonempty")
exp.predict("P2", "density <= 1 - 2^-d, close to it; d=1,2,3 -> 1/2, 3/4, 7/8")
exp.predict("P3", "two cells sharing one structure give d=1 (1/2), not 3/4")
exp.predict("P4", "v4's slice structures span d=2 and give C34's empty quadrant")
exp.predict("P5", "affine invariance: the law survives a random linear change "
                  "of variables, so it covers arbitrary subspace partitions")
exp.must_fail("C1", "an ODD dependency among the w_u must destroy the cap")
exp.must_fail("C2", "an unplanted random function must show no cap")

N = 14                                   # bits in the planted experiments
RNG = np.random.default_rng(20260808)


# ---------------------------------------------------------------------------
def span_dim(ws, n):
    return rank_kernel(list(dict.fromkeys(ws)), n)[0]


def consistent(ws, n):
    """Is {w . z = 1 for all w in ws} solvable? Fails iff some dependency
    among the ws has odd support."""
    uniq = list(dict.fromkeys(ws))
    r = rank_kernel(uniq, n)[0]
    r_aug = rank_kernel([(w << 1) | 1 for w in uniq], n + 1)[0]
    return r_aug == r


def excluded(ws, n):
    """E = {z : w . z = 1 for all w in ws}, as a boolean mask over all z."""
    z = np.arange(1 << n, dtype=np.int64)
    keep = np.ones(1 << n, dtype=bool)
    for w in dict.fromkeys(ws):
        par = np.zeros(1 << n, dtype=np.int64)
        ww = int(w)
        while ww:
            b = ww & -ww
            par ^= (z >> (b.bit_length() - 1)) & 1
            ww ^= b
        keep &= (par == 1)
    return keep


def plant(n, cond_bits, ws, rng):
    """Random f whose restriction to cell u has linear structure ws[u]."""
    m = len(cond_bits)
    assert len(ws) == (1 << m)
    for w in ws:
        assert all(not ((w >> c) & 1) for c in cond_bits), \
            "structures must not move the condition coordinates"
    bits = rng.integers(0, 2, size=1 << n).astype(np.int64)
    idx = np.arange(1 << n, dtype=np.int64)
    uval = np.zeros(1 << n, dtype=np.int64)
    for i, c in enumerate(cond_bits):
        uval |= ((idx >> c) & 1) << i
    g = np.empty(1 << n, dtype=np.int64)
    for u in range(1 << m):
        sel = np.nonzero(uval == u)[0]
        rep = np.minimum(sel, sel ^ np.int64(ws[u]))
        g[sel] = bits[rep]
    return g


def fn_support_local(g):
    chi = np.where(g == 1, -1.0, 1.0)
    c = walsh.wht(chi) / g.size
    return np.nonzero(np.abs(c) > 1e-9)[0].astype(np.int64)


def run_case(name, cond_bits, ws):
    g = plant(N, cond_bits, ws, RNG)
    zs = fn_support_local(g)
    d = span_dim(ws, N)
    ok_sys = consistent(ws, N)
    cap = 1 - 2.0 ** (-d) if ok_sys else 1.0
    dens = zs.size / (1 << N)
    E = excluded(ws, N)
    hits = int(E[zs].sum())
    exp.log(f"{name:34s} cells={1 << len(cond_bits)} d={d} "
            f"consistent={str(ok_sys):5s} |E|={int(E.sum()):6d} "
            f"cap={cap:.4f} density={dens:.4f} violations={hits}")
    return d, ok_sys, cap, dens, hits


# ---------------------------------------------------------------------------
exp.section("P1/P2  planted conditional structures: does the support avoid E?")
cases = [
    ("d=1, no conditioning (C30 form)", [], [0b0000000000110]),
    ("d=1, 2 cells SHARING one w", [0], [0b0000000000110, 0b0000000000110]),
    ("d=2, 2 cells, independent", [0], [0b0000000000110, 0b0000000001010]),
    ("d=3, 4 cells, independent", [0, 1],
     [0b0000000000100, 0b0000000001000, 0b0000000010000, 0b0000000011100]),
]
res = {}
p1_ok = p2_ok = True
for name, cb, ws in cases:
    d, ok_sys, cap, dens, hits = run_case(name, cb, ws)
    res[name] = (d, cap, dens)
    if ok_sys:
        p1_ok &= hits == 0
        p2_ok &= dens <= cap + 1e-12
exp.check("P1", p1_ok, "zero support elements inside E, all consistent cases")
exp.check("P2", p2_ok, "density never exceeds 1 - 2^-d")
exp.log(f"  caps hit: " + ", ".join(
    f"d={res[n][0]} cap={res[n][1]:.4f} got={res[n][2]:.4f}" for n in res))

exp.section("P3  is it the span dimension or the conditioning depth?")
d_share = res["d=1, 2 cells SHARING one w"][0]
dens_share = res["d=1, 2 cells SHARING one w"][2]
d_indep = res["d=2, 2 cells, independent"][0]
dens_indep = res["d=2, 2 cells, independent"][2]
exp.check("P3", d_share == 1 and dens_share <= 0.5 + 1e-12
          and d_indep == 2 and 0.5 < dens_indep <= 0.75 + 1e-12,
          f"same 2 cells: shared w -> d=1 density {dens_share:.4f} (<=1/2); "
          f"independent -> d=2 density {dens_indep:.4f} (in (1/2, 3/4])")

exp.section("P5  affine invariance: does the law survive a change of basis?")


def rows_from_cols(cols, n):
    return [sum(((cols[j] >> i) & 1) << j for j in range(n)) for i in range(n)]


def invert_rows(rows, n):
    m = [[rows[i], 1 << i] for i in range(n)]
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if (m[i][0] >> c) & 1), None)
        if p is None:
            return None
        m[r], m[p] = m[p], m[r]
        for i in range(n):
            if i != r and ((m[i][0] >> c) & 1):
                m[i][0] ^= m[r][0]
                m[i][1] ^= m[r][1]
        r += 1
    inv = [0] * n
    for i in range(n):
        inv[m[i][0].bit_length() - 1] = m[i][1]
    return inv


def apply_cols(cols, idx):
    out = np.zeros_like(idx)
    for j, c in enumerate(cols):
        out ^= np.where(((idx >> j) & 1) == 1, np.int64(c), np.int64(0))
    return out


def apply_rows(rows, w):
    return sum((int(r & w).bit_count() & 1) << i for i, r in enumerate(rows))


cond_bits, ws = [0], [0b0000000000110, 0b0000000001010]      # the d=2 case
g0 = plant(N, cond_bits, ws, RNG)
zs0 = fn_support_local(g0)

while True:
    cols = [int(RNG.integers(0, 1 << N)) for _ in range(N)]
    inv_rows = invert_rows(rows_from_cols(cols, N), N)
    if inv_rows is not None:
        break

idx = np.arange(1 << N, dtype=np.int64)
h = g0[apply_cols(cols, idx)]                    # h(x) = g(Lx)
vs = [apply_rows(inv_rows, w) for w in ws]       # transformed structures
zsh = fn_support_local(h)

# cells of h are the preimages of the coordinate cells -- no longer coordinate
# aligned, which is the whole point.
Lx = apply_cols(cols, idx)
cellh = (Lx >> cond_bits[0]) & 1
pointwise = all(
    bool(np.array_equal(h[sel ^ np.int64(vs[u])], h[sel]))
    for u in (0, 1)
    for sel in [np.nonzero(cellh == u)[0]])

d_h = span_dim(vs, N)
Eh = excluded(vs, N)
hits_h = int(Eh[zsh].sum())
exp.log(f"transformed: structures {[hex(v) for v in vs]}  d={d_h}  "
        f"|E|={int(Eh.sum())}  density={zsh.size / (1 << N):.4f}  "
        f"violations={hits_h}  pointwise-structures={pointwise}")
exp.check("P5", pointwise and hits_h == 0 and d_h == 2
          and zsh.size == zs0.size,
          "law survives a random GF(2) change of basis; cells are now cosets "
          "of a generic subspace, cap and support size unchanged")

exp.section("C1  must-fail: an ODD dependency destroys the cap")
# w0 ^ w1 ^ w2 = 0 with three cells constrained -> 1^1^1 = 1 != 0, unsolvable.
w0, w1 = 0b0000000000100, 0b0000000001000
odd = [w0, w1, w0 ^ w1, w0]
d_odd, ok_odd, cap_odd, dens_odd, hits_odd = run_case(
    "ODD dependency, 4 cells", [0, 1], odd)
naive_cap = 1 - 2.0 ** (-span_dim(odd, N))
exp.fail_check("C1", (not ok_odd) and dens_odd > naive_cap,
               f"span d={d_odd} would predict cap {naive_cap:.4f}, "
               f"measured density {dens_odd:.4f} -- cap gone, as derived")

# and the even-dependency partner: same cell count, cap survives
w2 = 0b0000000010000
even = [w0, w1, w2, w0 ^ w1 ^ w2]
d_ev, ok_ev, cap_ev, dens_ev, hits_ev = run_case(
    "EVEN dependency, 4 cells", [0, 1], even)
exp.log(f"  -> parity is the whole difference: same 4 cells, "
        f"odd dependency density {dens_odd:.4f} vs even {dens_ev:.4f} "
        f"(cap {cap_ev:.4f})")

exp.section("C2  must-fail: unplanted random function")
g = RNG.integers(0, 2, size=1 << N).astype(np.int64)
zs = fn_support_local(g)
st = structures(zs, N)
dens_rand = zs.size / (1 << N)
exp.fail_check("C2", dens_rand > 0.9 and st == [],
               f"random density {dens_rand:.4f}, structures {st}")

# ---------------------------------------------------------------------------
exp.section("P4  the real modexp instances are cases of the same law")
me = build_modexp(N=5, a=2, n_exp=2, wraps=("msb_t0_anc",))
assert verify_correctness(me)
qc = me.build()
g = ((walsh.classical_permutation(qc) >> me.x[0]) & 1).astype(np.int64)
msb, anc, t0 = me.b[me.m - 1], me.anc, me.t[0]

# per-cell structures, read off the two t0 slices
slices = {v: slice_fn(g, me.n_qubits, t0, v) for v in (0, 1)}
sl_st = {v: structures(fn_support_local(s), me.n_qubits - 1)
         for v, s in slices.items()}
exp.log(f"t0=0 slice structures: {[(hex(w), e) for w, e in sl_st[0]]}")
exp.log(f"t0=1 slice structures: {[(hex(w), e) for w, e in sl_st[1]]}")

zs = support(qc, me.x[0])
quad = quadrant_counts(zs, msb, anc)
dens = zs.size / (1 << me.n_qubits)
# E predicted from the derivation: w_0 = msb^anc, w_1 = msb  =>  d = 2
ws_real = [(1 << msb) | (1 << anc), (1 << msb)]
d_real = span_dim(ws_real, me.n_qubits)
exp.log(f"predicted d={d_real} -> cap {1 - 2.0 ** (-d_real):.4f}; "
        f"measured density {dens:.4f}; quadrant counts {quad}")
exp.check("P4", d_real == 2 and quad[(1, 0)] == 0 and dens <= 0.75 + 1e-12,
          "C34's empty quadrant IS E for the two slice structures, cap 3/4")

exp.finish()
