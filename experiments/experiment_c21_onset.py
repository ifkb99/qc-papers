"""TODO 12e -- C21's onset at alpha = 3 and alpha = 4, now that the GPU makes it
affordable.

C21 says the Walsh support of a beta=1 modexp pullback LOCKS at
n_exp = v2(r) + 1, because that is the first width carrying an identity block
(a^(2^i) = 1 iff i >= alpha). It is verified exactly at alpha = 1, 2. At
alpha = 3, 4 the record says only "still growing at the largest width we can
reach, as predicted" -- i.e. the sharpest consequence of the onset rule has
never been measured, purely for want of compute. The GPU backend (LAB_GPU=1)
puts q = 27 within reach, so this closes it.

Everything below follows from the C23/C24 theorem (experiment_c15_proof5.py):
the identity tail composes to V^p with p = XOR of the tail exponent bits, and
averaging the Walsh character over the tail confines the support to
z_I in {0, 1_I}. That proof is width-agnostic, so it predicts the onset AND the
structure at every alpha, not just the two already measured.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  alpha = 3 (N=17, a=2, r=8): |support| STRICTLY GROWS across
      n_exp = 1,2,3 and then |S(4)| == |S(5)| exactly. First lock at
      n_exp = alpha + 1 = 4.

  P2  alpha = 4 (N=17, a=3, r=16): strictly grows across n_exp = 1..4, then
      |S(5)| == |S(6)| exactly. First lock at n_exp = alpha + 1 = 5.
      This is the deepest onset ever measured here.

  P3  alpha = 3 at a SECOND modulus (N=41, a=3, r=8): same onset, n_exp = 4.
      Guards against the rule being an artifact of N = 17.

  P4  Structure, not just counts (C24 (iv)): for every locked width, every
      support element z restricted to the identity tail I is either all-zeros
      or all-ones. Zero violations.
      SCOPE CORRECTION, made after C2 refused to fail (see the note there and
      the OUTCOME below): this has content only at |I| >= 2. At |I| = 1 the
      tail restriction is a single bit, which is "all-zeros or all-ones" by
      definition, so the n_exp = alpha+1 rows are vacuous and are excluded
      from the check. The prediction is unchanged; the set of widths that can
      testify to it is smaller than first written.

  P5  Sharper than C24 as previously tested, and DERIVED rather than fitted:
      proof5 checked that the two halves are individually constant in SIZE.
      The proof gives more -- the coefficient c(z_rest, branch) does not
      mention |I| at all -- so the two halves must be identical as SETS. Write
      each z as (z_rest, tailflag) with tailflag in {0,1}; then the signature
      set is bit-for-bit identical at n_exp = alpha+1 and alpha+2.

  C1  MUST FAIL -- vacuity. If the support were flat from n_exp = 1 the whole
      test would confirm nothing. So growth BELOW the onset is required:
      |S(k)| != |S(k+1)| for every k < alpha.

  C2  MUST FAIL -- matched beta control. N=41, a=6 has r = 40 = 5 * 2^3: the
      SAME alpha = 3 and the SAME modulus, width and gate count as P3, with
      beta = 5 the only thing that differs. There is no identity block at any
      width, so it must NOT lock at n_exp = 4, and P4's tail confinement must
      also fail for it.

OUTCOME (2026-08-08): 7/7. C21's onset is now measured, not merely predicted,
at alpha = 3 (two moduli) and alpha = 4; P4/P5 confirm the C24 geometry and
strengthen it from equal counts to identical sets. **The must-fail control did
its job on the first pass**: C2 initially PASSED, which flagged that the tail
confinement is vacuous at |I| = 1 -- and P4 was reading three such rows as
evidence. Both were re-scoped to |I| >= 2 and the numbers are unchanged; what
changed is how much they are entitled to say.

Run:  LAB_GPU=1 uv run python -m experiments.experiment_c21_onset

Cost note: q = 3*ceil(log2 N) + 4 + n_exp, so the N=41 rows run to q = 27
(134M basis states, ~14k gates). CPU-only this is a multi-hour job; on one
A4500 it is minutes. Measurement uses the EXACT integer FWHT path
(accel.pullback_support_exact), so the support test is `!= 0` rather than a
magnitude threshold.
"""
from __future__ import annotations

import time

import numpy as np

from lab import Experiment, build_modexp, order, v2_split, support
from windowed_arith import replay, read

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "alpha=3 (N=17,a=2): grows through n_exp=3, locks at n_exp=4")
exp.predict("P2", "alpha=4 (N=17,a=3): grows through n_exp=4, locks at n_exp=5")
exp.predict("P3", "alpha=3 at N=41,a=3: same onset, n_exp=4")
exp.predict("P4", "locked widths: every z has z_I in {0, 1_I}, 0 violations")
exp.predict("P5", "the (z_rest, tailflag) signature SETS coincide across widths")
exp.must_fail("C1", "vacuity: support must NOT be flat below the onset")
exp.must_fail("C2", "matched beta=5 control (N=41,a=6, same alpha=3, same q): "
                    "must not lock, and tail confinement must fail")


# -- correctness gate, cheap enough to run at every width --------------------

def verify(me) -> bool:
    """x holds a^e mod N and every scratch qubit returns to 0, for all e.

    Single-state replay rather than the state vector: verify_correctness()
    evolves a 2^q amplitude array, hopeless at q = 27, while this is O(gates)
    per exponent value.
    """
    qc = me.build()
    scratch = list(me.b) + list(me.t) + [me.c0, me.anc]
    for e in range(1 << me.n_exp):
        y0 = 0
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1:
                y0 |= 1 << eq
        y = replay(qc, y0)
        if read(y, me.x) != pow(me.a, e, me.N):
            return False
        if read(y, me.exp) != e:
            return False
        if any((y >> q_) & 1 for q_ in scratch):
            return False
    return True


# -- one instance: sweep n_exp by 1, record support size and tail structure --

def sweep(N: int, a: int, widths, label: str):
    """Returns {n_exp: (q, |S|, tail_violations, signature, |I|)}.

    `tail_violations` is meaningless at |I| = 1: a one-bit tail is 0 or
    all-ones by definition, so the confinement test passes for free. Callers
    must filter on |I| >= 2 -- see the note at C2, where the missing filter
    made the control pass when it had to fail.
    """
    r = order(a, N)
    al, beta = v2_split(r)
    pred = f"predicted constant for n_exp >= {al + 1}" if beta == 1 \
        else "predicted: never constant"
    exp.log(f"{label}: N={N} a={a} r={r} -> alpha={al} beta={beta}, {pred}")
    print(f"    {'n_exp':>5} {'q':>3} {'gates':>7} {'|support|':>11} "
          f"{'density':>8} {'vs prev':>8} {'|I|':>4} {'tail z_I=0':>11} "
          f"{'z_I=1_I':>9} {'other':>9} {'secs':>7}")
    out = {}
    for ne in widths:
        t0 = time.time()
        me = build_modexp(N=N, a=a, n_exp=ne)
        qc = me.build()
        assert verify(me), f"{label} n_exp={ne} does not compute a^e mod N"
        zs = support(qc, me.x[0], exact=True)

        # tail = exponent bits alpha .. n_exp-1 (empty when ne <= alpha)
        tail_bits = [me.exp[i] for i in range(al, ne)]
        tmask = 0
        for q_ in tail_bits:
            tmask |= 1 << q_
        if tmask:
            sub = zs & tmask
            n0 = int(np.count_nonzero(sub == 0))
            n1 = int(np.count_nonzero(sub == tmask))
            other = int(zs.size - n0 - n1)
            # Canonical form of the pair (z_rest, tailflag). z_rest occupies
            # the bits below the tail and the flag goes in the first tail bit,
            # which z_rest cannot reach -- so the encoding is injective and two
            # widths agree as SETS iff their sorted key arrays are equal. (A
            # frozenset of tuples says the same thing and costs ~40x the
            # memory: at 33.5M support elements that is the difference between
            # 268 MB and an OOM.)
            rest_mask = (1 << me.exp[al]) - 1     # everything below the tail
            sig = np.sort((zs & rest_mask)
                          | ((sub == tmask).astype(np.int64) << me.exp[al]))
        else:
            n0 = n1 = other = -1
            sig = None

        prev = out.get(ne - 1)
        cmp = "-" if prev is None else ("SAME" if prev[1] == zs.size else "grew")
        dens = zs.size / (1 << qc.n)
        shown = other if len(tail_bits) >= 2 else ("vacuous" if tmask else "-")
        print(f"    {ne:5d} {qc.n:3d} {len(qc.logical):7d} {zs.size:11,d} "
              f"{dens:8.4f} {cmp:>8} {len(tail_bits):4d} "
              f"{n0 if tmask else '-':>11} {n1 if tmask else '-':>9} "
              f"{shown:>9} {time.time() - t0:7.1f}",
              flush=True)
        out[ne] = (qc.n, int(zs.size), other, sig, len(tail_bits))
    print()
    return al, beta, out


def locks_at(res, onset):
    """Grew at every step below `onset`, identical at every step from it on."""
    ks = sorted(res)
    grew = all(res[k][1] != res[k + 1][1]
               for k in ks if k + 1 in res and k + 1 <= onset)
    same = all(res[k][1] == res[k + 1][1]
               for k in ks if k + 1 in res and k >= onset)
    return grew, same


# ---------------------------------------------------------------------------
exp.section("P1  alpha = 3, N = 17, a = 2  (r = 8)")
al1, b1, r17a2 = sweep(17, 2, [1, 2, 3, 4, 5], "P1")
grew1, same1 = locks_at(r17a2, al1 + 1)
exp.check("P1", same1 and grew1,
          f"sizes {[r17a2[k][1] for k in sorted(r17a2)]}, "
          f"lock at n_exp={al1 + 1}")

exp.section("P2  alpha = 4, N = 17, a = 3  (r = 16)")
al2, b2, r17a3 = sweep(17, 3, [1, 2, 3, 4, 5, 6], "P2")
grew2, same2 = locks_at(r17a3, al2 + 1)
exp.check("P2", same2 and grew2,
          f"sizes {[r17a3[k][1] for k in sorted(r17a3)]}, "
          f"lock at n_exp={al2 + 1}")

exp.section("P3  alpha = 3 at a second modulus, N = 41, a = 3  (r = 8)")
al3, b3, r41a3 = sweep(41, 3, [1, 2, 3, 4, 5], "P3")
grew3, same3 = locks_at(r41a3, al3 + 1)
exp.check("P3", same3 and grew3,
          f"sizes {[r41a3[k][1] for k in sorted(r41a3)]}, "
          f"lock at n_exp={al3 + 1}")

exp.section("P4  tail confinement z_I in {0, 1_I}, at widths where it has content")
# |I| = 1 is excluded: a one-bit tail satisfies "0 or all-ones" by definition,
# so n_exp = alpha+1 rows carry no evidence. Only |I| >= 2 counts.
viol = []
for (al, res, tag) in ((al1, r17a2, "N=17,a=2"), (al2, r17a3, "N=17,a=3"),
                       (al3, r41a3, "N=41,a=3")):
    for k in sorted(res):
        if res[k][4] >= 2:
            viol.append((tag, k, res[k][2]))
exp.log("violations per non-vacuous width:",
        ", ".join(f"{t}/n_exp={k}: {v}" for t, k, v in viol))
exp.check("P4", bool(viol) and all(v == 0 for _, _, v in viol),
          f"{len(viol)} widths with |I| >= 2 checked, "
          f"{sum(1 for _, _, v in viol if v == 0)} clean")

exp.section("P5  the two halves coincide as SETS, not merely in size")
sig_ok, sig_detail = [], []
for (al, res, tag) in ((al1, r17a2, "N=17,a=2"), (al2, r17a3, "N=17,a=3"),
                       (al3, r41a3, "N=41,a=3")):
    ks = [k for k in sorted(res) if k >= al + 1 and res[k][3] is not None]
    for k in ks[1:]:
        same = bool(np.array_equal(res[ks[0]][3], res[k][3]))
        sig_ok.append(same)
        sig_detail.append(f"{tag}: n_exp={ks[0]} vs {k} "
                          f"{'identical' if same else 'DIFFER'}")
for d in sig_detail:
    exp.log(d)
exp.check("P5", bool(sig_ok) and all(sig_ok),
          f"{sum(sig_ok)}/{len(sig_ok)} width pairs identical")

exp.section("C1  vacuity control -- the support must grow below the onset")
below = []
for (al, res, tag) in ((al1, r17a2, "N=17,a=2"), (al2, r17a3, "N=17,a=3"),
                       (al3, r41a3, "N=41,a=3")):
    for k in sorted(res):
        if k + 1 in res and k + 1 <= al:
            below.append((tag, k, res[k][1] != res[k + 1][1]))
exp.log("strict growth below onset:",
        ", ".join(f"{t}/{k}->{k+1}:{'grew' if g else 'FLAT'}"
                  for t, k, g in below))
exp.fail_check("C1", all(g for _, _, g in below),
               f"{sum(g for _, _, g in below)}/{len(below)} steps below the "
               f"onset are strict increases (flat would make P1-P3 vacuous)")

exp.section("C2  matched control -- N = 41, a = 6 (r = 40 = 5*2^3), beta = 5")
al4, b4, r41a6 = sweep(41, 6, [1, 2, 3, 4, 5], "C2")
grew4, same4 = locks_at(r41a6, al4 + 1)
ctrl_tail = [r41a6[k][2] for k in sorted(r41a6) if r41a6[k][4] >= 2]
exp.log(f"beta={b4}, same alpha={al4} and same width as P3; "
        f"nominal-tail violations at |I| >= 2: {ctrl_tail}")
# NOTE, and it is the reason this control exists. The first version of this
# check included the |I| = 1 row, where the control reported 0 violations and
# so PASSED a test it was required to fail. It was not a discovery about
# beta = 5: a one-bit tail is trivially "0 or all-ones". The must-fail control
# caught a vacuous measurement in the *same* line of P4 -- exactly the failure
# mode METHOD.md records from step 9.
exp.fail_check("C2", (not same4) and bool(ctrl_tail)
               and all(v > 0 for v in ctrl_tail),
               f"sizes {[r41a6[k][1] for k in sorted(r41a6)]} -- "
               f"locked={same4} (must be False), "
               f"tail violations {ctrl_tail} (must all be > 0)")

exp.finish()
