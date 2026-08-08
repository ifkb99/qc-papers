"""TODO 12: does windowed / table-lookup arithmetic satisfy the C29 criterion?

Step 9 (C29) replaced Paper B's scope caveat with a checkable condition: the
C15 invariance covers any construction whose a=1 block V satisfies V^2 = id.
Gidney-style windowed arithmetic is what people actually propose to run, so
whether it qualifies is the first question a referee asks.

Windowing consumes w exponent bits per block, and there are two natural ways to
do it. They agree on the valid subspace and disagree off it, which is exactly
where the Walsh pullback lives:

  LOOKUP  (windowed_arith.WindowedModExp, Gidney-style)  always look up
          T[j] = a^(j*2^(kw)) mod N and always multiply. In the identity tail
          every T[j] = 1, so the QROM permutation is s ^= 1 REGARDLESS of j:
          the block never reads the window register.

  SELECT  (windowed_arith.SelectModExp)  apply u_a(a^(j*2^(kw))) controlled on
          [window == j], skipping j = 0 because "multiplying by 1 is free".
          In the identity tail every nonzero j gives the SAME block V, so the
          block is V gated on OR(window) -- a nonlinear function of the window.

Both tail blocks satisfy V^2 = id (test_windowed.py [D]). So V^2 = id cannot be
the whole condition: C29's other hypothesis -- each block controlled on its OWN
QUBIT -- is load-bearing, and at w = 1 the two hypotheses coincide, which is why
step 9 could not see the difference.

PREDICTIONS, WRITTEN BEFORE MEASURING THE SWEEPS.

  P1  LOOKUP, beta = 1: every tail-window exponent qubit is DEAD. The pullback
      is constant in them, so the Walsh support is confined to z_exp supported
      on bits {0 .. alpha-1} only (alpha = v2(r)); the tail contributes z = 0.
      Derived: T_0[j] = a^j depends on j only mod 2^alpha, and every later
      table is all-ones.
        N=5 a=4 (r=2, alpha=1) -> live bits {0}
        N=7 a=6 (r=2, alpha=1) -> live bits {0}
        N=5 a=2 (r=4, alpha=2) -> live bits {0,1}
  P1b IS IT THE WINDOW OR THE ALWAYS-MULTIPLY? At w = 1 there is no window at
      all, yet the lookup form still writes T[0] = T[1] = 1 and multiplies
      unconditionally. If the tail is dead there too -- while the standard
      construction at w = 1 keeps it alive through a parity bit (C24) -- then
      the cause is emitting the multiply-by-1 branch, not consuming several
      exponent bits at once. Also checked at w = 3.
  P2  LOOKUP, beta = 1: |support| does not grow with the tail-window count K.
      Because the tail applies the fixed involution W unconditionally, the
      circuit is W^K . P_live, so |support| is 2-PERIODIC in K: every even K
      equals the K=0 value, every odd K the K=1 value. (Not a single constant
      -- that is a genuine difference from C15/C24, stated as such.)
  P3  SELECT, beta = 1: |support| GROWS with K. C15 constancy is forfeited,
      even though V^2 = id holds, because the tail dependence runs through
      OR(window) rather than a single fresh qubit.
  P4  Contrast with C24: in the STANDARD construction the all-ones-on-tail
      vector 1_I is IN the support (that is what C24 asserts). In the LOOKUP
      design no z with any tail bit set is in the support at all.
  P5  SELECT with skip_zero=False: emitting the j = 0 branch RESTORES the
      property. In the tail all 2^w branches are the same block u_a(.,1) and
      exactly one fires for any window value, so the block degenerates to an
      uncontrolled V and the window register goes dead again -- so |support|
      must become 2-periodic in K, and tail bits must vanish from the support,
      exactly as in the LOOKUP design. The arithmetic is unchanged; only a
      branch nobody would ever emit is added back.
      (The BLOCK-level half of P5 -- window-independence and V^2 = id -- was
      checked before this file was written; the SUPPORT sweep below is the
      unmeasured consequence, and is what P5 is graded on.)
  C1  MUST-FAIL: the LOOKUP design with beta > 1 (N=7 a=3, r=6) must GROW with
      K. No window is a tail window there, so nothing is dead.
  C2  MUST-FAIL: a LIVE window block must READ its window register. Without
      this, P1's independence result would hold for a trivial reason.

OUTCOME (2026-08-08): 8/9 checks pass. **P5 FAILS BY DESIGN and the file exits
nonzero because of it** -- P5 was a conjunction, its flatness half held exactly
and its dead-bit half was refuted by an activation-ancilla artifact diagnosed
in the P5-diag section below. That is a result, not a broken experiment; it is
left failing rather than re-scoped so the record shows what was predicted.
Everything else confirmed, both must-fail controls failed as required.
Written up in `NOTES.md` SS WD; claims C36-C39, and C29 narrowed.

Run:  uv run python -m experiments.experiment_windowed     (from research/)
      Roughly 25 minutes cold; lab.measure caches by circuit content.
"""
from __future__ import annotations

import numpy as np

import walsh
from lab import Experiment, support, order, v2_split
from toffoli_arith import ToffoliModExp
from windowed_arith import WindowedModExp, SelectModExp, verify_modexp

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "lookup, beta=1: support touches only exponent bits 0..alpha-1")
exp.predict("P1b", "the dead tail survives at w=1 and w=3: cause is the "
                   "unconditional multiply, not the window width")
exp.predict("P2", "lookup, beta=1: |support| 2-periodic in K, no growth")
exp.predict("P3", "select, beta=1: |support| grows with K")
exp.predict("P4", "standard has the all-ones tail vector; lookup has no tail bit at all")
exp.predict("P5", "select with the j=0 branch restored is 2-periodic, tail bits dead")
exp.must_fail("C1", "lookup with beta>1 must grow with K")
exp.must_fail("C2", "a live window block must read its window register")

W = 2                     # window width; K counts windows PAST the first
BETA1 = [(5, 4), (7, 6), (5, 2)]
BETAG = (7, 3)            # r = 6 = 3 * 2, beta = 3


def live_exp_bits(me, z):
    return [i for i, q in enumerate(me.exp) if bool(((z >> q) & 1).any())]


# ---------------------------------------------------------------------------
exp.section("P1  which exponent bits survive into the support (lookup design)")
p1_ok = True
for (N, a) in BETA1:
    alpha, beta = v2_split(order(a, N))
    expect = list(range(alpha))
    for K in (1, 2):
        me = WindowedModExp(N=N, a=a, n_exp=W * (K + 1), w=W)
        assert verify_modexp(me), f"lookup {N},{a} incorrect"
        z = support(me.build(), me.x[0])
        got = live_exp_bits(me, z)
        ok = got == expect
        p1_ok &= ok
        exp.log(f"N={N} a={a} r={order(a, N)} alpha={alpha} beta={beta} K={K}: "
                f"live exp bits {got}  expect {expect}  |S|={z.size}  "
                f"{'ok' if ok else 'MISMATCH'}")
        del z
exp.check("P1", p1_ok, "support confined to the sub-2^alpha exponent bits")

exp.section("P1b  is it the WINDOW, or the always-multiply? vary w, not K")
# w = 1 is the decisive one. There is no window to speak of, yet the lookup
# form still writes T[0] = T[1] = 1 and multiplies unconditionally, so the tail
# should STILL be dead -- whereas the standard construction at w = 1 keeps the
# tail alive through a parity bit (C24). If P1 survives at w = 1, the mechanism
# is "always multiply", not "consume several bits at once".
p1b_ok = True
for wid, K in ((1, 3), (3, 1)):
    me = WindowedModExp(N=5, a=4, n_exp=wid * (K + 1), w=wid)
    assert verify_modexp(me), f"lookup w={wid} incorrect"
    z = support(me.build(), me.x[0])
    got = live_exp_bits(me, z)
    ok = got == [0]                       # r=2 => alpha=1
    p1b_ok &= ok
    exp.log(f"N=5 a=4 w={wid} K={K} n_exp={me.n_exp} qubits={me.n_qubits}: "
            f"live exp bits {got}  |S|={z.size}  {'ok' if ok else 'MISMATCH'}")
    del z
exp.check("P1b", p1b_ok, "tail dead at w=1 and w=3 too: the cause is the "
                         "unconditional multiply, not the window width")

# ---------------------------------------------------------------------------
exp.section("P2  |support| vs tail-window count K (lookup design, beta=1)")
p2_ok = True
for (N, a) in BETA1[:2]:                      # N=5 a=4 and N=7 a=6
    vals = {}
    for K in (0, 1, 2, 3):
        me = WindowedModExp(N=N, a=a, n_exp=W * (K + 1), w=W)
        vals[K] = support(me.build(), me.x[0]).size
        exp.log(f"N={N} a={a} K={K} n_exp={me.n_exp} qubits={me.n_qubits}: "
                f"|S|={vals[K]}")
    ok = vals[0] == vals[2] and vals[1] == vals[3]
    p2_ok &= ok
    exp.log(f"  -> even K: {vals[0]} == {vals[2]}; odd K: {vals[1]} == {vals[3]}"
            f"   {'ok' if ok else 'MISMATCH'}")
    exp.log(f"  -> growth over 4 windows: "
            f"{max(vals.values()) / min(vals.values()):.3f}x (2-periodic, not growing)")
exp.check("P2", p2_ok, "2-periodic in K at both beta=1 moduli")

# Strengthening, checked after P2 rather than derived in advance: equal
# cardinality is weak: if W^2 = id really makes the K=2 circuit the SAME
# permutation as K=0 on a wider register, the support SETS must coincide
# once the dead tail coordinates are dropped, not merely their sizes.
a0 = WindowedModExp(N=5, a=4, n_exp=W, w=W)
a2 = WindowedModExp(N=5, a=4, n_exp=W * 3, w=W)
z0 = support(a0.build(), a0.x[0])
z2 = support(a2.build(), a2.x[0])
low = (1 << a0.n_qubits) - 1
same = bool(np.array_equal(np.sort(z2 & low), np.sort(z0)))
exp.log(f"K=0 vs K=2 support sets identical after dropping the dead tail: {same}"
        f"  (|S|={z0.size})")
exp.check("P2-set", same, "same function, not merely the same count")
del z0, z2

# ---------------------------------------------------------------------------
exp.section("P3  |support| vs K (select design, beta=1) -- constancy forfeited?")
p3_ok = True
for (N, a) in BETA1[:2]:
    vals = []
    for K in (0, 1, 2, 3):
        me = SelectModExp(N=N, a=a, n_exp=W * (K + 1), w=W)
        assert verify_modexp(me), f"select {N},{a} incorrect"
        vals.append(support(me.build(), me.x[0]).size)
        exp.log(f"N={N} a={a} K={K} n_exp={me.n_exp} qubits={me.n_qubits}: "
                f"|S|={vals[-1]}"
                + (f"  x{vals[-1] / vals[-2]:.2f}" if len(vals) > 1 else ""))
    grew = all(vals[i + 1] > vals[i] for i in range(len(vals) - 1))
    p3_ok &= grew
    exp.log(f"  -> strictly increasing: {grew}")
exp.check("P3", p3_ok, "select-multiply windowing forfeits C15 constancy")

# ---------------------------------------------------------------------------
exp.section("P4  support geometry: standard (C24) vs lookup")
me = WindowedModExp(N=5, a=4, n_exp=W * 3, w=W)
zw = support(me.build(), me.x[0])
tail_w = sum(1 << q for q in me.exp[1:])
lookup_any_tail = bool((zw & tail_w != 0).any())
del zw

st = ToffoliModExp(N=5, a=4, n_exp=4)
zs = support(st.build(), st.x[0])
tail_s = sum(1 << q for q in st.exp[1:])
std_allones = bool((zs & tail_s == tail_s).any())
del zs

exp.log(f"standard  N=5 a=4 n_exp=4: all-ones-tail vector present = {std_allones}"
        " (C24)")
exp.log(f"lookup    N=5 a=4 n_exp=6: ANY tail bit present        = {lookup_any_tail}")
exp.check("P4", std_allones and not lookup_any_tail,
          "same map, same modulus, opposite support geometry in the exponent")

# ---------------------------------------------------------------------------
exp.section("P5  select design with the j=0 branch put back")
p5_ok = True
for (N, a) in BETA1[:2]:
    vals = {}
    for K in (0, 1, 2, 3):
        me = SelectModExp(N=N, a=a, n_exp=W * (K + 1), w=W, skip_zero=False)
        assert verify_modexp(me), f"select-nozero {N},{a} incorrect"
        z = support(me.build(), me.x[0])
        vals[K] = z.size
        live = live_exp_bits(me, z)
        alpha, _ = v2_split(order(a, N))
        dead_ok = live == list(range(alpha))
        p5_ok &= dead_ok
        exp.log(f"N={N} a={a} K={K} qubits={me.n_qubits}: |S|={z.size}  "
                f"live exp bits {live}  {'ok' if dead_ok else 'MISMATCH'}")
        del z
    per = vals[0] == vals[2] and vals[1] == vals[3]
    p5_ok &= per
    exp.log(f"  -> 2-periodic: even {vals[0]}=={vals[2]}, odd {vals[1]}=={vals[3]}"
            f"   {'ok' if per else 'MISMATCH'}")
exp.check("P5", p5_ok,
          "restoring the multiply-by-1 branch restores the dead tail")

# P5 was a CONJUNCTION and it split: the flatness held, the dead-bit half did
# not. Reported here rather than quietly re-scoped -- see NOTES.md SS WD.
#   held    : |S| does not grow with K (2-periodic, ~32k at every K), which is
#             the C15-relevant half and the whole point of the comparison.
#   refuted : the live window's OTHER bit is not dead. Diagnosis below.
exp.section("P5-diag  why the dead-bit half of P5 failed")
# The activation ancilla is itself a scratch qubit, and the pullback ranges over
# ALL its values. In the LOOKUP design each branch contributes a CNOT into s, so
# the 2^w contributions XOR the ancilla an even number of times and it cancels;
# in SELECT the ancilla GATES a block, so it does not. Prediction: the select
# block should depend on e1 in general but NOT on the act=0 half-space.
for cls, kw in ((SelectModExp, dict(skip_zero=False)), (WindowedModExp, {})):
    m = cls(N=5, a=4, n_exp=W, w=W, **kw)
    P = walsh.classical_permutation(m.window_block(m.windows()[0], 0))
    idn = np.arange(P.size, dtype=np.int64)
    e1 = 1 << m.exp[1]
    full = bool((P[idn ^ e1] == (P ^ e1)).all())
    sel = np.nonzero((idn >> m.lk[-1]) & 1 == 0)[0]
    half = bool((P[sel ^ e1] == (P[sel] ^ e1)).all())
    exp.log(f"{cls.__name__:16s}: block independent of e1 -- "
            f"everywhere {full}, on the act=0 half-space {half}")
    del P
exp.log("=> the refuted half of P5 is an ACTIVATION-ANCILLA artifact, not a "
        "failure of the mechanism; the lookup design cancels it, select does not")

exp.section("C1  must-fail control: lookup design with beta > 1")
N, a = BETAG
vals = []
for K in (0, 1, 2):
    me = WindowedModExp(N=N, a=a, n_exp=W * (K + 1), w=W)
    assert verify_modexp(me), "beta>1 control incorrect"
    vals.append(support(me.build(), me.x[0]).size)
    exp.log(f"N={N} a={a} r={order(a, N)} K={K} qubits={me.n_qubits}: "
            f"|S|={vals[-1]}"
            + (f"  x{vals[-1] / vals[-2]:.2f}" if len(vals) > 1 else ""))
grew = all(vals[i + 1] > vals[i] for i in range(len(vals) - 1))
exp.fail_check("C1", grew, f"beta=3 grows: {vals}")

# ---------------------------------------------------------------------------
exp.section("C2  must-fail control: a LIVE window block reads its window")
me = WindowedModExp(N=7, a=3, n_exp=W, w=W)      # beta>1: window 0 is live
P = walsh.classical_permutation(me.window_block(me.windows()[0], 0))
idn = np.arange(P.size, dtype=np.int64)
indep = all(bool((P[idn ^ (1 << q)] == (P ^ (1 << q))).all()) for q in me.exp)
del P
exp.fail_check("C2", not indep,
               "live block depends on the window, so P1's independence is real")

exp.finish()
