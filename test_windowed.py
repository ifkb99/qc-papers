"""Correctness gate for the windowed / table-lookup constructions (TODO 12).

Two things are being verified, and they are different in kind:

  [A]-[C] the ARITHMETIC -- the circuits really do compute a^e mod N, with
          every scratch register returned to |0>. Checked against an
          independent route (statevector at the smallest size, and the
          unrelated ToffoliModExp construction on the valid subspace), not
          only against the replay engine they are built on.

  [D]-[F] the STRUCTURE the experiment then reasons about -- the identity-tail
          block's permutation, its involutivity, and what it reads. These are
          exhaustive over the whole state space, so they are categorical.

Run:  uv run python test_windowed.py
"""
from __future__ import annotations
import numpy as np

import walsh
import statevec as sv
from circuits import Circuit
from toffoli_arith import ToffoliModExp
from windowed_arith import (WindowedModExp, SelectModExp, verify_modexp,
                            replay, read)

FAILED = []


def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra else ''}")
    if not cond:
        FAILED.append(name)


# ---------------------------------------------------------------------------
print("\n[A] replay() agrees with the vectorised permutation replay")
# The new single-state helper must agree with walsh.classical_permutation,
# which is itself gated against the state vector by test_core / test_walsh.
me = SelectModExp(N=5, a=4, n_exp=2, w=2)
qc = me.build()
perm = walsh.classical_permutation(qc)
rng = np.random.default_rng(0)
probe = rng.integers(0, 1 << me.n_qubits, size=200)
check("replay == classical_permutation",
      all(replay(qc, int(y)) == int(perm[y]) for y in probe),
      f"{probe.size} random basis states, {me.n_qubits} qubits")
del perm

# ---------------------------------------------------------------------------
print("\n[B] lookup tables are the intended classical tables")
me = WindowedModExp(N=7, a=3, n_exp=4, w=2)
ok = True
for k, win in enumerate(me.windows()):
    base = me.window_base(k)
    ok &= base == pow(me.a, 1 << (k * me.w), me.N)
    ok &= all(pow(base, j, me.N) == pow(me.a, j << (k * me.w), me.N)
              for j in range(1 << len(win)))
check("T_k[j] == a^(j*2^(kw)) mod N", ok, f"N=7 a=3 w=2, {len(me.windows())} windows")
check("tail windows are where a^(2^(kw)) == 1",
      [WindowedModExp(N=5, a=4, n_exp=8, w=2).is_tail_window(k) for k in range(4)]
      == [False, True, True, True],
      "N=5 a=4 (r=2, alpha=1): window 0 live, rest tail")

# ---------------------------------------------------------------------------
print("\n[C] the circuits compute a^e mod N, scratch clean")
for cls in (WindowedModExp, SelectModExp):
    for (N, a) in ((5, 4), (5, 2), (7, 3), (7, 6), (9, 2)):
        for w in (1, 2, 3):
            m = cls(N=N, a=a, n_exp=2 * w, w=w)
            check(f"{cls.__name__} N={N} a={a} w={w}", verify_modexp(m),
                  f"{m.n_qubits} qubits, all {1 << m.n_exp} exponents")

print("\n[C2] independent route: statevector, and vs ToffoliModExp")
# Statevector at the smallest reachable size -- shares no code with replay().
m = SelectModExp(N=5, a=4, n_exp=1, w=1)
ok = True
for e in range(2):
    qc = Circuit(m.n_qubits)
    if e:
        qc.x(m.exp[0])
    qc.extend(m.build())
    j, amp = sv.peak(sv.run(qc))
    ok &= sv.read_register(j, m.x) == pow(m.a, e, m.N) and abs(amp - 1) < 1e-8
check("statevector confirms SelectModExp N=5 a=4 w=1", ok, f"{m.n_qubits} qubits")

# The two constructions are structurally unrelated; agreeing on the valid
# subspace is a real cross-check of the arithmetic, not of the replay engine.
ok = True
for (N, a) in ((5, 4), (7, 3)):
    ref = ToffoliModExp(N=N, a=a, n_exp=2)
    for cls in (WindowedModExp, SelectModExp):
        m = cls(N=N, a=a, n_exp=2, w=2)
        for e in range(4):
            y0 = sum(1 << m.exp[i] for i in range(2) if (e >> i) & 1)
            z0 = sum(1 << ref.exp[i] for i in range(2) if (e >> i) & 1)
            ok &= read(replay(m.build(), y0), m.x) == \
                read(replay(ref.build(), z0), ref.x)
check("windowed constructions agree with ToffoliModExp on the valid subspace",
      ok, "N=5 a=4 and N=7 a=3, all 4 exponents")


# ---------------------------------------------------------------------------
def tail_block_perm(me, k=1):
    """Permutation of ONE window block, built with window k's multiplier but
    placed on the first window's qubits so the register stays small."""
    win = me.windows()[0]
    return walsh.classical_permutation(me.window_block(win, k)), win


def window_subperms(perm, win):
    """The block's action for each fixed window value, projected onto the
    other qubits. Valid because a block never changes the window register."""
    idn = np.arange(perm.size, dtype=np.int64)
    wval = sum(((idn >> q) & 1) << i for i, q in enumerate(win))
    out = []
    for j in range(1 << len(win)):
        sel = np.nonzero(wval == j)[0]
        out.append(np.searchsorted(sel, perm[sel]))
    return out


print("\n[D] identity-tail block: nontrivial, and an involution (the C29 test)")
for cls in (WindowedModExp, SelectModExp):
    m = cls(N=5, a=4, n_exp=2, w=2)          # r=2 => window 1 is a tail window
    P, win = tail_block_perm(m, k=1)
    idn = np.arange(P.size, dtype=np.int64)
    moved = int((P != idn).sum())
    check(f"{cls.__name__}: tail block != id (vacuity)", moved > 0,
          f"{moved} of {P.size} states moved")
    check(f"{cls.__name__}: tail block V^2 == id", bool((P[P] == idn).all()),
          "exhaustive over the full state space")
    del P

print("\n[E] what the tail block READS -- where the two designs diverge")
m = WindowedModExp(N=5, a=4, n_exp=2, w=2)
P, win = tail_block_perm(m, k=1)
idn = np.arange(P.size, dtype=np.int64)
indep = all(bool((P[idn ^ (1 << q)] == (P ^ (1 << q))).all()) for q in win)
check("lookup design: tail block is INDEPENDENT of the window register", indep,
      "all table entries equal 1, so the QROM permutation does not read j")
subs = window_subperms(P, win)
check("lookup design: V_j identical for every j including j=0",
      all(bool((s == subs[0]).all()) for s in subs), "4 window values")
del P, subs

m = SelectModExp(N=5, a=4, n_exp=2, w=2)
P, win = tail_block_perm(m, k=1)
idn = np.arange(P.size, dtype=np.int64)
indep = all(bool((P[idn ^ (1 << q)] == (P ^ (1 << q))).all()) for q in win)
check("select design: tail block DOES read the window register", not indep,
      "must-fail control for [E]: the designs must differ here")
subs = window_subperms(P, win)
check("select design: V_j identical for all j != 0 but different at j = 0",
      all(bool((s == subs[1]).all()) for s in subs[1:])
      and not bool((subs[0] == subs[1]).all()),
      "i.e. the block is gated on OR(window), a nonlinear function")
del P, subs

print("\n[F] control: a LIVE window block must read its window register")
m = WindowedModExp(N=7, a=3, n_exp=2, w=2)   # r=6, beta=3: no window is a tail
P, win = tail_block_perm(m, k=0)
idn = np.arange(P.size, dtype=np.int64)
indep = all(bool((P[idn ^ (1 << q)] == (P ^ (1 << q))).all()) for q in win)
check("lookup design, live window: block reads the window", not indep,
      "otherwise [E]'s independence result would be vacuous")
del P

# ---------------------------------------------------------------------------
print(f"\n{'ALL TESTS PASSED' if not FAILED else 'FAILURES: ' + ', '.join(FAILED)}")
raise SystemExit(1 if FAILED else 0)
