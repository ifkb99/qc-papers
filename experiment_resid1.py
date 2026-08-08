"""TODO 11.4/11.5 -- does the 0.716 residue persist, stack, and travel?

experiment_affstruct settled 11.1: the broken variants (density > 1/2) can
carry NO linear or affine structure, by counting. But 0.716 is still wildly
non-generic: a random function at q=15 has ~0.4% exact Walsh zeros (exact
balance of g^ell_z has probability ~ sqrt(2/(pi 2^q))), these variants have
~28%. Three questions decide what the residue IS before we anatomise it:

  (1) Does C15 constancy survive the breaking?  (beta=1 branch)
  (2) Does the excess-zero fraction persist as the circuit grows, or wash
      out?  (beta>1 branch, the C7-style trajectory)
  (3) Does stacking MORE independent nonlinear wraps push density to 1?

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  (derived, must hold) v4 with N=7 a=6 (r=2, beta=1): support size is
      CONSTANT across n_exp = 2,3,4. Reason: the wrap sits inside
      cc_add_mod, so the a=1 block becomes A'^-1 . S . A' with A' the wrapped
      multiply -- still a conjugate of the involutive cswap layer S, so
      V'^2 = id and theorem C29 applies verbatim. Breaking the linear
      structure and the 2-adic constancy are INDEPENDENT effects; the ~51%
      penalty does not forfeit the free exponent register.
      (If this fails, C29's reading of the construction is wrong -- that
      would be the bigger news.)

  P2  (open, the discriminating measurement) v4 with N=7 a=3 (r=6, beta>1):
      density at n_exp = 2..6. H-wash: density climbs toward 1.000 the way
      the baseline climbs toward 0.500 (residue is finite-size). H-persist:
      density plateaus near 0.72 (residue is intrinsic structure). No
      commitment; this is the fork in the road.

  P3  (open) stacking wraps -- v45 = both toffoli(msb,t0,anc) and
      toffoli(msb,t0,t1); v3x adds toffoli(msb,t1,anc). All must still
      compute a^e mod N (wrap controls are 0 on the valid subspace).
      If density saturates at ~0.72 under stacking, the residue does not
      live in the reduction's linear coupling at all.

  P4  (closes 11.5) The exact zero fractions across instances and widths.
      Already known: 23464/32768 (N=5) vs 23488/32768 (N=7) -- unequal, so
      "0.716" is not a universal rational. P2's sweep settles whether it is
      even width-stable for a fixed instance.
"""
from __future__ import annotations
import numpy as np
import time
from toffoli_arith import ToffoliModExp
from circuits import Circuit
import statevec as sv
import walsh

TOL = 1e-12


class Variant(ToffoliModExp):
    MODES = {"v0": (), "v4": ("A",), "v5": ("B",), "v45": ("A", "B"),
             "v3x": ("A", "B", "C")}

    def __init__(self, *a, mode="v0", **k):
        self.mode = mode
        super().__init__(*a, **k)

    def _wrap(self, qc):
        msb = self.b[self.m - 1]
        for tag in self.MODES[self.mode]:
            if tag == "A":
                qc.toffoli(msb, self.t[0], self.anc)
            elif tag == "B":
                qc.toffoli(msb, self.t[0], self.t[1])
            elif tag == "C":
                qc.toffoli(msb, self.t[1], self.anc)

    def cc_add_mod(self, qc, c1, c2, c):
        self._wrap(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._wrap(qc)


def support(me):
    qc = me.build()
    c = walsh.pullback_coefficients(qc, me.x[0])
    S = int((np.abs(c) > TOL).sum())
    return S, qc.n


def correct(me):
    for e in range(1 << me.n_exp):
        qc = Circuit(me.n_qubits)
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1:
                qc.x(eq)
        qc.extend(me.build())
        j, amp = sv.peak(sv.run(qc))
        if sv.read_register(j, me.x) != pow(me.a, e, me.N) or abs(amp - 1) > 1e-8:
            return False
    return True


print("=" * 80)
print("P3-pre  Correctness gate for the stacked variants (must-pass)")
print("=" * 80)
for mode in ("v45", "v3x"):
    me = Variant(N=5, a=2, n_exp=2, mode=mode)
    print(f"  {mode}: computes a^e mod N correctly = {correct(me)}")

print()
print("=" * 80)
print("P1  beta=1 branch: does constancy survive the breaking?  (derived: yes)")
print("=" * 80)
print(f"  {'mode':>5} {'n_exp':>6} {'q':>4} {'|supp|':>9} {'density':>9}  note")
for mode in ("v0", "v4"):
    base = None
    for ne in (2, 3, 4):
        me = Variant(N=7, a=6, n_exp=ne, mode=mode)
        S, q = support(me)
        tag = "-" if base is None else ("SAME" if S == base else f"{S-base:+d}")
        if base is None:
            base = S
        print(f"  {mode:>5} {ne:6d} {q:4d} {S:9d} {S/(1<<q):9.6f}  {tag}",
              flush=True)

print()
print("=" * 80)
print("P2  beta>1 branch: trajectory of the broken density (the fork)")
print("=" * 80)
print(f"  {'mode':>5} {'n_exp':>6} {'q':>4} {'|supp|':>9} {'density':>9} "
      f"{'zeros%':>8} {'time':>7}")
for mode in ("v0", "v4"):
    for ne in (2, 3, 4, 5, 6):
        t0 = time.time()
        me = Variant(N=7, a=3, n_exp=ne, mode=mode)
        S, q = support(me)
        el = time.time() - t0
        print(f"  {mode:>5} {ne:6d} {q:4d} {S:9d} {S/(1<<q):9.6f} "
              f"{100*(1-S/(1<<q)):7.2f}% {el:6.1f}s", flush=True)

print()
print("=" * 80)
print("P3/P4  Stacking, and the exact fractions across instances")
print("=" * 80)
print(f"  {'mode':>5} {'N':>3} {'a':>3} {'q':>4} {'|supp|':>9} {'density':>9}  "
      f"exact fraction")
for mode in ("v0", "v4", "v5", "v45", "v3x"):
    for N_, a_ in ((5, 2), (7, 6), (7, 3)):
        me = Variant(N=N_, a=a_, n_exp=2, mode=mode)
        S, q = support(me)
        print(f"  {mode:>5} {N_:3d} {a_:3d} {q:4d} {S:9d} {S/(1<<q):9.6f}  "
              f"{S}/{1<<q}", flush=True)
