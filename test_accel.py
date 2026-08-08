"""Gate for the CUDA backend: the GPU path must agree with the CPU reference.

Suite 9. Skips cleanly (exit 0) when no CUDA device is present, so the gate
still runs on machines without a card.

The standard here is deliberately strict. A permutation replay is integer
arithmetic with no rounding, so *exact* equality is required and anything less
is a bug. The transform is floating point, so allclose is the right test, but
we also require the *support sets* to match exactly, since that -- not the
coefficient values -- is the quantity every headline number is built on.
"""
from __future__ import annotations
import numpy as np

import accel
import walsh
from circuits import ripple_adder
from toffoli_arith import ToffoliModExp
from windowed_arith import WindowedModExp

FAILED = []


def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra else ''}")
    if not cond:
        FAILED.append(name)


if not accel.HAVE_GPU:
    print("\n[skip] no CUDA device -- accel backend not exercised")
    raise SystemExit(0)

print(f"\n[device] {accel.device_info()}")

# ---------------------------------------------------------------------------
print("\n[A] permutation replay: EXACT equality with the CPU reference")
cases = []
for nb in (3, 4):
    add, lay = ripple_adder(nb)
    cases.append((f"{nb}-bit adder", add, lay["b"][0]))
for N, a, ne in ((5, 2, 2), (7, 3, 2), (15, 7, 1)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    cases.append((f"modexp N={N} a={a} n_exp={ne}", me.build(), me.x[0]))
me = WindowedModExp(N=5, a=4, n_exp=2, w=2)
cases.append(("windowed modexp N=5 a=4", me.build(), me.x[0]))

for label, qc, tq in cases:
    pc = walsh.classical_permutation(qc)
    pg = accel.classical_permutation(qc)
    check(f"{label} (q={qc.n})", np.array_equal(pc, pg),
          f"{pc.size:,} entries, exact")
    del pc, pg

# The replay computes in int32 for n <= 30 (halves the working set, which is
# what makes q = 30 fit). Bit twiddling is only safe if every image really is
# below 2^31, and if the high-index gates survive the narrower type -- so check
# both directly rather than trusting the argument.
qc_hi, tq_hi = cases[2][1], cases[2][2]
check("replay index dtype is int32 for n <= 30",
      accel._index_dtype(30) == np.int32 and accel._index_dtype(31) == np.int64)
_pg = accel.classical_permutation(qc_hi)
check("int32 replay stays in range", int(_pg.max()) < (1 << qc_hi.n),
      f"max image {int(_pg.max()):,} < 2^{qc_hi.n}")
check("int32 replay is a permutation",
      np.array_equal(np.sort(_pg), np.arange(1 << qc_hi.n)))
del _pg

print("\n[B] FWHT: allclose on values, EXACT on the support set")
rng = np.random.default_rng(0)
for n in (12, 16, 20):
    v = rng.integers(0, 2, size=1 << n).astype(np.float64)
    wc, wg = walsh.wht(v), accel.wht(v)
    err = float(np.max(np.abs(wc - wg)))
    check(f"FWHT 2^{n} values", err < 1e-9, f"maxerr {err:.2e}")
    check(f"FWHT 2^{n} support set",
          np.array_equal(np.nonzero(np.abs(wc) > 1e-9)[0],
                         np.nonzero(np.abs(wg) > 1e-9)[0]))
    del v, wc, wg

print("\n[C] end-to-end pullback_support vs the CPU pipeline")
for label, qc, tq in cases:
    zc = np.nonzero(np.abs(walsh.pullback_coefficients(qc, tq)) > 1e-12)[0]
    zg = accel.pullback_support(qc, tq, tol=1e-12)
    check(f"{label}: support sets identical",
          np.array_equal(zc.astype(np.int64), zg),
          f"{zc.size:,} terms")
    del zc, zg

print("\n[F] exact integer path: identical support, no tolerance involved")
for label, qc, tq in cases:
    if qc.n > 30:
        continue
    ref = np.nonzero(np.abs(walsh.pullback_coefficients(qc, tq)) > 1e-12)[0]
    ex = accel.pullback_support_exact(qc, tq)
    check(f"{label}: exact-int support == float64 support",
          np.array_equal(ref.astype(np.int64), ex), f"{ref.size:,} terms")
    del ref, ex

# the property the exact path relies on: the FWHT of +/-1 data is integral
rng2 = np.random.default_rng(1)
chi = np.where(rng2.integers(0, 2, size=1 << 16) == 1, -1.0, 1.0)
w = walsh.wht(chi)
check("FWHT of +/-1 data is exactly integer-valued",
      np.array_equal(w, np.round(w)),
      f"max|W| = {int(np.abs(w).max()):,}, int32 bound 2^31")
check("int32 FWHT reproduces float64 bit-for-bit",
      np.array_equal(accel.wht_exact(chi.astype(np.int32)).astype(np.int64),
                     w.astype(np.int64)))

print("\n[G] pullback_stats: aggregates match counting the support by hand")
me_s = ToffoliModExp(N=5, a=2, n_exp=1)
qc_s, tq_s = me_s.build(), me_s.x[0]
w_lin = (1 << me_s.b[me_s.m - 1]) | (1 << me_s.anc)      # the C30 structure
w_bad = 1 << me_s.x[1]                                   # not a structure
zs_s = accel.pullback_support_exact(qc_s, tq_s)
st = accel.pullback_stats(qc_s, tq_s, masks=(w_lin, w_bad))
check("stats count == |support|", st["count"] == zs_s.size, f"{st['count']:,}")
check("stats density == count / 2^n",
      abs(st["density"] - zs_s.size / (1 << qc_s.n)) < 1e-15)
for w in (w_lin, w_bad):
    ref = sum((int(z) & w).bit_count() & 1 for z in zs_s.tolist())
    check(f"stats odd[{w}] matches host popcount", st["odd"][w] == ref,
          f"{st['odd'][w]:,} vs {ref:,}")
check("C30 structure has odd-count 0 (pinned), control does not",
      st["odd"][w_lin] == 0 and st["odd"][w_bad] > 0)
del zs_s

print("\n[D] control: a non-permutation circuit must be REJECTED, not silently wrong")
from circuits import Circuit
bad = Circuit(4)
bad.h(0)
rejected = False
try:
    accel.classical_permutation(bad)
except ValueError:
    rejected = True
check("GPU path rejects a non-classical circuit", rejected,
      "matches walsh.classical_permutation's contract")

print("\n[E] control: oversized register must raise rather than thrash")
class _Fake:
    n = accel.MAX_QUBITS + 1
    logical = []
    def is_classical(self): return True
raised = False
try:
    accel.classical_permutation(_Fake())
except MemoryError:
    raised = True
check("GPU path refuses n > MAX_QUBITS", raised, f"MAX_QUBITS={accel.MAX_QUBITS}")

# ---------------------------------------------------------------------------
print(f"\n{'ALL TESTS PASSED' if not FAILED else 'FAILURES: ' + ', '.join(FAILED)}")
raise SystemExit(1 if FAILED else 0)
