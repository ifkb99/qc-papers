"""C8: does Walsh-Hadamard sparsity predict the PPS term count for permutation
circuits, and is the prediction exact?

Key numerical fact used throughout: for an n-qubit permutation circuit the
coefficients c_z are exact multiples of 2^-n. For n=15 the smallest possible
nonzero coefficient is 2^-15 = 3.05e-5. Anything PPS reports between its
1e-13 working floor and that value is therefore accumulated floating-point
noise, not signal. That gives a clean cut for separating the two.
"""
from __future__ import annotations
import numpy as np, time
from toffoli_arith import ToffoliModExp
from circuits import ripple_adder
from pps import propagate
import walsh

SIGNAL_CUT = 1e-6          # noise floor << this << 2^-15


def compare(qc, target, label):
    n = qc.n
    quantum = 1.0 / (1 << n)

    t0 = time.time()
    c_exact = walsh.pullback_coefficients(qc, target)
    t_walsh = time.time() - t0

    t0 = time.time()
    r = propagate(qc, {(0, 1 << target): 1.0}, delta=0.0, max_terms=2_000_000)
    t_pps = time.time() - t0

    # rebuild the PPS coefficient vector, and check it stayed Z-type
    pps = np.zeros(1 << n)
    non_z = 0
    for (x, z), v in r.final_terms.items():
        if x != 0:
            non_z += 1
        else:
            pps[z] = v

    exact_nz = np.abs(c_exact) > quantum / 2
    pps_sig = np.abs(pps) > SIGNAL_CUT
    err = np.abs(pps - c_exact)

    print(f"  {label}")
    print(f"    qubits={n}  gates={len(qc.gates)}  walsh={t_walsh:.2f}s  pps={t_pps:.1f}s")
    print(f"    coefficients are multiples of 2^-n = {quantum:.3e}")
    print(f"    non-Z-type terms in PPS output      : {non_z}")
    print(f"    Walsh sparsity (exact)              : {int(exact_nz.sum())}")
    print(f"    PPS terms above {SIGNAL_CUT:.0e}            : {int(pps_sig.sum())}")
    print(f"    PPS raw term count (floor 1e-13)    : {r.n_terms[-1]}")
    print(f"    -> signal counts match              : {int(pps_sig.sum()) == int(exact_nz.sum())}")
    print(f"    -> supports identical               : {np.array_equal(pps_sig, exact_nz)}")
    print(f"    max |PPS - Walsh| over all 2^n coeffs: {err.max():.3e}")
    print(f"    max error on true-nonzero coeffs     : {err[exact_nz].max():.3e}")
    print(f"    max spurious |c| on true-zero coeffs : {np.abs(pps[~exact_nz]).max():.3e}")
    return int(exact_nz.sum()), int(pps_sig.sum())


print("=" * 78)
print("C8 -- Walsh sparsity vs exact PPS, Toffoli-compiled circuits")
print("=" * 78)
results = []
for N, a, ne in ((5, 2, 2), (7, 3, 2)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    for oi in (0, 1):
        results.append(compare(qc, me.x[oi], f"modexp N={N} a={a} obs=Z_x{oi}"))
        print()

add, lay = ripple_adder(4)
for i in (0, 2):
    results.append(compare(add, lay["b"][i], f"ripple adder 4-bit obs=Z_b{i}"))
    print()

print("=" * 78)
print(f"SUMMARY: {sum(1 for w, p in results if w == p)}/{len(results)} instances "
      f"where Walsh sparsity exactly predicts the PPS signal-term count")
print("=" * 78)
