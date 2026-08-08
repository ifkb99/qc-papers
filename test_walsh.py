"""Verify the logical trace and Walsh machinery before using them for C8."""
from __future__ import annotations
import numpy as np
from circuits import Circuit, ripple_adder
from toffoli_arith import ToffoliModExp
from modexp import ModExp
import statevec as sv
import walsh

FAILED = []
def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra else ''}")
    if not cond: FAILED.append(name)


print("\n[A] logical trace records atomically")
qc = Circuit(3); qc.toffoli(0, 1, 2)
check("toffoli logs one op, not its decomposition",
      qc.logical == [("toffoli", 0, 1, 2)], f"{qc.logical[:3]}")
qc = Circuit(2); qc.cnot(0, 1)
check("cnot logs one op", qc.logical == [("cnot", 0, 1)], f"{qc.logical}")
qc = Circuit(2); qc.h(0)
check("h logs nonclassical", qc.logical and qc.logical[0][0] == "nonclassical")
qc = Circuit(3); qc.swap(0, 1); qc.cswap(0, 1, 2)
check("swap/cswap decompose to logged primitives", qc.is_classical(),
      f"{[o[0] for o in qc.logical]}")

print("\n[B] classical replay == statevector simulation")
for name, qc in (("3-bit adder", ripple_adder(3)[0]),
                 ("toffoli modexp N=5 a=2", ToffoliModExp(5, 2, n_exp=2).build())):
    if not qc.is_classical():
        check(f"{name} is classical", False); continue
    perm = walsh.classical_permutation(qc)
    rng = np.random.default_rng(0)
    probe = rng.choice(1 << qc.n, size=min(12, 1 << qc.n), replace=False)
    ok = True
    for y in probe:
        psi = sv.run(qc, sv.basis(qc.n, int(y)))
        j, amp = sv.peak(psi)
        if j != perm[y] or abs(amp - 1) > 1e-8:
            ok = False; break
    check(f"{name}: replay matches statevec on {len(probe)} basis states", ok)

print("\n[C] Fourier compilation is correctly rejected")
fou = ModExp(5, 2, n_exp=2).build()
check("Fourier modexp is_classical() == False", not fou.is_classical())

print("\n[D] Walsh transform sanity")
qc = Circuit(3); qc.cnot(0, 1)                    # bit1 -> b1 XOR b0 : linear
c = walsh.pullback_coefficients(qc, 1)
nz = np.nonzero(np.abs(c) > 1e-12)[0]
check("CNOT: Z1 pulls back to a single term Z0Z1",
      len(nz) == 1 and nz[0] == 0b011 and abs(c[nz[0]] - 1) < 1e-12,
      f"nz={nz} c={c[nz]}")

qc = Circuit(3); qc.toffoli(0, 1, 2)              # bit2 -> b2 XOR (b0 AND b1)
c = walsh.pullback_coefficients(qc, 2)
nz = set(np.nonzero(np.abs(c) > 1e-12)[0].tolist())
check("Toffoli: Z2 pulls back to 4 terms (AND has 4 Walsh coeffs)",
      len(nz) == 4, f"{sorted(nz)}")
check("Toffoli is not affine", not walsh.is_affine(qc, 2))

print("\n[E] adder output bit is XOR-affine (explains the F1 collapse)")
add, lay = ripple_adder(3)
s = walsh.walsh_sparsity(add, lay["b"][0])
check("adder low bit: Walsh sparsity == 1 (affine)", s == 1, f"got {s}")
check("adder low bit: is_affine", walsh.is_affine(add, lay["b"][0]))
print("    (sparsity 1 <=> affine; the support point has weight 3 because the"
      "\n     bit is a0 XOR b0 XOR c0 -- weight is not algebraic degree)")
print("    higher bits: " + ", ".join(
    f"b{i}: sparsity={walsh.walsh_sparsity(add, lay['b'][i])}"
    f" wt={walsh.walsh_degree(add, lay['b'][i])}" for i in range(3)))

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
