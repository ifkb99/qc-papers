"""C15.3 -- is the Pauli support EXACTLY constant when the order r is a power
of two, as the exponent register (and hence the whole Hilbert space) grows?

If so, PPS cost for such instances is independent of the precision of the
period-finding register, which is a sharp and checkable statement.
"""
from __future__ import annotations
import numpy as np, time
from toffoli_arith import ToffoliModExp
import walsh

print("C15.3  Is sparsity EXACTLY constant when r is a power of two?")
print("N=7, a=6 (r=2). Exponent register grows; Hilbert space grows as 2^n_exp.\n")
print(f"  {'n_exp':>6} {'qubits':>7} {'2^q':>10} {'sparsity':>10} {'density':>9} "
      f"{'vs n_exp=2':>11} {'time':>7}")
base = None
for ne in (2, 4, 6, 8, 10):
    me = ToffoliModExp(N=7, a=6, n_exp=ne)
    qc = me.build()
    t0 = time.time()
    c = walsh.pullback_coefficients(qc, me.x[0])
    el = time.time() - t0
    sp = int((np.abs(c) > 1e-12).sum())
    dim = 1 << me.n_qubits
    if base is None:
        base = sp
    tag = "SAME" if sp == base else f"{sp - base:+d}"
    print(f"  {ne:6d} {me.n_qubits:7d} {dim:10d} {sp:10d} {sp/dim:9.6f} "
          f"{tag:>11} {el:6.1f}s", flush=True)
    del c

print("\n  CONTROL: a=3 (r=6, has an odd factor) at the same sizes")
print(f"  {'n_exp':>6} {'qubits':>7} {'sparsity':>10} {'density':>9} {'growth':>9}")
prev = None
for ne in (2, 4, 6, 8, 10):
    me = ToffoliModExp(N=7, a=3, n_exp=ne)
    qc = me.build()
    c = walsh.pullback_coefficients(qc, me.x[0])
    sp = int((np.abs(c) > 1e-12).sum())
    dim = 1 << me.n_qubits
    print(f"  {ne:6d} {me.n_qubits:7d} {sp:10d} {sp/dim:9.6f} "
          f"{(f'{sp/prev:.2f}x' if prev else '-'):>9}", flush=True)
    prev = sp
    del c

print("""
  Each +2 in n_exp multiplies the Hilbert space by 4. A power-of-two order
  should hold sparsity fixed (density falling 4x each row); an order with an
  odd factor should grow ~4x (density flat near 1/2).""")
