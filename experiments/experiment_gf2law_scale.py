"""C40/C41 at scale: the conditional-structure law on up to 2^30 points.

Answers the "toy sizes" objection to PAPER_A SS6.2. The law is PROVED, so scale
demonstrates rather than establishes it -- but n = 14 invites the question and
n = 30 (1.07e9 points) closes it cheaply.

Memory-lean by necessity: int8 plant built in chunks, an IN-PLACE FWHT (one
half-size temporary instead of two full copies), and a chunked support/E scan.
The in-place transform is gated against walsh.wht on a small case before use --
a hand-rolled replacement for a verified primitive is exactly the kind of thing
that silently corrupts a headline number.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Zero support elements inside E at every n, in every consistent case.
      This is the decisive claim and it is THRESHOLD-ROBUST: E-membership is
      combinatorial, not a magnitude test.
  P2  Densities sit slightly BELOW their caps, and the deficit is a
      MEASUREMENT artifact rather than structure: with an absolute cut of 1e-6
      and random-function coefficients ~N(0, 2^-n/2), the expected fraction
      erased is erf(tol / (sigma*sqrt2)) -- 0.7% at n=26, 1.3% at n=28, 2.6% at
      n=30. The no-cap (odd-dependency) row is the clean test since its true
      density is 1: it should read ~0.9935 / 0.9869 / 0.9739.
  C1  MUST-FAIL: the odd-dependency case must show NO cap at any n.

Run:  PYTHONPATH=. uv run python -m experiments.experiment_gf2law_scale 30
      (~80 s per case and ~15 GB at n=30; drop to 26 for a 5 s version.)
"""
import sys, time
import numpy as np
import walsh
from lab import rank_kernel

N = int(sys.argv[1]); CH = 1 << 24
RNG = np.random.default_rng(20260808)

def wht_inplace(a):
    n = a.size; h = 1
    while h < n:
        a = a.reshape(-1, 2, h)
        t = a[:, 0, :] - a[:, 1, :]
        a[:, 0, :] += a[:, 1, :]
        a[:, 1, :] = t
        del t
        a = a.reshape(-1)
        h *= 2
    return a

# gate the in-place FWHT against the reference on a small case
_v = RNG.integers(0, 2, size=1 << 12).astype(np.float64)
assert np.allclose(wht_inplace(_v.copy()), walsh.wht(_v)), "in-place FWHT mismatch"

def span_dim(ws, n): return rank_kernel(list(dict.fromkeys(ws)), n)[0]
def consistent(ws, n):
    u = list(dict.fromkeys(ws))
    return rank_kernel([(w << 1) | 1 for w in u], n + 1)[0] == rank_kernel(u, n)[0]

def run(name, cond_bits, ws):
    t0 = time.perf_counter()
    m = len(cond_bits)
    bits = RNG.integers(0, 2, size=1 << N, dtype=np.int8)
    chi = np.empty(1 << N, dtype=np.float64)
    leads = [w.bit_length() - 1 for w in ws]
    for lo in range(0, 1 << N, CH):                      # chunked plant
        y = np.arange(lo, lo + CH, dtype=np.int64)
        u = np.zeros(CH, dtype=np.int64)
        for i, c in enumerate(cond_bits):
            u |= ((y >> c) & 1) << i
        rep = y.copy()
        for uu in range(1 << m):
            sel = (u == uu)
            flip = ((y >> leads[uu]) & 1).astype(bool) & sel
            rep[flip] ^= np.int64(ws[uu])
        chi[lo:lo + CH] = np.where(bits[rep] == 1, -1.0, 1.0)
        del y, u, rep
    del bits
    c = wht_inplace(chi); c /= (1 << N)

    viol = 0; Esz = 0; nnz = 0
    for lo in range(0, 1 << N, CH):                      # chunked support/E scan
        nz = np.abs(c[lo:lo + CH]) > 1e-6
        nnz += int(nz.sum())
        z = np.arange(lo, lo + CH, dtype=np.int64)
        keep = np.ones(CH, dtype=bool)
        for w in dict.fromkeys(ws):
            par = np.zeros(CH, dtype=np.int8); ww = int(w)
            while ww:
                b = ww & -ww
                par ^= ((z >> (b.bit_length() - 1)) & 1).astype(np.int8)
                ww ^= b
            keep &= (par == 1)
        Esz += int(keep.sum()); viol += int((nz & keep).sum())
        del nz, z, keep
    del c
    d = span_dim(ws, N); ok = consistent(ws, N)
    cap = 1 - 2.0 ** (-d) if ok else 1.0
    print(f"  {name:26s} d={d} consistent={str(ok):5s} |E|={Esz:>14,} cap={cap:.4f} "
          f"density={nnz/(1<<N):.4f} violations={viol}  [{time.perf_counter()-t0:.0f}s]",
          flush=True)

print(f"== C40/C41 at n={N} (2^{N} = {1<<N:,} points) ==")
run("d=1 (C30 form)",        [],     [0b110])
run("d=2, 2 cells indep",    [0],    [0b110, 0b1010])
run("d=3, 4 cells indep",    [0,1],  [0b100, 0b1000, 0b10000, 0b11100])
w0, w1, w2 = 0b100, 0b1000, 0b10000
run("ODD dependency (C41)",  [0,1],  [w0, w1, w0 ^ w1, w0])
run("EVEN dependency (C41)", [0,1],  [w0, w1, w2, w0 ^ w1 ^ w2])

import math
print()
print("P2 check -- is the sub-cap deficit a thresholding artifact?")
sig = 2 ** (-N / 2); tol = 1e-6
print(f"  n={N}: sigma={sig:.2e}, tol={tol:.0e} -> expected erased fraction "
      f"{math.erf(tol / (sig * math.sqrt(2))):.4f}, so a TRUE density of 1.0 "
      f"should read ~{1 - math.erf(tol / (sig * math.sqrt(2))):.4f}")
print("  Compare against the ODD-dependency row above, whose true density is 1.")
