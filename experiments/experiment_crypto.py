"""TODO step 5 -- importing cryptanalytic results into PPS resource estimation.

C12 established a qualitative bridge (PPS-hard <-> resists linear
approximation). This makes it a usable inequality and tests it on functions
whose nonlinearity is PUBLISHED rather than measured here.

THE BOUND. With normalised Walsh coefficients (sum c_z^2 = 1 by Parseval),

    NL = 2^(n-1) (1 - max_z |c_z|)          =>   max|c| = 1 - NL/2^(n-1)
    S * max|c|^2  >=  sum c_z^2  =  1       =>   S >= 1 / max|c|^2

hence

    **S  >=  (1 - NL/2^(n-1))^(-2)**

Since (Paper A) S is exactly the Pauli support PPS carries for a reversible
circuit computing that bit, this converts any published nonlinearity into a
LOWER BOUND on Pauli-path simulation cost, for every circuit computing the
function -- compilation-independent.

Test functions with known values:
  * affine            NL = 0                      -> bound S >= 1   (tight)
  * inner product     bent, NL = 2^(n-1)-2^(n/2-1) -> bound S >= 2^n (tight)
  * AES S-box         NL = 112 at n = 8            -> bound S >= 64
  * modexp bits, random f                          -> for comparison
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm
import walsh

TOL = 1e-12


def spectrum_of(truth: np.ndarray, n: int):
    """Normalised Walsh coefficients of (-1)^truth."""
    return walsh.wht(np.where(truth == 1, -1.0, 1.0)) / (1 << n)


def stats(c: np.ndarray, n: int):
    a = np.abs(c)
    S = int((a > TOL).sum())
    mx = float(a.max())
    NL = 2 ** (n - 1) * (1 - mx)
    bound = (1.0 / mx ** 2) if mx > 0 else float("inf")
    NL_bent = 2 ** (n - 1) - 2 ** (n / 2 - 1)
    return S, mx, NL, bound, NL_bent


# ---------------------------------------------------------------- GF(2^8)
def gf_mul(a, b):
    p = 0
    for _ in range(8):
        if b & 1:
            p ^= a
        hi = a & 0x80
        a = (a << 1) & 0xFF
        if hi:
            a ^= 0x1B
        b >>= 1
    return p


INV = [0] * 256
for x in range(1, 256):
    for y in range(1, 256):
        if gf_mul(x, y) == 1:
            INV[x] = y
            break


def aes_sbox_table():
    """AES S-box = affine( inverse(x) )."""
    out = []
    for x in range(256):
        b = INV[x]
        s = 0
        for i in range(8):
            bit = ((b >> i) & 1) ^ ((b >> ((i + 4) % 8)) & 1) ^ \
                  ((b >> ((i + 5) % 8)) & 1) ^ ((b >> ((i + 6) % 8)) & 1) ^ \
                  ((b >> ((i + 7) % 8)) & 1) ^ ((0x63 >> i) & 1)
            s |= bit << i
        out.append(s)
    return out


SBOX = aes_sbox_table()

print("=" * 84)
print("1. THE BOUND  S >= (1 - NL/2^(n-1))^(-2)  on functions with KNOWN NL")
print("=" * 84)
print(f"  {'function':>26} {'n':>3} {'S (actual)':>11} {'bound':>10} "
      f"{'holds':>6} {'tight':>6} {'NL':>8} {'NL/NLbent':>10}")

rows = []

# affine
n = 10
y = np.arange(1 << n)
aff = (np.bitwise_count(y & 0b1011) & 1)
c = spectrum_of(aff, n)
S, mx, NL, bd, NLb = stats(c, n)
print(f"  {'affine (x0^x1^x3)':>26} {n:3d} {S:11d} {bd:10.1f} "
      f"{str(S >= bd - 1e-9):>6} {str(abs(S - bd) < 1e-6):>6} {NL:8.0f} {NL/NLb:10.4f}")

# inner product = bent, on n = 2m variables
for m in (4, 6):
    n = 2 * m
    y = np.arange(1 << n)
    ip = np.bitwise_count((y & ((1 << m) - 1)) & (y >> m)) & 1
    c = spectrum_of(ip, n)
    S, mx, NL, bd, NLb = stats(c, n)
    print(f"  {f'inner product (bent) m={m}':>26} {n:3d} {S:11d} {bd:10.1f} "
          f"{str(S >= bd - 1e-9):>6} {str(abs(S - bd) < 1e-6):>6} {NL:8.0f} "
          f"{NL/NLb:10.4f}")

# AES S-box coordinate bits
n = 8
tt = np.array(SBOX, dtype=np.int64)
nls = []
for bit in range(8):
    c = spectrum_of((tt >> bit) & 1, n)
    S, mx, NL, bd, NLb = stats(c, n)
    nls.append(NL)
    if bit < 3:
        print(f"  {f'AES S-box bit {bit}':>26} {n:3d} {S:11d} {bd:10.1f} "
              f"{str(S >= bd - 1e-9):>6} {str(abs(S - bd) < 1e-6):>6} {NL:8.0f} "
              f"{NL/NLb:10.4f}")
# min over all nonzero linear combinations -- the published figure
worst = min(
    stats(spectrum_of(np.bitwise_count(tt & mask) & 1, n), n)[2]
    for mask in range(1, 256)
)
print(f"  {'AES min over combos':>26} {n:3d} {'-':>11} {'-':>10} {'-':>6} "
      f"{'-':>6} {worst:8.0f}   <- published value is 112")

# modexp + random for comparison
me = ToffoliModExp(N=5, a=2, n_exp=1)
qc = me.build()
c = walsh.pullback_coefficients(qc, me.x[0])
S, mx, NL, bd, NLb = stats(c, qc.n)
print(f"  {'modexp N=5 (Z_x0)':>26} {qc.n:3d} {S:11d} {bd:10.1f} "
      f"{str(S >= bd - 1e-9):>6} {str(abs(S - bd) < 1e-6):>6} {NL:8.0f} {NL/NLb:10.4f}")

rng = np.random.default_rng(0)
n = 12
c = spectrum_of(rng.integers(0, 2, size=1 << n), n)
S, mx, NL, bd, NLb = stats(c, n)
print(f"  {'random f':>26} {n:3d} {S:11d} {bd:10.1f} "
      f"{str(S >= bd - 1e-9):>6} {str(abs(S - bd) < 1e-6):>6} {NL:8.0f} {NL/NLb:10.4f}")

print("""
  'tight' means the bound is attained. It is tight exactly at the two extremes
  -- affine (one spike) and bent (flat spectrum) -- because those are the cases
  where Parseval is saturated by a single magnitude.""")

print()
print("=" * 84)
print("2. CLOSING THE LOOP: a reversible circuit whose output bit is BENT")
print("=" * 84)
print("  Toffoli(x_i, y_i, out) for each i computes out ^= <x,y>, which is bent.")
print("  Published fact: bent => flat spectrum => PPS must carry EVERY Z-string.\n")
print(f"  {'m':>3} {'qubits':>7} {'perm_pps |supp|':>16} {'predicted':>10} {'match':>7}")
for m in (3, 4, 5):
    nq = 2 * m + 1
    qc = Circuit(nq)
    out = 2 * m
    for i in range(m):
        qc.toffoli(i, m + i, out)
    r = propagate_perm(qc, 1 << out, delta=0.0)
    pred = 1 << (2 * m)          # z_out=1 forced; bent => all 2^(2m) others
    got = len(r.final_terms)
    print(f"  {m:3d} {nq:7d} {got:16d} {pred:10d} {str(got == pred):>7}", flush=True)

print("""
  The support is exactly 2^(2m): the output qubit is forced into every term
  (z_out = 1) and the bent structure fills the remaining 2m coordinates
  completely. No truncation can help -- every coefficient has the same
  magnitude, so there is no 'large' subset to keep.""")

print()
print("=" * 84)
print("3. WHAT THIS BUYS")
print("=" * 84)
print(f"""  For AES's S-box, nonlinearity 112 at n = 8 is a published constant. The
  bound then says any reversible circuit computing an S-box output bit forces
  PPS to carry at least {int((1-112/128)**-2)} Pauli terms, whatever the compilation --
  without running any simulation, and without even building the circuit.

  Read the other way: a reversible circuit is cheap for Pauli propagation only
  if its output bits admit a good linear approximation. Cryptographic design
  deliberately maximises nonlinearity, so cryptographic primitives are
  systematically the worst case for this simulation method. Decades of published
  nonlinearity bounds transfer directly, and the transfer is one line.""")
