"""TODO step 3, attempt 3 -- identify WHY the branch magnitudes are equal.

Established in attempt 2: for beta=1, every Walsh coefficient has exactly one of
its two branch partners nonzero, and the surviving magnitudes reproduce the
smaller circuit's spectrum bit-for-bit (max err 0.00e+00). So

    h^(z) = eps(z) * g^(z),   eps(z) in {+1, -1}

where g and h are the two branches (new exponent qubit off / on). Read off from
the spectrum: eps(z) = +1 when the surviving term sits on the e=0 side
(g^+h^ != 0, g^-h^ = 0 => h^ = g^), and -1 when it sits on the e=1 side.

The question is what eps is. If eps is a LINEAR CHARACTER,

    eps(z) = (-1)^<z, v>   for some fixed v,

then h(y) = g(y XOR v): the block acts, after conjugation, as a TRANSLATION of
the input. Translations preserve Walsh magnitudes exactly and multiply
coefficients by a character -- which is precisely the observed structure, and
would complete the mechanism.

Tests:
  E1  Is eps multiplicative?  eps(z1)eps(z2) = eps(z1 XOR z2) whenever all three
      are in the support. (Necessary for a character.)
  E2  Solve for v from the support and verify eps(z) = (-1)^<z,v> on ALL of it.
  E3  What is v?  Which qubits does it touch?
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


def eps_from_spectrum(N, a, ne):
    """Return (zs, eps) for the branch sign pattern at the top exponent qubit."""
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    c = walsh.pullback_coefficients(qc, me.x[0])
    top = me.exp[-1]
    half = 1 << top
    lo, hi = c[:half], c[half:]
    nz_lo, nz_hi = np.abs(lo) > TOL, np.abs(hi) > TOL
    assert not np.any(nz_lo & nz_hi), "branch symmetry violated"
    zs = np.nonzero(nz_lo | nz_hi)[0]
    eps = np.where(nz_lo[zs], 1, -1)
    return me, qc, zs, eps


def analyse(N, a, ne, label):
    me, qc, zs, eps = eps_from_spectrum(N, a, ne)
    n = me.exp[-1]                      # number of non-top qubits
    emap = dict(zip(zs.tolist(), eps.tolist()))
    print(f"  {label}: n_exp={ne}, q={qc.n}, |support|={zs.size}")
    print(f"    eps = +1 on {int((eps>0).sum())}, -1 on {int((eps<0).sum())}")

    # E1 -- multiplicativity on sampled triples
    rng = np.random.default_rng(0)
    tried = ok = 0
    for _ in range(200000):
        i, j = rng.integers(0, zs.size, 2)
        z1, z2 = int(zs[i]), int(zs[j])
        z3 = z1 ^ z2
        if z3 in emap:
            tried += 1
            if emap[z1] * emap[z2] == emap[z3]:
                ok += 1
            if tried >= 5000:
                break
    print(f"    E1 multiplicative on {ok}/{tried} sampled triples "
          f"-> {'CHARACTER' if tried and ok == tried else 'NOT a character'}")

    # E2 -- solve for v over GF(2) from eps on the support, then verify on all
    # eps(z) = (-1)^<z,v>  =>  <z,v> = 0 if eps=+1, 1 if eps=-1
    rows, rhs = [], []
    for z, e in zip(zs.tolist(), eps.tolist()):
        rows.append(z)
        rhs.append(0 if e > 0 else 1)
    # Gaussian elimination over GF(2) on the bit-columns
    v = 0
    pivots = {}
    for z, b in zip(rows, rhs):
        cur, curb = z, b
        while cur:
            hb = cur.bit_length() - 1
            if hb in pivots:
                cur ^= pivots[hb][0]
                curb ^= pivots[hb][1]
            else:
                pivots[hb] = (cur, curb)
                break
        else:
            if curb:
                print("    E2 INCONSISTENT -- no v exists (eps is not a character)")
                return
    # back-substitute
    v = 0
    for hb in sorted(pivots):
        row, b = pivots[hb]
        if (bin(row & v).count("1") + b) % 2 == 1:
            v |= (1 << hb)
    pred = np.array([(-1) ** (bin(int(z) & v).count("1") % 2) for z in zs])
    match = int(np.count_nonzero(pred == eps))
    print(f"    E2 solved v = {v} (0b{v:b});  eps(z)=(-1)^<z,v> holds on "
          f"{match}/{zs.size}  -> {'CONFIRMED' if match == zs.size else 'FAILS'}")

    # E3 -- decode v by register
    if match == zs.size and v:
        reg = {}
        for q_ in me.b: reg[q_] = "b"
        for q_ in me.t: reg[q_] = "t"
        for q_ in me.x: reg[q_] = "x"
        reg[me.c0] = "c0"; reg[me.anc] = "anc"
        for q_ in me.exp: reg[q_] = "e"
        bits = [i for i in range(qc.n) if (v >> i) & 1]
        print(f"    E3 v touches qubits {bits} = "
              f"{[f'{reg.get(i,'?')}{i}' for i in bits]}")
    elif match == zs.size:
        print("    E3 v = 0, i.e. h = g on the support (no shift needed)")


print("=" * 78)
print("Is the branch sign pattern a linear character?  (beta = 1 cases)")
print("=" * 78)
analyse(7, 6, 3, "N=7 a=6 (r=2)")
print()
analyse(7, 6, 4, "N=7 a=6 (r=2), wider")
print()
analyse(5, 2, 4, "N=5 a=2 (r=4)")
print()
analyse(5, 4, 3, "N=5 a=4 (r=2)")
