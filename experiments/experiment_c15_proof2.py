"""TODO step 3, attempt 2 -- the branch-symmetry mechanism.

Attempt 1 (affineness of u_a(ctrl,1)) is REFUTED: the block is not affine
(14336 violations of 32768) and conjugating Z_j by it gives 2504 terms, not 1.

New route, derived rather than guessed. Adding exponent qubit e appends a block
CONTROLLED on e, so the pulled-back bit function splits into two branches:

    f(y, e) = g(y)  if e = 0        (block acts trivially, control off)
              h(y)  if e = 1        (block acts)

The Walsh coefficients of f in terms of those of (-1)^g and (-1)^h are

    c_(z,0) = ( g^(z) + h^(z) ) / 2
    c_(z,1) = ( g^(z) - h^(z) ) / 2

Three regimes:
  h = g            -> c_(z,1) = 0 for all z: size preserved but the new qubit
                      would be DEAD. Contradicted by measurement (it is live).
  supp disjoint    -> both survive everywhere: size DOUBLES. Not observed.
  |h^(z)|=|g^(z)|  -> exactly ONE of each pair survives: size preserved AND the
                      new qubit live. This is the observed behaviour.

So the mechanism to verify is **branch symmetry**: the two branches have Walsh
spectra equal in magnitude, differing only in sign pattern.

  M1  For each z, is exactly one of (z,0), (z,1) in the support?
  M2  Do the surviving magnitudes reproduce the n_exp = k spectrum exactly?
  M3  Does branch symmetry hold for beta=1 and FAIL for beta>1? (It must fail
      for beta>1, or the argument would prove constancy for every r.)
"""
from __future__ import annotations
import numpy as np
from toffoli_arith import ToffoliModExp
import walsh

TOL = 1e-12


def spectrum(N, a, ne):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    qc = me.build()
    return me, qc, walsh.pullback_coefficients(qc, me.x[0])


def branch_analysis(N, a, ne, label):
    """Split the n_exp=ne spectrum by its TOP exponent bit and test M1/M2."""
    me, qc, c = spectrum(N, a, ne)
    top = me.exp[-1]                      # highest exponent qubit
    assert top == qc.n - 1, (top, qc.n)
    half = 1 << top
    lo = c[:half]                          # z with top bit 0
    hi = c[half:]                          # z with top bit 1
    nz_lo = np.abs(lo) > TOL
    nz_hi = np.abs(hi) > TOL

    both = int(np.count_nonzero(nz_lo & nz_hi))
    exactly_one = int(np.count_nonzero(nz_lo ^ nz_hi))
    neither = int(np.count_nonzero(~nz_lo & ~nz_hi))
    total = int(nz_lo.sum() + nz_hi.sum())

    # M2: surviving magnitude vs the smaller circuit's spectrum
    _, qc_prev, c_prev = spectrum(N, a, ne - 1)
    surv = np.where(nz_lo, np.abs(lo), np.abs(hi))
    prev_abs = np.abs(c_prev)
    m2_ok = prev_abs.size == surv.size and np.allclose(surv, prev_abs, atol=1e-12)
    m2_err = float(np.max(np.abs(surv - prev_abs))) if prev_abs.size == surv.size else float("nan")

    print(f"  {label}")
    print(f"    n_exp={ne}, q={qc.n}, |support|={total} "
          f"(n_exp={ne-1} had {int((prev_abs>TOL).sum())})")
    print(f"    pairs with BOTH members present : {both}")
    print(f"    pairs with EXACTLY ONE          : {exactly_one}")
    print(f"    pairs with neither              : {neither}")
    print(f"    M1 (exactly-one everywhere)     : {both == 0 and neither + exactly_one == half}")
    print(f"    M2 (magnitudes match n_exp-1)   : {m2_ok}   max err {m2_err:.2e}")
    return both, exactly_one, m2_ok


print("=" * 78)
print("M1/M2  Branch symmetry where the invariance HOLDS (beta = 1)")
print("=" * 78)
branch_analysis(7, 6, 3, "N=7 a=6 (r=2, alpha=1): locked from n_exp=2")
print()
branch_analysis(7, 6, 4, "N=7 a=6 (r=2, alpha=1): another width")
print()
branch_analysis(5, 2, 4, "N=5 a=2 (r=4, alpha=2): locked from n_exp=3")
print()

print("=" * 78)
print("M3  Does it FAIL where the invariance fails (beta > 1)?")
print("=" * 78)
branch_analysis(7, 3, 3, "N=7 a=3 (r=6, beta=3): support grows 4x/step")
print()
branch_analysis(21, 4, 3, "N=21 a=4 (r=3, beta=3): pure odd order")

print("""
=" * 78
  If M1 holds exactly for beta=1 and fails for beta>1, the invariance is
  explained: the identity block acts on the observable so as to preserve the
  magnitude of every Walsh coefficient while redistributing signs across the
  new qubit. Size is then conserved term-for-term, which is exactly what the
  measurements show -- and the new qubit is live precisely because the sign
  pattern, not the magnitude, is what changes.""")
