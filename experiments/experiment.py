"""PARTIALLY STALE -- section 3 is invalid, sections 1/2/4/5 are fine.

Section 3 predates the 1e-13 noise floor in pps.py and was reading pure
floating-point residue (it reports coefficients at 2^-105 etc.). Ignore it.
Sections 1 and 2 (branch weights, coefficient distributions) still stand.
Not affected by the theta=pi bug: no X/Y/Z gates appear in these circuits.

---

Does the PPS resource-prediction framework transfer to arithmetic circuits?

The framework in Gharibyan et al. (arXiv:2507.10771) rests on the empirical
claim that the Pauli-coefficient distribution follows a power law
rho(t) ~ A/t^(m+1), which is what makes N_max extrapolation from cheap test
runs possible (their Eqs. 11, 17).

That was measured on brickwork circuits with *generic, uncorrelated* rotation
angles. Clifford+T arithmetic has every non-Clifford angle at exactly pi/4.
This script measures both families under identical conditions.
"""
from __future__ import annotations
import numpy as np, collections
from circuits import Circuit, ripple_adder, brickwork
from pps import propagate, exact_expectation, fit_power_law


def branch_weight_profile(circuit):
    """|sin(theta)| for every branching (non-Clifford) gate."""
    w = []
    for _, th in circuit.gates:
        s = abs(np.sin(th))
        if 1e-9 < s < 1 - 1e-9:
            w.append(s)
    return np.array(w)


def describe(coeffs, delta):
    u = np.unique(np.round(coeffs, 12))
    m, nfit = fit_power_law(coeffs, delta)
    return {
        "n_terms": len(coeffs),
        "n_distinct_magnitudes": len(u),
        "ratio": len(u) / max(len(coeffs), 1),
        "m": m,
        "max": coeffs.max() if len(coeffs) else 0,
        "min": coeffs.min() if len(coeffs) else 0,
    }


print("=" * 74)
print("1. BRANCH WEIGHTS: what each family feeds the truncator")
print("=" * 74)
rng = np.random.default_rng(1)
bw = brickwork(8, 12, rng)
add, lay = ripple_adder(3)
for name, qc in (("brickwork (random angles)", bw), ("3-bit ripple adder", add)):
    w = branch_weight_profile(qc)
    print(f"  {name:28s} branching gates={len(w):5d}  "
          f"|sin| mean={w.mean():.4f} min={w.min():.4f} max={w.max():.4f} "
          f"distinct={len(np.unique(np.round(w, 9)))}")
print("""
  Truncation works by discarding small coefficients. A branch multiplies a
  coefficient by cos(t) on one side and sin(t) on the other. Small |sin| means
  one child is negligible and can be dropped for free. |sin| = 1/sqrt(2) means
  both children are equal -- nothing is safe to drop.""")

print()
print("=" * 74)
print("2. COEFFICIENT DISTRIBUTION at fixed delta")
print("=" * 74)
delta = 1e-4
obs_bw = {(0, 1 << 4): 1.0}
obs_ad = {(0, 1 << lay["b"][0]): 1.0}
rows = []
for name, qc, obs in (("brickwork", bw, obs_bw), ("adder", add, obs_ad)):
    r = propagate(qc, obs, delta=delta, max_terms=400_000)
    d = describe(r.final_coeffs, delta)
    rows.append((name, d, r))
    mm = f"{d['m']:.3f}" if d["m"] is not None else "n/a"
    print(f"  {name:12s} terms={d['n_terms']:7d}  distinct |c| values={d['n_distinct_magnitudes']:7d}"
          f"  (ratio {d['ratio']:.4f})  power-law m={mm}")
print("""
  'distinct |c| values / terms' near 1.0 means a smooth continuum of
  coefficients -- a power law can describe it. Near 0 means the coefficients
  pile onto a few discrete magnitudes: the Dirac-spike limit the paper
  identifies (Sec 3.1) as what happens when every angle is identical.""")

print()
print("=" * 74)
print("3. WHERE THE ADDER'S COEFFICIENTS ACTUALLY SIT")
print("=" * 74)
r = propagate(add, obs_ad, delta=0.0, max_terms=400_000)
c = r.final_coeffs
hist = collections.Counter(np.round(-2 * np.log2(np.maximum(c, 1e-300)), 6))
print("  exact run (no truncation): every |c| should be 2^(-k/2) for integer k")
print(f"  {'k':>6}  {'|c| = 2^(-k/2)':>16}  {'count':>8}")
for k in sorted(hist)[:12]:
    print(f"  {k:6.2f}  {2**(-k/2):16.8f}  {hist[k]:8d}")
frac = sum(v for k, v in hist.items() if abs(k - round(k)) < 1e-6) / max(len(c), 1)
print(f"  fraction of coefficients at exactly 2^(-k/2): {frac:.4f}")

print()
print("=" * 74)
print("4. DOES TRUNCATION CONVERGE? (adder, exact ground truth available)")
print("=" * 74)
small, lay2 = ripple_adder(2)
obs2 = {(0, 1 << lay2["b"][0]): 1.0}
truth = exact_expectation(small, obs2)
print(f"  exact <O> = {truth:.10f}")
print(f"  {'delta':>10} {'N_max':>8} {'<O> est':>14} {'abs err':>12} {'dropped wt':>12}")
for d in (0.0, 1e-6, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
    rr = propagate(small, obs2, delta=d, max_terms=400_000)
    print(f"  {d:10.0e} {rr.n_max:8d} {rr.expectation:14.8f} "
          f"{abs(rr.expectation - truth):12.2e} {rr.truncated_weight:12.2e}")

print()
print("=" * 74)
print("5. N_max SCALING WITH delta  (paper predicts N_max ~ delta^-m)")
print("=" * 74)
for name, qc, obs in (("brickwork", bw, obs_bw), ("adder", add, obs_ad)):
    print(f"  {name}:")
    prev = None
    for d in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        rr = propagate(qc, obs, delta=d, max_terms=400_000)
        slope = ""
        if prev is not None and rr.n_max > 0 and prev[1] > 0:
            slope = f"  d(log N)/d(log 1/delta) = {np.log(rr.n_max/prev[1])/np.log(prev[0]/d):6.2f}"
        print(f"    delta={d:8.0e}  N_max={rr.n_max:7d}{slope}")
        prev = (d, rr.n_max)
print("""
  A clean power law gives a constant slope (= m). A slope that collapses to ~0
  means delta has stopped controlling cost: below the spike, nothing is
  discarded; above it, everything is.""")
