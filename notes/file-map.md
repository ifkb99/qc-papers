---
code: file-map
title: "File map"
date: 2026-08-08
outcome: retracted
claims: []
todo: []
---
# File map

Run everything with `uv run python` from `research/` (see trap 5).

**2026-08-08 restructure:** experiment scripts moved to `experiments/`
(filenames unchanged; run as `uv run python -m experiments.<name>`), fully
retracted ones to `archive/`, and the shared machinery was extracted into the
`lab/` package (`lab.harness`, `lab.measure`, `lab.gf2`, `lab.modarith`,
`lab.variants`, `lab.nulls`) with `test_lab.py` pinning it to the numbers
logged here and `test_claims.py` re-verifying the headline claims. Bare
`experiment_*.py` names below predate the move.

```
perm_pps.py      permutation-native PPS: X/CNOT/Toffoli as atomic gates,
                 stays Z-type throughout, 2x lower peak, ~10x faster
walsh.py         classical_permutation, FWHT, pullback_coefficients,
                 walsh_sparsity, is_affine, permutation_via_statevector
toffoli_arith.py Toffoli-compiled modexp (Cuccaro + add/sub-N reduction)
windowed_arith.py windowed modexp (§WD): WindowedModExp = table lookup +
                 register multiply; SelectModExp = select-multiply with the
                 skip_zero knob; replay()/verify_modexp() gate correctness on
                 circuits far too wide for a permutation array
pauli.py         symplectic algebra, rotation rule
circuits.py      Circuit (σ,θ list), gate decomps, Cuccaro adder, QFT+swaps,
                 cswap/ccphase/inverse, brickwork
statevec.py      O(2^n) per-gate state-vector sim (to_unitary dies past ~12q)
pps.py           propagate(), exact_expectation(), fit_power_law()
modexp.py        Beauregard modexp: phi_add -> cc_phi_add_mod -> cmult_mod
                 -> u_a -> build()/build_shor().  VERIFIED.
test_core.py     correctness gate — run first, always
test_modexp.py   modexp layers A-G — run second
experiment.py    branch weights + coeff distribution  (§3 STALE, trap 1)
experiment2.py   F1 mechanism + gate verification     (F2 section RETRACTED)
experiment3.py   old toy scaling                      (RETRACTED, see STATUS)
```

Order-finding sanity anchor: `ModExp(15, 7, n_exp=4).build_shor()` must give
exponent-register peaks at y = 0, 4, 8, 12 with p = 0.25 each (r = 4).
