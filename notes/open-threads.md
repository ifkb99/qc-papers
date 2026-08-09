---
code: open-threads
title: "Open threads"
outcome: record
claims: [C15]
todo: []
---
# Open threads

1. ~~Real modular exponentiation~~ **DONE** — `modexp.py`, verified end to end
   (`test_modexp.py` A–G all pass; N=5,7,15; period peaks land exactly).
2. ~~Re-run F2/F3 on real modexp~~ **DONE** — F2 did not survive. See STATUS block.
3. **[NOW THE MAIN THREAD] Toffoli-arithmetic vs Fourier-arithmetic A/B.**
   F4 says the compilation choice, not the logical function, sets PPS cost.
   Needs: a Toffoli-based modular multiplier to put beside the Beauregard one,
   so the comparison is like-for-like on the *same* N, a, and observable.
   Currently only have a Toffoli *adder* (non-modular) — that is the next build.
   Prediction to test: Toffoli modexp gives small `N_max` (permutation ⟹ Z-type
   preserved), Fourier modexp gives large `N_max`, identical logical circuit.
   If that holds cleanly it is a genuine, checkable statement about when PPS is
   applicable, and it is actionable (compile arithmetic to Toffoli before PPS).
4. Accuracy study now actually possible: non-degenerate `<O>` exist on the real
   circuit. Sweep δ, plot error vs `N_max`, use the |<O>| ≤ 1 check as a
   validity gate. Expect a cliff, but measure it properly this time.
5. Scale test for the pre-asymptotic worry (needs bit-packed PPS, App. B).
   Current Python dict PPS: ~15s for 8000 gates at 11 qubits, N_max ~50k.

6. **[HIGHEST VALUE, not started] Weight-truncation as Fourier tail mass.**
   The Pauli weight of `Z^z` is `popcount(z)`, which is exactly the Fourier
   degree of that coefficient. So for permutation circuits, **weight-truncation
   error is literally the Fourier tail mass** of the pulled-back bit function.
   The current work covers only coefficient-truncation (δ); weight-truncation is
   the other standard PPS knob and is completely untouched. Measuring
   mass-by-weight profiles (adder vs modexp vs the binomial profile of a random
   function) would extend the exact cost model to it. Logged in
   `experiment_review.py` R3.
7. **Prove the circuit-level C15 invariance.** The valid-subspace argument does
   not cover it (added qubits are live). Most promising route: is the support
   confined to an affine subspace? That would also explain density → ½ exactly.
   Paper B open problem 1.
8. Decode *which* qubits carry the 4 dominant Walsh coefficients that reproduce
   ⟨O⟩ exactly (F13) — would turn the aggressive-truncation rule from empirical
   to structural.
