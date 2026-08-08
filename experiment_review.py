"""UNFINISHED STUB -- planned checks only, no code. Kept for the ideas.

Started during the external review pass; the file was left truncated. The
proposals below are live open threads and are also logged in NOTES.md under
"Open threads" so they are findable without reading this file.

R1  Support structure across n_exp for r=2 (the "exactly 15549 = SAME" claim).
    Is the support literally nested as n_exp grows, or is it relabelled? How
    much spectral mass sits on the new exponent bits? Is it confined to an
    affine subspace -- which would explain density -> 1/2 exactly?
    STATUS: partially answered during review. The added exponent qubits are
    LIVE (each set in ~7779/15549 terms), so the support is relabelled rather
    than nested. The affine-subspace question is still open and is the most
    promising route to PROVING the circuit-level invariance (Paper B, open
    problem 1).

R2  Prediction test: N=5, a=2 has r=4=2^2, so u_a blocks with i>=2 multiply by
    1. Sparsity should therefore be constant for n_exp >= 2.
    STATUS: the equivalent test was run for N=7 a=6 (r=2) and N=21 a=8 (r=2)
    and confirmed. N=5 a=2 specifically was not swept -- cheap to add as a
    third modulus for Paper B.

R3  Fourier-mass-by-weight profiles: adder vs modexp vs the binomial profile of
    a random function. The Pauli weight of Z^z is popcount(z), which is exactly
    the Fourier degree, so **weight-truncation error for PPS on permutation
    circuits IS the Fourier tail mass**. This is the sharpest unexplored idea
    here: it would extend the exact cost model from coefficient-truncation
    (delta) to weight-truncation (the other standard PPS knob), which the
    current work does not cover at all.
    STATUS: not started. Highest value of the four.

R4  Decode the dominant coefficients by register -- which qubits carry the few
    large Walsh coefficients that reproduce <O> exactly (4 of 15493, see F13).
    Would turn the aggressive-truncation rule from empirical to structural.
    STATUS: not started.
"""

raise SystemExit(__doc__)
