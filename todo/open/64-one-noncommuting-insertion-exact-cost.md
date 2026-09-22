---
id: 64
state: open
title: "Queued: exact Heisenberg term count of the modexp pullback with one Rx(theta) between blocks, and its growth factor R(q)"
claims: [C8, C52, C53, C56, C99, C104]
---
# One non-permutation gate between blocks (queued)

**Triage (2026-09-22, METHOD.md barrier check 4):** an exact term count at small sizes; computes nothing known to be hard. Passes.

From the SN slate (2026-09-18), deriver candidate C1, which sits beside PAPER_B §13
item 2. The deriver's diagnostics (exact, CPU, N = 3, 5, 7, t = 2) show three things.
The exact support of the Z_x0 pullback is the same for every θ with irrational cos θ,
and a rational cos θ loses only diagonal terms. The off-diagonal branch has 27×, 110×
and 220× the permutation count, as about R(q) X-classes of about 0.41·2^(q−1) terms
each, with R a DDT-row count of the prefix permutation. No exponent bit appears in any
X-mask, so C104's slices survive one insertion.

**Premise check first.** The deriver cites none of C52–C57, which already treat single
defects, including Rx(θ) kicks, in known-order circuits, and their detectable
coherences. Before a derivation, establish which of those claims concern this object
(the pre-QFT Heisenberg term count) and which concern output sampling. The next
derivation's content would be R(q) derived from the carry butterfly, tested at
N = 11, q = 18.
