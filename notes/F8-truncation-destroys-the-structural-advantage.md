---
code: F8
title: "truncation destroys the structural advantage"
outcome: record
claims: []
todo: []
---
# F8 — truncation destroys the structural advantage

At δ=1e-3 the ordering **reverses**: Toffoli N_max 24692, Fourier 11581. The
compilation that is better exactly is worse under truncation, because δ-truncation
discards terms that were going to cancel, so the exact Z-closure cancellation
never completes. Norm violations appear on both (⟨O⟩ up to +1.109).

Practical reading: telling someone "compile to Toffoli before running PPS" is
**only** sound advice at δ=0, which is not a regime anyone runs in. This
substantially weakens the applied claim.

---
