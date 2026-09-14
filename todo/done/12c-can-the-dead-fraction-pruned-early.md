---
id: 12c
state: done
title: Can dead exponent terms be pruned early?
outcome: "Solved by exact last-reverse-use contraction in both order branches; repeated β=1 tails additionally compress to one block"
claims: [C18, C23, C24, C42, C45, C46]
---
# Can dead exponent terms be pruned early?

**Solved 2026-09-09.** See C45/C46 and `notes/CT-finished-control-contraction.md`
for the safe contraction schedule, proof, implementation and measurements.
The original anticipated useful fractions were approximate, not exact (C42).
The proposed β>1 negative control below was also wrongly scoped: finished
controls can be contracted in either branch. β>1 instead defeats the separate
identity-tail compression. The original proposal follows as historical context.

Opened by §BI/C42. For the real Shor input state a fraction 1 − 2^−(α+1) of the
Walsh support contributes **exactly zero** to ⟨O⟩ — 75% at α=1, 87.5% at α=2,
verified 3/3. But PPS propagates backwards and only meets the input state at
the end, so those terms are carried at full cost and *then* discarded: peak
memory (C18) is unchanged. **The question is whether the dead set can be
recognised early.** If a cheap invariant identifies "this branch will end with
exponent support" partway through propagation, that is a real constant-factor
win of 4× (α=1) to 8× (α=2) in peak memory, on top of everything else — and
unlike δ-truncation it is **exact**, not approximate.

> **The naive version is DEAD — answered the same day it was written, by
> derivation, before any effort went into it.** Exponent support is **not
> monotone** under back-propagation, so a term cannot be pruned when it
> acquires it. Through `CCX(a,b,c)` the pullback of `Z_c` is
> `½(Z_c + Z_aZ_c + Z_bZ_c − Z_aZ_bZ_c)`, and the pullback of `Z_aZ_c` is
> `½(Z_aZ_c + Z_c + Z_aZ_bZ_c − Z_bZ_c)` — because `Z_a·Z_a = I`. Verified
> exactly: `Z_aZ_c` produces the term `Z_c` with coefficient +0.50, carrying
> **no** support on the control qubit a. A term with exponent support can
> therefore lose it and go on to contribute, so zeroing it early is simply
> wrong, not merely suboptimal.

What survives the above: is there a *certificate* weaker than "has exponent
support now" that predicts the final exponent support of a whole branch? The
C23/C24 mechanism is known analytically (the tail is entered through a single
parity bit), so the dead set is characterised at the END — the question is
whether that characterisation can be pushed backwards through the propagation.
Unclear, and much less likely to be cheap than it looked. A must-fail control
is easy: the β>1 case should show no such structure.

Related: the same framing applies to any input state via ⟨O⟩ = Σ_z c_z ∏ δ_i,
so a partially-biased register interpolates between "free" and "carried".
