---
id: 12d
state: done
title: Why is the peak ratio exactly rot = 2·perm − 2 for modexp?
claims: [C17, C44]
---
# Why is the peak ratio exactly `rot = 2·perm − 2` for modexp?

**Review correction, 2026-09-09.** The universal derivation in the historical
resolution below is false. C44 retains only the measured arithmetic relation;
C17 retains the exact suffix-Walsh characterization, not a factor-two bound.
See `experiments/experiment_independent_review.py` for the four-qubit
counterexample and the corrected §PK note. The earlier Proposition 2 has
been removed from Paper A and subsequent statements renumbered.

**SOLVED 2026-08-08 — and it resolved in an afternoon, as the note below
guessed.** See `NOTES.md` §PK; file `experiments/experiment_c17_deficit.py`;
claim **C44**; C17 regraded from "the factor is empirical" to derived; Paper A
§5 rewritten as Proposition 2.

**The identity.** Only Z_c-carrying strings meet the gadget's H and leave the
diagonal, and the four T gates on c are rotations *about* Z_c, which act inside
span{X_c, Y_c} — a closed 2-dimensional space. So they branch **once between
them**, not 2⁴ times, and

```
    N_max^rot = 2·N_max^perm − |B|,    B = { z ∈ S : z_c = 0 }
```

with c the target of the gadget containing the rotation-level peak. **The
deficit is a countable set, not a constant.** Verified 9/9, and sharper than a
count: at the peak every Pauli has X-support ∅ or exactly {c}, and folding the
X_c/Y_c partners recovers the perm peak set *set for set*.

**Why 2.** B = {Z_x0, Z_x0·Z_e0} in 6/6 modexp instances — and those are
exactly the **|coefficient| = ½ Walsh terms**, i.e. the dominant Fourier modes
that §W3/W4 had already singled out from the *final* spectrum by an unrelated
route. They live off the scratch register, never acquire Z_c, never double.

**The correction it forces.** The constant belongs to the **observable**, not
to modexp: same circuit, observable moved to b0 → deficit 4004, anc → 4014,
t0 → 0 (ratio exactly 2, like the adders).

**Left open, deliberately ungraded:** `|B| = 2` itself is measured, not proved.
Proving it needs a characterisation of which peak-time strings avoid the
scratch register — a statement about the peak moment, not about the final
spectrum. The identity above is the theorem; the membership is evidence.

**Method note worth keeping.** Pass 1 had the identity right and the *pairing*
wrong: it compared the peak against the Z-set at the gadget **boundary**,
reconstructed with ~100 lines of gate-index bookkeeping. 0/9. The peak actually
sits mid-gadget, and the gadget target never needed reconstructing — it is the
single bit set in the peak's common x-mask. The checks that only counted B
passed 6/6 through both passes, which localised the failure immediately.

<details><summary>Original note, kept for the record</summary>

Opened 2026-08-08 by the C17 correction. The permutation-native peak reduction
is **2.000000 exactly** for ripple-carry adders (128/64, 512/256, 2048/1024) but
**1.9997** for modular exponentiation — and every modexp instance satisfies

```
  rot = 2·perm − 2       13666 = 2·6834 − 2
                         16386 = 2·8194 − 2
                        128138 = 2·64070 − 2
```

exactly, 3 of 3, across N = 5, 7, 15 and wildly different peak magnitudes. A
constant additive deficit of exactly 2 that is independent of instance size is
not a coincidence and is currently **unexplained**.

**The half that is understood.** At the gadget's Hadamard on target qubit c, a
Z-type string containing Z_c is carried to a mirrored X-sector partner while a
string without Z_c is untouched. So the doubling applies only to the
Z_c-containing subset, which bounds the ratio above by 2 and explains why the
adders attain it (there, at the peak, every string apparently contains Z_c) and
modexp does not.

**The half that is not.** Why the shortfall is exactly 2 rather than
instance-dependent. Two obvious candidates, neither checked: the identity string
Z^0 (which cannot acquire Z_c and so never doubles) plus one partner; or a
parity constraint pinning a second string. Cheapest first move is to dump the
peak-time Pauli set for the smallest modexp instance and simply *look at* which
two strings fail to double — this is a 2-line diagnostic, not a research
programme, and it either resolves in an afternoon or reveals something.

If it resolves, `PAPER_A.md` §5 upgrades from "upper bound of 2, empirical" to a
theorem, which is worth having: it is currently the only quantitative claim in
Paper A that is measured rather than derived.

</details>

*The "cheapest first move" above — dump the peak-time Pauli set and simply look
at it — was exactly right, and is what produced the answer. Both named
candidates (the identity string Z^0; a parity constraint) were wrong.*
