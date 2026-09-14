---
code: PK
title: "TODO 12d: measured peak deficit set and the two dominant Fourier modes (universal derivation withdrawn)"
outcome: record
claims: [C2, C17, C44]
todo: [12d]
---
# PK — TODO 12d SOLVED: the peak deficit is a COUNTABLE SET, and it is the two dominant Fourier modes

## two dominant Fourier modes

> **Review correction (2026-09-09).** The global peak formula below was
> over-promoted. The local Toffoli observations explain possible doubling but do
> not determine the global maxima after intervening CNOTs and interference. A
> four-qubit counterexample has atomic peak 4 and rotation peak 10, disproving
> the universal derivation. Retain the formula only as a measured regularity in
> the nine arithmetic rows; C44 is narrowed accordingly. The atomic
> peak/Walsh-suffix characterization remains valid.

`experiments/experiment_c17_deficit.py`. Claim **C44**; the arithmetic
regularity is retained, while the universal derivation is withdrawn. Paper A §5
is rewritten around the narrower result.

**The question.** Permutation-native propagation halves peak memory —
2.000000 exactly for ripple-carry adders, but **1.9997** for modexp, where all
three logged instances satisfied `rot = 2·perm − 2` exactly. A constant
additive deficit of 2, independent of instance size, was unexplained.

### The local gadget observation (not a universal peak derivation)

Two target-local facts motivated the original, invalid global inference.

1. A Z-type string commutes with every Z-rotation, so it **cannot branch until
   an H turns a Z into an X**. The gadget's only H acts on the target c, so
   only strings carrying Z_c ever leave the diagonal.
2. Once a string is X_c-type, the gadget's four T gates on c are rotations
   about Z_c, and those act *within* span{X_c, Y_c}:
   `X_c → cos·X_c − sin·Y_c`, `Y_c → cos·Y_c + sin·X_c`. That target-local
   space is closed only for a fixed rest-label component. Interleaved CNOTs
   change other labels, so it does **not** follow that full strings branch once.

Hence a Z_c-carrying string can contribute two Paulis within the gadget while
other strings contribute one at that local step. Intervening CNOTs and
interference prevent this local count from proving a relation between global
maxima. The arithmetic rows nevertheless satisfy the following measured
regularity:

```
    N_max^rot = 2·N_max^perm − |B|,     B = { z ∈ S : z_c = 0 }
```

with S the permutation-native peak set and c the target of the gadget in which
the rotation-level peak falls. **The deficit is not a constant; it is the
number of live strings that miss that gadget's own target qubit.**

### Verified 9/9, and sharper than a count

At the peak, every Pauli has X-support either empty or exactly {c} — the
x-masks are literally `{0, 128}` for the q=14 instances and `{0, 512}` for
q=17 (and `{8}`/`{16}`/`{32}` for the adders, where the identity sector is
absent precisely because |B| = 0). Folding the X_c/Y_c partners back onto
their parent Z-string recovers S **set for set**, with multiplicity 2 on
z_c = 1 and 1 on z_c = 0, in all nine instances.

### Why 2 — and this is the pretty part

```
  instance     perm      rot   deficit   B
  N=5  a=2    6,834   13,666      2      {x0}, {x0,e0}
  N=7  a=3    8,194   16,386      2      {x0}, {x0,e0}
  N=15 a=7   64,070  128,138      2      {x0}, {x0,e0}
  N=5  a=4    8,194   16,386      2      {x0}, {x0,e0}
  N=7  a=6    6,910   13,818      2      {x0}, {x0,e0}
  N=15 a=2   65,538  131,074      2      {x0}, {x0,e0}
  adders          -        -      0      empty
```

B is the same pair in every instance — and **those two strings are exactly the
Walsh coefficients of magnitude ½**, checked against the full final spectrum,
6/6. §W3/W4 found that same pair from a completely unrelated direction
(coefficient-magnitude analysis of the final spectrum: "the two large
coefficients are supported almost entirely on the x-register bit being measured
plus one exponent qubit"). So:

> **The deficit is 2 because the modexp bit function has exactly two dominant
> Fourier modes, and those modes live off the scratch register — they never
> acquire Z_c, so they never double.**

### The correction this forces: the constant belongs to the OBSERVABLE

Must-fail control C2, and it fires hard. Same circuit (N=5, a=2), observable
moved:

```
  observable    perm      rot   deficit
  x0  (standard) 6,834   13,666        2
  b0            12,196   20,388     4004
  anc           12,206   20,398     4014
  t0             1,024    2,048        0     <- ratio exactly 2, like the adders
```

"Modexp has deficit 2" is a statement about the **standard computational-basis
observable**, not about modular exponentiation. Paper A now says so.

### Grading, and what is still open

- `N_max^rot = 2·N_max^perm − |B|` — **measured 9/9 in the arithmetic rows.**
  It is not a universal derivation; the counterexample above disproves that
  broader claim. The exact atomic-peak/Walsh-suffix characterization remains.
- `|B| = 2` for the standard observable, with B the two |c| = ½ modes —
  **verified 6/6, mechanism identified, NOT proved.** Proving it needs a
  characterisation of which peak-time strings avoid the scratch register, which
  is a statement about the peak moment rather than about the final spectrum.
  Do not write it as a theorem.

### Method: the first pass was wrong, and the way it was wrong is worth keeping

Pass 1 (P1/P2) had the identity right but paired the peak against the Z-set at
the **gadget boundary**, reconstructed by re-emitting the circuit and mapping
gate indices to logical ops. 0/9. Diagnosis from dumping the peak set: the peak
sits deep *inside* the gadget, after the Toffoli's own 4-way expansion has run
at rotation level, so no boundary carries the right set — the partner is the
perm-level **peak**. And the gadget target c never needed reconstructing at
all: it is the single bit set in the peak's common x-mask. **A hundred lines of
index bookkeeping were solving a problem that did not exist.** Meanwhile the
checks that only *counted* B (P4/P5) passed 6/6 through both passes, which is
what made the diagnosis quick — the failure was localised to the boundary
machinery, not to the idea.
