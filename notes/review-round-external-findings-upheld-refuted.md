---
code: review-round-external-findings-upheld-refuted
title: "REVIEW ROUND (external) — 3 of 4 findings upheld, 1 refuted"
date: 2026-08-08
outcome: record
claims: [C3, C7, C8, C11, C15, C17]
todo: []
---
# REVIEW ROUND (external) — 3 of 4 findings upheld, 1 refuted

An independent review raised four objections. Verified each rather than
accepting; results below. **Net: the core (F9/C8) is unharmed, C15 is
strengthened, one real gap found and fixed, one review claim was wrong.**

### R1 UPHELD (important) — Walsh predicts N_final, not N_max (peak memory)

PPS memory cost is the *peak* term count during propagation, not the final one.
Walsh sparsity equals the **final** count. These differ:

```
        circuit    walsh/final    rot N_max   ratio
 3-bit adder Z_b0            1          128    128x
 4-bit adder Z_b0            1          512    512x
 5-bit adder Z_b0            1         2048   2048x
      modexp N=5         3086        13666    4.4x
     modexp N=15        31176       128138    4.1x
```

**So F1's "the adder is FREE for PPS" is misleading** — it peaks at the full
lightcone (2^(q-1)) before collapsing to 1. And C11's "exact cost model" predicts
*result complexity*, not peak memory. Both need restating. For modexp the ratio
is a bounded ~4x, so **C7's Θ(2ⁿ) asymptotic is unaffected**.

### R1 FIX — permutation-native PPS (`perm_pps.py`), a real improvement

Root cause: standard PPS decomposes Toffoli into Clifford+T and propagates
through the *rotations*, which include Hadamards — so the operator leaves the
diagonal mid-circuit even though the Toffoli as a whole is a permutation.

Treating X/CNOT/Toffoli as **atomic permutation gates** keeps everything Z-type
at every step. Conjugation rules (from `(-1)^{popcount(perm(y) & z)}`):

```
X(q)           : Z^z -> (-1)^{z_q} Z^z                        1 term
CNOT(c,t)      : Z^z -> Z^{z XOR (z_t << c)}                  1 term
Toffoli(a,b,c) : z_c=0 -> Z^z                                 1 term
                 z_c=1 -> 1/2[Z^z + Z^{z^a} + Z^{z^b} - Z^{z^a^b}]   4 terms
```

Measured (verified against Walsh and against rotation-level PPS,
`test_perm_pps.py` all pass):

```
        circuit    final   rot N_max   perm N_max   saving
 3-bit adder Z_b0      1         128           64     2.0x
      modexp N=5    3086       13666         6834     2.0x
     modexp N=15   31176      128138        64070     2.0x
```

> **CORRECTED 2026-08-08 — "exactly 2.0x everywhere" is a ROUNDING ARTIFACT of
> the `%.1f` in the table above.** The true ratios are **2.000000** for the
> adders (128/64, 512/256, 2048/1024) but **1.9997** for modexp:
> 13666/6834, 16386/8194, 128138/64070. Every modexp instance satisfies
> `rot = 2·perm − 2` exactly, 3/3 — a clean regularity we have **not**
> explained. The upper bound of 2 is arguable from the gadget's H exchanging
> the Z- and X-sectors (strings containing Z_c get a mirrored partner, strings
> without Z_c do not), which also explains why the adders attain it and modexp
> does not; the constant deficit of 2 is unexplained. Claim C17 regraded: the
> factor is **empirical, not proved**. Caught by Fable reviewing `PAPER_A.md`
> — the abstract said "halves peak memory exactly" while §5 supported only
> "2.0× on all 6 instances", which are different epistemic states.

At the peak the rotation-level run holds Z and non-Z terms in near-equal
measure, so dropping the non-Z part is close to a factor of two. Modest, but two
things matter more than the constant:

1. The peak becomes a **well-defined Walsh quantity**: after k gates the term
   count *is* the Walsh sparsity of the k-gate suffix's pullback, so
   `N_max = max over suffixes of Walsh sparsity`. Exact, not extrapolated.
2. It is **much faster** — logical ops instead of rotations (1699 vs 21247 for
   N=5) and half the terms. Did q=21 in 7.4s where rotation-level PPS failed to
   finish q=17 in 20 minutes.

Actionable claim for a PPS engine: *implement reversible arithmetic blocks as
atomic permutation primitives; the propagation never leaves the diagonal.*

### R1 consequence — C15 SURVIVES for peak memory, not just final count

The critical test the review did not run. Using `perm_pps`:

```
   a   r  n_exp    q   N_final  N_max(peak)  peak growth
   6   2      2   15     15549        24369            -
   6   2      4   17     15549        24369        1.00x
   6   2      6   19     15549        24369        1.00x
   6   2      8   21     15549        24369        1.00x

   3   6      2   15     15539        24412            -
   3   6      4   17     64353        98018        4.02x
```

Peak is **exactly** constant at 24369 across a 64x growth in dimension for r=2,
against 4.02x/step for r=6. **C15's headline holds for memory cost, not merely
result size.** (Also confirmed at rotation level: N_max 40753 at both n_exp=2
and n_exp=3.)

### R2 UPHELD — the C15 *mechanism* stated in ABSTRACT.md is wrong

Checked directly: for r=2, n_exp=6, the added exponent qubits (15–18) are each
set in ~7779/15549 support terms — **genuinely live**, not confined. So the
abstract's "the Walsh support is confined to the low bits of the exponent
register" is false at circuit level. It is true for the *idealised* function
(F12), where the domain is e alone. Cause: `u_a(ctrl, 1)` is identity only on
the valid subspace; as a full-space unitary it entangles the control with ancilla
garbage. The count is constant; the support is not confined.

NOTES already carried this caveat; **the abstract did not, and has been fixed.**

### R3 REFUTED — the "2-adic vs factors of r" objection was based on the abstract only

The review fetched only Dang et al.'s abstract ("depends on its factors") and
concluded MPS keys on general smoothness while PPS keys on powers of two,
offering r=6 as a counterexample. **The full text says otherwise.** Section 4:

> "α is the number of trailing zeroes in the binary representation of r" …
> "due to the **odd factor β ≡ r/2^α** of r which cannot be localised to
> specific qubits"

and §5.2: MPS matrices for the qubits in A "reduced by a factor of **β²**".

So Dang et al. decompose r = β·2^α *exactly* as 2-adic valuation times odd part,
and it is the **odd part β that costs memory** — β=1 (r a power of two) is the
free case. r=6 has β=3, so it is in the expensive branch for MPS too, and is
**not** a counterexample. The correspondence is real and now citable with a
specific mechanism. Claim strengthened rather than weakened.

### R4 UPHELD in substance, but its numbers are confounded

True and useful: at rotation level the Toffoli decomposition passes through H
gates, so the operator does leave the diagonal mid-circuit — "Z-closure" is a
property of the *endpoint*, not the trajectory. That is exactly R1, and
`perm_pps.py` fixes it.

But the specific figures offered (Toffoli peak non-Z 13664 vs Fourier 848) compare
**q=14 against q=9** — the same ancilla/qubit-count confound already logged as
killing C3. Not a valid compilation comparison.

---
