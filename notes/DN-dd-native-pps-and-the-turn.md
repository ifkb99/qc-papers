---
code: DN
title: "DD-NATIVE PPS (TODO 13 rounds, 2026-09-15/17). A decision diagram of the propagated operator: ~2x fewer nodes than a matched null, 12-41x more bytes, and C103's crossover refuted as the mechanism for the node-count turn."
outcome: mixed
claims: [C102, C103, C104]
todo: [13, 55]
---
# DN — DD-native PPS propagation and the peak-node turn

Two board rounds (2026-09-15/16) plus coordinator experiments (2026-09-16/17) on
whether a decision diagram of the *propagated operator* beats the sparse Walsh
dictionary. C102/C103 own the canonical diagram of the *final function*; this note owns
the propagation-native route and the negative result about its one interesting feature.
The route is superseded for the peak-memory goal by C104 (note ES), which removes the
peak above a streamed output without diagrams.

## The object and what is proved

CNOT is memory-free in Heisenberg form: `perm_pps` maps a key by z ^ (z_t << c), a
bijection, so term count cannot change; under C82's lazy frame the dictionary is not
touched at all. Toffoli is the only gate that grows memory, and its branch is an
XOR-convolution on a 2-dimensional space in frame coordinates. Two derived rules were
verified against the engine on 300/300 random fixtures with a 300/300 control (gating:
only Claim 1, the predicate-invariance rule, is harness-gated; see the promotion bullet
below):
the predicate parity is frame-invariant, and the scatter rule has a gather form. They
are implemented in the frozen pilot `out/agent-board/artifacts/dd-pilot-20260915/`
(ddprop.py sha256 d00cc123…), whose MTBDD is the propagator all later work used.

## Measured, accepted on the board

**Node counts** (task `T563881b120674641`, submission `S18ba76c1588945c1`, accepted
`Vb2484fd40baa44bf` with four binding corrections I1–I4). ToffoliModExp(N=7, a=3),
plain convention, peak DD nodes over peak dictionary terms, n_exp = 1..6:

    dict peak   8,194  24,412  48,855  98,018  196,060  392,458
    plain DD    2,301  12,657  30,989  63,673  113,648  185,608
    ratio      0.2808  0.5185  0.6343  0.6496   0.5797   0.4729
    weighted   0.2235  0.4173  0.3902  0.4399   0.3734   0.2896

The plain ratio rises then falls, with its maximum at n_exp = 4; that is the "turn".
Against a matched null (same qubit count, support cardinality and value multiset,
random placement) the real/null node ratios are 0.351 / 0.534 / 0.648 / 0.667 / 0.606 /
0.499, so the diagram is about 2x smaller than a structureless function of the same
size, not the 3.5x a comparison against the dictionary suggests. (A separate
node-over-support figure of 0.440 on a random support is not the matched null and must
not be quoted as one.) Correction I1 of the review corrects a mis-denominated
retained-byte row in the submission itself, and separately the coordinator's relay of
it: like for like the retained fractions are 0.423 / 0.736 / 0.868 at n_exp = 1..3
[ESTIMATE, from a synthetic dictionary validated to 1.25x on bytes per term], with the
rigorous bracket 0.181–0.886 / 0.322–1.260 / 0.376–1.440.

**Bytes** (same submission; two instruments, kept apart). Peak-to-peak the DD route
costs 11.69x / 40.69x / 36.68x the dictionary route at n_exp = 1..3, from tracemalloc
traced peaks under the default collector. Its never-collected node table alone is
2.77x / 14.92x / 19.90x the dictionary route's whole traced peak, from the gc-forced
retained instrument. The collector matters; three readings at n_exp = 2, paired like
with like: 795,750,300 B peak with the collector off against 234,670,240 B with it
forced (3.4x), both referee-dd-asym-1's; current against current is 9.05x; and the
author's default-collector peak at the same point is 239,740,820 B. So the node-count
advantage does not survive as bytes in this implementation.

**Promotion** (task `T163d291c2cdd4cd0`, submission `S4279919efc104d37`, accepted
`V0c95e8735e50477a`, 11/11 checks) produced a proposed `lab` module and experiment that
were never integrated. Its registered prediction P5 PASSED on its two registered points
(0.281 then 0.518, recorded "deteriorates=True"); what the later six-point series
contradicts is P5's trend clause, and the reason for non-integration is the
coordinator's integration hold (`Mff82e0d959c0493c` item 1, decision
`M889892e48cad42b9`), not a refuted prediction. That round also established, and these
are its findings' only ledger home:

* **C103's variable order does not transfer.** At n_exp = 2 the exponent-first order
  gives a ratio of 0.920 against index-ascending 0.518 on the same object (P6), so the
  order that makes C103's bound work is not the order this store uses.
* **Order invariance where the corrected prefix-saturation condition licenses it.** The
  condition, in the form the earlier review corrected (a worker's first version had an
  invalid induction): if the level widths saturate, w_i = 2^i for i = 0..m inclusive,
  then permuting positions 0..m-1 changes nothing. At n_exp = 1 the profile gives
  w_4 = 16 but w_5 = 31 ≠ 32, so m = 4 and positions 0..3 are licensed: all 24
  permutations of them give total 1,581 nodes with byte-identical level profiles
  (exhaustive). An order permuting positions 0..7, which the condition does not license
  there, gives 1,585 — the must-fail control, earned.
* **Guard and mutants.** The module and `perm_pps` both refuse the same 4/4 malformed
  inputs, and five derived mutants are detected 40/40 each with 0/40 on the unmutated
  run. In the accepted version only Claim 1's 300/300 is harness-gated; the
  gather-equals-scatter 300/300 comes from the pilot's print-only check and superseded
  submissions.

The pilot's validator was print-only and exited 0 regardless (`M841349bedd124e52`);
nothing print-only may be promoted.

## The turn: C103's crossover refuted as its mechanism

C103 puts the first non-trivial exponent level of the *final function's* ROBDD at
k = α + μ(β) + 1, so the first t with any non-trivial level is α + μ(β) + 2. That
equals 4 for N=7, a=3, which is where the node-count ratio turns. Two coordinator
experiments (`experiments/experiment_dd_ordering.py`,
`experiments/experiment_dd_alpha_replicate.py`; rows and verdict logs under
`out/dd-ordering/`, gitignored) tested the prediction across bases. **Measured by the
coordinator and never independently reviewed**, unlike the two board rounds above:

| family | α | predicted turn | measured turn |
|---|---|---|---|
| N=7, a=2 | 0 | 3 | 5 (0.3067, 0.3376, 0.3819, 0.5119, 0.5488, 0.4567, 0.3442) |
| N=7, a=3 | 1 | 4 | 4 |
| N=7, a=4 | 0 | 3 | 5 |
| N=7, a=5 | 1 | 4 | 4 |
| N=5, a=2 | 2 | 5 | degenerate, see below |

So the prediction fails on two bases, and not by an offset: it says a=2 turns *earlier*
than a=3 and the measurement says later. Where the turn sits can be restated as the
point at which the diagram's growth per exponent bit falls through the dictionary's
steady 2.00x, but that is a restatement of where a ratio's argmax lies, not a mechanism;
no mechanism for the turn is known. Two differences make the transfer unsurprising in
hindsight: C103 bounds the exponent-first ROBDD of the Boolean pullback, while this is
an index-ascending multi-terminal diagram of the Walsh COEFFICIENT function with the
exponent variables at the bottom, and the promotion round measured that the order itself
costs 0.920 against 0.518 at n_exp = 2. Note also, when comparing the two texts, that
C103 indexes exponent levels by k while this note indexes by n_exp = k + 1; that is a
convention difference, not a reason the transfer fails. Grouping by α fits the four
non-degenerate families, but the α-partners at N=7 are necessarily inverse bases
(2·4 ≡ 3·5 ≡ 1 mod 7), so that grouping cannot be separated from inversion symmetry;
the test would need same-order non-inverse NON-DEGENERATE bases. Other small moduli have
same-order non-inverse bases (N=8 and N=12) but only at order 2, which forces c_k = 1
for k ≥ 1 and reproduces the N=5 degeneracy below; the first usable modulus is N=11.

**N=5, a=2 is degenerate.** c_k = a^(2^k) mod N is 1 for k ≥ 2, so those blocks apply
C23's involution V rather than growing the support (they are NOT identities: at N=7 the
multiplier-1 block moves 7,168 of 8,192 work states and squares to the identity, and at
N=11 it moves 61,440 of 65,536 (both measured, `V947329feabd540ad`); nothing measures it
at N=5). The dictionary peak freezes at
48,972 from n_exp = 3, so its ratio cannot turn. A must-fail control caught this.
An exact GPU route (dense replay plus integer WHT, `out/dd-ordering/gpu_saturation.py`)
confirms it independently: the final Walsh support is pinned at 32,143 for
n_exp = 3..14, while the non-degenerate families approach exactly half the index space
(0.4938 and 0.4960 of 2^27 at q = 27).

## Withdrawn along the way (all the coordinator's)

* "0.281 is a DD win": random ±1 functions on a RANDOM support of the same cardinality
  need 3,606–3,633 nodes at n = 14, i.e. a node-over-support figure of 0.440–0.443
(`Mb64ff4535cea4767`). That is the
  free baseline, not the matched null; the matched-null ratios above (0.351–0.667) are
  what carry "about 2x".
* "The ratio deteriorates": refuted as a trend by the six-point series.
* An over-retraction: the coordinator withdrew "it peaks at n_exp = 4 in both
  conventions", which its own rows support (the weighted argmax is 0.4399 at n_exp = 4).
  What was actually withdrawn (`Mff82e0d959c0493c` 2b) is the rise-then-fall SHAPE in the
  weighted convention, which dips at n_exp = 3.
* A live-diagram byte measurement (~8 B per node) that measured pointers, and a
  tracemalloc window that charged only allocations made inside it.
* "c_k = 1 blocks are identity permutations" (they are the involution V).
* A relayed claim that ord(3 mod 7) = 6 governs the block sequence: the multipliers are
  3, 2, 4, 2, 4, 2 with period 2 in the tail.

## Where this leaves TODO 13 and the direction

TODO 13 asked for a third simulation method keying on r = β·2^α. C102/C103 answered it
for the canonical diagram of the final function. This note answers it negatively for the
propagation-native diagram: the one feature that looked like an arithmetic signature is
not one. TODO 13 stays open for its original target, a real DD simulator's state
diagram, which nothing here touches. For the peak-memory goal the route is superseded:
under a materialized output no method gains more than about 1.5x (note SL), and under a
streamed output C104's exponent slices remove the peak without diagrams.
