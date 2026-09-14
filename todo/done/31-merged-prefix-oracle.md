---
id: 31
state: done
title: "Can coherent merging supply every gate-by-gate prefix without rejection?"
outcome: "C72 proves and implements all merged prefixes with no-rejection sampling; complete tiny laws pass, with a runtime/storage tradeoff against C71"
claims: [C56, C59, C60, C61, C68, C70, C71, C72]
---
# Use the small backward walk as an amplitude oracle

**Completed 2026-09-11 at mathematical/float scope.** C72 owns the proof,
batched memory/cost contract and implementation; §MP owns the all-prefix,
complete-law and matched comparison evidence, including corrected controls.
The no-rejection sampler is not a new finite-TV certificate, and its batching
cost does not give a universal runtime win over C71. TODO 32 owns the next
fixed-route-alphabet question. No default or manuscript change.

## Original plan (historical)

TODO 30 / C71 / §CM complete sparse coherent reverse rejection at
mathematical/float scope. The next stronger simple baseline is to reuse the
same merging in C59's EXISTING no-rejection gate-by-gate sampler. Removing
history enumeration need not require paying C71's b*S rejection envelope.
This is an oracle substitution in an established sampler, not a new generic
propagator or a claim that no-rejection is automatically faster.

## Initial algebra audit, not yet an implementation

For every allowed C59 `prefix_vector` request, reverse-contract its actual
prefix from each target fine boundary j and coherently sum the initial
coarse components. Conjugate to match the forward-amplitude convention and
return sqrt(M) times the b-vector, preserving all relative phases. Reuse
shared shifts, work/reflection gates and QFT factors. Batch columns if useful;
do not create a second independent implementation of circuit semantics.

Fixed unmeasured exponent bit x_i contributes T^(-2^i*x_i)/sqrt(2).
A measured bit contributes [I+(-1)^z*conj(phi_i)*T^(-2^i)]/2. In reverse
order, apply any included R_(i+1)^dag, then W_(i+1)^dag, before the control
factor. Unvisited prepared controls contribute 2^(-(t-stop)/2).

Match `_contract_route_component` exactly at EVERY boundary: arithmetic
excludes insertion-stop W/R, background includes W but not R, reflection
includes both. At stop=0, the adjoint order is W0^dag R0^dag when both
are included (apply R0^dag first). The all-prepared initial state is not
the same prefix as one with insertion-0 gates already applied.

The C71 bound applies to the reflections actually included in each prefix,
including intermediate peak support. Naively computing all b fine columns
costs O((stop+1)*S*b^3+k_prefix*S*b^2) per b-vector query. A complete BGL
draw then has a candidate O((t+k+1)*(t+1)*S*b^3) arithmetic bound. Compare
with C71's expected O(b*S*((t+1)*S*b^2+k*S*b)), including setup, dictionary
work, all queries/attempts, matrices and bit costs. The ratio suggests a
tradeoff in S versus t+k, not a universal winner. Audit before promotion.

## Smallest discriminating experiment

Use lower-cost agents for independent initial prefix/law tests while main
implements and audits. Start from the harness template and run core first.

1. Freeze a tiny b=3 circuit with noncommuting W, initial/terminal reflections,
   and collisions. Enumerate ALL valid stop/boundary/exponent/measured/output
   labels within a pre-budgeted cap. Compare COMPLEX vectors against C59 and
   existing independent full-r prefix helpers. Cover t=0/1, fixed points,
   exact zero amplitudes and generic-label fallback.
2. Feed the audited oracle into the existing block sampler. Enumerate the
   full tiny transition law, not just selected prefixes or a sample histogram.
   Forbid `_histories` and verify zero rejection, exact query counts and
   claimed peak/storage bounds. Count actual internal sparse work.
3. Must-fail controls: wrong boundary inclusion, incoherent merge, discarded
   phase or tail normalization. The controls must exercise intermediate
   prefixes; final-output agreement alone cannot validate this oracle.
4. Only after correctness, compare at fixed b,t,r and frozen route labels
   while varying k within the existing cap. Use BOTH adjacent and separated
   coprime routes as separately frozen fixtures from §CM. Include C71 reverse,
   C59 history methods and cheap static/input simplifications where applicable.
   Charge returned samples and setup; separate float diagnostics from any
   requested-TV comparison. Do not infer hardness from a huge exact group.

Consult the actual BGL theorem/algorithm and sparse-amplitude prior work before
claiming a new contribution. Reuse C60's oracle budget only after verifying
its assumptions for the new coordinate oracle; C70's different certificate
does not transfer by renaming the state. Finite-bit certification, broader
route alphabets and backend-wide precision complexity are deferred, not
implicitly solved. Avoid another large synthetic size record or manuscript
change until this baseline question is settled.
