---
code: OB
date: 2026-09-11
title: "Certified finite-work sampling extends beyond the binary scalar tail"
outcome: confirmed
claims: [C56, C60, C61, C62, C63, C64, C65]
todo: [14, 24, 25, 26]
---
# OB — Observable fine-work dynamics need not break certification

TODO 25's bounded implementation/comparison is complete. C61 owns the extended
exact-input promise, C62 the dimension-aware norm transport, C63 the accepted
law, and C64 the forward/backward mass budget. This extends an existing proved
construction; it is not a new simulation principle or a factoring breakthrough.

The qsim-research skill kept the comparison on a fixture that fails the binary
scalar shortcut, required separate coordinate/CDF/matrix dimension accounting,
and led main to repair verifier predicates rather than accept unanimous pass
reports at face value. Lower-cost agents supplied full-law/edge tests and a
read-only proof audit; main implemented production changes, strengthened the
tests and reran them. Default b=2 behavior remains available and tested.

## Implementation boundary

`VerifiedReflectionCircuit(..., block_size=3)` is opt-in. Angles remain exact
rational multiples of pi; backgrounds rotate fine labels 0/1 and leave label
2 unchanged. This is a supplied indexed orbit promise, not a new physical
three-qubit gate interface or arbitrary rounded-unitary input. The b=2 scalar
sampler explicitly rejects b=3. The optional arithmetic backend is unchanged.

Main generalized vector lengths, twisted-shift divisors, exact label updates,
proposal snapshots and coherent acceptance sums. The existing finite
contraction and forward-state/backward-effect identities are reused. The
amplitude planner charges all six real coordinates and three-bin categorical
updates; norm mode uses the corresponding conservative integer radius factor.
The finite-work tree remains binary, so its probability error budget does not
change merely because the work block grows. Matrix costs still do.

## Complete-law and edge evidence

The primary fixture is the frozen r=9,b=3,t=4 circuit from §FW. For target
TV 1e-3 and 1e-6, the experiment covers both prefix and finite-work proposals,
both rectangular and norm enclosure modes, all four histories, three initial
sectors and all sixteen outputs. Each component law and complete accepted law
normalizes exactly. Per-component P192/P256 outward TV bounds meet their own
plans; proposal TV, acceptance success and final TV meet the common plan.
Independent full-r products are normalized and agree with ideal interval
midpoints within the declared 2e-14 tolerance. That independent floating
comparison is not itself the numerical certificate.

At target 1e-6, the accepted-law TV upper bounds are about 1.844e-8 for the
prefix proposal and 1.808e-8 for finite work. Both enclosure modes satisfy
the contract; matching reported summary numbers is not proof that their
entire approximate laws are bit-for-bit identical. The comparison also gives
an outward late-background-removal TV lower bound about 0.05095, confirming
that its non-scalar work dynamics are observable.

The edge experiment checks invalid dimensions/divisibility/rounded angles and
scalar specialization, exact route endpoints, width zero, and complete forced
laws. An insertion-0 angle pi*2^-128 supplies a strictly positive tiny prefix
coordinate mass, with an outward upper bound below 2^-250; the known prefix
normalization factor is explicitly checked. This is a coordinate mass, NOT
proof of an exact-zero binary output branch. No all-zero binary block was
observed at the requested full-law accuracy in that fixture. A separate p=0
mass-oracle test rounds an exact half-half root to zero weights and checks
the declared finite-bit fallback; it does not claim that coarse test meets
the primary target accuracy.

C56's early span-preserving cancellation is tested on its exponent MARGINAL,
not silently strengthened to a joint-sector assertion. Moving the same fixed
rotation to the next insertion supplies a visible reference control. Both
compared laws normalize. The fixed-initial verified output laws differ as
well, which is distinct from a single binary conditional ratio.

For actual adaptive replay, main widens a valid backward effect only AFTER
selecting a bit. The production cursor rebuilds twice in total and replays
that selected bit once; the resulting weights exactly equal a fresh run at
the final precision. Resetting E to identity at the SAME depth/output prefix
changes its normalized binary law by about 0.2200. This is a genuine
conditional-ratio control, unlike §FW's binary mass-only example. All six
real/imaginary coordinates are exercised in the enclosure comparison, with
endpoints extracted at their evaluation precision. Higher-precision midpoint
containment is labeled a diagnostic, not a new containment theorem.

Authoritative main-reviewed reports:

| Experiment | Checks | Report | Log |
|---|---:|---|---|
| `experiment_verified_odd_block` | 4/4 | `out/verified_odd_block_20260911T051252893314Z.json` | `out/verified_odd_block_audited.log` |
| `experiment_odd_block_edges` | 6/6 | `out/odd_block_edges_20260911T051833615443Z.json` | `out/odd_block_edges_audited.log` |
| `experiment_odd_block_comparison` | 7/7 | `out/odd_block_comparison_20260911T051014870500Z.json` | `out/odd_block_comparison.log` |

## Bounded matched timing

Only after the tiny full-law gates pass, the comparison runs a SEPARATE
supplied r=3*(2^60-1),b=3,t=63 circuit with work rotations at insertions
16,32,48 and two coherent reflections at 21,42. This is not a width-scaling
fit or evidence that the tiny visibility magnitude persists unchanged.

At target TV 1e-6, median whole-construction-and-sampling times are about
0.1152 seconds for prefix-component rejection and 0.01480 seconds for
finite-work rejection, a ratio about 7.78. The two methods consume different
random words: their three seeded attempt counts are respectively 4,2,1 and
2,1,2, each repeated twice with alternating execution order. There are three
random traces per method, not six independent traces. This is not a universal
speedup or an optimality claim against all tensor/history decompositions.

All setup, rejected attempts, CDF words and refinement are charged. The wide
finite-work runs rebuild their forward matrices at higher precision before
the first bit, so their recorded replay counts are zero; the edge experiment
separately exercises replay AFTER a selected bit. Maximum working precisions
in these runs are 202 for prefix rejection and 212 for finite work: higher
precision does not imply more total work when repeated contractions disappear.
Retained scalar counts are not RSS; no wide orbit/output table is allocated.
The backend-wide primitive-error constant remains deferred in TODO 24.

## Preserved failures and corrections

- The full-law agent first supplied an incorrectly shaped float comparison
  gate; `out/verified_odd_block_20260911T051005477475Z.json` records the rejected
  input. Its repaired report `out/verified_odd_block_20260911T051025646760Z.json`
  remains available. Main then added the omitted independent TARGET-error
  predicate, guarded accumulated report rows, and used the standard harness
  schema. Per-component agreement alone would not check the coherent target.
- The first edge invocation omitted the optional backend; its failure is
  `out/odd_block_edges_failure_20260911T050958854802Z.json`. A subsequent
  false invalid-mask control and insufficient replay predicate are retained
  in `out/odd_block_edges_20260911T051031155739Z.json` and the matching
  `out/odd_block_edges_failure_20260911T051031160350Z.json`.
- Further agent reports at `051356958085Z`, `051445451623Z` and
  `051453674461Z` under the `out/odd_block_edges_` prefix are preserved.
  The final agent report `out/odd_block_edges_20260911T051528897662Z.json`
  explicitly disclosed that it did not trigger production replay. Its
  reset-effect comparison also used a root versus a later node. Main did not
  treat this as the requested same-node adaptive test; the authoritative
  correction above now exercises both actual replay and the correct control.
- Main also checked exact route destinations rather than only determinism,
  aligned the near-zero forced circuit with its tiny-coordinate diagnostic,
  tested the coarse fallback separately, required all reference masses to
  normalize, and selected a nontrivial six-coordinate enclosure label.
  These are verifier corrections, not evidence of a failed production theorem.

## Validation and follow-up

Fresh core passed before science: `out/odd_verified_core.log`. Affected
suites pass in `out/odd_verified_lab.log` and `out/odd_verified_claims.log`.
Without the optional backend, lab passes with explicit arithmetic skips in
`out/odd_verified_lab_no_flint.log`. Existing binary complete-law regressions
pass in `out/odd_verified_sampling_regression.log`,
`out/odd_verified_finite_work_regression.log` and
`out/odd_verified_rejection_regression.log`. The other six science suites
were not rerun. The documentation gate is `out/odd_verified_docs.log`.
Neither paper nor either abstract workshop was edited; no commit was made.

Reproduce, substituting any experiment name from the table:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_odd_block_comparison
```

An initial algebra audit suggests a connection between forward channel memory
loss and true-mass-weighted sampling error, without discarding backward
measurement effects. Main and a lower-cost agent checked the candidate POVM
argument; no mixing/compression experiment or new sampler is supplied here.
The primary-source check found directly relevant product-contraction/MPS work:
main read Definition 2.1, Propositions 2.2/2.4 and the good-block condition
and proof of Proposition 3.9 in
[Pathirana, arXiv:2605.00157v1](https://arxiv.org/html/2605.00157v1).
It bounds memory loss using contractive channel products, including cases
where a one-step criterion is insufficient. This is relevant prior work,
not our discovery of a general replacement principle. Its random-cocycle
and MPS application sections have not been audited here. TODO 26 owns the
specific next question and its required counterexamples.
