---
code: FM
date: 2026-09-11
title: "Forward-state approximation can preserve the output law without resetting measurement memory"
outcome: confirmed
claims: [C56, C64, C65, C66]
todo: [14, 26, 27]
---
# FM — Forgetting one side of the contraction is not forgetting both

C66 owns the mathematical approximation contract, its proof, the finite
arithmetic budget and the conditional mixing consequence. This note records
the experiment and its limits. No production sampler or manuscript changed.

The qsim-research skill required an observable odd-block fixture, full output
laws, a stronger binary null and false shortcuts that genuinely fail. That
kept a fixed-point observation from being mistaken for decoupling. Lower-cost
agents audited the algebra and tested channel products; main implemented the
full-law diagnostic, repaired verifier issues, and reran the evidence.

## Complete-law experiment

`experiment_forward_memory.py` fixes §OB's r9,b3,t4 input and varies only
epsilon in rho_tilde_i=(1-epsilon)rho_i+epsilon I/3. It checks all four route
histories and three initial sectors at P192 and P256. There are 96 complete
laws across the four epsilon values, not sampled histograms. Their finite
CDF probabilities are exact Fractions, with p=L=80. The tiny experiment
retains its whole diagnostic tree; it does not claim to save memory itself.

Every branch comes from the existing verified finite-work builder. Every
backward child comes from the shared QFT branch helper. At each depth, summing
all child effects contains the identity, including feedback phases. Forward
states remain normalized PSD by construction, not by clipping eigenvalues.
The nonzero-epsilon rows separately check aggregate effect-weighted mass
errors against outward Frobenius-to-trace bounds. Ideal conditional output
intervals come from the distinct shared prefix contraction; bijective routes
allow multiplying a joint final-sector cell by the number of sectors.

The complete finite laws normalize exactly and obey the state-plus-numerical
budget. The independent full-r floating reference agrees within about
4.2e-17; this diagnostic is not the interval certificate. At epsilon=1/16,
maximum outward output TV is about 0.00836 against a maximum bound about
0.166. At epsilon=0, the integer kernels match the production cursor with
matched p,L, and maximum outward TV is about 2.5e-24. No precision/gap scaling
claim follows from these two tested precisions.

Controls matter here:

- Resetting E on the original pure-input fixture changes its full law by
  outward TV at least about 0.562. Both compared chains normalize; this is
  not a mass-loss artifact or a different-depth comparison.
- A separately specified mixed initial I/3 is exactly preserved by its
  unital work channels. The resulting t=2 law nevertheless has correlations:
  the analytic reference and finite calculation agree exactly, while resetting
  E gives TV 1/8. This mixed-input diagnostic tests C56's general identity,
  not a new initial-state option in the production exact-input API.
- A five-dimensional diagonal POVM demonstrates why unrelated per-prefix
  states cannot use one maximum state-error budget. Its exact arithmetic is
  regressed in `test_claims.py`. This is explicitly not claimed to be a
  realization of our three-state QFT circuit.
- The binary late-gate null gives zero change under the same PSD mixture,
  confirming that it would hide sensitivity in the main experiment.

Authoritative report:
`out/forward_memory_20260911T053532681815Z.json`, 7/7 checks;
log `out/forward_memory_audited.log`.

## Initial mixing probe and its scope

The separate repeated-gate fixture uses r9,b3, no routes and the same embedded
Rx(pi/7) after each of sixteen controls. Only the channel-product length
varies. All three initial sectors are checked. The SVD is explicitly an
uncertified diagnostic of the eight-dimensional traceless Hermitian channel
representation, not an interval trace-norm certificate.

Each one-step Hilbert–Schmidt norm is numerically one; the length-sixteen
products have norms about 0.599 in sector zero and 0.595 in the other sectors.
The no-background control retains an explicit nonzero conserved mode and
norm one. The implication from no single-step contraction to no product
contraction fails. This is a familiar channel-product mechanism, not a newly
discovered general principle; C66 gives the primary-source positioning.

Main reused `_density_forward_step`, added branch-unitarity checks, rebuilt
branches at P192, corrected a minimum-versus-maximum eigenvalue check in the
null, preflighted aggregate live storage and retained the harness schema.
Authoritative report: `out/forward_mixing_probe_20260911T053605193970Z.json`,
3/3 checks; log `out/forward_mixing_probe_audited.log`.

## A certified finite-order rate, with a poor worst-case warm-up

`experiment_forward_mixing_certificate.py` evaluates an exact orthonormal
traceless Hermitian basis with Acb enclosures. For each of three r9 sectors
and six starting offsets, it forms a twelve-channel product. All eighteen
outward Frobenius bounds are below 4/5 at both P192 and P256; the reported
bounds range from about 0.7643 to 0.7659. This is an induced
Hilbert–Schmidt upper bound, NOT the trace-norm coefficient. The phase
identity and periodicity turn the finite certificate into a repeated-block
rate for that supplied r9 family; they do not extend the numerical constant
to untested sector phases.

Main audited the complete code, added invariant-subspace predicates and
zero-containment of the preserved mode, and corrected the phase identity's
denominator to M rather than b (equal only on this fixture). No-background
spectral projector differences supply an explicit conserved nonzero mode.
All branch-unitarity and basis checks pass. The conserved null mode refutes
inferring contraction merely from a maximally mixed fixed point.

Authoritative report:
`out/forward_mixing_certificate_20260911T054348729950Z.json`, 4/4 checks;
log `out/forward_mixing_certificate_audited.log`. The initial successful agent
report is `out/forward_mixing_certificate_20260911T054211166701Z.json`.

The sufficient bound is not yet practical. Using 3/2 as a rational upper
bound on the pure-state factor sqrt(2), the rate 4/5 per twelve controls
gives sufficient warm-ups of 624 and 996 checkpoints for approximation error
targets 1e-3 and 1e-6 respectively. These are conservative sufficient counts,
not observed mixing times or lower bounds. Both exceed the production API's
width limit. There is no measured memory saving, and raising that limit just
to display one is not the next priority.

## Why the next step is no longer an exceptional-sector search

While auditing the bounded product probe, main noticed that scalar arithmetic
phases cancel in FORWARD conjugation. Its schedule has period two at every
sector phase, even when raw branch matrices have a longer period. Main then
derived a norm-equality/eigendirection argument excluding a preserved
traceless Hermitian mode for either two-step product. A lower-cost agent
independently audited it. C66 owns the resulting qualitative uniform-gap
proof. Compactness gives a gap for every sector, not its numerical size.

This is the strongest structural outcome of this investigation, but it is
still an application of established channel-contraction geometry. It concerns
one fixed mixer, not arbitrary perturbations or a factoring algorithm.
The exact-input production sampler and its width limit remain unchanged.

## Preserved failures and audit corrections

The first main full-law run failed in its new Frobenius bound, before accepting
any laws: squaring an interval containing zero and taking its square root
could straddle the domain. The repair uses nonnegative outward upper endpoints
before the root, preserving an upper bound. The failed report remains at
`out/forward_memory_20260911T053147246115Z.json` and log at
`out/forward_memory_initial.log`. The intermediate six-check successful report
`out/forward_memory_20260911T053227154844Z.json` precedes the additional binary
null and aggregate-error predicates; it is not the final coverage record.

The initial lower-cost mixing probe had TWO scientific predicate failures
before a JSON serialization failure, not only a reporting problem. It
incorrectly demanded real entries in Hermitian matrices, which may have
imaginary off-diagonals. Correct tests check Hermiticity and real coefficients
in a Hermitian basis. NumPy boolean serialization was also fixed. The initial
log is `out/forward_mixing_probe_test_20260911T053119Z.log`; the agent's repaired
report is `out/forward_mixing_probe_20260911T053201074130Z.json`. Its repaired
log was accidentally placed outside `out/`, at
`/home/djneko/Workspace/qsim-test/20260911T053200Z.log`; it is left untouched.

The algebra auditor initially omitted the QFT feedback phase in a proposed
two-control reset example. Main corrected it before implementation: the full
law is the four-entry one in C66, not two identical conditional pairs. The
independently checked complete-law test uses the corrected feedback phases.

The first certificate attempts had harness/API failures, retained as
`out/forward_mixing_certificate_failure_20260911T053745712825Z.json`,
`out/forward_mixing_certificate_failure_20260911T053755595992Z.json`,
`out/forward_mixing_certificate_failure_20260911T053804355324Z.json`, and
`out/forward_mixing_certificate_failure_20260911T054103623359Z.json`.
They respectively used a complex value as a real endpoint, missed a helper
argument, passed a branch pair where the all-branch audit expected a list,
and called an unsupported inverse method. Main also required shared density
updates, an aggregate allocation allowance and preservation of complete
per-window evidence in the report. No failed contraction prediction was
hidden by changing its window length or threshold.

## Validation and remaining scope

Core passed before science (`out/forward_memory_core.log`). Affected lab and
claims suites pass (`out/forward_memory_lab.log`,
`out/forward_memory_claims.log`), with the optional verified backend present.
The other six science suites were not rerun. No production helper, either
paper or either abstract workshop was changed; no commit was made.
The documentation gate is `out/forward_memory_docs.log`; `git diff --check`
is clean. Enclosure and row-storage allowances are conservative preflight
counts, not measured native memory or RSS.

Reproduce any of the experiments with this environment:

```bash
OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' --with 'python-flint==0.9.0' python -m experiments.experiment_forward_memory
```

The current next question belongs only in TODO 27. A law bound is useful
progress, not a simulation breakthrough: current APIs still store their full
forward history, and general numerical bit costs remain in TODO 24.
