---
code: GS
date: 2026-09-11
title: "A cheap physical phase couples every coarse sector, but this is not sampling hardness"
outcome: confirmed
claims: [C53, C56, C58, C74, C75, C76]
todo: [14, 34, 35, 36]
---
# GS — Additive characters test the representation, not just gate complexity

C76 owns the physical rotation identity, exact prime-field nonzero proof,
one-initial-sector truncation bound and their scope. This note records the
discovery and finite evidence. The prime-field Gauss magnitudes are known
mathematics, not a numerical discovery or a new phase-estimation algorithm.
Main derived the subgroup and cyclotomic arguments; lower-cost agents supplied
initial kernel/output tests and independent proof audits. Main strengthened
the verifier and reran all three experiments using existing Circuit/statevec
and sequential_path. No shared production helper or manuscript changed.

## Frozen physical and algebraic checks

The elementary kernel report is
`out/gauss_kernels_20260911T093101284276Z.json` (4/4 PASS). It checks the
previously frozen primes/bases (7,3), (13,2), (17,3), zero/nonzero phases,
and every complex column of the binary phase gate, including its global
phase. Physical compositions G_1 G_-1 and G_1 squared agree with identity
and G_2 respectively. Removing the kernel phases fails the inverse identity
with error .952855. Simple coefficient magnitudes are not classical transition
probabilities for a coherent circuit. The separately projected N13,b3 subgroup
matrix agrees with the direct formula; the full-orbit flat-magnitude law was
not silently applied to a subgroup.

All 1,496 direct summands reconcile with the preflight. The physical column
checks charge 19,200 gate-entry updates, distinct from gate-column applications.
The simultaneous numerical payload bound is 409,600 bytes under 16 MiB.
These counters are neither native FLOPs nor bit complexity nor total RSS;
the tiny diagnostic still enumerates the orbit and reference matrices.

The exact structure report is
`out/additive_phase_structure_20260911T093102955518Z.json` (6/6 PASS).
For those same primes, k=1,2 and all divisors b<=6, all 264 subgroup
coefficients are nonzero by integer cyclotomic-polynomial reduction. This
is independent of the floating cutoff used to display magnitudes. The
Gauss decomposition agrees within 6.8e-16; one-sector mass and every tested
best-L truncation bound pass for individual p inputs and two fine-coordinate
weightings. Different fine amplitudes are orthogonal, so their phases cannot
cancel sector probabilities for a one-sector input.

The composite N9,a2,b2,r6,k1 control has exact zeros: p=0 retains frequency
2, p=1 frequency 1. This is one frequency PER p, not a single block-wide
route for a superposition of fine labels. The report's shorter control
description should be read with that qualification. Polynomial input slots
(24,718), reductions (270) and floating summands (10,744) reconcile exactly;
the maximum root order is 272. The numeric payload bound is again 409,600
bytes. Integer polynomial coefficient storage is separately capped by input
slots/root degree; these are not a bound on backend bit-runtime or RSS.

## Actual exponent output and strongest tested shortcuts

`out/additive_phase_output_20260911T093102634931Z.json` passes 4/4 checks.
It freezes C75's N13,a2,b3,r12,t3 initial/post-control mixer schedule, removes
the earlier cell reflections, and inserts G_k after the second post-control
mixer. Only k varies over 0,1,2. A fourth row places G_1 at the end as a
trace-out null. The physical register contains residues, with the same ten
qubits and clean unused flag as the earlier fixture; this is not an indexed
gate being described as a physical implementation.

The complete eight-output compiled law agrees with independently assembled
full-r branch matrices and sequential_path to at most 7.683e-14 per outcome.
Maximum scratch leakage is below 3.73e-30 and norm error below 3.75e-13.
Omitting the interior phase changes the law by TV .0220140067 (k1) and
.0355523339 (k2). Initial coarse dephasing incurs TV .0176795970 and
.0195881265 respectively. Zero kick and terminal insertion are nulls, and
coarse dephasing works for the no-kick background. These are meaningful
failures of two particular exact shortcuts, not of every approximation.
Indeed, an allowed TV of .05 would already tolerate omission on these rows.
The ideal-order-finding law is a DIFFERENT baseline: its TV from the no-kick
mixed-work background is .2086007086. Neither is assumed to equal the other.

Each physical circuit uses 15,213 Pauli rotations; the four runs jointly
charge 62,312,448 gate-entry updates under the original 100,000,000 cap.
The independent reference makes 144 sequential_path calls and charges
746,496 units of t*d^3, with reference payload bound 185,856 bytes. It uses
finite full-r arrays, not an orbit-free algorithm. Some final diagnostic runs
overlapped briefly; no matched timing or host-health conclusion is drawn.

The structure report also directly tests C53's scalar phase-transfer premise
on the four actual low-control histories. The best complex scalar leaves
relative work-vector residuals between .243 and .921 for the eight nonzero
kick/history pairs. The k0 cases are identity NULLS, not evidence for a
nontrivial transferable phase. This excludes diagonal scalar replacement at
that insertion on these histories, not a general early-exponent unitary or
an output-specific rewriting. TODO 36 explicitly tests that stronger baseline.

## Proof audit and primary context

After the direct projection, main noticed a sharper exact fact than dense
floating support: a vanished coefficient would contradict the minimal
polynomial of a prime root of unity over the coprime cyclotomic field.
Independent lower-cost audits checked the proof, signs, subgroup decomposition
and one-sector-only probability bound. C76 states the proof and primary
sources; it does not claim novelty for this elementary application.

Main read van Dam/Seroussi's finite-field definitions and Gauss norm facts
with their proof in [quant-ph/0207131](https://arxiv.org/pdf/quant-ph/0207131),
and Evans's irreducibility theorem/proof and coprime compositum argument in
[the cyclotomic notes](https://maths.dur.ac.uk/users/daniel.evans/GaloisTheory/Notes/cyclotomic-extensions.html).
The paper's Gauss-phase/discrete-log oracle reduction is not a lower bound
for this fixed-output task. C74 separately records the earlier reading of
its powering reduction. No assertion relies on the unread quantum-algorithm
or finite-ring sections.

## Preserved failures and verifier corrections

The first Gauss attempt failed an API guard on k=N-1, outside the physical
builder's frozen signed range, after earlier b1/b3 measurements had run:
`out/gauss_kernels_20260911T092044205720Z.json`. Using the equivalent k=-1
within the documented gate contract fixed the call; the failure JSON remains.
The initial log was subsequently reused, so it is not the failed-run log.
The corrected lower-cost report is
`out/gauss_kernels_20260911T092201282162Z.json`. Main then distinguished
gate-column applications from gate-entry updates and checked their exact
total rather than accepting a nonnegative cost counter.

The lower-cost output and structure reports are
`out/additive_phase_output_20260911T092353002435Z.json` and
`out/additive_phase_structure_20260911T092840301919Z.json`. Before the first
output execution, the builder signature and reference-counter reset were
repaired; no failed-run artifact exists for those pre-run defects. Main then
removed import-time experiment logging, corrected a false orbit-free-storage
description, and added physical gate/layout/cumulative-cost predicates.
The main reports above include those fixes. Tolerances, fixtures and controls
were not changed to select a favorable scientific outcome.

## Validation and continuation

Core passed first in `out/additive_phase_core.log`. The affected claims suite
passes in `out/additive_phase_claims.log`, including independent tiny monic
division regressions over Gaussian integers and integers for C76's prime and
composite controls, plus a one-sector Gauss mass check. The other seven science
suites were not rerun; shared production helpers are unchanged. System Python
3.12.3, NumPy 2.4.6 and single-threaded BLAS were used; the exact structure and
claim checks also use python-flint 0.9.0. These small passes do not settle
TODO 34's runtime/host problem. No host configuration, firmware or long timing
benchmark was changed or launched. Documentation validation is recorded in
`out/additive_phase_docs.log` after index regeneration.

Reproduce from research/:

```bash
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_gauss_kernels
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' python -u -m experiments.experiment_additive_phase_output
PYTHONFAULTHANDLER=1 OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12 --with 'numpy==2.4.6' --with 'python-flint==0.9.0' python -u -m experiments.experiment_additive_phase_structure
```

The skill's stronger-baseline rule prompted exact support and actual-output
checks before interpreting dense sectors. C57's fixed-history support cone
does not by itself remove the coherent history sum, and its localized-kick
omission bound is trivial for an extensive additive phase. The next bounded
question, including the stronger exponent-unitary transfer test, lives only
in TODO 36. Neither a hardness claim nor a new efficient sampler follows;
the open-ended research goal remains active.
