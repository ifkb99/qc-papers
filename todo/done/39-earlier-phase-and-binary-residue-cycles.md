---
id: 39
state: done
title: "One earlier phase: complementary progression covers and sampling"
outcome: "C79/ER implement and validate complementary earlier-phase covers; C80/CW add an opt-in root-mass proposal with independent law/RNG checks and charged bounded comparison"
claims: [C57, C76, C77, C78, C79, C80]
---
# Test a real schedule boundary, not a larger known-order size record

**Coefficient and complete-output checkpoints complete:** C79 owns the
proofs, including cancellation revival, the complementary remaining-history
cover and the two-history coherence identity. ER owns independently checked
output comparisons, controls, audit failures and resource limits. The new
opt-in EarlierPhaseProgressions helper implements this schedule; the original
C78 class still does not. The derivations/test plans below are retained as
the pre-measurement record.

**Sampler checkpoint complete:** the bounded opt-in prototype now compares
the same joint/output law and accuracy with C78 omission, feedback-aware
dephasing and C79's exact v=0 phase-into-W0 reduction. ER owns the final
evidence and limits. A full-r sequential method
remains a relevant width-linear baseline; no timing superiority was tested.
The omitted-phase law is not the strongest baseline at the initial endpoint.

**Weighted-proposal checkpoint complete:** C80 owns the envelope/cost-model
proof; CW owns the independently reproduced formula, actual RNG and physical
same-law comparison, including failed verifier/preflight history. The opt-in
root_mass mode retains the original default. No numerical certificate or
runtime superiority is implied. TODO40 owns the next multiple-diagonal-phase
question. The unimplemented-state wording in the frozen plan below describes
its pre-measurement checkpoint, not the current API.

**Separate deferred question:**
Charge proposal generation AND ratio evaluation if pursuing the separate
C59 factor-two envelope; its cheap acceptance
rate does not make those oracles free. Reuse the existing progression kernel,
retain the sampled work label across retries, and gate the actual RNG path
before a production integration claim.

The completed first width follow-up exposes long progressions rather than another
tiny law whose candidate classes are nearly all singletons. A fixed family
with one remaining low control (v=s-1) keeps K=2 while the cycle description
can grow. Select bounded forced rows/outputs and preflight their references;
do not allocate full Q-by-r matrices merely to exhibit this elementary bound.
The frozen fixture and scoped evidence are retained below and in ER. Numerical
certification, arbitrary multiple earlier mixers and TODO 32/34 host-stress
timings remain separate. No generic factoring or hardness claim is intended.

C78's one-late-insertion sampler is implemented; WF owns its validation.
Its compact rows rely on the diagonal phase being evaluated at an index fixed
by the final work label and the small late history. The original question changed
that premise: add ONE earlier pointwise phase after v low ascending controls,
before the unchanged late mixer/phase at split s. Do not feed this new schedule
to the original LateWorkProgressions class and call its output exact.

## Derivation and frozen prediction, before measurement

Keep e=l+2^s*h, rho=(c+q-u) mod r and each old component's
l=rho+n*r. The extra phase is evaluated at

    g_early((u + (l mod 2^v)) mod r).

Along this component the argument repeats after

    P = 2^v / gcd(r,2^v)

steps in n. Splitting by n mod P therefore restores constant coefficients,
with stride r*P=lcm(r,2^v), at most H*min(r,2*b-1)*P components (truncate
to the actual candidate progression lengths). IMPORTANT: recompute the sum
over initial u/intermediate q inside each cycle class. The early phase may
revive an old exact cancellation, so multiplying the old scalar coefficient
or refining only the old NONZERO components is invalid. This is a sufficient bound, not a claimed
minimal period or hardness lower bound. If 2^v divides r, P=1 and the old
geometric component bound survives even though the insertion is earlier. The endpoint
v=s has an additional simplification: u+l=c+q mod r, so the phase is fixed
locally and constant coefficients survive even when the period bound is loose.
No new general multi-mixer claim or novelty assertion is intended.

Freeze the r60,b3,t8,s7,H2 fixture, W0/W1 and late G_1 unchanged; add an
identical pointwise G_1 at v=0,1,...,7. Only its insertion position varies.
Lower-cost initial testers independently audit the formula and compare all
literal columns/conditional row coefficients against this cyclic refinement.
Use pointwise modular powers, no supplied orbit table in the refined formula.
Bound reference storage at 16 MiB and actual local/root/query work before loops.
Expect P=1 at v<=2 and the independent endpoint identity at v=7. For interior
v=3..6, whether the unrefined coefficients actually vary is an explicit
empirical branch: exhibit a witness or explain cancellation; do not assume
an upper bound is tight. A wrong unrefined constant-row or omitted-phase
control must fail somewhere before calling the probe discriminating.

These are initial coefficient/amplitude tests, not a new production sampler,
full-output comparison or speedup. If they pass, next compare complete output
laws and the strongest approximation before deciding whether the P factor
is useful. Carry/finite-state compression of the cycle is a possible later
angle, not a result. Numerical certification and host-stress timings remain
separate; do not revive TODO 32/34's long runs.

## Frozen complete-output test (2026-09-11, after the coefficient checkpoint)

Retain the SAME t8,s7,r60,b3,H2 matrices and two G_1 phases. Vary only the
earlier insertion v=0..7. Before measurement, predict that the geometric-
progression Fourier formula, a tiny FFT of the literal columns, and selected
existing Circuit/statevec inverse-QFT checks agree as COMPLETE output laws.
Whether any inexpensive candidate reaches TV 1e-3 or 1e-2 is an empirical
branch, not a prediction that it must fail. No numerical certificate/timing
claim or production sampler extension is implied.

Two baselines are required: (i) omit only the added earlier phase, giving
the implemented C78 schedule; (ii) FB-style dephasing of the single high input
history while retaining the early phase and exact Fourier feedback. For the
latter, set chi_l=U^(l-a) D_early U^a W0|0>, a=l mod 2^v. The common late
W1/G and fixed high shift cancel in each diagonal history. Its full law is

    p_deph(y) = ||sum_(l<L) chi_l exp(-2*pi*i*y*l/Q)||^2 / (Q*L).

Equivalently draw the low output bit z uniformly, then Fourier-sample the
length-L state with feedback exp(-2*pi*i*z*l/Q). This is a valid normalized
candidate, not assumed equal to the target. Support offsets lie in [-2,4];
strip a common unit work translation to center the cone at radius 3. The
high-history separation dist_60(128,0)=8>6 then certifies an exactly uniform
first output bit for EVERY target v. This marginal equality does not assert
independence or justify removing the high-history coherences.

Lower-cost initial testers own (1) progression/full-FFT laws, (2) independently
assembled selected Circuit-QFT references, (3) the feedback-aware baseline
and an independent contraction/cancellation audit. Main audits predicates,
phases, normalization, counters and scope. All arrays/temporary state vectors
must be preflighted under 16 MiB BEFORE allocation; budget reference/FFT/root/
QFT work separately, including actual pointwise modular powers. No uncharged
orbit tables or sampler-support scans. Keep failed reports. Must-fail controls
include dropping coherent Fourier cross terms and a genuine missing-feedback
counterexample. Gate the actual output, not only the conditional coefficient.

## Candidate dual cover, frozen before a pilot (2026-09-11)

The remaining low controls suggest a complementary partition. Put A=2^v,
K=L/A and l=a+A*k with a<A and k<K. For fixed (j,h,k), retain the previous
w,c,p and merge local pairs by rho'=(c+q-A*k-u) mod r. The early phase is
now fixed at (c+q-A*k) mod r. The proposed component is

    start=L*h+A*k+rho', stride=r,
    count=max(0,1+floor((A-1-rho')/r)),
    gamma=g(w) sum_(u,q: c+q-A*k-u = rho' mod r)
                    W1[p,q] W0[u,0] f((c+q-A*k) mod r).

Different k occupy disjoint length-A intervals. This suggests a dual bound
H*K*min(r,2*b-1), with local construction O(H*K*b^2), and hence choosing the
smaller of the cycle factor P and remaining-history factor K. The v=s
endpoint would become the K=1 case of a general construction. This is a
candidate elementary regrouping, not a new sampler or a generic polynomial
bound: worst-case min(P,K) can still grow exponentially with width.

First obtain an independent algebra audit. If sound, use the same eight
insertion positions and tiny literal columns to test every reconstructed row,
disjointness and Fourier law, retaining both component counts. Freeze the
prediction that both covers agree, not that one is always smaller. A wrong
phase that drops the A*k shift must fail on an interior insertion. Keep
16 MiB aggregate numeric payload and separately bounded actual local/root
terms before loops. Do not change the production helper during this pilot.

## Frozen sampling implementation test (2026-09-11)

The preceding coefficient/output pilots are complete. Implement an opt-in
EarlierPhaseProgressions specialization of the C78 helper: add an early
insertion/phase, construct either proven cover, and expose its uniform row
stride. Reuse the original work draw, progression primitive and coherent
acceptance loop, preserving the work label across retries. Existing C78
row tuples remain (start,count,coefficient); add a stride field to the row
record without changing defaults. Choose the cover BEFORE matrix copies or
row construction from conservative local-pair-visit bounds, including finite
progression truncation. Explicit cycle/dual modes are diagnostic comparators,
not permission to bypass caps. No prefix/output/orbit arrays in the sampler.

Predict full joint work/output equality with normalized literal-column FFTs
for the established r60,b3,t8,s7,v0..7 fixture. Verify actual RNG flow and
zero rows, cancellation revival, aliases, invalid phases/caps and the new
stride (including its gcd lift). Include the exact v=0 folded-W0 C78 baseline;
omission/dephasing comparisons retain their stated accuracy scope. Dropping
coherent acceptance or resampling work on rejection must fail discriminating
controls. A finite attempt cap must raise without returning a substitute.

For the untruncated follow-up freeze N61,a2,r60,b3, the SAME W0/W1 and two
G1 phases, s=7..12 in unit steps, t=s+1 and early insertion v=s-1. Only s
varies under this explicit schedule rule; H=K=2 remain fixed. Predict the
dual description stays within its fixed component bound while both covers
give the SAME selected complete conditional Fourier laws. Test work labels
0,1,59. Literal references may enumerate the union of their known global
support residues and call the existing tiny column engine; they must charge
those calls and prove omitted entries zero from the support cone. No full
Q-by-r array. Full-r sequential sampling remains a relevant width-linear
baseline: a cycle-versus-dual gap is not a speedup over every classical method.

Main owns the helper and regression integration. Lower-cost initial testers
own independent joint-law/strong-baseline checks, RNG/edge checks, and the
long-progression comparison. Preflight aggregate numeric storage below16MiB
including old/new arrays, FFT magnitudes, scalar lists and reference setup;
record actual cumulative named counters before loops/calls with separate
caps for local/oracle/transform/proposal work. Freeze each executable budget
before its run. No numerical certificate, broad timing or novelty claim.

## Completed bounded test: component masses, not just counts (frozen 2026-09-11)

The sampler/long-progression checkpoint above is complete at float diagnostic
scope. A read-only independent algebra audit confirmed a direct application
of C59's weighted Cauchy-Schwarz argument. For an active disjoint row let
lambda_c=count_c*|gamma_c|^2/Z, B=sum_c sqrt(lambda_c). Select component c
with sqrt(lambda_c)/B and retain the existing normalized progression kernel.
If F_c(y) is its UNNORMALIZED Fourier amplitude, put

    S(y)=sum_c |F_c(y)|^2/sqrt(lambda_c),
    q(y)=S(y)/(Q*Z*B),
    accept(y)=|sum_c F_c(y)|^2/(B*S(y)).

Predict the SAME accepted conditional law, mean attempts B^2<=m, and equality
for equal masses. Ignore exact zero components only; at S=0 both proposal
and target vanish. Keep work fixed. This is an established weighted-envelope
application, not a breakthrough or a precision certificate. The implemented
sampler has NOT yet changed to these proposal weights.

First run a tiny formula/actual-kernel pilot before integration. Freeze Q=8,
component 1 singleton e=0 with mass9/10, component 2 progression e=2,4 with
mass1/10 and constant amplitude sqrt(1/20). Their Fourier laws differ, so
mixing the new component weights with the OLD acceptance rule is a meaningful
must-fail control. Predict B^2=8/5 versus the old mean2, full-law normalization
and FFT agreement. Also retain an exact-zero component and an equal-mass
case. Splitting one progression into disjoint pieces must preserve the TARGET
law but generally changes the proposal and increases B; don't confuse a
representation-dependent envelope with a property of the circuit alone.

Then reuse the unchanged physical r60 fixture and both covers, charging row
construction, both phase oracles, progression proposals and every coherent
ratio query. Whether reduced retries lower TOTAL named arithmetic work is an
empirical branch: weights add divisions/square roots and cover construction
still matters. No uncharged complete law as a sampler oracle, no broad native
timing sweep on the unresolved host, no arbitrary multi-defect extension.
Budget and freeze a bounded actual-RNG decision test before changing defaults.

### Implementation gate for that pilot

After the formula/progression-primitive pilot passes, add explicit
`proposal="root_mass"` to sample/forced_joint while retaining `"mass"` as
the default. Use raw omega and R=sum sqrt(omega), avoiding division by Z to
form tiny lambda values. Prepare component weights once per selected row,
then reuse the SAME conditional rejection loop. Extracting that loop into a
private row method permits an isolated synthetic-row RNG test; it must not
be mislabeled a physically consistent full-joint circuit fixture.

Freeze separate counters for component-weight preparation, square roots and
weighted-ratio divisions, in addition to the existing local/phase/Fourier/
progression/marginal counters. New weights add m square roots per row and m
divisions per coherent ratio. Common row stride means the interval marginal
count per proposal is the same for all its components. Return the mathematical
row-local expected attempts alongside actual attempts, never conflate them.
Empty rows, exact zeros, invalid modes and finite exhaustion remain explicit;
nonzero underflow must raise rather than silently drop a component.

Test the isolated Q8 conditional row by enumerating actual component selection,
interval bits, gcd lift and accept/reject branches with independent path
weights. Test the full physical joint law separately on r60 using all v0..7,
all work rows for auto and selected outputs for explicit cycle/dual. Use32
actual draws per v per proposal mode, cap512, and report expected arithmetic
categories plus observed counts; neither a seedwise win nor a timing win is
a frozen prediction. Preserve default-mode RNG paths against the prior report.
