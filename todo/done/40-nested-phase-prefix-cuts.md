---
id: 40
state: done
title: "Multiple diagonal insertions: does a binary cut avoid multiplying phase periods?"
outcome: "C81 proves nested-prefix covers and exact finite-truncated construction counts; bounded amplitudes, edge controls and longer-row construction savings independently reproduced"
claims: [C78, C79, C80, C81]
---
# Use nested prefixes before enumerating phase histories

**Complete at mathematical / bounded amplitude scope:** C81 owns the proofs
and selector; NC owns authoritative integer, amplitude, edge and selected
long-row reports, controls, audit corrections and validation. The multiple-
phase row constructor is experiment-local, not a production sampler. TODO41
owns integration and end-to-end work accounting. The plans and unimplemented-
state wording below are retained as their pre-measurement historical record;
do not rerun these completed discriminators as new research.

This follows the bounded C80 component-weighting checkpoint. The potentially
more structural question is whether several earlier diagonal phases require
independent period factors. Main derived the following candidate; an
independent lower-cost read-only audit confirmed the algebra. It is not yet
a multi-phase sampler, measured construction improvement or novelty claim.

## Candidate construction, frozen before testing

Retain C79's known/indexed r=b*M and single late mixer at split s. Insert d
pointwise unit-modulus phases f_i after v_i ascending low controls, v_i<=s,
with NO intervening mixers. For fixed work j, high history h and local source
x=c+q, choose a cut w<=s and write l=a+2^w*k. Here 0<=a<2^w,
0<=k<K=2^(s-w), and u+l=x mod r on a compatible local path.

For v_i>=w the phase argument becomes

    x - 2^v_i * floor(k / 2^(v_i-w)) mod r.

It is fixed within (h,k,q), even for phases exactly AT the cut. For v_i<w,
the argument is u+(a mod 2^v_i) mod r. Along a=rho+n*r, all left-prefix
residues repeat with ONE period

    P_left = 2^V/gcd(r,2^V), V=max{v_i:v_i<w},

or P_left=1 if there is no left phase. Nested powers of two make this the
maximum period, not their product. For rho=(x-2^w*k-u) mod r, refine n mod
P_left; each nonempty class is a progression with common stride r*P_left,
start L*h+2^w*k+rho+r*z, and the usual finite count below 2^w. Multiply ALL
phase factors into each local pair BEFORE coherent merging. Old canceled
candidates cannot be dropped. Different h,k,rho,z have disjoint supports.

The sufficient local construction scale is O(H*b^2*K*P_left), with truncation
to actual nonempty classes, plus the d-dependent phase-query/product work.
Candidate cuts at phase positions and s suffice for the untruncated factor:
between insertions P_left is constant and K decreases. This suggests a link
to the largest phase-free interval in the ascending binary controls, not an
independent exponential penalty per phase. Derive the precise bound before
using that phrasing as a theorem. The result can still be exponential in
width; phase-oracle cost and known-order access remain assumptions.

Main's sharper derived prediction, frozen before its follow-up integer run: let
alpha=nu_2(r), B=max(0,s-alpha), and sort the distinct effective insertion
positions max(0,v_i-alpha), including endpoints0,B. If Delta is their largest
adjacent gap (zero when B=0), predict

    min_w K*P_left = 2^(B-Delta).

For w>=alpha the exponent s-w+max(0,V_left-alpha) is exactly B minus the
effective gap ending at that cut. Thus d insertions give the sufficient
worst-case upper factor 2^floor(d*B/(d+1)); d=1 recovers C79's square-root
factor. Endpoints and repeated positions add no effective gaps. This is an
untruncated structural bound, not a runtime or realized component count.
An independent read-only proof audit confirmed this identity, including the
at-cut/endpoint and duplicate-position cases. Its finite integer regression
is a separately frozen follow-up to the initial partition pilot.

## First discriminator: exact integer pilot only

Freeze r60,s7,u=0,1,2, all source labels x=0..59, and all cuts w=0..7.
Keep one phase at v=3; move only the second through v=3,4,5,6,7. Enumerate
compatible l in [0,128) and compare the original phase-argument word with
the cut formula, then the refined progression partition and constant word
within each class. Include equal insertion positions and endpoint v=s.
These are exact residue/partition tests, NOT amplitude or quantum-law tests.

Predict complete coverage without collisions and exact word agreement. For
positions (3,6), cut6 predicts K*P_left=4, versus endpoint factors16 and16.
The must-fail control drops the left refinement: u=x=k=rho=0, w6, the
compatible a=0,60 give different v3 arguments (0 and4), although both right
v6 arguments agree. Require this explicit witness; do not use mere unequal
upper bounds as the control.

Budget at most 5,000,000 actual integer relation/argument/partition visits,
charged before loops, and 1 MiB aggregate numeric payload. Stream one
(positions,cut,u,x) case at a time with at most128 labels; do not retain an
orbit/Q-by-r table. Preserve failed reports. Main owns interpretation and
integration; the lower-cost initial tester may implement this bounded pilot
from TEMPLATE, using lab.harness, with no production helper edits.

## If the integer mechanism survives

Next freeze a literal multi-phase amplitude test against existing local
block/shift/diagonal primitives, before generalizing the sampler. All-left
cycle and earliest-cut dual are exact competing covers of the SAME state;
full-r sequential propagation remains a width-linear baseline at fixed r.
Short progressions may hide the benefit in the r60,s7 fixture. The derived
s12, positions5,11 example predicts factors16 versus512 and128, but has NOT
been measured: use bounded selected rows, not a full dense-law size record.
Numerical certification, intervening mixers, native timings and manuscript
promotion remain separate. Stop or redirect if this only redescribes the
existing finite-support cap without reducing charged construction work.

## Frozen coherent-amplitude discriminator (2026-09-11)

The previous goal turn made concrete progress: C80's weighted sampler and
the NC integer mechanism were independently checked. This next stage tests
the unmeasured coherent coefficients, not a repeat of that partition audit.

Keep r60,b3,t8,s7,H2,N61,base2 and the ER work blocks/late G1 phase fixed.
Insert TWO identical pointwise G1 phases, at (3,v) with only v varying through
3,4,5,6,7. For each v, construct all work rows at each distinct cut in
{3,v,7}. The existing ER literal direct-column engine will accept a sequence
of insertion positions (same phase oracle) while keeping its one-insertion
default unchanged. Do not write a new vector propagator.

An experiment-local cut_row constructor groups local pairs for every h,k,
then refines each nonempty rho class by n mod P_left. Use the derived cut
arguments and multiply both phases BEFORE coherent summation. Its supplied
input contract is (r,b,t,s,W0,W1,g,[(v_i,f_i)],cut); no production sampler
integration at this stage. Return common-stride (start,count,coefficient)
components and the positive disjoint norm. Omit exact-zero coefficients only.

Predict coordinatewise agreement with every literal column, no overlapping
components, correct row/work norms, and selected geometric Fourier queries
y=0,1,Q/2,Q-1 agreeing with independent FFTs. Check both geometric amplitudes
and NORMALIZED conditional laws; normalization alone is not a discriminator.
The existing C79 helper is an exact cheap baseline when both insertions share
v=3 (combine their pointwise phases) and when v=s (absorb the source phase
into the relevant local coefficient, not the unchanged helper blindly).
Use all-left cut7 and earliest dual cut3 as exact competing covers for every
row. Whether the hybrid cut saves ACTUAL construction terms in this short
fixture is an empirical branch; the finite-support cap may erase the gain.

Must-fail control: freeze the earlier coefficient at z=0 while retaining the
correct refined supports for (3,6), cut6. Predict a nonzero amplitude error
somewhere against the literal columns; require a named work/exponent witness.
An independent r6,b2 Hadamard edge check will verify revival of an originally
canceled rho=0 candidate with f1(j)=(-1)^j at v1 and f2(j)=i^(j mod2) at v=s.
Its exact equivalent C79 reference absorbs diag(1,i) into W1. A wrong builder
that drops the original canceled candidate must fail on that same physical
row. Also test duplicate phases, initial/end insertions, empty phase lists,
and r=b aliases through exact reductions to existing literal/C79 machinery.

Freeze preflight <=16 MiB aggregate numeric payload for the full r60 pilot,
including simultaneous reference/FFT/magnitude arrays, retained summaries
and conversion copies. Stream one reference matrix per v and one expanded
row per cut/work; never supply the reference as a sampler oracle. Charge
actual setup, group-pair visits, coefficient-pair visits, phase calls/products,
expansion labels and Fourier component terms separately, BEFORE loops/calls.
Reserve each complete workload before execution; cap each named scalar
category at 5,000,000 and total reserved categories at 15,000,000. Counts
sharing work must not be presented as additive native FLOPs. The small
independent edge test gets <=4 MiB and <=1,000,000 named scalar visits.

Only after tiny/edge agreement, freeze a selected-row long-progression check
at s12,t13, positions(5,11), same r60 inputs and work labels0,1,59. Compare
cuts5,11,12, deriving component/local/phase budgets first. Reference only the
support-cone-compatible columns through the SAME literal engine, with exact
omitted-zero justification; do not allocate a full Q-by-r reference. The
untruncated factors are the earlier derived prediction, not measured timing.
Full-r sequential methods remain relevant fixed-r width-linear baselines.

### Long-row discriminator: frozen detailed budget

After the tiny amplitude and edge gates pass their audits, keep r60,b3 and
the same physical phase/work blocks; fix s12,t13,H2 and phase positions(5,11).
Vary only the representation cut through 5,11,12 for work labels0,1,59.
No new phase family or circuit parameter sweep is included. Compare every
amplitude and complete normalized conditional FFT law, plus geometric queries
y=0,1,Q/2,Q-1 with Q8192 explicitly supplied. The prederived K*P factors are
128,16,512. Whether the middle cut saves actual local/phase work is an empirical
branch; report losing as well as winning rows. Retain the wrong-left-z0
coefficient control on cut11 with a named amplitude witness.

For cuts5,11,12 respectively, per-row bounds are grouping2304,36,18;
coefficient pairs2304,288,1242; components690,160,690; phase calls4610,578,2486;
phase products7602,1024,4416. All are conservative geometric bounds, not
observed costs. Include an additional cut11 control row, expansion/comparison/
normalization visits, setup/unitarity work, and FFT input entries separately.
Reserve <=5,000,000 named units per category and <=15,000,000 total overlapping
ledger units, and <=4MiB aggregate numeric payload before execution. This is
not an additive FLOP count or RSS measurement.

An independent support-cone argument makes sparse reference queries complete:
initial label u and the late mixer input/output q,p lie in0..2, so for the
FULL exponent e=l+Lh, j-e=u+p-q lies in[-2,4] modulo60. Thus query only the
seven residues e=(j-delta) mod60, delta=-2..4, using the existing ER
direct_column(..., split=12, early_insertions=(5,11)); all omitted columns
are algebraic zeros for that work row. Stream one 60-entry column at a time,
at most959 calls per work and2877 total (not Q-by-r storage). Reserve at most
517860 reference block products and517860 reference phase queries (the latter
also bounds modular-power queries), plus four shifts per call. Reuse existing
cut_row; do not introduce a second vector propagator. A full-r sequential
method still provides a fixed-r width-linear baseline, so any measured gain
is only a comparison between these exact row constructions.

### Follow-up conjecture: exact finite-truncated selector (frozen before test)

Main noticed a simplification after the independent long-row geometry counts.
For A=2^w and canonical rho, N_rho is either floor(A/r) or ceil(A/r).
The integer P therefore either clips ALL candidate lengths to P, or clips
NONE. Let T be the total compatible local-pair count over the full low
interval, summed over h,u,q, and S the analogous count over distinct geometric
rho classes per h. Both can be computed with at most H*b^2 local terms,
without enumerating K or P. Predict for any cut:

    coefficient_pair_visits = min(H*K*b^2*P, T),
    geometric_candidates_before_cancellation = min(H*K*R*P, S),
    group_pair_visits = H*K*b^2,  R=min(r,2*b-1).

In the unclipped case the k-slices partition the original compatible low
exponents, so their total is T or S; in the clipped case all pairs/residues
contribute P per slice. This candidate identity would give an exact selector
for THIS naive direct construction's named counts, including phase calls
H+d*C and phase products(d+1)*C+G, without evaluating a phase oracle. It
does not predict cancellations, caching, Fourier/rejection cost or runtime.
Between phase positions P is fixed and K decreases, so phase positions and
s still suffice for minimizing these nonnegative named construction costs.

Before promoting, require an independent proof audit and <=100,000 exact
integer pair/class visits including non-divisor/aliased r=b, empty and
duplicate/endpoint phases, short slices and both clipping branches. Compare
direct h,k,rho,z enumeration against independently computed T/S minima.
Retain the C81 short-row largest-gap cost inversion as the must-fail control
against using K*P alone. No amplitude or phase-oracle queries in this test.
