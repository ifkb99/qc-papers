---
id: 38
state: done
title: "Can measuring work first replace dense sectors by short compatible-history rows?"
outcome: "C78/WF implement opt-in work-first progression draws; independent tiny laws, actual decision paths, edge guards and same-output float comparisons pass; numerical certification is not implied"
claims: [C54, C55, C56, C57, C76, C77, C78]
---
# Change the latent measurement, not the physical circuit

**Complete at bounded scope:** C78/WF own the proof and opt-in sampler,
independent complete-law comparisons, actual RNG transition/retry checks and
bounded draws. The stronger FB approximation still misses the frozen tiny-law
accuracy; the new helper agrees with the exact sequential reference. This is
not a finite-bit accuracy certificate or timing superiority result. Everything
below is the retained pre-measurement plan. TODO 39 owns the next earlier-phase
schedule boundary; do not rerun these completed formulas as that investigation.

This arose while TODO 37's stronger feedback-aware baseline was being tested.
Measuring work in an indexed computational basis commutes with the terminating
exponent QFT. Unlike final-eigenphase conditioning in C54, a finite-displacement
promise makes the compatible exponent histories short when Q/r is bounded.
Dense sector coupling need not imply a large conditional row in this basis.
This is a promising restricted sampler route, not a claimed new general theorem
or factoring breakthrough. Check against the actual same-accuracy baselines.

## Derivation before the first pilot

Write a_(e,j)=<j|phi_e>, with each work column phi_e normalized and all Q=2^t
exponent histories initially uniform. Suppose phi_e has support in
e+[-R,R] mod r, in a supplied known/indexed orbit. A common final work unitary
may be stripped before defining this basis and R.

1. Draw an auxiliary exponent e uniformly, evaluate its coherent work column,
   and sample j from |a_(e,j)|^2. This gives the exact work marginal
   p(j)=Z_j/Q, where Z_j=sum_e |a_(e,j)|^2. The sampled e is only a classical
   auxiliary variable; DISCARD it rather than using a measured-exponent state.
2. Enumerate all allowed e for this j via e=(j-delta mod r)+n*r in [0,Q),
   -R<=delta<=R, deduplicating residues if needed. Evaluate all their complex
   amplitudes a_(e,j), not just their squared magnitudes. The normalized
   conditional exponent state has these amplitudes divided by sqrt(Z_j).
3. Fourier-sample that sparse state. A uniform-output rejection envelope has
   mean S_j proposals, where S_j is the number of candidate labels (zeros may
   be retained). A stronger baseline uses partial-Fourier orthogonality:

       P_d(z|j) = 1/(2^d Z_j) sum_(c mod L)
          |sum_(e=c mod L) a_(e,j) exp(-2*pi*i*z*floor(e/L)/2^d)|^2,
       L=2^(t-d).

   Each prefix groups only the supplied sparse labels. Two child masses
   give the next output bit, without an output array or rejection loop.
   Charge integer grouping/sorting and phase evaluation, not just amplitudes.

For 2R+1<=r, S_j <= min(Q,(2R+1)*ceil(Q/r)); a wider cone can use min(r,2R+1)
distinct residues but generally loses the useful small-row bound. Thus this
does NOT remove the Q/r factor at high precision. The usual large-Q regime
and supplied order/index discovery remain substantive limitations.

The support promise alone does not supply cheap amplitudes. State the cost
T_col of evaluating a coherent fixed-history column explicitly. Finite-
displacement matrices plus pointwise computable diagonal phases can be applied
within a small intermediate cone, but that implementation and cost must be
validated separately; sampling a classical path through mixers is invalid.
Arithmetic-operation claims do not certify float rare branches or bit cost.

## Frozen initial test and controls

Use the SAME N61,a2,r60,b3,t6,k1 fixture from UP/FB; strip only its common
terminal W so R=4. Lower-cost initial pilot owns
experiment_work_first_rows.py and uses the existing closed tiny dense boundary
as a CHARGED oracle, not a scalable implementation. Main audits independently.
TEMPLATE/harness, predictions before measurement, core already passed at the
C77 science checkpoint, 16 MiB aggregate numeric cap and before-loop work caps.

Enumerate the tiny complete joint law only for verification. Check the work
weights from uniform-history sampling, modular row enumeration, all prefix
normalizations and the complete 64-point exponent law against independent UP
QFT/sequential references. Count queried columns, candidate entries and actual
prefix/Fourier terms; no uncharged orbit/output tables. A must-fail control
dephases each conditional exponent row: every such diagonal state gives
uniform Fourier output, which the frozen full target does not have.
Do not assume uniform work weights; classify any such extra control honestly.

If this pilot passes, the next separate step is a no-orbit-array fixed-column
implementation for the literal sparse schedule, bounded samples and a charged
comparison to omission/feedback baselines at the same requested accuracy.
Do not rewrite existing propagators merely to check them. General gate APIs,
large timing sweeps and numerical certification remain deferred until useful.

Relevant primary context already read by main: Van den Nest, arXiv:0911.1624,
Definition 1 and Section 4.3/Theorem 3 (CT/sample-query and sparse operators);
Schwarz and Van den Nest, arXiv:1310.6749, Theorem 1 (sparse OUTPUT promise).
Neither citation alone establishes priority for this conditional-row
specialization. The proposed Fourier prefix step is elementary orthogonality,
not a new QFT identity. TODO 37 still owns its earlier bounded comparison.

## Stronger follow-up: few late histories, long progression components

**Formula checks complete.** WF owns the primary and nontrivial width-eight
checks, corrected pointwise phase access, main reproduction and retained
failures. The next action is the opt-in sampler described below. The intervening
diagnostic plan is historical; do not rerun it in place of implementing a draw.

C78 derives the one-middle-mixer/pointwise-diagonal-phase formula. Keep the
literal schedule: W_init, the first t-1 ascending controls, W_pre then G_1,
the final high control, and a common terminal W that can be stripped. Its
conditional work rows have a bounded number of progression components,
independent of Q/r. The old Q/r-sized row is not the strongest comparator for
THIS more restricted schedule. Multiple extra noncommuting insertions remain
outside the proof.

The lower-cost initial progression pilot is experiment_work_first_progressions.
Main found two limits: its phase/label queries read the supplied orbit table,
and the UP width-six fixture has only singleton progression terms. Correct
the first using pointwise modular powers with charged queries. Add a separately
frozen width-eight check at the SAME N,a,r,b,angles,k with the same
penultimate-control-relative insertion. It is a specified one-parameter family,
not a silent retuning of UP. Now L exceeds r and actual multi-term progressions
must occur. Compare expanded rows against literal sparse columns, and geometric
Fourier sums against direct finite-root sums over those rows. This secondary
check is not an independently compiled physical-circuit reference.

Before its loops, cap numeric storage at 16 MiB and preflight actual column,
row, expansion, component and direct-root term counts. No amplitude threshold
in an acceptance denominator; zero proposal gets zero acceptance. Retain
the original bounded reports and all failures. No random sampling or timing
is implied by verifying an enumerated accepted law.

After those formula tests, implement an opt-in sampler using existing
progression_sample, with coherent local-column sampling of work j, merged row
components, component-weight proposal and coherent acceptance. Require actual
bounded draws, zero-weight/zero-proposal tests, safe exhaustion and invalid-
input gates. Never keep the seed exponent after conditioning.
Keep the sampled work label fixed across output-proposal rejection retries;
redrawing it would bias its marginal by the row-dependent acceptance rate.
No orbit table, work/output vector or list of all progression entries may be
built in a draw.
Counter audit must include pointwise modular powers, phase precision and setup.

Compare to FB's strongest inexpensive candidate and G omission at stated
accuracy, and to the existing small exact reference where feasible. Do not
claim a large-size speedup from removing deliberately dense reference arrays.
Scalar exactness, float empirical accuracy and certified finite-bit sampling
are distinct outcomes. Wide/high-precision samples are a separate preflighted
test, not permission to resume TODO 32/34's crash-prone timing work.

## Frozen implementation test (2026-09-11)

Add opt-in `lab/work_first.py:LateWorkProgressions`, taking known period,
block size, exponent width/split, two supplied small unitary ndarrays and a
deterministic pointwise unit-modulus phase callable. Preflight H*b^2 local
terms, H*min(r,2*b-1) components and numeric storage before copying/looping.
Expose bounded column/row diagnostics, forced joint probabilities and actual
RNG sampling; reuse `progression_sample` unchanged. Count actual column/row
local terms, phase calls, Fourier components, progression proposals and scalar
prefix queries. No order discovery or global phase-oracle certification is
implied. Invalid inputs, numerical nonfinite/zero anomalies and exhaustion
raise; zero mathematical proposal has acceptance zero, without a cutoff.

Before measurement: predict all tiny forced joint laws agree with the earlier
independent reference, mean conditional attempts equal m, and actual bounded
draws exercise both work sampling and coherent rejection. Freeze primary
t=6 and relative-insertion t=8 fixtures unchanged; use modest fixed-seed
samples as smoke tests, NOT a 1% distribution certificate. A separately tiny
deterministic RNG traversal checks the implemented progression/acceptance
decisions. Controls must detect dropping coherence and redrawing the work
label during retries. Compare same full-output laws against FB and the exact
small-work sequential reference; any accuracy conclusion concerns enumerated
tiny float laws, not certified arbitrary-width sampling or timing superiority.
Lower-cost agents own independent initial tests/audits; main owns the helper,
regressions and final evidence review. All runs stay bounded and single-BLAS.
