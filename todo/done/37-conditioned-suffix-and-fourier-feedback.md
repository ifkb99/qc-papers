---
id: 37
state: done
title: "Does Fourier feedback explain the conditional suffix, or is coherent history needed?"
outcome: "C77/FB validate the exact conditioned boundary; the feedback-aware approximation improves on omission but misses the frozen 1% TV target"
claims: [C49, C56, C57, C76, C77]
---
# Retain the known Fourier feedback before testing a coherence approximation

**Completed at bounded scope, 2026-09-11.** FB owns independent boundary and
sector/full-r checks, the retained failures and same-output comparison.
The stronger approximation is not exact and misses the requested accuracy.
The plan below is preserved as history. C78/WF own the work-first connection
that arose here; TODO 38 owns its implementation and comparison.

C77/UP establish exact uniformity of the first two output bits on the frozen
N=61,a=2,r=60,b=3,t=6 sparse-mixer fixture, but not their independence from
the suffix. A product-of-marginals comparison throws away even ordinary
inverse-QFT feedback. Test a stronger inexpensive comparator on the SAME
fixture before pursuing a general contraction algorithm.

## Derivation before measurement

Put H=4,L=16,Q=LH, e=l+Lh and y=z+Hw. Condition on z=y mod H. Because
p(z)=1/H, the normalized low-input/work boundary before the remaining F_L is

    |eta_z> = (LH)^(-1/2) sum_(l,h)
       exp(-2*pi*i*z*h/H) exp(-2*pi*i*z*l/Q) |l>|phi_(l+Lh)>.

The second phase is the known feedback; it must not be omitted. Check the
normalization, sign and bit convention against the existing full QFT.

For this sparse schedule, with h=h0+2h1,

    phi_(l+16h) = T_h U^l W_init |0>,
    T_h = W_final U^(32h1) G_1 W_pre U^(16h0).

Dephase the high INPUT histories h, dropping only h!=h' terms while keeping
the feedback phase on l. Each diagonal h component has a common work unitary
T_h, which cancels under the work trace. The candidate conditional law is
therefore the first four controls from W_init|0>, with product exponent
phases exp(-2*pi*i*z*2^i/Q). It is independent of the late G_1 and mixers.
This is an approximation to test, not a consequence of uniform prefix alone.

The existing periodic-sector representation has a supplied uniform coarse
sector and three-state initial block. Use its existing shift construction
and the existing sequential_path instrument with phased branch-one matrices.
No new propagator or production sampler/API is justified at this stage.
Charge supplied order, sector setup, enumeration of sectors for the tiny
complete reference, matrix products, QFT updates and all retained arrays.

## Frozen bounded comparison

Lower-cost agents perform initial independent tests; main audits formulas,
predicates, counters and failure retention. Use TEMPLATE/harness; core first
(the preceding C77 checkpoint passed). Before allocation, enforce 16 MiB
aggregate numerical storage and explicit finite call/work caps; this is not
a host stress test, timing comparison or RSS bound.

1. Reproduce the full target and exact conditioned-boundary law through the
   existing QFT engine. Record all four suffix laws, not just their average.
2. Compare the high-history-dephased candidate with that SAME complete law.
   Record TV and whether it meets fixed tolerances 1e-3 and 1e-2; either outcome
   is informative. Independently reproduce the candidate using small sectors.
3. Include omission of G_1 as a full-task baseline. Vary ONLY k=0,1 in this
   new experiment, not the frozen angles, modulus, width or insertion time.
   Do not silently expand the completed UP experiment's call budget.
4. A synthetic normalized boundary with one low-input Fourier basis state
   must expose omission of the feedback phase for a nonzero z. Derive its
   shifted finite-Fourier law first. This checks convention, not physical
   circuit difficulty. Additional physical no-feedback/sign controls may be
   reported without assuming they must be distinguishable.

If the candidate succeeds only at coarse accuracy, say so. If it fails,
inspect the surviving h!=h' overlap terms and their finite-support geometry;
do not respond with another static dense-sector count. A useful bound must
control the actual full-output error and its cost without enumerating all
histories. No efficient full sampler, asymptotic speedup or novelty is implied
by one small law or by known terminating-QFT feedback.

TODO 34 host/runtime reliability, TODO 32's unfinished timing sweep and
TODO 24's working-precision complexity remain separate and deferred here.
