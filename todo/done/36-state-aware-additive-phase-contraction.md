---
id: 36
state: done
title: "Can a state-aware rewrite bypass dense additive-phase sector coupling?"
outcome: "C53/UT close early-register transfer; C77/UP prove and validate the uniform-prefix support certificate; conditional suffix work continues in TODO 37"
claims: [C53, C54, C56, C57, C76, C77]
---
# Test the strongest small rewrite before scaling the dense representation

**Completed at bounded scope, 2026-09-11.** C53/UT own the first diagnostic;
C77/UP own the uniform-prefix theorem, exact integer helper, independent
complete-law checks and retained failures. The prefix does not factor from
the suffix. TODO 37 owns the stronger Fourier-feedback-aware comparison.
The plans below are historical, not instructions to rerun completed probes.

C76/GS complete TODO 35's physical kernels and actual-output question.
Dense prime-sector coupling closes a specific static grouping route, not
simulation of the chosen input and output. Omitting the kick remains a good
low-accuracy baseline on the tiny rows. Keep comparisons at a stated accuracy.

## First bounded test: general exponent-unitary transfer

**Completed.** C53 owns the general criterion and output-parity obstruction;
UT owns the failed normalized-control attempt, corrected physical laws and
independent formula audit. The state-optimal candidate loses to omission
on the measured output, and the parity witness excludes exact early-register-
only replacement in the fixed numerical cases. The plan below is preserved
as history; do not rerun or optimize the state-fidelity objective by default.

The C53 test in GS only rules out a diagonal scalar phase per early control
history. A general early-exponent unitary could still transfer the work gate.
Freeze the SAME N13,a2,b3,r12,t3 background, insertion after s=2 controls,
and k=0,1,2. Do not tune the circuit to increase the effect.

Let A be the r-by-L matrix of pre-kick conditional work columns divided by
sqrt(L), L=2^s, and B=G_k A. An early-exponent E acts as A E^T. With
C=A^dag B=U Sigma V^dag, choose E=(U V^dag)^T. Before measuring, derive

    F_root = ||C||_1,
    min_E ||B-A E^T||_F^2 = 2-2 F_root,
    min_E pure_joint_trace_distance = sqrt(1-F_root^2).

Exact transfer is possible iff G_k rho G_k^dag=rho, rho=A A^dag.
The best diagonal-phase overlap is sum_l |C_ll|, while omission has overlap
|Tr C|. Neither output TV nor gate complexity is optimized by this SVD.
This is standard purification fidelity/unitary optimization, not a new
quantum information theorem. Read the relevant primary statement before
crediting it. E is placed at the insertion; it cannot generally be moved
before the earlier controlled arithmetic.

Use the existing branch matrices and a tiny matrix diagnostic, not another
propagator. k0 is the exact null. A normalized embedded four-column Fourier
purification has diagonal work rho and must permit nontrivial exact general
transfer while defeating diagonal-only transfer: this is a deliberately
synthetic algebra control, not the physical initial state. Classify the
actual k1/k2 result without assuming exact transfer or its failure.

If feasible, embed target and candidate states at the insertion into the
existing physical Circuit/statevec suffix, then compare COMPLETE output
laws. Include omission and diagonal-only candidates. The SVD optimum is
for the joint state; a worse output law than omission is allowed and must
be retained. Compare the target with the already independent full-r law.
Charge the physical suffix and finite reference arrays explicitly.

Lower-cost agents do initial tests and proof audits; main verifies predicates,
normalization, resource guards and retained failures. Core first (already
passed for the current C76 science turn), TEMPLATE/harness, 16 MiB aggregate
numeric cap and no long timing. No production API, manuscript or default
changes are needed for this diagnostic. A length-2^s matrix is not a scalable
exponent-only sampler merely because it is small for s=2.

## After that diagnostic

**Next discriminating question: exact uniform-prefix certificates from
physical support, not a sector-sparse state.** Derive this candidate first.
For t controls, let d be the requested LOW-output prefix length, H=2^d and
L=2^(t-d). Trace the unmeasured low input controls l. Conditional on l and
high input history h, C57's co-moving argument gives final work support
within l+L*h+[-R,R] modulo r, where R is the sum of all repeated-mixer
displacement radii. Orbit-diagonal phases do not enlarge this support.
If

    min_(1<=q<H) circular_distance(L*q mod r, 0) > 2R,

then different high histories have disjoint work support. Their reduced
high-input state is maximally mixed, suggesting that the first d inverse-QFT
output bits are exactly uniform even with dense sector coupling. Verify the
partial-QFT/no-signalling identity and normalizations explicitly. This is a
sufficient criterion; a failed separation test does not imply nonuniformity.
It does not by itself finish conditional sampling after the certified prefix.

Freeze a genuinely nontrivial bounded pilot: supplied (N,a,r,b)=(61,2,60,3),
t=6; verify the order before any claim. Use the existing three-state W with
initial angle pi/4, one W after control i=4 with angle -pi/10, and one after
i=5 with angle pi/11; omit other W gates. Put G_1 after the i=4 W. This is
an explicitly NEW sparse-mixer fixture, not a retuned row in GS/UT. Here the
conservative radius is R=6. Vary only requested prefix length d=1,2,3 on
this SAME circuit. The d1 high shift is 32 and has circular distance 28,
so the proposed certificate predicts an exactly uniform first bit despite
the nontrivial mixers and phase. d2 reaches the boundary distance 12 and
is not certified; retain either output outcome. No dense-sector truncation
or high-order phase oracle is supplied.

Use the existing full-r sequential_path within its dimension cap for the
independent small reference; no new propagator or large physical truth-table
compiler. Preflight reference temporaries and full-output enumeration if
needed, with 16 MiB numeric cap and explicit call/matrix-work budgets.
Prefer prefix marginals without complete enumeration if the EXISTING API
supports them correctly. A must-fail control sets the last controlled work
shift to identity: the high bit then stays |+>, yielding deterministic even
parity instead of uniformity. This intentionally altered control violates
separation; verify common terminal-W cancellation. Add k0 as a phase null
only if it answers a separate predicate, not as a replacement for nonzero-k1.
Lower-cost agents should initially check the certificate and reference,
main audits the partial-output target and actual allocation/counters.

If the certificate survives, ask whether its modular near-collision test
can be evaluated without enumerating H, and how much it saves in a SAME-
accuracy sampler. Do not promote a uniform-prefix shortcut into a full
efficient sampler; exact intervals after conditioning are a different task.

Seek a formula for the actual output/low-bit marginals, not another dense
support count. C57's short orbit-displacement cone survives a diagonal G_k
for each FIXED control history, but it does not remove the coherent sum over
histories. Moving G_k through controlled multipliers makes its additive
frequency history-dependent. Any proposed additive-character/finite-offset
contraction must explicitly evaluate or bound that remaining sum, charging
modular arithmetic, order information and precision.

Reject a synthetic scaling test in which N grows with a fixed tiny prefix
and bounded physical residues: G_1 can then tend to identity for a trivial
reason. Derive a discriminating fixed-input or controlled-parameter test
first. Known Gauss norms, a low-rank matrix or easy pointwise phase evaluation
do not supply coherent Gauss phases or a sampling algorithm for free.

TODO 34's host/runtime question and TODO 32's unfinished timing sweep remain
separate; do not launch stress tests or infer hardware repair from small passes.
