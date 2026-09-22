---
id: TX51
field: "Low-space generic-group order algorithms and spectral sampling"
status: open
effect: unknown
one_line: "Two bars: for the factor output, ECM and Pollard rho already factor with polynomial memory (heuristic); for C52's exponent law, a fully charged low-space exact-order routine (audit deprioritized)"
source: "Sutherland Order Computations in Generic Groups, sections 2.3 and 3.1 (Algorithms 2.1, 3.1; Propositions 2.5, 3.1) and 5.2; Lenstra 1987 (ECM); Pollard 1975, BIT 15; Brent 1980, BIT 20 (ECM/rho cited as standard, not body-audited); C52"
claims: [C52, C124]
notes: [FG, GR]
todo: [74]
---
# TX51 — Move the memory bill into explicit order discovery

## Dictionary

C52 already obtains the exact ideal exponent-output law from a uniformly
sampled latent eigenphase once r=ord_N(a) is known. Its simple traversal already
uses small space at O(r) modular-multiplication time. The new proposed transfer
replaces that deliberately weak discovery baseline with a known low-space
order routine, charging its whole pipeline before claiming a better tradeoff.
No quantum state vector is required for the conditional extraction step.

## Hypotheses

[Sutherland's thesis](https://math.mit.edu/~drew/sutherland-phd.pdf),
section 3.1, Algorithm 3.1 and Proposition 3.1, computes a collision exponent
E, factors it, and extracts the exact order. Its square-root walk analysis
assumes a random-function approximation. Section 5.2 offers a low-space
alternative to stored searches. Do not transfer median-order distributions
to RSA instances or count a collision multiple as the exact order.

Algorithm 3.1 uses a distinguished-point table; specify its bounded storage
or an actual alternative detector rather than calling the printed routine
constant-space. Include walk setup, collision detection, exponent sizes, restarts, E
factorization with its own workspace, certified prime-divisor reduction,
and scalar sampler precision. Target total space is discovery workspace plus
polynomial scalar-sampling workspace and retained outputs. Expected walk
behavior is not a worst-case N-only runtime theorem. C52's scalar operation
count is not already a numerical bit-cost certificate.

The contract is ideal clean-arithmetic exponent outputs, not full-state
fidelity or joint work measurements. A joint TV-epsilon implementation must
bound all phase/conditional probability errors and random-bit draws. Compare
complete order algorithms and factoring-first reconstruction at the same
output accuracy; compare GNFS separately when the requested answer is a factor.

## Consequence for the goal

This may improve an honest memory/time baseline, while discovery time remains
exponential in input bit length. It is not a new Shor speedup and does not
remove TX15's barrier. The deciding audit is the complete discovery and
extraction bill, especially E factorization. Since 2026-09-22 this is bar (ii)
below, and TODO74 deprioritizes it: it matters only for C52's exponent law.

## Split into two bars (2026-09-22)

**(i) Factor-output bar.** ECM (Lenstra 1987; heuristic L_p[1/2, sqrt 2],
polynomial memory) and Pollard rho (heuristic expected O(p^(1/2)) <= O(N^(1/4))
with Floyd/Brent detection) already factor with polynomial memory. A rho-style
walk on powers of a is kept only as a construction. Its collision gives a
multiple E of r (Sutherland 3.1, with the restart strategy of Proposition 3.1(1)
bounding E). Algorithm 2.1, given the multiple E, runs its squaring chain from
a^m with m the odd part of E, and succeeds exactly when Shor's
condition holds (C124), so E never needs factoring. For RSA N and typical a its
sqrt(r), about N^(1/2), is weaker than both methods, so it is not the bar. For the factor
output, memory is not the exponential that matters.

**(ii) C52 exponent-law bar.** The exact r, and so this row's original audit:
E factorization, certified reduction and the scalar sampler's precision. The
audit is deprioritized, not dropped (TODO74). It matters only for the
exponent-output law.
