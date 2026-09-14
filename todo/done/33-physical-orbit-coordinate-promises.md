---
id: 33
state: done
title: "Which indexed reflection gates have an inexpensive physical-orbit realization?"
outcome: "C74/C75 prove coordinate promises and clean physical constructions; tiny complex-gate and complete-output tests pass against the small static baseline"
claims: [C52, C56, C58, C59, C73, C74, C75]
---
# Audit the input promise before further synthetic scaling

**Completed at bounded first-pass scope, 2026-09-11: C74/PC and C75/CG.**
Item 4 now has a constructive clean-unitary proof, independently audited
tiny physical gate compilation and complete-output comparison. The static
baseline remains small; this is not a simulation speedup. A scalable
arithmetic compiler was not implemented: the test compiler explicitly
enumerates tiny truth tables. Do not scale it as a performance experiment.
TODO 35 owns the next cheap-physical-phase question; TODO 34 owns reliability.
The original scope and predictions below are retained as history.

**Algebraic stage completed: C74/PC.** The independent integer/rational tests
and claim regression pass. Main pursued this cheap algebraic direction while
TODO 32's long timing was paused for the runtime issue now owned by TODO 34.
Items 1–3 and the constant-error reduction below now have direct proofs;
the original predictions are retained as history, not open conjectures.
At that checkpoint the next step was item 4, now completed above. Its
distinction between basis-action verification and gate synthesis remains
important: C75 separates the existence proof from the finite compiler.

No new propagator, unbounded orbit scan or manuscript expansion.
Use lower-cost initial testers, with main auditing the proofs and predicates.
Keep actual modular encodings x=a^j mod N distinct from an explicit j register.
The order r=bM, generator, orbit membership and small b are supplied here.

## Original candidate predictions (before measuring)

1. The fine coordinate p=j mod b can be recovered from x^M using a b-entry
   table of powers of a^M. Thus cell reflection j=bm+p -> -bm+p can be
   evaluated as x -> a^(2p)/x mod N without learning m. Check on prime and
   composite moduli, non-coprime b,M, endpoints and every tiny orbit label.
   Plain inversion is a must-fail replacement when p is not zero.
2. A CLASSICAL evaluator returning D_q's phase in turns to circular error
   strictly below 1/(2s), s=M/gcd(q,M), reveals m modulo s by rounding and
   inverting q/gcd(q,M). With p this yields j modulo bs. If s=M it recovers
   the full discrete logarithm. Test exact rational error intervals, including
   wraparound and a deliberately insufficient-precision counterexample.
   This is an oracle reduction, not a new discrete-log algorithm: a reference
   phase generated from a known j must never be described as computing it
   from x. It does NOT prove that every quantum gate synthesis requires a
   classical phase evaluator or that quantum implementation is impossible.
   A stronger candidate is a UNIFORM constant-error classical evaluator:
   query D_q on (x*a^(-p))^(2^i), then unwrap inverse doublings from highest
   i to lowest. For circular error at most 1/16, nearest-half selection
   should contract error at each step, allowing root rounding after
   ceil(log2 M) steps when q is coprime to M. Audit wraparound, uniqueness
   and adversarial rational errors before claiming this reduction. It is not
   constant-cost access to high powers of a black-box quantum gate.
3. If s is small, extracting j mod bs instead uses a bs-entry subgroup table
   and x^(r/(bs)). This evaluates D_q without a full orbit table. Compare with
   C59's exact static regrouping: several small-order route labels may already
   give a small static group via their least common multiple. Fixed alphabet
   size alone is not evidence that the finite-time method beats this baseline.
4. Check whether repeated fine-work gates can be realized by extracting p,
   removing a^p, applying W on a small p register and restoring the modular
   encoding with clean ancillas. A basis-action sketch is not a compiled
   circuit proof; leave synthesis open unless independently validated.

Freeze tiny fixture sizes, setup/query counters and allocations before tests.
No runtime-hardness claim from these rows. An efficient special physical family
would be interesting; discovering hidden coordinate work is also a useful
negative result and should redirect the next experiment.

## Prior-art lead already inspected

Main read Theorem 6's statement and Appendix B's setup, especially B.1.1 and
the rational phase-evaluation premise in B.2 of
[Bermejo-Vega, Lin and Van den Nest](https://arxiv.org/pdf/1409.4800).
It explicitly distinguishes forward encoding from inverse-coordinate
recovery and relates the latter to generalized discrete logarithms. This is
positioning for the question, not a claim that our arbitrary W or coherent
reflection rotations belong to their normalizer gate class, nor an audit of
every theorem dependency in that paper.
