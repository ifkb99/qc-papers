---
id: 73
state: open
title: "Factoring memory after the cache and recurrence audits"
outcome: "Cache lifetime pilot measured; exact recurrence audit complete; production applicability and constructive lemmas remain open"
claims: [C118, C119, C120, C121]
---
# 73: Factoring memory and aggregation

**Triage (2026-09-22, METHOD.md barrier check 4):** none of the three parts computes anything complexity-forbidden (they concern memory constants and certificates inside known factoring algorithms), so it passes. It concerns factoring software rather than these simulations: kept as the user's preserved direction, not scheduled.

Notes FD and FE own discovery/provenance; C120/C121 own the completed cache
pilot and exact symbolic audit; TX44-TX47 own transfer scope. TODO72 (the
original factor-state direction) closed at complexity triage on 2026-09-22. The selected user-approved measurements
are complete; the items below are future decisions, not active experiments.
The user accepts moderate slowdown for memory; no numerical tolerance was given.

## Remaining practical question

C120 establishes a construction-memory saving on the tested actual CADO path.
Before expanding the implementation, identify the dominant peak on a
representative filtered GNFS matrix: dispatch, cache preparation, final cache,
or another stage. Pin layout, vector width, allocator and actual backend path.
A later comparison must include complete read/build/save cost and account for
production dispatch and concurrent submatrices. Compare existing sequential
build and finer partitions when applicable. No percentage from synthetic data
is a production guarantee.

If temporary lifetime overlap dominates, validate the small existing B patch on
that workload with a new bounded reviewed design. If steady-state cache/solver
storage dominates, this patch does not address it; audit that allocation before
proposing compression or recomputation. The standalone reader change did not
demonstrate a material independent peak benefit here. Smooth-part scheduling
remains secondary, conditional on that stage being the relevant peak.

## Remaining algebraic question

C121 completes the selected ratio/seed audit and refutes the one-class/one-run
shortcut for a legitimate published-generator implementation. Do not repeat a
numerical recurrence sweep: the missing object is a constructive short cover
or small shift-boundary lemma for the actual required candidates.

A new derivation must preserve coverage, multiplicities, exact-match processing,
per-prime order conditions and witness recovery, while charging cover creation,
seeds and exceptional gcd recovery. Discovering a favorable cover only after
storing the old full list retains its construction peak. If a useful lemma
closes, compare its full time-space bound with structured multipoint evaluation,
blocking and streaming products. C121 is not a generic factoring obstruction.
TX45's later source leads remain unaudited for new memory guarantees.

## Speculative interval-count question

TX46 supplies an exact interface but no efficient constructor. Specify one
shared dissection rule for N-1 < pq <= N and derive exact clipped integer-region
counts, including boundaries. Decide whether common regions can be discarded
without already solving their divisor question; bound retained regions,
coordinate bit lengths and repeated adaptive queries.

If this merely repackages separate summatory calculations or visits every
denominator, stop that variant. If a cheaper certificate exists, prove its scope
before designing a test. A counterfamily is a scoped obstruction, not a generic
factoring bound. Squares, empty/singleton intervals and one-divisor intervals
are future correctness cases. Approximation would need absolute-error control
of the small count, not only relative error in large sums.

## Interleaving and completion

The candidate pipeline emits sound empty-region certificates and disjoint
unresolved intervals, then runs bounded-memory arithmetic search. Bound total
remaining work and charge all preprocessing/rejection. Product gcds are not
branch counts (TX47); geometric formulas do not automatically apply to
arithmetic progressions. Compare against factoring first and constructing the
state afterward. GNFS is the practical output-level baseline; published
collision algorithms are matched arithmetic baselines.

Use local reading/derivation first, no paid compute or challenge-scale search.
Later execution needs a specific reviewed design and explicit cap. Completion
may be a concrete memory tradeoff, a recurrence with hypotheses or a useful
obstruction. C120 is a scoped measured preparation saving; no generic RSA breakthrough
or full Shor-simulation improvement follows. C121 supplies no additional
measured saving to multiply into it.
