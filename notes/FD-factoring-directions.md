---
code: FD
date: 2026-09-21
title: "Factoring directions: practical memory, arithmetic aggregation and interval counts"
outcome: record
claims: [C118, C119]
todo: [72, 73]
---
# FD — Three directions and their interfaces

The user requested preservation of the factor-state direction on arb, deeper
brainstorming of three alternatives and possible interleaving. Decision
M8c90bbdf49a94689 in topic:rsa72 preserves TODO72 and its C118/C119 boundaries.
This round produced reviewed proposals, with no solver, experiment, measured
saving, factoring attempt or claim promotion. A moderate slowdown for lower
memory is acceptable; no numerical tradeoff has been selected.

## Independent slate

Two fresh contexts generated candidates without seeing each other's output.
The coordinator sealed a candidate privately (Me5a4c8b3dfdc48d1), revealing it
only after both froze (M6cbe829568884f55). Exact text, hash convention and merged
ranking are in out/rsa73/coordinator_sealed_reveal_v1.md and
out/rsa73/merged_slate_v1.md. This is procedural independence, not sandboxing.

| Candidate | Attribution | Assessment and owning record |
|---|---|---|
| Build one GNFS matrix cache with less temporary storage | surveyor G; deriver A2 | first practical source audit; TX44 |
| Schedule smooth-part batches under a memory cap | deriver A1; surveyor G+A; coordinator overlap | secondary target if this stage is the relevant peak; TX44 |
| Block arithmetic collision aggregation | surveyor A; deriver B baseline | established-mathematics memory baseline; TX45 |
| Generate actual giant steps through short recurrences | deriver B | concrete algebraic bet; run-cover/seed lemmas missing; TX45 |
| Jointly dissect neighboring hyperbolas to count factor intervals | surveyor C | speculative upside; cancellation/region bound missing; TX46 |
| Use progression-product gcd certificates | deriver C; surveyor C+A independently | sound interface, not a branch-count oracle; TX47 |

The coordinator independently emphasized typed arithmetic services, construction
cost and practical baselines. The workers supplied the specific cache-build,
shift-recurrence and paired-strip proposals. Keep all three broad directions
alive. Their ranking was shown to the user before any follow-on derivation or
experiment assignment; none was dispatched. TODO73 owns the deciding steps.

## Combinations and rejected shortcuts

Practical interleaving shares memory/recompute accounting and, where operations
match, product/remainder-tree scheduling. Smooth parts consume explicit integer
norms and a prime set; collision search consumes polynomials and evaluation
points modulo N. Neither replaces GF(2) matrix-cache construction.

The strongest speculative combination is certified geometric rejection followed
by bounded-memory arithmetic search in disjoint unresolved intervals. Charge
total unresolved work, repeated queries and gcd=N recovery. Positive product
gcd can find a global factor without certifying a nonempty local factor branch
(TX47). Arithmetic progressions are not geometric runs: TX45's identity cannot
be copied to the interval interface without a new derivation.

A complete arithmetic or GNFS factorization can construct the tiny final factor
state afterward; C119 owns the discovery/access bill. Partial GNFS relations
are not conditional factor counts. Combining methods gives no automatic product
of speedups, and a speculative prefilter needs a cutoff charging failed work.

Generic compression, ordinary product trees and sequential cache building are
already available. Checkpoint deletion does not automatically save RAM. Compact
final states, easy pointwise multiplication and cheap modular-inverse sampling
do not supply efficient interval-conditioned construction. C118's settled
low-bit frontier is not scheduled again. These are scoped rejections, leaving
new recurrences, certificates and allocation schedules open.

## Evidence and review

| Work | Frozen submission | Review |
|---|---|---|
| Structural slate | S7e2bdedf5f784fed | V9975e0f543ff4ce1 |
| Source slate | S3f12499d0ee74a57 | Vfbc3ad231e824150 |

The referee was fresh for the first slate and disclosed that exposure for the
second. Acceptance is for proposals, including explicitly open lemmas, not
evidence of savings. Coordinator integration prose is a separate write task,
Td9a65c91cfe94de8, with a fresh integration referee.

Frozen reports contain primary-source body locations; TX44-TX47 retain the
dictionaries and hypotheses. The coordinator additionally inspected CADO tree
70354d7a8d54e985e46ca0fb6fb64d10716ba8bc; snapshots/hashes and limitations are
in out/rsa73/primary_manifest_v1.json and coordinator_source_audit_v1.md.
Upstream README estimates are not measurements. Newer balanced-semiprime and
order-finding preprint leads remain in TX45; the 2021 deterministic baseline
is not called the latest result or a practical GNFS competitor.

Only evidence lint, generated indexes and documentation checks apply. The
initial transfer index was frozen byte-for-byte before adding rows; existing
scientific input records remain unchanged. Unrelated dirty files are checked
against out/rsa73/preintegration_hashes_v1.json. Board/out evidence is local
and gitignored; canonical records retain the substance independently of board
resolution. No science suite, commit or publication occurs in this round.
