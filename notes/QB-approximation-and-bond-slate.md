---
code: QB
title: "Approximation versus exact representation cost: ideation and reviewed fixed-cut MPS obstruction"
outcome: "proved within the selected scope; C113 records the result"
claims: [C19, C107, C108, C109, C110, C111, C112, C113]
todo: [66, 69, 70]
---
# QB — From the swarm API check to an approximation obstruction

The user asked to verify familiarity with the current arb API and ideate a
useful, feasible research task. At the ideation stage this was a single-session proposal, not an
independent phase-0 swarm slate or an execution round: no claims or board
assignments were created then. The subsequently authorized round is recorded
below; no scientific code was run in either stage.

## Candidates retained

| Candidate and generator | Goal and consequence if successful | Strongest obstruction or limitation | Cheapest deciding step |
|---|---|---|---|
| Approximate MPS bond from arithmetic collisions; change the question plus obstruction first (TX34, TODO 70) | Determine whether approximation changes representation cost, using the actual physical state | Exact scalar rank and Pauli-LIMDD node counts do not bound the relevant Schmidt tail; general efficient output sampling faces TX15 | Derive the purity rectangle identity at one specified mixed cut and identify the arithmetic estimate it needs |
| Clean-block compression with interference; subspace structure (TX9, TODO 66 and TODO 64) | Remove the scratch-related exponential from a useful reduced representation | Already proposed in VR; a noncommuting kick can leave the x<N domain on which clean return holds | Check invariance/leakage of the exact kick fixture before adopting a reduced block |
| Interleaved Pauli-LIMDD obstruction; extend the known classification (TX29, TODO 69) | Close the gap between the two extreme register orders | C112's last-exponent count can become small; another level must carry the cost | Pair its invariant with a second cut and seek a counterexample to a proposed tradeoff |

The first is recommended for its combination of a new output-accuracy question
and an exact, bounded mathematical entry point. The second is a practical
fallback but should not be repackaged as a new proposal. The third remains
valuable, but an all-order classification is a less bounded first target.
No candidate is declared refuted by this ranking. TODO 70 owns the recommended
contract and stopping rules; TX34 owns the proposed external transfer.

The distinction between exact rank and approximation was checked against the
source bodies linked in TX34, rather than inferred from the rank headlines.
The key scope choice is the normalized clean state: an average error over the
truth table of a free work input would not certify Shor's fixed-input output.
Even a full-state approximation obstruction would cover only that simulation
route, not every possible output sampler.

## arb API verification

The installed executable reports 0.1.0 and is an editable installation of
`/home/djneko/Workspace/agent-research-board`. Live operation schemas and its
`docs/AGENT_GUIDE.md` were checked for task creation, claiming, heartbeat,
closure, submissions, review, run start/finish and inbox/resume/ack. Source
inspection confirmed recursive dependency freshness, frozen reads with
`candidate_changed`, and the stronger checks on accepting a review.

`arb call` currently supports `--file` (including `--file -`), inline `--data`,
`key=@file` and `key:=@file`; `--data @file` is not file loading. Resume and inbox
create delivery receipts, while acknowledgment advances the consumed cursor.
`run.start` records intent and does not execute a command. These are observations
of the installed API, not new protocol rules.

`uv run python tools/test_swarm.py` passed all ten tests against that executable,
including the disposable-board lifecycle and read-only-state checks. Live
board observation showed no open assignments, questions or held resources;
this dated observation does not replace the board as the owner of live state.

## Authorized derivation round (2026-09-21)

After the user approved the proposal, two independent workers received its
fixed contract: a derivation and a primary-source/hypothesis audit. Neither worker
received the other's results before both reports were frozen. Both independently
found the same simplification. Maximum row and column degrees bound the largest Schmidt weight,
so the proposed fourth-order estimate was unnecessary. A separate interval
Fourier argument connected the existing full-frequency character-sum hypothesis
to this state. C113 owns the theorem, proof, constants and all scope limits.

The derivation additionally explained the arithmetic rectangle structure using
a multiplicative determinant and the sizes of the low/high work coordinates.
The source audit checked approximation conventions and identified the printed
source-equation caveat now recorded in TX34. Neither worker ran numerical
science. Independent fresh referees reviewed the frozen mathematical reports;
the canonical claim then received a separate integration review.

The useful choice was to ask about approximation of the actual normalized
state at one explicit mixed cut. The initially anticipated character-sum
calculation and optional spectrum sweep were unnecessary once the simpler proof
settled the bounded question. An all-order theorem, or an inference from exact
rank to approximation cost, would have outrun the evidence. No claim of novelty,
a general sampling obstruction, or measured simulation performance is made.

TODO 70 is complete at its proof stopping point. The optional diagnostic was
not launched because it would not change the mathematical decision. Further
Schmidt-tail optimization or different orders would be separate research, not
unfinished work in this round.

Evidence and independent review:

| Deliverable | Frozen submission | Accepting fresh review |
|---|---|---|
| Derivation, including rectangle structure and degree/Fourier proofs | `Sefc89aee03564e66` | `Vf26584b6f448426a` |
| Independent primary-source audit and degree/Fourier proofs | `S11c34596ca6840e6` | `V5da1f5062b2b4693` |
| Canonical C113 integration | `Sc4aa448e3fc04871` | `V6e531a1c9f834029` |

The claim contains the self-contained proof and report hashes. These three
reviews required no scientific corrections. Board acceptance and canonical
promotion were separate steps, and the proposal records were updated only
after the producer and dependency-bound integration tasks had closed.
