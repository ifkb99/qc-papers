---
code: DS
date: 2026-09-12
title: "Dirty-prefix reversal symmetry redirects the scalar; the disabled sector has a polynomial construction"
outcome: exact three-macro full-space discriminator found; growing disabled sector reduced; enabled closure remains open
claims: [C92]
todo: [50]
---
# DS — A full-space discriminator after a symmetry exclusion

C92 owns the exact symmetries, all-width three-macro map witness, Mersenne
prefix identities and constructive symbolic bounds. TODO50 owns the next
enabled-sector question. This investigation makes progress, without a general
simulation or CNOT-memory breakthrough; the broad user goal remains active.

## Derivation and changes of direction

The authorized board swarm retained Astra for algebra and Sol for source
and independent test work. Root remained the only canonical writer. Starting
from C91, root and Astra derived a small correction count for Mersenne
prefixes. That count leaves correlations with the scratch, enable word and
flag history, so it was not treated as a compact scalar state.

Astra derived a common reversing involution. Root gave a shorter direct
proof and the workers independently checked it. This rules out complete
reversal as a discriminator for equal low-bit/flag endpoints at any depth,
even before scratch averaging. It does not rule out adjacent swaps. An
intermediate suggestion that odd scratch/carry parity might make such a
full coefficient reversal-sensitive was eliminated by the independent
global-complement symmetry; no experiment was run on that suggestion.

Root and Astra then constructed a polynomial affine partition for the
disabled common-control sector, with full dirty-scratch averaging. Root's
canonical transfer uses open vertical strips plus separate integer critical
slices to make threshold boundaries explicit. This stronger growing-depth
baseline remains unimplemented and intentionally unoptimized. It does not
remove the enabled sector's nested dyadic prefix terms.

Sol's source audit constructed a native bit-column carry baseline and checked
the strongest applicable published constructions. Root inspected the five
actual add/subtract calls and the cited primary theorem/algorithm bodies.
C92 charges state construction, coefficient bits and extraction, and compares
with C89's polynomial-memory interval enumeration. Minimal automata or compact
word expressions are not assumed to be freely supplied representations.

The other Sol worker supplied an explicit dirty input at every Mersenne
width n>=3 whose three-macro output character changes under reversal and
omitted incoming carry. Walsh invertibility and Parseval prove that separate
nonzero chronological coefficients distinguish the two references. This
justified a bounded MASK DISCOVERY before a separate frozen-mask verification;
it did not predict the discovered mask or coefficient magnitudes.

## Bounded discovery and exact values

Source: `experiments/experiment_dirty_prefix_walsh.py`. Reproduce with:

```bash
timeout 60s uv run python -u -X faulthandler -m experiments.experiment_dirty_prefix_walsh
```

Task `T6f4e677c084545a9`, attempt `A2119f77916944163`, prediction
`Mb738ced3a45741a1` preceded run `Rb2dceaef08e04fe9`. It exited zero;
all eight checks passed without warnings. The complete report is
`out/dirty_prefix_walsh_report.json`, with raw log
`out/dirty_prefix_walsh_main.log`.

The fixture is N=7,a=1,n_exp=1,q=3, constants (1,2,4), with fourteen
physical wires and 16,384 uniformly averaged input labels. Both chronological
and reversed circuits have 426 logical operations. Reversal retains each
constant's own multiplicand control. C89's independent word map agrees with
both gate maps on all 32,768 compared labels. Both maps satisfy every P/J
orbit check; the two selection rules leave 4,096 candidate input masks.

For the output character C=b0+f, the ascending-mask search found:

| input character | physical mask | chronological | reversed | incoming h omitted |
|---|---:|---:|---:|---:|
| b1+b3 | 10 | -19/512 | -15/512 | -19/512 |
| b3+h | 2056 | -15/512 | -13/512 | 0 |

The output physical mask is 4097. There are 3,138 chronological-nonzero
order witnesses, 1,804 carry witnesses and 1,570 common witnesses within
the allowed set. The second row is the first carry and first common mask.
A common mask was an open discovery branch, not guaranteed by the proof.
Message `M35a79c5f21434060` froze both rows before independent verification.

The full average includes every b,t,x,h,f,u assignment. Omitting h changes
only the arithmetic, while retaining its physical wire and input phase.
The chosen common mask has no x character. In u=0 all macros are the same
map, so its forward/reverse difference comes entirely from u=1; its scalar
value can still contain a disabled-sector contribution.

The reference uses existing `walsh.classical_permutation` on logical
X/CNOT/Toffoli operations. Existing `walsh.wht` stores float64 butterflies,
but this use is exact: all inputs are integer signs and each intermediate
is an integer of magnitude at most 16,384, below 2^53. Integer-conversion,
Parseval and balance checks passed, and the selected sums agree with direct
integer character sums. This is not a floating-point tolerance result or an
independent unitary verification of the stored Clifford+T decomposition.
No statevector, dense unitary, size sweep or comparative memory/timing
measurement was made. Harness elapsed times are descriptive only.

## Independent interval contraction and retained verifier correction

Task `T0596b7520d904530`, attempt `Af53b117e11f341f5`, uses only the
existing C89 `compile_intervals` and `interval_coefficient` helpers. It
contracts all 32 b,f states through exact range-character sums for each of
512 outside fibers. It uses no gate replay or WHT and retains every fiber
in its JSON report. It is an independent contraction route, not an independent
derivation of the already-reviewed C89 helpers.

The initial prediction `Me2a13c1257eb4e1c` preceded run
`R610eb83d46024315`, terminal exit zero, six checks passed. Its submission
`S52fc7675ddec4ca6` reported both frozen triples exactly. Root then found a
latent verifier helper error: `outside_sign` added bit counts and took the
bit parity of that COUNT rather than the count modulo two. For these frozen
masks the count is only zero or h, so both reported triples remain correct.
Counts two and three would have been wrong. Review `V87ef380ebaa34da8`
requested preservation, repair and an independent physical-mask sign check
before a same-fixture rerun. This is a verifier-generalization defect, not
a failed scientific prediction or a reason to discard the initial evidence.

Revision 1 source/report/log/review are preserved in the attempt's
`revision1/` directory. The corrected prediction `M84586123f81d4148`
preceded run `Rc59d17f61025445c`, terminal exit zero, seven checks passed
without warnings. New count-two/count-three cases reject the old helper;
the unchanged frozen triples and both wrong-reference controls pass.
All three variants have a maximum of 26 interval pieces in this fixture.
The accepted revision is `S282ecdab7d9f42d0`. Its authoritative files are
`verify_dirty_prefix_intervals.py`, `verify_dirty_prefix_intervals_report.json`,
and `verify_dirty_prefix_intervals.log` in
`out/agent-board/workers/Af53b117e11f341f5/`.

Root read the full original source, exact correction and full logs, checked
the archived hashes, then parsed every final fiber row. The fiber set is
complete without duplication, all enabled controls are correctly paired,
disabled-fiber forward/reverse results agree, and root's re-summed Fraction
totals match both the report and the independently constructed gate sums.
No additional size experiment was needed after these checks.

## Proof and source evidence

| board task | accepted submission | accepting review |
|---|---|---|
| T5157b49f7df64a4c, Astra algebra | Sc7d13ba772164ac6 | Vc0be156ce73743b5 |
| Td13dad1c7f324a9d, Sol symbolic baseline | Seeeac35b283e4023 | V0ebb930df3ef462d |
| Tf4de407e93bd44f0, Sol map discriminator | S73b4832414f14243 | V977c7261e1ac4566 |

Proof: `out/agent-board/workers/Ad4bfaca3d7054c0f/mersenne_prefix_audit.md`.
Source audit: `out/agent-board/workers/A045d055712424a18/dirty_mersenne_symbolic_baseline_audit.md`.
Discriminator: `out/agent-board/workers/A693bffaa51a7437b/q3_dirty_prefix_symmetry_proposal.md`.
Root read each complete artifact, inspected its mathematical and native-code
reasoning, and checked current files against the archived SHA256 identities.
The proposal's redundant old two-macro fixture was not repeated: the new
experiment has its own exact normalization, direct-sum and wrong-reference
checks. None of these proof tasks executed scientific code.

Astra and Sol independently reviewed the complete canonical C92 transfer
after integration, recorded as `Mcd4e617b366b4c99` and `Ma067f4ca2cc54058`.
Neither requested a correction. The direct reversor proof and more explicit
integer-boundary partition retain the submitted proof's scope. Root's write
task includes the new experiment and final claim/note/backlog/handoff records
for a separate source-and-evidence acceptance review.

## Validation scope

Root ran core first, observed terminal exit zero, and inspected the full log
`out/mersenne_core_initial.log`. The independent worker also ran core before
its verifier, with its separate log in the attempt directory. The new
bounded experiment is the affected scientific check; production helpers
were unchanged. Other eight science suites were not rerun. No timeout or
native crash occurred in these bounded runs; TODO34 retains the prior issue.
No manuscript, abstract, dependency, default or host setting changed; nothing
was committed or published.

Index generation and all ten documentation checks pass, recorded in
`out/dirty_prefix_reindex.log` and `out/dirty_prefix_docs_check.log`.
New experiment syntax and scoped whitespace checks pass. The board correctly
rejected initial closure of three proof tasks after root updated their
snapshotted TODO50 input. Those tasks were reopened for fresh input snapshots
and read-only revalidation of their unchanged evidence; no scientific rerun
or retraction was required. The independent verification task is accepted
and closed. Final proof/write-task acceptance and closure remain in the board
history rather than duplicating version identifiers here.
