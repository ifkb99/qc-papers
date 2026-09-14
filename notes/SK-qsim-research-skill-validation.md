---
code: SK
date: 2026-09-10
title: "QSim research skill: task-specific workflow and isolated behavioral smoke tests"
outcome: validated
claims: []
todo: []
---
# SK — QSim research skill validation

The user requested a small research-alignment skill and a practical test.
This is tooling/workflow validation, not a new scientific claim or experiment
supporting either paper. Skill Creator guided authoring and the independent
forward tests; OpenAI Docs supplied the discovery/invocation conventions.

## Installation and scope

Source: `~/.codex/skills/qsim-research/SKILL.md`, with
`agents/openai.yaml` providing its UI name, invocation hint, and enabled
implicit-matching policy. The user-level discovery path
`~/.agents/skills/qsim-research` is a symlink to that single source, not a
second maintained skill copy. This follows the current
[skill-discovery documentation](https://learn.chatgpt.com/docs/build-skills).
The skill is installed outside the repository; it is not included in a commit
of this research tree.

The skill points to existing project documents and claim files. It does not
embed result numbers, claim statuses, or a fixed next research direction.
It distinguishes discussion, design, execution and write-up, asks for a precise
simulation task and matched baselines, and preserves scope for exploration.
There is no new simulator, general-purpose plugin, or additional research
framework hidden in the package.

`AGENTS.md` is the repository entrypoint to the existing `CLAUDE.md` and skill.
The test-gate wording in `CLAUDE.md` was clarified: scientific code work starts
with the core science gate; conceptual discussion and skill/docs-only changes
do not automatically execute science suites. This removes conflicting workflow
instructions without weakening the checks for scientific code.

## Tests

The supplied Skill Creator validator passed. Separate semantic checks verified
UI metadata constraints, the explicit invocation hint, enabled implicit policy,
and the discovery symlink target. The installed skill and test snapshot have
identical SHA-256:

`6aba6e5ce0b67c51869281cf4779a705f7b06f917f39413b0e052b9f862e5eb4`

Independent evaluators started without this conversation's history. Each read
the same frozen skill, then addressed a realistic request using a project
snapshot at `/tmp/qsim-skill-eval.H5QuBM/research`. They were given the task and
raw repository artifacts, not expected answers. The discussion/review trials
were read-only; the execution trial could create only its prefixed test files
and output artifacts. No evaluation changes were requested in the live papers,
claim ledger or scientific code.

| Trial | User task | Observed behavior | Assessment |
|---|---|---|---|
| Conceptual | Explain within 120 words whether scalar tensor rank one implies cheap Shor output sampling; no tests or edits | Distinguished scalar/joint ranks and sampling cost; charged order discovery; cited current claims; reported only read-only inspection | Pass |
| Claim review | Assess whether the recorded small-amplitude run establishes a memory advantage over classical Shor simulators; suggest wording without editing | Distinguished retained amplitudes from process memory, compared orbit and spectral baselines, charged discovery, checked primary literature, and noted absent raw JSON in the snapshot instead of claiming to reverify it | Pass |
| Reproduction | Test the proposed universal factor-two peak bound on the supplied tiny circuit, using the existing harness and propagators, with a control and saved report | Reproduced the existing counterexample with a same-layout single-Toffoli control; retained the false candidate prediction, saved JSON with `ok: false`, and reported exit 1; used existing engines and checked dense/Walsh references plus higher precision | Pass |

The first two trials cite the existing ledger and do not assert a new finding.
The reproduction is likewise a test of workflow behavior on an existing
research question, not a promotion of its hypothesis. A failed scientific
prediction can be the correct outcome of a passing skill-behavior trial.

The main agent inspected the reproduction source and report and reran it in
the isolated snapshot: exit 1, four successful checks and the one refuted
candidate bound, with no harness warnings. The negative-control condition
was resolved as expected. This is a successful failure-preservation test,
not an all-green scientific experiment. Existing witness results remain in
their original claims; no new scientific claim was created.

The live papers, claim directory, and inspected propagation/harness source
files were compared against the pretrial snapshot and were unchanged. Review
evaluators reported read-only operations; these reports are not a complete
host-level audit of every possible side effect. The third evaluator's source,
JSON and audit, plus the frozen skill, are retained under ignored
`out/skill_validation/`. The experiment source is an archival fixture: copy it
to a fresh snapshot's `experiments/` before rerunning its module command.

## Limits and reproduction

These are qualitative smoke tests, not a controlled with/without-skill study.
They cannot establish that the skill caused better decisions than an unassisted
agent, guarantee future correctness or discovery, or establish automatic
selection by the running UI. The trials explicitly supplied the skill path;
automatic matching is enabled in metadata, not behaviorally measured here.

To rerun, use fresh independent contexts with the three task descriptions
above and the installed skill. Use an isolated project snapshot for any
execution and inspect both the answer and persisted report. Do not include
this validation note or an answer rubric in the evaluator's task context.
The temporary snapshot is disposable and is not the authoritative repository.

The core science gate was run before the instruction clarification and passed;
no scientific implementation was changed for this skill. Documentation and
skill validation are the relevant final gates for this task; the supplied
skill validator, metadata/discovery checks and all ten documentation checks
passed. No full nine-suite science rerun was necessary for instruction-only
changes. No commits or publication actions were taken.
