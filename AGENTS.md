# Codex entrypoint for the research workspace

Read [CLAUDE.md](CLAUDE.md) for repository conventions and validation rules;
it remains the shared instruction source. Read [METHOD.md](METHOD.md) before
designing or running research, and [HANDOFF.md](HANDOFF.md) when resuming it.
Current claim status and evidence live in the individual `claims/` records.

For authorized swarm work, read [SWARM.md](SWARM.md), the shared Codex/Claude
protocol for future assignments, independent references/referees, dependencies
and evidence lint. Installed host skills are adapters to that protocol.

Use the installed `qsim-research` skill for this project's research discussion,
claim review, experiments, and paper refinement. It is a workflow guide, not
a second ledger or an instruction to continue research during a review-only
request. If skill discovery is unavailable, its source is
`$CODEX_HOME/skills/qsim-research/SKILL.md` (default
`~/.codex/skills/qsim-research/SKILL.md`); read it directly when relevant.

Match validation to the requested work as specified in `CLAUDE.md`. Do not
start science suites for a short conceptual answer. Preserve the current dirty
worktree; do not commit or publish unless requested.
