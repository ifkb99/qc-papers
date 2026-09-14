# CLAUDE.md — read this before touching anything

Research repo: Pauli-path simulation (PPS) of reversible quantum arithmetic.
Two papers are drafted; the experimental backlog is largely closed.

**One coordinating session and one canonical source writer at a time.**
Authorized delegated workers coordinate through the local `arb` board configured
by `arb.toml`. The coordinator edits shared source files; workers put drafts and
evidence in their assigned `out/agent-board/workers/` directory. Board ownership
is advisory: this host does not enforce filesystem isolation.

On research resumption, run `arb --actor coordinator resume` from this directory;
workers use their assigned actor ID. Read the inbox and acknowledge its receipt
after consuming it. `arb schema OPERATION` gives the agent interface; `arb board`
and `arb thread topic:general` provide a human-readable entry. The installed
research skill supplies the coordination workflow. Model choice remains with
the existing session, and board use does not imply permission for new delegation.

---

## Running things

```bash
uv run python test_core.py        # science gate — run first for science work
uv run python tools/check.py      # documentation gate — before any doc commit
```

The nine `test_*.py` suites gate the science. `tools/check.py` is separate and
gates the documentation; it is deliberately **not** one of the nine, so that
"nine-suite correctness gate" in the papers keeps meaning nine science suites.

Run the core gate before changing or executing scientific code/experiments,
then select affected suites in proportion to the change. Read-only discussion
and document/skill-only edits do not require science execution. Documentation
changes use the documentation gate below; skill changes also need their own
validation. Do not interpret a passing science suite as a skill-behavior test.

1. **`uv run python` from `research/`.** Not bare `python3`.
2. **`LAB_GPU=1` or it silently runs on CPU at a tenth the speed.**
   `accel.enabled()` says whether it took effect. Capacity n ≤ 30 on one 20 GiB
   card; it raises rather than thrashing past that.
3. **Use Python 3.12/3.13 for long `perm_pps` jobs.** CPython 3.14 (the pinned
   version) dies intermittently with `Fatal Python error: _TAIL_CALL_CACHE` —
   an interpreter bug, not this code. Presents as a hang or a death with no
   traceback.
4. **Don't pipe a long experiment through `tail`** — the pipe buffers the whole
   run and it looks hung. Redirect to a file.
5. **Don't `pkill` in the same command as a heredoc write** — it kills the
   write. Two scripts were lost this way.
6. **Don't build large `frozenset`s of tuples over a support.** At 33.5M
   elements that is multi-GB and an OOM. Encode as int64 and compare sorted
   numpy arrays.

Some experiments **exit nonzero on purpose** because their predictions were
refuted and they were left failing rather than re-scoped
(`experiment_c7_deficit`, `experiment_windowed` P5). Read the header before
concluding something is broken.

---

## Where facts live — one fact, one home

| kind of fact | home | never |
|---|---|---|
| what is known, and its status | `claims/<ID>.md` | duplicated into prose elsewhere |
| how something was found | `notes/<CODE>-*.md` | restated in a claim body |
| what to do next | `todo/open/*.md` | tracked anywhere else |
| current state, pointers | `HANDOFF.md` | facts that belong in a claim |
| how the work is conducted | `METHOD.md` | — |

Board assignments link to the existing research questions and record temporary
execution, ownership, discussion, submissions, and review. They do not duplicate
the scientific ledger or replace the research backlog. A submission accepted on
the board does not itself promote a scientific claim.

**`INDEX.md` files and `CLAIMS.md` are GENERATED. Never hand-edit them.**
Run `uv run python tools/reindex.py` after changing any claim, note or todo.

A hand-maintained mirror is the single failure mode this repo has actually
suffered: a 2026-08-08 audit found sixteen defects, essentially all of them a
fact updated in one place and not its copy. If you find yourself typing a
number that already exists somewhere else, stop and cite it instead.

---

## Research discipline

These are not decoration. They have caught real errors that had already
produced confident-looking numbers.

- **Every experiment carries a control that must FAIL.** It has caught two
  vacuous measurements. `lab.harness` warns when none is registered.
- **Write the prediction before measuring.** `lab.harness` enforces ordering.
- **A unanimous pass is a tell, not a triumph.** Three vacuous tests in this
  project passed unanimously and unearned.
- **Vary exactly one parameter.** A sweep that drew a fresh random table per
  width confounded width with table variance.
- **In research on your own code, a surprising result is almost always your own
  bug** — and a result *agreeing* with your hypothesis is more dangerous than
  one contradicting it, because you will not look.
- **Derive before measuring**, and promote regularities to proofs by reading
  the construction rather than measuring more.

Start new experiments from `experiments/TEMPLATE.py`. Reach for
`lab.measure`, `lab.gf2`, `lab.variants`, `lab.nulls` before writing new
machinery. Never write a second propagator to inspect the first — that has
bitten this project before.

The installed `qsim-research` skill guides Codex through these existing
instructions and the current ledger; `AGENTS.md` is its repository entrypoint.
It does not depend on a separate `computational-research` skill.

---

## Before committing documentation changes

```bash
uv run python tools/reindex.py    # regenerate indexes
uv run python tools/check.py      # must pass
```

`check.py` fails on: unresolvable claim IDs, papers citing retracted claims as
support, appendix maps that disagree with real headings, dangling `§`
references, lab notation in paper prose, duplicate or non-sequential
Theorem/Proposition numbering, reproducibility commands naming files that do
not exist, stale indexes, and a wrong suite count.

---

## Orientation

`HANDOFF.md` for current state · `claims/INDEX.md` for what is known ·
`todo/INDEX.md` for what is next · `METHOD.md` for how the work is conducted ·
`RESTRUCTURE.md` for why the repo is laid out this way.
