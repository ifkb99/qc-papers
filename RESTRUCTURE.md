# Restructure plan — context, drift, and the shape of the repo

*Written 2026-08-08, after a session that fixed sixteen consistency defects
across the two paper drafts and the ledger.*

> **STATUS: EXECUTED 2026-08-08.** All four phases are done; this document is
> kept as the rationale, not as a to-do. Outcomes against §10's criteria:
>
> | criterion | target | measured |
> |---|---|---|
> | cold-start orientation | < 5,000 tok | **4,912** (was 16,769) |
> | working one item | ≈ 9,000 tok | **≈ 9,000** |
> | `check.py` fails on a reintroduced defect | must | **verified on three** |
>
> Counts landed as predicted for claims (44 live + 19 retracted) but not for
> the other two: TODO had **29** items, not 20, because nine were loose
> Housekeeping bullets; NOTES had **42** sections, not 43, because one heading
> is wrapped across two lines. Both are recorded in the migration commits.
>
> **Two decisions taken during execution that this plan did not anticipate.**
> `CLAIMS.md` and `NOTES.md` were kept as *generated aggregates* rather than
> deleted — the papers cite `CLAIMS.md` as supplementary material, and roughly
> a hundred references to `NOTES.md` live in claim files and in historical
> experiment scripts that are records and should not be edited. Regenerating
> costs nothing and every existing pointer keeps working. And `check.py` is
> deliberately **not** one of the nine `test_*.py` suites: it gates
> documentation, not science, and folding it in would change what "nine-suite
> correctness gate" means in the papers.

**Operating assumption, decided 2026-08-08: sessions are SEQUENTIAL, one at a
time.** No two agents work on this repo simultaneously. That decision removes
a large amount of machinery from the design — there is no ownership protocol,
no locking, no claim/release convention, and no post-merge index regeneration
below. If that assumption ever changes, see "If parallelism becomes real" at
the end; the layout proposed here gives conflict-freedom as a side effect, so
the change would be additive rather than a redesign.

---

## 1. Diagnosis: the problem is duplication, not size

The TODO item that prompted this reads *"TODO, claims, and even handoff files
are becoming quite large."* They are. But size is the symptom.

Take four headline numbers and ask which files contain them:

| number | files containing it |
|---|---|
| `0.716` | HANDOFF, CLAIMS, ABSTRACT_SHOR_2ADIC, PAPER_A, PAPER_B, TODO, NOTES — **7** |
| `1.006` | HANDOFF, CLAIMS, PAPER_B, PAPER_A, NOTES, TODO — **6** |
| `24369` | CLAIMS, ABSTRACT_SHOR_2ADIC, PAPER_B, NOTES, TODO — **5** |
| `15549` | CLAIMS, ABSTRACT_SHOR_2ADIC, PAPER_B, NOTES — **4** |

Each is one fact stored 4–7 times, with **no mechanism keeping the copies
equal**. The files are large *because* facts are duplicated, and duplication is
also what makes them wrong.

### The evidence that this is not theoretical

The 2026-08-08 consistency pass found **sixteen defects**. Essentially all of
them were the same shape — *a fact was updated in one place and not its
mirror*:

| defect | copies that drifted |
|---|---|
| `rot = 2·perm` described as "an exact factor of two" after Proposition 2 superseded that framing | 4 (§1, §11.2, §12, and the README) |
| correctness gate described as "eight-suite" after `test_accel.py` made it nine | 4 (PAPER_A, PAPER_B, `test_lab.py`, `test_claims.py`) |
| C11's speedup recorded as `~10³×` when Table 1 measures 143–219× | ledger vs paper, and **the ledger is canonical** — so the wrong copy was the authoritative one |
| §OS4 marked "RETRACTED" when only its cliff sub-observation died, while CLAIMS C7 cites §OS4 for a live result | ledger pointed at a section labelled retracted |
| C44 defined in the ledger, cited in neither paper | the headline result of a whole session was orphaned |

Splitting these files into a wiki *without* addressing duplication would
produce the same sixteen defects spread across eighty files instead of six —
harder to find, not easier.

### Current orientation cost

| file | lines | bytes | ~tokens |
|---|---|---|---|
| `NOTES.md` | 2790 | 141,925 | ~35,500 |
| `TODO.md` | 799 | 47,020 | ~11,800 |
| `CLAIMS.md` | 159 | 30,779 | ~7,700 |
| `HANDOFF.md` | 509 | 28,767 | ~7,200 |
| `METHOD.md` | 157 | 7,536 | ~1,900 |
| `README.md` | 112 | 6,331 | ~1,600 |
| **total** | | **262,358** | **~65,600** |

A cold session must read HANDOFF + CLAIMS + METHOD (~17K tokens) before it can
do anything at all, and usually several NOTES sections on top.

Two further measurements:

- **81% of `TODO.md` is closed items** — 644 lines of 790. Valuable (it is what
  stops re-derivation) but not hot-path.
- **A third of `CLAIMS.md` is the retracted section** — 52 lines of 159.

### The missing file

**There is no `CLAUDE.md` anywhere in this repository.** Every session starts
with zero project instructions and reconstructs the conventions from
`HANDOFF.md` — which is exactly why HANDOFF has grown to 509 lines carrying
traps, GPU capacity notes, reading order and a historical thread. Roughly half
of HANDOFF is content that should be auto-loaded instruction, not prose a
session has to find and read.

Separately: `.claude/settings.local.json` (in the parent directory) contains
permissions for an unrelated Java project — `java T.java`, `probe.Runner`.
Leftovers. Nothing there helps this project.

---

## 2. Principles

1. **One fact, one home.** Every number, status and claim lives in exactly one
   file. Everything else points at it or is generated from it.
2. **Indexes are generated, never hand-written.** A hand-maintained index is
   just another mirror, and mirrors drift. If it can be derived, derive it.
3. **Consistency checks are executable.** The six checks written ad hoc during
   the 2026-08-08 pass become a permanent suite. Cleaning up drift once is
   worth little; preventing its recurrence is the actual deliverable.
4. **The historical record is preserved but leaves the hot path.** Retracted
   claims, closed TODOs and superseded notes are what stop a fresh session
   re-deriving dead ends. They must stay greppable and must stop being resident.

---

## 3. Target layout

```
research/
  CLAUDE.md                  NEW — auto-loaded conventions, traps, invariants
  HANDOFF.md                 slimmed to ~80 lines: current state + pointers
  METHOD.md                  unchanged (157 lines, always relevant, no duplication)
  README.md                  unchanged in role; trimmed of duplicated status

  claims/
    C8.md  C15.md  …  C44.md      one claim per file, YAML frontmatter
    retracted/C3.md  F2.md  …     dead claims, same schema, out of hot path
    INDEX.md                      GENERATED — id | one-line | status

  notes/
    PK-peak-deficit.md            one investigation per file
    GF-gf2-pattern.md             keyed to the existing §codes so §PK, §GF still resolve
    …
    INDEX.md                      GENERATED

  todo/
    open/12c-early-pruning.md     one item per file
    open/14-inverse-qft.md
    done/12e-gpu-sweeps.md        state = which directory
    INDEX.md                      GENERATED

  tools/
    reindex.py                    regenerates every INDEX.md
    check.py                      the consistency gate — suite 10

  PAPER_A.md  PAPER_B.md          unchanged, whole
  ABSTRACT.md  ABSTRACT_SHOR_2ADIC.md   unchanged (abstract workshops / prior-art dossiers)
  lab/  experiments/  archive/  out/    unchanged
```

Counts, measured rather than estimated:

| directory | files | source |
|---|---|---|
| `claims/` | **44** | 42 identified live rows (C1…C44 with gaps, plus F12) + 2 unnumbered |
| `claims/retracted/` | **19** | 12 identified (C3, C4, C9, F1–F8, H1, H2) + 7 unnumbered |
| `notes/` | **43** | `## ` sections in `NOTES.md` |
| `todo/` | **20** | `## N.` items in `TODO.md` |
| | **126 total** | replacing three large files |

**Migration wrinkle worth knowing before starting: nine rows have no
identifier.** Two live (the biased-input "dictionary" row, the superseded
residual-ladder row) and seven retracted (the "earlier guess" rows, the
original Paper A thesis, the weight-truncation row). They need names assigned
before they can become files — `claims/dictionary-biased-inputs.md` and so on.
Do not let the migration script silently drop them; they are cited in prose.

---

## 4. Frontmatter schemas

### `claims/C44.md`

```yaml
---
id: C44
status: derived            # proven | derived | established | narrowed | retracted
paper: A                   # A | B | none
one_line: >
  The peak ratio is exact, not empirical: N_max^rot = 2·N_max^perm − |B|.
experiments: [experiment_c17_deficit]
notes: [PK-peak-deficit]
supersedes: []
unproved: ["|B| = 2 itself"]
---

## Statement
…

## Evidence
…

## Mechanism
…
```

`one_line` is what the generated index shows. `unproved` is deliberately a
field rather than prose: it is the thing most likely to be lost in a rewrite,
and it is what keeps the honesty record honest.

### `notes/PK-peak-deficit.md`

```yaml
---
code: PK                   # the existing §code, so old cross-references resolve
date: 2026-08-08
todo: 12d
outcome: solved            # solved | negative | retracted | partial
claims: [C44]
supersedes_in_this_file: []
---
```

### `todo/open/14-inverse-qft.md`

```yaml
---
id: 14
state: open                # open | done | dead
rank: 1
title: Through the inverse QFT
outcome: null              # filled in on close
claims: []
---
```

No `owner` field. Sessions are sequential; ownership machinery would be
ceremony with no counterparty.

---

## 5. `tools/check.py` — the gate

This is the part that matters most. Each check below is listed with the actual
defect from the 2026-08-08 pass it would have caught, because a check that has
never caught anything is decoration.

| # | check | would have caught |
|---|---|---|
| 1 | every claim ID referenced anywhere resolves to a claim file | C44 orphaned — defined, cited nowhere |
| 2 | no paper cites a retracted claim as support | preventive |
| 3 | appendix claim maps match real section headings | stale appendix rows after §9 was cut |
| 4 | every `§N` cross-reference resolves to a heading | — (clean, but cheap insurance) |
| 5 | no lab notation (`§PK`, bare `C15`, `TODO 12`) in paper prose | the `§I` leak; the `**On C1/C6.**` heading; the `(TODO 12)` comment |
| 6 | numbered statements are unique and sequential | **two different results both numbered "Proposition 2"** |
| 7 | reproducibility commands resolve to files that exist | Paper A had no repro block at all |
| 8 | `INDEX.md` files match a fresh regeneration | index drift |
| 9 | suite count in prose matches `ls test_*.py` | "eight-suite" after the ninth was added |
| 10 | a headline number appearing in two files must appear in a claim file | C11's `~10³×` against Table 1's 143–219× |

Checks 6 and 10 are the highest-value ones and neither is obvious in advance:
6 catches a class of error that is invisible on a section-by-section read, and
10 is the direct answer to the duplication table in §1.

`check.py` runs as **suite 10** of the correctness gate, alongside the nine
existing suites. It must exit nonzero on failure so it cannot be ignored.

---

## 6. `CLAUDE.md` — what goes in it

Auto-loaded, so it must be short and must contain only what *every* session
needs. Target: under 120 lines.

- **Run everything with `uv run python` from `research/`.** Not bare `python3`.
- **Use 3.12/3.13 for long `perm_pps` jobs** — CPython 3.14 dies intermittently
  with `_TAIL_CALL_CACHE`.
- **`LAB_GPU=1` or it silently runs on CPU at a tenth the speed.** Capacity
  n ≤ 30 on one 20 GiB card.
- **The one-fact-one-home rule**, and where each kind of fact lives.
- **Indexes are generated. Never hand-edit an `INDEX.md`.**
- **Every experiment carries a must-fail control**, and predictions are written
  before measurement. This rule has caught two vacuous results.
- **A unanimous pass is a tell, not a triumph.**
- **Run `tools/check.py` before committing documentation changes.**
- Pointers: `HANDOFF.md` for state, `METHOD.md` for how the work is conducted,
  `claims/INDEX.md` for what is known.

Everything else currently in HANDOFF — GPU benchmark tables, the historical
resolved thread, reading order, the do-not-redo list — moves to
`notes/` or stays in a slimmed HANDOFF as pointers.

---

## 7. Migration phases

Ordered so that value lands early and the risky part comes last.

### Phase 1 — cheap, high value, no migration risk

1. Write `CLAUDE.md`.
2. Write `tools/check.py` with checks 4, 5, 6, 7, 9 (the ones that need no new
   file layout), wire it in as suite 10.
3. Split `TODO.md` into `todo/open/` and `todo/done/`; write `tools/reindex.py`
   and generate `todo/INDEX.md`.

Effect: 644 lines leave the hot path; the drift class becomes detectable; every
future session starts oriented. **`CLAIMS.md` and `NOTES.md` are untouched.**

### Phase 2 — claims

4. Split `CLAIMS.md` into `claims/*.md` + `claims/retracted/*.md`, generate
   `claims/INDEX.md`.
5. Enable checks 1, 2, 8, 10.
6. Update the papers' Appendix A tables to point at claim files.

Effect: the canonical ledger becomes selectively loadable, and check 10 —
the direct answer to the duplication problem — switches on.

### Phase 3 — notes

7. Split `NOTES.md` into `notes/*.md` keyed by existing §codes.
8. Generate `notes/INDEX.md`.

Largest file, lowest urgency: NOTES is rarely read whole, and the §codes
already function as an index of sorts.

### Phase 4 — handoff

9. Slim `HANDOFF.md` to current state plus pointers, once there is somewhere
   for the rest to go.

---

## 8. Deliberately not doing

- **Not splitting the papers.** They are single documents with a single
  authorial voice; a paper assembled from fragments reads like one.
- **Not splitting `METHOD.md`.** 157 lines, always relevant, contains no
  duplicated facts. It is already the right size.
- **Not running a wiki server.** A directory of markdown plus generated indexes
  gives every benefit that matters here — searchable, diffable, reviewable,
  testable by the existing gate — with no infrastructure and no new failure
  mode. "Wiki" was the right instinct; a server is not the right implementation.
- **No ownership or locking protocol.** Sessions are sequential.
- **Not deleting anything.** Retracted claims, closed TODOs and dead ends move;
  they do not disappear. They are the most expensively acquired knowledge here.

---

## 9. Risks

| risk | mitigation |
|---|---|
| The migration itself introduces the drift it prevents | Migrate mechanically with a script, not by hand. Assert counts before and after (44 live claims, 19 retracted, 43 notes, 20 todos). Run `check.py` immediately. |
| The nine unnumbered ledger rows are silently dropped | Assert the retracted-row count explicitly; they are the easiest thing to lose and are cited in paper prose. |
| `git log` on `NOTES.md` becomes harder to follow | Accept it. The history stays in the repository; the file stops being the working surface. |
| 105 files is more `ls` noise than 3 | Generated `INDEX.md` in each directory is the entry point, not the file listing. |
| Frontmatter rots the way prose does | It is machine-read by `reindex.py` and `check.py`, so rot fails a test rather than sitting unnoticed. |
| Effort spent on tooling instead of research | Phase 1 is the only part that must happen. Phases 2–4 are each independently useful and can stop at any point. |

---

## 10. Success criteria

- A cold session reaches "I know what is true and what I am doing" in **under
  5K tokens** (`CLAUDE.md` + slim `HANDOFF.md` + `claims/INDEX.md`), against
  ~17K today.
- Working one TODO item costs roughly **9K tokens** all in — the item file,
  three to five claim files, one or two note files.
- `tools/check.py` passes, and **fails when a fact is edited in one place and
  not its mirror.** This is the real test: if it cannot be made to fail by
  reintroducing one of the sixteen defects, it is not doing its job.

---

## 11. If parallelism becomes real

The layout above is already conflict-free — two sessions adding `C45.md` and
`C46.md` touch different files and git merges them cleanly. What would need
adding, and nothing more:

- an `owner` field in TODO frontmatter, set in its own commit before work
  starts;
- regenerate indexes after merge rather than in the working session;
- run `check.py` before every commit rather than before documentation commits.

No part of the layout would change. That is the reason to prefer one-fact-per-
file even under the sequential assumption: the property is free, and acquiring
it later would be a migration.
