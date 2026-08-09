---
id: hk9
state: done
title: Structure files to reduce their context
outcome: DONE 2026-08-08 — one fact one home, generated indexes, executable checks
claims: []
---
# Structure files to decreace their context

**Structure files to decreace their context.** TODO, claims, and even handoff files are becoming quite large. We need to create a system to store information in a format that will allow [multiptle] agents to have only what they need in their context window, and ideally allow for parallel work. Perhaps a wiki would help.

---

## Outcome — DONE 2026-08-08

Executed as the four phases of `RESTRUCTURE.md`.

**The diagnosis moved before the work started.** Size was the symptom;
duplication was the disease. `0.716` lived in seven files, `1.006` in six, with
nothing keeping the copies equal — the same failure that produced sixteen
defects in the consistency audit earlier that day. A wiki alone would have
spread those sixteen across 126 files instead of six.

So the answer is **one fact, one home; indexes generated; checks executable** —
not a wiki server. `tools/check.py` is the part that matters: it was verified
adversarially by reintroducing three real defects and confirming it fails.

**Delivered:**

- `CLAUDE.md`, which did not exist. Every session had been rebuilding the
  conventions from `HANDOFF.md`, which is why that file reached 509 lines.
- `tools/check.py` — 10 checks, each mapped to a defect that got through
  without it. Two of its own bugs were caught on first run.
- `tools/reindex.py` — generates every `INDEX.md`, plus `CLAIMS.md` and
  `NOTES.md` as aggregates so ~100 existing references keep resolving.
- `todo/` 29 files · `claims/` 63 files · `notes/` 44 files.

**Measured:** cold-start orientation **16,769 → 4,912 tokens** (3.4×); working
one item ≈ 9,000 all in. `HANDOFF.md` 509 → 113 lines.

**Not done, deliberately:** no ownership or locking protocol. Sessions are
sequential (decided 2026-08-08). The layout is conflict-free anyway, so
parallelism would be additive — see `RESTRUCTURE.md` §11.
