---
id: _preamble
skip_index: true
---

# Ledger preamble — conventions carried over from CLAIMS.md

# Claims ledger — single source of truth

**This file is the canonical record of claim statuses for both papers.**
`ABSTRACT.md` (Paper A, the Walsh cost model) and `ABSTRACT_SHOR_2ADIC.md`
(Paper B, Shor / 2-adic) no longer carry ledger tables of their own; they point
here. `NOTES.md` is the **historical working narrative** — it records what was
believed at each point in time and therefore contains **superseded statuses**.
Where `NOTES.md` prose disagrees with a row below, **this file wins**.

Each row must be independently reproducible before the abstract goes anywhere.

**All figures below are post-bugfix** (the θ=π propagator bug, `NOTES.md`
STATUS 2). Anything citing pre-fix numbers is void.

Every claim has its full row in exactly **one** section. Claims that both papers
rest on are cross-referenced by ID rather than duplicated — duplicate rows are
what caused the ledgers to drift once already (a stale C7 row, removed in
`3606e9e`).

The "Where established" column cites only pointers stated in the source
material (`ABSTRACT.md`, `ABSTRACT_SHOR_2ADIC.md`, `todo/INDEX.md`, `HANDOFF.md`,
`README.md`, `NOTES.md` section headers). A blank cell means no pointer was
recorded, not that none exists.

**Paths note (2026-08-08 restructure).** Experiment scripts moved into
`experiments/` with filenames unchanged; a bare `experiment_*.py` reference
below means `experiments/experiment_*.py`. The headline rows are additionally
re-verified by `test_claims.py`, which runs as part of the correctness gate.

---

## Cross-paper notes

**Also relied on by Paper A, full rows in the Paper B section:** C7 (Θ(2ⁿ) for
generic r — the original "results are pre-asymptotic" framing is retracted, see
the Retracted section), C15 and C18 (both moved wholesale to Paper B), C30, C31
and C32 (the density-½ linear structure and the one compilation-dependent cost
effect), **C40 and C41** (the conditional-structure law and its parity
condition, presented in `PAPER_A.md` §6.2 — filed below only because they were
found during the Paper B residue work, and a candidate to move if the two
ledgers are ever split), and **C42** (the inert-support fraction, cited in
`PAPER_A.md` §1.1 to state what the cost model does *not* claim).

**Depends on Paper A, full rows in the Paper A section:** C8 (the Walsh identity,
used as the measurement instrument) and C17 (`perm_pps.py`, which supplies the
peak-memory numbers and makes the n_exp = 8 measurement tractable).

 mapping
Early findings were relabelled when they were promoted to claims. The surviving
ones:
| old | now | note |
|---|---|---|
| F5 | C5 | truncation can violate \|⟨O⟩\| ≤ 1 |
| F9 | C8 | PPS term count = Walsh sparsity |
| F10 | C13 | scope of the identity: Z-type yes, X/Y no |
| F11 | C14 | δ non-monotonicity |
| F13 | C16 | heavy-tailed spectrum explains it |
| F12 | — | kept as F12 in the Paper B section |
