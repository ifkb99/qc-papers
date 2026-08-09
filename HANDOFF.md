# Handoff — pick up here

Written for a session with no prior context. **Read `CLAUDE.md` first** for the
conventions and traps; this file is only *where things stand*.

**On `main`** as of 2026-08-08 — `todo-12e-gpu-sweeps` was fast-forwarded in
and both now point at the same commit, so the branch can be deleted whenever
wanted. **Nothing is pushed**; `origin` is
[`ifkb99/qc-papers`](https://github.com/ifkb99/qc-papers.git) and is behind.
All nine science suites and the documentation gate pass.

---

## Where things stand

**Both papers are drafted in full and have had a consistency pass and a
structural edit.** `PAPER_A.md` is the exact cost model, `PAPER_B.md` the
2-adic result for modular exponentiation. `ABSTRACT.md` and
`ABSTRACT_SHOR_2ADIC.md` are abstract workshops and prior-art dossiers, **not**
superseded copies to edit in parallel.

Paper A: for any circuit implementing a permutation of the computational basis,
the Pauli support carried by PPS with a computational-basis observable is
*exactly* the Walsh–Hadamard spectrum of the corresponding output-bit Boolean
function — an exact cost model where the field uses empirical extrapolation.
Paper B: applied to modular exponentiation, cost is governed by the 2-adic
structure of the order r = β·2^α — independent of exponent-register width when
β = 1, Θ(2ⁿ) otherwise, with the mechanism proved.

**No live research thread is blocked.** `todo/INDEX.md` is the ranked list;
7 open, 22 done. The best remaining candidate is **14, through the inverse
QFT** — highest risk, highest reach, a Paper C candidate.

---

## → The current task: PAPER REFINEMENT (the user is now editing personally)

Decisions taken 2026-08-08, so they need not be re-asked:

1. **Venue.** Shape both papers **as if for a journal** — focused, every
   paragraph earning its place — but the likely first destination is **arXiv**.
   The user has not published before.
2. **Order.** Consistency pass across both (**done**), then Paper A end to end,
   then Paper B.
3. **The third paper** (compilation-dependence) stays **scattered as caveats
   inside A and B**; revisit after the first two are submitted.

**Done so far:** a consistency pass that fixed sixteen defects, and a
structural pass on Paper A that fixed a duplicate Proposition number, cut §9,
and removed a note-to-self from the shipped text. See commits `184f49c`,
`d68eca3`, `0f2bd6e`.

**Still open on Paper A**, and deliberately left for the user:

- **§10.1 / §10.2** (limitations, retractions) are load-bearing and currently
  in-body. Whether they stay there at journal length is an unmade decision.
- The §3.2 / §4.1 / §6.3 trio explains the same invariance relationship three
  times. Defensible signposting, but it is where a "focused" edit would bite.
- **|B| = 2 is measured, not proved** (Paper A §5). The paper says so. Keep it
  that way unless someone proves it.
- **C33, C34 and C35 are cited in neither paper** — three established claims
  with no home. That is the Paper C material.

---

## Do NOT redo these

Listed here so a fresh session does not burn time; each has a retracted claim
file or a note carrying the full reasoning.

- **The compilation thesis** ("Toffoli vs Fourier compilation determines
  cost"). Dead twice over — a propagator bug plus an ancilla-count confound.
  Both compilations have exact Z-closure and the identity is
  compilation-invariant.
- **"The QFT is the PPS bottleneck."** Dead — bit-reversal bug plus a toy that
  was not a proxy for Shor.
- **The toy "Shor sandwich"** (H → adder → QFT). Not a proxy. Use `modexp.py`
  or `toffoli_arith.py`.
- **Affineness of the identity block** as a route to the C15 proof. Refuted
  (14336/32768 violations). The proof goes via the parity reduction instead.
- **The (z, z⊕e) pairing** as the C15 mechanism. An n_exp=1 artifact.
- **An intermediate 2-adic law in α.** There is none; the β=1/β>1 split is the
  real structure and the apparent intermediate case was sampling aliasing.
- **The msb↔anc CNOT pairing as the cause of the linear structure.** Refuted;
  see §L2 for the actual three-part mechanism.
- **A "cliff" in the hyperplane deficit between n_exp = 2 and 3.** Refuted the
  same day it was conjectured (§DF); D = 1 − 2·density is monotone decreasing
  in every β > 1 instance.
- **The §I period-ord₂(β) oscillation at CIRCUIT level.** Function-level only;
  it does not transfer. That import is what made 12g look worth a sweep.
- **The identity string Z⁰, or a parity constraint, as the peak deficit.** Both
  were the named candidates in TODO 12d and both are wrong; the answer is the
  two |c| = ½ Walsh modes (§PK).
- **Cîrstoiu for TODO 14.** Pulled and read; their group indexes *circuit
  parameters* and their results concern *ensembles*, while Shor's circuit is
  fixed. Useful for framing, not for the technical problem.

---

## Where everything is

| you want | read |
|---|---|
| conventions, traps, discipline | `CLAUDE.md` |
| what is known | `claims/INDEX.md` → `claims/<ID>.md` |
| what is next | `todo/INDEX.md` → `todo/open/*.md` |
| how something was found | `notes/INDEX.md` → `notes/<CODE>-*.md` |
| how the work is conducted | `METHOD.md`, and the `computational-research` skill |
| why the repo is laid out this way | `RESTRUCTURE.md` |
| GPU capacity and measurements | `notes/GPU-*.md` |
| the resolved 0.716 thread | `notes/HIST-*.md` |

`CLAIMS.md`, `NOTES.md` and every `INDEX.md` are **generated** — edit the
individual files, then run `tools/reindex.py`. `tools/check.py` must pass
before any documentation commit.
