---
code: ST
title: "SLATE ROUND 3. The transfer register's open rows: noise is trivial or barrier-bound, gate-uniform holography is #P-hard, and OBDDs of modexp are exponential in every order."
outcome: record
claims: [C52, C86, C102, C103, C104, C106]
todo: [66, 68]
---
# ST — Slate round 3 (2026-09-18): the register's open rows TX14, TX13, TX8

First round under the goal written into METHOD.md on 2026-09-18 ("The goal and the
barrier check"). Three mutually blind fresh contexts, one per open row, each asked for
a premise verdict, a slate and proposed register rows: surveyor 3a on TX14
(`S9be0b75bf1f64536`, v2 of `S3502a0f10bdc43c8`), surveyor 3b on TX13
(`S9c60f8dce5194cf7`), deriver on TX8 (`S0434232124d846c6`). The coordinator sealed a
candidate before dispatch (`Me8495d4c440247ea`, sha256 f780d725…ad50).

## Premise verdicts (all three rows were mis-framed)

* **TX14.** Its obstacle (anticoncentration) is not the binding one; the fixed-input,
  fixed-circuit contract is. Every published polynomial-time noisy Pauli-path result
  read in the body needs random gates or an input ensemble (sources and locations in
  3a §2). Row → `not-applicable`.
* **TX13.** The dictionary holds, but "#P-hard closes every holographic route" is
  wrong: the hardness is worst-case over the gate language and closes only
  gate-uniform methods; C86 is itself a gadget-level escape. Row → `obstruction` with
  that scope.
* **TX8.** Sound question, wrong tool. The deriver's route needs no embedding of
  Bryant's middle bit: bit x0 of z mod N (N odd) is a half-interval test on z·2⁻¹,
  and a direct fooling-set count follows. Row stays `open` pending review (below).

## Results offered (graded by their authors; none reviewed yet)

* **Dephasing lemma (3a S2, proved given Pauli noise with a dephasing part).** Every
  operation touching an exponent qubit is block-diagonal on it (C104 Roles; checked at
  rotation level in run `R7ace5f08b1564841`), so a dephasing event commutes back to the
  start; if every exponent qubit has one, the output is exactly uniform.
  TV ≤ 1 − Π_j(1 − (1−λ)^{g_j}), 13n ≤ g_j ≤ 12n²+n. The coordinator re-read the
  argument and found no gap; it is not refereed and likely folklore (decoherence in
  Shor, Miquel–Paz–Zurek, not searched).
* **Low-noise edge (3a S3).** With a fault-free run of probability ≥ n^−c, a classical
  sampler would give classical order finding: TX15 binds for γ ≤ c ln n / Θ(t n²).
  Only a window of width Θ(t·n) (worst case) remains between the two, where published
  run times are superpolynomial.
* **Gate-language hardness (3b §1.2–1.3).** T4/CCZ lie outside affine, product and
  Hadamard-matchgate classes (Cai–Fu Thm 1.1, 2.31; CFGW Thm 8.1; membership checked
  exactly in run `Rc5b7bf01443648f7`); exact full-space Walsh coefficients of general
  Toffoli circuits are GapP-hard by an explicit dirty-register compiler.
* **Treewidth (3b K3, derived by reading code).** One `cmult_mod` has an m × 5n grid
  minor, so every contraction order of the native network costs 2^Θ(n).
* **Theorem 1 (deriver; proved modulo the Erdős–Turán inequality).** If
  κ(⟨a⟩) ≤ N^(−2η), η ≤ 1/6, every OBDD agreeing with the clean x0 bit on clean inputs
  has ≥ N^η/(C log N) nodes in every variable order; this covers the full-space
  pullback. Restricting e after the order is fixed selects the constant. For prime N,
  r ≥ N^(1/2+2η) suffices (Gauss sums), and typical bases qualify. The coordinator
  re-read the proof (Lemmas R, C, S, G, steps 1–4) and found no gap. Composite N = pq,
  Shor's case, is conjecture.
* **Theorem 2 (deriver; proved, checked n = 3..7).** N = 2^n − 1, a = 2 escapes:
  Σ_{k<t} min(2^k, n) + n + 2 nodes. **The ledger's N = 7 fixtures with a = 2, 4 are in
  this family**, so N = 7 diagram data cannot speak to growth in n. A registered bound
  failed first (run `R43578ea35a36407b`, exit 1, kept) because the script reduced the
  invalid input mod N; v2 (`Rabfa584eb4674c5a`) held.
* **Single block (deriver K3, proved for c = 2).** On the clean code the x0 bit of
  u_a(ctrl, 2) is a comparison, O(n) nodes: C102's single-block growth is all dirty
  inputs, matching note VR.

## Candidates, ranked for the next phase

| rank | candidate | kind | cheapest deciding step |
|---|---|---|---|
| 1 | TX8 Theorem 1, and its composite-N extension | obstruction | referee the proof; surveyor verifies Erdős–Turán, Gauss sums, BGK, and composite-modulus sums |
| 2 | Deriver K4: rank at x-internal cuts ≥ distinct rows | obstruction (TX17) | seconds: exact rank of the arc matrices, n ≤ 11 |
| 3 | Deriver K5: word-level diagrams (*BMD) escape Theorem 1? | the one escape in the diagram family | survey modular-reduction lower bounds for *BMD |
| 4 | 3a S2 dephasing lemma | noise threshold | referee plus decoherence-in-Shor literature check |
| 5 | 3b K3 grid minor | obstruction | seconds: mechanical minor and min-fill bound at N = 7, 11 |
| 6 | 3a S3 fault window | noise | minutes: single-fault enumeration at N = 15, t = 8 |
| 7 | 3b K2 gadget-level holography for `cc_add_mod` | escape (TODO 50) | as TODO 50 |
| 8 | Deriver K2: which ⟨a⟩ escape (N = 2^n + 1, a = 2 next) | escape | derivation |

Killed, with reasons in the submissions: embedding Bryant's MUL (not needed; ⟨a⟩ need
not contain its constants); every published noisy algorithm for fixed-input Shor;
GCT's worst-case theorem (dominated, p ≳ 0.25 is already trivial); gate-uniform
holographic methods; modular counting (3b K6); Montanaro's hitting set (duplicates TX5).

## The sealed candidate

Coordinator K0: under depolarizing noise, weight concentration of the reduced
observable near q/2 would put the truncation threshold at p* ~ 1/(n² t). Surveyor 3a
found a different and sharper mechanism (exponent dephasing, independent of Pauli
weight, per-qubit g_j) and proved it; K0's weight argument is unnecessary and its t
dependence does not appear. The coordinator's sealed candidate lost for the third round
running.

## Reading

The round produced obstructions, which METHOD.md now counts as results. Under growth in
n: no published noisy algorithm helps; gate-uniform holography is dead; every graph
contraction order is 2^Θ(n); and, if Theorem 1 survives review, no ordered bit-level
diagram of the modexp bit is polynomial for typical prime N, even on the clean code.
What escapes: word-level diagrams (K5), gadget/value-level algebra (TX6, K2), and
explicit small-order structure (Theorem 2), which TX15 says must be named.

## Process

* **Output cap.** The deriver's first turn ended on an API error: one response over the
  64k output-token cap, nothing saved after 45 minutes. Resumed with the instruction
  to append files in sections of ≤ 150 lines; it finished in 63 minutes more. Briefs
  for long derivations should say this up front.
* **Lint gating.** 3a listed its submission JSON as evidence and ran lint in the same
  command as `submission.create`; a provenance-only `changes_requested` and v2 (same
  scientific hashes, verified) fixed it. Gate creation on the lint exit.
* **"Changed on disk."** 3b saw its slate.md altered after a write; no formatting hook
  exists and the frozen hash equals the file, so the reviewed bytes are intact.
* **Coordinator ordering slip.** The coordinator edited six register files that were
  snapshotted task inputs before closing the slate tasks, against SWARM.md ("close
  read/proposal tasks before editing their inputs"); the board refused acceptance with
  `stale_inputs`. Recovery: the merged versions were stashed, the archived input bytes
  restored, the reviews and closures recorded (closure notes name the edited rows), and
  the merged versions reapplied; the regenerated index is byte-identical to the stash.
  Detector for next time: merge into a draft outside `transfers/`, close, then apply.
