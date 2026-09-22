---
code: TL
title: "PHASE 1 AFTER ST. Ordered representations of the modexp output bit are exponential in n in every variable order: OBDD (refereed), word-level diagrams and tensor-train bond (via a new rank identity), for typical prime N; Shor's N = pq for typical bases."
outcome: record
claims: [C102, C103, C106, C48]
todo: [68]
---
# TL — Phase 1 after slate round ST (2026-09-18): lower bounds for ordered representations

The user chose items 1–3 of note ST's ranking. Four fresh, mutually blind contexts:
a referee of TX8's derivation (task `T1e74d6a66ec34037`, review filed as submission
`S86f11aef594648bd`, since the author's task was closed and its register inputs renamed);
surveyor 4a on the exponential-sum inputs (`Sf8eb4afd15a84c9e`, v2 of `S67f4791a54ad4dac`);
surveyor 4b on word-level diagrams (`S20cdcbcfcbcf4547`, v2 of `S6f97e76d496b417f`);
deriver 4k on rank at x-internal cuts (`S11b9708e34d14abe`). No claim is promoted here.

## Object

f(e, x) = bit x0 of `ToffoliModExp(N, a, n_exp=t).build()`: on clean inputs,
bit_0(a^e (x XOR 1) mod N); the full-space pullback restricts to it. Every result below
is about representing this Boolean function, not about sampling Shor's output.

## Results and their status

* **Theorem 1 (TX8): refereed, accept with corrections R1–R8.** Every OBDD of f, in every
  variable order, has ≥ N^η/(C log N) nodes when κ(⟨a⟩) ≤ N^(−2η), η ≤ 1/6. The referee
  re-derived every step and found nothing refuted. Required before integration: R1 add
  the proof that the circuit's clean bit equals the formula for all N (the referee
  proved it; the author verified 7 instances); R2–R3 constants (N0 justification; the
  1/K term is 2N^(−η)); R4 Theorem 2's "exactly" needs 2^t ≥ n (n = 7, t = 2: 9, not
  12); R5–R6 scope wording (v1's mod-N statement refuted; clean code vs full space for
  the N = 7 fixtures); R7 the single-block "O(n)" holds in MSB-first order only (7 → 37
  over n = 4..14 in evens-then-odds); R8 K3's proposed step is not to be scheduled as
  written.
* **External inputs (surveyor 4a): verified.** Erdős–Turán with C0 = 6
  (Kuipers–Niederreiter Thm 2.5, p. 112, read from the page image). The √p bound
  follows from Parseval, |Σ_{c∈H} e(cu/p)| ≤ √(p − r), so Gauss sums are not needed.
  Bourgain–Konyagin 2003 Thm 2.1 read in the body; BGK 2006 only as restated.
* **Shor's moduli (surveyor 4a, proved there, unrefereed).** For N = pq,
  κ ≤ max(√N/r, √p/r_p, √q/r_q), and κ ≥ (p−1)^(−1/2) unless a is a primitive root mod p.
  So for balanced N, η ≤ 1/8 unless a is a primitive root modulo both factors; typical a
  get η = 1/8 − θ/2 when gcd(p−1, q−1) ≤ N^(1/4−θ); primitive roots of both get η → 1/6 only
  when gcd(p−1, q−1) ≤ N^(1/6) (corrected 2026-09-18; C109 owns the statement).
  The hypothesis fails exactly when some order r_p or r_q is small.
* **Rank identity (deriver 4k, proved, unrefereed).** For any 0/1 matrix whose rows are
  cyclic translates of one arc sampled on point sets S, T ⊂ Z_N:
  rank_Q = E − β₁(G) + ε, with G the gap graph (one edge per distinct nonconstant row).
  Hence (D − 2)/4 ≤ rank ≤ D, and Theorem 1's fooling-set rows have rank ≥ ⌈K'/2⌉.
  1_I has no DFT zero for odd N, so rank loss comes from sampling (corrected: a general
  arc of length ℓ has one iff gcd(ℓ, N) > 1; C108).
  Checked at 148,678 instances (corrected 2026-09-18: the evidence brackets the rank where a
  full row and a winding cycle coexist; C108 owns the statement).
* **Consequence (4k, modulo Theorem 1).** In every order, the cut after half of x's low
  bits has rank ≥ N^η/(870 ln N) (corrected; C108): every tensor train/MPS, MPO or weighted automaton
  reading that order has at least that bond there. This is TODO 68's missing lower
  bound, in every order, for prime N.
* **Word-level diagrams (surveyor 4b).** Scholl–Becker–Weis ICCAD'98 Thm 1, Lemma 1,
  Thm 2 (pp. 674–675, read in the body): every ordered MTBDD, EVBDD, *BMD, HDD, K*BMD or
  *PHDD maps into a WLCD with no more nodes, and WLCD size ≥ Q-rank at each prefix cut.
  With 4k's consequence, every ordered word-level diagram of f is exponential in n in
  every order (prime N, κ regime). Independently, exponent-first or exponent-last
  orders need ≥ 1 + r/2 nodes (derived modulo two recalled L-function facts; measured
  for primes 11..113). The known division lower bound does not apply (variable divisor).
  EVBDD/MTBDD also inherit Theorem 1 via Becker–Drechsler–Enders ASP-DAC'97 p. 463.

## Reading

For typical prime moduli, and for typical bases of Shor's N = pq, **no ordered linear
representation of the modexp output bit is polynomial in n in any variable order**
(corrected 2026-09-18: Pauli-LIMDDs and GF(2)-frame diagrams are outside this; see C107, C108):
OBDD (refereed), and, pending review of 4k's identity, every ordered word-level
diagram and every tensor-network bond at some cut. This covers the diagram and
tensor-network directions the register listed as escapes (TX7, TX9, TX17, TX26). It
respects TX15 (a lower bound). What it does not cover: unordered or free diagrams
(FBDD, free BMDs, SDDs), representations in a GF(2)-transformed frame, non-linear
representations, sampling contracts, and small-order families (Theorem 2's
N = 2^n − 1, a = 2, which contains the ledger's N = 7 fixtures).

## Process

* **Lint false positive, twice.** 4a and 4b printed "FAILED as intended/required" for a
  control that behaved correctly; submission lint reads FAILED in an exit-0 log as a
  contradiction (EXIT-CONTRADICTS-LOG). Both fixed by provenance-only v2 runs. The
  working-file lint does not compare logs with exits, so only the submission lint
  catches it. Protocol fix: report a correct must-fail outcome without that token (as
  `lab.harness.fail_check` does).
* **Seal location.** 4b listed the session scratchpad and saw the sealed-candidate file
  name (not its content; it was already revealed). The scratchpad is shared with
  workers; seals belong somewhere workers cannot list.
* **Closed-task review.** A review of a closed submission whose inputs were renamed goes
  through a separate referee task, as here.
