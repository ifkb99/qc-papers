---
id: 75
state: done
title: "Short straight-line programs for multiples of a hidden factor: what is known, what is conditional, and the strongest obstruction"
outcome: "Mapped; barrier in both directions (TX53). Closed without a round: survey v3 accepted in V31dbc84c4b964cf4; project steps in note NI, section SLP survey"
claims: []
---
# 75: Straight-line programs and factoring (lead (a) of note NI)

**Triage (2026-09-22, METHOD.md barrier check 4):** an explicit family of polylog-length programs producing multiples of p would be nonuniform fast factoring, a TX15 barrier instance. The lead names its structure (compressible integers with hidden divisibility), so the user directed a bounded survey rather than a round. The tractable output is a map of the literature: exact statements, which results are conditional, and the strongest obstruction.

Note NI owns the motivation; TX53 owns the transfer row.

## Question

State precisely, from the primary sources (bodies, not abstracts):

1. Lipton's reduction from short straight-line programs to factoring: the model, the quantifiers (uniform or nonuniform, which integers), and the program length needed.
2. The Shub–Smale τ-conjecture and its known consequences. Bürgisser's results linking τ(n!) to algebraic complexity (VP vs VNP), and in which direction they go.
3. The best known upper bounds on τ(n!) and on programs for products of many small primes, and how Pollard–Strassen and Harvey's deterministic N^{1/5} algorithm fit the model.
4. Any unconditional lower bound, or a barrier (such as a VP ≠ VNP consequence), showing that short programs for these multiples are unlikely or hard to prove impossible.

Deliver an obstruction-first reading: what a success would imply, and whether any restricted program family (for example exponentiation chains, as in p−1 and ECM) has a proved limit. Seconds-scale exact computation on toy integers is allowed only to check a stated lemma. No factoring runs, and nothing at challenge scale.
