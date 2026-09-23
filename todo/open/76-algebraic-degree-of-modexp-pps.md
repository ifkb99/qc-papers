---
id: 76
state: open
title: "Algebraic degree of the full-space u_a pulled-back bit: a certified PPS lower bound via monomial trails (TX52)"
outcome: ""
claims: [C8, C102, C126]
---
# 76: Degree lower bounds for the modexp PPS object

**Triage (2026-09-22, METHOD.md barrier check 4):** a lower bound on PPS cost computes nothing hard, so it passes.

TX52 owns the dictionary and lemma (a degree-k ANF monomial gives PPS term count ≥ 2^k); note NI owns the motivation.

## Question

For C102's full-space object (x0 bit of u_a(ctrl, 2) as built by `toffoli_arith.py`), derive from the construction a monomial of provably odd monomial-trail count, and so a degree lower bound, as a function of n. Compare the resulting 2^deg bound with C102's measured Walsh support at n = 4..8, before any new measurement.

## Steps

1. Referee the TX52 lemma, including the ±1 vs 0/1 normalization and the k = 1 control.
2. Derivation: read the actual gate list, including how constants enter. Identify a candidate top monomial and argue that its trail count is odd, for example through a unique trail as in the adder carry.
3. Only after an accepted derivation: a bounded check. Either an exact ANF at small n, or a monomial-prediction MILP/SAT at larger n, with a must-fail control. The candidate is a monomial the derivation says is absent.

## Also open

Whether the Witt-vector isobaric grading (TX52, second section) constrains multiplier monomials. This is reading and derivation only.

## Status (2026-09-22; superseded by the next section)

Steps 1 and 2 are done. The TX52 lemma was accepted in V2b8a1d9005e24dd4. The derivation is claim C126, which states its grades and what was not reached. What remains is step 3: the bounded experiment in section 9 of the accepted plan, derivation v2 (Scb013aeebb964e0f). It must apply review V16acbf38c34b4414's execution-time corrections C1–C6. It answers Q1 (does deg f reach 3n + 2?) and, conditionally, Q2 (n = 8..10 by the interval route).

## Status (2026-09-22, after the experiment)

Step 3 is done. The experiment (T066df2c305284ceb, accepted in V69a199ae5c5a4984) is recorded in C126's attainment section. What remains:

- the top coefficient, or the exact degree, for general n (C126, not reached); its data are the experiment tables (`out/agent-board/workers/Adbf05967c8a5413a/tables_v1.txt`), whose moduli C126 lists;
- the Witt-vector item above.

