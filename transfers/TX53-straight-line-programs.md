---
id: TX53
field: "Arithmetic complexity: straight-line programs and the tau-conjecture"
status: barrier
effect: barrier
one_line: "Short +,-,x program families that split every k-bit N would be nonuniform factoring; proving none exist would prove a constant-free permanent lower bound"
source: "Buergisser, ECCC TR06-113 (2006), Thm 1.1, Thm 4.1, Cor. 3.9, Lemma 2.12, Cor. 4.2 (ECCC numbering; the journal version, Comput. Complexity 18 (2009) 81-103, numbers differently: journal Cor 3.12 = ECCC Cor 3.9 per Buergisser arXiv:2606.25121 pp.7, 12, 15; the journal version was not read); Dutta, Jindal, Pandey, Sinhababu, CCC 2021 / ECCC TR21-072, section 1 p.2 and Remark 1.4; Cheng, TCS 326 (2004), Thm 1 and section 5; Smale 1998, Problem 4; Markstroem arXiv:1306.3091; Harvey arXiv:2010.05450, Prop. 2.5, Algorithms 4.1-4.2; Lipton ANTS 1994, Shub-Smale 1995, Strassen 1976, BCSS 1998 and Shamir 1979 unread (credited as stated below); project steps in note NI, section SLP survey"
claims: []
notes: [NI]
todo: [75]
---
# TX53 — Arithmetic complexity: straight-line programs and the tau-conjecture

Survey by ni75_surveyor (task T3b864af4007948b2). Reviewed three times; v3
(Sf971590cb261423f) was accepted in V31dbc84c4b964cf4. Its integration
corrections E1–E2 are applied here, and E3–E5 in note NI.

## Dictionary

Note NI lead (a).
- A straight-line program (+, −, × from the constant 1; Smale 1998 Problem 4)
  corresponds to an integer whose residue mod N is tested by a gcd with N.
- A program with one input corresponds to a polynomial evaluated mod N.
  Pollard–Strassen evaluates (x+1)…(x+d) at fixed points (Harvey
  arXiv:2010.05450, Prop. 2.5).
- Harvey's N^{1/5} uses exponents that depend on N (Algorithms 4.1–4.2), which
  places it in TX45.

The family model the project uses, and its findings on it, are in note NI,
section SLP survey.

## Hypotheses

- **Program length for one integer.** τ(m) ≤ 2 log₂ m (Buergisser TR06-113
  §2.2 p.6; Markstroem eq. (1)), and τ(m) ≥ log₂ log₂ m + 1 (Markstroem
  eq. (1)). The second is the only explicit lower bound found.
- **Success direction (secondary sources only):**
  - Dutta et al. §1 p.2: "A candidate hard family is the Pochhammer-Wilkinson
    polynomial f_d := ∏_{i∈[d]}(x+i), for if it turns out to be easy, it
    would imply that integer factorization is also easy [Lip94; Bür09]."
  - Buergisser TR06-113 §1 p.2 credits Lipton 1994 with the implication that
    factoring hard on average gives a weaker form of the τ-conjecture. This is
    a paraphrase, not a quotation.
  - "n! easy ⇒ factoring easy" is credited to Strassen 1976 and BCSS 1998
    p.126 (Buergisser TR06-113 §1 p.2), and to Shamir 1979 (Dutta et al.
    Remark 1.4). All three are unread.
- **Impossibility direction (proved in the source; ECCC numbering):**
  - Thm 4.1(2) with Cor. 3.9: if τ(Per_n) = n^{O(1)}, then
    2^{e(n)} ∏_{k≤n}(X − k) has polylog constant-free, division-free programs,
    with e(n) = (log n)^{O(1)} (proof of Thm 4.1(2), p.15).
  - Thm 1.1: n! hard, or the τ-conjecture, implies τ(Per_n) superpolynomial,
    hence VP⁰ ≠ VNP⁰.
  - Lemma 2.12: τ(Per_n) = n^{O(1)} implies PP ⊆ P/poly.
  - Over C with constants, only under GRH and only for L-complexity
    (Cor. 4.2).
- **Upper bounds:**
  - τ(n!) = O(√n log² n) (Strassen 1976, unread; stated in Cheng 2004 §1 p.2).
  - Some nonzero multiple of n! has length L_n[c] under a conjecture on
    smooth numbers in short intervals (Cheng 2004 Thm 1).
  - Short multiples of n! are not known to factor (Cheng 2004 §5 p.11).
- **Project steps** built on these (the family model, the reading of the n!
  quantifier, the permanent-barrier assembly, a p−1-family bound) are in note
  NI, section SLP survey, with their grades.

## Consequence for the goal

None for the simulation exponentials.
- A success is a TX15 instance: factoring in P/poly.
- A proof that no such family exists would be a constant-free permanent lower
  bound, of VP⁰ ≠ VNP⁰ strength. That step assembles Thm 4.1(2) above and is
  recorded in note NI, section SLP survey.
