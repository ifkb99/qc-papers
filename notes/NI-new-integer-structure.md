---
code: NI
date: 2026-09-22
title: "What a new structure in the integers would have to do: Joux's anatomy, the transplant failures, and four leads"
outcome: record
claims: []
todo: [75]
---
# NI — The missing structure in the integers (2026-09-22)

## Why

The user asked whether simulating Shor or Grover could factor RSA challenge
numbers. Complexity triage (METHOD.md, check 4) marks that as a TX15 barrier
instance: it computes factoring. The user then asked a different question. Joux
and others found new structure for small-characteristic discrete logs, so what
would a comparable structure in the integers have to do? This note records that
discussion and the leads it produced. It is conceptual and coordinator-authored.
Nothing here is refereed, and no claim rests on it. Literature statements marked
**(unverified)** were written from memory and must be read in their sources
before anything cites them.

## Anatomy of the small-characteristic collapse

The route ran from Joux 2013 (L[1/4]) to Barbulescu–Gaudry–Joux–Thomé 2014
(heuristic quasi-polynomial), then to Kleinjung–Wesolowski, who proved
quasi-polynomial time for fixed characteristic **(unverified)**. Four ingredients:

1. **Representation freedom.** All fields of size q^k are isomorphic, so one
   can choose a defining polynomial that makes X^q cheap (X^q = h0/h1 with low
   degree).
2. **An identity that splits completely.** X^q − X = ∏_{a∈F_q}(X − a) is
   smooth by construction, and PGL₂(F_q) turns one relation into many.
3. **Smoothness is constructed, not waited for.** Every index calculus pays
   ρ(u) ≈ u^(−u) for luck. Here one side of each relation is free.
4. **Self-similar descent.** Degree D goes to about D/2 with polynomial
   branching, which gives n^{O(log n)}.

## Two templates for fast factoring

* **A. Annihilate a group** (p−1, ECM, Shor): a group whose order mod p
  divides a cheaply computable exponent M.
* **B. Collect relations** (quadratic sieve, NFS, Schnorr's lattice attempt):
  an auxiliary ring R with a map φ: R → Z/N, elements smooth on both sides,
  then linear algebra to a congruence of squares. Testing numbers of size L[β]
  against a factor base of size L[γ] gives L[α] with α ≈ β/2. QS tests numbers
  of size L[1] and runs in L[1/2]; NFS tests L[2/3] and runs in L[1/3]. L[1/4]
  needs smooth parts of size L[1/2]; quasi-polynomial time needs relations
  with no size bound.

## Why each direct transplant fails

| Joux's ingredient | What happens in Z/N | Requirement it leaves |
|---|---|---|
| Representation freedom | Auxiliary rings reach Z/N only through a known root mod N, built from N's digits (polynomial selection). A volume bound gives coefficients ≈ N^{1/(d+1)}. | Relations involving N that escape the volume bound, without being SNFS in disguise |
| Computable symmetries | Other roots mod N are unknown; square roots mod N are equivalent to factoring (Rabin) | Every symmetry used must be computable without p |
| A small constant field | Z/N ≅ F_p × F_q with both large. Prime-field DL, with p known, still sits at L[1/3]. | Something playing a small base field for Z |
| Frobenius | x ↦ x^p is hidden. AKS: x ↦ x^N acts as a Frobenius only when N is prime. | Partial Frobenius access without p |
| Forced order structure (template A) | Lang–Weil gives order ≈ p^{dim}, random-looking; forcing structure (CM) needs p. | A group family with smooth order mod p, built from N alone |

The pattern: every structure that obviously helps is a function of p, not N.

## Test for any proposed structure

1. Is it computable from N alone?
2. Does it produce smoothness or relations **by construction**, and by which identity?
3. Does it involve N, or only an auxiliary ring?
4. Why does it not also break prime-field discrete logs, which remain at L[1/3]?

An idea with no answer to question 4 probably works nowhere. The operational
form of this test is the `integer-structure` skill.

## Four leads

Each carries METHOD.md's triage line. All four concern factoring, so each is a
TX15 barrier instance **unless** it names its structure; each lead below names
one, which is why they are recorded rather than dismissed. None is scheduled
as compute.

* **(a) Short straight-line programs (TX53, now `barrier`; TODO 75 closed).**
  Find a compressible integer that is secretly divisible by p: a short +,−,×
  program from 1 whose value mod N shares a factor with N. The from-memory
  citations first written here were partly wrong: Lipton 1994 concerns
  polynomials with many roots, not n!. The corrected record is TX53 plus the
  section "SLP survey" below.
  *If this succeeds, it computes:* a factor of N nonuniformly, which is a TX15
  instance. Proving that it fails would be a constant-free permanent lower
  bound. It is barrier-bound in both directions.
* **(b) Free relations mod N (TX54).** A polynomial-size set b_i with a supply
  of multiplicative relations mod N that are not consequences of relations in
  Z. With one element this is order finding. Schnorr's short-lattice-vector
  relations were not smooth often enough.
  *If this succeeds, it computes:* factoring. Kept only as a statement of
  template B's requirement.
* **(c) A constant field for Z (TX55).** The F₁ program, Borger's λ-rings
  (commuting Frobenius lifts), and Deninger's proposed cohomology
  **(unverified in detail)**. As far as the coordinator knows, none has
  computational content.
  *If this succeeds, it computes:* nothing directly. It is the place the
  missing Frobenius is sought, a reading lead.
* **(d) Carries as Witt vectors (TX52, second section).** F₂[x] and Z differ
  only by carries. Factoring in F₂[x] is polynomial (Berlekamp: squaring is
  F₂-linear), while in Z squaring carries. Z₂ = W(F₂), and the binary digits of
  x + y are the 2-typical Witt sum polynomials S_n mod 2.
  - Verified for bits 0–3 on all 256 input pairs: `out/ni_witt/witt_v2.py`
    (`uv run --with sympy`; v1 superseded after review V2b8a1d9005e24dd4), exit 0
  - Must-fail control: carry-free XOR disagrees, as it must.
  - This is exploration, not a claim.
  *If this succeeds, it computes:* degree and sparsity bounds on carry
  polynomials, a property of the simulator's cost. It passes triage.

## Calibration

The right shape would make smoothness happen by construction for relations
involving N, the role Frobenius plays for F_q. No candidate is known, and the
field has looked since 1990. Leads (a)–(c) are reading directions, not rounds.
Lead (d) connects the structure question to the carry costs this repository
already measures.

## SLP survey (lead (a), TODO 75)

Author: ni75_surveyor (task T3b864af4007948b2), v3, accepted in review
V31dbc84c4b964cf4. The integration corrections E3–E5 are applied here.
Runs cited: R4580260e0c024ce9 and Rc5a434b1aa804b54.

This section records the steps this project took. TX53 holds the imported
statements the row needs, and they are cited from here, not restated. The one
exception is D2's number-theoretic input, which only D2 uses; it is recorded
below with its (secondhand) source.

Grades:
- **derived**: argument written out here, not refereed;
- **reading**: an interpretation of sources, not a theorem;
- **proved**: a published theorem, held in TX53.

Buergisser theorem numbers use ECCC TR06-113 numbering, as in TX53.

**Model (project definition).**
- For a single N, a short program for a factor always exists (TX53, first hypothesis). A single
  integer M cannot split every semiprime either: of any three primes, two lie on the same side of
  "divides M" (derived, elementary). So the object with content is a family.
- For each bit length k, the family is a list of programs from the constant 1, each optionally
  with one input. A procedure running in time polynomial in k and the total program length L(k) evaluates them mod N, may multiply outputs mod N, and
  takes gcds with N.
- The procedure chooses programs and evaluation points using only k and the outcomes of earlier
  gcds.
- A family of poly(k) total length that splits every k-bit N puts factoring in P/poly.
- Pollard-Strassen is a family in this model, but of length exponential in k (its source is in
  TX53). Methods with N-dependent exponents are outside the model (TX53, TX45).

**C-d, the quantifier in "n! short => factoring" (reading).** TX53 records who is credited with the
statement; those primary sources are unread.
- Three bodies were read that sketch the reduction as a binary search for the least m with
  gcd(m! mod N, N) > 1:
  - Cheng TCS 326 (2004), section 1 pp.1-2, which states it uniformly;
  - Dutta et al., proof of Thm 6.2 and Algorithm 1 (ECCC TR21-072 pp.18-19, 39);
  - the Lipton-Regan blog post of 2009-02-23, whose model also allows division.
- The search needs m! at N-dependent points. For balanced N = pq with p < q close together, some
  m in [p, q) is needed, so a poly(k)-size list of per-integer programs cannot serve every N.
- The statement is safe in two forms:
  1. uniformly (programs for m! constructible in poly(log m) time); and
  2. nonuniformly, via one-input programs for prod_{i<=2^j}(x+i), j <= k/2, evaluated at
     gcd-chosen points.
- The unread primary sources may avoid the gap. Short multiples of n! are a separate matter; see
  TX53.

**D3, impossibility is as hard as a permanent lower bound (proved theorem from TX53; the added
step is derived).**
- Imported: TX53's Buergisser Thm 4.1(2) with Cor. 3.9.
- Added step: evaluating the program that TX53's Thm 4.1(2) supplies at X = a + n + 1 mod odd N
  gives a unit times prod_{t=a+1}^{a+n} t (the power of 2 is a unit because N is odd). With
  n = 2^j, j <= k/2, a gcd-driven binary search over products of such blocks finds the least
  prime factor of N. That is a poly(k)
  family in the model above.
- Hence a proof that no poly(k) family exists would give tau(Per_n) superpolynomial, and so
  VP^0 != VNP^0. The Boolean route through TX53's Lemma 2.12 does not give this, because the
  model forbids N-dependent choices.

**D2, p-1-type families (derived; unrefereed; the input is secondhand).**
- *Input* (held only here, since only D2 uses it): Goldfeld 1969 (some theta > 1/2) and
  Baker-Harman (theta = 0.677), as stated in Ding-Wang arXiv:2510.04026v2, section 1, eq. (1.2):
  #{p <= x : P+(p-1) >= p^theta} >= delta pi(x). Under Elliott-Halberstam every theta < 1 works
  (eq. (1.3)). The primary bodies were not read.
- *Statement.* Take integers M_i = a_i^{E_i} - 1 with |a_i| >= 2 and total base size
  sum_i log2(|a_i|+1) <= x^{2 theta - 1 - eps}. This caps the number of programs below
  x^{2 theta - 1 - eps}: about x^{0.354} = X^{0.177} at theta = 0.677, approaching x = X^{1/2}
  as theta -> 1 under EH. Suppose some M_i is divisible by p or p' for every pair of distinct
  primes p, p' <= x = X^{1/2}. Then sum_i log2 E_i >= x^{theta - o(1)}.
- *Measures.*
  - M1 := sum_i tau(M_i), one program per output (no input, no sharing). Here the bound is
    X^{theta/2 - o(1)}, which is X^{0.338} (X^{1/2-o(1)} under EH). For gcd-adaptive procedures it holds along the
    worst-case path, meaning the path taken for some semiprime N <= X; it does not hold along every path.
  - M2, one shared multi-output program: X^{theta/4 - o(1)}, about X^{0.169}. A heuristic
    Fibonacci-exponent-chain estimate (about 0.21 L^2 prime bits from L multiplications; not
    run) suggests the M2 counting cannot do better.
  - One-input programs reused at several points (for example x^E - 1 at many bases) are not
    covered by D2.
- *No comparison with Pollard-Strassen is established* in either measure.
- *Mechanism.*
  - p | a^E - 1 iff ord_p(a) | E.
  - Primes with ord_p(a) <= x^{1-theta} number at most x^{2-2theta} log2(|a|+1).
  - For the other primes, P+(p-1) >= p^theta divides every ord_p(a_i), so it must divide
    prod E_i for all but one of them.
  - Each such prime factor serves at most x^{1-theta+o(1)} primes p.
- *Scope limits.* The output form a^E - 1 and the base-size cap are essential; a base 1 + M'
  makes a^E - 1 a multiple of any M'. The ECM analogue was not checked. The weakness of p-1 when
  p-1 has a large prime factor is folklore; no source read states this bound, and no novelty is
  claimed.

**Consequence.** These steps leave TX53 a two-sided barrier (see its Consequence). They carry
nothing for the simulation exponentials.
