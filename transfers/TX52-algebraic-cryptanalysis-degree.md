---
id: TX52
field: "Symmetric cryptography: algebraic cryptanalysis (division property, monomial prediction, cube attacks)"
status: imported
effect: lower-bound
one_line: "An ANF monomial of degree k >= 2 in the pulled-back bit forces PPS term count >= 2^k (lemma accepted); monomial prediction certifies such monomials gate by gate through the X/CNOT/Toffoli u_a blocks, where TX2's Parseval bound is weakest (C26)"
source: "Tsang-Xie-Zhang arXiv:1508.02158, section 1 p.2 (deg2(f) <= log spar(f), uncited there, normalisation unstated) and section 1.1 p.2 (full-degree bound, attainment example), Lemma 5 (restriction); Hu-Sun-Wang-Wang ASIACRYPT 2020, eprint 2020/1048: Definition 1 p.6, Definition 2 p.7, Proposition 1 p.8, COPY/AND/XOR trail rules p.13 (MILP models). All read in the body 2026-09-22"
claims: [C8, C26, C102, C126]
notes: [NI]
todo: [76]
---
# TX52 — Algebraic cryptanalysis: algebraic degree and monomial trails

Coordinator-authored on 2026-09-22 at the user's direction.
- The lemma and the gate check were accepted in review V2b8a1d9005e24dd4
  (submission S1711eb7e5c204612, task Tde7986be1e184ff6). Corrections C1–C7
  are applied here.
- The lemma was applied to the modexp object in C126; the resulting bounds
  and their grades are there.
- C126 certified its monomials by subcube (Möbius) sums, not by counting
  trails. The monomial-prediction half of this row therefore remains unused;
  see Open.

## Dictionary

C8 is the object: the PPS term count for Z_j through a basis permutation U is
the Walsh sparsity of the ±1 form of f = bit j of U. TX2 imported linear
cryptanalysis (the Walsh side). This row imports the *algebraic* side, the
algebraic normal form (ANF) of the same f.

* **ANF monomial ↔ lower bound on PPS support.** If x^S appears in the ANF of
  f with |S| = k ≥ 2, then the ±1 form has ≥ 2^k nonzero Walsh coefficients.
  So PPS term count ≥ 2^{deg₂ f} **when deg₂ f ≥ 2**. At degree 1 the count
  is 1: an affine nonconstant f is a single character.
* **Monomial trail ↔ certificate.** Hu et al. (Proposition 1, p.8): a monomial
  is present iff the number of monomial trails through the chosen
  decomposition is odd. Their models are MILP.

## Hypotheses (checked against the construction)

**Lemma (proved; accepted in V2b8a1d9005e24dd4).** Let f: {0,1}^m → {0,1}
have x^S in its ANF with |S| = k ≥ 2. Then spar(±f) ≥ 2^k.

1. Let g be f with the variables outside S set to 0. Its ANF contains x^S,
   and its top coefficient is Σ_x g(x) mod 2, so wt(g) is odd.
2. For ±g, W(u) = 2^k − 2·wt(g ⊕ u·x).
   - u = 0: W(0) = 2^k − 2·(odd) ≠ 0, because 2^{k−1} is even.
   - u ≠ 0: wt(u·x) = 2^{k−1} is even, so wt(g ⊕ u·x) is odd, and W(u) ≠ 0.
   So ±g has full support.
3. ĝ(u) = Σ_v f̂(u, v). A nonzero ĝ(u) needs a nonzero f̂(u, v), and distinct
   u give distinct pairs.

Step 3 is self-contained; Tsang et al. Lemma 5 (restriction does not increase
sparsity, for real-valued f) agrees with it.
- **Edge cases.** k = 0 is trivially true. k = 1 fails exactly on affine
  nonconstant f.
- **Final vs peak.** The bound is on the final support, so it also bounds the
  propagation peak from below.

**Source statements.** Tsang et al. §1 (p.2) states deg₂ f ≤ log spar(f)
"for every f", without citation. Its normalisation is not stated there, and
reading it as 0/1 is an inference. §1.1 (p.2) states the full-degree case
(spar ≥ 2^n − 1). It also gives the example of a single maxonomial x₁…x_d
whose sparsity can be as small as 2^d, so the lemma's bound is attained.

**Gates and constants (checked by the referee).**
- The logical op set of a u_a block is exactly {x, cnot, toffoli}; every
  Toffoli acts on three distinct wires.
- cswap (circuits.py:96–97) is a correct CNOT–Toffoli–CNOT decomposition.
- The constants a and N enter only through gate placement in `_load`
  (toffoli_arith.py:65–79) and through uncontrolled X gates. Hu's Proposition 1
  holds with X treated as XOR with the constant 1: the referee checked it on
  verbatim toffoli_arith sub-circuits with q ≤ 9, with three mutants detected.
- H appears only in build_shor at the logical level. The rotation-level
  Toffoli decomposition contains H and T; the lemma applies to the logical
  permutation.

**Scope.** It is C8's object: permutation blocks, computational-basis
observables, full space (matching C102's reading of C8).
- **Full degree is impossible on C8's object.** A coordinate of a bijection on
  q bits is balanced, so its weight is even and its top ANF coefficient is 0.
  That gives deg ≤ q − 1. For the single-block u_a x0 bit, C126 sharpens this
  bound and lists the moduli at which the sharper bound is attained.

**Exploration, not a test of a registered prediction:**
- `out/tx52/degree_vs_sparsity_v1.py`, exit 0.
  - 22 functions, 0 violations of the bound at deg ≥ 2. The affine control
    violates the naive k = 1 bound, as required.
  - Modexp rows use the function-level map e ↦ a^e mod N, not C8's object.
    There the bound is exact for N = 21 and 91 only because those functions
    have full degree, which forces it. It is exact without being forced for
    N = 33 and N = 15.
- `out/tx52/parseval_vs_degree_v2.py`, exit 0; supersedes v1.
  - On adder carry-out with c0 = 0 (n = 2..7), the measured sparsity is
    3·2^n − 2. The degree bound is 2^{n+1}, while TX2's Parseval bound
    4^m / max W² is 4 at every width. That is C26's weakness, and the degree
    bound escapes it.
  - Must-fail control: on a bent function (m = 8) Parseval is tight, and the
    claim "degree bound > Parseval bound" is false, as required.
  - This carry-out (c0 = 0, carry kept) is not C102's object. C102 proves
    3·2^(m−1) − 2 for the Cuccaro top sum bit with c0 free and carry-out
    dropped.

**Open:**
- The exact degree of the full-space u_a x0 bit for general n, and a lower
  bound valid for every N (C126, not reached). TODO 76's bounded experiment
  (n ≤ 10) can inform these but cannot settle either one.
- Whether monomial-prediction MILP scales to u_a blocks.
- Whether the linear-rank conjecture (Tsang et al., Conjecture 1) gives more;
  it touches TX30.

## Consequence for the goal

This is an obstruction-side tool. It certifies exponential PPS cost from
algebraic structure the circuit carries, gate by gate, where TX2 is weak. It
removes no exponential. Degree alone does not bound sparsity from above.
Upper-bound uses would need granularity results (Tsang et al. §6), which are
unread here.

## Second section: carries as Witt vectors (note NI, lead (d))

Z₂ = W(F₂), and bit n of binary x + y is the 2-typical Witt sum polynomial
S_n reduced mod 2.
- Checked for n ≤ 3 on all 256 input pairs by `out/ni_witt/witt_v2.py`
  (`uv run --with sympy`; exits 1 on mismatch; supersedes v1): 0 mismatches.
- Must-fail control: carry-free XOR, checked inside the loop, disagrees 272
  times, as required.
- This is exploration.

The ANF degree of the carry into bit n is n + 1.

Open question: does the isobaric grading of S_n (weight p^i on X_i; standard,
not re-read) survive multilinearization over F₂ and constrain which
monomials, and so which Walsh characters, can appear in multiplier bits?
