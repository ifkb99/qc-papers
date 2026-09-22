---
id: TX14
field: "Noisy circuit simulation by Pauli-path truncation"
status: not-applicable
effect: none
one_line: "Published polynomial-time noisy Pauli-path algorithms need random gates or input ensembles; none applies to fixed-input Shor"
source: "Aharonov-Gao-Landau-Liu-Vazirani arXiv:2211.03999 Def 3-4, Lemma 3, eq 23 pp.9-10; Schuster-Yin-Gao-Yao arXiv:2407.12768v2 Def 1, Thms 1-2 p.3; Gonzalez-Garcia-Cirac-Trivedi arXiv:2407.16068 (read in slate 3a)"
claims: [C5, C14, C16, C99]
notes: [ST]
todo: []
---
# TX14 — Noisy circuit simulation by Pauli-path truncation

## Dictionary

Pauli noise rescales paths without creating new ones: a diagonal pullback keeps the
C8 support with each path damped; exponent (post-QFT) observables are X/Y on the
exponent and damped per location (TX20).

## Hypotheses

Random gates and anticoncentration (AGLLV) fail: the ideal Shor output has 2^t Σp² ≈
2^t/r ≥ 2^n. The input-ensemble average (SYGY) fails for one fixed input, which that
paper calls fundamental because of error correction. Mele et al. and Angrisani et
al. are average-case (abstracts only). GCT covers fixed circuits but needs 2D and p
≳ 0.25, inside TX20's trivial region (TX22).

## Consequence for the goal

None of these removes the exponential for Shor's output. The noise question for this
circuit is answered by TX20 (upper edge) and TX15's low-noise amendment (lower
edge).
