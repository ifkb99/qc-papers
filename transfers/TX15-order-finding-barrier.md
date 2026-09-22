---
id: TX15
field: "Computational number theory: order finding and factoring"
status: barrier
effect: barrier
one_line: "Polynomial simulation of Shor's output distribution is a classical order-finding algorithm"
source: "Shor's reduction; C52 (known r gives scalar sampling)"
claims: [C52, C45, C104]
notes: [SN, ST, EF]
todo: [14]
---
# TX15 — Computational number theory: order finding and factoring

## Dictionary

Efficient classical sampler for Shor's output ↔ classical polynomial-time order
finding (hence factoring).

## Hypotheses

Holds for any exact or approximate sampler accurate enough for Shor's post-
processing.

Low-noise amendment (slate 3a S3, derived, not refereed; uses Shor's standard success bound, not re-read): a fault-free run has probability e^(-γL), L = Θ(t n²) noise locations, and r is checkable, so a classical sampler of the noisy output for γ ≤ c ln n / L is classical order finding. The barrier covers noisy sampling below that rate.

Order-class amendment (slate 4 deriver S2, proved there, unrefereed): every exponent-first ordered diagram of the modexp bit has >= r nodes and every x-first one >= the largest prime power dividing r, so a polynomial diagram in those order classes would itself give order finding (via smoothness of r). Interleaved orders are open.

## Consequence for the goal

Not a prohibition: it is the barrier check of METHOD.md. A candidate claiming to
remove the exponential on the Shor family must name the number-theoretic structure
it exploits, or it is a bug or a changed contract.
