---
id: TX43
field: "Difference-of-squares factoring and close-factor promises"
status: imported
effect: classification
one_line: "Symmetric factor coordinates give a known efficient close-factor subfamily; no qualifying factor-gap promise is established for RSA challenge inputs"
source: "Elementary Fermat difference-of-squares derivation; phase-0 structural slate Se58d9b2a8c4b45ee, independently reviewed in Vefe60c16abb4432b"
claims: []
notes: [FC]
todo: [72]
---
# TX43 — Symmetry gives a useful subfamily, not a general challenge method

## Dictionary

For odd p<=q use s=(p+q)/2 and d=(q-p)/2. Then s^2-d^2=N and
p=s-d,q=s+d. Testing s from ceil(sqrt(N)) is ordinary difference-of-squares
factoring, with exact square tests and increment z=(s^2-N) by 2s+1.

## Hypotheses

Let Delta=q-p for the first factor pair reached. Since
s-sqrt(N)=Delta^2/(4(s+sqrt(N))), the inclusive number of square tests is at
most 1+Delta^2/(8sqrt(N)). A promise Delta<=N^(1/4)(log N)^c therefore gives
polynomially many tests for fixed c, including polynomial-time integer work
per test. This is a first-discovery bound, not a bound for enumerating all
factor pairs. Square pairs count once; unequal ordered orientations twice.
Balanced bit lengths alone do not establish the required small gap.

## Consequence for the goal

Keep this known restricted-family method as a baseline or bounded prefilter.
Do not select it as a generic RSA challenge route without a justified gap
condition. No new factoring algorithm or factor-gap evidence was obtained.
