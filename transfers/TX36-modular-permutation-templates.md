---
id: TX36
field: "Group actions and coherent modular collision sums"
status: imported
effect: classification
one_line: "Multiplicative work permutations share a low-exponent template; full Fourier output still requires coherent modular-power collisions"
source: "Direct clean-state decomposition and collision expansion in C116; toffoli_arith.py u_a"
claims: [C116]
notes: [MC]
todo: []
---
# TX36 — Modular permutation labels and their evaluation cost

## Dictionary

Write the exponent e=s+Lh. A work permutation P_c maps |x> to |cx mod N>
on valid clean inputs. Each high-exponent branch is a transformed copy of
one low-exponent template, with label c=(a^L)^h. These are modular permutation
labels, not the tensor-product Pauli labels of TX29 or the GF(2) frames of TX30.

## Hypotheses

N is odd, a is a unit modulo N, scratch is zero, work starts at |1>,
and the exponent register is uniform; L=2^ell. The normalized orbit state
of C116 stays within x<N at arithmetic block boundaries. A unitary acting
only on low exponent bits commutes with every
work label; one transformed template can therefore be reused. Mixing high
bits forms sums of labels. Full inverse-QFT probabilities involve pairs
with a^e=a^f, equivalently ord_N(a) dividing e-f. The direct derivation and
tiny exact controls are in C116; no imported general lower bound is used.

## Consequence for the goal

Sharing is valid before this collision evaluation. A short symbolic expression
does not establish low evaluated memory or fast sampling. Selected-branch
streaming already avoids duplicate templates. The order-finding barrier TX15
remains; C113's specific MPS cut does not prohibit these nonlinear labels.
No performance advantage was measured. The optional classical order-discovery
benchmark was retired during design because its comparison set and RSS
accounting were incomplete; no measured loss or win is inferred.
