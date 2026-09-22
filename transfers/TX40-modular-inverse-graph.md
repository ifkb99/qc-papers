---
id: TX40
field: "Modular inverse graphs and Schmidt approximation"
status: obstruction
effect: lower-bound
one_line: "Complete modular inverse graphs have flat Schmidt spectra at the whole-factor-register cut; interval-conditioned symbolic access remains open"
source: "Elementary inverse-bijection and Schmidt-projection arguments in C119; TX34 supplies the general Schmidt-tail framework"
claims: [C119]
notes: [FC]
todo: [72]
---
# TX40 — Modular filtering before exact product constraints

## Dictionary

The congruence uv=N mod M, with N a unit and both registers covering all
residues, is a permutation graph on the units. Its coefficient matrix gives
the whole-register Schmidt decomposition directly. C119 owns the exact bound
and the rejection cost after imposing factor intervals.

## Hypotheses

gcd(N,M)=1; complete residue domains; coherent equal amplitudes; the cut
separates the entire factor registers. The power-of-two case assumes k>=1.
CRT relabeling within each whole register preserves that cut's spectrum.
Restricted domains, interleaved physical bits, symbolic inverse maps and the
final equality state are different objects. C113 is not invoked.

## Consequence for the goal

Low-rank truncation of this complete intermediate state cannot provide fixed
fidelity at small register-cut rank. Cheap modular output sampling is still
possible. The unresolved arithmetic step is efficient interval conditioning
with all preparation and rejected candidates charged.
