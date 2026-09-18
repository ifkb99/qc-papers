---
id: 58
state: open
title: "Promote the K2 exponent-slice prototype to a registered lab module and experiment so C104's numbers are reproducible from git"
claims: [C104]
---
# Promote the K2 prototype into the tracked tree

Opened 2026-09-17 from the C104 integration review (`V99e34f8a27704055`, RC-7). Every
artifact behind C104 — the K2 prototype, the per-arm readings, the dense pullback, the
survey's query list — lives under the gitignored `out/` tree. So C104 cites no
`experiments/` script, nothing in git reproduces its numbers, and no board record will
flag the claim if those files change. Their sha256 are in note ES.

What this needs: the accepted prototype (k2_proto_v3.py, sha fbca01e4…) as a `lab`
module with the repository's conventions, an `experiments/` script that reproduces the
measured points through `lab.harness` with its own preregistered predictions and
must-fail controls, and the claim updated to cite it. Compare `T163d291c2cdd4cd0`, which
holds the same promotion for the DD-native pilot and is itself unresolved.

Not a correctness gap: the numbers were independently reviewed (`V0e29250446b9478a`).
It is a reproducibility gap.

Prediction to register in that experiment (note SN, deriver candidate C3, proved for the
multiplier count): distinct multipliers T_L(t) = min(t, alpha + nu(beta)), nu = ord_beta(2),
so F2 (N=11, a=2) at t = 7 has 5 tables. Whether distinct multipliers always give distinct
K2 table keys is a property of k2_proto_v3's canonicalisation and must be checked there.
Must-fail: N = 23, a = 5 (beta = 11, nu = 10) gives T_L = t for every t <= 10.
