---
id: 15b
state: done
title: "2/4: tensor-rank diagnostics and finite-state representations of dense Walsh spectra"
claims: [C8, C45, C48]
outcome: "C48 derived; 72-check rank/spectrum diagnostic and residue realization; full scratch ranks substantially exceed ideal-function ranks"
---
# 2/4 — Tensor rank and finite-state memory

Completed first-pass diagnostics 2026-09-09: C48 and §TR. No general
compressed propagator or minimal-automaton solver is claimed.

After 15a, add reusable real/complex tensor-cut diagnostics, distinct from the
existing GF(2) rank tools. Verify normalized Walsh transforms preserve cut
singular values. Measure final and selected intermediate arithmetic suffixes;
compare the clean-input scalar function with the full scratch-space function.
Use a deterministic random-function negative control and distinguish observed
numerical ranks from exact algebraic ranks. Investigate finite-state/running-
residue realizations, but do not infer a cheap propagation algorithm solely
from a small final rank.

Primary starting points: Kiefer, arXiv:2009.01217; Li, Precup and Rabusseau,
arXiv:2010.10029. Then proceed to 15c.
