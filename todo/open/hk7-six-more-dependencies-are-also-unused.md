---
id: hk7
state: open
title: Six more dependencies are also unused
outcome: decide
claims: []
---
# Six more dependencies are also unused

**Six more dependencies are also unused — decide.** An import scan
found *only* numpy is imported anywhere: `click`, `matplotlib`, `numba`,
`qiskit`, `quimb`, `scipy`, `stim` all have **zero** imports. Not removed
unilaterally because `qiskit` and `stim` were the *source-reading* tools for
the prior-art sweeps (TODO 1) and may be wanted again; the rest look like
genuine leftovers. Recommendation: drop `click`, `matplotlib`, `quimb`,
`scipy`; keep `numba` only if the CPU fallback is ever parallelised; keep or
drop `qiskit`/`stim` on how likely another source-level sweep is.
