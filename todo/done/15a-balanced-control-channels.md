---
id: 15a
state: done
title: "1/4: balanced control functions, reduced channels, and certified tail approximation"
claims: [C45, C46, C47]
outcome: "C47 derived; toy and arithmetic block checks pass; exact biased-tail error and support formula recorded in BC"
---
# 1/4 — Balanced control functions

Completed first pass 2026-09-09: C47 and §BC. Low-level synthesis/performance
of the alternative Boolean controls remains a possible extension, not part
of the block-level result.

Requested 2026-09-09, first in the agreed sequence. For a fresh independent
uniform control window, derive the reduced channel from p=Pr[phi=1], not its
Walsh sparsity. Test parity, majority and OR on a nonvacuous toy and on an
actual arithmetic identity-tail involution. Check full coefficient vectors,
the predicted geometric residual, and a biased-control counterexample.
Record theorem scope separately from circuit implementation and measured cost.

Then proceed to 15b; do not wait for a general optimized compiler.
