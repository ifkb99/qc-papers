---
id: hk6
state: done
title: Check if julia and dependencies are still needed
claims: []
---
# Check if julia and dependencies are still needed

~~Check if julia and dependencies are still needed~~ — **DONE
2026-08-08. Removed.** No `.py` file imported `juliacall`; the prior-art work
on PauliPropagation.jl was done by *reading* its source, so nothing
reproducible depended on it. Dropped `juliacall`, deleted the orphaned 1.1 GB
depot and `.env`. Venv 1.6 GB → 534 MB, all suites pass, and trap 2 in
`HANDOFF.md` (the `LD_PRELOAD`/`longdouble` segfault) is retired.
