---
id: hk1
state: done
title: Yao.jl backend
claims: []
---
# Yao.jl backend

Yao.jl backend — closed as far as possible without installing it. Its
docs list Toffoli under **"Clifford Gates: Two-qubit gates"**, wrong on both
counts (Toffoli is neither two-qubit nor Clifford), so the doc is unreliable;
and PauliPropagation.jl, which it most plausibly wraps, has **zero** mentions
of Toffoli/CCX/CCZ and no permutation gate type. No library documents or
exploits diagonal closure. A definitive check would need Yao installed.
