---
id: 16
state: done
title: Reachable compiled order-finding versus charged classical baselines
outcome: "scratch tables removed; orbit-sized support, no advantage over tested classical baselines; latent-eigenphase sampler isolates order-discovery cost"
claims: [C51, C52]
---
# Reachable compiled order-finding versus charged classical baselines

Completed 2026-09-10 after the two-paper consolidation. Proofs and scope:
C51/C52. Experimental record and reproducibility: §RO. Implementation:
`lab/reachable.py`, the selected-input entry point in `walsh.py`, and
`lab/semiclassical.eigenphase_path`. Existing propagators are reused.

The broader frontier remains TODO 14. Do not interpret removal of a full
scratch allocation as an advantage over orbit enumeration, or a latent-phase
sampler supplied with r as efficient order discovery.
