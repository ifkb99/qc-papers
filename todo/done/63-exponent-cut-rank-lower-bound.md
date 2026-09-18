---
id: 63
state: done
outcome: "C106/CR: the cut rank lower-bounds every representation linear across the exponent-prefix cut and is at most |S_k|; for beta=3, alpha=0, where the letters commute, rank = min(k+1, d), and at N=7, a=2 d = 1936 exactly (certified over Z), so C103's widths are the minimal linear bond at every reachable t (<= 1935) and exceed it only beyond; the non-commuting case is not reached (verified at listed sizes only)"
title: "Does the exponent-prefix cut rank of the full-space modexp pullback equal C103's |S_k|, making C103's widths a lower bound for every representation linear across that cut?"
claims: [C48, C103, C104, C106]
---
# Exponent-cut rank as a lower bound

Chosen from the SN slate (2026-09-18), surveyor candidate R1. Let F_k be the ±1 matrix
of f(e, w) = x0 of `ToffoliModExp(N, a, t).build()` on the full dirty space, with rows
indexed by the first k exponent bits and columns by the rest. By C103 Lemma 0 its rows
depend only on the prefix product, so rank F_k ≤ |S_k|. Any representation linear
across that cut (tensor train or MPO over the exponent bits, a weighted automaton, the
class form of note SN, K2's leaf sums) has bond at least rank F_k. By C48 the same rank
governs the Walsh spectrum at that cut.

If rank F_k = |S_k|: C103's widths are optimal among linear representations, and the
ledger gets its first lower bound that is not an output floor. If rank F_k < |S_k|: an
exact TT below the class form exists.

A prior worth stating: a few dozen distinct pseudo-random ±1 rows over thousands of
columns are almost surely independent, so the measured rank will probably equal |S_k|,
and measuring it alone decides little. The content is a proof of independence, or an
honest "not reached". The rank computation is the cheap test that proof predicts. It
needs a must-fail control: the ideal function, whose cut rank is 3 against r = 6 at
N = 7, a = 3 (note TR).
