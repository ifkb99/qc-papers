---
id: 56
state: done
outcome: "Closed without a round (note SN): both slate-2 workers generalised the class form to every beta via the prefix automaton and found its value conditional -- no ledger consumer of coefficient queries, and C103's exponent-first ROBDD is a smaller exact O3 object at the one benchmark point (8,517 nodes vs 58,279 class terms); minimality among linear representations moves to TODO 63"
title: "Can the exact final modexp operator be stored in exponent-class form for beta in {1, 3}, below the ~2^(q-1) output floor?"
claims: [C24, C101, C103]
---
# Exponent-class form of the output (queued behind TODO 55)

Chosen from the SL slate (2026-09-16) under output contract **O3**: an exact
compressed representation of the final operator that answers coefficient
queries, instead of an enumerated output.

For β = 3, α = 0 every block letter is M^{±1} with M = M_{c_0} (C101), so the
product depends on the exponent only through the signed count
s(e) = Σ (−1)^k e_k, and the Walsh transform factors through Krawtchouk
polynomials over the even and odd exponent positions (deriver slate lemma L3,
unreviewed). Storing one work spectrum per reachable class gives a proxy of
58,279 class terms against 265,232,710 Walsh terms at t = 16, N = 7, a = 2
[consistency check against C103's recorded value, not a blind test].

Open before any derivation task: the class-form PEAK during propagation (not
only the final form), β ≥ 5 where letters do not commute, and the cost of
answering queries. Queued: start after TODO 55's design review.

## Closed (2026-09-18)

See note SN, "TODO 56, closed". Reopen only with a named consumer of coefficient
queries, stated against the (S_store, T_build, T_query) contract there.
