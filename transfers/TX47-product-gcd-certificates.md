---
id: TX47
field: "Exact arithmetic certificates for constrained factor search"
status: open
effect: unknown
one_line: "Progression-product gcds give sound rejection or global factors, not conditional counts or a demonstrated fast constructor"
source: "Elementary product/gcd argument in RSA73 structural slate; TX45 factorial baseline"
claims: [C118, C119]
notes: [FD]
todo: [72, 73]
---
# TX47 — A certificate with a limited output

## Dictionary and hypotheses

Let a prefix constrain p to a finite arithmetic progression A inside [2,N-1].
Compute its product modulo N and g=gcd(product,N). If g=1, no member divides N,
so the factor branch is empty. A proper gcd yields a global factor; the bounds
still need checking. For distinct-prime semiprime N, a g=N branch can be split
until isolating a proper gcd, charging recomputation. A nonunit leaf in [2,N-1]
cannot have gcd=N.

Positive gcd does not prove a nonempty factor branch. For N=35 and
A={10,11,12}, gcd(1320,35)=5, but none of those integers divides 35. The product
detects multiples too. Recovering 5 solves the global problem and only then
resolves local bounds; this is not TX46's count oracle.

## Hypotheses

Streaming products already keep small live accumulators while visiting all
members. TX45 supplies stronger factorial/multipoint baselines. New progression
products or sibling-query reuse must improve the complete bill, including
construction and splitting. Low memory alone can leave exponential time.

## Consequence for the goal

This lies outside C118's permitted simple pruning, but escaping that scope
proves no cost improvement. TX41's compact pointwise multiplier does not
supply global products for free. The transfer remains open: a sound interface
is available, but no competitive constructor or additional compression benefit
is established. TODO73 owns its comparison role, not a newly scheduled search.
