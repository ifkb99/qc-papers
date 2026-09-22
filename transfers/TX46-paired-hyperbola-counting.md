---
id: TX46
field: "Divisor sums and exact lattice-point dissection"
status: open
effect: unknown
one_line: "Paired hyperbola dissection supplies a factor-count interface only with a missing constructive cancellation and region bound"
source: "Elementary floor identity; Sladkey arXiv:1206.3369 sections 6-8; RSA73 surveyor slate"
claims: [C118, C119]
notes: [FD]
todo: [72, 73]
---
# TX46 — Count the thin strip rather than compile every candidate

## Dictionary

For positive integer bounds a<=b,

    D_N([a,b]) = sum_(u=a)^b (floor(N/u)-floor((N-1)/u))

counts divisors of N in that interval: division with remainder makes each
summand one exactly when u divides N. For positive factor bounds, clip the
p interval further to [ceil(N/U_q), floor(N/L_q)]. An MSB prefix of p specifies
a contiguous interval; this primitive could therefore give its branch count,
with q=N/p determined after selecting p. It names the missing arithmetic, not
an efficient routine. C119 owns the access consequence; C118 treats a different
specified low-bit constructor.

## Hypotheses

Dissect uv<=N and uv<=N-1 in shared rational/unimodular coordinates, cancel
certified common integer regions before expansion, and recurse on unresolved
parts of the thin strip. The missing lemma bounds construction, residual
regions, coordinate bits and repeated interval queries. Certifying a region
must not silently call the divisor-count oracle itself.

[Sladkey](https://arxiv.org/pdf/1206.3369), sections 6-8, gives a geometric
divisor-summatory baseline and discusses partial sums. Its Theorem 1 counts
operations/regions and recursion depth, not full bit time or allocated bytes.
It does not prove this paired, interval-conditioned speedup. Scans, ordinary
hyperbola grouping and factoring first are additional mandatory baselines.

## Consequence for the goal

No factor state is needed during this arithmetic. This neither constructs the
complete modular-inverse state nor gets zero extraction from compact word
multiplication (TX40-TX42). A small difference requires exact arithmetic or a
certified total absolute error; relative accuracy of large terms is inadequate.

No useful cancellation bound or removal of an exponential is proved. A
restricted-family result or scoped obstruction would also help. TODO73 owns
the first rule to prove/refute. TODO72's constructor question closed at
complexity triage (2026-09-22) as a barrier instance.
