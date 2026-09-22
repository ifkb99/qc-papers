---
id: TX44
field: "Exact arithmetic and GNFS memory scheduling"
status: open
effect: unknown
one_line: "C120 measures a bounded CADO cache-construction saving; production peaks and batch schedules remain open"
source: "CADO-NFS BWC and bucket sources; Bernstein smooth-parts Algorithm 2.1"
claims: [C120]
notes: [FD, FE]
todo: [73]
---
# TX44 — Reduce live intermediate storage in existing factoring stages

## Dictionary and baseline

A GNFS filtered matrix is an exact GF(2) linear map, not a quantum factor-state
coefficient matrix. Keep its layout, vector width and consumer fixed. The
audited builder opportunity concerns temporary allocation lifetimes within one
submatrix, preserving the existing encoded matrix; C120 owns the tested change.

The [CADO BWC README](https://github.com/cado-nfs/cado-nfs/blob/70354d7a8d54e985e46ca0fb6fb64d10716ba8bc/linalg/bwc/README),
`sequential_cache_build` and `build`, already discusses construction peaks and
sequential builds. Its pinned
[bucket implementation](https://github.com/cado-nfs/cado-nfs/blob/70354d7a8d54e985e46ca0fb6fb64d10716ba8bc/linalg/bwc/matmul-bucket.cpp)
already packs selected entries and uses cache-sized processing. The open delta
is temporary storage within one submatrix, compared with sequential building
and finer partitions. C120 owns the completed source audit and bounded actual-
CADO pilot. Production dispatch, concurrency and the full-factorization peak
remain separate applicability questions; upstream documentation is not itself
a measured RSS result.

An alternative uses [Bernstein](https://cr.yp.to/factorization/smoothparts-20040510.pdf),
section 2, on explicit positive norms/cofactors and a prime list. It returns
smooth parts, not complete relation factorizations or algebraic ideal labels.
Its sieving/batch combination is established, and CADO's README.batch already
supports precomputed batch files. Only memory/recompute scheduling is proposed.

## Hypotheses

The cache route matters if steady-state cache/workspace fits and construction
is the obstacle. Bound metadata, buffers, scans and output residency; it cannot
solve an oversized final cache by assertion. Checkpoint deletion primarily
concerns disk and reconstruction, not automatically RAM.

The batch route charges prime-product construction, repeated scans, rejected
survivors, factor recovery and relation metadata. A local saving need not reduce
the full factorization peak.

## Consequence for the goal

C120 supplies a scoped measured construction-memory result, without a new
exponent, smaller final cache or demonstrated whole-factorization peak saving.
The row remains open for production applicability and the separate batch route.
TODO73 owns the next deciding steps; FE records the completed audit and pilot.
