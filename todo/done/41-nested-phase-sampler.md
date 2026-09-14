---
id: 41
state: done
title: "Do nested-phase construction savings survive actual conditional sampling?"
outcome: "Opt-in C81 sampler integrated; tiny laws, actual weighted RNG paths and bounded returned-sample comparison validated, including exhaustion costs"
claims: [C78, C79, C80, C81]
---
# Turn a cheaper row into a cheaper returned sample

**Completed at bounded float scope, 2026-09-12.** C81 owns the implemented
contract and preflight bounds; NS owns authoritative laws, RNG paths, matched
costs, failed verifiers and validation. Existing defaults remain unchanged.
The selected construction saves named sampling work in the frozen fixture;
it is not a general sampling-cost optimizer or a finite-bit certificate.
TODO42 owns the next selector question. The original frozen plan follows.

TODO40 is complete. C81 owns the multi-phase cover and exact finite-truncated
selector proof; NC owns the independently reproduced amplitudes, geometric
counts and selected long-row advantage. Do not repeat those as the next
experiment. The row builder still lives in an experiment. Existing production
classes implement C78/C79 schedules only.

**Historical pause at the user's usage-limit request (2026-09-11).** Only context/design
review and the fresh core gate were completed; no integration code or new
sampler experiment exists. The following implementation sketch is provisional,
not an independently audited or measured result. Resume on user request.

Candidate implementation: subclass LateWorkProgressions, reusing sample,
_sample_row and forced_joint unchanged. Minimize grouping plus coefficient
visits per work row using C81; cache modes must not be advertised as a phase-
cost optimum. Before matrix copies, bound selector setup by H*b^2. A possible
uniform auto-construction bound is the minimum, over candidate cuts, of their
safe grouping-plus-coefficient upper bounds; the selected actual local count
is no larger. A conservative selected component cap follows from G<=C, also
bounded by the full-L support cap. Derive phase/product and optional memo-cache
storage bounds explicitly before coding, and audit row-dependent strides in
the shared loop. The independent audit of this sketch was stopped, not passed.

## First discriminating implementation

Integrate the restricted multiple-diagonal-phase schedule as an opt-in helper,
reusing the actual existing work-first/progression sampling loop and C80's
mass/root_mass modes. Preserve existing defaults and RNG streams. No new
generic vector propagator and no intervening work mixers. Derive the local
column with all d phases while its initial b-label support remains small.
Specify indexed known order, small block, pointwise unit phases and high-history
access explicitly; charge phase internals and integer costs separately.

Use C81's T/S formula to choose a cut per work row without first constructing
every competing row. State the nonnegative cost weights and include selector
setup. The formula is exact for the naive constructor, not optimized phase
caching or cancellation counts. If caching right phases changes the ranking,
derive the amended count before presenting an automatic cost choice. Bound
all candidate and selected allocations/counters before construction, including
per-row stride changes and work retained across rejection retries.

## Freeze before testing

Start with complete tiny literal-column/FFT laws using the already validated
reference engines, then inspect the actual RNG path, retries, gcd lifts,
exhaustion and nonfinite/invalid phase inputs. Include empty, duplicate,
initial/end phases, r=b, non-divisor periods, revived cancellations and exact
zero rows. Preserve NC's near-canceled conditional-law failure: do not declare
tiny positive masses zero or call a float prototype a finite-bit certificate.
When assessing a sampler, distinguish conditional instability from its much
smaller work-probability-weighted joint error.

Only then freeze a bounded returned-sample comparison on the SAME physical
state, with explicit cut alternatives and both proposal modes. Charge column
selection, T/S selection, phase calls/products, construction, Fourier and
progression marginals, every rejected attempt, and failed/exhausted calls.
The empirical branch is whether cheaper construction survives the resulting
rejection envelope and per-proposal query cost. A losing comparison is useful.
Include the strongest inexpensive specialized/cached cut implementation;
full-r sequential or table-based methods remain relevant fixed-r baselines.

Use lower-cost agents for independent initial tests/proof audits while main
owns implementation and integration. Freeze numeric-payload and named-work
caps before execution; stream references instead of scaling a Q-by-r archive.
No native timing sweep until TODO34's host issue is resolved. Working-precision
certification, arbitrary intervening mixers, broad scaling and manuscript /
novelty promotion remain separate. Check primary literature before any claim
of novelty or a competing simulator's limitation.
