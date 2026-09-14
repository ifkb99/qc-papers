---
id: 42
state: open
title: "Can cut selection account for total expected sampling work without constructing every row?"
outcome: "TODO41 confirms a favorable bounded case, but C81 selects only grouping-plus-coefficient work"
claims: [C80, C81]
---
# Choose the representation for the requested sample

TODO41/NS implement and test nested-phase sampling. Do not repeat those
integer, amplitude or RNG pilots. The current selector minimizes a precise
construction metric; actual component masses, cancellations, stride-dependent
progression costs and C80 rejection envelopes are not part of its objective.

First derive what an inexpensive pre-construction selector can know from
C81's T/S counts. Distinguish a bound on returned work from an exact optimizer.
Root masses require phase-weighted amplitudes; assuming those values for free
would hide the row construction being avoided. The implemented right-cache
and supplied-table options are relevant inexpensive baselines.

Freeze one bounded discriminator: keep r,b,width, phase positions and the
physical mixer fixed, vary one phase parameter through a cancellation regime,
and compare construction selection with actual costs of candidate cuts on
the same state. The all-cut calculation is a diagnostic reference, not a
free production oracle. Include selection/setup cost, both proposal modes and
stride-dependent queries. A losing selector or a proof that a bound is too
loose is useful; do not promise a general improvement in advance.

Use existing reference engines and a nonvacuous incorrect-selection or
omitted-phase control. Budget numeric payload and named work before runs.
Keep known order/index, phase internals, precision and exhaustion explicit.
Use lower-cost initial audits with distinct file ownership; main reviews and
integrates. No native timing sweep until TODO34 is resolved, and no generic
simulation or novelty claim follows from a fixed-r comparison.
