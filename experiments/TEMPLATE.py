"""<one-line question this experiment answers>

<Context: what is known, which claim/section this extends, why now.>

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  <a consequence DERIVED from the conjectured mechanism, not yet measured>
  P2  <...>
  C1  <the must-fail control: where the effect must vanish, and why>

<If a prediction is open rather than derived, say so and state both branches:
"H-wash: ...; H-persist: ...; no commitment".>

Run:  uv run python -m experiments.<name>     (from research/)

Conventions (see METHOD.md; the harness warns on violations):
  * predictions declared before any measurement;
  * exactly one parameter varies per sweep, stepped by 1 near thresholds;
  * every new circuit variant passes verify_correctness() before being
    measured;
  * a surprising result is treated as a bug until the bug check (precision
    sweep, independent route, smallest failing case) is discharged;
  * on completion: log the outcome in NOTES.md, grade any new claim in
    CLAIMS.md (proved / verified-not-proved / empirical / conjecture), and
    mark this file's header STALE/RETRACTED if later work kills it.
"""
from __future__ import annotations

import numpy as np

from lab import (Experiment, build_modexp, verify_correctness,
                 support, sparsity, density, peak_pps,
                 find_structures, coset_split, quadrant_counts,
                 order, v2_split, bit_table, tile,
                 random_boolean, random_table)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__)

exp.predict("P1", "<prediction>")
exp.must_fail("C1", "<where the effect must vanish>")

exp.section("P1  <what is being measured>")
# me = build_modexp(N=7, a=6, n_exp=2, wraps=())
# assert verify_correctness(me)          # gate every new variant
# s = sparsity(me.build(), me.x[0])      # cached; LAB_NO_CACHE=1 to bypass
# exp.check("P1", s == <derived value>, f"measured {s}")

exp.section("C1  control")
# exp.fail_check("C1", <effect vanished>, "<detail>")

exp.finish()
