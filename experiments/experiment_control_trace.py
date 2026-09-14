"""Exact contraction of finished |+> exponent controls, 2026-09-09.

DERIVATION BEFORE MEASUREMENT.
After the first forward occurrence of q has been back-propagated, all remaining
gates avoid q. Contracting q with its independent |+> input therefore commutes
with the remaining propagation. Since this representation is diagonal, this
simply drops Z-terms carrying q. This argument holds for BOTH order branches.
For consecutive independently controlled identical involutions, the reduced
map E=(Id+Ad_V)/2 is idempotent, so all but one identity-tail block can be omitted.
That second simplification relies on beta=1 in this arithmetic construction.

The observable is a work-register Z BEFORE the inverse QFT. This is a reduced
observable computation, not a simulation of Shor output sampling. Peaks count
retained terms before contractions, not bytes or simultaneous dictionaries.

PREDICTIONS.
  P1  Every final reduced coefficient agrees with terminal projection of the
      uncontracted operator, and small instances agree independently with Walsh.
  P2  With m work qubits and one exponent control per contiguous block, peak
      support is <= 2^(m+1), for beta=1 AND beta>1, independently of width.
  P3  Keeping alpha+1 blocks preserves the reduced operator when beta=1.
  P4  All reduced expectations agree with direct classical enumeration of the
      legitimate modular-exponentiation input distribution.
  H1  Measured peak and runtime improve; sizes are measured, not predicted.
  C1  Contracting a reused control before its final use must give a wrong result.
  C2  Applying the tail-collapse prescription to beta>1 must fail at some width.

Run (Python 3.12 avoids the documented 3.14 long-propagation crash):
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      -m experiments.experiment_control_trace --max-width 6
Results: out/control_trace.json (generated, ignored by git).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from circuits import Circuit
from lab import Experiment
from perm_pps import propagate_perm
from toffoli_arith import ToffoliModExp
import walsh


def max_error(left, right):
    return max((abs(left.get(z, 0.) - right.get(z, 0.))
                for z in left.keys() | right.keys()), default=0.)


def prefix(me, block_count):
    qc = Circuit(me.n_qubits)
    qc.x(me.x[0])
    for i, eq in enumerate(me.exp[:block_count]):
        qc.extend(me.u_a(eq, pow(me.a, 1 << i, me.N)))
    return qc


def measured(qc, target, trace=()):
    start = time.perf_counter()
    result = propagate_perm(qc, 1 << target, trace_plus=trace, max_terms=600_000)
    return result, time.perf_counter() - start


def extended_checks(width, output):
    """Beyond the baseline sweep: bound, actual input, and tail-control checks.

    No full uncontracted operator or Walsh vector is constructed at this width.
    """
    exp = Experiment("control_trace_extended", doc=__doc__)
    exp.predict("P2", "peak <= 2^(work+1) in both order branches")
    exp.predict("P3", "beta=1 reduced operator equals its compressed prefix")
    exp.predict("P4", "expectation matches all legitimate classical inputs")
    exp.must_fail("C2", "beta>1 reduced operator must differ from compressed prefix")
    rows = []
    for base in (6, 3):
        me = ToffoliModExp(N=7, a=base, n_exp=width)
        result, seconds = measured(me.build(), me.x[0], me.exp)
        collapsed, collapsed_seconds = measured(prefix(me, 2), me.x[0], me.exp)
        assert not result.hit_cap and not collapsed.hit_cap
        error = max_error(result.final_terms, collapsed.final_terms)
        truth = sum(1 - 2 * (pow(base, e, 7) & 1)
                    for e in range(1 << width)) / (1 << width)
        exp.check("P2", result.n_max <= 1 << (me.exp[0] + 1),
                  f"a={base}, t={width}: peak={result.n_max}")
        exp.check("P4", abs(result.expectation - truth) < 1e-10,
                  f"a={base}: expectation={result.expectation}, true={truth}")
        if base == 6:
            exp.check("P3", error < 1e-11, f"every compressed coefficient: error={error}")
        else:
            exp.fail_check("C2", error > 1e-6, f"invalid compression error={error}")
        rows.append(dict(N=7, a=base, t=width, qubits=me.n_qubits,
                         trace_peak=result.n_max, trace_final=len(result.final_terms),
                         trace_seconds=seconds, expectation=result.expectation,
                         classical_expectation=truth, collapsed_error=error,
                         collapsed_seconds=collapsed_seconds))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(rows, indent=2) + "\n")
    exp.finish()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-width", type=int, default=6)
    parser.add_argument("--baseline-width", type=int, default=6)
    parser.add_argument("--output", default=None)
    parser.add_argument("--extended-width", type=int, default=0,
                        help="run only extended checks without the full baseline (4..20)")
    args = parser.parse_args()
    if args.extended_width:
        if not 4 <= args.extended_width <= 20:
            parser.error("--extended-width must be between 4 and 20")
        extended_checks(args.extended_width, args.output or "out/control_trace_extended.json")
        return
    args.output = args.output or "out/control_trace.json"
    if args.max_width < 4:
        parser.error("--max-width must be >=4 to test beyond the parity onset")

    exp = Experiment("control_trace", doc=__doc__)
    exp.predict("P1", "reduced coefficient vectors match terminal projection and independent Walsh")
    exp.predict("P2", "peak <= 2^(work+1) for either order branch")
    exp.predict("P3", "one identity-tail block suffices after |+> contraction")
    exp.predict("P4", "expectations match the actual uniform-exponent input")
    exp.predict("H1", "measure improvement in at least one nontrivial case")
    exp.must_fail("C1", "premature contraction of a reused control must be wrong")
    exp.must_fail("C2", "the beta>1 branch must defeat identity-tail compression")

    reuse = Circuit(2).cnot(0, 1).cnot(0, 1)
    correct = propagate_perm(reuse, 2, trace_plus=[0])
    half = propagate_perm(Circuit(2).cnot(0, 1), 2)
    wrong = {z: c for z, c in half.final_terms.items() if not (z & 1)}
    exp.fail_check("C1", correct.final_terms == {2: 1.} and wrong == {},
                   "the full identity circuit gives Z1; premature pruning gives zero")

    rows = []
    control_errors = []
    improvement = False
    for modulus, base in [(7, 6), (7, 3), (5, 2)]:
        exp.section(f"N={modulus}, a={base}: vary only exponent width")
        for width in range(2, args.max_width + 1):
            me = ToffoliModExp(N=modulus, a=base, n_exp=width)
            order = me.order()
            alpha = (order & -order).bit_length() - 1
            beta = order >> alpha
            qc = me.build()
            mask = sum(1 << q for q in me.exp)
            traced, trace_seconds = measured(qc, me.x[0], me.exp)
            assert not traced.hit_cap, "contraction run hit the term cap"
            baseline = None
            baseline_seconds = None
            exact_error = None
            if width <= args.baseline_width:
                baseline, baseline_seconds = measured(qc, me.x[0])
                if not baseline.hit_cap:
                    terminal = {z: c for z, c in baseline.final_terms.items()
                                if not (z & mask)}
                    exact_error = max_error(traced.final_terms, terminal)
                    exp.check("P1", exact_error < 1e-11,
                              f"t={width}: every projected coefficient, error={exact_error:.2e}")
                    improvement |= traced.n_max < baseline.n_max and trace_seconds < baseline_seconds
            if width <= 4:
                coeffs = walsh.pullback_coefficients(qc, me.x[0])
                ids = np.flatnonzero(coeffs)
                terminal = {int(z): float(coeffs[z]) for z in ids if not (int(z) & mask)}
                independent_error = max_error(traced.final_terms, terminal)
                exp.check("P1", independent_error < 1e-11,
                          f"t={width}: independent Walsh error={independent_error:.2e}")
            bound = 1 << (me.exp[0] + 1)
            exp.check("P2", traced.n_max <= bound,
                      f"t={width}: peak={traced.n_max}, work+control bound={bound}")
            truth = sum(1 - 2 * (pow(base, e, modulus) & 1)
                        for e in range(1 << width)) / (1 << width)
            exp.check("P4", abs(traced.expectation - truth) < 1e-10,
                      f"t={width}: actual input expectation={traced.expectation:.12g}, true={truth:.12g}")

            collapsed, collapsed_seconds = measured(prefix(me, min(width, alpha + 1)),
                                                     me.x[0], me.exp)
            collapsed_error = max_error(traced.final_terms, collapsed.final_terms)
            if beta == 1:
                exp.check("P3", collapsed_error < 1e-11,
                          f"t={width}: keeping {min(width, alpha+1)} blocks, error={collapsed_error:.2e}")
            elif width > alpha + 1:
                control_errors.append(collapsed_error)

            row = dict(N=modulus, a=base, r=order, alpha=alpha, beta=beta,
                       t=width, qubits=me.n_qubits, logical_gates=len(qc.logical),
                       baseline_peak=baseline.n_max if baseline else None,
                       baseline_capped=baseline.hit_cap if baseline else None,
                       baseline_seconds=baseline_seconds, trace_peak=traced.n_max,
                       trace_final=len(traced.final_terms), trace_seconds=trace_seconds,
                       trace_events=traced.trace_events, exact_error=exact_error,
                       expectation=traced.expectation, classical_expectation=truth,
                       collapsed_peak=collapsed.n_max, collapsed_seconds=collapsed_seconds,
                       collapsed_error=collapsed_error)
            rows.append(row)
            exp.log(json.dumps({k: v for k, v in row.items() if k != "trace_events"}))
            output = Path(args.output)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(rows, indent=2) + "\n")

    exp.check("H1", improvement, "at least one row reduces both peak retained terms and wall time")
    exp.fail_check("C2", any(error > 1e-6 for error in control_errors),
                   f"odd-order tail-compression errors={control_errors}")
    exp.finish()


if __name__ == "__main__":
    main()
