"""Working-precision growth for one verified exact-input prefix expression.

The fixed wide fixture and seed are inherited from ``experiment_verified_sampling``.
An instrumented subclass records the first sampled prefix query that required
an Arb refinement; selecting that label is diagnostic and is not a worst-case
claim.  That label is then frozen while only Acb working precision P varies.
Radii are extracted as exact Fractions, and midpoint differences are compared
with a P=512 diagnostic ball.  A bounded r=10,t=4 label set is also checked.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The fixed-label radius/error curve exposes precision growth and the
      requested dyadic accuracy is not guaranteed by the initial P=p+16 value.
  P2  The tiny reference set has bounded labels and converges under the same
      exact-input contraction, without implying a uniform worst-case P bound.

  C1  The initial P=max(64,p+16) precision must fail the selected p=61
      radius target before adaptive refinement.

The P=512 comparison is diagnostic, not a machine-checked exact truth oracle.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from random import Random

import numpy as np

from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction
from experiments.experiment_verified_sampling import fixtures


MAX_BYTES = 32 << 20
PRECISIONS = tuple(range(64, 161))
TARGET_BITS = 61

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "fixed-label Arb radius and midpoint error are recorded across P=64..160")
exp.predict("P2", "bounded r10,t4 reference labels converge without a worst-case-P claim")
exp.must_fail("C1", "initial P=max(64,p+16)=77 is insufficient for the selected p=61 label")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def scalar_mid_rad(value):
    return binary_fraction(value.mid()), binary_fraction(value.rad().upper())


def ball_stats(ball):
    mids, rads = [], []
    for z in ball:
        mr, rr = scalar_mid_rad(z.real)
        mi, ri = scalar_mid_rad(z.imag)
        mids.extend((mr, mi))
        rads.extend((rr, ri))
    return tuple(mids), tuple(rads)


class TracingVerifiedReflectionCircuit(VerifiedReflectionCircuit):
    """Test-only tracer; production helper is not modified."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.trace = []

    def prefix_dyadic(self, sector, exponent, stop, *, accuracy_bits,
                      boundary="reflection", measured=0, output=0,
                      max_precision=None):
        result = super().prefix_dyadic(
            sector, exponent, stop, accuracy_bits=accuracy_bits,
            boundary=boundary, measured=measured, output=output,
            max_precision=max_precision)
        self.trace.append(dict(sector=sector, exponent=exponent, stop=stop,
                               boundary=boundary, measured=measured, output=output,
                               accuracy_bits=accuracy_bits,
                               working_precision=result.working_precision,
                               refinements=result.refinements))
        return result


def tiny_labels():
    labels = []
    for stop in range(5):
        for boundary in ("arithmetic", "background", "reflection"):
            for sector in range(5):
                for exponent in range(16):
                    labels.append((sector, exponent, stop, boundary, 0, 0))
    for measured in range(1, 5):
        for sector in range(5):
            for exponent in range(1 << (4-measured)):
                for output in range(1 << measured):
                    labels.append((sector, exponent, 4, "reflection", measured, output))
    return labels


def tiny_fixture():
    return VerifiedReflectionCircuit(
        10, 4, {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
        {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))})


def wide_fixture():
    return TracingVerifiedReflectionCircuit(
        (1 << 61)-2, 63,
        {1: ("x", Fraction(1, 7)), 31: ("z", Fraction(1, 5)),
         63: ("x", Fraction(1, 7))},
        {2: (0, Fraction(1, 5)), 32: (1, Fraction(1, 5)),
         63: (2, Fraction(1, 5))})


def main():
    started = time.time()
    report = {"status":"PASS", "rows":[], "controls":{},
              "selection_rule":"first requiring retry in fixed seed=624 trace; diagnostic only"}
    try:
        wide = wide_fixture()
        sample = wide.sample(Random(624), target_tv=Fraction(1, 1_000_000))
        first_retry = next((x for x in wide.trace if x["refinements"] > 0), None)
        if first_retry is None:
            raise AssertionError("fixed seed trace had no retry to select")
        report["wide_trace"] = dict(sample={k: str(v) if isinstance(v, Fraction) else v
                                             for k,v in sample.items()},
                                    calls=len(wide.trace), first_retry=first_retry)
        label = tuple(first_retry[key] for key in
                      ("sector", "exponent", "stop", "boundary", "measured", "output"))
        high_ball = wide.prefix_enclosure(*label[:3], boundary=label[3],
                                          measured=label[4], output=label[5],
                                          working_precision=512)
        high_mid, high_rad = ball_stats(high_ball)
        selected_min_midpoint_abs_squared = min(
            mid*mid + imag*imag for mid, imag in zip(high_mid[::2], high_mid[1::2]))
        curve = []
        for precision in PRECISIONS:
            ball = wide.prefix_enclosure(*label[:3], boundary=label[3],
                                         measured=label[4], output=label[5],
                                         working_precision=precision)
            mids, rads = ball_stats(ball)
            midpoint_error = max(abs(a-b) for a,b in zip(mids, high_mid))
            high_uncertainty = max(high_rad)
            curve.append(dict(working_precision=precision,
                              max_radius=str(max(rads)),
                              max_midpoint_error_vs_P512=str(midpoint_error),
                              P512_max_radius=str(high_uncertainty)))
        initial = wide.prefix_enclosure(*label[:3], boundary=label[3],
                                        measured=label[4], output=label[5],
                                        working_precision=max(64, TARGET_BITS+16))
        _, initial_rads = ball_stats(initial)
        initial_failed = max(initial_rads) > Fraction(1, 1 << (TARGET_BITS+1))
        report["selected_label"] = dict(zip(("sector","exponent","stop","boundary","measured","output"), label))
        report["radius_curve"] = curve
        report["curve_summary"] = dict(initial_radius=curve[0]["max_radius"],
                                       final_radius=curve[-1]["max_radius"],
                                       initial_midpoint_error=curve[0]["max_midpoint_error_vs_P512"],
                                       final_midpoint_error=curve[-1]["max_midpoint_error_vs_P512"],
                                       radius_decreases=Fraction(curve[-1]["max_radius"]) < Fraction(curve[0]["max_radius"]),
                                       selected_min_P512_midpoint_abs_squared=str(selected_min_midpoint_abs_squared),
                                       no_amplitude_denominator=True)
        report["curve_summary"]["first_success_P_in_scanned_range"] = next(
            (row["working_precision"] for row in curve
             if Fraction(row["max_radius"]) <= Fraction(1, 1 << (TARGET_BITS+1))), None)
        # The doubling endpoint is NOT the smallest sufficient mantissa.
        # Scan each retried query from initial P upwards, retaining exact
        # radii; this is a FIXED sampled-label diagnostic, not all-label proof.
        retry_thresholds = []
        for query in wide.trace:
            if query["refinements"] == 0:
                continue
            for candidate in range(77, query["working_precision"] + 1):
                ball = wide.prefix_enclosure(query["sector"], query["exponent"], query["stop"],
                    boundary=query["boundary"], measured=query["measured"], output=query["output"],
                    working_precision=candidate)
                radius = max(ball_stats(ball)[1])
                if radius <= Fraction(1, 1 << (TARGET_BITS+1)):
                    retry_thresholds.append(dict(query=query, first_success_P_from_77=candidate,
                                                 success_radius=str(radius)))
                    break
        report["retry_thresholds"] = retry_thresholds
        report["initial_P"] = dict(P=max(64, TARGET_BITS+16),
                                    max_radius=str(max(initial_rads)),
                                    target_radius=str(Fraction(1, 1 << (TARGET_BITS+1))),
                                    failed=initial_failed)
        # Fixed tiny set: modest exhaustive reference labels, with only two
        # diagnostic precisions and a P=512 midpoint reference.
        tiny = tiny_fixture()
        labs = tiny_labels()
        if len(labs) > 1520:
            raise MemoryError("tiny diagnostic exceeds its Python label-count cap")
        tiny_rows = []
        for precision in (64, 128):
            max_radius = Fraction(0)
            max_error = Fraction(0)
            for lab in labs:
                ball = tiny.prefix_enclosure(*lab[:3], boundary=lab[3],
                                             measured=lab[4], output=lab[5],
                                             working_precision=precision)
                ref = tiny.prefix_enclosure(*lab[:3], boundary=lab[3],
                                            measured=lab[4], output=lab[5],
                                            working_precision=512)
                mids, rads = ball_stats(ball)
                ref_mids, _ = ball_stats(ref)
                max_radius = max(max_radius, max(rads))
                max_error = max(max_error, max(abs(a-b) for a,b in zip(mids, ref_mids)))
            tiny_rows.append(dict(P=precision, labels=len(labs),
                                  max_radius=str(max_radius),
                                  max_midpoint_error_vs_P512=str(max_error)))
        report["tiny_reference"] = tiny_rows
        report["controls"] = dict(initial_P_claim_failed=initial_failed,
                                   initial_P=max(64, TARGET_BITS+16),
                                   selected_accuracy_bits=TARGET_BITS,
                                   selected_trace_refinements=first_retry["refinements"])
        p1 = (len(curve) == len(PRECISIONS)
              and report["curve_summary"]["radius_decreases"]
              and Fraction(curve[-1]["max_midpoint_error_vs_P512"])
              < Fraction(curve[0]["max_midpoint_error_vs_P512"]))
        p2 = (len(labs) <= 1520 and all(row["labels"] == len(labs) for row in tiny_rows)
              and Fraction(tiny_rows[1]["max_radius"]) < Fraction(tiny_rows[0]["max_radius"])
              and Fraction(tiny_rows[1]["max_midpoint_error_vs_P512"])
                  < Fraction(tiny_rows[0]["max_midpoint_error_vs_P512"]))
        exp.check("P1", p1, f"selected label={label}, curve points={len(curve)}")
        exp.check("P2", p2, f"tiny labels={len(labs)}, rows={tiny_rows}")
        exp.fail_check("C1", initial_failed,
                       f"initial P=77 radius {max(initial_rads)} > target {Fraction(1, 1 << (TARGET_BITS+1))}")
        if not exp.finish():
            raise AssertionError("Experiment harness failed")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type":type(exc).__name__, "message":str(exc),
                            "traceback":traceback.format_exc()}
    report["elapsed_seconds"] = time.time()-started
    path = Path("out") / f"precision_growth_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"status":report["status"], "report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
