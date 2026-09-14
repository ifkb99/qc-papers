"""Audit the opt-in norm-ball verified prefix enclosure.

The fixed exact rational-pi fixture is r=10,b=2,t=4 with Rx(pi/7) at s=1,
Rz(pi/5) at s=3, and reflections K_2(q=0,pi/5), K_3(q=1,pi/5).  Both the
unchanged rectangular enclosure and the opt-in norm enclosure are checked on
all 1,520 distinct prefix labels at requested accuracy bits 4, 12 and 24.
The existing full-r product is only a float diagnostic; P=512 Acb midpoints
and radii are recorded separately and interval overlap is not treated as a
proof.  The production default is not modified here.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Both enclosure modes meet the coordinatewise dyadic promise against the
      float diagnostic at all requested accuracies; P=512 same-oracle bounds
      provide a higher-precision check of the same contraction.
  P2  Analytical exact-zero and Bell/H amplitudes remain exact in both modes,
      with normalized resource/label accounting.

  C1  A low-P norm-ball midpoint-only claim must fail when its midpoint lies
      outside the high-P reference ball.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction
from lab.coherent_routes import CoherentReflectionCircuit
from experiments.experiment_coherent_route_sampling import direct_prefix


PERIOD, BLOCK, WIDTH = 10, 2, 4
MAX_BYTES = 32 << 20
ACCURACIES = (4, 12, 24)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "rectangular and norm modes meet coordinatewise dyadic checks with high-P diagnostics")
exp.predict("P2", "analytical zero/Bell controls and bounded label accounting pass in both modes")
exp.must_fail("C1", "low-P norm midpoint is excluded by a high-P Acb diagnostic ball")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def rx(theta):
    return np.array([[np.cos(theta/2), -1j*np.sin(theta/2)],
                     [-1j*np.sin(theta/2), np.cos(theta/2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j*theta/2), np.exp(1j*theta/2)]).astype(complex)


def labels():
    result = []
    for stop in range(WIDTH+1):
        for boundary in ("arithmetic", "background", "reflection"):
            for sector in range(PERIOD//BLOCK):
                for exponent in range(1 << WIDTH):
                    result.append((sector, exponent, stop, boundary, 0, 0))
    for measured in range(1, WIDTH+1):
        for sector in range(PERIOD//BLOCK):
            for exponent in range(1 << (WIDTH-measured)):
                for output in range(1 << measured):
                    result.append((sector, exponent, WIDTH, "reflection",
                                   measured, output))
    return result


def float_reference():
    return CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH,
        {1: rx(np.pi/7), 3: rz(np.pi/5)},
        {2: (0, np.pi/5), 3: (1, np.pi/5)})


def verified(mode):
    return VerifiedReflectionCircuit(
        PERIOD, WIDTH,
        {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
        {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))},
        enclosure_mode=mode)


def dyadic(item):
    scale = 1 << item.accuracy_bits
    return np.asarray([complex(re/scale, im/scale)
                       for re, im in item.coordinates], complex)


def ball_components(ball):
    mids, rads = [], []
    for z in ball:
        for component in (z.real, z.imag):
            mids.append(binary_fraction(component.mid()))
            rads.append(binary_fraction(component.rad().upper()))
    return tuple(mids), tuple(rads)


def main():
    started = time.time()
    report = {"status":"PASS", "labels":0, "modes":{}, "controls":{},
              "float_reference":"diagnostic only; no interval containment claim"}
    try:
        labs = labels()
        if len(labs) > 1520:
            raise MemoryError("distinct-label cap exceeded")
        # Conservative cache-object allowance, not a NumPy payload or RSS
        # claim: two 4-coordinate Fraction records per finite label.
        if len(labs) * 8192 > MAX_BYTES:
            raise MemoryError("tiny high-precision cache allowance exceeded")
        report["labels"] = len(labs)
        ref = float_reference()
        mode_results = {}
        all_checks = []
        for mode in ("rectangular", "norm"):
            exact = verified(mode)
            rows = []
            high_mid_cache = {}
            high_rad_cache = {}
            for label in labs:
                ball = exact.prefix_enclosure(*label[:3], boundary=label[3],
                                              measured=label[4], output=label[5],
                                              working_precision=512)
                high_mid_cache[label], high_rad_cache[label] = ball_components(ball)
            for bits in ACCURACIES:
                max_complex_error = 0.
                max_coordinate_error = 0.
                max_declared_radius = Fraction(0)
                max_high_mid_error = Fraction(0)
                max_high_error_upper = Fraction(0)
                high_radius_max = Fraction(0)
                for label in labs:
                    item = exact.prefix_dyadic(*label[:3], accuracy_bits=bits,
                                               boundary=label[3], measured=label[4],
                                               output=label[5], max_precision=512)
                    got = dyadic(item)
                    float_value = np.asarray(direct_prefix(
                        ref, label[0], label[1], label[2], label[3],
                        label[4], label[5]), complex)
                    max_complex_error = max(max_complex_error,
                                            float(np.max(np.abs(got-float_value))))
                    max_coordinate_error = max(max_coordinate_error,
                                               float(np.max(np.abs(np.real(got-float_value)))),
                                               float(np.max(np.abs(np.imag(got-float_value)))))
                    max_declared_radius = max(max_declared_radius,
                                              Fraction(item.max_ball_radius))
                    high_mid, high_rad = high_mid_cache[label], high_rad_cache[label]
                    got_mid = tuple(Fraction(x, 1 << bits)
                                    for pair in item.coordinates for x in pair)
                    max_high_mid_error = max(max_high_mid_error,
                                             max(abs(a-b) for a,b in zip(got_mid, high_mid)))
                    high_radius_max = max(high_radius_max, max(high_rad))
                    max_high_error_upper = max(max_high_error_upper,
                        max(abs(a-b)+r for a,b,r in zip(got_mid, high_mid, high_rad)))
                rows.append(dict(accuracy_bits=bits,
                                 max_complex_error=max_complex_error,
                                 max_coordinate_error=max_coordinate_error,
                                 max_declared_radius=str(max_declared_radius),
                                 max_midpoint_error_vs_P512=str(max_high_mid_error),
                                 max_P512_radius=str(high_radius_max),
                                 outward_coordinate_error_upper=str(max_high_error_upper),
                                 coordinate_bound_pass=(max_coordinate_error <= 2**(-bits)+1e-12
                                     and max_high_error_upper <= Fraction(1, 1 << bits))))
            mode_results[mode] = rows
            all_checks.extend(row["coordinate_bound_pass"] for row in rows)
        report["modes"] = mode_results
        # Analytic exact-zero and Bell/H ground truths in both enclosure modes.
        analytic = {}
        for mode in ("rectangular", "norm"):
            zero = VerifiedReflectionCircuit(2, 0, {}, {}, enclosure_mode=mode)
            zball = zero.prefix_enclosure(0, 0, 0, boundary="arithmetic", working_precision=128)
            zdy = zero.prefix_dyadic(0, 0, 0, accuracy_bits=24,
                                     boundary="arithmetic", max_precision=512)
            bell = VerifiedReflectionCircuit(2, 1, {}, {}, enclosure_mode=mode)
            bell_values = [bell.prefix_dyadic(0, 0, 1, accuracy_bits=20,
                                              measured=1, output=y).coordinates
                           for y in range(2)]
            analytic[mode] = dict(exact_zero_ball=zball[1].is_zero() and zball[1].is_exact(),
                                  exact_zero_dyadic=zdy.coordinates[1] == (0, 0),
                                  bell_exact=bell_values == [((1 << 19, 0), (1 << 19, 0)),
                                                             ((1 << 19, 0), (-(1 << 19), 0))])
        report["analytic"] = analytic
        low = verified("norm")
        chosen = None
        low_mid = high_mid = high_rad = None
        for candidate in labs:
            low_ball = low.prefix_enclosure(*candidate[:3], boundary=candidate[3],
                                            measured=candidate[4], output=candidate[5],
                                            working_precision=16)
            high_ball = low.prefix_enclosure(*candidate[:3], boundary=candidate[3],
                                             measured=candidate[4], output=candidate[5],
                                             working_precision=512)
            candidate_low_mid, _ = ball_components(low_ball)
            candidate_high_mid, candidate_high_rad = ball_components(high_ball)
            if any(abs(a-b) > radius for a,b,radius in
                   zip(candidate_low_mid, candidate_high_mid, candidate_high_rad)):
                chosen = candidate
                low_mid, high_mid, high_rad = (candidate_low_mid,
                                               candidate_high_mid,
                                               candidate_high_rad)
                break
        if chosen is None:
            raise AssertionError("no low-P norm midpoint excluded by P512 diagnostic ball")
        midpoint_outside = True
        report["controls"] = dict(low_P=16, low_label=chosen,
                                   midpoint_only_outside_high_ball=midpoint_outside,
                                   low_midpoint_vs_high_midpoint=str(max(abs(a-b) for a,b in zip(low_mid, high_mid))),
                                   high_ball_max_radius=str(max(high_rad)),
                                   analytic=analytic)
        p1 = all(all_checks)
        p2 = (all(all(v.values()) for v in analytic.values())
              and len(labs) == 1520 and report["labels"] * 8192 <= MAX_BYTES)
        exp.check("P1", p1, f"coordinatewise rows pass={p1}; interval overlap not used")
        exp.check("P2", p2, f"labels={len(labs)}, analytic={analytic}")
        exp.fail_check("C1", midpoint_outside,
                       f"low P16 midpoint outside P512 diagnostic ball={midpoint_outside}")
        if not exp.finish():
            raise AssertionError("Experiment harness failed")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type":type(exc).__name__, "message":str(exc),
                            "traceback":traceback.format_exc()}
    report["elapsed_seconds"] = time.time()-started
    path = Path("out") / f"norm_prefix_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"status":report["status"], "report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
