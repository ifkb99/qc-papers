"""Verified exact-input reference audit for coherent-route prefix amplitudes.

The exact-input oracle uses rational multiples of pi and Arb/Acb.  The
independent comparison is the existing bounded full-r product reference from
``experiment_coherent_route_sampling``; this file does not add a propagator.
The fixed fixture is r=10,b=2,t=4 with Rx(pi/7) at s=1, Rz(pi/5) at s=3,
and K_2(q=0,pi/5), K_3(q=1,pi/5).  Every requested prefix label is checked at
accuracy bits 4, 12 and 24, including all partial-QFT output prefixes.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Dyadic prefixes converge with the declared absolute coordinate bound;
      increasing requested bits reduces the independent-reference discrepancy.
  P2  128-bit Acb midpoints agree with the independent float diagnostic at the
      stated tolerance, exact zeros are handled absolutely, and query
      order/context restoration do not change deterministic results.
  P3  Exact rational input and explicit working precision remain within the
      declared tiny reference-array and label-count caps; no RSS claim is made.

  C1  Wrong reflection sign and a moved-insertion schedule must disagree with
      the verified prefixes (the latter is not a same-insertion gate-order test).
  C2  Discarding all imaginary coordinates must disagree on this complex
      fixture.

This is an input/oracle reference audit, not a certificate for complex128.
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
from lab.verified_prefix import VerifiedReflectionCircuit
from experiments.experiment_coherent_route_sampling import direct_prefix
from lab.coherent_routes import CoherentReflectionCircuit


PERIOD, BLOCK, WIDTH = 10, 2, 4
MAX_BYTES = 32 << 20
ACCURACIES = (4, 12, 24)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "dyadic absolute-error prefixes converge across requested accuracy bits")
exp.predict("P2", "Acb midpoints meet the stated float diagnostic tolerance and deterministic labels are stable")
exp.predict("P3", "exact-input/reference setup stays within bounded arrays and labels; no RSS claim")
exp.must_fail("C1", "wrong reflection sign and moved-insertion schedule disagree with verified prefixes")
exp.must_fail("C2", "discarding imaginary coordinates disagrees on complex prefixes")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def rx(theta):
    return np.array([[np.cos(theta/2), -1j*np.sin(theta/2)],
                     [-1j*np.sin(theta/2), np.cos(theta/2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j*theta/2), np.exp(1j*theta/2)]).astype(complex)


def circuit(sign=1, moved_insertion=False):
    # The direct reference's normal order is W then K.  The sign variant is
    # the deliberately wrong q=+/-1 route control.
    backgrounds = {1: rx(np.pi/7), 3: rz(np.pi/5)}
    reflections = {2: (0, np.pi/5), 3: (sign, np.pi/5)}
    if moved_insertion:
        # Deliberately move Rz from insertion 3 to insertion 2.  This is an
        # insertion-schedule control, not a claim that W/K order at one
        # insertion fails (those supplied operations commute in this model).
        backgrounds = {1: rx(np.pi/7), 2: rz(np.pi/5)}
        reflections = {2: (0, np.pi/5), 3: (sign, np.pi/5)}
    return CoherentReflectionCircuit(PERIOD, BLOCK, WIDTH, backgrounds, reflections)


def verified():
    return VerifiedReflectionCircuit(
        PERIOD, WIDTH,
        {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
        {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))})


def labels():
    result = []
    for stop in range(WIDTH+1):
        for boundary in ("arithmetic", "background", "reflection"):
            for sector in range(PERIOD//BLOCK):
                for exponent in range(1 << WIDTH):
                    result.append((sector, exponent, stop, boundary, 0, 0))
    for measured in range(1, WIDTH+1):
        for sector in range(PERIOD//BLOCK):
            # Once measured low output bits are selected, only the remaining
            # low exponent labels are distinct; ignored high labels are tested
            # separately as an invariance control.
            for exponent in range(1 << (WIDTH-measured)):
                for output in range(1 << measured):
                    result.append((sector, exponent, WIDTH, "reflection",
                                   measured, output))
    return result


def direct_value(ref, label):
    sector, exponent, stop, boundary, measured, output = label
    return np.asarray(direct_prefix(ref, sector, exponent, stop, boundary,
                                    measured, output), dtype=complex)


def dyadic_value(item):
    scale = 1 << item.accuracy_bits
    return np.asarray([complex(re/scale, im/scale)
                       for re, im in item.coordinates], complex)


def midpoint_float_diagnostic(ball, reference, tolerance=2e-12):
    worst = 0.
    over_tolerance = 0
    for z, target in zip(ball, reference):
        for component, exact in ((z.real, target.real), (z.imag, target.imag)):
            midpoint = float(component.mid())
            worst = max(worst, abs(midpoint-exact))
            if abs(midpoint-exact) > tolerance:
                over_tolerance += 1
    return worst, over_tolerance


def main():
    started = time.time()
    report = {"status":"PASS", "labels":0, "rows":[], "controls":{},
              "fixture":"r10,b2,t4,Rx(pi/7)@1,Rz(pi/5)@3,K2(q0)@2,K3(q1)@3"}
    try:
        labs = labels()
        if len(labs) > 2000:
            raise MemoryError("reference Python label-count cap exceeded")
        guard((PERIOD, PERIOD), label="largest live full-r matrix")
        report["labels"] = len(labs)
        exact = verified()
        float_ref = circuit()
        wrong_sign = circuit(sign=-1)
        moved_insertion = circuit(sign=1, moved_insertion=True)
        rows = []
        by_accuracy = {}
        exact_zero_count = 0
        wrong_sign_error = 0.
        wrong_order_error = 0.
        imaginary_drop_error = 0.
        midpoint_bad = 0
        max_midpoint_error = 0.
        near_zero_labels = set()
        max_coordinate_errors = {}
        # Exhaustive finite labels; each exact oracle query is rebuilt from
        # exact rational input, so order-dependent caches cannot hide defects.
        for bits in ACCURACIES:
            max_error = 0.
            max_coordinate_error = 0.
            max_declared = 0.
            zero_count = 0
            refinement_max = 0
            for label in labs:
                sector, exponent, stop, boundary, measured, output = label
                item = exact.prefix_dyadic(sector, exponent, stop,
                                           accuracy_bits=bits,
                                           boundary=boundary, measured=measured,
                                           output=output, max_precision=512)
                got = dyadic_value(item)
                ref = direct_value(float_ref, label)
                max_error = max(max_error, float(np.max(np.abs(got-ref))))
                max_coordinate_error = max(
                    max_coordinate_error,
                    float(np.max(np.abs(np.real(got-ref)))),
                    float(np.max(np.abs(np.imag(got-ref)))))
                max_declared = max(max_declared, float(item.max_ball_radius))
                refinement_max = max(refinement_max, item.refinements)
                if all(re == 0 and im == 0 for re, im in item.coordinates):
                    zero_count += 1
                if float(np.max(np.abs(ref))) < 1e-14:
                    near_zero_labels.add(label)
                    if float(np.max(np.abs(got))) > 2**(-bits) + 1e-12:
                        raise AssertionError("near-zero diagnostic exceeded absolute tolerance")
                if bits == ACCURACIES[-1]:
                    wrong_sign_error = max(wrong_sign_error,
                                           float(np.max(np.abs(direct_value(wrong_sign, label)-ref))))
                    wrong_order_error = max(wrong_order_error,
                                            float(np.max(np.abs(direct_value(moved_insertion, label)-ref))))
                    imaginary_drop_error = max(imaginary_drop_error,
                                               float(np.max(np.abs(got.real-ref))))
                if bits == 12:
                    ball = exact.prefix_enclosure(sector, exponent, stop,
                                                  boundary=boundary,
                                                  measured=measured, output=output,
                                                  working_precision=128)
                    excess, bad = midpoint_float_diagnostic(ball, ref)
                    max_midpoint_error = max(max_midpoint_error, excess)
                    midpoint_bad += bad
            max_coordinate_errors[bits] = max_coordinate_error
            by_accuracy[bits] = dict(max_abs_error=max_error,
                                     max_coordinate_error=max_coordinate_error,
                                     max_ball_radius=max_declared,
                                     all_zero_quantized_vectors=zero_count,
                                     max_refinements=refinement_max)
        # Determinism in reverse label order and context restoration are checked
        # on a representative exact label set at the middle requested accuracy.
        order_labels = labs[::17]
        forward = [exact.prefix_dyadic(*label[:3], accuracy_bits=12,
                                       boundary=label[3], measured=label[4],
                                       output=label[5], max_precision=512).coordinates
                    for label in order_labels]
        reverse = [exact.prefix_dyadic(*label[:3], accuracy_bits=12,
                                       boundary=label[3], measured=label[4],
                                       output=label[5], max_precision=512).coordinates
                   for label in reversed(order_labels)]
        deterministic = all(a == b for a, b in zip(forward, reversed(reverse)))
        ignored_highbit_invariance = True
        for measured in range(1, WIDTH):
            for sector in range(PERIOD//BLOCK):
                for exponent in range(1 << (WIDTH-measured)):
                    label = (sector, exponent, WIDTH, "reflection", measured, 0)
                    low = exact.prefix_dyadic(sector, exponent, WIDTH,
                                              accuracy_bits=12, boundary="reflection",
                                              measured=measured, output=0,
                                              max_precision=512).coordinates
                    high = exact.prefix_dyadic(sector, exponent + (1 << (WIDTH-measured)),
                                               WIDTH, accuracy_bits=12,
                                               boundary="reflection", measured=measured,
                                               output=0, max_precision=512).coordinates
                    ignored_highbit_invariance &= low == high
        from flint import ctx
        context_before = ctx.prec
        with ctx.workprec(77):
            context_inside_before = ctx.prec
            exact.prefix_enclosure(0, 0, 1, boundary="arithmetic", working_precision=128)
            context_inside_after = ctx.prec
            try:
                exact.prefix_dyadic(0, 0, 1, accuracy_bits=24,
                                    boundary="arithmetic", max_precision=16)
            except ArithmeticError:
                pass
            context_after_exception = ctx.prec
        context_restored = (ctx.prec == context_before
                            and context_inside_before == 77
                            and context_inside_after == 77
                            and context_after_exception == 77)
        empty = VerifiedReflectionCircuit(10, 0, {}, {})
        zero_ball = empty.prefix_enclosure(0, 0, 0, boundary="arithmetic", working_precision=128)
        zero_dyadic = empty.prefix_dyadic(0, 0, 0, accuracy_bits=24,
                                          boundary="arithmetic", max_precision=512)
        exact_zero_gate = (zero_ball[1].is_zero() and zero_ball[1].is_exact()
                           and zero_dyadic.coordinates[1] == (0, 0))
        report["rows"] = [{"accuracy_bits": bits, **row}
                           for bits, row in by_accuracy.items()]
        report["midpoint_float_diagnostic"] = dict(working_precision=128,
                                                     tolerance=2e-12,
                                                     max_midpoint_minus_float=max_midpoint_error,
                                                     over_tolerance=midpoint_bad)
        report["resources"] = dict(label_count=len(labs),
                                    largest_full_r_matrix_payload_bytes=PERIOD**2*16,
                                    per_array_payload_cap_bytes=MAX_BYTES,
                                    python_label_count_cap=2000,
                                    no_rss_claim=True)
        report["controls"] = dict(wrong_sign_max_error=wrong_sign_error,
                                   moved_insertion_max_error=wrong_order_error,
                                   discarded_imaginary_max_error=imaginary_drop_error,
                                   deterministic_reverse_order=deterministic,
                                   near_zero_prefixes=len(near_zero_labels),
                                   ignored_highbit_invariance=ignored_highbit_invariance,
                                   context_restored=context_restored,
                                   analytical_exact_zero_gate=exact_zero_gate)
        p1 = (max_coordinate_errors[4] >= max_coordinate_errors[12]
              and max_coordinate_errors[12] >= max_coordinate_errors[24]
              and all(row["max_coordinate_error"] <= 2**(-bits)+1e-12
                      for bits, row in by_accuracy.items()))
        p2 = (midpoint_bad == 0 and deterministic and ignored_highbit_invariance
              and context_restored and exact_zero_gate
              and max_midpoint_error <= 2e-12)
        p3 = (PERIOD**2*16 <= MAX_BYTES and report["labels"] <= 2000
              and all(row["max_refinements"] < 20 for row in by_accuracy.values()))
        exp.check("P1", p1, f"dyadic errors={by_accuracy}")
        exp.check("P2", p2, f"midpoint over-tolerance={midpoint_bad}, context={context_restored}, exact-zero={exact_zero_gate}")
        exp.check("P3", p3, f"labels={report['labels']}, full-r matrix bytes={PERIOD**2*16}, no RSS claim")
        exp.fail_check("C1", wrong_sign_error > 1e-8 and wrong_order_error > 1e-8,
                       f"wrong sign/moved-insertion errors={wrong_sign_error:.6g}/{wrong_order_error:.6g}")
        exp.fail_check("C2", imaginary_drop_error > 1e-8,
                       f"discarded-imaginary error={imaginary_drop_error:.6g}")
        if not exp.finish():
            raise AssertionError("Experiment harness failed")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time()-started
    path = Path("out") / f"verified_prefix_reference_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"status":report["status"], "report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
