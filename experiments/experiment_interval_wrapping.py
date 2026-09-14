"""Isolated Arb/Acb probe of rectangular interval wrapping.

The input is the rectangular set
    (1 +/- eps) + i(0 +/- eps),  eps=2^-20,
and each step multiplies it by the exact unit phase exp(+/- pi*i/4), with
signs alternating.  The exact combined phase is the identity after every
two steps, but sequential rectangular arithmetic forgets the correlation
between real and imaginary coordinates.  At pi/4, equal axis radii therefore
grow by sqrt(2) per multiplication; an alternating pair grows by two even
though the mathematical product is unitary and exactly the identity.

P1: high-precision Acb radii are input-width dominated, and one-shot combined
    phases track the exact projected rectangle bounds.
P2: sequential axis radii inflate across steps 1..64 while the even-step
    combined phase remains the original identity enclosure.
P3: exact rational-radius comparisons and the analytic projected-rectangle
    bounds agree without float-based certification.
C1: unitary multiplication implies a nonincreasing rectangular axis radius.
C2: an alternating +pi/4,-pi/4 pair must reset independently propagated
    rectangular boxes to the initial radius.

This is standard interval-analysis wrapping, not new quantum behavior and not
a circuit propagator.  Arb/Acb enclosures use fixed high precision; all pass
predicates convert finite binary upper radii to exact Fraction values.
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


STEPS = 64
WORKING_PRECISION = 256
EPS_BITS = 20
MAX_BYTES = 16 << 20


def guard(shape, itemsize: int, label: str) -> None:
    if math.prod(int(x) for x in shape) * int(itemsize) > MAX_BYTES:
        raise MemoryError(f"{label} exceeds 16 MiB")


def report_path(prefix: str = "interval_wrapping") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def binary_upper_fraction(value) -> Fraction:
    """Exact Fraction from a finite nonnegative Arb upper endpoint."""
    upper = value.upper()
    if not upper.is_finite():
        raise ArithmeticError("nonfinite Arb radius")
    mantissa, exponent = (int(v) for v in upper.man_exp())
    return (Fraction(mantissa << exponent) if exponent >= 0
            else Fraction(mantissa, 1 << -exponent))


def radius_pair(z) -> tuple[Fraction, Fraction]:
    return binary_upper_fraction(z.real.rad()), binary_upper_fraction(z.imag.rad())


def main() -> None:
    from flint import acb, arb, ctx, fmpq

    exp = Experiment("interval_wrapping", doc=__doc__)
    exp.predict("P1", "input-width dominates high-precision Acb and one-shot bounds track the true rectangle")
    exp.predict("P2", "sequential rectangular radii inflate through 64 alternating steps")
    exp.predict("P3", "exact Fraction radius comparisons agree with analytic identity/rotation bounds")
    exp.must_fail("C1", "unitarity makes rectangular axis radii nonincreasing")
    exp.must_fail("C2", "an alternating inverse pair resets independently propagated boxes")

    with ctx.workprec(WORKING_PRECISION):
        if STEPS > 64:
            raise ValueError("isolated diagnostic is capped at 64 steps; no byte/RSS claim")
        eps = Fraction(1, 1 << EPS_BITS)
        eps_arb = arb(1) / arb(1 << EPS_BITS)
        initial = acb(arb(arb(1), eps_arb), arb(arb(0), eps_arb))
        input_radii = radius_pair(initial)
        plus = acb(fmpq(1, 4)).exp_pi_i()
        minus = acb(fmpq(-1, 4)).exp_pi_i()
        phase_radius = max(radius_pair(plus) + radius_pair(minus))

        sequential = []
        current = initial
        for step in range(1, STEPS + 1):
            current = (plus if step % 2 else minus) * current
            sequential.append(radius_pair(current))

        combined = []
        for step in range(1, STEPS + 1):
            net_turn = 1 if step % 2 else 0
            phase = acb(fmpq(net_turn, 4)).exp_pi_i()
            combined.append(radius_pair(phase * initial))

        expected_squared = []
        for step in range(1, STEPS + 1):
            axis_radius = eps if step % 2 == 0 else None
            # A pi/4 projection of an eps-by-eps rectangle has radius
            # sqrt(2)*eps; retain its exact squared value instead of sqrt.
            expected_squared.append(eps * eps if axis_radius is not None
                                   else 2 * eps * eps)

        rows = []
        sequential_pass = True
        combined_pass = True
        for step, (seq, one_shot, expected_sq) in enumerate(
                zip(sequential, combined, expected_squared), start=1):
            seq_sq = seq[0] * seq[0]
            combined_sq = one_shot[0] * one_shot[0]
            ideal_seq_sq = (1 << step) * eps * eps
            # Fixed 30-bit radius rounding compounds slightly. This allowance
            # is a frozen diagnostic predicate, not a backend-wide theorem.
            sequential_pass &= all(ideal_seq_sq <= radius**2 <=
                ideal_seq_sq * Fraction(10001, 10000)**2 for radius in seq)
            # One-shot radii have only high-precision transcendental rounding
            # beyond the exact projected rectangle, bounded by a loose factor.
            combined_pass &= all(expected_sq <= radius**2 <=
                expected_sq * Fraction(10001, 10000)**2 for radius in one_shot)
            rows.append(dict(
                series="step", step=step,
                sequential_real_radius=str(seq[0]),
                sequential_imag_radius=str(seq[1]),
                combined_real_radius=str(one_shot[0]),
                combined_imag_radius=str(one_shot[1]),
                sequential_radius_over_eps=float(seq[0] / eps),
                combined_radius_over_eps=float(one_shot[0] / eps),
                expected_axis_radius_squared=str(expected_sq),
                ideal_sequential_axis_radius_squared=str(ideal_seq_sq),
                sequential_radius_squared=str(seq_sq),
                combined_radius_squared=str(combined_sq),
            ))

        first = sequential[0][0]
        second = sequential[1][0]
        final = sequential[-1][0]
        even_combined = combined[1][0]
        odd_combined = combined[0][0]
        input_radius = input_radii[0]
        symmetry_pass = all(
            abs(seq[0] - seq[1]) <= max(seq[0], seq[1]) / 100_000_000
            for seq in sequential)
        p1 = (input_radii[0] == input_radii[1]
              and phase_radius < eps / 100
              and first > Fraction(6, 5) * input_radius
              and odd_combined <= 2 * input_radius)
        p2 = final > (1 << 20) * eps and even_combined <= 2 * eps
        p3 = combined_pass and symmetry_pass and sequential_pass
        exp.section("P1 input-width and one-shot combined phase")
        exp.check("P1", p1,
                  f"input_radii={input_radii}, phase_radius={phase_radius}, "
                  f"first/eps={float(first/eps):.6g}, odd_combined/eps={float(odd_combined/eps):.6g}")
        exp.section("P2 sequential wrapping growth")
        exp.check("P2", p2,
                  f"step64/eps={float(final/eps):.6g}, "
                  f"even_combined/eps={float(even_combined/eps):.6g}")
        exp.section("P3 exact radius comparisons")
        exp.check("P3", p3,
                  f"analytic_bound_checks={combined_pass}, sequential_growth={sequential_pass}, sequential_axis_symmetry={symmetry_pass}")
        exp.section("must-fail controls")
        exp.fail_check("C1", first > eps,
                       f"one-step rectangular radius={float(first/eps):.6g}*eps")
        exp.fail_check("C2", second > Fraction(3, 2) * eps,
                       f"two-step inverse-pair radius={float(second/eps):.6g}*eps")

        report = report_path()
        exp.finish(report_path=report, rows=rows,
                   metadata=dict(steps=STEPS, working_precision=WORKING_PRECISION,
                                 epsilon=str(eps), input_radii=tuple(str(x) for x in input_radii),
                                 phase_radius=str(phase_radius),
                                 true_even_axis_radius=str(eps),
                                 true_odd_axis_radius_squared=str(2 * eps * eps),
                                 arithmetic="exact Fraction predicates from Arb binary upper radii",
                                 interpretation="rectangular interval wrapping; not quantum dynamics"))
        print(f"report: {report}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("interval_wrapping_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
