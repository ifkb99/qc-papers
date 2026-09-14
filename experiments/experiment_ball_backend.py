"""Bounded python-flint/Arb probe for the C60 verified-oracle boundary.

The probe uses exact fmpq inputs, acb.exp_pi_i for rational multiples of pi,
and direct man_exp -> Fraction extraction.  It checks absolute ball radii,
including cancellation, rather than relative errors that divide by a small
or zero answer.  This is an arithmetic-backend probe only: it does not add a
propagator or certify the existing complex128 sampler.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Rational Arb/Acb inputs, exact binary midpoint/radius extraction, and
      workprec restoration behave as specified without float/string conversion.
  P2  Independent rational-phase cancellation radii decrease with precision
      and meet a fixed absolute target at a bounded precision.
  P3  A normalized two-phase input remains enclosed by Arb norm bounds, while
      the report retains the explicit cap and no relative-accuracy test.
  C1  Converting a high-precision midpoint through Python float leaves the
      claimed tiny Arb ball, so that shortcut must fail.
  C2  A deliberately low explicit precision cap must fail the absolute target
      even for normalized phase inputs.

Run:
  uv run --no-project --python 3.12 --with 'numpy<2.5' --with
  'python-flint==0.9.0' python -m experiments.experiment_ball_backend
"""
from __future__ import annotations

import json
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from lab import Experiment


ABS_TARGET = Fraction(1, 10**30)
PRECISIONS = (32, 64, 128, 256)
MAX_PRECISION = max(PRECISIONS)


def report_path(prefix="ball_backend"):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def binary_fraction(x):
    """Lossless exact binary Arb value -> Fraction, with no float round-trip."""
    if not x.is_exact() or not x.is_finite():
        raise ValueError("man_exp requires an exact finite Arb value")
    mantissa, exponent = x.man_exp()
    mantissa, exponent = int(mantissa), int(exponent)
    if exponent >= 0:
        return Fraction(mantissa * (1 << exponent), 1)
    return Fraction(mantissa, 1 << (-exponent))


def midpoint_fraction(x):
    return binary_fraction(x.mid())


def outward_radius_fraction(x):
    # x.rad() is itself an Arb ball.  Its upper endpoint, not its midpoint,
    # is the certified absolute radius to carry into a rational budget.
    return binary_fraction(x.rad().upper())


def rational_phase(numerator, denominator):
    return acb(fmpq(int(numerator), int(denominator))).exp_pi_i()


def cancellation_row(precision):
    old = ctx.prec
    with ctx.workprec(precision):
        # exp(pi*i*(4/3)) = -exp(pi*i*(1/3)); each input is a unit phase.
        value = rational_phase(1, 3) + rational_phase(4, 3)
        real_mid = midpoint_fraction(value.real)
        imag_mid = midpoint_fraction(value.imag)
        real_rad = outward_radius_fraction(value.real)
        imag_rad = outward_radius_fraction(value.imag)
        contains_zero = value.contains(0)
    restored = ctx.prec == old
    return dict(precision=precision, real_mid=str(real_mid), imag_mid=str(imag_mid),
                real_radius=str(real_rad), imag_radius=str(imag_rad),
                max_radius=str(max(real_rad, imag_rad)), contains_zero=bool(contains_zero),
                restored=restored)


def normalized_phase_row(precision):
    old = ctx.prec
    with ctx.workprec(precision):
        root_two = arb(2).sqrt()
        phase_a = rational_phase(1, 7)
        phase_b = rational_phase(2, 5)
        amp_a = phase_a / root_two
        amp_b = phase_b / root_two
        norm_sq = (amp_a.real**2 + amp_a.imag**2
                   + amp_b.real**2 + amp_b.imag**2)
        lower = binary_fraction(norm_sq.lower())
        upper = binary_fraction(norm_sq.upper())
    return dict(precision=precision, norm_lower=str(lower), norm_upper=str(upper),
                contains_one=lower <= 1 <= upper, restored=ctx.prec == old)


def main():
    exp = Experiment("ball_backend", doc=__doc__)
    exp.predict("P1", "exact rational Arb/Acb conversion and workprec restoration")
    exp.predict("P2", "absolute cancellation radii converge under bounded precision")
    exp.predict("P3", "normalized phase inputs remain enclosed and cap is explicit")
    exp.must_fail("C1", "naive float midpoint conversion leaves a tiny claimed ball")
    exp.must_fail("C2", "low explicit precision cap misses the absolute target")

    old_precision = ctx.prec
    rows = []
    try:
        # Exact rational input and binary extraction checks.
        with ctx.workprec(192):
            rational = arb(fmpq(1, 3))
            midpoint = midpoint_fraction(rational)
            radius = outward_radius_fraction(rational)
            phase = rational_phase(1, 3)
            phase_midpoints = (midpoint_fraction(phase.real), midpoint_fraction(phase.imag))
        exp.check("P1", isinstance(midpoint, Fraction) and isinstance(radius, Fraction)
                  and midpoint.denominator & (midpoint.denominator - 1) == 0
                  and radius >= 0 and all(isinstance(x, Fraction) for x in phase_midpoints)
                  and ctx.prec == old_precision,
                  f"binary midpoint={midpoint}, radius={radius}, "
                  f"phase_midpoints={phase_midpoints}, restored={ctx.prec == old_precision}")

        for precision in PRECISIONS:
            row = cancellation_row(precision)
            rows.append(row)
        radii = [Fraction(row["max_radius"]) for row in rows]
        exp.check("P2", all(rows[i]["contains_zero"] and rows[i]["restored"]
                             for i in range(len(rows)))
                  and all(radii[i+1] < radii[i] for i in range(len(radii)-1))
                  and radii[-1] < ABS_TARGET <= radii[1],
                  f"absolute_target={ABS_TARGET}, radii={[str(x) for x in radii]}")

        normalized = normalized_phase_row(128)
        exp.check("P3", normalized["contains_one"] and normalized["restored"]
                  and MAX_PRECISION == 256,
                  f"normalized={normalized}, explicit_cap={MAX_PRECISION}")

        # Deliberate shortcut failure: at 256-bit Arb precision, float64 cannot
        # preserve the midpoint inside the tiny exact Arb enclosure.
        with ctx.workprec(MAX_PRECISION):
            exact_ball = arb(fmpq(1, 3))
            naive = arb(float(exact_ball.mid()))
            inside = exact_ball.lower() <= naive and naive <= exact_ball.upper()
        exp.fail_check("C1", not bool(inside),
                       f"naive_float_inside_claimed_ball={bool(inside)}, "
                       f"radius={outward_radius_fraction(exact_ball)}")

        # Normalized phases still require the explicit absolute-accuracy cap:
        # precision 32 does not meet the fixed 1e-30 cancellation target.
        low = normalized_phase_row(32)
        low_cancel = rows[0]
        exp.fail_check("C2", low["contains_one"] and
                       Fraction(low_cancel["max_radius"]) > ABS_TARGET
                       and MAX_PRECISION == 256,
                       f"normalized_low_precision={low}, low_cancel_radius="
                       f"{low_cancel['max_radius']}, target={ABS_TARGET}, cap={MAX_PRECISION}")
    finally:
        ctx.prec = old_precision

    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(
        package="python-flint==0.9.0", absolute_target=str(ABS_TARGET),
        precisions=PRECISIONS, explicit_max_precision=MAX_PRECISION,
        conversion="fmpq -> Arb/Acb; man_exp -> Fraction; no float/string midpoint path",
        phase_api="acb.exp_pi_i", radius_api="arb.rad().upper()",
        context_api="ctx.workprec", certifies_float_sampler=False))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        path = report_path("ball_backend_failure")
        path.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                        traceback=traceback.format_exc()), indent=2) + "\n")
        raise
