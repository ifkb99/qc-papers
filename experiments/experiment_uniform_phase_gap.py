"""A finite phase-mesh certificate for the b=3 repeated-mixer gap.

TODO27 asks whether the qualitative C66 two-step gap can be given a useful
uniform number for the fixed embedded Rx(pi/7) mixer.  This probe samples the
periodic phase grid theta=j/N with N=2048 using the existing exact-input
finite-work builder.  For each sample it forms the two 8-by-8 traceless
Hilbert--Schmidt channel matrices and bounds both alternating 12-step products
by their outward Frobenius norm.  The channel phase is Lipschitz between mesh
points after the exact scalar phase has been removed:

    ||P_theta-P_phi||_(HS->HS) <= 24*pi*d_circle(theta,phi).

Here P is either alternating twelve-step product, not one homogeneous
channel power. The nearest circular grid point is within 1/(2N), so the
padding is 12*pi/N. Frobenius bounds only the SAMPLED matrix norm; the
between-mesh estimate is in induced norm and has no extra sqrt(8) factor.

The rational upper pi <= 355/113 is checked against Arb's pi enclosure.  The
reported uniform certificate is therefore max(mesh bound)+12*pi_upper/N; the
floating mesh values alone are not treated as a proof.  This is an HS bound,
with only ||Delta||_1 <= sqrt(3)||Delta||_F available downstream.
"""
from __future__ import annotations

import json
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


GRID = 2048
PERIOD_MULTIPLE = 3
BLOCK = 3
WIDTH = 2
WINDOW = 12
PRECISIONS = (192, 256)
THRESHOLD = Fraction(4, 5)
MAX_BYTES = 32 << 20


def report_path(prefix="uniform_phase_gap"):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def guard_entries(entries, bytes_per_entry=512):
    if type(entries) is not int or entries < 0:
        raise ValueError("entry count must be a nonnegative integer")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("uniform phase probe exceeds 32 MiB preallocation budget")


def spectral_mode(shift, alpha, sectors, acb, acb_mat, rational_phase):
    """Difference of two shift spectral projectors, traceless and conserved."""
    unit = acb_mat(BLOCK, BLOCK)
    for i in range(BLOCK):
        unit[i, i] = acb(1)
    omega = rational_phase(1, BLOCK)
    lam = rational_phase(alpha, BLOCK*sectors)
    inv0 = 1 / lam
    inv1 = 1 / (lam * omega)
    return ((unit + inv0 * shift + (inv0 ** 2) * (shift * shift))
            - (unit + inv1 * shift + (inv1 ** 2) * (shift * shift))) / BLOCK


def null_mode_audit(branches, alpha, sectors, binary_fraction, acb,
                    acb_mat, rational_phase, arb, density_forward_step):
    from experiments.experiment_forward_mixing_certificate import (
        adjoint, frobenius_upper, interval_contains, trace,
        matrix_difference_contains_zero)

    shift = branches[0][1]
    mode = spectral_mode(shift, alpha, sectors, acb, acb_mat, rational_phase)
    current = mode
    for b0, b1 in branches:
        current = density_forward_step(
            b0, b1, current, multiply=lambda left, right: left * right,
            adjoint=adjoint)
    difference = current - mode
    trace_mode = trace(mode, acb)
    lower_sq = arb(0)
    for row in range(BLOCK):
        for column in range(BLOCK):
            lower_sq += current[row, column].abs_lower() ** 2
    hermitian = True
    hermitian_difference = mode - adjoint(mode)
    for row in range(BLOCK):
        for column in range(BLOCK):
            entry = hermitian_difference[row, column]
            hermitian &= interval_contains(entry.real, 0, binary_fraction)
            hermitian &= interval_contains(entry.imag, 0, binary_fraction)
    return dict(
        difference_upper=str(frobenius_upper(difference, binary_fraction, arb)),
        difference_contains_zero=matrix_difference_contains_zero(current,mode,binary_fraction),
        norm_lower=str(binary_fraction(lower_sq.sqrt().lower())),
        trace_contains_zero=(interval_contains(trace_mode.real, 0, binary_fraction)
                             and interval_contains(trace_mode.imag, 0, binary_fraction)),
        hermitian=hermitian)


def alternating_product(first, second, identity, repeats):
    product = identity
    for _ in range(repeats):
        product = first * product
        product = second * product
    return product


def main():
    exp = Experiment("uniform_phase_gap", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "N=2048 phase mesh plus Lipschitz padding certifies both 12-step products below 4/5 at P192/P256")
    exp.predict("P2", "the no-background b=3 spectral traceless null survives the width-2 product")
    exp.must_fail("C1", "N=1 has a contractive sample but its charged Lipschitz bound cannot certify a sub-unit norm")
    report = {"status": "PASS", "precisions": {}, "controls": {}}
    started = time.perf_counter()
    try:
        from flint import acb, acb_mat, arb, ctx
        from lab.semiclassical import _density_forward_step
        from lab.verified_finite_work import VerifiedFiniteWork
        from lab.verified_prefix import (VerifiedReflectionCircuit,
                                         _rational_phase, binary_fraction)
        from experiments.experiment_forward_mixing_certificate import (
            DIM, basis_factory, basis_audit, channel_matrix, frobenius_upper,
            identity, matrix_unitarity_audit)

        # Two 8x8 channels, two products, and all retained mesh rows are
        # charged before any exact matrix is built.
        per_sample_entries = 2 * DIM * DIM + 2 * DIM * DIM + 2 * 2 * BLOCK * BLOCK
        # Both precision reports remain live, plus a serialization copy.
        retained_row_entries = 2 * len(PRECISIONS) * GRID * 6
        guard_entries(per_sample_entries + retained_row_entries + DIM * BLOCK * BLOCK + 4096)
        period = PERIOD_MULTIPLE * GRID
        backgrounds = {1: ("x", Fraction(1, 7)), 2: ("x", Fraction(1, 7))}
        circuit = VerifiedReflectionCircuit(period, WIDTH, backgrounds, {}, block_size=BLOCK)
        pi_upper = Fraction(355, 113)
        p1 = True
        p2 = True
        for precision in PRECISIONS:
            with ctx.workprec(precision):
                # Verify the classical rational pi upper bound in the same
                # Arb context used for every enclosure below.
                pi_ball = arb.pi()
                pi_upper_ball = arb(355) / arb(113)
                pi_verified = (binary_fraction(pi_upper_ball.lower())
                               >= binary_fraction(pi_ball.upper()))
                basis = basis_factory(acb, arb, acb_mat)
                basis_row = basis_audit(basis, binary_fraction, acb)
                worker = VerifiedFiniteWork(circuit, 0)
                mesh_rows = []
                mesh_max = Fraction(0)
                mesh_arg = None
                for alpha in range(GRID):
                    branches, _states, _final = worker._build(alpha)
                    if not matrix_unitarity_audit(branches,binary_fraction,acb)["ok"]:
                        raise AssertionError("mesh branch unitarity failed")
                    channels = []
                    for pair in branches:
                        channel, _audit = channel_matrix(
                            pair, basis, acb, acb_mat, binary_fraction,
                            _density_forward_step)
                        channels.append(channel)
                    product10 = alternating_product(
                        channels[0], channels[1], identity(DIM, acb_mat, acb), 6)
                    product01 = alternating_product(
                        channels[1], channels[0], identity(DIM, acb_mat, acb), 6)
                    upper10 = frobenius_upper(product10, binary_fraction, arb)
                    upper01 = frobenius_upper(product01, binary_fraction, arb)
                    local_max = max(upper10, upper01)
                    mesh_max = max(mesh_max, local_max)
                    if mesh_arg is None or local_max > mesh_arg[0]:
                        mesh_arg = (local_max, alpha, upper10, upper01)
                    mesh_rows.append(dict(alpha=alpha, theta=f"{alpha}/{GRID}",
                                         product10=str(upper10),
                                         product01=str(upper01),
                                         maximum=str(local_max)))
                    if alpha % 256 == 0:
                        print(f"  P{precision}: phase mesh {alpha}/{GRID}", flush=True)
                mesh_padding = Fraction(12) * pi_upper / GRID
                uniform_bound = mesh_max + mesh_padding
                # C65/C66 null: no backgrounds preserve a nonzero spectral
                # projector difference even though the work channel mixes
                # other directions.
                null_circuit = VerifiedReflectionCircuit(
                    period, WIDTH, {}, {}, block_size=BLOCK)
                null_worker = VerifiedFiniteWork(null_circuit, 0)
                null_branches, _states, _final = null_worker._build(1)
                null_row = null_mode_audit(
                    null_branches, 1, GRID, binary_fraction, acb, acb_mat,
                    _rational_phase, arb, _density_forward_step)
                p1 &= (pi_verified and basis_row["ok"]
                       and uniform_bound < THRESHOLD)
                p2 &= (Fraction(null_row["difference_upper"]) < Fraction(1, 1 << (precision // 2))
                       and Fraction(null_row["norm_lower"]) > Fraction(1, 2)
                       and null_row["difference_contains_zero"]
                       and null_row["trace_contains_zero"] and null_row["hermitian"])
                report["precisions"][str(precision)] = dict(
                    grid=GRID, mesh_max=str(mesh_max), mesh_arg=mesh_arg,
                    mesh_padding=str(mesh_padding), uniform_bound=str(uniform_bound),
                    threshold=str(THRESHOLD), pi_upper=str(pi_upper),
                    pi_verified=pi_verified, basis=basis_row,
                    null_mode=null_row, rows=mesh_rows)

        first = report["precisions"][str(PRECISIONS[0])]
        one_sample_bound = Fraction(first["rows"][0]["maximum"])
        small_grid_bound = one_sample_bound + 12 * pi_upper
        control_failed = one_sample_bound < 1 and small_grid_bound >= 1
        report["controls"] = dict(
            undersized_grid_N1_bound=str(small_grid_bound),
            undersized_grid_cannot_certify_subunit=control_failed,
            note="N=1 sample contraction is not itself a uniform certificate")
        exp.check("P1", p1, "mesh plus 12*pi/N padding is below 4/5 at both precisions")
        exp.check("P2", p2, "no-background conserved spectral mode")
        exp.fail_check("C1", control_failed,
                       f"N=1 charged bound {small_grid_bound} is not below one")
        report["status"] = "PASS" if p1 and p2 and control_failed else "FAIL"
        path = report_path()
        if not exp.finish(
                report_path=path,
                rows=[dict(series="precision", precision=precision,
                           **json_safe(report["precisions"][str(precision)]))
                      for precision in PRECISIONS],
                metadata=dict(period=period, block_size=BLOCK, width=WIDTH,
                              grid=GRID, threshold=str(THRESHOLD),
                              lipschitz="24*pi*d_circle; nearest distance <=1/(2N); pi<=355/113",
                              norm="induced HS norm bounded by sampled Frobenius plus induced-norm padding",
                              trace_conversion="||Delta||_1 <= sqrt(3)||Delta||_F",
                              allocation_bound_bytes=MAX_BYTES,
                              controls=json_safe(report["controls"]),
                              elapsed_seconds=time.perf_counter() - started)):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = report_path("uniform_phase_gap_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
