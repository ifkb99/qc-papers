"""Certified b=3 forward-channel window bounds for TODO26.

For r=9, b=3 and repeated exact-input Rx(pi/7) blocks, this probe takes the
existing ``VerifiedFiniteWork._build`` branch matrices and forms the induced
real channel on the eight-dimensional traceless Hermitian subspace.  The
basis is an Acb enclosure of the normalized Gell-Mann basis.  Each 12-channel
window is bounded by its outward Hilbert--Schmidt Frobenius matrix norm,
which is an upper bound on the induced Hilbert--Schmidt norm.  It is not a
trace-norm contraction certificate; the consequence conversion is
||Delta||_1 <= sqrt(3)||Delta||_F.

P1: at both working precisions P=192 and P=256, every one of the 18
    (three sectors by six period-offset) repeated-Rx windows has certified
    Frobenius upper bound below 4/5.
P2: the Acb basis is traceless and orthonormal, every branch is unitary, and
    the shift/channel schedule agrees between offsets six apart.
P3: the no-background null has an explicit conserved traceless spectral-mode
    difference through a 12-channel product.
C1: a mixed fixed point by itself implies traceless contraction.

The arithmetic is bounded before allocation and uses no sampler, histogram,
new propagator, or floating eigensolver.  Exact rational-pi inputs are
reconstructed by the existing verified builder at each working precision.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


PERIOD = 9
BLOCK = 3
WINDOW = 12
OFFSETS = 6
WIDTH = WINDOW + OFFSETS
PRECISIONS = (192, 256)
THRESHOLD = Fraction(4, 5)
MAX_BYTES = 16 << 20
DIM = BLOCK * BLOCK - 1


def report_path(prefix: str = "forward_mixing_certificate") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_entries(entries: int, bytes_per_entry: int = 512, label: str = "entries"):
    if type(entries) is not int or entries < 0:
        raise ValueError("entry count must be a nonnegative integer")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError(f"{label} exceeds 16 MiB preallocation budget")


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def interval_contains(value, target, binary_fraction):
    lo = binary_fraction(value.lower())
    hi = binary_fraction(value.upper())
    return lo <= target <= hi


def adjoint(matrix):
    return matrix.conjugate().transpose()


def trace(matrix, acb_zero):
    result = acb_zero(0)
    for index in range(matrix.nrows()):
        result += matrix[index, index]
    return result


def identity(size, acb_mat, acb_one):
    result = acb_mat(size, size)
    for index in range(size):
        result[index, index] = acb_one(1)
    return result


def basis_factory(acb, arb, acb_mat):
    """Return eight orthonormal Herm_0(C^3) Acb basis matrices."""
    inv_sqrt2 = arb(1) / arb(2).sqrt()
    inv_sqrt6 = arb(1) / arb(6).sqrt()
    result = []

    diagonal = acb_mat(BLOCK, BLOCK)
    diagonal[0, 0] = inv_sqrt2
    diagonal[1, 1] = -inv_sqrt2
    result.append(diagonal)

    diagonal = acb_mat(BLOCK, BLOCK)
    diagonal[0, 0] = inv_sqrt6
    diagonal[1, 1] = inv_sqrt6
    diagonal[2, 2] = -2 * inv_sqrt6
    result.append(diagonal)

    for row in range(BLOCK):
        for column in range(row + 1, BLOCK):
            symmetric = acb_mat(BLOCK, BLOCK)
            symmetric[row, column] = inv_sqrt2
            symmetric[column, row] = inv_sqrt2
            result.append(symmetric)

            antisymmetric = acb_mat(BLOCK, BLOCK)
            antisymmetric[row, column] = acb(0, -inv_sqrt2)
            antisymmetric[column, row] = acb(0, inv_sqrt2)
            result.append(antisymmetric)
    if len(result) != DIM:
        raise AssertionError("Hermitian traceless basis dimension mismatch")
    return tuple(result)


def basis_audit(basis, binary_fraction, acb_zero):
    failures = []
    max_real_radius = Fraction(0)
    max_imag_radius = Fraction(0)
    for i, left in enumerate(basis):
        trace_value = trace(left, acb_zero)
        if not interval_contains(trace_value.real, 0, binary_fraction) \
                or not interval_contains(trace_value.imag, 0, binary_fraction):
            failures.append(("trace", i))
        for j, right in enumerate(basis):
            value = trace(adjoint(left) * right, acb_zero)
            max_real_radius = max(max_real_radius,
                                  binary_fraction(value.real.rad().upper()))
            max_imag_radius = max(max_imag_radius,
                                  binary_fraction(value.imag.rad().upper()))
            target = 1 if i == j else 0
            if not interval_contains(value.real, target, binary_fraction) \
                    or not interval_contains(value.imag, 0, binary_fraction):
                failures.append(("gram", i, j))
    return dict(ok=not failures, failure_count=len(failures),
                failures=failures[:12], max_real_radius=str(max_real_radius),
                max_imag_radius=str(max_imag_radius))


def matrix_unitarity_audit(branches, binary_fraction, acb_zero):
    failures = []
    for branch_index, pair in enumerate(branches):
        for member_index, branch in enumerate(pair):
            product = adjoint(branch) * branch
            for row in range(BLOCK):
                for column in range(BLOCK):
                    target = 1 if row == column else 0
                    value = product[row, column]
                    if not interval_contains(value.real, target, binary_fraction) \
                            or not interval_contains(value.imag, 0, binary_fraction):
                        failures.append((branch_index, member_index, row, column))
    return dict(ok=not failures, failure_count=len(failures), failures=failures[:12])


def channel_matrix(branches, basis, acb_zero, acb_mat, binary_fraction,
                   density_forward_step):
    b0, b1 = branches
    result = acb_mat(DIM, DIM)
    max_imaginary = Fraction(0)
    max_trace = Fraction(0)
    for column, source in enumerate(basis):
        image = density_forward_step(
            b0, b1, source, multiply=lambda left, right: left * right,
            adjoint=adjoint)
        image_trace = trace(image, acb_zero)
        if (not interval_contains(image_trace.real,0,binary_fraction)
                or not interval_contains(image_trace.imag,0,binary_fraction)
                or not matrix_difference_contains_zero(image,adjoint(image),binary_fraction)):
            raise AssertionError("channel failed traceless Hermitian invariant-subspace audit")
        max_trace = max(max_trace, abs(binary_fraction(image_trace.real.rad().upper())))
        for row, target in enumerate(basis):
            coefficient = trace(adjoint(target) * image, acb_zero)
            if not interval_contains(coefficient.imag,0,binary_fraction):
                raise AssertionError("Hermitian channel coordinate does not enclose a real number")
            max_imaginary = max(max_imaginary,
                                abs(binary_fraction(coefficient.imag.rad().upper())))
            result[row, column] = coefficient
    return result, dict(max_imaginary_radius=str(max_imaginary),
                        max_trace_radius=str(max_trace))


def frobenius_upper(matrix, binary_fraction, arb):
    squared = arb(0)
    for row in range(matrix.nrows()):
        for column in range(matrix.ncols()):
            squared += matrix[row, column].abs_upper() ** 2
    return binary_fraction(squared.sqrt().upper())


def matrix_difference_upper(left, right, binary_fraction):
    maximum = Fraction(0)
    for row in range(left.nrows()):
        for column in range(left.ncols()):
            maximum = max(maximum,
                          binary_fraction((left[row, column] - right[row, column]).abs_upper()))
    return maximum


def matrix_difference_contains_zero(left, right, binary_fraction):
    for row in range(left.nrows()):
        for column in range(left.ncols()):
            difference = left[row, column] - right[row, column]
            if (not interval_contains(difference.real, 0, binary_fraction)
                    or not interval_contains(difference.imag, 0, binary_fraction)):
                return False
    return True


def phase_identity_audit(branches, offsets, period, block, alpha,
                         binary_fraction, rational_phase):
    """Check the exact scalar identity behind channel period two/six.

    In a b=3 sector, T^3 is a scalar.  Thus the powers separated by
    3*2**i differ by a scalar, and their conjugation channels are identical.
    The divisibility check is exact integer algebra; the matrix residual is
    only used to certify that the enclosed branch matrices contain that
    identity.
    """
    rows = []
    for separation in (2, 6):
        algebra_ok = True
        residuals = []
        contains_zero_rows = []
        for offset in range(offsets):
            left_index = offset
            right_index = offset + separation
            difference = (1 << right_index) - (1 << left_index)
            divisible = (difference % block) == 0
            algebra_ok &= divisible
            if not divisible:
                residuals.append(None)
                contains_zero_rows.append(False)
                continue
            scalar = rational_phase(alpha * (difference // block), period // block)
            left = branches[left_index]
            right = branches[right_index]
            left_shift = adjoint(left[0]) * left[1]
            right_shift = adjoint(right[0]) * right[1]
            contains_zero = matrix_difference_contains_zero(
                right_shift, scalar * left_shift, binary_fraction)
            contains_zero_rows.append(contains_zero)
            residuals.append(matrix_difference_upper(
                right_shift, scalar * left_shift, binary_fraction))
        rows.append(dict(separation=separation, algebra_ok=algebra_ok,
                         residuals=residuals,
                         residual_contains_zero=all(contains_zero_rows),
                         max_residual=(max(residuals) if residuals else None)))
    return rows


def mode_projector(shift, alpha, acb, acb_mat, acb_zero, rational_phase):
    """Exact spectral projector polynomial for the no-background shift."""
    unit = identity(BLOCK, acb_mat, acb)
    omega = rational_phase(1, 3)
    lam = rational_phase(alpha, PERIOD)
    lam_inv = 1 / lam
    projector0 = (unit + (lam_inv * shift)
                  + (lam_inv ** 2) * (shift * shift)) / 3
    mode1_inv = 1 / (lam * omega)
    projector1 = (unit + (mode1_inv * shift)
                  + ((mode1_inv ** 2) * (shift * shift))) / 3
    return projector0 - projector1


def conserved_mode_audit(branches, alpha, basis, binary_fraction,
                         acb, acb_mat, acb_zero, rational_phase, arb,
                         density_forward_step):
    shift = branches[0][1]
    mode = mode_projector(shift, alpha, acb, acb_mat, acb_zero, rational_phase)
    current = mode
    for pair in branches[:WINDOW]:
        current = density_forward_step(
            pair[0], pair[1], current, multiply=lambda left, right: left * right,
            adjoint=adjoint)
    difference = current - mode
    difference_upper = frobenius_upper(difference, binary_fraction, arb)
    difference_contains_zero = matrix_difference_contains_zero(current,mode,binary_fraction)
    lower_sq = arb(0)
    for row in range(BLOCK):
        for column in range(BLOCK):
            lower_sq += current[row, column].abs_lower() ** 2
    norm_lower = binary_fraction(lower_sq.sqrt().lower())
    trace_mode = trace(mode, acb_zero)
    hermitian_failures = []
    hermitian_difference = mode - adjoint(mode)
    for row in range(BLOCK):
        for column in range(BLOCK):
            if not interval_contains(hermitian_difference[row, column].real, 0, binary_fraction) \
                    or not interval_contains(hermitian_difference[row, column].imag, 0, binary_fraction):
                hermitian_failures.append((row, column))
    return dict(channel_mode_difference_upper=str(difference_upper),
                channel_mode_difference_contains_zero=difference_contains_zero,
                conserved_mode_norm_lower=str(norm_lower),
                trace_contains_zero=interval_contains(trace_mode.real, 0, binary_fraction)
                and interval_contains(trace_mode.imag, 0, binary_fraction),
                hermitian=not hermitian_failures,
                hermitian_failures=hermitian_failures[:12])


def main() -> None:
    exp = Experiment("forward_mixing_certificate", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "all 18 repeated-Rx channel windows have certified Frobenius upper < 4/5")
    exp.predict("P2", "Acb basis, branch unitarity, and exact period-6 channel identities pass")
    exp.predict("P3", "the no-background spectral traceless mode survives the 12-channel product")
    exp.must_fail("C1", "a mixed fixed point alone implies traceless contraction")
    started = time.perf_counter()
    report = {"status": "PASS", "precisions": {}, "controls": {}}
    try:
        from flint import acb, acb_mat, arb, ctx
        from lab.verified_finite_work import VerifiedFiniteWork
        from lab.semiclassical import _density_forward_step
        from lab.verified_prefix import VerifiedReflectionCircuit, _rational_phase, binary_fraction

        # Count the simultaneously live branch/state/channel/basis/report
        # payloads before constructing any Acb matrices.  This is deliberately
        # conservative and includes both the repeated and null builds.
        repeated_branch_entries = 3 * WIDTH * 2 * BLOCK * BLOCK
        repeated_state_entries = 3 * (WIDTH + 1) * BLOCK * BLOCK
        repeated_channel_entries = 3 * WIDTH * DIM * DIM
        repeated_product_entries = 3 * OFFSETS * DIM * DIM
        null_branch_entries = 3 * WINDOW * 2 * BLOCK * BLOCK
        null_state_entries = 3 * (WINDOW + 1) * BLOCK * BLOCK
        basis_entries = DIM * BLOCK * BLOCK
        report_entries = 3 * (OFFSETS * 2 + 40)
        allocation_entry_bound = (repeated_branch_entries + repeated_state_entries
            + repeated_channel_entries + repeated_product_entries
            + null_branch_entries + null_state_entries + basis_entries
            + report_entries + 2048)
        guard_entries(allocation_entry_bound,
                      label="aggregate branch/state/channel/basis/report objects")
        backgrounds = {s: ("x", Fraction(1, 7)) for s in range(1, WIDTH + 1)}
        repeated = VerifiedReflectionCircuit(PERIOD, WIDTH, backgrounds, {},
                                              block_size=BLOCK)
        null = VerifiedReflectionCircuit(PERIOD, WINDOW, {}, {}, block_size=BLOCK)
        p1 = True
        p2 = True
        p3 = True
        for precision in PRECISIONS:
            with ctx.workprec(precision):
                basis = basis_factory(acb, arb, acb_mat)
                basis_row = basis_audit(basis, binary_fraction, acb)
                repeated_worker = VerifiedFiniteWork(repeated, 0)
                branches_by_sector = {}
                for sector in range(repeated.sectors):
                    branches, _states, _final = repeated_worker._build(sector)
                    branches_by_sector[sector] = branches
                window_rows = []
                period_diffs = []
                unitary_rows = []
                phase_identity_rows = []
                period2_all_diffs = []
                for sector in range(repeated.sectors):
                    branches = branches_by_sector[sector]
                    channels = []
                    channel_audits = []
                    for pair in branches:
                        channel, audit = channel_matrix(
                            pair, basis, acb, acb_mat, binary_fraction,
                            _density_forward_step)
                        channels.append(channel)
                        channel_audits.append(audit)
                    for offset in range(OFFSETS):
                        product = identity(DIM, acb_mat, acb)
                        for channel in channels[offset:offset + WINDOW]:
                            product = channel * product
                        upper = frobenius_upper(product, binary_fraction, arb)
                        row = dict(sector=sector, offset=offset,
                                   frobenius_upper=str(upper),
                                   below_four_fifths=upper < THRESHOLD,
                                   channel_imaginary_radius=max(
                                       Fraction(a["max_imaginary_radius"])
                                       for a in channel_audits[offset:offset + WINDOW]),
                                   channel_trace_radius=max(
                                       Fraction(a["max_trace_radius"])
                                       for a in channel_audits[offset:offset + WINDOW]))
                        window_rows.append(row)
                        p1 &= row["below_four_fifths"]
                    for offset in range(OFFSETS):
                        period_diffs.append(matrix_difference_upper(
                            channels[offset], channels[offset + 6], binary_fraction))
                    unitary_rows.append(matrix_unitarity_audit(
                        branches, binary_fraction, acb))
                    phase_rows = phase_identity_audit(
                        branches, OFFSETS, PERIOD, BLOCK, sector,
                        binary_fraction, _rational_phase)
                    # Channel period two is the exact scalar-conjugation
                    # identity; period six is retained as the raw r=9 window
                    # schedule requested by the frozen experiment.
                    period2_diffs = [matrix_difference_upper(
                        channels[offset], channels[offset + 2], binary_fraction)
                        for offset in range(WIDTH - 2)]
                    phase_identity_rows.append(dict(sector=sector,
                                                    rows=phase_rows))
                    period2_all_diffs.extend(period2_diffs)
                period_bound = Fraction(1, 1 << (precision // 2))
                period_ok = all(value < period_bound for value in period_diffs)
                phase_ok = all(row["rows"][j]["algebra_ok"]
                               and row["rows"][j]["max_residual"] is not None
                               and row["rows"][j]["residual_contains_zero"]
                               and row["rows"][j]["max_residual"] < period_bound
                               for row in phase_identity_rows
                               for j in range(2))
                period2_ok = all(value < period_bound for value in period2_all_diffs)
                unitary_ok = all(row["ok"] for row in unitary_rows)
                p2 &= (basis_row["ok"] and period_ok and period2_ok
                       and phase_ok and unitary_ok)
                report["precisions"][str(precision)] = dict(
                    basis=basis_row, windows=window_rows,
                    period6_channel_max_difference=str(max(period_diffs)),
                    period6_threshold=str(period_bound), period6_ok=period_ok,
                    period2_channel_max_difference=str(max(period2_all_diffs)),
                    period2_ok=period2_ok,
                    phase_identity=phase_identity_rows,
                    phase_identity_ok=phase_ok,
                    branch_unitarity=unitary_rows, basis_and_unitarity_ok=unitary_ok,
                    all_windows_below_four_fifths=all(row["below_four_fifths"]
                                                     for row in window_rows),
                    allocation_entry_bound=allocation_entry_bound)

                null_worker = VerifiedFiniteWork(null, 0)
                null_modes = {}
                for sector in range(null.sectors):
                    null_branches, _states, _final = null_worker._build(sector)
                    null_modes[str(sector)] = conserved_mode_audit(
                        null_branches, sector, basis, binary_fraction, acb, acb_mat,
                        acb, _rational_phase, arb, _density_forward_step)
                report["precisions"][str(precision)]["null_modes"] = null_modes
                p3 &= all(item["channel_mode_difference_upper"]
                          and Fraction(item["channel_mode_difference_upper"]) < period_bound
                          and item["channel_mode_difference_contains_zero"]
                          and item["conserved_mode_norm_lower"]
                          and Fraction(item["conserved_mode_norm_lower"]) > Fraction(1, 2)
                          and item["trace_contains_zero"] and item["hermitian"]
                          for item in null_modes.values())

        report["controls"] = dict(
            mixed_fixed_point_is_not_contraction=bool(p3),
            null_mode_norm_lower={precision: {
                sector: report["precisions"][str(precision)]["null_modes"][sector]["conserved_mode_norm_lower"]
                for sector in report["precisions"][str(precision)]["null_modes"]}
                for precision in PRECISIONS},
            consequence_trace_norm_note="Only ||Delta||_1 <= sqrt(3)||Delta||_F is available; no trace contraction is claimed.")
        exp.check("P1", p1, "all 3 sectors x 6 offsets have outward Frobenius upper < 4/5")
        exp.check("P2", p2, "Acb basis/branch unitarity/period-6 channel checks")
        exp.check("P3", p3, "explicit no-background conserved traceless spectral mode")
        exp.fail_check("C1", bool(p3),
                       "mixed fixed point control: null spectral mode remains noncontracting")
        report["status"] = "PASS" if p1 and p2 and p3 else "FAIL"
        path = report_path()
        if not exp.finish(report_path=path, rows=[
                dict(series="precision", precision=precision,
                     **json_safe(report["precisions"][str(precision)]))
                for precision in PRECISIONS],
                metadata=dict(period=PERIOD, block_size=BLOCK, window=WINDOW,
                              offsets=OFFSETS, precisions=PRECISIONS,
                              threshold=str(THRESHOLD), norm="Hilbert-Schmidt Frobenius upper",
                              trace_conversion="||Delta||_1 <= sqrt(3)||Delta||_F",
                              allocation_bound_bytes=MAX_BYTES,
                              allocation_entry_bound=allocation_entry_bound,
                              arithmetic="Acb/Arb outward enclosures")):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = report_path("forward_mixing_certificate_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
