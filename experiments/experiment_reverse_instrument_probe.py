"""Bounded reverse-instrument identity audit for the b=3 finite-work law.

For fixed unitary branch pairs B_i0,B_i1, let

    K_i,z = (B_i0 + (-1)^z exp(-i*pi*prefix/2^(j+1)) B_i1)/2.

The reverse instrument uses A_i,z=K_i,z^dag.  Since
sum_z K_i,z K_i,z^dag=I, starting sigma=I/3 gives a normalized proposal
q_y=Tr(sigma_y).  With the actual initial rho_0 (including the fixed W_0),
the accepted submeasure is Tr(rho_0 sigma_y)=p_y/3, so the expected number
of proposals is b=3.  This probe checks that identity for all four route
histories and three initial sectors of the tiny r=9,b=3,t=4 fixture.

The CONDITIONAL target p_y is BLOCK times the existing verified joint
``target_interval`` helper; an independent full-r float ``direct_joint`` law
is only a diagnostic.  No float reverse transition is treated as certified,
and no division by a zero effect is performed. Raw effects below are
unnormalized; they diagnose the exact rejection identity, not an implemented
finite-bit proposal/acceptance algorithm. Interval overlaps are numerical
checks of the independently proved identity, not proofs of equality.
"""
from __future__ import annotations

import json
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

from lab import Experiment


PERIOD = 9
BLOCK = 3
WIDTH = 4
HISTORY_COUNT = 4
INITIAL_SECTORS = 3
OUTPUTS = 1 << WIDTH
PRECISIONS = (192, 256)
MAX_BYTES = 16 << 20


def report_path(prefix="reverse_instrument_probe"):
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
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def guard_entries(entries, bytes_per_entry=512):
    if type(entries) is not int or entries < 0:
        raise ValueError("entry count must be a nonnegative integer")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("reverse-instrument retained payload exceeds 16 MiB")


def adjoint(matrix):
    return matrix.conjugate().transpose()


def trace(matrix, acb):
    result = acb(0)
    for index in range(matrix.nrows()):
        result += matrix[index, index]
    return result


def interval(value, binary_fraction):
    return (binary_fraction(value.lower()), binary_fraction(value.upper()))


def real_interval(value, binary_fraction):
    lo, hi = interval(value.real, binary_fraction)
    if hi < 0:
        raise AssertionError("real mass enclosure is strictly negative")
    return max(Fraction(0), lo), max(Fraction(0), hi)


def contains_zero(value, binary_fraction):
    lo, hi = interval(value, binary_fraction)
    return lo <= 0 <= hi


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2),
                    np.exp(1j * theta / 2)]).astype(complex)


def fixture():
    """W0 plus the original W1/W3 and both supplied route declarations."""
    backgrounds = {0: ("x", Fraction(1, 4)),
                   1: ("x", Fraction(1, 7)),
                   3: ("z", Fraction(1, 5))}
    all_reflections = {2: (0, Fraction(1, 5)),
                       3: (1, Fraction(1, 5))}
    from lab.verified_prefix import VerifiedReflectionCircuit
    return VerifiedReflectionCircuit(PERIOD, WIDTH, backgrounds, all_reflections,
                                     block_size=BLOCK)


def float_reference(mask):
    from lab.coherent_routes import CoherentReflectionCircuit
    all_reflections = {2: (0, np.pi), 3: (1, np.pi)}
    selected = {s: value for bit, (s, value) in enumerate(all_reflections.items())
                if mask & (1 << bit)}
    w0 = np.eye(BLOCK, dtype=complex)
    w0[:2, :2] = rx(np.pi / 4)
    w1 = np.eye(BLOCK, dtype=complex)
    w1[:2, :2] = rx(np.pi / 7)
    w3 = np.eye(BLOCK, dtype=complex)
    w3[:2, :2] = rz(np.pi / 5)
    return CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH,
        {0: w0, 1: w1, 3: w3},
        selected)


def reverse_leaf(branches, rho0, output, acb, acb_mat, rational_phase,
                 binary_fraction, qft_branch_operators):
    """Return q_y and accepted submeasure intervals for one output label."""
    sigma = acb_mat([[int(i == j) for j in range(BLOCK)]
                     for i in range(BLOCK)]) / BLOCK
    prefix = 0
    for depth in range(WIDTH):
        index = WIDTH - 1 - depth
        phase = rational_phase(-prefix, 1 << (depth + 1))
        K0, K1 = qft_branch_operators(
            branches[index][0], branches[index][1], phase)
        # The terminating inverse-QFT measures controls in descending
        # arithmetic order, while its returned output integer records those
        # measured bits in ascending (measurement-prefix) positions.
        bit = (output >> depth) & 1
        K = K1 if bit else K0
        sigma = adjoint(K) * sigma * K
        # Feedback prefix and returned output both use measurement order.
        prefix |= bit << depth
    q_value = trace(sigma, acb)
    accepted_value = trace(rho0 * sigma, acb)
    mixed = acb_mat([[int(i == j) for j in range(BLOCK)]
                     for i in range(BLOCK)]) / BLOCK
    null_value = trace(mixed * sigma, acb)
    q_lo, q_hi = real_interval(q_value, binary_fraction)
    a_lo, a_hi = real_interval(accepted_value, binary_fraction)
    n_lo, n_hi = real_interval(null_value, binary_fraction)
    return dict(q=(q_lo, q_hi), accepted=(a_lo, a_hi),
                null_accepted=(n_lo, n_hi),
                q_imag_zero=contains_zero(q_value.imag, binary_fraction),
                accepted_imag_zero=contains_zero(accepted_value.imag, binary_fraction),
                null_imag_zero=contains_zero(null_value.imag, binary_fraction))


def overlap(left, right):
    return max(left[0], right[0]) <= min(left[1], right[1])


def midpoint(value):
    return (value[0] + value[1]) / 2


def tv_lower_intervals(left, right):
    """Outward lower bound: charge BOTH laws' enclosure radii."""
    total = Fraction(0)
    for key, (lo, hi) in right.items():
        other_lo, other_hi = left[key]
        total += max(Fraction(0), lo - other_hi, other_lo - hi)
    return total / 2


def main():
    exp = Experiment("reverse_instrument_probe", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "reverse proposal laws normalize and accepted submeasures agree with target intervals at P192/P256")
    exp.predict("P2", "the expected proposal count is b=3 and independent full-r references match verified target intervals")
    exp.must_fail("C1", "skipping terminal acceptance leaves the reverse proposal equal to the target law")
    exp.must_fail("C2", "a normalized maximally mixed initial state has nonconstant terminal acceptance")
    report = {"status": "PASS", "rows": [], "controls": {}}
    started = time.perf_counter()
    try:
        from flint import acb, acb_mat, ctx
        from lab.semiclassical import _qft_branch_operators
        from lab.verified_finite_work import VerifiedFiniteWork
        from lab.verified_prefix import _rational_phase, binary_fraction
        from experiments.experiment_verified_rejection import target_interval
        from experiments.experiment_coherent_route_sampling import direct_joint

        # Both precision reports are retained, along with scalar strings,
        # serialization copies and the live tiny direct-reference arrays.
        # This conservative structural preflight is NOT measured native RSS.
        retained = HISTORY_COUNT * INITIAL_SECTORS * len(PRECISIONS) * OUTPUTS * 12
        guard_entries(3 * retained + 4096)
        p1 = True
        p2 = True
        proposal_tv_lowers = []
        null_constant_rows = []
        for mask in range(HISTORY_COUNT):
            circuit = fixture()
            comparator = float_reference(mask)
            independent = direct_joint(comparator)
            worker = VerifiedFiniteWork(circuit, mask)
            from lab.verified_rejection import VerifiedRejectionSampler
            target_component = VerifiedRejectionSampler(circuit).component(mask)
            for precision in PRECISIONS:
                with ctx.workprec(precision):
                    row_set = []
                    for initial in range(INITIAL_SECTORS):
                        branches, states, final = worker._build(initial)
                        rho0 = states[0]
                        target = {y: target_interval(target_component, final, y, precision)
                                  for y in range(OUTPUTS)}
                        reverse = {y: reverse_leaf(
                            branches, rho0, y, acb, acb_mat, _rational_phase,
                            binary_fraction, _qft_branch_operators)
                                   for y in range(OUTPUTS)}
                        q_intervals = {y: reverse[y]["q"] for y in range(OUTPUTS)}
                        accepted_intervals = {y: reverse[y]["accepted"]
                                              for y in range(OUTPUTS)}
                        null_intervals = {y: reverse[y]["null_accepted"]
                                          for y in range(OUTPUTS)}
                        q_mass = (sum(v[0] for v in q_intervals.values()),
                                  sum(v[1] for v in q_intervals.values()))
                        accepted_mass = (sum(v[0] for v in accepted_intervals.values()),
                                         sum(v[1] for v in accepted_intervals.values()))
                        normalized_target = {y: (3 * target[y][0], 3 * target[y][1])
                                             for y in range(OUTPUTS)}
                        normalized_accepted = {y: (3 * accepted_intervals[y][0],
                                                   3 * accepted_intervals[y][1])
                                               for y in range(OUTPUTS)}
                        enclosure_overlap = all(
                            overlap(accepted_intervals[y], target[y])
                            and overlap(normalized_accepted[y], normalized_target[y])
                            for y in range(OUTPUTS))
                        proposal_imag_ok = all(reverse[y]["q_imag_zero"]
                                               and reverse[y]["accepted_imag_zero"]
                                               and reverse[y]["null_imag_zero"]
                                               for y in range(OUTPUTS))
                        target_mass = (sum(v[0] for v in target.values()),
                                       sum(v[1] for v in target.values()))
                        full_ref_error = max(
                            abs(float(midpoint(target[y])) - independent[final, y])
                            for y in range(OUTPUTS))
                        max_interval_width = max(
                            hi-lo for law in (q_intervals, accepted_intervals,
                                              null_intervals, target)
                            for lo, hi in law.values())
                        proposal_tv_lower = tv_lower_intervals(
                            q_intervals, normalized_target)
                        proposal_tv_lowers.append(proposal_tv_lower)
                        # Null acceptance with rho0=I/3: compare the raw
                        # numerator Tr((I/3)sigma_y) to q_y/3, without a
                        # ratio and hence without dividing a zero effect.
                        null_constant = all(
                            overlap(null_intervals[y],
                                    (q_intervals[y][0] / BLOCK,
                                     q_intervals[y][1] / BLOCK))
                            for y in range(OUTPUTS))
                        null_constant_rows.append(null_constant)
                        p1 &= (q_mass[0] <= 1 <= q_mass[1]
                               and accepted_mass[0] <= Fraction(1, BLOCK) <= accepted_mass[1]
                               and target_mass[0] <= Fraction(1, BLOCK) <= target_mass[1]
                               and enclosure_overlap and proposal_imag_ok
                               and max_interval_width < Fraction(1, 1 << 120))
                        p2 &= full_ref_error < 2e-14
                        row_set.append(dict(initial=initial, final=final,
                                            q_mass=[str(x) for x in q_mass],
                                            accepted_mass=[str(x) for x in accepted_mass],
                                            target_mass=[str(x) for x in target_mass],
                                            enclosure_overlap=enclosure_overlap,
                                            max_interval_width=str(max_interval_width),
                                            full_reference_max_error=full_ref_error,
                                            proposal_tv_lower=str(proposal_tv_lower),
                                            outputs=[dict(y=y, q=list(map(str, q_intervals[y])),
                                                          accepted=list(map(str, accepted_intervals[y])),
                                                          target=list(map(str, target[y])),
                                                          normalized_accepted=list(map(str, normalized_accepted[y])),
                                                          normalized_target=list(map(str, normalized_target[y])))
                                                      for y in range(OUTPUTS)]))
                    report["rows"].append(dict(mask=mask, precision=precision,
                                               rows=row_set))
        report["controls"] = dict(
            skip_terminal_acceptance_tv_lower=str(max(proposal_tv_lowers)),
            skip_terminal_acceptance_changes_target=max(proposal_tv_lowers) > Fraction(1, 1000),
            maxmixed_null_constant=all(null_constant_rows),
            expected_attempts=BLOCK,
            target_scope="joint final-sector/output target mass; accepted mass is 1/3 per fixed initial")
        exp.check("P1", p1, "all 24 initial/history/precision rows normalize and overlap target intervals")
        exp.check("P2", p2, "full-r direct_joint midpoint diagnostic error below 2e-14")
        exp.fail_check("C1", max(proposal_tv_lowers) > Fraction(1, 1000),
                       f"proposal-vs-normalized-target TV lower={float(max(proposal_tv_lowers)):.6g}")
        exp.fail_check("C2", all(null_constant_rows),
                       "maximally mixed normalized initial state gives constant acceptance 1/3")
        report["status"] = "PASS" if p1 and p2 and max(proposal_tv_lowers) > Fraction(1, 1000) and all(null_constant_rows) else "FAIL"
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]),
                          metadata=dict(period=PERIOD, block_size=BLOCK, width=WIDTH,
                                        histories=HISTORY_COUNT, initial_sectors=INITIAL_SECTORS,
                                        precisions=PRECISIONS, allocation_bound_bytes=MAX_BYTES,
                                        arithmetic="Acb/Arb outward intervals; overlap checks are not equality proofs",
                                        controls=json_safe(report["controls"]),
                                        retained_entry_allowance=3 * retained + 4096,
                                        accepted_identity="q_y*a_y=Tr(rho0 sigma_y)=p_y/3",
                                        float_reference="direct_joint diagnostic only",
                                        elapsed_seconds=time.perf_counter() - started)):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = report_path("reverse_instrument_probe_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
