"""Tiny complete-law verifier for the finite-bit reverse rejection sampler.

This is the first finite-bit audit of ``VerifiedReverseWork``.  For the
fixed r=9,b=3,t=4 fixture, every deterministic route history, initial coarse
sector, boundary label j and output y is enumerated.  The API returns a
proposal probability conditional on j and a finite acceptance probability;
the experiment averages both over the uniform j label, checks the accepted
submeasure, and normalizes it by its measured accepted mass before comparing
with the existing outward target intervals.

The independent full-r ``direct_joint`` calculation is a float diagnostic
only.  No histogram, generic propagator, forward checkpoint or acceptance
cap is used.  The width-zero identity fixture explicitly exercises exact
zero acceptance NUMERATORS (its proposal denominator is one).
"""
from __future__ import annotations

import json
import random
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
HISTORIES = 4
INITIALS = 3
OUTPUTS = 1 << WIDTH
TARGETS = (Fraction(1, 1000), Fraction(1, 1_000_000))
PRECISIONS = (192, 256)
MAX_BYTES = 16 << 20


def report_path(prefix="verified_reverse_work"):
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
        raise MemoryError("verified reverse retained payload exceeds 16 MiB")


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2),
                    np.exp(1j * theta / 2)]).astype(complex)


def base_fixture():
    from lab.verified_prefix import VerifiedReflectionCircuit
    return VerifiedReflectionCircuit(
        PERIOD, WIDTH,
        {0: ("x", Fraction(1, 4)),
         1: ("x", Fraction(1, 7)),
         3: ("z", Fraction(1, 5))},
        {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))},
        block_size=BLOCK)


def float_reference(mask):
    from lab.coherent_routes import CoherentReflectionCircuit
    selected = {2: (0, np.pi), 3: (1, np.pi)}
    selected = {s: value for bit, (s, value) in enumerate(selected.items())
                if mask & (1 << bit)}
    w0 = np.eye(BLOCK, dtype=complex)
    w0[:2, :2] = rx(np.pi / 4)
    w1 = np.eye(BLOCK, dtype=complex)
    w1[:2, :2] = rx(np.pi / 7)
    w3 = np.eye(BLOCK, dtype=complex)
    w3[:2, :2] = rz(np.pi / 5)
    return CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH, {0: w0, 1: w1, 3: w3}, selected)


def tv_from_intervals(law, intervals):
    total = Fraction(0)
    for key, (lo, hi) in intervals.items():
        value = law.get(key, Fraction(0))
        total += max(abs(value - lo), abs(value - hi))
    return total / 2


def tv_lower_from_intervals(law, intervals):
    total = Fraction(0)
    for key, (lo, hi) in intervals.items():
        value = law.get(key, Fraction(0))
        total += max(Fraction(0), lo - value, value - hi)
    return total / 2


def run_zero_width(target_tv):
    from lab.verified_prefix import VerifiedReflectionCircuit
    from lab.verified_reverse_work import VerifiedReverseWork
    circuit = VerifiedReflectionCircuit(PERIOD, 0, {}, {}, block_size=BLOCK)
    worker = VerifiedReverseWork(circuit, 0)
    rows = []
    for initial in range(INITIALS):
        for boundary in range(BLOCK):
            result = worker.attempt(initial, boundary, output=0,
                                    target_tv=target_tv)
            rows.append(dict(initial=initial, boundary=boundary,
                             proposal=str(result["proposal_path_probability"]),
                             acceptance=str(result["acceptance_probability"]),
                             final=result["final_coarse_sector"],
                             zero_fallbacks=result["zero_approximate_block_fallbacks"]))
    return rows


def main():
    exp = Experiment("verified_reverse_work", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "all finite-bit proposal and accepted laws normalize and meet their declared target budgets")
    exp.predict("P2", "accepted laws agree with target intervals and independent full-r diagnostics")
    exp.predict("P3", "zero-width exact-zero boundary labels are handled without division or deletion")
    exp.must_fail("C1", "skipping terminal acceptance leaves the proposal equal to the target law")
    report = {"status": "PASS", "rows": [], "controls": {}}
    started = time.perf_counter()
    try:
        from lab.verified_reverse_work import VerifiedReverseWork
        from lab.verified_rejection import VerifiedRejectionSampler
        from experiments.experiment_verified_rejection import target_interval
        from experiments.experiment_coherent_route_sampling import direct_joint

        # Retain summarized complete laws, not every full attempt dictionary.
        # Charge rows plus serialization and current small j-conditioned laws.
        retained = HISTORIES * INITIALS * len(TARGETS) * len(PRECISIONS) * OUTPUTS * 6
        guard_entries(3 * retained + 4096)
        circuit = base_fixture()
        p1 = True
        p2 = True
        skip_acceptance_lowers = []
        attempt_rows = 0
        for mask in range(HISTORIES):
            worker = VerifiedReverseWork(circuit, mask)
            component = VerifiedRejectionSampler(circuit).component(mask)
            independent = direct_joint(float_reference(mask))
            for target_tv in TARGETS:
                plan = worker.plan(target_tv)
                sample_rng = random.Random(28000 + mask * 17 + target_tv.denominator)
                sample = worker.sample(sample_rng, target_tv=target_tv)
                sample_ok = (sample["attempts"] >= 1
                             and sample["forward_steps"] == 0
                             and not sample["orbit_or_output_tables"]
                             and sample["attempt_cap"] is None
                             and sample["precision_cap"] is None)
                for precision in PRECISIONS:
                    row_set = []
                    for initial in range(INITIALS):
                        proposal = {y: Fraction(0) for y in range(OUTPUTS)}
                        accepted = {y: Fraction(0) for y in range(OUTPUTS)}
                        cost_sum = {key: 0 for key in (
                            "reverse_vector_matvecs", "replayed_vector_steps",
                            "branch_matrix_pair_constructions", "refinement_retries",
                            "acceptance_enclosure_evaluations")}
                        final = None
                        proposal_rows_ok = True
                        boundary_tv_bounds = []
                        for boundary in range(BLOCK):
                            boundary_proposal = {}
                            boundary_intervals = {}
                            for y in range(OUTPUTS):
                                result = worker.attempt(
                                    initial, boundary, output=y,
                                    target_tv=target_tv)
                                prop = Fraction(result["proposal_path_probability"])
                                acc = Fraction(result["acceptance_probability"])
                                proposal[y] += prop / BLOCK
                                accepted[y] += prop * acc / BLOCK
                                final = result["final_coarse_sector"]
                                for key in cost_sum:
                                    cost_sum[key] += result[key]
                                attempt_rows += 1
                                boundary_proposal[y] = prop
                                boundary_intervals[y] = result["acceptance_denominator_bounds"]
                                proposal_rows_ok &= (0 <= acc <= 1
                                    and result["stored_forward_matrix_count"] == 0
                                    and result["stored_branch_matrix_count"] == 0
                                    and result["acceptance_width_sum"] <= plan["acceptance_interval_width_tolerance"])
                            proposal_bound = tv_from_intervals(boundary_proposal,boundary_intervals)
                            boundary_tv_bounds.append(proposal_bound)
                            proposal_rows_ok &= (sum(boundary_proposal.values()) == 1
                                and proposal_bound <= plan["proposal_tv_upper_bound"])
                        if final is None:
                            raise AssertionError("empty finite reverse law")
                        proposal_mass = sum(proposal.values())
                        accepted_mass = sum(accepted.values())
                        target_intervals = {
                            y: target_interval(component, final, y, precision)
                            for y in range(OUTPUTS)}
                        conditional_target = {
                            y: (BLOCK * lo, BLOCK * hi)
                            for y, (lo, hi) in target_intervals.items()}
                        if accepted_mass <= 0:
                            raise AssertionError("accepted mass vanished")
                        normalized_accepted = {
                            y: accepted[y] / accepted_mass for y in range(OUTPUTS)}
                        normalized_intervals = {
                            y: (lo, hi) for y, (lo, hi) in conditional_target.items()}
                        tv = tv_from_intervals(normalized_accepted,
                                               normalized_intervals)
                        tv_lower = tv_lower_from_intervals(
                            proposal, normalized_intervals)
                        skip_acceptance_lowers.append(tv_lower)
                        target_mass = (sum(lo for lo, _ in target_intervals.values()),
                                       sum(hi for _, hi in target_intervals.values()))
                        full_error = max(
                            abs(float(sum(target_intervals[y]) / 2)
                                - independent[final, y])
                            for y in range(OUTPUTS))
                        success_lower = Fraction(plan["success_probability_lower_bound"])
                        expected_upper = Fraction(plan["expected_attempts_upper_bound"])
                        strict_success = (accepted_mass >= success_lower
                                          and accepted_mass > 0
                                          and Fraction(1, accepted_mass) <= expected_upper)
                        law_ok = (proposal_mass == 1
                                  and target_mass[0] <= Fraction(1, BLOCK) <= target_mass[1]
                                  and accepted_mass >= success_lower
                                  and tv <= Fraction(plan["total_tv_upper_bound"])
                                  and strict_success and sample_ok and proposal_rows_ok)
                        p1 &= law_ok
                        p2 &= (full_error < 2e-14
                               and tv <= Fraction(plan["total_tv_upper_bound"])
                               and target_mass[0] <= Fraction(1, BLOCK) <= target_mass[1])
                        row_set.append(dict(
                            initial=initial, final=final,
                            proposal_mass=str(proposal_mass),
                            accepted_mass=str(accepted_mass),
                            success_lower=str(success_lower),
                            expected_attempts_upper=str(expected_upper),
                            normalized_tv_upper=str(tv),
                            skip_acceptance_tv_lower=str(tv_lower),
                            boundary_proposal_tv_uppers=list(map(str,boundary_tv_bounds)),
                            plan=json_safe(plan), costs=cost_sum,
                            full_reference_max_error=full_error,
                            proposal=[str(proposal[y]) for y in range(OUTPUTS)],
                            accepted=[str(accepted[y]) for y in range(OUTPUTS)],
                            normalized_accepted=[str(normalized_accepted[y])
                                                 for y in range(OUTPUTS)],
                            conditional_target=[list(map(str, conditional_target[y]))
                                                for y in range(OUTPUTS)],
                            law_ok=law_ok))
                    report["rows"].append(dict(mask=mask, target_tv=str(target_tv),
                                               precision=precision, rows=row_set))
        zero_rows = run_zero_width(TARGETS[-1])
        zero_proposals = [Fraction(row["proposal"]) for row in zero_rows]
        zero_acceptances = [Fraction(row["acceptance"]) for row in zero_rows]
        zero_ok = (all(value == 1 for value in zero_proposals)
                   and all(Fraction(row["acceptance"]) == int(row["boundary"] == 0)
                           and row["final"] == row["initial"] for row in zero_rows))
        p3 = zero_ok
        report["controls"] = dict(
            skip_acceptance_tv_lower=str(max(skip_acceptance_lowers)),
            skip_acceptance_changes_target=max(skip_acceptance_lowers) > Fraction(1, 1000),
            zero_width=zero_rows,
            zero_width_exact_zero_boundary=zero_ok,
            forced_attempt_rows=attempt_rows,
            arithmetic="exact Fraction finite-bit laws with Acb/Arb target intervals")
        exp.check("P1", p1, "all 48 rows pass per-boundary proposal and accepted-law normalization/success/TV checks")
        exp.check("P2", p2, "target intervals and independent direct_joint diagnostics")
        exp.check("P3", p3, "zero-width zero-acceptance boundary labels retained")
        exp.fail_check("C1", max(skip_acceptance_lowers) > Fraction(1, 1000),
                       f"proposal-vs-target TV lower={float(max(skip_acceptance_lowers)):.6g}")
        report["status"] = "PASS" if p1 and p2 and p3 and max(skip_acceptance_lowers) > Fraction(1, 1000) else "FAIL"
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]),
                          metadata=dict(period=PERIOD, block_size=BLOCK, width=WIDTH,
                                        histories=HISTORIES, initials=INITIALS,
                                        targets=[str(x) for x in TARGETS],
                                        precisions=PRECISIONS,
                                        allocation_bound_bytes=MAX_BYTES,
                                        method="VerifiedReverseWork exact finite-bit enumeration",
                                        no_forward_checkpoints=True,
                                        controls=json_safe(report["controls"]),
                                        elapsed_seconds=time.perf_counter() - started)):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = report_path("verified_reverse_work_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
