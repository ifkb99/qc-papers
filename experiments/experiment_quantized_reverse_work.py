"""Bounded finite-law audit of the locally quantized reverse trajectory.

The quantized sampler keeps unnormalized Gaussian-integer reverse vectors and
rounds only their projective state after each selected child.  Its proposal is
not assumed to equal the ideal reverse proposal: every forced quantized
proposal mass is compared with the tight ``VerifiedReverseWork`` proposal at
the same boundary label.  Accepted laws are compared separately with the
existing outward ``target_interval`` reference and an independent full-r
float diagnostic.

The fixture is the fixed r=9,b=3,t=4 route family with W0 plus the original
W1/W3 blocks.  All four histories, three initial sectors, three boundary
labels, sixteen outputs, both requested target TVs and P192/P256 are covered.
Only summarized laws/costs are retained; no attempt dictionary or output
histogram is used.
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
REFERENCE_TARGET = Fraction(1, 10**12)
PRECISIONS = (192, 256)
MAX_BYTES = 16 << 20


def report_path(prefix="quantized_reverse_work"):
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
        raise MemoryError("quantized reverse retained payload exceeds 16 MiB")


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
    all_routes = {2: (0, np.pi), 3: (1, np.pi)}
    selected = {s: value for bit, (s, value) in enumerate(all_routes.items())
                if mask & (1 << bit)}
    w0 = np.eye(BLOCK, dtype=complex)
    w0[:2, :2] = rx(np.pi / 4)
    w1 = np.eye(BLOCK, dtype=complex)
    w1[:2, :2] = rx(np.pi / 7)
    w3 = np.eye(BLOCK, dtype=complex)
    w3[:2, :2] = rz(np.pi / 5)
    return CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH, {0: w0, 1: w1, 3: w3}, selected)


def tv(a, b):
    return sum(abs(a.get(k, 0) - b.get(k, 0))
               for k in set(a) | set(b)) / 2


def tv_from_intervals(law, intervals):
    return sum(max(abs(law.get(k, 0) - lo), abs(law.get(k, 0) - hi))
               for k, (lo, hi) in intervals.items()) / 2


def tv_lower_from_intervals(law, intervals):
    return sum(max(Fraction(0), lo - law.get(k, 0), law.get(k, 0) - hi)
               for k, (lo, hi) in intervals.items()) / 2


def quantization_controls():
    from lab.verified_quantized_reverse import (canonical_quantize,
                                                integer_matvec, integer_norm)
    zero, zero_flag = canonical_quantize(((0, 0), (0, 0), (0, 0)), 4)
    tied, tied_flag = canonical_quantize(((-3, 3), (3, -3), (0, 0)), 4)
    matrix = (((1, 0), (0, 1)), ((0, 1), (1, 0)))
    vector = ((2, 1), (-1, 3))
    product = integer_matvec(matrix, vector)
    expected = ((-1, 0), (-2, 5))
    return dict(zero_fallback=zero, zero_flag=zero_flag,
                tied=tied, tied_flag=tied_flag,
                tied_max_coordinate=max(abs(x) for z in tied for x in z),
                integer_norm=integer_norm(vector), matvec=product,
                helper_ok=(zero_flag and zero[0] == (16, 0)
                           and not tied_flag and tied[0] == (-16, 16)
                           and tied[1] == (16, -16)
                           and max(abs(x) for z in tied for x in z) == 16
                           and product == expected
                           and integer_norm(vector) == 15))


def main():
    exp = Experiment("quantized_reverse_work", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "quantized proposals stay within the local-vs-reference proposal budget and accepted laws meet plan bounds")
    exp.predict("P2", "accepted laws match target intervals and independent full-r diagnostics")
    exp.predict("P3", "fixed boundary-label proposals differ from the required uniform-j mixture")
    exp.predict("P4", "zero/tied/negative-scale Gaussian-integer helper cases satisfy canonical invariants")
    exp.must_fail("C1", "omitting terminal acceptance leaves the quantized proposal equal to the target law")
    report = {"status": "PASS", "rows": [], "controls": {}}
    started = time.perf_counter()
    try:
        from lab.verified_quantized_reverse import VerifiedQuantizedReverseWork
        from lab.verified_reverse_work import VerifiedReverseWork
        from lab.verified_rejection import VerifiedRejectionSampler
        from experiments.experiment_verified_rejection import target_interval
        from experiments.experiment_coherent_route_sampling import direct_joint

        guard_entries(HISTORIES * INITIALS * len(TARGETS) * len(PRECISIONS)
                      * BLOCK * OUTPUTS * 8 + 8192)
        helper = quantization_controls()
        circuit = base_fixture()
        p1 = True
        p2 = True
        p3 = True
        fixed_boundary_tvs = []
        omitted_acceptance_lowers = []
        for mask in range(HISTORIES):
            quantized = VerifiedQuantizedReverseWork(circuit, mask)
            reference = VerifiedReverseWork(circuit, mask)
            component = VerifiedRejectionSampler(circuit).component(mask)
            independent = direct_joint(float_reference(mask))
            for target_tv in TARGETS:
                qplan = quantized.plan(target_tv)
                refplan = reference.plan(REFERENCE_TARGET)
                sample = quantized.sample(
                    random.Random(37000 + mask * 31 + target_tv.denominator),
                    target_tv=target_tv)
                sample_ok = (sample["attempts"] >= 1
                             and sample["forward_steps"] == 0
                             and not sample["orbit_or_output_tables"]
                             and sample["attempt_cap"] is None
                             and sample["precision_cap"] is None
                             and sample["vector_compressions"]
                                 == sample["resampled_blocks"]
                             and sample["reverse_vector_matvecs"]
                                 == 2 * sample["resampled_blocks"]
                             and sample["stored_forward_matrix_count"] == 0
                             and sample["stored_branch_matrix_count"] == 0
                             and sample["max_state_coordinate_bits"]
                                 <= sample["state_coordinate_bits_upper_bound"]
                             and sample["max_initial_coordinate_bits"]
                                 <= sample["initial_coordinate_bits_upper_bound"]
                             and sample["max_child_coordinate_bits"]
                                 <= sample["child_coordinate_bits_upper_bound"]
                             and sample["working_scalar_units"]
                                 == "complex-coordinate slots; two integers per Gaussian-integer slot; excludes variable bit storage and backend scratch")
                for precision in PRECISIONS:
                    for initial in range(INITIALS):
                        q_by_boundary = []
                        reference_by_boundary = []
                        reference_denominator_intervals = []
                        accepted = {y: Fraction(0) for y in range(OUTPUTS)}
                        final = None
                        costs = {key: 0 for key in (
                            "operator_enclosure_evaluations",
                            "branch_matrix_pair_constructions",
                            "reverse_vector_matvecs", "vector_compressions",
                            "refinement_retries")}
                        for boundary in range(BLOCK):
                            qlaw = {y: Fraction(0) for y in range(OUTPUTS)}
                            rlaw = {y: Fraction(0) for y in range(OUTPUTS)}
                            dlaw = {}
                            for y in range(OUTPUTS):
                                qr = quantized.attempt(
                                    initial, boundary, output=y,
                                    target_tv=target_tv)
                                rr = reference.attempt(
                                    initial, boundary, output=y,
                                    target_tv=REFERENCE_TARGET)
                                qprob = Fraction(qr["proposal_path_probability"])
                                rprob = Fraction(rr["proposal_path_probability"])
                                dlaw[y] = tuple(Fraction(v) for v in
                                                rr["acceptance_denominator_bounds"])
                                qlaw[y] = qprob
                                rlaw[y] = rprob
                                accepted[y] += (
                                    qprob * Fraction(qr["acceptance_probability"])
                                    / BLOCK)
                                final = qr["final_coarse_sector"]
                                for key in costs:
                                    costs[key] += qr[key]
                                plan_bounds_ok = (
                                    qr["vector_compressions"] == WIDTH
                                    and qr["reverse_vector_matvecs"] == 2 * WIDTH
                                    and qr["replayed_vector_steps"] == 0
                                    and qr["forward_steps"] == 0
                                    and qr["stored_forward_matrix_count"] == 0
                                    and qr["stored_branch_matrix_count"] == 0
                                    and qr["max_state_coordinate_bits"]
                                        <= qr["state_coordinate_bits_upper_bound"]
                                    and qr["max_initial_coordinate_bits"]
                                        <= qr["initial_coordinate_bits_upper_bound"]
                                    and qr["max_child_coordinate_bits"]
                                        <= qr["child_coordinate_bits_upper_bound"]
                                    and not qr["orbit_or_output_tables"])
                                p1 &= plan_bounds_ok
                            q_by_boundary.append(qlaw)
                            reference_by_boundary.append(rlaw)
                            reference_denominator_intervals.append(dlaw)
                            denominator_tv = sum(
                                max(abs(qlaw[y] - dlaw[y][0]),
                                    abs(qlaw[y] - dlaw[y][1]))
                                for y in range(OUTPUTS)) / 2
                            p1 &= (sum(qlaw.values()) == 1
                                   and denominator_tv
                                       <= Fraction(qplan["proposal_tv_upper_bound"]))
                        proposal = {y: sum(q[y] for q in q_by_boundary) / BLOCK
                                    for y in range(OUTPUTS)}
                        reference_proposal = {
                            y: sum(q[y] for q in reference_by_boundary) / BLOCK
                            for y in range(OUTPUTS)}
                        proposal_ref_tv = tv(proposal, reference_proposal)
                        qbound = Fraction(qplan["proposal_tv_upper_bound"])
                        rbound = Fraction(refplan["proposal_tv_upper_bound"])
                        proposal_budget = qbound + rbound
                        denominator_tvs = [
                            sum(max(abs(q_by_boundary[j][y]
                                        - reference_denominator_intervals[j][y][0]),
                                    abs(q_by_boundary[j][y]
                                        - reference_denominator_intervals[j][y][1]))
                                for y in range(OUTPUTS)) / 2
                            for j in range(BLOCK)]
                        fixed_tv = tv(q_by_boundary[0], proposal)
                        fixed_boundary_tvs.append(fixed_tv)
                        target_intervals = {
                            y: target_interval(component, final, y, precision)
                            for y in range(OUTPUTS)}
                        conditional_target = {
                            y: (BLOCK * lo, BLOCK * hi)
                            for y, (lo, hi) in target_intervals.items()}
                        accepted_mass = sum(accepted.values())
                        if accepted_mass <= 0:
                            raise AssertionError("quantized accepted mass vanished")
                        normalized = {y: accepted[y] / accepted_mass
                                       for y in range(OUTPUTS)}
                        law_tv = tv_from_intervals(normalized, conditional_target)
                        law_tv_lower = tv_lower_from_intervals(
                            proposal, conditional_target)
                        omitted_acceptance_lowers.append(law_tv_lower)
                        target_mass = (sum(lo for lo, _ in target_intervals.values()),
                                       sum(hi for _, hi in target_intervals.values()))
                        full_error = max(
                            abs(float(sum(target_intervals[y]) / 2)
                                - independent[final, y])
                            for y in range(OUTPUTS))
                        success_lower = Fraction(qplan["success_probability_lower_bound"])
                        expected_upper = Fraction(qplan["expected_attempts_upper_bound"])
                        p1 &= (sum(proposal.values()) == 1
                               and proposal_ref_tv <= proposal_budget
                               and all(v <= qbound for v in denominator_tvs)
                               and accepted_mass >= success_lower
                               and Fraction(1, accepted_mass) <= expected_upper
                               and law_tv <= Fraction(qplan["total_tv_upper_bound"])
                               and sample_ok)
                        p2 &= (proposal_ref_tv <= proposal_budget
                               and all(v <= qbound for v in denominator_tvs)
                               and law_tv <= Fraction(qplan["total_tv_upper_bound"])
                               and target_mass[0] <= Fraction(1, BLOCK) <= target_mass[1]
                               and full_error < 2e-14)
                        row = dict(initial=initial, final=final,
                                   proposal_mass=str(sum(proposal.values())),
                                   accepted_mass=str(accepted_mass),
                                   proposal_vs_reference_tv=str(proposal_ref_tv),
                                   proposal_budget=str(proposal_budget),
                                   per_boundary_interval_tv=[str(v)
                                                             for v in denominator_tvs],
                                   fixed_boundary_tv=str(fixed_tv),
                                   normalized_tv_upper=str(law_tv),
                                   normalized_tv_lower=str(law_tv_lower),
                                   success_lower=str(success_lower),
                                   expected_attempts_upper=str(expected_upper),
                                   full_reference_max_error=full_error,
                                   costs=costs,
                                   sample_cost_check=sample_ok,
                                   plan=json_safe(qplan),
                                   proposal=[str(proposal[y]) for y in range(OUTPUTS)],
                                   reference_proposal=[str(reference_proposal[y])
                                                       for y in range(OUTPUTS)],
                                   normalized=[str(normalized[y])
                                              for y in range(OUTPUTS)],
                                   law_ok=(proposal_ref_tv <= proposal_budget
                                           and law_tv <= Fraction(qplan["total_tv_upper_bound"])))
                        report["rows"].append(dict(mask=mask,
                                                   target_tv=str(target_tv),
                                                   precision=precision, row=row))
        p3 = max(fixed_boundary_tvs) > Fraction(1, 1000)
        p4 = helper["helper_ok"]
        report["controls"] = dict(
            omitted_acceptance_tv_lower=str(max(omitted_acceptance_lowers)),
            omitted_acceptance_changes_target=max(omitted_acceptance_lowers) > Fraction(1, 1000),
            max_fixed_boundary_tv=str(max(fixed_boundary_tvs)),
            helper=helper,
            reference_target=str(REFERENCE_TARGET),
            report_rows=len(report["rows"]))
        exp.check("P1", p1, "proposal/reference/local budgets, accepted success and costs")
        exp.check("P2", p2, "target intervals and independent direct_joint diagnostics")
        exp.check("P3", p3, "fixed j differs from required uniform boundary mixture")
        exp.check("P4", p4, "canonical quantization helper edge cases")
        exp.fail_check("C1", max(omitted_acceptance_lowers) > Fraction(1, 1000),
                       f"proposal-vs-target TV lower={float(max(omitted_acceptance_lowers)):.6g}")
        report["status"] = "PASS" if p1 and p2 and p3 and p4 and max(omitted_acceptance_lowers) > Fraction(1, 1000) else "FAIL"
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]),
                          metadata=dict(period=PERIOD, block_size=BLOCK, width=WIDTH,
                                        histories=HISTORIES, initials=INITIALS,
                                        targets=[str(x) for x in TARGETS],
                                        precisions=PRECISIONS,
                                        allocation_bound_bytes=MAX_BYTES,
                                        method="VerifiedQuantizedReverseWork vs VerifiedReverseWork",
                                        no_forward_checkpoints=True,
                                        controls=json_safe(report["controls"]),
                                        elapsed_seconds=time.perf_counter() - started)):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = report_path("quantized_reverse_work_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
