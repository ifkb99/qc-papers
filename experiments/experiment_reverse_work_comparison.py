"""Bounded comparison of reverse-vector and forward finite-work samplers.

This is an implementation/cost audit for C68, not a claim of a universal
speedup.  The reverse sampler keeps an unnormalised work vector and uses a
fresh terminal boundary label on each rejection attempt; the forward baselines
are ``VerifiedFiniteWork`` with all forward matrices or spacing eight.  The
same target is used throughout, but seeded outputs are deliberately not
expected to agree because the decompositions draw different random variables.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Wide b=3,r=3*(2^60-1),t=63 samples complete with no forward/branch
      tables, and report the reverse/forward costs rather than silently
      conditioning on a capped run.
  P2  The reverse and k=8 workers preserve the declared acceptance/proposal
      invariants on the tiny b=2 scalar baseline and b=3 edge fixtures;
      endpoint W0/route and a width-32 no-background zero case complete.
  P3  Reverse attempts retain O(b) work-vector storage counters while the
      all-forward baseline retains a width-dependent matrix bound; timing is
      reported as a scoped diagnostic only.
  C1  Omitting terminal acceptance must not equal the independent full-r
      target law on a routed tiny fixture.
  C2  Resetting a reverse cursor at an already measured prefix must change a
      subsequent conditional weight vector for at least one tiny case.

No pathwise equality is predicted between the samplers.  Matrix scalar bounds
are conservative algorithmic allowances, not native RSS measurements.
"""
from __future__ import annotations

import json
import math
import statistics
import time
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from random import Random

import numpy as np

from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_reverse_work import VerifiedReverseWork
from lab.verified_finite_work import VerifiedFiniteWork, VerifiedScalarWork


TARGET = Fraction(1, 1_000_000)
MAX_BYTES = 32 << 20
WIDE_PERIOD = 3 * ((1 << 60) - 1)
WIDE_WIDTH = 63
BLOCK = 3
SEEDS = (624, 625, 626, 627, 628)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "the wide reverse/vector and forward baselines finish with explicit work and storage counters")
exp.predict("P2", "tiny scalar, endpoint, and zero-prefix cases preserve valid proposal/acceptance ranges")
exp.predict("P3", "reverse retained-vector counters stay O(b), while all-forward storage is width dependent")
exp.must_fail("C1", "dropping terminal acceptance leaves the routed reverse proposal equal to the independent target")
exp.must_fail("C2", "resetting a reverse vector at a measured prefix leaves its next conditional weights unchanged")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"reverse_work_comparison_{stamp}.json"


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
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def guard_entries(entries, bytes_per_entry=512):
    if type(entries) is not int or entries < 0:
        raise ValueError("retained entry count must be a nonnegative integer")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("retained diagnostic entries exceed 32 MiB allowance")


def circuit3(width=4, *, backgrounds=True, reflections=True):
    bg = ({s: ("x", Fraction(1, 7)) for s in range(1, width + 1)}
          if backgrounds else {})
    refs = ({2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}
            if reflections and width >= 3 else {})
    return VerifiedReflectionCircuit(9, width, bg, refs, block_size=3)


def wide_circuit():
    bg = {s: ("x", Fraction(1, 7)) for s in range(1, WIDE_WIDTH + 1)}
    return VerifiedReflectionCircuit(WIDE_PERIOD, WIDE_WIDTH, bg, {}, block_size=3)


def b2_scalar_circuit():
    bg = {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))}
    refs = {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}
    return VerifiedReflectionCircuit(6, 3, bg, refs, block_size=2)


def summarize_samples(worker, *, seeds=SEEDS):
    rows = []
    for seed in seeds:
        started = time.perf_counter()
        instance = worker() if callable(worker) else worker
        setup_elapsed = time.perf_counter() - started
        sample_started = time.perf_counter()
        result = instance.sample(Random(seed), target_tv=TARGET)
        sample_elapsed = time.perf_counter() - sample_started
        elapsed = time.perf_counter() - started
        rows.append({
            "seed": seed, "elapsed_seconds": elapsed,
            "setup_seconds": setup_elapsed, "sample_seconds": sample_elapsed,
            "output": result["output"],
            "initial_coarse_sector": result["initial_coarse_sector"],
            "final_coarse_sector": result["final_coarse_sector"],
            "attempts": result.get("attempts", 1),
            "max_working_precision": result.get("max_working_precision"),
            "total_tv_upper_bound": result["total_tv_upper_bound"],
            "orbit_or_output_tables": result["orbit_or_output_tables"],
            "initial_vector_builds": result.get("initial_vector_builds",0),
            "acceptance_inner_products": result.get("acceptance_inner_products",0),
            "reverse_child_evaluations": result.get("reverse_child_evaluations", 0),
            "reverse_vector_matvecs": result.get("reverse_vector_matvecs", 0),
            "replayed_vector_steps": result.get("replayed_vector_steps", 0),
            "branch_matrix_pair_constructions": result.get("branch_matrix_pair_constructions", 0),
            "forward_steps": result.get("forward_steps", result.get("forward_builds", 0)),
            "recomputed_forward_steps": result.get("recomputed_forward_steps", 0),
            "stored_forward_matrix_count": result.get("stored_forward_matrix_count"),
            "peak_retained_forward_matrix_count": result.get("peak_retained_forward_matrix_count"),
            "stored_branch_matrix_count": result.get("stored_branch_matrix_count"),
            "stored_matrix_count_upper_bound": result.get("stored_matrix_count_upper_bound"),
            "retained_work_vector_count_upper_bound": result.get("retained_work_vector_count_upper_bound"),
            "working_matrix_scalar_count_upper_bound": result.get("working_matrix_scalar_count_upper_bound"),
            "zero_approximate_block_fallbacks": result.get("zero_approximate_block_fallbacks", 0),
        })
    elapsed = [row["elapsed_seconds"] for row in rows]
    return rows, {
        "elapsed_median_seconds": statistics.median(elapsed),
        "elapsed_min_seconds": min(elapsed), "elapsed_max_seconds": max(elapsed),
        "attempts": [row["attempts"] for row in rows],
        "max_working_precision": max(row["max_working_precision"] for row in rows),
    }


def valid_attempt(result):
    acceptance = result["acceptance_probability"]
    return (Fraction(0) <= acceptance <= Fraction(1)
            and result["proposal_path_probability"] >= 0
            and result["proposal_path_probability"] <= 1
            and result["stored_forward_matrix_count"] == 0
            and result["stored_branch_matrix_count"] == 0
            and result["retained_work_vector_count_upper_bound"] <= 4)


def edge_rows():
    rows = []
    cases = [
        ("b3_width0_W0_route", VerifiedReflectionCircuit(
            9, 0, {0: ("x", Fraction(1, 4))}, {0: (0, Fraction(1, 5))}, block_size=3), 1),
        ("b3_width1_endpoint", VerifiedReflectionCircuit(
            9, 1, {0: ("x", Fraction(1, 4)), 1: ("z", Fraction(1, 7))},
            {0: (1, Fraction(1, 5)), 1: (0, Fraction(1, 5))}, block_size=3), 3),
        ("b3_width3_endpoint", VerifiedReflectionCircuit(
            9, 3, {0: ("x", Fraction(1, 4)), 1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
            {0: (1, Fraction(1, 5)), 3: (0, Fraction(1, 5))}, block_size=3), 3),
        ("b3_width32_no_background", VerifiedReflectionCircuit(9, 32, {}, {}, block_size=3), 0),
    ]
    guard_entries(4096, bytes_per_entry=512)
    for name, circuit, mask in cases:
        guard((circuit.sectors, 1 << min(circuit.width, 5)), label=name)
        reverse = VerifiedReverseWork(circuit, mask)
        finite = VerifiedFiniteWork(circuit, mask)
        forced = []
        # All boundaries/sectors, complete tiny laws; only two outputs at t=32.
        initials = range(circuit.sectors)
        boundaries = range(circuit.b)
        outputs = range(1 << circuit.width) if circuit.width <= 3 else (0, 1)
        for initial in initials:
            for boundary in boundaries:
                for output in outputs:
                    result = reverse.attempt(initial, boundary, output=output, target_tv=TARGET)
                    forced.append({
                        "initial": initial, "boundary": boundary, "output": output,
                        "proposal_path_probability": result["proposal_path_probability"],
                        "acceptance_probability": result["acceptance_probability"],
                        "zero_approximate_block_fallbacks": result["zero_approximate_block_fallbacks"],
                        "max_working_precision": result["max_working_precision"],
                        "final": result["final_coarse_sector"],
                        "valid": valid_attempt(result),
                    })
        edge_laws = []
        if circuit.width <= 3:
            for initial in initials:
                accepted = {y:sum((row["proposal_path_probability"]*row["acceptance_probability"]/circuit.b
                    for row in forced if row["initial"] == initial and row["output"] == y),Fraction(0)) for y in outputs}
                mass = sum(accepted.values())
                baseline = {y:finite.path(initial,output=y,target_tv=TARGET) for y in outputs}
                joint_ok = all(row["final"] == baseline[row["output"]]["final_coarse_sector"]
                               for row in forced if row["initial"] == initial)
                tv = sum(abs(accepted[y]/mass-baseline[y]["conditional_path_probability"]) for y in outputs)/2
                bound = reverse.plan(TARGET)["total_tv_upper_bound"]+max(row["total_tv_upper_bound"] for row in baseline.values())
                edge_laws.append(dict(initial=initial,accepted_mass=mass,finite_baseline_tv=tv,
                                      comparison_bound=bound,ok=joint_ok and tv<=bound))
        finite_probe = finite.path(0, output=0, target_tv=TARGET)
        rows.append({"name": name, "mask": mask, "forced": forced,
                     "all_forced_valid": all(x["valid"] for x in forced),
                     "reverse_zero_fallbacks": sum(x["zero_approximate_block_fallbacks"] for x in forced),
                     "complete_edge_laws":edge_laws,
                     "finite_probe": {k: finite_probe.get(k) for k in (
                         "conditional_path_probability", "zero_approximate_block_fallbacks",
                         "max_working_precision", "forward_steps")}})
    return rows


def reverse_proposal_law(worker, outputs):
    """Mix forced reverse proposal paths over uniform initial sector and boundary."""
    M, b = worker.circuit.sectors, worker.circuit.b
    law = {}
    for initial in range(M):
        for boundary in range(b):
            for output in outputs:
                result = worker.attempt(initial, boundary, output=output, target_tv=TARGET)
                key = (result["final_coarse_sector"], output)
                law[key] = law.get(key, Fraction(0)) + result["proposal_path_probability"] / (M * b)
    return law


def routed_float_target():
    from lab.coherent_routes import CoherentReflectionCircuit
    from experiments.experiment_coherent_route_sampling import direct_joint
    def rx(theta):
        return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                         [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)
    bg = {}
    for s, matrix in ((0, rx(np.pi / 4)), (1, rx(np.pi / 7)), (3, np.diag(
            [np.exp(-1j * np.pi / 10), np.exp(1j * np.pi / 10)]).astype(complex))):
        embedded = np.eye(3, dtype=complex); embedded[:2, :2] = matrix
        bg[s] = embedded
    comparator = CoherentReflectionCircuit(9, 3, 3, bg, {2: (0, np.pi), 3: (1, np.pi)})
    return direct_joint(comparator)


def omission_control():
    # Keep the independent full-r comparator at the same tiny width.
    circuit = VerifiedReflectionCircuit(
        9, 3,
        {0: ("x", Fraction(1, 4)), 1: ("x", Fraction(1, 7)),
         3: ("z", Fraction(1, 5))},
        {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}, block_size=3)
    worker = VerifiedReverseWork(circuit, 3)
    proposal = reverse_proposal_law(worker, range(1 << circuit.width))
    target = routed_float_target()
    tv = 0.0
    for gamma in range(circuit.sectors):
        for output in range(1 << circuit.width):
            tv += abs(float(proposal.get((gamma, output), 0)) - float(target[gamma, output]))
    return tv / 2, sum(proposal.values()), len(proposal)


def reset_control():
    circuit = circuit3(width=3, backgrounds=True, reflections=True)
    worker = VerifiedReverseWork(circuit, 3)
    bits = worker.plan(TARGET)["proposal_plan"]["weight_accuracy_bits"]
    for initial in range(circuit.sectors):
        for boundary in range(circuit.b):
            cursor = worker.cursor(initial, boundary, accuracy_bits=bits)
            cursor.weights(); cursor.advance(0)
            correct = cursor.weights()
            # Reset to the measured boundary basis e_j, rather than the
            # physical initial vector, to model the forbidden checkpoint.
            from flint import acb_mat
            cursor.vector = acb_mat([[int(i == boundary)] for i in range(circuit.b)])
            cursor._pending = None
            reset = cursor.weights()
            def normalized(weights):
                total = sum(weights)
                return tuple(Fraction(w, total) for w in weights) if total else None
            if normalized(correct) != normalized(reset):
                return True, {"initial": initial, "boundary": boundary,
                              "correct": correct, "reset": reset,
                              "correct_normalized": normalized(correct),
                              "reset_normalized": normalized(reset)}
    return False, {}


def main():
    report = {"status": "PASS", "wide": {}, "scalar_null": {}, "edges": [], "controls": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        # These are only small reference allocations; the wide workers allocate no
        # orbit/output arrays.  Keep the explicit product guard before tiny law work.
        guard((9, 16), label="tiny diagnostic law")
        # Frozen five-seed samples, all retained edge/control rows, current
        # exact-law dictionaries and serialization copies, not native RSS.
        guard_entries(16000)
        wide = wide_circuit()
        report["wide"]["reverse"], report["wide"]["reverse_summary"] = summarize_samples(
            lambda: VerifiedReverseWork(wide, 0))
        report["wide"]["full"], report["wide"]["full_summary"] = summarize_samples(
            lambda: VerifiedFiniteWork(wide, 0))
        report["wide"]["checkpoint8"], report["wide"]["checkpoint8_summary"] = summarize_samples(
            lambda: VerifiedFiniteWork(wide, 0, checkpoint_spacing=8))
        for kind in ("reverse", "full", "checkpoint8"):
            rows = report["wide"][kind]
            p1 &= all(row["stored_forward_matrix_count"] == 0 for row in rows if kind == "reverse")
            p1 &= all(row["stored_branch_matrix_count"] == 0 for row in rows if kind == "reverse")
            p1 &= all(row["max_working_precision"] >= 64 for row in rows)
            p1 &= all(row["total_tv_upper_bound"] <= TARGET and not row["orbit_or_output_tables"]
                      and 0 <= row["output"] < 1 << WIDE_WIDTH for row in rows)
        p3 &= all(row["retained_work_vector_count_upper_bound"] == 4
                  for row in report["wide"]["reverse"])
        # Compare like-scoped working allowances; do not compare reverse
        # vectors alone against forward matrices PLUS temporary allowances.
        reverse_scalars = max(row["working_matrix_scalar_count_upper_bound"] for row in report["wide"]["reverse"])
        full_scalars = report["wide"]["full"][0]["stored_matrix_count_upper_bound"] * BLOCK**2
        checkpoint_scalars = report["wide"]["checkpoint8"][0]["stored_matrix_count_upper_bound"] * BLOCK**2
        p3 &= reverse_scalars < checkpoint_scalars < full_scalars
        p3 &= all(row["stored_matrix_count_upper_bound"] > 0
                  for row in report["wide"]["full"] + report["wide"]["checkpoint8"])

        scalar = b2_scalar_circuit()
        scalar_rows = {}
        for name, worker in (("reverse", lambda: VerifiedReverseWork(scalar, 3)),
                             ("finite", lambda: VerifiedFiniteWork(scalar, 3)),
                             ("scalar", lambda: VerifiedScalarWork(scalar, 3))):
            scalar_rows[name], scalar_rows[name + "_summary"] = summarize_samples(worker)
        report["scalar_null"] = scalar_rows
        report["edges"] = edge_rows()
        p2 &= all(case["all_forced_valid"] for case in report["edges"])
        p2 &= all(row["ok"] for case in report["edges"] for row in case["complete_edge_laws"])
        p2 &= all(Fraction(case["finite_probe"]["conditional_path_probability"]) >= 0
                  for case in report["edges"])
        p2 &= any(case["reverse_zero_fallbacks"] > 0
                  for case in report["edges"]
                  if case["name"] == "b3_width32_no_background")

        omission_tv, proposal_mass, support = omission_control()
        report["controls"]["omission"] = {"proposal_tv_vs_float": omission_tv,
                                              "proposal_mass": proposal_mass, "support": support}
        reset_changed, reset_detail = reset_control()
        report["controls"]["reset"] = {"changed": reset_changed, "detail": reset_detail}
        exp.check("P1", p1, "wide workers completed with explicit counters")
        exp.check("P2", p2, "tiny edge proposal/acceptance checks")
        exp.check("P3", p3, "retained/recomputed storage counters reported")
        exp.fail_check("C1", omission_tv > 1e-4 and proposal_mass == 1,
                       f"omitting acceptance TV={omission_tv:.6g}, proposal mass={proposal_mass}")
        exp.fail_check("C2", reset_changed, "reset conditional weights differ in a tiny routed case")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        # Resolve every prediction/control as failed, preserving the raw report.
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
        exp.fail_check("C2", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "PASS" if all(ok for _, _, ok, _ in exp._results) else "FAIL"
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)],
                    metadata={"target_tv": str(TARGET), "seeds": SEEDS,
                              "retained_entry_allowance":16000,
                              "setup_scope":"worker snapshot plus sample; supplied input circuit built once outside timer",
                              "storage_scope":"structural scalar-entry allowances, not measured native memory or RSS",
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
