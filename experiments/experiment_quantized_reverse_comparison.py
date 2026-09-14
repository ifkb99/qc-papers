"""Bounded comparison of locally quantized and unnormalized reverse work.

TODO29 asks whether projective dyadic integer coordinates can reduce the
precision burden of C69's reverse trajectory.  This experiment compares the
new quantized implementation with C69, all-forward finite work, checkpoint-8,
and the b=2 scalar baseline at the same requested TV.  Seeded outputs are not
expected to agree: these methods draw different proposal variables.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The quantized wide sampler completes without orbit/output tables and
      reports every attempt's aggregate rounding, operator, integer-bit,
      matvec, refinement, and bounded-working-state costs.
  P2  Quantized accepted laws on the matched width 0/1/3 endpoint fixtures
      normalize and agree with the independent finite-work law within the
      declared tiny diagnostic tolerance; the b=2 scalar comparison remains a
      separate strongest baseline.
  P3  Quantization exposes bounded coordinate/weight integer bitlengths and a
      fixed O(b^2) working scalar allowance.  Timing is reported descriptively;
      no faster-method prediction is made.
  C1  Omitting projective compression causes integer coordinates to grow on a
      bounded repeated scalar update using the same Gaussian-integer arithmetic.
  C2  Replacing physical child weights by norms after compression changes a
      tiny conditional law; weights must be taken before compression.

The scalar growth control is an arithmetic control, not a second propagator.
Integer/array bounds are structural payload guards, not native RSS claims.
"""
from __future__ import annotations

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
from lab.verified_quantized_reverse import (
    VerifiedQuantizedReverseWork, canonical_quantize, integer_matvec, integer_norm,
)
from lab.verified_reverse_work import VerifiedReverseWork
from lab.verified_finite_work import VerifiedFiniteWork, VerifiedScalarWork
from lab.sampling_error import plan_prefix_mass_accuracy


TARGET = Fraction(1, 1_000_000)
MAX_BYTES = 32 << 20
WIDE_PERIOD = 3 * ((1 << 60) - 1)
WIDE_WIDTH = 63
BLOCK = 3
SEEDS = (624, 625, 626, 627, 628)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "the quantized and comparison samplers complete the wide run with all requested costs")
exp.predict("P2", "matched tiny accepted laws normalize and agree within their declared diagnostic tolerance")
exp.predict("P3", "quantized coordinate/working-state bounds are explicit while timings remain descriptive")
exp.must_fail("C1", "omitting projective compression leaves bounded repeated integer coordinates unchanged")
exp.must_fail("C2", "using child norms after compression preserves the physical pre-compression weights")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"quantized_reverse_comparison_{stamp}.json"


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
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("retained experiment entries exceed 32 MiB allowance")


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


def summarize_samples(factory, *, seeds=SEEDS):
    rows = []
    for seed in seeds:
        started = time.perf_counter()
        worker = factory()
        setup_seconds = time.perf_counter() - started
        sample_started = time.perf_counter()
        result = worker.sample(Random(seed), target_tv=TARGET)
        sample_seconds = time.perf_counter() - sample_started
        rows.append({
            "seed": seed, "elapsed_seconds": time.perf_counter() - started,
            "setup_seconds": setup_seconds, "sample_seconds": sample_seconds,
            "output": result["output"],
            "initial_coarse_sector": result["initial_coarse_sector"],
            "final_coarse_sector": result["final_coarse_sector"],
            "attempts": result.get("attempts", 1),
            "total_tv_upper_bound": result.get("total_tv_upper_bound"),
            "orbit_or_output_tables": result.get("orbit_or_output_tables"),
            "max_working_precision": result.get("max_working_precision"),
            "grid_bits": result.get("grid_bits"),
            "operator_grid_bits": result.get("operator_grid_bits"),
            "initial_coordinate_bits_upper_bound": result.get("initial_coordinate_bits_upper_bound"),
            "operator_enclosure_evaluations": result.get("operator_enclosure_evaluations", 0),
            "vector_compressions": result.get("vector_compressions", 0),
            "reverse_child_evaluations": result.get("reverse_child_evaluations", 0),
            "reverse_vector_matvecs": result.get("reverse_vector_matvecs", 0),
            "replayed_vector_steps": result.get("replayed_vector_steps", 0),
            "branch_matrix_pair_constructions": result.get("branch_matrix_pair_constructions", 0),
            "refinement_retries": result.get("refinement_retries", result.get("refinements", 0)),
            "max_state_coordinate_bits": result.get("max_state_coordinate_bits"),
            "max_initial_coordinate_bits": result.get("max_initial_coordinate_bits"),
            "max_child_coordinate_bits": result.get("max_child_coordinate_bits"),
            "child_coordinate_bits_upper_bound": result.get("child_coordinate_bits_upper_bound"),
            "max_norm_weight_bits": result.get("max_norm_weight_bits"),
            "max_acceptance_integer_bits": result.get("max_acceptance_integer_bits"),
            "working_scalar_count_upper_bound": result.get("working_scalar_count_upper_bound"),
            "working_scalar_units": result.get("working_scalar_units"),
            "forward_steps": result.get("forward_steps", result.get("forward_builds", 0)),
            "recomputed_forward_steps": result.get("recomputed_forward_steps", 0),
            "peak_retained_forward_matrix_count": result.get("peak_retained_forward_matrix_count"),
            "stored_forward_matrix_count": result.get("stored_forward_matrix_count"),
            "stored_branch_matrix_count": result.get("stored_branch_matrix_count"),
            "stored_matrix_count_upper_bound": result.get("stored_matrix_count_upper_bound"),
            "retained_work_vector_count_upper_bound": result.get("retained_work_vector_count_upper_bound"),
            "working_matrix_scalar_count_upper_bound": result.get("working_matrix_scalar_count_upper_bound"),
            "diagnostic_zero_child_fallbacks": result.get("diagnostic_zero_child_fallbacks", 0),
        })
    elapsed = [r["elapsed_seconds"] for r in rows]
    return rows, {"elapsed_median_seconds": statistics.median(elapsed),
                  "elapsed_min_seconds": min(elapsed), "elapsed_max_seconds": max(elapsed),
                  "attempts": [r["attempts"] for r in rows],
                  "max_working_precision": max(r["max_working_precision"] for r in rows)}


def edge_fixture_rows():
    cases = [
        ("binary_width1_endpoint", VerifiedReflectionCircuit(
            10, 1, {0: ("x", Fraction(1, 4))},
            {0: (1, Fraction(1, 5)), 1: (0, Fraction(1, 5))}), 3),
        ("binary_width3_endpoint", VerifiedReflectionCircuit(
            10, 3, {0: ("x", Fraction(1, 4)), 1: ("z", Fraction(1, 7))},
            {0: (1, Fraction(1, 5)), 3: (0, Fraction(1, 5))}), 3),
        ("width0_endpoint", VerifiedReflectionCircuit(
            9, 0, {0: ("x", Fraction(1, 4))}, {0: (0, Fraction(1, 5))}, block_size=3), 1),
        ("width1_endpoint", VerifiedReflectionCircuit(
            9, 1, {0: ("x", Fraction(1, 4)), 1: ("z", Fraction(1, 7))},
            {0: (1, Fraction(1, 5)), 1: (0, Fraction(1, 5))}, block_size=3), 3),
        ("width3_endpoint", VerifiedReflectionCircuit(
            9, 3, {0: ("x", Fraction(1, 4)), 1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
            {0: (1, Fraction(1, 5)), 3: (0, Fraction(1, 5))}, block_size=3), 3),
    ]
    guard_entries(2048)
    rows = []
    for name, circuit, mask in cases:
        quantized = VerifiedQuantizedReverseWork(circuit, mask)
        finite = VerifiedFiniteWork(circuit, mask)
        M, Q, b = circuit.sectors, 1 << circuit.width, circuit.b
        accepted_by_initial = []
        for initial in range(M):
            accepted = {}
            for boundary in range(b):
                for output in range(Q):
                    result = quantized.attempt(initial, boundary, output=output, target_tv=TARGET)
                    key = (result["final_coarse_sector"], output)
                    accepted[key] = accepted.get(key, Fraction(0)) + (
                        result["proposal_path_probability"] * result["acceptance_probability"] / b)
            success = sum(accepted.values())
            if success <= 0:
                raise AssertionError(f"zero accepted mass in fixed initial sector {initial}")
            accepted_by_initial.append((success, {
                key: value / success for key, value in accepted.items()}))
        # The production sampler holds initial fixed through rejection.  Only
        # after each sector's accepted submeasure is normalized do we mix the
        # sectors uniformly.
        normalized = {}
        for success, sector_law in accepted_by_initial:
            for key, value in sector_law.items():
                normalized[key] = normalized.get(key, Fraction(0)) + value / M
        ideal = {}
        for initial in range(M):
            for output in range(Q):
                result = finite.path(initial, output=output, target_tv=TARGET)
                key = (result["final_coarse_sector"], output)
                ideal[key] = ideal.get(key, Fraction(0)) + result["conditional_path_probability"] / M
        tv = sum(abs(normalized.get(key, 0) - ideal.get(key, 0))
                 for key in set(normalized) | set(ideal)) / 2
        planned_sum = (quantized.plan(TARGET)["total_tv_upper_bound"]
                       + plan_prefix_mass_accuracy(circuit.width, TARGET)["total_tv_upper_bound"])
        rows.append({"name": name, "sector_success_masses": [mass for mass, _ in accepted_by_initial],
                     "sector_normalized_masses": [sum(law.values()) for _, law in accepted_by_initial],
                     "success_mass": sum(mass for mass, _ in accepted_by_initial) / M,
                     "normalized_mass": sum(normalized.values()),
                     "ideal_mass": sum(ideal.values()), "tv_vs_finite": tv,
                     "planned_tv_bound_sum": planned_sum,
                     "entries": len(normalized)})
    return rows


def no_compression_control():
    # Same Gaussian-integer matvec arithmetic as the quantized implementation;
    # this toy only isolates coordinate growth, not a new circuit propagator.
    matrix = (((1, 0), (1, 0)), ((1, 0), (1, 0)))
    vector = ((1, 0), (0, 0))
    initial_bits = max(abs(x).bit_length() for z in vector for x in z)
    trace = []
    for _ in range(10):
        vector = integer_matvec(matrix, vector)
        trace.append(max(abs(x).bit_length() for z in vector for x in z))
    return trace[-1] > initial_bits, {"initial_bits": initial_bits,
                                      "final_bits": trace[-1], "trace": trace}


def weights_after_compression_control():
    circuit = circuit3(width=3, backgrounds=True, reflections=True)
    worker = VerifiedQuantizedReverseWork(circuit, 3)
    plan = worker.plan(TARGET)
    for initial in range(circuit.sectors):
        for boundary in range(circuit.b):
            cursor = worker.cursor(initial, boundary, target_tv=TARGET)
            for depth in range(circuit.width):
                physical = cursor.weights()
                children = cursor._pending[1]
                compressed = [canonical_quantize(child, plan["grid_bits"])[0]
                              for child in children]
                wrong = tuple(integer_norm(child) for child in compressed)
                if sum(physical) and sum(wrong):
                    left = tuple(Fraction(x, sum(physical)) for x in physical)
                    right = tuple(Fraction(x, sum(wrong)) for x in wrong)
                    if left != right:
                        return True, {"initial": initial, "boundary": boundary,
                                      "depth": depth, "physical": physical,
                                      "post_compression": wrong,
                                      "physical_normalized": left,
                                      "post_normalized": right}
                cursor.advance(0)
    return False, {}


def main():
    report = {"status": "PASS", "wide": {}, "scalar_null": {},
              "edges": [], "controls": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        guard((9, 16), label="tiny diagnostic law")
        wide = wide_circuit()
        report["wide"]["quantized"], report["wide"]["quantized_summary"] = summarize_samples(
            lambda: VerifiedQuantizedReverseWork(wide, 0))
        report["wide"]["raw_reverse"], report["wide"]["raw_reverse_summary"] = summarize_samples(
            lambda: VerifiedReverseWork(wide, 0))
        report["wide"]["full"], report["wide"]["full_summary"] = summarize_samples(
            lambda: VerifiedFiniteWork(wide, 0))
        report["wide"]["checkpoint8"], report["wide"]["checkpoint8_summary"] = summarize_samples(
            lambda: VerifiedFiniteWork(wide, 0, checkpoint_spacing=8))
        qrows = report["wide"]["quantized"]
        p1 &= all(not row["orbit_or_output_tables"] and row["attempts"] >= 1 for row in qrows)
        for kind in ("quantized", "raw_reverse", "full", "checkpoint8"):
            p1 &= all(Fraction(row["total_tv_upper_bound"]) <= TARGET
                      for row in report["wide"][kind])
        p1 &= all(row["vector_compressions"] > 0 and row["reverse_vector_matvecs"] > 0
                  for row in qrows)
        p1 &= all(row["vector_compressions"] == WIDE_WIDTH * row["attempts"]
                  and row["reverse_vector_matvecs"] == 2 * WIDE_WIDTH * row["attempts"]
                  and row["branch_matrix_pair_constructions"] == WIDE_WIDTH * row["attempts"]
                  and row["operator_enclosure_evaluations"] == (WIDE_WIDTH + 1) * row["attempts"]
                  for row in qrows)
        p1 &= all(row["working_scalar_count_upper_bound"] == 240 for row in qrows)
        p1 &= all(row["max_state_coordinate_bits"] >= 1 and row["max_child_coordinate_bits"] >= 1
                  and row["max_norm_weight_bits"] >= 1 and row["max_acceptance_integer_bits"] >= 1
                  for row in qrows)
        p3 &= all(row["max_state_coordinate_bits"] <= row["grid_bits"] + 1 for row in qrows)
        p3 &= all(row["max_initial_coordinate_bits"] <= row["initial_coordinate_bits_upper_bound"]
                  and row["working_scalar_units"] for row in qrows)
        p3 &= all(row["max_child_coordinate_bits"] <= row["child_coordinate_bits_upper_bound"]
                  and row["replayed_vector_steps"] == row["forward_steps"] == 0
                  and row["stored_forward_matrix_count"] == row["stored_branch_matrix_count"] == 0
                  and row["diagnostic_zero_child_fallbacks"] == 0 for row in qrows)
        p3 &= all(row["working_scalar_count_upper_bound"] == 24 * BLOCK * BLOCK + 8 * BLOCK
                  for row in qrows)

        scalar = b2_scalar_circuit()
        for name, factory in (
                ("quantized", lambda: VerifiedQuantizedReverseWork(scalar, 3)),
                ("raw_reverse", lambda: VerifiedReverseWork(scalar, 3)),
                ("finite", lambda: VerifiedFiniteWork(scalar, 3)),
                ("scalar", lambda: VerifiedScalarWork(scalar, 3))):
            report["scalar_null"][name], report["scalar_null"][name + "_summary"] = summarize_samples(factory)
        report["edges"] = edge_fixture_rows()
        p2 &= all(0 < Fraction(row["success_mass"]) <= 1
                  and all(mass > 0 for mass in row["sector_success_masses"])
                  and all(mass == 1 for mass in row["sector_normalized_masses"])
                  and row["normalized_mass"] == 1
                  and row["ideal_mass"] == 1
                  and row["tv_vs_finite"] <= row["planned_tv_bound_sum"]
                  for row in report["edges"])
        growth_failed, growth_detail = no_compression_control()
        wrong_failed, wrong_detail = weights_after_compression_control()
        report["controls"] = {"no_compression": growth_detail,
                              "weights_after_compression": wrong_detail}
        exp.check("P1", p1, "wide quantized/reverse/forward comparison completed")
        exp.check("P2", p2, "endpoint accepted laws normalized against finite work")
        exp.check("P3", p3, "quantized bits and working scalar bounds recorded")
        exp.fail_check("C1", growth_failed, "uncompressed integer coordinates grow on the bounded scalar update")
        exp.fail_check("C2", wrong_failed, "post-compression weights differ from physical pre-compression weights")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
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
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
