"""Matched bounded cost comparison for sparse coherent reverse routing.

The two-label support audit established a linear reached-sector envelope.  This
experiment measures the resulting float-only ``SparseCoherentReverse`` against
the existing coherent-history rejection sampler and the stronger gate-by-gate
prefix sampler.  It varies only the number of supplied coherent reflections
and reports setup, attempts, history work, sparse support, and sampler costs.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every k in {0,2,4,6,8} completes on the fixed wide input (or records a
      genuine capped resource failure), with no orbit/history/output table;
      the sparse reverse peak sector support obeys its declared two-label bound.
  P2  The tiny r=9,b=3,t=4 full-r sequential_path marginal agrees with the
      existing coherent-history contraction for the same output marginal.
  P3  Reports charge actual rejection attempts, sparse work, and streamed
      history work; no timing order is predicted.
  C1  Ignoring rejection multiplicity must fail when a history sampler makes
      more than one proposal: actual history work then exceeds one-history work.
  C2  Applying the two-label support bound to the known generic eight-insertion
      word must fail.

The wide fixture is r=3*(2^40-1), b=3, t=15, with W0 and alternating Rx/Rz
block mixers, and q=0/1 alternating reflections at fixed insertions.  The
float samplers are diagnostics only; no finite-TV certificate or speedup claim
is made.  Matrix payload and report-entry guards are explicit and not RSS.
"""
from __future__ import annotations

import argparse
import itertools
import math
import statistics
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.coherent_reverse import SparseCoherentReverse
from lab.semiclassical import sequential_path


MAX_BYTES = 32 << 20
PERIOD = 3 * ((1 << 40) - 1)
BLOCK = 3
SECTORS = PERIOD // BLOCK
WIDTH = 15
POSITIONS = (0, 2, 4, 6, 8, 10, 12, 15)
KS = (0, 2, 4, 6, 8)
SEEDS = (7130, 7131, 7132)
MAX_PROPOSALS = 2048

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all requested k rows finish or preserve explicit capped failures without tables")
exp.predict("P2", "tiny full-r sequential_path and coherent-history output marginals agree")
exp.predict("P3", "actual attempts, streamed histories, support and work counters are charged")
exp.must_fail("C1", "ignoring rejection multiplicity gives the actual history-work cost")
exp.must_fail("C2", "the two-label support bound applies to the generic eight-insertion word")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"coherent_merge_comparison_{stamp}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
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


def guard_entries(entries, bytes_per_entry=1024):
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("comparison report/retained payload exceeds 32 MiB")


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)]).astype(complex)


def embed(matrix):
    result = np.eye(BLOCK, dtype=complex)
    result[:2, :2] = matrix
    return result


def wide_fixture(k, second_route=1):
    defects = {0: embed(rx(np.pi / 4))}
    for s in range(1, WIDTH + 1):
        defects[s] = embed(rx(np.pi / 7) if s % 2 else rz(np.pi / 5))
    reflections = {s: (0 if i % 2 == 0 else second_route, np.pi / 4)
                   for i, s in enumerate(POSITIONS[:k])}
    guard((BLOCK, BLOCK), label="wide block mixer")
    return CoherentReflectionCircuit(PERIOD, BLOCK, WIDTH, defects, reflections)


def detail(result, *, method, k, seed, elapsed, setup_seconds, counters):
    # Preserve the raw field (BGL reports zero) separately from the logical
    # completed-call count used for timing/accounting.
    raw_rejection_proposals = int(result.get("rejection_proposals", 0))
    attempts = raw_rejection_proposals if raw_rejection_proposals > 0 else 1
    history_count = int(result.get("history_count_upper_bound", 0))
    return {
        "method": method, "k": k, "seed": seed,
        "elapsed_seconds": elapsed, "setup_seconds": setup_seconds,
        "sample_seconds": elapsed - setup_seconds,
        "attempts": attempts, "completed_calls": 1,
        "rejection_proposals": raw_rejection_proposals,
        "history_count_upper_bound": history_count,
        "history_component_evaluations": counters["history_components"],
        "history_work_units": counters["history_components"],
        "history_generator_calls": counters["history_calls"],
        "component_sample_calls": counters["component_sample_calls"],
        "prefix_vector_queries": counters["prefix_queries"],
        "history_enumerations": result.get("history_enumerations"),
        "support_bound": result.get("support_bound"),
        "peak_sector_count": result.get("peak_sector_count"),
        "reflection_vector_contributions": result.get("reflection_vector_contributions"),
        "work_matvecs": result.get("work_matvecs"),
        "qft_matrix_pair_constructions": result.get("qft_matrix_pair_constructions"),
        "prefix_vector_evaluations": result.get("prefix_vector_evaluations"),
        "resampled_blocks": result.get("resampled_blocks"),
        "output": result.get("output"),
        "final_coarse_sector": result.get("final_coarse_sector"),
        "orbit_table_entries": result.get("orbit_table_entries"),
        "sector_table_entries": result.get("sector_table_entries"),
        "output_table_entries": result.get("output_table_entries"),
        "full_sector_table_entries": result.get("full_sector_table_entries"),
        "history_storage": result.get("history_storage"),
        "working_complex_coordinate_slots_upper_bound": result.get("working_complex_coordinate_slots_upper_bound"),
        "supplied_gate_complex_entries": result.get("supplied_gate_complex_entries"),
        "owned_matrix_payload_bytes": result.get("owned_matrix_payload_bytes"),
        "method_reported": result.get("method"),
        "support_bound_formula": result.get("support_bound"),
    }


def run_method(k, seed, method, second_route=1):
    started = time.perf_counter()
    circuit = wide_fixture(k, second_route)
    if circuit.sectors != SECTORS or any(
            circuit.reflections[s][0] != (0 if i % 2 == 0 else second_route)
            for i, s in enumerate(POSITIONS[:k])):
        raise AssertionError("reported route labels differ from realized circuit labels")
    worker = SparseCoherentReverse(circuit) if method == "sparse_reverse" else circuit
    target = worker.circuit if method == "sparse_reverse" else circuit
    counters = {"history_components": 0, "history_calls": 0,
                "component_sample_calls": 0, "prefix_queries": 0}
    original_histories = target._histories
    def counted_histories(*args, **kwargs):
        counters["history_calls"] += 1
        for item in original_histories(*args, **kwargs):
            counters["history_components"] += 1
            yield item
    target._histories = counted_histories
    original_component_sample = target._component_sample
    def counted_component_sample(*args, **kwargs):
        counters["component_sample_calls"] += 1
        return original_component_sample(*args, **kwargs)
    target._component_sample = counted_component_sample
    original_prefix_vector = target.prefix_vector
    def counted_prefix_vector(*args, **kwargs):
        counters["prefix_queries"] += 1
        return original_prefix_vector(*args, **kwargs)
    target.prefix_vector = counted_prefix_vector
    setup_seconds = time.perf_counter() - started
    rng = np.random.default_rng(seed)
    if method == "sparse_reverse":
        result = worker.sample(rng, max_proposals=MAX_PROPOSALS)
    elif method == "history_rejection":
        result = circuit.sample_rejection(rng, max_proposals=MAX_PROPOSALS)
    elif method == "gate_by_gate":
        result = circuit.sample(rng)
    else:
        raise ValueError(method)
    elapsed = time.perf_counter() - started
    return detail(result, method=method, k=k, seed=seed, elapsed=elapsed,
                  setup_seconds=setup_seconds, counters=counters)


def tiny_static_comparator():
    circuit = CoherentReflectionCircuit(
        9, 3, 4,
        {0: embed(rx(np.pi / 4)), 1: embed(rx(np.pi / 7)),
         2: embed(rz(np.pi / 5)), 3: embed(rx(np.pi / 7)), 4: embed(rz(np.pi / 5))},
        {0: (0, np.pi / 4), 2: (1, np.pi / 4), 4: (0, np.pi / 4)})
    guard((9, 9), label="tiny static full-r reference")
    from experiments.experiment_coherent_route_sampling import full_gates
    started = time.perf_counter()
    pairs, initial = full_gates(circuit)
    sequential = np.array([
        sequential_path(pairs, initial, output=y)["conditional_path_probability"]
        for y in range(1 << circuit.width)])
    coherent = np.array([
        sum(circuit.joint_probability(gamma, y) for gamma in range(circuit.sectors))
        for y in range(1 << circuit.width)])
    elapsed = time.perf_counter() - started
    return {"diagnostic_seconds": elapsed, "timing_scope": "full_gates plus both output-marginal references; not a sampler timing",
            "sequential_mass": float(sequential.sum()),
            "coherent_mass": float(coherent.sum()),
            "max_abs_error": float(np.max(np.abs(sequential - coherent))),
            "l1_tv_diagnostic": float(np.sum(np.abs(sequential - coherent)) / 2),
            "outputs": len(sequential), "full_r_pair_count": len(pairs),
            "setup_included": True}


def generic_control():
    M, gamma = 1009, 17
    word = (0, 1, 0, 3, 0, 9, 0, 27)
    reached = {gamma}
    for q in reversed(word):
        reached |= {(-x - q) % M for x in reached}
    return len(reached), min(M, 1 + 2 * len(word))


def main(second_route=1):
    report = {"status": "PASS", "preflight": [], "rows": [], "summaries": {},
              "tiny_static": {}, "controls": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        if not (0 < int(second_route) < SECTORS) or math.gcd(int(second_route), SECTORS) != 1:
            raise ValueError("second route must satisfy 0<q<M and gcd(q,M)=1")
        second_route = int(second_route)
        report["benchmark"] = {"q0": 0, "q1": second_route,
                                "period": PERIOD, "sectors": SECTORS,
                                "gcd_q1_sectors": math.gcd(second_route, SECTORS),
                                "label": "fixed q0=0 and supplied q1 across all k; tiny_static remains q0/q1=(0,1)"
                                if second_route != 1 else "default adjacent q1=1 fixture"}
        guard_entries(3 * len(KS) * 3 + 32)
        # Required budget preflight: run one k=8 sample per method before the
        # full sweep, retaining failures rather than silently reducing work.
        for method in ("sparse_reverse", "history_rejection", "gate_by_gate"):
            try:
                row = run_method(8, SEEDS[0], method, second_route)
                report["preflight"].append(row)
            except Exception as exc:
                report["preflight"].append({"method": method, "k": 8,
                                             "error": repr(exc)})
                p1 = False
        for k in KS:
            for seed in SEEDS:
                for method in ("sparse_reverse", "history_rejection", "gate_by_gate"):
                    try:
                        row = run_method(k, seed, method, second_route)
                    except Exception as exc:
                        row = {"method": method, "k": k, "seed": seed,
                               "error": repr(exc)}
                        p1 = False
                    report["rows"].append(row)
        for k in KS:
            for method in ("sparse_reverse", "history_rejection", "gate_by_gate"):
                group = [r for r in report["rows"]
                         if r.get("k") == k and r.get("method") == method
                         and "error" not in r]
                if group:
                    elapsed = [r["elapsed_seconds"] for r in group]
                    report["summaries"][f"k{k}_{method}"] = {
                        "k": k, "method": method, "count": len(group),
                        "elapsed_median_seconds": statistics.median(elapsed),
                        "elapsed_min_seconds": min(elapsed),
                        "elapsed_max_seconds": max(elapsed),
                        "attempts": [r["attempts"] for r in group],
                        "history_component_evaluations": [r["history_component_evaluations"] for r in group],
                        "prefix_vector_queries": [r["prefix_vector_queries"] for r in group],
                    }
        sparse_rows = [r for r in report["rows"] if r["method"] == "sparse_reverse" and "error" not in r]
        p1 &= len(sparse_rows) == len(KS) * len(SEEDS)
        p1 &= all(r["peak_sector_count"] <= r["support_bound"]
                  and r["orbit_table_entries"] == 0 for r in sparse_rows)
        p1 &= all(r["history_enumerations"] == 0 for r in sparse_rows)
        p1 &= all(0 <= r["output"] < (1 << WIDTH)
                  and 0 <= r["final_coarse_sector"] < PERIOD // BLOCK
                  for r in report["rows"] if "error" not in r)
        p1 &= all(16 * (r["working_complex_coordinate_slots_upper_bound"]
                       + r["supplied_gate_complex_entries"]) < MAX_BYTES
                  for r in sparse_rows)
        p1 &= all(all(r.get(name) in (0, None) for name in (
            "orbit_table_entries", "sector_table_entries", "output_table_entries",
            "full_sector_table_entries")) for r in report["rows"] if "error" not in r)
        p3 &= all("error" not in r for r in report["rows"])
        p3 &= all(r["attempts"] >= 1 and r["elapsed_seconds"] > 0 for r in report["rows"])
        p3 &= all(r["history_work_units"] == r["attempts"] * r["history_count_upper_bound"]
                  and r["component_sample_calls"] == r["attempts"]
                  for r in report["rows"] if r["method"] == "history_rejection")
        p3 &= all(r["history_work_units"] == 0 and r["prefix_vector_queries"] == 0
                  for r in report["rows"] if r["method"] == "sparse_reverse")
        p3 &= all(r["prefix_vector_queries"] > 0 and r["history_work_units"] > 0
                  and r["rejection_proposals"] == 0
                  and r["prefix_vector_queries"] == r["prefix_vector_evaluations"]
                  for r in report["rows"] if r["method"] == "gate_by_gate")

        report["tiny_static"] = tiny_static_comparator()
        p2 &= abs(report["tiny_static"]["sequential_mass"] - 1.0) < 1e-12
        p2 &= abs(report["tiny_static"]["coherent_mass"] - 1.0) < 1e-12
        p2 &= report["tiny_static"]["max_abs_error"] < 2e-12

        actual, bound = generic_control()
        rejection_failure = any(r["attempts"] > 1 and
                                r["history_work_units"] > r["history_count_upper_bound"]
                                for r in report["rows"]
                                if r["method"] == "history_rejection" and "error" not in r)
        report["controls"] = {"rejection_cost": {"failed": rejection_failure},
                              "generic_two_label_bound": {"support": actual, "bound": bound}}
        exp.check("P1", p1, "wide k/seed/method rows and sparse support bounds")
        exp.check("P2", p2, "tiny full-r sequential_path marginal")
        exp.check("P3", p3, "attempt/history/sparse work accounting")
        exp.fail_check("C1", rejection_failure,
                       "history work includes attempts rather than one uncharged history pass")
        exp.fail_check("C2", actual > bound,
                       f"generic support {actual} exceeds two-label bound {bound}")
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
                    metadata={"k_values": KS, "seeds": SEEDS,
                              "q0": 0, "q1": second_route, "period": PERIOD,
                              "sectors": SECTORS,
                              "gcd_q1_sectors": math.gcd(second_route, SECTORS),
                              "benchmark_label": report.get("benchmark", {}).get("label"),
                              "max_proposals": MAX_PROPOSALS,
                              "native_memory_measured": False,
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="bounded coherent-route merge comparison")
    parser.add_argument("--second-route", type=int, default=1,
                        help="fixed q1 route label in (0,M), coprime to the sector modulus M=r/b")
    main(parser.parse_args().second_route)
