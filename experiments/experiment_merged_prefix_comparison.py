"""Matched wide cost comparison for the merged-prefix baseline.

This compares four supplied-input samplers on exactly the same wide indexed
orbit: MergedCoherentPrefixes, SparseCoherentReverse, the existing C59
gate-by-gate sampler, and C59 coherent-history rejection.  The two q1 values
and each k row are frozen independently; no timing result is pooled across
route fixtures.  Setup, adapter installation, returned-sample work, actual
history components, rejection proposals, and merged prefix queries are kept
as separate fields.  Native RSS and a universal speed claim are out of scope.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  all requested q1/k/method rows complete within the explicit proposal
      cap without orbit/output/history tables; every method reports its own
      actual work counters.
  P2  per-(q1,k,method) summaries remain separately attributable, with no
      history work charged to merged/reverse methods and no raw rejection work
      charged to C59 gate-by-gate sampling.
  P3  the tiny static full-r comparator remains normalized; it is a correctness
      diagnostic, not a sampler timing.

  C1  ignoring prefix queries, history components, or rejected proposals must
      fail to represent the charged baseline work.

This is a supplied-input complex128 diagnostic comparison.  The adjacent and
large-q1 rows are separate frozen fixtures; the C71 two-label shortcut is not
silently generalized into a generic-route claim.
"""
from __future__ import annotations

import json
import math
import statistics
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_reverse import SparseCoherentReverse
from lab.merged_prefix import MergedCoherentPrefixes

from experiments.experiment_coherent_merge_comparison import (
    BLOCK, PERIOD, POSITIONS, SEEDS, WIDTH, KS, MAX_PROPOSALS,
    wide_fixture, tiny_static_comparator)


MAX_BYTES = 16 << 20
Q1_VALUES = (1, 679535556937)
PROPOSAL_CAP = MAX_PROPOSALS

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all frozen wide q1/k/method rows complete with bounded table-free work")
exp.predict("P2", "method-specific setup, adapter, history, rejection, and prefix costs remain separated")
exp.predict("P3", "tiny static full-r correctness diagnostic remains normalized")
exp.must_fail("C1", "ignoring query/history/rejection work misstates the comparison")


def stamp(prefix="merged_prefix_comparison"):
    now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{now}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"{prefix}_{now}_{serial}.json"
    return path


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


def guard_entries(entries, bytes_per_entry=1024):
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("comparison report/retained payload exceeds 16 MiB")


def install_counters(circuit, *, merged=None):
    """Count only work belonging to the selected implementation."""
    counts = dict(history_calls=0, history_components=0,
                  component_sample_calls=0, prefix_queries=0,
                  first_prefix_block_products=None)
    target = circuit
    original_histories = target._histories

    def counted_histories(*args, **kwargs):
        counts["history_calls"] += 1
        for item in original_histories(*args, **kwargs):
            counts["history_components"] += 1
            yield item

    target._histories = counted_histories
    original_component = target._component_sample

    def counted_component(*args, **kwargs):
        counts["component_sample_calls"] += 1
        return original_component(*args, **kwargs)

    target._component_sample = counted_component
    prefix_owner = merged if merged is not None else target
    original_prefix = prefix_owner.prefix_vector

    def counted_prefix(*args, **kwargs):
        counts["prefix_queries"] += 1
        value = original_prefix(*args, **kwargs)
        if merged is not None and counts["first_prefix_block_products"] is None:
            local = merged.last_prefix_stats or {}
            counts["first_prefix_block_products"] = (
                int(local.get("background_block_products", 0))
                + int(local.get("control_block_products", 0)))
        return value

    prefix_owner.prefix_vector = counted_prefix
    return counts


def detail(result, *, method, q1, k, seed, elapsed, setup_seconds,
           adapter_seconds, counts, worker_stats):
    attempts = int(result.get("rejection_proposals", 0))
    return {
        "method": method, "q1": q1, "k": k, "seed": seed,
        "elapsed_seconds": elapsed,
        "setup_seconds": setup_seconds,
        "input_setup_seconds": setup_seconds - adapter_seconds,
        "adapter_seconds": adapter_seconds,
        "sample_seconds": elapsed - setup_seconds,
        "charged_wall_scope": "setup (including adapter) + returned sample; not native RSS",
        "attempts": attempts if attempts else 1,
        "rejection_proposals": attempts,
        "history_calls": counts["history_calls"],
        "history_components": counts["history_components"],
        "history_work_units": counts["history_components"],
        "component_sample_calls": counts["component_sample_calls"],
        "prefix_queries": counts["prefix_queries"],
        "first_prefix_block_products": counts["first_prefix_block_products"],
        "total_block_products": (
            int(result.get("background_block_products", 0))
            + int(result.get("control_block_products", 0))),
        "history_enumerations": result.get("history_enumerations"),
        "support_bound": result.get("support_bound"),
        "peak_sector_count": result.get("peak_sector_count"),
        "reflection_vector_contributions": result.get("reflection_vector_contributions"),
        "work_matvecs": result.get("work_matvecs"),
        "qft_matrix_pair_constructions": result.get("qft_matrix_pair_constructions"),
        "prefix_vector_evaluations": result.get("prefix_vector_evaluations"),
        "resampled_blocks": result.get("resampled_blocks"),
        "reverse_steps": result.get("reverse_steps"),
        "output": result.get("output"),
        "final_coarse_sector": result.get("final_coarse_sector"),
        "orbit_table_entries": result.get("orbit_table_entries"),
        "sector_table_entries": result.get("sector_table_entries"),
        "output_table_entries": result.get("output_table_entries"),
        "full_sector_table_entries": result.get("full_sector_table_entries"),
        "history_storage": result.get("history_storage"),
        "working_complex_coordinate_slots_upper_bound": result.get(
            "working_complex_coordinate_slots_upper_bound"),
        "supplied_gate_complex_entries": result.get("supplied_gate_complex_entries"),
        "worker_stats": worker_stats,
        "method_reported": result.get("method"),
    }


def run_method(k, q1, seed, method):
    started = time.perf_counter()
    circuit = wide_fixture(k, q1)
    if circuit.sectors != PERIOD // BLOCK:
        raise AssertionError("wide fixture sector count changed")
    expected_labels = tuple(0 if i % 2 == 0 else q1
                            for i in range(k))
    realized_labels = tuple(circuit.reflections[s][0]
                            for s in POSITIONS[:k])
    if realized_labels != expected_labels:
        raise AssertionError(f"route labels differ: {realized_labels} != {expected_labels}")

    adapter_started = time.perf_counter()
    merged = None
    if method == "merged_prefix":
        merged = MergedCoherentPrefixes(circuit)
        worker = merged
        counter_target = worker.circuit
    elif method == "sparse_reverse":
        worker = SparseCoherentReverse(circuit)
        counter_target = worker.circuit
    else:
        worker = circuit
        counter_target = circuit
    counts = install_counters(counter_target, merged=merged)
    adapter_seconds = time.perf_counter() - adapter_started
    setup_seconds = time.perf_counter() - started
    rng = np.random.default_rng(seed)
    if method == "merged_prefix":
        result = worker.sample(rng)
        worker_stats = worker.stats()
    elif method == "sparse_reverse":
        result = worker.sample(rng, max_proposals=PROPOSAL_CAP)
        worker_stats = worker._stats()
    elif method == "gate_by_gate":
        result = worker.sample(rng)
        worker_stats = worker.stats()
    elif method == "history_rejection":
        result = worker.sample_rejection(rng, max_proposals=PROPOSAL_CAP)
        worker_stats = worker.stats()
    else:
        raise ValueError(method)
    elapsed = time.perf_counter() - started
    return detail(result, method=method, q1=q1, k=k, seed=seed,
                  elapsed=elapsed, setup_seconds=setup_seconds,
                  adapter_seconds=adapter_seconds, counts=counts,
                  worker_stats=worker_stats)


def summarize(rows):
    summaries = {}
    for q1 in Q1_VALUES:
        for k in KS:
            for method in ("merged_prefix", "sparse_reverse", "gate_by_gate",
                           "history_rejection"):
                group = [row for row in rows
                         if row.get("q1") == q1 and row.get("k") == k
                         and row.get("method") == method and "error" not in row]
                if not group:
                    continue
                key = f"q1_{q1}_k{k}_{method}"
                summaries[key] = {
                    "q1": q1, "k": k, "method": method,
                    "count": len(group),
                    "elapsed_seconds_median": statistics.median(
                        row["elapsed_seconds"] for row in group),
                    "elapsed_seconds_range": [
                        min(row["elapsed_seconds"] for row in group),
                        max(row["elapsed_seconds"] for row in group)],
                    "setup_seconds_median": statistics.median(
                        row["setup_seconds"] for row in group),
                    "sample_seconds_median": statistics.median(
                        row["sample_seconds"] for row in group),
                    "attempts": [row["attempts"] for row in group],
                    "history_components": [row["history_components"] for row in group],
                    "prefix_queries": [row["prefix_queries"] for row in group],
                    "rejection_proposals": [row["rejection_proposals"] for row in group],
                }
    return summaries


def main():
    started = time.perf_counter()
    report = {"status": "PASS", "preflight": [], "rows": [],
              "summaries": {}, "controls": {}}
    methods = ("merged_prefix", "sparse_reverse", "gate_by_gate",
               "history_rejection")
    p1 = p2 = p3 = True
    try:
        if math.gcd(PERIOD // BLOCK, Q1_VALUES[1]) != 1:
            raise AssertionError("large q1 fixture is not coprime to M")
        guard_entries(2 * len(Q1_VALUES) * len(KS) * len(SEEDS) * len(methods) + 128)
        # Required preflight: each frozen route fixture gets one k=8 run per
        # method before the complete timing sweep.
        for q1 in Q1_VALUES:
            for method in methods:
                try:
                    report["preflight"].append(run_method(8, q1, SEEDS[0], method))
                except Exception as exc:
                    report["preflight"].append({"q1": q1, "k": 8,
                                                 "method": method,
                                                 "error": repr(exc),
                                                 "traceback": traceback.format_exc()})
                    p1 = False
        for q1 in Q1_VALUES:
            for k in KS:
                for seed in SEEDS:
                    for method in methods:
                        try:
                            row = run_method(k, q1, seed, method)
                        except Exception as exc:
                            row = {"q1": q1, "k": k, "seed": seed,
                                   "method": method, "error": repr(exc),
                                   "traceback": traceback.format_exc()}
                            p1 = False
                        report["rows"].append(row)
        report["summaries"] = summarize(report["rows"])
        expected_rows = len(Q1_VALUES) * len(KS) * len(SEEDS) * len(methods)
        p1 &= len(report["rows"]) == expected_rows
        valid = [row for row in report["rows"] if "error" not in row]
        p1 &= len(valid) == expected_rows
        p1 &= all(row["history_enumerations"] in (0, None)
                  and row["orbit_table_entries"] in (0, None)
                  and row["sector_table_entries"] in (0, None)
                  and row["output_table_entries"] in (0, None)
                  and row["full_sector_table_entries"] in (0, None)
                  for row in valid)
        p1 &= all(0 <= row["output"] < (1 << WIDTH)
                  and 0 <= row["final_coarse_sector"] < PERIOD // BLOCK
                  for row in valid)
        p1 &= all(
            row["peak_sector_count"] <= row["support_bound"]
            and 16 * (row["working_complex_coordinate_slots_upper_bound"]
                      + row["supplied_gate_complex_entries"]) < MAX_BYTES
            for row in valid if row["method"] in ("merged_prefix", "sparse_reverse"))
        # The merged and reverse methods must not be charged with C59 history
        # work; the C59 gate-by-gate method must not be charged with rejection.
        merged_rows = [row for row in valid if row["method"] == "merged_prefix"]
        sparse_rows = [row for row in valid if row["method"] == "sparse_reverse"]
        gate_rows = [row for row in valid if row["method"] == "gate_by_gate"]
        rejection_rows = [row for row in valid if row["method"] == "history_rejection"]
        p2 &= all(row["history_components"] == 0
                  and row["rejection_proposals"] == 0
                  and row["prefix_queries"] > 0
                  and row["prefix_queries"] == row["prefix_vector_evaluations"]
                  and row["total_block_products"] > row["first_prefix_block_products"]
                  for row in merged_rows)
        p2 &= all(row["history_components"] == 0
                  and row["prefix_queries"] == 0
                  and row["rejection_proposals"] >= 1 for row in sparse_rows)
        p2 &= all(row["history_components"] > 0
                  and row["rejection_proposals"] == 0
                  and row["prefix_queries"] > 0
                  and row["prefix_queries"] == row["prefix_vector_evaluations"]
                  for row in gate_rows)
        p2 &= all(row["history_components"] > 0
                  and row["component_sample_calls"] == row["attempts"]
                  and row["rejection_proposals"] == row["attempts"]
                  and row["history_components"] == row["attempts"] * (1 << row["k"])
                  for row in rejection_rows)
        tiny = tiny_static_comparator()
        report["tiny_static"] = tiny
        p3 &= abs(tiny["sequential_mass"] - 1.0) < 1e-12
        p3 &= abs(tiny["coherent_mass"] - 1.0) < 1e-12
        p3 &= tiny["max_abs_error"] < 2e-12
        # Must-fail accounting control: omitted retries understate a history
        # sampler, and charging one merged prefix understates its total block
        # products.  Both are required to expose the accounting mistake.
        retry_understatement = any(
            row["method"] == "history_rejection"
            and row["history_components"] > (1 << row["k"])
            for row in valid)
        merged_understatement = any(
            row["method"] == "merged_prefix"
            and row["total_block_products"] > row["first_prefix_block_products"]
            for row in valid)
        c1 = retry_understatement and merged_understatement
        report["controls"] = {
            "ignored_work_differs": c1,
            "retry_understatement": retry_understatement,
            "merged_single_prefix_understatement": merged_understatement,
            "actual_nonzero_history_rows": sum(row["history_components"] > 0
                                                for row in valid),
            "actual_nonzero_prefix_rows": sum(row["prefix_queries"] > 0
                                               for row in valid),
            "actual_retried_rows": sum(row["rejection_proposals"] > 1
                                        for row in valid),
        }
        exp.check("P1", p1, "wide rows complete, bounded, and table-free")
        exp.check("P2", p2, "method-specific work attribution")
        exp.check("P3", p3, "tiny static full-r diagnostic")
        exp.fail_check("C1", c1, "omitting retries and charging only one prefix both understate actual work")
        report["status"] = "PASS" if all((p1, p2, p3, c1)) else "FAIL"
        path = stamp()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]),
                          metadata=json_safe({
                              "period": PERIOD, "block_size": BLOCK,
                              "sectors": PERIOD // BLOCK, "width": WIDTH,
                              "positions": POSITIONS, "k_values": KS,
                              "seeds": SEEDS, "q1_values": Q1_VALUES,
                              "proposal_cap": PROPOSAL_CAP,
                              "allocation_bound_bytes": MAX_BYTES,
                              "native_memory_measured": False,
                              "timing_scope": "setup (including adapter)+returned sample",
                              "c71_caveat": "adjacent q1=1 admits a cheap smooth-state approximation; exact task timings do not establish hardness or superiority over approximate grouping",
                              "preflight": report["preflight"],
                              "tiny_static": report.get("tiny_static"),
                              "summaries": report["summaries"],
                              "controls": report["controls"],
                              "elapsed_seconds": time.perf_counter() - started})):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = stamp("merged_prefix_comparison_failure")
        path.write_text(json.dumps(json_safe(report), indent=2) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
