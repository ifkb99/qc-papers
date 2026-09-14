"""Bounded same-output comparison for the C78 work-first sampler.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  LateWorkProgressions' enumerated k=1 joint law agrees with an independent
      existing work-60 sequential_path reference; the separately computed k=0
      omit-G reference agrees with the cached FB omit-G row.
  P2  The feedback-aware FB comparator and explicit k=0 omit-G baseline are
      reported against the same complete output law at TV 1e-3 and 1e-2.
      These are diagnostic float laws, not a finite-bit sampling certificate.
  P3  A few bounded helper samples exercise work selection and coherent
      rejection, with returned counters and component bounds retained.

C1  Replacing coherent conditional rows by the fully dephased uniform output
    law must disagree with the target.
C2  Redrawing the work label after every rejected output must produce the
    component-count-weighted wrong law when row component counts vary.

The physical fixture is the TODO37/C78 N=61, a=2, r=60, b=3, t=6 schedule,
with split=5 (L=32,H=2).  The common terminal W is omitted because it is a
work-only unitary traced by the output measurement.  The independent reference
uses the existing finite-work sequential_path instrument with supplied dense
60-by-60 branches; this is a tiny exact reference, not a new propagator.
The helper's actual forced/sample counters and setup/reference work are charged
separately.  No timing, scaling, or certified arbitrary-precision conclusion
is made.
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_uniform_prefix_output import (
    _repeated_block as repeated_block_60,
    _shift as shift_matrix_60,
)
from lab import Experiment
from lab.fourier_sampling import unit_phase
from lab.semiclassical import sequential_path
from lab.work_first import LateWorkProgressions


N = 61
BASE = 2
PERIOD = 60
BLOCK = 3
WIDTH = 6
SPLIT = 5
Q = 1 << WIDTH
L = 1 << SPLIT
H = 1 << (WIDTH - SPLIT)
BYTE_CAP = 16 * 1024 * 1024
MAX_FORCED_CALLS = PERIOD * Q
MAX_SAMPLE_CALLS = 8
MAX_ATTEMPTS = 10_000
MAX_REFERENCE_CALLS = 2 * Q
MAX_REFERENCE_WORK = MAX_REFERENCE_CALLS * WIDTH * PERIOD**3
MAX_REFERENCE_SAMPLE_CALLS = 4
MAX_SAMPLE_ATTEMPT_TOTAL = 256
MAX_REFERENCE_DENSE_TERMS = 32 * PERIOD**3
MAX_REFERENCE_PHASE_QUERIES = 2 * PERIOD
TOL = 3e-10

FB_REPORT = Path("out/conditional_feedback_boundary_20260911T104152447926Z.json")


def report_path() -> Path:
    path = Path("out") / (
        "work_first_sampler_comparison_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        + ".json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(payload: int, label: str) -> int:
    payload = int(payload)
    if payload < 0 or payload > BYTE_CAP:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def guard_shape(shape, dtype=np.complex128, label="array") -> int:
    return guard_bytes(math.prod(int(x) for x in shape)
                       * np.dtype(dtype).itemsize, label)


def tv(left: np.ndarray, right: np.ndarray) -> float:
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != right.shape:
        raise ValueError("TV shape mismatch")
    return float(np.sum(np.abs(left - right)) / 2)


def valid_law(values: np.ndarray, size: int = Q) -> bool:
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (size,) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) < TOL)


def phase(index: int) -> complex:
    """The supplied pointwise G_1 character, queried by orbit index."""
    return unit_phase(pow(BASE, int(index), N), N)


def load_fb_laws() -> dict[str, np.ndarray]:
    if not FB_REPORT.exists():
        raise FileNotFoundError(f"required FB provenance report missing: {FB_REPORT}")
    payload = json.loads(FB_REPORT.read_text())
    row = payload["rows"][0]
    full = {str(k): np.asarray(v, dtype=float)
            for k, v in row["full_laws"].items()}
    candidate_conditionals = {
        int(z): np.asarray(v, dtype=float)
        for z, v in row["candidate_conditionals"].items()
    }
    fb_h = len(candidate_conditionals)
    if fb_h < 1 or Q % fb_h:
        raise ValueError("invalid cached FB candidate prefix dimensions")
    candidate = np.empty(Q, dtype=float)
    for z in range(fb_h):
        candidate[z::fb_h] = candidate_conditionals[z] / fb_h
    return {"full_k0": full["0"], "full_k1": full["1"],
            "candidate_k1": candidate}


def reference_pairs(k: int, counters: dict) -> tuple[list[tuple[np.ndarray, np.ndarray]], np.ndarray]:
    """Build supplied work-60 sequential branches for k=0 or k=1."""
    if k not in (0, 1):
        raise ValueError("reference character k must be 0 or 1")
    guard_shape((PERIOD, PERIOD), label="reference identity")
    guard_shape((PERIOD, PERIOD), label="reference G diagonal")
    # Count the matrices actually allocated below, not both pair references:
    # the branch-zero matrix is shared, whereas after@shift is a new matrix.
    counters["reference_setup_matrix_allocated_entries"] += 4 * PERIOD**2
    identity = np.eye(PERIOD, dtype=complex)
    W0 = repeated_block_60(work_block(math.pi / 4))
    W1 = repeated_block_60(work_block(-math.pi / 10))
    initial = W0[:, 0]
    if k:
        phase_values = []
        for j in range(PERIOD):
            if counters["reference_phase_queries"] + 1 > MAX_REFERENCE_PHASE_QUERIES:
                raise MemoryError("reference phase-query cap before query")
            counters["reference_phase_queries"] += 1
            phase_values.append(phase(j))
    else:
        phase_values = [1.+0j] * PERIOD
    g = np.diag(phase_values)
    if counters["reference_dense_matmul_terms"] + PERIOD**3 > MAX_REFERENCE_DENSE_TERMS:
        raise MemoryError("reference dense setup cap before matrix product")
    counters["reference_dense_matmul_terms"] += PERIOD**3
    counters["reference_setup_matrix_allocated_entries"] += PERIOD**2
    post = g @ W1
    pairs = []
    for i in range(WIDTH):
        counters["reference_shift_entry_writes"] += PERIOD
        counters["reference_setup_matrix_allocated_entries"] += PERIOD**2
        shift = shift_matrix_60(1 << i)
        after = post if i == WIDTH - 2 else identity
        if counters["reference_dense_matmul_terms"] + PERIOD**3 > MAX_REFERENCE_DENSE_TERMS:
            raise MemoryError("reference dense setup cap before matrix product")
        counters["reference_dense_matmul_terms"] += PERIOD**3
        counters["reference_pair_count"] += 2
        counters["reference_pair_matrix_allocated_entries"] += PERIOD**2
        pairs.append((after, after @ shift))
    return pairs, initial


def sequential_law(k: int, counters: dict) -> np.ndarray:
    pairs, initial = reference_pairs(k, counters)
    result = np.zeros(Q, dtype=float)
    for output in range(Q):
        if counters["reference_calls"] + 1 > MAX_REFERENCE_CALLS:
            raise MemoryError("sequential reference call cap before forced call")
        counters["reference_calls"] += 1
        counters["reference_work_units"] += WIDTH * PERIOD**3
        value = sequential_path(pairs, initial, output=output,
                                max_payload_bytes=4 << 20)
        probability = float(value["conditional_path_probability"])
        if not math.isfinite(probability):
            raise ArithmeticError("nonfinite sequential reference probability")
        result[output] = probability
    return result


def aggregate_counter(total: dict, current: dict) -> None:
    for key, value in current.items():
        if isinstance(value, (int, np.integer)):
            total[key] = total.get(key, 0) + int(value)


def main() -> None:
    exp = Experiment("work_first_sampler_comparison", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "forced helper law agrees with independent work-60 sequential references")
    exp.predict("P2", "FB candidate and omit-G baselines are compared at both TV thresholds")
    exp.predict("P3", "bounded helper samples exercise work selection and coherent rejection")
    exp.must_fail("C1", "fully dephased uniform output differs from the target")
    exp.must_fail("C2", "redrawing work labels after rejection biases varying component rows")

    counters = {
        "forced_calls": 0, "sample_calls": 0,
        "reference_calls": 0, "reference_work_units": 0,
        "reference_sample_calls": 0, "reference_sample_work_units": 0,
        "reference_pair_count": 0,
        "reference_pair_matrix_allocated_entries": 0,
        "reference_setup_matrix_allocated_entries": 0,
        "reference_shift_entry_writes": 0,
        "reference_dense_matmul_terms": 0,
        "reference_phase_queries": 0,
        "helper_forced": {}, "helper_samples": {},
    }
    report = {"status": "FAIL", "rows": [], "counters": counters}
    p1 = p2 = p3 = c1 = c2 = False
    try:
        # Before allocating any dense reference/helper arrays: retained 60x60
        # branch matrices, pair products, joint laws, and conservative QR/
        # sequential temporaries. Python JSON containers are not claimed as RSS.
        preflight_components = {
            "reference_setup_matrices": 8 * PERIOD * PERIOD * 16,
            "reference_pair_matrices": 2 * WIDTH * PERIOD * PERIOD * 16,
            "reference_sequential_temporaries": 46 * PERIOD * PERIOD * 16,
            "helper_joint_law": PERIOD * Q * 8,
            "output_law_arrays": 16 * Q * 8,
            "helper_numeric_payload_reserve": 16 * (
                16 * BLOCK * BLOCK + 12 * (H * min(PERIOD, 2 * BLOCK - 1)) + 128),
        }
        planned = sum(preflight_components.values())
        guard_bytes(planned, "aggregate sampler comparison preflight")
        reference_total_work = MAX_REFERENCE_WORK + MAX_REFERENCE_SAMPLE_CALLS * WIDTH * PERIOD**3
        if MAX_FORCED_CALLS > 4_000 or reference_total_work > 200_000_000:
            raise MemoryError("forced/reference work cap before allocation")

        fb = load_fb_laws()
        for key, law in fb.items():
            if not valid_law(law):
                raise ValueError(f"invalid cached FB law {key}")

        W0 = work_block(math.pi / 4)
        W1 = work_block(-math.pi / 10)
        worker = LateWorkProgressions(
            PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, phase,
            max_local_terms=1_000_000, max_components=4096,
            max_payload_bytes=BYTE_CAP)
        helper_stats = worker.stats()
        helper_numeric_reserve_ok = (
            helper_stats["numeric_payload_bound_bytes"]
            <= preflight_components["helper_numeric_payload_reserve"])
        if not helper_numeric_reserve_ok:
            raise AssertionError("helper numeric reserve omitted from preflight")
        guard_shape((PERIOD, Q), np.float64, "helper joint law")
        helper_joint = np.zeros((PERIOD, Q), dtype=float)
        work_probabilities = np.zeros(PERIOD, dtype=float)
        component_counts = np.zeros(PERIOD, dtype=np.int64)
        helper_counter_totals = {}
        forced_result_fields = 0
        for work in range(PERIOD):
            for output in range(Q):
                if counters["forced_calls"] + 1 > MAX_FORCED_CALLS:
                    raise MemoryError("forced helper call cap before call")
                counters["forced_calls"] += 1
                result = worker.forced_joint(work, output)
                if "joint_probability" not in result or "work_probability" not in result:
                    raise AssertionError("forced helper result schema changed")
                joint = float(result["joint_probability"])
                work_probability = float(result["work_probability"])
                if not math.isfinite(joint) or joint < -TOL:
                    raise ArithmeticError("invalid forced joint probability")
                if not math.isfinite(work_probability) or work_probability < -TOL:
                    raise ArithmeticError("invalid forced work probability")
                if output == 0:
                    work_probabilities[work] = work_probability
                    component_counts[work] = int(result["component_count"])
                    if work_probability <= 0 or component_counts[work] <= 0:
                        raise ArithmeticError("nonpositive forced work row")
                elif abs(work_probabilities[work] - work_probability) > TOL:
                    raise AssertionError("work probability changed across forced outputs")
                elif int(result["component_count"]) != int(component_counts[work]):
                    raise AssertionError("component count changed across forced outputs")
                helper_joint[work, output] = joint
                aggregate_counter(helper_counter_totals, result["counters"])
                forced_result_fields += 1

        helper_law = helper_joint.sum(axis=0)
        helper_work_law = helper_joint.sum(axis=1)
        exact_k1 = sequential_law(1, counters)
        exact_k0 = sequential_law(0, counters)
        reference_sample_rows = []
        for k in (0, 1):
            pairs, initial = reference_pairs(k, counters)
            for seed in (39_01 + k, 39_11 + k):
                if counters["reference_sample_calls"] + 1 > MAX_REFERENCE_SAMPLE_CALLS:
                    raise MemoryError("reference sample call cap before call")
                counters["reference_sample_calls"] += 1
                counters["reference_sample_work_units"] += WIDTH * PERIOD**3
                sampled = sequential_path(
                    pairs, initial, rng=np.random.default_rng(seed),
                    max_payload_bytes=4 << 20)
                reference_law = exact_k1 if k else exact_k0
                if (not 0 <= sampled["output"] < Q
                        or reference_law[sampled["output"]] <= 0
                        or abs(reference_law[sampled["output"]]
                               - sampled["conditional_path_probability"]) > TOL):
                    raise AssertionError("returned reference draw disagrees with forced law")
                reference_sample_rows.append({
                    "k": k, "seed": seed, "output": int(sampled["output"]),
                    "path_probability": float(sampled["conditional_path_probability"]),
                    "steps": int(sampled["steps"]),
                })
        candidate = fb["candidate_k1"]
        omit = exact_k0
        wrong_redraw_unnormalized = np.sum(
            helper_joint / component_counts[:, None], axis=0)
        wrong_redraw = wrong_redraw_unnormalized / float(wrong_redraw_unnormalized.sum())
        uniform = np.full(Q, 1. / Q)

        sample_rows = []
        sample_counter_totals = {}
        for seed in range(38_01, 38_01 + MAX_SAMPLE_CALLS):
            if counters["sample_calls"] >= MAX_SAMPLE_CALLS:
                raise MemoryError("sample call cap before call")
            rng = np.random.default_rng(seed)
            counters["sample_calls"] += 1
            remaining_attempt_budget = MAX_SAMPLE_ATTEMPT_TOTAL - sum(
                row["attempts"] for row in sample_rows)
            if remaining_attempt_budget <= 0:
                raise MemoryError("sample attempt budget exhausted before call")
            attempt_cap = min(MAX_ATTEMPTS, remaining_attempt_budget)
            result = worker.sample(rng, max_attempts=attempt_cap)
            if not all(key in result for key in ("work", "output", "attempts",
                                                  "component_count", "counters")):
                raise AssertionError("sample helper result schema changed")
            if not (0 <= int(result["work"]) < PERIOD
                    and 0 <= int(result["output"]) < Q
                    and 1 <= int(result["attempts"]) <= MAX_ATTEMPTS
                    and 1 <= int(result["component_count"]) <= helper_stats["component_bound"]):
                raise AssertionError("invalid bounded sample result")
            aggregate_counter(sample_counter_totals, result["counters"])
            sample_rows.append({
                "seed": seed, "work": int(result["work"]),
                "output": int(result["output"]),
                "attempts": int(result["attempts"]),
                "component_count": int(result["component_count"]),
                "attempt_cap": attempt_cap,
                "sample_counter_ok": bool(
                    result["counters"].get("work_draws") == 1
                    and result["counters"].get("component_draws") == int(result["attempts"])
                    and result["counters"].get("progression_proposals") == int(result["attempts"])
                    and result["counters"].get("acceptance_draws") == int(result["attempts"])
                    and result["counters"]["column_local_terms"] == BLOCK**2
                    and result["counters"]["row_local_terms"] == H*BLOCK**2
                    and result["counters"]["fourier_component_terms"]
                        == result["attempts"]*result["component_count"]
                    and result["counters"]["marginal_queries"] == 8*result["attempts"]
                    and helper_joint[result["work"],result["output"]] > 0),
                "counters": {str(k): int(v) for k, v in result["counters"].items()},
            })

        helper_forced_expected = {
            "row_local_terms": MAX_FORCED_CALLS * H * BLOCK**2,
            "phase_queries": MAX_FORCED_CALLS * H,
            "fourier_component_terms": int(Q * np.sum(component_counts)),
        }
        forced_counter_ok = all(
            helper_counter_totals.get(key, 0) == value
            for key, value in helper_forced_expected.items())
        reference_cost_expected = {
            "reference_calls": MAX_REFERENCE_CALLS,
            "reference_work_units": MAX_REFERENCE_CALLS * WIDTH * PERIOD**3,
            "reference_sample_calls": MAX_REFERENCE_SAMPLE_CALLS,
            "reference_sample_work_units": MAX_REFERENCE_SAMPLE_CALLS * WIDTH * PERIOD**3,
            "reference_pair_count": 4 * 2 * WIDTH,
            "reference_pair_matrix_allocated_entries": 4 * WIDTH * PERIOD**2,
            "reference_setup_matrix_allocated_entries": 4 * (5 + WIDTH) * PERIOD**2,
            "reference_shift_entry_writes": 4 * WIDTH * PERIOD,
            "reference_dense_matmul_terms": 4 * (WIDTH + 1) * PERIOD**3,
            "reference_phase_queries": 2 * PERIOD,
        }
        reference_cost_ok = all(
            counters[key] == value for key, value in reference_cost_expected.items())
        comparison_rows = {
            "helper_vs_exact_k1": tv(helper_law, exact_k1),
            "helper_vs_cached_fb_k1": tv(helper_law, fb["full_k1"]),
            "candidate_vs_exact_k1": tv(candidate, exact_k1),
            "omit_g_vs_exact_k1": tv(omit, exact_k1),
            "exact_k1_vs_cached_fb_k1": tv(exact_k1, fb["full_k1"]),
            "exact_k0_vs_cached_fb_k0": tv(exact_k0, fb["full_k0"]),
        }
        threshold_rows = {
            key: {"tv": value, "meets_1e-3": value <= 1e-3,
                  "meets_1e-2": value <= 1e-2}
            for key, value in comparison_rows.items()
        }
        work_marginal_error = tv(helper_work_law, work_probabilities)
        work_marginal_mass_error = abs(float(helper_work_law.sum()) - 1.)
        component_values = sorted(set(int(x) for x in component_counts))
        redraw_tv = tv(wrong_redraw, exact_k1)
        uniform_tv = tv(uniform, exact_k1)
        sample_attempt_total = sum(row["attempts"] for row in sample_rows)

        p1 = (valid_law(helper_law) and valid_law(exact_k1)
              and valid_law(exact_k0)
              and comparison_rows["helper_vs_exact_k1"] < TOL
              and comparison_rows["helper_vs_cached_fb_k1"] < TOL
              and comparison_rows["exact_k1_vs_cached_fb_k1"] < TOL
              and comparison_rows["exact_k0_vs_cached_fb_k0"] < TOL
              and work_marginal_mass_error < TOL
              and work_marginal_error < TOL
              and forced_counter_ok and reference_cost_ok)
        p2 = (valid_law(candidate) and valid_law(omit)
              and all(np.isfinite(v) for v in comparison_rows.values())
              and all("meets_1e-3" in row and "meets_1e-2" in row
                      for row in threshold_rows.values()))
        p3 = (len(sample_rows) == MAX_SAMPLE_CALLS
              and sum(row["counters"].get("work_draws", 0) for row in sample_rows)
              == MAX_SAMPLE_CALLS
              and sum(row["counters"].get("acceptance_draws", 0)
                      for row in sample_rows) > 0
              and sample_attempt_total <= MAX_SAMPLE_ATTEMPT_TOTAL
              and all(row["sample_counter_ok"] for row in sample_rows)
              and len(reference_sample_rows) == MAX_REFERENCE_SAMPLE_CALLS
              and all(1 <= row["steps"] <= WIDTH for row in reference_sample_rows)
              and all(1 <= row["component_count"] <= helper_stats["component_bound"]
                      for row in sample_rows)
              and helper_stats["orbit_table_entries"] == 0
              and helper_stats["exponent_table_entries"] == 0
              and helper_stats["output_table_entries"] == 0)
        c1 = uniform_tv > 1e-8
        c2 = len(component_values) > 1 and redraw_tv > 1e-8

        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "split": SPLIT, "L": L, "H": H},
            "preflight_payload_bytes": planned,
            "preflight_components": preflight_components,
            "fb_provenance_report": str(FB_REPORT),
            "helper_stats": helper_stats,
            "helper_joint_law": helper_joint.tolist(),
            "helper_output_law": helper_law.tolist(),
            "helper_work_law": helper_work_law.tolist(),
            "reference_k1_law": exact_k1.tolist(),
            "reference_k0_omit_g_law": exact_k0.tolist(),
            "fb_candidate_law": candidate.tolist(),
            "component_counts": component_counts.tolist(),
            "component_count_values": component_values,
            "work_probability_error": work_marginal_error,
            "work_marginal_mass_error": work_marginal_mass_error,
            "comparison_tvs": comparison_rows,
            "comparison_thresholds": threshold_rows,
            "wrong_redraw_law": wrong_redraw.tolist(),
            "wrong_redraw_tv": redraw_tv,
            "uniform_dephased_tv": uniform_tv,
            "sample_rows": sample_rows,
            "reference_sample_rows": reference_sample_rows,
            "sample_attempt_total": sample_attempt_total,
            "sample_attempt_budget": MAX_SAMPLE_ATTEMPT_TOTAL,
            "helper_forced_counters": helper_counter_totals,
            "helper_forced_expected_counters": helper_forced_expected,
            "forced_counter_reconciliation": forced_counter_ok,
            "reference_cost_expected": reference_cost_expected,
            "reference_cost_reconciliation": reference_cost_ok,
            "sample_counter_totals": sample_counter_totals,
            "reference_counters": {k: v for k, v in counters.items()
                                    if k not in ("helper_forced", "helper_samples")},
        })
        report["status"] = "PASS" if p1 and p2 and p3 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "forced helper agrees with independent sequential references")
        exp.check("P2", p2, "same-output candidate and omit-G TV rows are complete")
        exp.check("P3", p3, "bounded actual samples exercise helper work/rejection")
        exp.fail_check("C1", c1, "uniform dephasing differs from target")
        exp.fail_check("C2", c2, "redrawing work labels produces a biased law")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.check("P1", False, "exception before forced-law comparison")
        exp.check("P2", False, "exception before TV comparison")
        exp.check("P3", False, "exception before bounded samples")
        exp.fail_check("C1", False, "exception before dephasing control")
        exp.fail_check("C2", False, "exception before redraw control")

    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK, "t": WIDTH},
        "reference": "existing lab.semiclassical.sequential_path with supplied 60x60 work branches",
        "fb_reference": str(FB_REPORT),
        "same_output_law": True,
        "float_diagnostic": True,
        "sample_count_is_not_distribution_certificate": True,
        "no_timing_claim": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
