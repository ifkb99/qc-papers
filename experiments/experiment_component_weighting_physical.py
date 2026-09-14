"""Physical-law and actual-RNG audit of weighted progression proposals.

This is the bounded follow-up frozen in TODO 39.  It keeps the established
N=61, a=2, r=60, b=3, t=8, s=7 fixture and varies only the earlier phase
position v=0,...,7.  Existing C79 rows are compared against the independent
ER direct-column engine.  The released weighted sampler API is also exercised
with bounded draws after the algebraic checks.  This is a supplied-input
complex128 diagnostic, not a timing benchmark or finite-precision certificate.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Cycle, dual and auto rows match the independent full work/output law.
  P2  For every active row, 1 <= E=(sum sqrt(lambda))^2 <= m; weighted and
      old acceptance formulas give the same accepted law. Complete AUTO
      queries check q*A*E coordinatewise and normalize each row separately.
  P3  The weighted envelope changes expected proposal/ratio counts according
      to E; 512 actual draws finish and default seeded streams match the
      frozen prior report. No seedwise improvement or timing win is predicted.
  C1  Applying the old 1/m acceptance to the weighted proposal changes a
      physical unequal-mass row law somewhere.

The reference is built one v at a time as a 256-by-60 work-row matrix, using
the existing direct-column schedule. These dense verification arrays are
never used as a sampler oracle. No orbit lookup table is supplied.
"""
from __future__ import annotations

import math
import json
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from experiments import experiment_earlier_phase_cycles as er
from experiments.experiment_clean_orbit_output import work_block
from lab import Experiment
from lab.fourier_sampling import geometric_sum, unit_phase
from lab.work_first import EarlierPhaseProgressions


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT, Q = 8, 7, 1 << 8
INSERTIONS = tuple(range(8))
WORKS = tuple(range(PERIOD))
COVERS = ("cycle", "dual", "auto")
SELECTED_OUTPUTS = (0, 1, Q // 2, Q - 1)
PRIOR_DEFAULT_REPORT = Path("out/earlier_phase_sampler_20260911T123806945928Z.json")
BYTE_CAP = 16 * 1024 * 1024
TERM_CAP = 10_000_000
TOL = 5e-10


class AuditedEarlierPhaseProgressions(EarlierPhaseProgressions):
    """Keep the in-flight sample counters if a bounded call raises."""

    def _new_counters(self):
        counters = super()._new_counters()
        self._audit_last_counters = counters
        return counters


def report_path() -> Path:
    path = Path("out") / (
        "component_weighting_physical_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%fZ")
        + ".json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value: int, label: str) -> int:
    value = int(value)
    if value < 0 or value > BYTE_CAP:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def guard_shape(shape, dtype=np.complex128, label="array") -> int:
    return guard_bytes(math.prod(int(x) for x in shape)
                       * np.dtype(dtype).itemsize, label)


def add_count(counters: dict, key: str, amount: int,
              cap: int = TERM_CAP, label: str | None = None) -> None:
    amount = int(amount)
    if amount < 0 or counters[key] + amount > cap:
        raise MemoryError(f"{label or key} cap before loop/call")
    counters[key] += amount


def phase(index: int) -> complex:
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def preflight(prior_source_bytes: int = 0) -> dict:
    reference_calls = len(INSERTIONS) * Q
    reference_blocks = reference_calls * (PERIOD // BLOCK) * BLOCK**2
    reference_phases = reference_calls * 2 * PERIOD
    helper_rows = len(INSERTIONS) * len(COVERS) * len(WORKS)
    # Conservative structural bounds before constructing any worker.  The
    # explicit dual cover at v=0 has the largest frozen component bound.
    component_bounds = {}
    max_components = 0
    auto_components_sum = 0
    auto_components_rows_sum = 0
    auto_local_sum = 0
    auto_phase_sum = 0
    for v in INSERTIONS:
        A = 1 << v
        P = 1 if v == SPLIT else A // math.gcd(PERIOD, A)
        K = (1 << SPLIT) // A
        occupied = min(Q, 2 * 5 * math.ceil((1 << SPLIT) / PERIOD))
        cycle = min(Q, 2 * min(1 << SPLIT,
                    5 * min(P, math.ceil((1 << SPLIT) / PERIOD))), occupied)
        dual = min(Q, 2 * K * min(A, 5), occupied)
        cycle_pair = 18 * (1 + min(P, math.ceil((1 << SPLIT) / PERIOD)))
        dual_pair = 18 * K
        auto = cycle if cycle_pair <= dual_pair else dual
        auto_pair = min(cycle_pair, dual_pair)
        auto_cover = "cycle" if cycle_pair <= dual_pair else "dual"
        component_bounds[str(v)] = {"cycle": cycle, "dual": dual, "auto": auto,
                                    "cycle_pairs": cycle_pair, "dual_pairs": dual_pair,
                                    "auto_pairs": auto_pair, "auto_cover": auto_cover,
                                    "cycle_phase": 9 + 9 + 2 * 9 * min(P, math.ceil((1 << SPLIT) / PERIOD)),
                                    "dual_phase": 9 + 9 + 2 * 3 * K,
                                    "auto_phase": (9 + 9 + 2 * 9 * min(P, math.ceil((1 << SPLIT) / PERIOD))
                                                    if auto_cover == "cycle" else 9 + 9 + 2 * 3 * K)}
        max_components = max(max_components, cycle, dual)
        auto_components_sum += auto
        auto_components_rows_sum += auto * len(WORKS)
        auto_local_sum += auto_pair + 9
        auto_phase_sum += component_bounds[str(v)]["auto_phase"]
    expansion = helper_rows * Q
    selected_formula_terms = helper_rows * max_components * len(SELECTED_OUTPUTS)
    # Retain one direct reference matrix, one expanded row, FFT output, and
    # complex comparison temporaries.  Python containers/report JSON are not
    # represented as native numerical RSS.
    components = {
        "direct_work_rows": Q * PERIOD * 16,
        "direct_fft_output": Q * PERIOD * 16,
        "direct_laws_magnitudes_and_arithmetic_temps": 4 * Q * PERIOD * 8,
        "helper_expanded_row": Q * 16,
        "helper_fft_output": Q * 16,
        "work_blocks_and_temps": 12 * BLOCK * BLOCK * 16,
        "scalar_component_reserve": 2 * 1024 * 1024,
    }
    payload = guard_bytes(sum(components.values()), "component weighting preflight")
    # At most512KiB field chunks plus joined text/decoding and scalar records.
    # The full source file is streamed, not retained. Extracted numeric records
    # subsequently fit inside scalar_component_reserve above.
    prior_peak = guard_bytes(2 * 1024 * 1024, "prior sample extraction")
    operation_bounds = {
        "reference_direct_calls": reference_calls,
        "reference_block_products": reference_blocks,
        "reference_phase_queries": reference_phases,
        "reference_shift_ops": 3 * reference_calls,
        "setup_work_block_terms": 2 * BLOCK**3,
        "setup_helper_construction_terms": len(INSERTIONS) * len(COVERS) * 2 * BLOCK**3,
        "helper_rows": helper_rows,
        "helper_expansion_assignments": expansion,
        "helper_fft_calls": helper_rows,
        "helper_fft_input_entries": expansion,
        "selected_formula_component_terms": selected_formula_terms,
        "weighted_ratio_component_terms": selected_formula_terms,
        "wrong_acceptance_component_terms": selected_formula_terms,
        "wrong_full_component_terms": max_components * Q,
        "wrong_full_weighted_ratio_terms": max_components * Q,
        "mass_sqrt_terms": helper_rows * max_components,
        "helper_row_local_terms": sum(
            (component_bounds[str(v)]["cycle_pairs"] + 9
             if (cover == "cycle" or (cover == "auto"
                 and component_bounds[str(v)]["auto_cover"] == "cycle"))
             else component_bounds[str(v)]["dual_pairs"] + 9)
            * len(WORKS) for v in INSERTIONS for cover in COVERS),
        "helper_phase_queries": sum(
            (component_bounds[str(v)]["cycle_phase"] if cover == "cycle"
             else component_bounds[str(v)]["dual_phase"] if cover == "dual"
             else component_bounds[str(v)]["auto_phase"])
            * len(WORKS) for v in INSERTIONS for cover in COVERS),
        "full_auto_query_component_terms": 2 * Q * auto_components_rows_sum,
        "full_auto_query_calls": 2 * Q * PERIOD * len(INSERTIONS),
        "full_auto_query_weight_terms": 2 * auto_components_rows_sum,
        "full_auto_query_sqrt_terms": auto_components_rows_sum,
        "full_auto_query_ratio_divisions": Q * auto_components_rows_sum,
        "full_auto_query_weight_actual": 2 * auto_components_rows_sum,
        "full_auto_query_sqrt_actual": auto_components_rows_sum,
        "full_auto_query_ratio_actual": Q * auto_components_rows_sum,
        "sample_calls": 2 * len(INSERTIONS) * 32,
        "sample_attempts": 2 * len(INSERTIONS) * 32 * 512,
        "sample_component_weight_terms": 2 * 32 * auto_components_sum,
        "sample_component_sqrt_terms": 32 * auto_components_sum,
        "sample_weighted_ratio_divisions": 32 * 512 * auto_components_sum,
        "sample_fourier_component_terms": 2 * 32 * 512 * auto_components_sum,
        "sample_row_local_terms": 2 * 32 * auto_local_sum,
        "sample_column_local_terms": 2 * len(INSERTIONS) * 32 * BLOCK * BLOCK,
        # The worker's selected-cover phase bound is charged once per bounded
        # sampler call; sampled column/row work has separate named counters.
        "sample_phase_queries": 2 * 32 * auto_phase_sum,
        "sample_marginal_queries": 2 * len(INSERTIONS) * 32 * 512 * 2 * WIDTH,
        "max_components": max_components,
    }
    for key, value in operation_bounds.items():
        if value > TERM_CAP and key != "max_components":
            raise MemoryError(f"{key} preflight exceeds scalar cap")
    return {
        "planned_numeric_payload_bytes": payload,
        "prior_sample_source_bytes": int(prior_source_bytes),
        "prior_sample_extraction_peak_bytes": prior_peak,
        "peak_numeric_payload_bytes": max(payload, prior_peak),
        "payload_components": components,
        "operation_bounds": operation_bounds,
        "component_bounds_by_insertion": component_bounds,
        "max_numeric_payload_bytes": BYTE_CAP,
        "max_scalar_terms": TERM_CAP,
        "reference_is_live_direct_column": True,
        "no_Q_by_r_orbit_array": True,
    }


def extract_prior_samples(path: Path) -> dict:
    """Stream the frozen pretty-printed sample_rows field only.

    This intentionally accepts the frozen report layout rather than loading
    the 14-MiB JSON document or repeatedly reparsing a growing prefix.
    """
    if not path.is_file():
        raise FileNotFoundError(f"required frozen default-sample report missing: {path}")
    marker = '"sample_rows":'
    chunks = []
    found = False
    extracted_bytes = 0
    max_field_bytes = 512 * 1024
    with path.open("r") as source:
        for line in iter(lambda: source.readline(4097), ""):
            if len(line) > 4096:
                raise MemoryError("frozen report line exceeds bounded layout")
            if not found:
                start = line.find(marker)
                if start < 0:
                    continue
                indent = line[:start]
                if (not indent or any(char != " " for char in indent)
                        or not line.startswith(indent + marker)):
                    raise ValueError("frozen sample_rows indentation/layout mismatch")
                found = True
                suffix = line[start + len(marker):]
                opening = suffix.find("{")
                if opening < 0:
                    raise ValueError("frozen sample_rows opening object missing")
                line = suffix[opening:]
            elif line.rstrip("\r\n") in (indent + "}", indent + "},"):
                line = indent + "}"
                if extracted_bytes + len(line.encode("utf-8")) > max_field_bytes:
                    raise MemoryError("frozen sample_rows exceeds bounded field reserve")
                chunks.append(line)
                extracted_bytes += len(line.encode("utf-8"))
                break
            if extracted_bytes + len(line.encode("utf-8")) > max_field_bytes:
                raise MemoryError("frozen sample_rows extraction exceeds bounded field reserve")
            chunks.append(line)
            extracted_bytes += len(line.encode("utf-8"))
    if not found:
        raise ValueError("frozen report has no sample_rows field")
    if not chunks or not line.rstrip("\r\n").startswith(indent + "}"):
        raise ValueError("frozen sample_rows field is incomplete")
    samples = json.loads("".join(chunks))
    del chunks
    if not isinstance(samples, dict) or set(samples) != {str(v) for v in INSERTIONS}:
        raise ValueError("frozen sample_rows insertion set mismatch")
    if any(len(samples[str(v)]) != 32 for v in INSERTIONS):
        raise ValueError("frozen default sample count mismatch")
    return samples


def expand_row(row: dict, counters: dict, label: str) -> np.ndarray:
    guard_shape((Q,), label=f"{label} expanded row")
    expanded = np.zeros(Q, dtype=np.complex128)
    for start, count, gamma in row["components"]:
        for n in range(int(count)):
            add_count(counters, "helper_expansion_assignments", 1,
                      counters["helper_expansion_cap"], label="row expansion")
            exponent = int(start) + int(row["stride"]) * n
            if not 0 <= exponent < Q or expanded[exponent] != 0:
                raise AssertionError(f"invalid/colliding progression in {label}")
            expanded[exponent] = gamma
    add_count(counters, "helper_fft_calls", 1,
              counters["helper_fft_call_cap"], label="row FFT")
    add_count(counters, "helper_fft_input_entries", Q,
              counters["helper_fft_input_cap"], label="row FFT input")
    return expanded


def component_values(row: dict, output: int, counters: dict,
                     counter_key: str = "selected_formula_component_terms") -> list[complex]:
    values = []
    for start, count, gamma in row["components"]:
        cap_key = "selected_formula_cap" if counter_key == "selected_formula_component_terms" \
            else "wrong_full_cap"
        add_count(counters, counter_key, 1, counters[cap_key], label="formula component")
        values.append(complex(gamma) * unit_phase(-int(output) * int(start), Q)
                      * geometric_sum(int(count), -int(output) * int(row["stride"]), Q))
    return values


def reference_rows(insertion: int, W0: np.ndarray, W1: np.ndarray,
                   counters: dict) -> tuple[np.ndarray, np.ndarray]:
    guard_shape((Q, PERIOD), label=f"reference v={insertion} rows")
    rows = np.zeros((Q, PERIOD), dtype=np.complex128)
    for exponent in range(Q):
        if counters["reference_direct_calls"] + 1 > counters["reference_call_cap"]:
            raise MemoryError("reference call cap before direct_column")
        block_terms = (PERIOD // BLOCK) * BLOCK**2
        if counters["reference_block_products"] + block_terms > counters["reference_block_cap"]:
            raise MemoryError("reference block cap before direct_column")
        add_count(counters, "reference_direct_calls", 1,
                  counters["reference_call_cap"], label="direct column")
        add_count(counters, "reference_block_products", block_terms,
                  counters["reference_block_cap"], label="direct block")
        add_count(counters, "reference_shift_ops", 3,
                  counters["reference_shift_cap"], label="direct shift")
        # er.direct_column receives the same cumulative dictionary, so its
        # phase/block guards remain cumulative across all eight positions.
        if (counters["direct_phase_queries"] + 2 * PERIOD > er.MAX_TERMS
                or counters["direct_block_products"] + block_terms > er.MAX_TERMS):
            raise MemoryError("ER cumulative cap before direct_column")
        vector = er.direct_column(exponent, insertion, W0, W1, counters,
                                  split=SPLIT)
        rows[exponent, :] = vector
    add_count(counters, "reference_fft_calls", 1,
              counters["reference_fft_call_cap"], label="reference FFT")
    add_count(counters, "reference_fft_input_entries", Q * PERIOD,
              counters["reference_fft_input_cap"], label="reference FFT input")
    raw = np.fft.fft(rows, axis=0)
    norms = np.sum(np.abs(rows)**2, axis=0)
    laws = np.abs(raw)**2 / (Q * norms[None, :])
    return rows, laws


def main() -> None:
    exp = Experiment("component_weighting_physical", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "physical cycle/dual/auto rows match live direct-column laws")
    exp.predict("P2", "weighted envelope bounds and accepted-law identities hold")
    exp.predict("P3", "expected named proposal/ratio work follows E without timing claim")
    exp.must_fail("C1", "old acceptance on weighted proposals changes an unequal row")
    started = time.perf_counter()
    counters = {
        "reference_direct_calls": 0, "reference_block_products": 0,
        "reference_phase_queries": 0, "reference_fft_calls": 0,
        "reference_fft_input_entries": 0, "reference_shift_ops": 0,
        "direct_phase_queries": 0, "modular_pow_queries": 0,
        "direct_block_products": 0, "direct_shift_ops": 0,
        "setup_work_block_terms": 0, "setup_helper_construction_terms": 0,
        "helper_rows": 0, "helper_row_local_terms": 0,
        "helper_early_phase_queries": 0, "helper_late_phase_queries": 0,
        "helper_expansion_assignments": 0, "helper_fft_calls": 0,
        "helper_fft_input_entries": 0, "selected_formula_component_terms": 0,
        "mass_sqrt_terms": 0, "weighted_ratio_terms": 0,
        "wrong_ratio_terms": 0, "wrong_weighted_ratio_terms": 0,
        "wrong_full_component_terms": 0,
        "full_auto_query_calls": 0, "full_auto_query_component_terms": 0,
        "full_auto_query_weight_terms": 0, "full_auto_query_sqrt_terms": 0,
        "full_auto_query_ratio_divisions": 0,
        "full_auto_query_weight_actual": 0,
        "full_auto_query_sqrt_actual": 0,
        "full_auto_query_ratio_actual": 0,
        "sample_calls": 0, "sample_attempts": 0,
        "sample_fourier_component_terms": 0, "sample_column_local_terms": 0,
        "sample_row_local_terms": 0, "sample_phase_queries": 0,
        "sample_marginal_queries": 0, "sample_component_weight_terms": 0,
        "sample_component_sqrt_terms": 0, "sample_weighted_ratio_divisions": 0,
    }
    report = {"status": "FAIL", "rows": [], "checks": {},
              "predictions": {"insertions": list(INSERTIONS),
                              "covers": list(COVERS), "works": "all 60"}}
    p1 = p2 = p3 = c1 = False
    # Initialize partial-result containers before preflight so an early
    # resource failure still produces a truthful failure report.
    actual_sample_rows = {}
    full_auto_report = {}
    prior_default_samples = None
    try:
        report["fixture"] = {"N": N, "base": BASE, "period": PERIOD,
                             "block_size": BLOCK, "width": WIDTH,
                             "split": SPLIT, "insertions": list(INSERTIONS),
                             "covers": list(COVERS)}
        if not PRIOR_DEFAULT_REPORT.is_file():
            raise FileNotFoundError(f"required frozen default report missing: {PRIOR_DEFAULT_REPORT}")
        report["preflight"] = preflight(PRIOR_DEFAULT_REPORT.stat().st_size)
        bounds = report["preflight"]["operation_bounds"]
        counters.update({
            "reference_call_cap": bounds["reference_direct_calls"],
            "reference_block_cap": bounds["reference_block_products"],
            "reference_shift_cap": bounds["reference_shift_ops"],
            "setup_work_block_cap": bounds["setup_work_block_terms"],
            "setup_helper_cap": bounds["setup_helper_construction_terms"],
            "reference_fft_call_cap": len(INSERTIONS),
            "reference_fft_input_cap": len(INSERTIONS) * Q * PERIOD,
            "helper_row_cap": bounds["helper_rows"],
            "helper_expansion_cap": bounds["helper_expansion_assignments"],
            "helper_fft_call_cap": bounds["helper_fft_calls"],
            "helper_fft_input_cap": bounds["helper_fft_input_entries"],
            "selected_formula_cap": bounds["selected_formula_component_terms"],
            "weighted_ratio_cap": bounds["weighted_ratio_component_terms"],
            "wrong_ratio_cap": bounds["wrong_acceptance_component_terms"],
            "wrong_full_cap": bounds["wrong_full_component_terms"],
            "wrong_full_weighted_ratio_cap": bounds["wrong_full_weighted_ratio_terms"],
            "mass_sqrt_cap": bounds["mass_sqrt_terms"],
            "helper_row_local_cap": bounds["helper_row_local_terms"],
            "helper_phase_cap": bounds["helper_phase_queries"],
            "full_auto_query_call_cap": bounds["full_auto_query_calls"],
            "full_auto_query_component_cap": bounds["full_auto_query_component_terms"],
            "full_auto_query_weight_cap": bounds["full_auto_query_weight_terms"],
            "full_auto_query_sqrt_cap": bounds["full_auto_query_sqrt_terms"],
            "full_auto_query_division_cap": bounds["full_auto_query_ratio_divisions"],
            "full_auto_query_weight_actual_cap": bounds["full_auto_query_weight_actual"],
            "full_auto_query_sqrt_actual_cap": bounds["full_auto_query_sqrt_actual"],
            "full_auto_query_ratio_actual_cap": bounds["full_auto_query_ratio_actual"],
            "sample_call_cap": bounds["sample_calls"],
            "sample_attempt_cap": bounds["sample_attempts"],
            "sample_fourier_cap": bounds["sample_fourier_component_terms"],
            "sample_row_local_cap": bounds["sample_row_local_terms"],
            "sample_column_local_cap": bounds["sample_column_local_terms"],
            "sample_phase_cap": bounds["sample_phase_queries"],
            "sample_marginal_cap": bounds["sample_marginal_queries"],
            "sample_weight_cap": bounds["sample_component_weight_terms"],
            "sample_sqrt_cap": bounds["sample_component_sqrt_terms"],
            "sample_division_cap": bounds["sample_weighted_ratio_divisions"],
        })
        # Extract only the bounded frozen sample field before numerical setup.
        prior_default_samples = extract_prior_samples(PRIOR_DEFAULT_REPORT)
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        add_count(counters, "setup_work_block_terms", 2 * BLOCK**3,
                  counters["setup_work_block_cap"], label="work-block setup")
        if not (np.allclose(W0.conj().T @ W0, np.eye(BLOCK), atol=1e-12, rtol=0)
                and np.allclose(W1.conj().T @ W1, np.eye(BLOCK), atol=1e-12, rtol=0)):
            raise AssertionError("supplied work_block matrices are not unitary")
        rows_report = {}
        max_law_tv = max_amp_error = max_work_error = 0.
        max_weighted_error = max_mass_error = 0.
        max_wrong_error = 0.
        wrong_control_tv = 0.
        wrong_control_example = None
        unequal_rows = 0
        for insertion in INSERTIONS:
            reference, reference_laws = reference_rows(insertion, W0, W1, counters)
            insertion_rows = {}
            workers = {}
            stats = {}
            auto_rows = {}
            for cover in COVERS:
                add_count(counters, "setup_helper_construction_terms", 2 * BLOCK**3,
                          counters["setup_helper_cap"], label="helper construction")
                worker = AuditedEarlierPhaseProgressions(
                    PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, phase,
                    early_split=insertion, early_phase=phase, cover=cover,
                    max_local_terms=1_000_000, max_components=4096,
                    max_payload_bytes=BYTE_CAP)
                workers[cover] = worker
                stats[cover] = json_safe(worker.stats())
            for cover, worker in workers.items():
                insertion_rows[cover] = {}
                for work in WORKS:
                    add_count(counters, "helper_rows", 1, counters["helper_row_cap"],
                              label="helper row")
                    row_local_bound = int(stats[cover]["local_term_bound"])
                    row_phase_bound = int(stats[cover]["phase_queries_per_draw_bound"])
                    if (counters["helper_row_local_terms"] + row_local_bound
                            > counters["helper_row_local_cap"]
                            or counters["helper_early_phase_queries"]
                            + counters["helper_late_phase_queries"]
                            + row_phase_bound > counters["helper_phase_cap"]):
                        raise MemoryError("helper row local/phase cap before row()")
                    row = worker.row(work)
                    row_counters = row["counters"]
                    counters["helper_row_local_terms"] += int(row_counters.get("row_local_terms", 0))
                    counters["helper_early_phase_queries"] += int(row_counters.get("early_phase_queries", 0))
                    counters["helper_late_phase_queries"] += int(row_counters.get("late_phase_queries", 0))
                    expanded = expand_row(row, counters,
                                          f"v={insertion},cover={cover},work={work}")
                    helper_raw = np.fft.fft(expanded)
                    helper_norm = float(np.sum(np.abs(expanded)**2))
                    helper_law = np.abs(helper_raw)**2 / (Q * helper_norm)
                    target_law = reference_laws[:, work]
                    amp_error = float(np.max(np.abs(expanded - reference[:, work])))
                    law_error = float(np.sum(np.abs(helper_law - target_law)) / 2)
                    work_error = abs(helper_norm / Q
                                     - float(np.sum(np.abs(reference[:, work])**2) / Q))
                    max_law_tv = max(max_law_tv, law_error)
                    max_amp_error = max(max_amp_error, amp_error)
                    max_work_error = max(max_work_error, work_error)
                    lambdas = [float(count) * abs(complex(gamma))**2 / float(row["norm"])
                               for _, count, gamma in row["components"]]
                    add_count(counters, "mass_sqrt_terms", len(lambdas),
                              counters["mass_sqrt_cap"], label="mass square roots")
                    lambda_sum = math.fsum(lambdas)
                    roots = [math.sqrt(x) for x in lambdas]
                    Bsum = math.fsum(roots)
                    E = Bsum * Bsum
                    m = len(lambdas)
                    mass_error = abs(lambda_sum - 1.)
                    max_mass_error = max(max_mass_error, mass_error)
                    if m and not (1. - TOL <= E <= m + TOL):
                        raise AssertionError("weighted envelope outside [1,m]")
                    if len({round(x, 12) for x in lambdas}) > 1:
                        unequal_rows += 1
                    weighted_error = 0.
                    wrong_ratios = []
                    add_count(counters, "weighted_ratio_terms",
                              len(lambdas) * len(SELECTED_OUTPUTS),
                              counters["weighted_ratio_cap"],
                              label="weighted selected ratios")
                    for output in SELECTED_OUTPUTS:
                        values = component_values(row, output, counters)
                        add_count(counters, "wrong_ratio_terms", len(values),
                                  counters["wrong_ratio_cap"], label="wrong ratio")
                        coherent = abs(sum(values, 0j))**2
                        old_square = math.fsum(abs(x)**2 for x in values)
                        weighted_square = math.fsum(
                            abs(x)**2 / root
                            for x, root in zip(values, roots))
                        target = float(target_law[output])
                        accepted_wrong = 0.
                        if weighted_square == 0:
                            accepted_weighted = 0.
                        else:
                            # q*A has total mass 1/B^2.  Rescale by B^2
                            # before comparing the accepted conditional law.
                            accepted_weighted = (weighted_square / (Q * float(row["norm"]) * Bsum)
                                                 * coherent / (Bsum * weighted_square)
                                                 * Bsum * Bsum)
                            # Deliberately wrong: old 1/m envelope paired with
                            # the new weighted proposal.
                            accepted_wrong = (weighted_square / (Q * float(row["norm"]) * Bsum)
                                              * coherent / (m * old_square)) if old_square else 0.
                        weighted_error = max(weighted_error, abs(accepted_weighted - target))
                        if target > TOL:
                            # A wrong accepted sublaw can have a different
                            # total mass; compare its shape by ratios, which
                            # is invariant under the eventual normalization.
                            wrong_ratios.append(accepted_wrong / target)
                    wrong_error = (max(wrong_ratios) - min(wrong_ratios)
                                   if wrong_ratios else 0.)
                    if unequal_rows == 1 and wrong_control_example is None:
                        wrong_sublaw = np.zeros(Q, dtype=np.float64)
                        add_count(counters, "wrong_weighted_ratio_terms", len(lambdas) * Q,
                                  counters["wrong_full_weighted_ratio_cap"],
                                  label="wrong full weighted ratios")
                        for output in range(Q):
                            values = component_values(
                                row, output, counters,
                                counter_key="wrong_full_component_terms")
                            weighted_square = math.fsum(
                                abs(x)**2 / root
                                for x, root in zip(values, roots))
                            old_square = math.fsum(abs(x)**2 for x in values)
                            coherent = abs(sum(values, 0j))**2
                            if weighted_square and old_square:
                                wrong_sublaw[output] = (
                                    weighted_square / (Q * float(row["norm"]) * Bsum)
                                    * coherent / (m * old_square))
                        wrong_mass = float(np.sum(wrong_sublaw))
                        if wrong_mass > 0:
                            wrong_normalized = wrong_sublaw / wrong_mass
                            wrong_control_tv = float(
                                np.sum(np.abs(wrong_normalized - target_law)) / 2.)
                            wrong_control_example = {
                                "insertion": insertion, "work": work,
                                "component_count": m, "wrong_mass": wrong_mass,
                                "tv": wrong_control_tv,
                            }
                    max_weighted_error = max(max_weighted_error, weighted_error)
                    max_wrong_error = max(max_wrong_error, wrong_error)
                    insertion_rows[cover][str(work)] = {
                        "component_count": m, "norm": float(row["norm"]),
                        "stride": int(row["stride"]), "lambda_sum_error": mass_error,
                        "work_probability": float(helper_norm / Q),
                        "envelope_E": E, "old_mean_attempts": m,
                        "weighted_mean_attempts": E,
                        "expected_old_ratio_terms": float(m * m),
                        "expected_weighted_ratio_terms": float(E * m),
                        "law_tv_vs_direct": law_error,
                        "amplitude_max_error": amp_error,
                        "work_mass_error": work_error,
                        "weighted_selected_accepted_error": weighted_error,
                        "wrong_old_acceptance_error": wrong_error,
                        "law_normalized": abs(float(np.sum(helper_law)) - 1.) < TOL,
                    }
                    if cover == "auto":
                        auto_rows[str(work)] = row
            # Complete AUTO conditional-law query through the released
            # prepared-row API.  Both proposal modes are checked against the
            # same independent direct-column law; no probability table is
            # retained after this insertion.
            full_mode_checks = {}
            for mode in ("mass", "root_mass"):
                max_error = 0.
                conditional_mass = 0.
                proposal_mass = 0.
                expected_values = []
                row_checks = []
                for work in WORKS:
                    row = auto_rows[str(work)]
                    local = workers["auto"]._new_counters()
                    m = len(row["components"])
                    # Reserve all preparation and ratio operations before the
                    # corresponding helper calls.  The separate *_actual
                    # fields below reconcile returned counters afterward.
                    add_count(counters, "full_auto_query_weight_terms", m,
                              counters["full_auto_query_weight_cap"],
                              label="full query weights")
                    if mode == "root_mass":
                        add_count(counters, "full_auto_query_sqrt_terms", m,
                                  counters["full_auto_query_sqrt_cap"],
                                  label="full query square roots")
                    prepared = workers["auto"]._prepare_proposal(row, mode, local)
                    expected_values.append(float(prepared["expected_attempts"]))
                    if mode == "root_mass":
                        add_count(counters, "full_auto_query_ratio_divisions", Q * m,
                                  counters["full_auto_query_division_cap"],
                                  label="full query ratio divisions")
                    row_conditional_mass = 0.
                    row_proposal_mass = 0.
                    row_accepted_mass = 0.
                    row_coord_error = 0.
                    for output in range(Q):
                        add_count(counters, "full_auto_query_calls", 1,
                                  counters["full_auto_query_call_cap"], label="full query")
                        add_count(counters, "full_auto_query_component_terms", m,
                                  counters["full_auto_query_component_cap"], label="full query components")
                        result = workers["auto"]._fourier(
                            row, output, local, prepared=prepared)
                        target = float(reference_laws[output, work])
                        conditional = float(result["conditional_probability"])
                        proposal = float(result["proposal_probability"])
                        acceptance = float(result["acceptance"])
                        expected = float(result["expected_attempts"])
                        max_error = max(max_error,
                                        abs(conditional - target))
                        conditional_mass += conditional
                        proposal_mass += proposal
                        row_conditional_mass += conditional
                        row_proposal_mass += proposal
                        accepted_coordinate = proposal * acceptance * expected
                        row_accepted_mass += accepted_coordinate
                        row_coord_error = max(row_coord_error,
                                              abs(accepted_coordinate - target))
                    actual_weight = int(local.get("component_weight_terms", 0))
                    actual_sqrt = int(local.get("component_sqrt_terms", 0))
                    actual_ratio = int(local.get("weighted_ratio_divisions", 0))
                    if (actual_weight != m
                            or actual_sqrt != (m if mode == "root_mass" else 0)
                            or actual_ratio != (Q*m if mode == "root_mass" else 0)
                            or local["fourier_component_terms"] != Q*m):
                        raise AssertionError("full query observed counters do not reconcile")
                    counters["full_auto_query_weight_actual"] += actual_weight
                    counters["full_auto_query_sqrt_actual"] += actual_sqrt
                    counters["full_auto_query_ratio_actual"] += actual_ratio
                    if (counters["full_auto_query_weight_actual"]
                            > counters["full_auto_query_weight_actual_cap"]
                            or counters["full_auto_query_sqrt_actual"]
                            > counters["full_auto_query_sqrt_actual_cap"]
                            or counters["full_auto_query_ratio_actual"]
                            > counters["full_auto_query_ratio_actual_cap"]):
                        raise MemoryError("full query actual counter cap")
                    row_checks.append({
                        "work": work,
                        "conditional_mass_error": abs(row_conditional_mass - 1.),
                        "proposal_mass_error": abs(row_proposal_mass - 1.),
                        "accepted_mass_error": abs(row_accepted_mass - 1.),
                        "accepted_coordinate_error": row_coord_error,
                    })
                full_mode_checks[mode] = {
                    "max_conditional_error": max_error,
                    "conditional_mass_error": abs(conditional_mass / PERIOD - 1.),
                    "proposal_mass_error": abs(proposal_mass / PERIOD - 1.),
                    "expected_attempts_min": min(expected_values),
                    "expected_attempts_max": max(expected_values),
                    "rows": row_checks,
                }
            full_auto_report[str(insertion)] = full_mode_checks

            # Actual bounded RNG comparison: 32 draws per insertion and mode,
            # with the historical default mass stream retained exactly.
            actual_sample_rows[str(insertion)] = {}
            for mode in ("mass", "root_mass"):
                # Reset both proposal modes to the same frozen stream; mode
                # changes only the proposal computation, not RNG provenance.
                seed = 624 + insertion
                rng = np.random.default_rng(seed)
                draws = []
                actual_sample_rows[str(insertion)][mode] = draws
                for draw_index in range(32):
                    auto_bound = int(stats["auto"]["component_bound"])
                    row_local_bound = int(stats["auto"]["local_term_bound"])
                    phase_bound = int(stats["auto"]["phase_queries_per_draw_bound"])
                    if counters["sample_attempts"] + 512 > counters["sample_attempt_cap"]:
                        raise MemoryError("sample attempt cap before sample()")
                    if counters["sample_calls"] + 1 > counters["sample_call_cap"]:
                        raise MemoryError("sample call cap before sample()")
                    for key, amount, cap_key in (
                            ("sample_fourier_component_terms", 512 * auto_bound, "sample_fourier_cap"),
                            ("sample_component_weight_terms", auto_bound, "sample_weight_cap"),
                            ("sample_component_sqrt_terms",
                             auto_bound if mode == "root_mass" else 0, "sample_sqrt_cap"),
                            ("sample_weighted_ratio_divisions",
                             512 * auto_bound if mode == "root_mass" else 0, "sample_division_cap"),
                            ("sample_row_local_terms", row_local_bound, "sample_row_local_cap"),
                            ("sample_column_local_terms", BLOCK * BLOCK, "sample_column_local_cap"),
                            ("sample_phase_queries", phase_bound, "sample_phase_cap"),
                            ("sample_marginal_queries", 2 * WIDTH * 512, "sample_marginal_cap")):
                        if counters[key] + amount > counters[cap_key]:
                            raise MemoryError(f"{key} cap before sample()")
                    counters["sample_calls"] += 1
                    sampler = workers["auto"]
                    sampler._audit_last_counters = None
                    try:
                        sample = sampler.sample(rng, max_attempts=512, proposal=mode)
                    except Exception as sample_exc:
                        observed = getattr(sampler, "_audit_last_counters", {}) or {}
                        for key, counter_key in (
                                ("sample_attempts", "component_draws"),
                                ("sample_fourier_component_terms", "fourier_component_terms"),
                                ("sample_column_local_terms", "column_local_terms"),
                                ("sample_row_local_terms", "row_local_terms"),
                                ("sample_phase_queries", "early_phase_queries"),
                                ("sample_marginal_queries", "marginal_queries"),
                                ("sample_component_weight_terms", "component_weight_terms"),
                                ("sample_component_sqrt_terms", "component_sqrt_terms"),
                                ("sample_weighted_ratio_divisions", "weighted_ratio_divisions")):
                            counters[key] += int(observed.get(counter_key, 0))
                        counters["sample_phase_queries"] += int(
                            observed.get("late_phase_queries", 0))
                        draws.append({"draw": draw_index, "failed": True,
                                      "exception": repr(sample_exc),
                                      "counters": json_safe(observed)})
                        raise
                    observed = sample["counters"]
                    counters["sample_attempts"] += int(sample["attempts"])
                    counters["sample_fourier_component_terms"] += int(observed.get("fourier_component_terms", 0))
                    counters["sample_column_local_terms"] += int(observed.get("column_local_terms", 0))
                    counters["sample_row_local_terms"] += int(observed.get("row_local_terms", 0))
                    counters["sample_phase_queries"] += int(observed.get("early_phase_queries", 0))
                    counters["sample_phase_queries"] += int(observed.get("late_phase_queries", 0))
                    counters["sample_marginal_queries"] += int(observed.get("marginal_queries", 0))
                    counters["sample_component_weight_terms"] += int(observed.get("component_weight_terms", 0))
                    counters["sample_component_sqrt_terms"] += int(observed.get("component_sqrt_terms", 0))
                    counters["sample_weighted_ratio_divisions"] += int(observed.get("weighted_ratio_divisions", 0))
                    for key, cap_key in (
                            ("sample_attempts", "sample_attempt_cap"),
                            ("sample_fourier_component_terms", "sample_fourier_cap"),
                            ("sample_column_local_terms", "sample_column_local_cap"),
                            ("sample_row_local_terms", "sample_row_local_cap"),
                            ("sample_phase_queries", "sample_phase_cap"),
                            ("sample_marginal_queries", "sample_marginal_cap"),
                            ("sample_component_weight_terms", "sample_weight_cap"),
                            ("sample_component_sqrt_terms", "sample_sqrt_cap"),
                            ("sample_weighted_ratio_divisions", "sample_division_cap")):
                        if counters[key] > counters[cap_key]:
                            raise MemoryError(f"{key} cap after sample()")
                    prior = (prior_default_samples.get(str(insertion), [])[draw_index]
                             if mode == "mass" and prior_default_samples else None)
                    draws.append({
                        "draw": draw_index, "seed": seed, "work": int(sample["work"]),
                        "output": int(sample["output"]), "attempts": int(sample["attempts"]),
                        "component_count": int(sample["component_count"]),
                        "proposal_mode": sample.get("proposal_mode"),
                        "expected_attempts": float(sample.get("expected_attempts", 0.)),
                        "counters": json_safe(observed),
                        "default_seed_parity": (None if prior is None else
                            all(int(sample[key]) == int(prior[key])
                                for key in ("work", "output", "attempts", "component_count"))),
                    })

            rows_report[str(insertion)] = {
                "stats": stats, "rows": insertion_rows,
                "max_components": max(insertion_rows[c][str(w)]["component_count"]
                                       for c in COVERS for w in WORKS),
                "full_auto_queries": full_mode_checks,
                "samples": actual_sample_rows[str(insertion)],
            }
            # target_law is a view and would otherwise retain the old laws
            # while the next insertion's rows/FFT temporaries are allocated.
            del target_law, reference, reference_laws
        cost_by_insertion_cover = {}
        for insertion, insertion_data in rows_report.items():
            cost_by_insertion_cover[insertion] = {}
            for cover in COVERS:
                items = tuple(insertion_data["rows"][cover].values())
                cost_by_insertion_cover[insertion][cover] = {
                    "work_mass_sum": math.fsum(x["work_probability"] for x in items),
                    "expected_old_attempts": math.fsum(
                        x["work_probability"] * x["old_mean_attempts"] for x in items),
                    "expected_weighted_attempts": math.fsum(
                        x["work_probability"] * x["weighted_mean_attempts"] for x in items),
                    "expected_old_ratio_terms": math.fsum(
                        x["work_probability"] * x["expected_old_ratio_terms"] for x in items),
                    "expected_weighted_ratio_terms": math.fsum(
                        x["work_probability"] * x["expected_weighted_ratio_terms"] for x in items),
                    "row_mass_sqrt_setup_terms": math.fsum(
                        x["component_count"] for x in items),
                }
        p1 = (max_law_tv < TOL and max_amp_error < TOL
              and max_work_error < TOL
              and all(item["law_normalized"]
                      for ins in rows_report.values() for cov in ins["rows"].values()
                      for item in cov.values()))
        full_query_ok = all(
            item["max_conditional_error"] < TOL
            and item["conditional_mass_error"] < TOL
            and item["proposal_mass_error"] < TOL
            and all(row["conditional_mass_error"] < TOL
                    and row["proposal_mass_error"] < TOL
                    and row["accepted_mass_error"] < TOL
                    and row["accepted_coordinate_error"] < 2*TOL
                    for row in item["rows"])
            for insertion in full_auto_report.values()
            for item in insertion.values())
        p2 = (max_mass_error < TOL and max_weighted_error < 2*TOL
              and full_query_ok)
        # This is an analytical expected-work result, not a timing assertion.
        expected_cost_ok = all(1. - TOL <= item["weighted_mean_attempts"]
                 <= item["old_mean_attempts"] + TOL
                 for ins in rows_report.values() for cov in ins["rows"].values()
                 for item in cov.values())
        actual_flat = [item for modes in actual_sample_rows.values()
                       for draws in modes.values() for item in draws]
        actual_cost_by_mode = {}
        cost_fields = (
            "attempts", "component_count", "fourier_component_terms",
            "column_local_terms", "row_local_terms", "early_phase_queries",
            "late_phase_queries", "marginal_queries", "component_weight_terms",
            "component_sqrt_terms", "weighted_ratio_divisions")
        for mode in ("mass", "root_mass"):
            draws = [item for modes in actual_sample_rows.values()
                     for values in modes.items() if values[0] == mode
                     for item in values[1]]
            totals = {key: 0 for key in cost_fields}
            for item in draws:
                totals["attempts"] += int(item["attempts"])
                totals["component_count"] += int(item["component_count"])
                for key in cost_fields[2:]:
                    totals[key] += int(item["counters"].get(key, 0))
            actual_cost_by_mode[mode] = {
                "draw_count": len(draws), "totals": totals,
                "mean_attempts": (totals["attempts"] / len(draws) if draws else None),
            }
        parity_ok = all(item["default_seed_parity"] is True
                        for item in actual_flat if item["proposal_mode"] == "mass")
        actual_ok = (len(actual_flat) == 2 * len(INSERTIONS) * 32
                     and all(0 <= item["work"] < PERIOD and 0 <= item["output"] < Q
                             and 1 <= item["attempts"] <= 512
                             and item["proposal_mode"] in ("mass", "root_mass")
                             for item in actual_flat))
        p3 = expected_cost_ok and actual_ok and parity_ok
        c1 = unequal_rows > 0 and wrong_control_tv > 1e-8
        counters["reference_phase_queries"] = counters["direct_phase_queries"]
        report.update({
            "rows": rows_report,
            "full_auto_query_checks": full_auto_report,
            "actual_sample_rows": actual_sample_rows,
            "prior_default_report": str(PRIOR_DEFAULT_REPORT),
            "counters": counters,
            "cost_summary": {
                "max_law_tv": max_law_tv, "max_amplitude_error": max_amp_error,
                "max_work_mass_error": max_work_error,
                "max_weighted_selected_error": max_weighted_error,
                "max_wrong_acceptance_error": max_wrong_error,
                "wrong_control_tv": wrong_control_tv,
                "wrong_control_example": wrong_control_example,
                "max_lambda_sum_error": max_mass_error,
                "unequal_mass_rows": unequal_rows,
                "weighted_sampler_api_executed": True,
                "expected_work_only": False,
                "by_insertion_cover": cost_by_insertion_cover,
                "full_auto_query_checks": full_auto_report,
                "actual_samples_executed": len(actual_flat),
                "default_seed_parity": parity_ok,
                "actual_cost_by_mode": actual_cost_by_mode,
            },
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
            "python_version": sys.version, "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "all physical rows match direct-column laws")
        exp.check("P2", p2, "weighted envelope and accepted-law identities hold")
        exp.check("P3", p3, "expected named work follows weighted envelope")
        exp.fail_check("C1", c1, "wrong old acceptance fails on unequal masses")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        report["counters"] = counters
        report["actual_sample_rows_partial"] = actual_sample_rows
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before physical weighting audit")
        exp.fail_check("C1", False, "exception before unequal-mass control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": report.get("fixture"),
        "weighted_sampler_api_executed": True,
        "reference": "live ER direct_column, one 256x60 matrix per insertion",
        "no_timing_claim": True, "max_numeric_payload_bytes": BYTE_CAP,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
