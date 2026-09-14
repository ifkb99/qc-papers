"""Bounded long-progression discriminator for the nested phase cut.

This is the first post-tiny C81 amplitude check.  It keeps the physical
fixture fixed at r=60,b=3,t=13,s=12 and compares cuts 5, 11, and 12 on work
labels 0, 1, and 59.  The existing experiment-local ``cut_row`` is reused;
the literal reference is the existing ``direct_column`` engine, queried only
on the seven residues permitted by the [-2,4] support cone.  No generic
propagator, sampler, or Q-by-r reference matrix is constructed.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every expanded cut row agrees with the sparse literal reference at all
      8192 exponent coordinates, including its complete normalized FFT law.
  P2  The four explicit geometric queries (Q=8192) agree with the FFT law;
      all literal columns and positive row norms are normalized.
  P3  The frozen per-category reservations and actual per-cut ledgers remain
      bounded; component/cost ratios are descriptive, not timing claims.
  C1  Freezing the left phase coefficient at z=0 for cut 11 changes a named
      amplitude at the frozen named work-0 row.

The reported overlapping checksum is a ledger bound, not an additive FLOP
count or RSS measurement.  Exact zero is the only coefficient omission rule.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.fourier_sampling import geometric_sum, unit_phase
from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_earlier_phase_cycles import direct_column
from experiments.experiment_nested_phase_amplitudes import (
    charge, cut_row, expand_row, g_phase)


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT, Q, L, H = 13, 12, 1 << 13, 1 << 12, 2
POSITIONS = (5, 11)
CUTS = (5, 11, 12)
WORKS = (0, 1, 59)
DELTA_SUPPORT = tuple(range(-2, 5))
MAX_BYTES = 4 << 20
MAX_CATEGORY_TERMS = 5_000_000
MAX_TOTAL_RESERVED = 15_000_000
TOL = 3e-10

# Frozen conservative per-row formula bounds from TODO40.  These are
# separate from the support-reference and comparison ledgers below.
ROW_BOUNDS = {
    5: {"group_pair_visits": 2304, "coefficient_pair_visits": 2304,
        "components": 690, "phase_calls": 4610,
        "phase_products": 7602},
    11: {"group_pair_visits": 36, "coefficient_pair_visits": 288,
         "components": 160, "phase_calls": 578,
         "phase_products": 1024},
    12: {"group_pair_visits": 18, "coefficient_pair_visits": 1242,
         "components": 690, "phase_calls": 2486,
         "phase_products": 4416},
}


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_long_rows_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def json_safe(value):
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    return value


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value}>{MAX_BYTES}")


def reserve(ledger, key, amount=1, *, cap=MAX_CATEGORY_TERMS):
    amount = int(amount)
    if ledger.get(key, 0) + amount > cap:
        raise MemoryError(f"{key} reservation exceeds cap before operation")
    ledger[key] = ledger.get(key, 0) + amount


def support_count(work):
    """Count the seven arithmetic progressions without allocating a list."""
    total = 0
    for delta in DELTA_SUPPORT:
        residue = (int(work) - int(delta)) % PERIOD
        total += 0 if residue >= Q else 1 + (Q - 1 - residue) // PERIOD
    return total


def support_exponents(work, counters=None):
    """Enumerate only the seven allowed residue progressions."""
    result = []
    for delta in DELTA_SUPPORT:
        residue = (int(work) - int(delta)) % PERIOD
        count = 0 if residue >= Q else 1 + (Q - 1 - residue) // PERIOD
        for offset in range(count):
            if counters is not None:
                charge(counters, "support_generation_visits")
            result.append(residue + PERIOD * offset)
    return result


def selected_geometric_long(row, output, counters):
    total = 0j
    for start, count, gamma in row["components"]:
        charge(counters, "fourier_terms")
        total += (gamma * unit_phase(-int(output) * int(start), Q)
                  * geometric_sum(int(count), -int(output) * int(row["stride"]), Q))
    return total


def preflight():
    """Reserve all arrays and named work before the first science call."""
    support_counts = {str(work): support_count(work) for work in WORKS}
    reference_calls = sum(support_counts.values())
    if reference_calls > 2877:
        raise MemoryError("support-cone reference-call bound exceeded")

    formula_bounds = defaultdict(int)
    for cut in CUTS:
        for key, amount in ROW_BOUNDS[cut].items():
            formula_bounds[key] += len(WORKS) * amount
        # Every expanded row has at most H*L entries; comparisons and the
        # FFT input are separately charged from local formula construction.
        formula_bounds["expansion_terms"] += len(WORKS) * H * L
        formula_bounds["comparison_visits"] += len(WORKS) * (3 * Q + 4)
        formula_bounds["scatter_entries"] += len(WORKS) * Q
        # Row norm plus two unnormalized Born-law vectors are all charged.
        formula_bounds["normalization_visits"] += len(WORKS) * 3 * Q
        formula_bounds["formula_fft_calls"] += len(WORKS)
        formula_bounds["formula_fft_entries"] += len(WORKS) * Q
        formula_bounds["fourier_terms"] += len(WORKS) * 4 * ROW_BOUNDS[cut]["components"]
    # One extra cut-11 wrong-left row is retained as a genuine must-fail
    # control.  It is not silently folded into the positive formula ledger.
    for key, amount in ROW_BOUNDS[11].items():
        formula_bounds[f"wrong_{key}"] = amount
    formula_bounds["wrong_expansion_terms"] = H * L
    formula_bounds["wrong_comparison_visits"] = Q

    reference_bounds = {
        "reference_calls": reference_calls,
        "direct_block_products": reference_calls * (PERIOD // BLOCK) * BLOCK * BLOCK,
        "direct_phase_queries": reference_calls * 3 * PERIOD,
        "modular_pow_queries": reference_calls * 3 * PERIOD,
        "direct_shift_ops": reference_calls * 4,
        "column_norm_entries": reference_calls * PERIOD,
        "support_generation_visits": reference_calls,
        "reference_fft_calls": len(WORKS),
        "reference_fft_entries": len(WORKS) * Q,
    }
    setup_bounds = {"matrix_entries": 2 * BLOCK * BLOCK,
                    "unitary_check_terms": 2 * BLOCK ** 3}
    total_reserved = (sum(formula_bounds.values()) + sum(reference_bounds.values())
                      + sum(setup_bounds.values()))
    if any(value > MAX_CATEGORY_TERMS for value in formula_bounds.values()):
        raise MemoryError("formula category exceeds 5M preflight")
    if any(value > MAX_CATEGORY_TERMS for value in reference_bounds.values()):
        raise MemoryError("reference category exceeds 5M preflight")
    if total_reserved > MAX_TOTAL_RESERVED:
        raise MemoryError(f"long-row preflight {total_reserved}>15M")

    # Three selected reference rows and their FFTs are retained, plus formula
    # rows, rebindings and temporary vectors. This is not a Q-by-r archive.
    # Scalar reserves include expansion dictionaries and JSON numeric copies;
    # Python container headers and allocator overhead are not claimed as RSS.
    complex_vectors = 12 * Q * np.dtype(np.complex128).itemsize
    real_vectors = 8 * Q * np.dtype(np.float64).itemsize
    direct_vector = PERIOD * np.dtype(np.complex128).itemsize
    component_scalars = (max(ROW_BOUNDS[c]["components"] for c in CUTS)
                         * 6 * 16)
    call_records = reference_calls * 6 * 16 * 2
    support_lists = reference_calls * 8
    scalar_report = 512 * 1024
    payload = (complex_vectors + real_vectors + direct_vector
               + component_scalars + call_records + support_lists + scalar_report)
    guard_bytes(payload, "long-row aggregate numeric payload")
    return {
        "support_counts": support_counts,
        "reference_bounds": reference_bounds,
        "formula_bounds": dict(formula_bounds),
        "setup_bounds": setup_bounds,
        "total_reserved_terms": total_reserved,
        "planned_numeric_payload_bytes": payload,
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_category_terms": MAX_CATEGORY_TERMS,
        "max_total_reserved_terms": MAX_TOTAL_RESERVED,
    }


def expand_for_report(row, counters):
    """Use the audited cut-row expansion while retaining exact zero omission."""
    return expand_row(row, counters)


def main():
    exp = Experiment("nested_phase_long_rows", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "all long-row amplitudes and complete FFT laws agree")
    exp.predict("P2", "selected Q=8192 geometric laws and norms agree")
    exp.predict("P3", "frozen named reservations and per-cut ledgers are bounded")
    exp.must_fail("C1", "wrong frozen-left coefficient has a named witness")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    formula_actual = defaultdict(int)
    formula_budget = defaultdict(int)
    reference_actual = defaultdict(int)
    reference_budget = defaultdict(int)
    support_actual = defaultdict(int)
    support_budget = defaultdict(int)
    reference_call_records = []
    formula_call_records = []
    setup_actual = defaultdict(int)
    setup_budget = defaultdict(int)
    fft_actual = defaultdict(int)
    fft_budget = defaultdict(int)
    row_details = {}
    try:
        pre = preflight()
        report["preflight"] = pre
        reserve(setup_budget, "matrix_entries", 2 * BLOCK * BLOCK)
        reserve(setup_budget, "unitary_check_terms", 2 * BLOCK ** 3)
        setup_actual["matrix_entries"] = 2 * BLOCK * BLOCK
        setup_actual["unitary_check_terms"] = 2 * BLOCK ** 3
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        for block_matrix in (W0, W1):
            if (not np.all(np.isfinite(block_matrix))
                    or not np.allclose(block_matrix.conj().T @ block_matrix,
                                       np.eye(BLOCK), rtol=0, atol=1e-12)):
                raise AssertionError("work block is not unitary")

        # Sparse literal reference: retain three selected Q-entry rows, while
        # streaming one dense 60-entry literal column per oracle query.
        reference_rows = {}
        max_reference_column_norm_error = 0.
        for work in WORKS:
            ref_row = np.zeros(Q, dtype=np.complex128)
            support_n = pre["support_counts"][str(work)]
            reserve(support_budget, "support_generation_visits", support_n)
            for exponent in support_exponents(work, support_actual):
                reserve(reference_budget, "reference_calls")
                for key, amount in {
                    "direct_block_products": (PERIOD // BLOCK) * BLOCK * BLOCK,
                    "direct_phase_queries": 3 * PERIOD,
                    "modular_pow_queries": 3 * PERIOD,
                    "direct_shift_ops": 4,
                    "column_norm_entries": PERIOD}.items():
                    reserve(reference_budget, key, amount)
                call_counter = defaultdict(int)
                literal = None
                try:
                    literal = direct_column(
                        exponent, 5, W0, W1, call_counter,
                        split=SPLIT, early_insertions=POSITIONS)
                finally:
                    delta = dict(call_counter)
                    for key, value in delta.items():
                        reference_actual[key] += int(value)
                    reference_call_records.append({
                        "work": work, "exponent": exponent,
                        "counters": delta, "completed": literal is not None})
                if literal is None:
                    raise AssertionError("literal reference returned no vector")
                reference_actual["column_norm_entries"] += PERIOD
                column_norm = float(np.sum(np.abs(literal) ** 2))
                max_reference_column_norm_error = max(
                    max_reference_column_norm_error, abs(column_norm - 1.))
                if abs(column_norm - 1.) > TOL or not np.all(np.isfinite(literal)):
                    raise AssertionError("literal reference column is not normalized")
                ref_row[exponent] = literal[work]
            reference_rows[work] = ref_row

        reference_fft = {}
        for work in WORKS:
            reserve(fft_budget, "reference_fft_calls")
            reserve(fft_budget, "reference_fft_entries", Q)
            fft_actual["reference_fft_calls"] += 1
            fft_actual["reference_fft_entries"] += Q
            reference_fft[work] = np.fft.fft(reference_rows[work])

        early_phases = ((5, g_phase), (11, g_phase))
        max_amp_error = max_norm_error = 0.
        max_law_error = max_probability_error = 0.
        max_full_law_tv = 0.
        max_selected_probability_error = 0.
        max_geometric_error = 0.
        max_components = 0
        for cut in CUTS:
            for work in WORKS:
                bound = ROW_BOUNDS[cut]
                row_key = f"cut{cut}:work{work}"
                for key, amount in bound.items():
                    if key != "components":
                        reserve(formula_budget, key, amount)
                reserve(formula_budget, "expansion_terms", H * L)
                reserve(formula_budget, "comparison_visits", 3 * Q + 4)
                reserve(formula_budget, "scatter_entries", Q)
                reserve(formula_budget, "normalization_visits", 3 * Q)
                reserve(formula_budget, "fourier_terms",
                        4 * bound["components"])
                before = dict(formula_actual)
                row = None
                try:
                    row = cut_row(
                        work, period=PERIOD, block=BLOCK, width=WIDTH,
                        split=SPLIT, initial=W0, middle=W1, phase=g_phase,
                        early_phases=early_phases, cut=cut,
                        counters=formula_actual)
                finally:
                    formula_call_records.append({
                        "cut": cut, "work": work,
                        "counters": {key: formula_actual.get(key, 0)
                                     - before.get(key, 0)
                                     for key in formula_actual},
                        "completed": row is not None})
                if row is None:
                    raise AssertionError("cut row returned no result")
                expanded = expand_for_report(row, formula_actual)
                if len(row["components"]) > bound["components"]:
                    raise AssertionError(f"component bound exceeded {row_key}")
                formula_vec = np.zeros(Q, dtype=np.complex128)
                for exponent, value in expanded.items():
                    charge(formula_actual, "scatter_entries")
                    if not 0 <= int(exponent) < Q:
                        raise AssertionError("expanded exponent outside Q")
                    formula_vec[int(exponent)] = value
                ref_row = reference_rows[work]
                if (not math.isfinite(row["norm"]) or row["norm"] <= 0
                        or not np.all(np.isfinite(formula_vec))):
                    raise AssertionError("nonfinite/zero formula row")
                formula_actual["normalization_visits"] += Q
                norm = float(np.sum(np.abs(ref_row) ** 2))
                if norm <= 0 or not math.isfinite(norm):
                    raise AssertionError("zero/nonfinite retained work row")
                norm_error = abs(float(row["norm"]) - norm)
                max_norm_error = max(max_norm_error, norm_error)
                for exponent in range(Q):
                    formula_actual["comparison_visits"] += 1
                    error = abs(formula_vec[exponent] - ref_row[exponent])
                    max_amp_error = max(max_amp_error, error)
                    if error > TOL:
                        raise AssertionError(
                            f"long-row amplitude mismatch cut={cut} work={work} "
                            f"exponent={exponent}")
                reserve(fft_budget, "formula_fft_calls")
                reserve(fft_budget, "formula_fft_entries", Q)
                fft_actual["formula_fft_calls"] += 1
                fft_actual["formula_fft_entries"] += Q
                formula_fft = np.fft.fft(formula_vec)
                formula_actual["normalization_visits"] += 2 * Q
                ref_sq = np.abs(reference_fft[work]) ** 2
                formula_sq = np.abs(formula_fft) ** 2
                ref_total = float(ref_sq.sum())
                formula_total = float(formula_sq.sum())
                if (not math.isfinite(ref_total) or not math.isfinite(formula_total)
                        or ref_total <= 0 or formula_total <= 0
                        or abs(ref_total / (Q * norm) - 1.) > TOL
                        or abs(formula_total / (Q * row["norm"]) - 1.) > TOL):
                    raise AssertionError(f"unnormalized FFT law total mismatch {row_key}")
                # Normalize each law by its own positive total before TV.
                ref_law = ref_sq / ref_total
                formula_law = formula_sq / formula_total
                charge(formula_actual, "comparison_visits", 2 * Q)
                law_error = float(np.max(np.abs(formula_law - ref_law)))
                max_law_error = max(max_law_error, law_error)
                law_tv = float(np.sum(np.abs(formula_law - ref_law)) / 2.)
                max_probability_error = max(max_probability_error, law_tv)
                max_full_law_tv = max(max_full_law_tv, law_tv)
                if law_error > TOL or law_tv > TOL:
                    raise AssertionError(f"long-row FFT law mismatch {row_key}")
                for output in (0, 1, Q // 2, Q - 1):
                    charge(formula_actual, "comparison_visits")
                    geometric = selected_geometric_long(row, output, formula_actual)
                    ferr = abs(geometric - reference_fft[work][output])
                    max_geometric_error = max(max_geometric_error, ferr)
                    selected_probability_error = abs(
                        abs(geometric) ** 2 / formula_total
                        - ref_sq[output] / ref_total)
                    max_selected_probability_error = max(
                        max_selected_probability_error, selected_probability_error)
                    max_probability_error = max(max_probability_error,
                                                selected_probability_error)
                    if ferr > TOL:
                        raise AssertionError(
                            f"geometric query mismatch {row_key} y={output}")
                max_components = max(max_components, len(row["components"]))
                row_delta = {
                    key: formula_actual.get(key, 0) - before.get(key, 0)
                    for key in formula_actual}
                # Independently enumerated integer geometry was supplied
                # before the amplitude run. Nonzero components can be fewer.
                predicted_coeff = 288 if cut == 11 else (1224 if work == 59 else 1230)
                predicted_candidates = 160 if cut == 11 else (680 if work == 59 else 683)
                if (row_delta["coefficient_pair_visits"] != predicted_coeff
                        or row_delta["phase_products"] != 3 * predicted_coeff + predicted_candidates
                        or row_delta["group_pair_visits"] != bound["group_pair_visits"]):
                    raise AssertionError("actual counters disagree with frozen integer geometry")
                row_details[row_key] = {
                    "cut": cut, "work": work,
                    "component_count": len(row["components"]),
                    "expanded_count": len(expanded),
                    "stride": row["stride"],
                    "norm": row["norm"],
                    "K": row["K"], "P_left": row["P_left"],
                    "structural_factor": row["structural_factor"],
                    "geometric_candidates_before_cancellation": predicted_candidates,
                    "progression_count_min": min(
                        (int(count) for _start, count, _gamma in row["components"]),
                        default=0),
                    "progression_count_max": max(
                        (int(count) for _start, count, _gamma in row["components"]),
                        default=0),
                    "actual_counter_delta": row_delta,
                    "cost_ratios_to_frozen_bounds": {
                        key: (row_delta.get(key, 0) / bound[key])
                        for key in ("group_pair_visits", "coefficient_pair_visits",
                                    "phase_calls", "phase_products")},
                }

        comparisons = {}
        for work in WORKS:
            hybrid = row_details[f"cut11:work{work}"]
            hybrid_counts = hybrid["actual_counter_delta"]
            comparisons[str(work)] = {}
            for baseline_cut in (5, 12):
                baseline = row_details[f"cut{baseline_cut}:work{work}"]
                counts = baseline["actual_counter_delta"]
                comparisons[str(work)][str(baseline_cut)] = {
                    "baseline_to_hybrid_local_pair_ratio":
                        (counts["group_pair_visits"] + counts["coefficient_pair_visits"])
                        / (hybrid_counts["group_pair_visits"] + hybrid_counts["coefficient_pair_visits"]),
                    "baseline_to_hybrid_phase_call_ratio": counts["phase_calls"] / hybrid_counts["phase_calls"],
                    "baseline_to_hybrid_component_ratio": baseline["component_count"] / hybrid["component_count"],
                }

        # Wrong-left-z=0 control on the frozen named work-0 row only.
        wrong_witness = None
        wrong_scanned = []
        for work in (WORKS[0],):
            if wrong_witness is not None:
                break
            bound = ROW_BOUNDS[11]
            for key, amount in bound.items():
                if key != "components":
                    reserve(formula_budget, f"wrong_{key}", amount)
            reserve(formula_budget, "wrong_expansion_terms", H * L)
            reserve(formula_budget, "wrong_comparison_visits", Q)
            before = dict(formula_actual)
            wrong = None
            try:
                wrong = cut_row(
                    work, period=PERIOD, block=BLOCK, width=WIDTH,
                    split=SPLIT, initial=W0, middle=W1, phase=g_phase,
                    early_phases=early_phases, cut=11,
                    counters=formula_actual, wrong_left_z0=True)
            finally:
                formula_call_records.append({
                    "cut": 11, "work": work, "wrong_left_z0": True,
                    "counters": {key: formula_actual.get(key, 0)
                                 - before.get(key, 0)
                                 for key in formula_actual},
                    "completed": wrong is not None})
            if wrong is None:
                raise AssertionError("wrong-left row returned no result")
            expanded = expand_for_report(wrong, formula_actual)
            wrong_scanned.append(work)
            ref_row = reference_rows[work]
            for exponent in range(Q):
                formula_actual["wrong_comparison_visits"] += 1
                error = abs(expanded.get(exponent, 0j) - ref_row[exponent])
                if error > TOL:
                    wrong_witness = {"work": work, "exponent": exponent,
                                     "error": error}
                    break
        c1 = wrong_witness is not None
        p1 = bool(max_amp_error <= TOL and max_norm_error <= TOL)
        p2 = bool(max_law_error <= TOL and max_geometric_error <= TOL
                  and max_probability_error <= TOL
                  and max_full_law_tv <= TOL
                  and max_selected_probability_error <= TOL
                  and max_reference_column_norm_error <= TOL)
        combined_formula_bounds = {
            "group_pair_visits": formula_budget.get("group_pair_visits", 0)
            + formula_budget.get("wrong_group_pair_visits", 0),
            "coefficient_pair_visits": formula_budget.get("coefficient_pair_visits", 0)
            + formula_budget.get("wrong_coefficient_pair_visits", 0),
            "phase_calls": formula_budget.get("phase_calls", 0)
            + formula_budget.get("wrong_phase_calls", 0),
            "phase_products": formula_budget.get("phase_products", 0)
            + formula_budget.get("wrong_phase_products", 0),
            "expansion_terms": formula_budget.get("expansion_terms", 0)
            + formula_budget.get("wrong_expansion_terms", 0),
            "comparison_visits": formula_budget.get("comparison_visits", 0)
            + formula_budget.get("wrong_comparison_visits", 0),
            "wrong_comparison_visits": formula_budget.get("wrong_comparison_visits", 0),
            "fourier_terms": formula_budget.get("fourier_terms", 0),
            "normalization_visits": formula_budget.get("normalization_visits", 0),
            "scatter_entries": formula_budget.get("scatter_entries", 0),
        }
        formula_reconciled = all(
            formula_actual.get(key, 0) <= value
            for key, value in combined_formula_bounds.items())
        fft_reconciled = all(
            fft_actual.get(key, 0) <= fft_budget.get(key, 0)
            for key in fft_actual)
        reference_reconciled = all(
            reference_actual.get(key, 0) <= reference_budget.get(key, 0)
            for key in reference_actual)
        support_reconciled = all(
            support_actual.get(key, 0) <= support_budget.get(key, 0)
            for key in support_actual)
        setup_reconciled = all(
            setup_actual.get(key, 0) <= setup_budget.get(key, 0)
            for key in setup_actual)
        reservation_reconciled = (
            all(value <= pre["formula_bounds"].get(key, -1)
                for key, value in formula_budget.items())
            and all(value <= pre["reference_bounds"].get(key, -1)
                    for key, value in reference_budget.items())
            and all(value <= pre["reference_bounds"].get(key, -1)
                    for key, value in support_budget.items())
            and all(value <= (pre["formula_bounds"].get(
                key, pre["reference_bounds"].get(key, -1)))
                    for key, value in fft_budget.items())
            and all(value <= pre["setup_bounds"].get(key, -1)
                    for key, value in setup_budget.items()))
        actual_counter_checksum = (sum(formula_actual.values())
                                   + sum(reference_actual.values())
                                   + sum(support_actual.values())
                                   + sum(fft_actual.values())
                                   + sum(setup_actual.values()))
        reservation_checksum = (sum(formula_budget.values())
                                + sum(reference_budget.values())
                                + sum(support_budget.values())
                                + sum(fft_budget.values())
                                + sum(setup_budget.values()))
        combined_formula = sum(formula_actual.values())
        combined_reference = sum(reference_actual.values())
        all_bounds = list(pre["formula_bounds"].values()) + list(
            pre["reference_bounds"].values()) + list(pre["setup_bounds"].values())
        p3 = bool(
            all(value <= MAX_CATEGORY_TERMS for value in formula_actual.values())
            and all(value <= MAX_CATEGORY_TERMS for value in reference_actual.values())
            and combined_formula + combined_reference <= MAX_TOTAL_RESERVED
            and all(value <= MAX_CATEGORY_TERMS for value in formula_budget.values())
            and all(value <= MAX_CATEGORY_TERMS for value in reference_budget.values())
            and all(value <= MAX_CATEGORY_TERMS for value in support_budget.values())
            and all(value <= MAX_CATEGORY_TERMS for value in fft_budget.values())
            and all(value <= MAX_CATEGORY_TERMS for value in setup_budget.values())
            and all(value <= MAX_CATEGORY_TERMS for value in all_bounds)
            and formula_reconciled and fft_reconciled
            and reference_reconciled and support_reconciled
            and setup_reconciled and reservation_reconciled
            and actual_counter_checksum <= pre["total_reserved_terms"]
            and reservation_checksum <= pre["total_reserved_terms"])
        report.update({
            "fixture": {"N": N, "base": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L,
                        "phase_positions": list(POSITIONS),
                        "cuts": list(CUTS), "works": list(WORKS),
                        "support_deltas": list(DELTA_SUPPORT)},
            "row_details": row_details,
            "same_row_construction_comparisons": comparisons,
            "max_components": max_components,
            "max_amplitude_error": max_amp_error,
            "max_norm_error": max_norm_error,
            "max_full_law_error": max_law_error,
            "max_full_law_tv": max_full_law_tv,
            "max_selected_probability_error": max_selected_probability_error,
            "max_geometric_error": max_geometric_error,
            "max_reference_column_norm_error": max_reference_column_norm_error,
            "wrong_left_scanned_works": wrong_scanned,
            "wrong_left_witness": wrong_witness,
            "formula_actual": dict(formula_actual),
            "formula_budget": dict(formula_budget),
            "combined_formula_bounds": combined_formula_bounds,
            "formula_reconciled": formula_reconciled,
            "reference_actual": dict(reference_actual),
            "reference_budget": dict(reference_budget),
            "support_actual": dict(support_actual),
            "support_budget": dict(support_budget),
            "reference_call_records": reference_call_records,
            "formula_call_records": formula_call_records,
            "fft_actual": dict(fft_actual),
            "fft_budget": dict(fft_budget),
            "fft_reconciled": fft_reconciled,
            "reference_reconciled": reference_reconciled,
            "support_reconciled": support_reconciled,
            "setup_reconciled": setup_reconciled,
            "reservation_reconciled": reservation_reconciled,
            "setup_actual": dict(setup_actual),
            "setup_budget": dict(setup_budget),
            "overlapping_counter_checksum_not_physical_cost": actual_counter_checksum,
            "reservation_checksum_overlapping": reservation_checksum,
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
            "status": "PASS" if p1 and p2 and p3 and c1 else "FAIL",
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "scope": "bounded long-row coherent amplitudes; no sampler/timing claim",
        })
        exp.check("P1", p1, "all expanded amplitudes and complete FFT laws")
        exp.check("P2", p2, "selected geometric laws and normalization")
        exp.check("P3", p3, "frozen reservations and ledgers")
        exp.fail_check("C1", c1, "wrong frozen-left coefficient has witness")
    except Exception as exc:
        report.update({
            "exception": repr(exc),
            "traceback": __import__("traceback").format_exc(),
            "formula_actual": dict(formula_actual),
            "formula_budget": dict(formula_budget),
            "reference_actual": dict(reference_actual),
            "reference_budget": dict(reference_budget),
            "support_actual": dict(support_actual),
            "support_budget": dict(support_budget),
            "reference_call_records": reference_call_records,
            "formula_call_records": formula_call_records,
            "fft_actual": dict(fft_actual),
            "fft_budget": dict(fft_budget),
            "setup_actual": dict(setup_actual),
            "setup_budget": dict(setup_budget),
        })
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before long-row audit")
        exp.fail_check("C1", False, "exception before wrong-left control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "fixture": "r60,b3,t13,s12, positions(5,11), cuts(5,11,12)",
        "scope": "bounded long-row coherent amplitude pilot",
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_category_terms": MAX_CATEGORY_TERMS,
        "max_total_reserved_terms": MAX_TOTAL_RESERVED,
        "execution": "released bounded initial run",
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
