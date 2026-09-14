"""Coherent-amplitude pilot for the C81 nested binary phase cut.

The frozen family is N=61, a=2, r=60, b=3, t=8, s=7 with W0/W1 and the
late G1 phase from C79.  Two identical G1 phases occur at (3,v), v=3..7.
This experiment-local constructor implements the independently derived cut
row; the literal reference is the existing ER ``direct_column`` with its
multi-insertion option.  No sampler or generic propagator is introduced.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every cut row's expanded amplitudes and positive disjoint norm agree
      coordinatewise with the literal reference for all works and frozen cuts.
  P2  Selected geometric Fourier amplitudes and normalized conditional
      probabilities agree with independent FFT values; duplicate phases at
      v=3 and endpoint phases at v=7 remain covered.
  P3  The independently computed cut factors and row component construction
      remain within the frozen per-category budgets; actual construction
      savings are recorded as an empirical branch, not predicted.
  C1  Freezing the left phase coefficient at z=0 while retaining the correct
      refined supports fails against the literal reference at a named row and
      exponent witness.

This is a bounded complex128 coefficient/reference audit, not a numerical
certificate, sampler integration, timing claim, or novelty result.
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
from lab.work_first import EarlierPhaseProgressions
from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_earlier_phase_cycles import direct_column


N, BASE, PERIOD, BLOCK, WIDTH, SPLIT = 61, 2, 60, 3, 8, 7
Q, L, H = 1 << WIDTH, 1 << SPLIT, 1 << (WIDTH - SPLIT)
POSITIONS = tuple(range(3, 8))
MAX_BYTES = 16 << 20
MAX_CATEGORY_TERMS = 5_000_000
MAX_TOTAL_RESERVED = 15_000_000
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_amplitudes_{stamp}.json"
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
    if int(value) > MAX_BYTES:
        raise MemoryError(f"{label} {value}>{MAX_BYTES}")


def g_phase(index):
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def phase_factor(index, oracle, period):
    """Query an indexed oracle after canonicalizing to the supplied period."""
    period = int(period)
    if period <= 0:
        raise ValueError("period must be positive")
    value = complex(oracle(int(index) % period))
    if (not math.isfinite(value.real) or not math.isfinite(value.imag)
            or abs(abs(value) - 1.) > 1e-12):
        raise ValueError("phase oracle returned a non-unit value")
    return value


def charge(counters, key, amount=1, max_terms=MAX_CATEGORY_TERMS):
    amount = int(amount)
    if counters.get(key, 0) + amount > max_terms:
        raise MemoryError(f"{key} cap before operation")
    counters[key] = counters.get(key, 0) + amount


def cut_row(work, *, period, block, width, split, initial, middle, phase,
            early_phases, cut, counters, max_terms=MAX_CATEGORY_TERMS,
            wrong_left_z0=False):
    """Construct one cut row as (components,stride,norm,K,P_left,factor).

    Components are ``(start,count,gamma)`` with common ``stride``.  All
    candidate rho classes are formed before zero coefficients are omitted;
    only exact zero gamma values are dropped from the returned support.
    """
    period, block, width, split, cut = map(int, (period, block, width, split, cut))
    length = 1 << split
    A = 1 << cut
    K = length // A
    left_positions = [int(v) for v, _oracle in early_phases if int(v) < cut]
    if left_positions:
        V = max(left_positions)
        P_left = (1 << V) // math.gcd(period, 1 << V)
    else:
        P_left = 1
    components = []
    for high in range(1 << (width - split)):
        v_work = (int(work) - length * high) % period
        cell, p = divmod(v_work, block)
        charge(counters, "phase_calls", 1, max_terms)
        late = phase_factor(v_work, phase, period)
        for k in range(K):
            pairs = defaultdict(list)
            for u in range(block):
                for q in range(block):
                    charge(counters, "group_pair_visits", 1, max_terms)
                    rho = (cell * block + q - A * k - u) % period
                    x = cell * block + q
                    pairs[rho].append((u, q, x,
                                       complex(middle[p, q] * initial[u, 0])))
            for rho, local_pairs in sorted(pairs.items()):
                count_n = max(0, 1 + (A - 1 - rho) // period)
                for z in range(min(P_left, count_n)):
                    a = rho + period * z
                    gamma = 0j
                    for u, _q, x, amplitude in local_pairs:
                        charge(counters, "coefficient_pair_visits", 1, max_terms)
                        product = 1.+0j
                        for v, oracle in early_phases:
                            v = int(v)
                            if v < cut:
                                argument = u + (a % (1 << v))
                            else:
                                argument = x - (1 << v) * (k // (1 << (v - cut)))
                            if wrong_left_z0 and v < cut:
                                # Deliberate control: preserve supports but
                                # freeze only the left coefficient at z=0.
                                argument = u + (rho % (1 << v))
                            charge(counters, "phase_calls", 1, max_terms)
                            value = phase_factor(argument, oracle, period)
                            charge(counters, "phase_products", 1, max_terms)
                            product *= value
                        charge(counters, "phase_products", 1, max_terms)
                        gamma += amplitude * product
                    charge(counters, "phase_products", 1, max_terms)
                    gamma *= late
                    if gamma == 0j:
                        continue
                    count = 1 + (count_n - 1 - z) // P_left
                    start = length * high + A * k + rho + period * z
                    stride = period * P_left
                    components.append((start, count, gamma))
    norm = math.fsum(count * abs(gamma)**2
                     for _start, count, gamma in components)
    return {"components": tuple(components), "stride": period * P_left,
            "norm": norm, "K": K, "P_left": P_left,
            "structural_factor": K * P_left,
            "wrong_left_z0": bool(wrong_left_z0)}


def expand_row(row, counters):
    expanded = {}
    for start, count, gamma in row["components"]:
        for index in range(int(count)):
            charge(counters, "expansion_terms")
            exponent = int(start) + int(row["stride"]) * index
            if exponent in expanded:
                raise AssertionError("overlapping cut components")
            expanded[exponent] = complex(gamma)
    return expanded


def selected_geometric(row, output, counters):
    total = 0j
    for start, count, gamma in row["components"]:
        charge(counters, "fourier_terms")
        total += (gamma * unit_phase(-int(output) * int(start), Q)
                  * geometric_sum(int(count), -int(output) * int(row["stride"]), Q))
    return total


def row_work_bounds(v, cut):
    """Conservative per-work bounds valid for arbitrary complex b-blocks.

    R bounds the number of possible residue classes geometrically; F bounds
    the number of refined z values in one class.  The two component-count
    bounds are the independent ``K`` and full-row bounds from the cut proof.
    """
    A = 1 << int(cut)
    K = L // A
    left = [3, int(v)] if int(v) != 3 else [3, 3]
    left = [x for x in left if x < int(cut)]
    P = ((1 << max(left)) // math.gcd(PERIOD, 1 << max(left))) if left else 1
    R = min(PERIOD, 2 * BLOCK - 1)
    F = min(P, (A + PERIOD - 1) // PERIOD)
    m_bound = min(H * K * min(A, R * F),
                  H * R * ((L + PERIOD - 1) // PERIOD))
    coefficient_bound = H * K * BLOCK * BLOCK * F
    return {
        "group_pair_visits": H * K * BLOCK * BLOCK,
        "coefficient_pair_visits": coefficient_bound,
        "phase_calls": coefficient_bound * 2 + H,
        "phase_products": coefficient_bound * 3 + m_bound,
        "expansion_terms": H * L,
        "fourier_terms": 4 * m_bound,
        "components": m_bound,
        "K": K, "P_left": P, "F": F, "R": R,
    }


def preflight():
    cuts_by_position = {v: tuple(sorted({3, v, 7})) for v in POSITIONS}
    bounds = defaultdict(int)
    for v in POSITIONS:
        for cut in cuts_by_position[v]:
            one = row_work_bounds(v, cut)
            for key in ("group_pair_visits", "coefficient_pair_visits",
                        "phase_calls", "phase_products", "expansion_terms",
                        "fourier_terms"):
                bounds[key] += PERIOD * one[key]
    # The wrong-left control repeats one complete (3,6), cut-6 row family.
    wrong = row_work_bounds(6, 6)
    for key in ("group_pair_visits", "coefficient_pair_visits", "phase_calls",
                "phase_products", "expansion_terms"):
        bounds[f"wrong_{key}"] = wrong[key] * PERIOD
    bounds["reference_calls"] = len(POSITIONS) * Q
    # These are the actual counter names used by direct_column.  The loose
    # upper bounds charge every physical label, even when a vector entry is
    # zero, before the reference call is made.
    bounds["direct_block_products"] = len(POSITIONS) * Q * (PERIOD // BLOCK) * BLOCK * BLOCK
    bounds["direct_phase_queries"] = len(POSITIONS) * Q * 3 * PERIOD
    bounds["direct_shift_ops"] = len(POSITIONS) * Q * 5
    bounds["modular_pow_queries"] = bounds["direct_phase_queries"]
    bounds["unitary_check_terms"] = 4 * BLOCK ** 3
    bounds["matrix_entries"] = 2 * BLOCK * BLOCK
    # One FFT per (phase position, work), with a conservative entry charge.
    bounds["fft_calls"] = len(POSITIONS) * PERIOD
    bounds["fft_entries"] = bounds["fft_calls"] * Q
    # Duplicate-v=3 EarlierPhaseProgressions baseline, all r rows.
    duplicate_P = (1 << 3) // math.gcd(PERIOD, 1 << 3)
    duplicate_classes = min(duplicate_P, (L + PERIOD - 1) // PERIOD)
    bounds["duplicate_baseline_row_terms"] = PERIOD * H * BLOCK * BLOCK * (1 + duplicate_classes)
    bounds["duplicate_baseline_phase_queries"] = PERIOD * (
        BLOCK + BLOCK * BLOCK + H + H * BLOCK * BLOCK * duplicate_classes)
    bounds["duplicate_baseline_expansion_terms"] = PERIOD * H * L
    bounds["duplicate_baseline_phase_squares"] = bounds["duplicate_baseline_phase_queries"]
    row_count = PERIOD * sum(len(cuts) for cuts in cuts_by_position.values())
    bounds["comparison_entries"] = (row_count + 2 * PERIOD) * Q
    bounds["norm_entries"] = (2 * len(POSITIONS) * Q * PERIOD
                              + (row_count + PERIOD) * Q)
    total = sum(bounds.values())
    if any(value > MAX_CATEGORY_TERMS for value in bounds.values()):
        raise MemoryError("named category exceeds 5M preflight")
    if total > MAX_TOTAL_RESERVED:
        raise MemoryError(f"total preflight {total}>{MAX_TOTAL_RESERVED}")
    ref = Q * PERIOD * np.dtype(np.complex128).itemsize
    payload = 4 * ref + 8 * PERIOD * np.dtype(np.complex128).itemsize + 2 * ref + 256 * 1024
    guard_bytes(payload, "nested amplitude aggregate payload")
    return {"cuts_by_position": {str(k): list(v) for k, v in cuts_by_position.items()},
            "category_bounds": dict(bounds), "total_reserved_terms": total,
            "planned_numeric_payload_bytes": payload,
            "max_numeric_payload_bytes": MAX_BYTES,
            "max_category_terms": MAX_CATEGORY_TERMS,
            "max_total_reserved_terms": MAX_TOTAL_RESERVED}


def main():
    exp = Experiment("nested_phase_amplitudes", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "cut amplitudes and norms match literal columns")
    exp.predict("P2", "selected geometric amplitudes and normalized laws match FFT")
    exp.predict("P3", "all named cut/reference workloads remain bounded")
    exp.must_fail("C1", "freezing left z coefficient fails on an interior row")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    counters = defaultdict(int)
    reference_counters = defaultdict(int)
    reference_budget = defaultdict(int)
    row_budget = defaultdict(int)
    baseline_budget = defaultdict(int)
    fft_budget = defaultdict(int)
    setup_counters = defaultdict(int)
    baseline_counters = defaultdict(int)
    baseline_live = []
    try:
        pre = preflight()
        report["preflight"] = pre
        # Reserve and perform the two small-block unitarity checks before any
        # reference or formula products.  This is setup accounting, not a
        # claim about native BLAS work.
        charge(setup_counters, "matrix_entries", 2 * BLOCK * BLOCK)
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        charge(setup_counters, "unitary_check_terms", 2 * BLOCK ** 3)
        for block_matrix in (W0, W1):
            if (not np.all(np.isfinite(block_matrix))
                    or not np.allclose(block_matrix.conj().T @ block_matrix,
                                       np.eye(BLOCK), rtol=0, atol=1e-12)):
                raise AssertionError("frozen work block is not unitary")
        early_cache = {v: ((3, g_phase), (v, g_phase)) for v in POSITIONS}
        summaries = {}
        max_amp_error = max_norm_error = max_fourier_error = 0.
        max_probability_error = 0.
        max_literal_column_norm_error = 0.
        work_masses = {}
        # Keep only the v=6 reference needed by the explicit wrong-left
        # control; streaming the other four matrices avoids retaining a
        # needless five-matrix archive.
        direct_for_wrong = None
        direct_for_duplicate = None
        fft_cache = {}
        cut_counter_totals = defaultdict(lambda: defaultdict(int))
        work0_counter_examples = {}
        for v in POSITIONS:
            positions = (3, v)
            direct = np.zeros((Q, PERIOD), dtype=np.complex128)
            for exponent in range(Q):
                # Reserve the conservative full-label work before entering
                # direct_column; the returned counters below record actual
                # operations and are not precharged.
                charge(reference_budget, "reference_calls")
                for key, amount in {
                    "direct_block_products": (PERIOD // BLOCK) * BLOCK * BLOCK,
                    "direct_phase_queries": 3 * PERIOD,
                    "direct_shift_ops": 5,
                    "modular_pow_queries": 3 * PERIOD}.items():
                    charge(reference_budget, key, amount)
                literal = direct_column(
                    exponent, 3, W0, W1, reference_counters,
                    early_insertions=positions)
                direct[exponent] = literal
                if not np.all(np.isfinite(literal)):
                    raise AssertionError("nonfinite literal reference column")
                charge(counters, "norm_entries", PERIOD)
                column_norm = float(np.sum(np.abs(literal) ** 2))
                max_literal_column_norm_error = max(
                    max_literal_column_norm_error, abs(column_norm - 1.))
                if abs(column_norm - 1.) > TOL:
                    raise AssertionError(f"literal column norm mismatch: {column_norm}")
            charge(counters, "norm_entries", Q * PERIOD)
            work_masses[v] = float(np.sum(np.abs(direct) ** 2)) / Q
            if abs(work_masses[v] - 1.) > TOL:
                raise AssertionError("work marginal does not normalize")
            if v == 6:
                direct_for_wrong = direct.copy()
            if v == 3:
                direct_for_duplicate = direct.copy()
            for work in range(PERIOD):
                charge(fft_budget, "fft_calls")
                charge(fft_budget, "fft_entries", Q)
                fft_cache[work] = np.fft.fft(direct[:, work])
                counters["fft_calls"] += 1
                counters["fft_entries"] += Q
            for cut in pre["cuts_by_position"][str(v)]:
                cut = int(cut)
                one_bound = row_work_bounds(v, cut)
                for work in range(PERIOD):
                    # Reserve this row's complete formula/expansion/Fourier
                    # workload before invoking any local operation.
                    for key in ("group_pair_visits", "coefficient_pair_visits",
                                "phase_calls", "phase_products",
                                "expansion_terms", "fourier_terms"):
                        charge(row_budget, key, one_bound[key])
                    before_row = dict(counters)
                    row = cut_row(
                        work, period=PERIOD, block=BLOCK, width=WIDTH,
                        split=SPLIT, initial=W0, middle=W1, phase=g_phase,
                        early_phases=early_cache[v], cut=cut,
                        counters=counters)
                    if len(row["components"]) > one_bound["components"]:
                        raise AssertionError("component bound exceeded")
                    expanded = expand_row(row, counters)
                    # expand_row already rejects overlapping components; a
                    # component can of course contain several stride entries.
                    if not expanded and row["norm"] > 0:
                        raise AssertionError("positive row norm has empty expansion")
                    charge(counters, "norm_entries", Q)
                    expected_norm = float(np.sum(np.abs(direct[:, work])**2))
                    norm_error = abs(float(row["norm"]) - expected_norm)
                    max_norm_error = max(max_norm_error, norm_error)
                    charge(counters, "comparison_entries", Q)
                    for exponent in range(Q):
                        observed = expanded.get(exponent, 0j)
                        error = abs(observed - direct[exponent, work])
                        max_amp_error = max(max_amp_error, error)
                        if error > TOL:
                            raise AssertionError(
                                f"amplitude mismatch v={v} cut={cut} "
                                f"work={work} exponent={exponent}")
                    fft = fft_cache[work]
                    for output in (0, 1, Q // 2, Q - 1):
                        geometric = selected_geometric(row, output, counters)
                        ferr = abs(geometric - fft[output])
                        max_fourier_error = max(max_fourier_error, ferr)
                        probability_error = abs(
                            (abs(geometric)**2 / (Q * row["norm"]))
                            - (abs(fft[output])**2 / (Q * expected_norm)))
                        max_probability_error = max(max_probability_error,
                                                    probability_error)
                        if ferr > TOL or probability_error > TOL:
                            raise AssertionError(
                                f"selected Fourier mismatch v={v} cut={cut} "
                                f"work={work} output={output}")
                    row_delta = {key: counters.get(key, 0) - before_row.get(key, 0)
                                 for key in ("group_pair_visits",
                                             "coefficient_pair_visits",
                                             "phase_calls", "phase_products",
                                             "expansion_terms", "fourier_terms")}
                    cut_key = f"v{v}:cut{cut}"
                    for key, amount in row_delta.items():
                        cut_counter_totals[cut_key][key] += amount
                    if work == 0:
                        work0_counter_examples[cut_key] = row_delta
                    key = f"v{v}:cut{cut}"
                    entry = summaries.setdefault(key, {
                        "positions": positions, "cut": cut,
                        "factor": row["structural_factor"],
                        "P_left": row["P_left"], "K": row["K"],
                        "max_components": 0, "max_norm_error": 0.,
                    })
                    entry["max_components"] = max(
                        entry["max_components"], len(row["components"]))
                    entry["max_norm_error"] = max(
                        entry["max_norm_error"], norm_error)
        # Genuine wrong-left control on the frozen (3,6), cut-6 family.
        wrong_witness = None
        wrong_bound = row_work_bounds(6, 6)
        for work in range(PERIOD):
            for key in ("group_pair_visits", "coefficient_pair_visits",
                        "phase_calls", "phase_products", "expansion_terms"):
                charge(row_budget, key, wrong_bound[key])
            row = cut_row(
                work, period=PERIOD, block=BLOCK, width=WIDTH, split=SPLIT,
                initial=W0, middle=W1, phase=g_phase,
                early_phases=early_cache[6], cut=6, counters=counters,
                wrong_left_z0=True)
            expanded = expand_row(row, counters)
            direct = direct_for_wrong
            if direct is None:
                raise AssertionError("missing streamed v=6 reference")
            charge(counters, "comparison_entries", Q)
            for exponent in range(Q):
                error = abs(expanded.get(exponent, 0j) - direct[exponent, work])
                if error > TOL and wrong_witness is None:
                    wrong_witness = {"work": work, "exponent": exponent,
                                     "error": error}
        # Strong C79 positive baseline: two identical phases at the same
        # insertion are represented by one EarlierPhaseProgressions helper
        # with early oracle g^2.  Compare every expanded row to the streamed
        # literal v=3 reference, rather than only checking a marginal.
        if direct_for_duplicate is None:
            raise AssertionError("missing streamed duplicate-v=3 reference")
        def phase_squared(index):
            value = g_phase(int(index))
            charge(baseline_counters, "phase_squares")
            return value * value
        class AuditedEarlier(EarlierPhaseProgressions):
            def _new_counters(self):
                live = super()._new_counters()
                baseline_live.append(live)
                return live
        charge(setup_counters, "unitary_check_terms", 2 * BLOCK ** 3)
        duplicate = AuditedEarlier(
            PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, g_phase,
            early_split=3, early_phase=phase_squared, cover="auto",
            max_local_terms=1_000_000, max_components=4096,
            max_payload_bytes=MAX_BYTES)
        duplicate_stats = duplicate.stats()
        duplicate_row_pair_bound = int(duplicate_stats[
            f"{duplicate_stats['cover']}_row_pair_visit_bound"])
        duplicate_phase_bound = int(duplicate_stats["phase_queries_per_draw_bound"])
        duplicate_max_amp_error = 0.
        duplicate_max_norm_error = 0.
        duplicate_max_components = 0
        for work in range(PERIOD):
            charge(baseline_budget, "row_local_terms", duplicate_row_pair_bound)
            charge(baseline_budget, "phase_queries", duplicate_phase_bound)
            charge(baseline_budget, "expansion_terms", H * L)
            live_start = len(baseline_live)
            try:
                row = duplicate.row(work)
            finally:
                for live in baseline_live[live_start:]:
                    for key, value in live.items():
                        baseline_counters[key] += int(value)
            expanded = {}
            for start, count, gamma in row["components"]:
                for n in range(int(count)):
                    charge(baseline_counters, "expansion_terms")
                    exponent = int(start) + int(row["stride"]) * n
                    if exponent in expanded:
                        raise AssertionError("duplicate baseline row overlaps")
                    expanded[exponent] = complex(gamma)
            duplicate_max_components = max(duplicate_max_components,
                                           len(row["components"]))
            charge(counters, "norm_entries", Q)
            expected_norm = float(np.sum(np.abs(direct_for_duplicate[:, work]) ** 2))
            duplicate_max_norm_error = max(
                duplicate_max_norm_error, abs(float(row["norm"]) - expected_norm))
            charge(counters, "comparison_entries", Q)
            for exponent in range(Q):
                duplicate_max_amp_error = max(
                    duplicate_max_amp_error,
                    abs(expanded.get(exponent, 0j)
                        - direct_for_duplicate[exponent, work]))
        duplicate_ok = (duplicate_max_amp_error <= TOL
                        and duplicate_max_norm_error <= TOL)
        c1 = wrong_witness is not None
        p1 = bool(max_amp_error <= TOL and max_norm_error <= TOL
                  and all(entry["max_norm_error"] <= TOL
                          for entry in summaries.values())
                  and duplicate_ok)
        p2 = bool(max_fourier_error <= TOL and max_probability_error <= TOL)
        # The wrong-left control is deliberately charged in separate
        # preflight categories.  Recombine those categories before checking
        # the cumulative actual counters, since the live counters aggregate
        # ordinary and control rows under the same names.
        category_bounds = pre["category_bounds"]
        combined_bounds = {
            "group_pair_visits": category_bounds["group_pair_visits"]
            + category_bounds["wrong_group_pair_visits"],
            "coefficient_pair_visits": category_bounds["coefficient_pair_visits"]
            + category_bounds["wrong_coefficient_pair_visits"],
            "phase_calls": category_bounds["phase_calls"]
            + category_bounds["wrong_phase_calls"],
            "phase_products": category_bounds["phase_products"]
            + category_bounds["wrong_phase_products"],
            "expansion_terms": category_bounds["expansion_terms"]
            + category_bounds["wrong_expansion_terms"],
            "fourier_terms": category_bounds["fourier_terms"],
            "comparison_entries": category_bounds["comparison_entries"],
            "norm_entries": category_bounds["norm_entries"],
        }
        category_reconciled = all(
            counters.get(name, 0) <= limit
            for name, limit in combined_bounds.items())
        row_budget_reconciled = all(
            row_budget.get(name, 0) <= combined_bounds[name]
            for name in combined_bounds)
        fft_reconciled = (
            fft_budget.get("fft_calls", 0) <= category_bounds["fft_calls"]
            and fft_budget.get("fft_entries", 0) <= category_bounds["fft_entries"]
            and counters.get("fft_calls", 0) <= category_bounds["fft_calls"]
            and counters.get("fft_entries", 0) <= category_bounds["fft_entries"])
        reference_reconciled = all(
            reference_budget.get(name, 0) <= category_bounds[name]
            for name in ("reference_calls", "direct_block_products",
                         "direct_phase_queries", "direct_shift_ops",
                         "modular_pow_queries")) and all(
            reference_counters.get(name, 0) <= reference_budget.get(name, 0)
            for name in ("direct_block_products", "direct_phase_queries",
                         "direct_shift_ops", "modular_pow_queries"))
        # ``_checked_phase`` is the shared phase-query counter and includes
        # early calls; early_phase_queries is a named subset, not an additive
        # cost category.
        baseline_phase_actual = baseline_counters.get("phase_queries", 0)
        baseline_reconciled = (
            baseline_counters.get("row_local_terms", 0)
            <= category_bounds["duplicate_baseline_row_terms"]
            and baseline_phase_actual
            <= category_bounds["duplicate_baseline_phase_queries"]
            and baseline_counters.get("early_phase_queries", 0)
            <= category_bounds["duplicate_baseline_phase_queries"]
            and baseline_counters.get("expansion_terms", 0)
            <= category_bounds["duplicate_baseline_expansion_terms"]
            and baseline_counters.get("phase_squares", 0)
            <= category_bounds["duplicate_baseline_phase_squares"]
            and baseline_budget.get("row_local_terms", 0)
            <= category_bounds["duplicate_baseline_row_terms"]
            and baseline_budget.get("phase_queries", 0)
            <= category_bounds["duplicate_baseline_phase_queries"]
            and baseline_budget.get("expansion_terms", 0)
            <= category_bounds["duplicate_baseline_expansion_terms"])
        setup_reconciled = (
            setup_counters.get("unitary_check_terms", 0)
            <= category_bounds["unitary_check_terms"]
            and setup_counters.get("matrix_entries", 0)
            <= category_bounds["matrix_entries"])
        # This checksum intentionally overlaps modular_pow_queries with the
        # phase-query counters; it is a conservative ledger bound, not a
        # claim of distinct physical cost.
        actual_counter_checksum = (sum(counters.values())
                                   + sum(reference_counters.values())
                                   + sum(baseline_counters.values())
                                   + sum(setup_counters.values()))
        p3 = bool(all(value <= MAX_CATEGORY_TERMS
                      for value in counters.values())
                  and all(value <= MAX_CATEGORY_TERMS
                          for value in reference_counters.values())
                  and all(value <= MAX_CATEGORY_TERMS
                          for value in reference_budget.values())
                  and all(value <= MAX_CATEGORY_TERMS
                          for value in baseline_counters.values())
                  and row_budget_reconciled and fft_reconciled
                  and reference_reconciled and baseline_reconciled
                  and setup_reconciled
                  and category_reconciled
                  and actual_counter_checksum <= pre["total_reserved_terms"]
                  and actual_counter_checksum <= MAX_TOTAL_RESERVED)
        report.update({
            "fixture": {"N": N, "base": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L,
                        "phase_positions": list(POSITIONS)},
            "summaries": summaries,
            "max_amplitude_error": max_amp_error,
            "max_norm_error": max_norm_error,
            "max_literal_column_norm_error": max_literal_column_norm_error,
            "work_marginal_totals": work_masses,
            "max_selected_fourier_error": max_fourier_error,
            "max_selected_probability_error": max_probability_error,
            "wrong_left_witness": wrong_witness,
            "counters": dict(counters),
            "reference_counters": dict(reference_counters),
            "reference_budget": dict(reference_budget),
            "row_budget": dict(row_budget),
            "fft_budget": dict(fft_budget),
            "setup_counters": dict(setup_counters),
            "baseline_counters": dict(baseline_counters),
            "baseline_budget": dict(baseline_budget),
            "duplicate_baseline_stats": duplicate_stats,
            "duplicate_baseline_max_amplitude_error": duplicate_max_amp_error,
            "duplicate_baseline_max_norm_error": duplicate_max_norm_error,
            "duplicate_baseline_max_components": duplicate_max_components,
            "duplicate_baseline_ok": duplicate_ok,
            "cut_counter_totals": {key: dict(value)
                                    for key, value in cut_counter_totals.items()},
            "work0_counter_examples": work0_counter_examples,
            "actual_counter_checksum_overlapping": actual_counter_checksum,
            "combined_category_bounds": combined_bounds,
            "category_reconciled": category_reconciled,
            "row_budget_reconciled": row_budget_reconciled,
            "fft_reconciled": fft_reconciled,
            "reference_reconciled": reference_reconciled,
            "baseline_reconciled": baseline_reconciled,
            "setup_reconciled": setup_reconciled,
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
            "status": "PASS" if p1 and p2 and p3 and c1 else "FAIL",
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "scope": "coherent amplitudes/selected normalized laws; no sampler",
        })
        exp.check("P1", p1, "cut amplitudes, disjoint expansion, and norms")
        exp.check("P2", p2, "selected geometric amplitudes/probabilities")
        exp.check("P3", p3, "named category and total workload caps")
        exp.fail_check("C1", c1, "wrong frozen-left coefficient has witness")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        report.update({
            "counters": dict(counters),
            "reference_counters": dict(reference_counters),
            "reference_budget": dict(reference_budget),
            "row_budget": dict(row_budget),
            "fft_budget": dict(fft_budget),
            "setup_counters": dict(setup_counters),
            "baseline_counters": dict(baseline_counters),
            "baseline_budget": dict(baseline_budget),
        })
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before amplitude audit")
        exp.fail_check("C1", False, "exception before wrong-left control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "fixture": "r60,b3,t8,s7, two G1 phases at (3,v)",
        "scope": "bounded coherent amplitude pilot",
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_category_terms": MAX_CATEGORY_TERMS,
        "max_total_reserved_terms": MAX_TOTAL_RESERVED,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
