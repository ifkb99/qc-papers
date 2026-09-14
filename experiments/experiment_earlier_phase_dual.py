"""Independent dual-cover pilot for the TODO39 earlier-phase schedule.

The frozen indexed fixture is N=61, a=2, r=60, b=3, t=8, s=7 with the
same W0/W1 and pointwise phases used by experiment_earlier_phase_cycles.  For
each v=0..7, this compares the existing cyclic-refined row with the
complementary partition l=a+A*k, A=2**v, K=L/A.  Literal columns are built
through the existing ER direct_column reference; no generic propagator is
introduced.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Both row covers reconstruct every literal column amplitude, preserve
      disjoint nonwrapping exponent supports, and reproduce row norms.
  P2  Full Fourier laws from the direct columns and both covers agree.  The
      measured cycle/dual component counts are recorded without assuming one
      cover is always smaller.
  P3  The proposed conservative component-count bound is checked explicitly;
      if it is too small for the dual partition, retain that counterexample
      rather than weakening the fixture silently.
  C1  Dropping the dual partition's -A*k in the early phase must change an
      explicit interior coefficient witness (v=3, j=0), if that control is
      genuinely discriminating.

The dual coefficient for rho'=(c+q-A*k-u) mod r uses
f((c+q-A*k) mod r), then has start L*h+A*k+rho', stride r, and the
nonwrapping count inside the length-A interval.  All dense arrays and named
local/phase/expansion work are preflighted before allocation.  This is a
bounded float amplitude/Fourier diagnostic, not a sampler or efficiency claim.
"""
from __future__ import annotations

import math
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


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT = 8, 7
Q, L, H = 1 << WIDTH, 1 << SPLIT, 1 << (WIDTH - SPLIT)
INSERTIONS = tuple(range(SPLIT + 1))
R = min(PERIOD, 2 * BLOCK - 1)
BYTE_CAP = 16 * 1024 * 1024
TERM_CAP = 5_000_000
TOL = 3e-10


def report_path() -> Path:
    path = Path("out") / (
        "earlier_phase_dual_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
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


def add_count(counters: dict, key: str, amount: int, cap: int = TERM_CAP,
              label: str | None = None) -> None:
    amount = int(amount)
    if amount < 0 or counters[key] + amount > cap:
        raise MemoryError(f"{label or key} cap before loop/call")
    counters[key] += amount


def dual_phase(index: int, counters: dict) -> complex:
    index = int(index) % PERIOD
    add_count(counters, "dual_phase_queries", 1, label="dual phase")
    add_count(counters, "dual_modular_power_queries", 1,
              label="dual modular power")
    return unit_phase(pow(BASE, index, N), N)


def dual_row(work: int, insertion: int, W0: np.ndarray, W1: np.ndarray,
             counters: dict, *, wrong_shift: bool = False):
    """Complementary l=a+A*k row, including zero/canceled candidates."""
    A, K = 1 << insertion, L >> insertion
    components = []
    for high in range(H):
        v = (int(work) - L * high) % PERIOD
        cell, p = divmod(v, BLOCK)
        late = dual_phase(v, counters)
        for k in range(K):
            shift = A * k
            pairs = {}
            for u in range(BLOCK):
                for q in range(BLOCK):
                    add_count(counters, "dual_local_terms", 1,
                              label="dual local pair")
                    rho = (cell * BLOCK + q - shift - u) % PERIOD
                    pairs.setdefault(rho, []).append(
                        (u, q, W1[p, q] * W0[u, 0]))
            for rho, local_pairs in sorted(pairs.items()):
                count = max(0, 1 + (A - 1 - rho) // PERIOD)
                if count == 0:
                    continue
                gamma = 0j
                for _u, q, amplitude in local_pairs:
                    # The wrong control deliberately omits -A*k only in the
                    # phase argument; rho' and support remain correct.
                    phase_index = (cell * BLOCK + q
                                   - (0 if wrong_shift else shift)) % PERIOD
                    gamma += amplitude * dual_phase(phase_index, counters)
                gamma *= late
                if gamma != 0:
                    components.append((L * high + shift + rho,
                                       PERIOD, count, gamma))
    return components


def expand_components(components: list[tuple], counters: dict,
                      key: str = "dual_expansion_terms") -> dict[int, complex]:
    result = {}
    for start, stride, count, gamma in components:
        for n in range(int(count)):
            add_count(counters, key, 1, label=key)
            exponent = int(start) + int(stride) * n
            if exponent in result:
                raise AssertionError("duplicate component exponent")
            if not 0 <= exponent < Q:
                raise AssertionError("component leaves exponent register")
            result[exponent] = gamma
    return result


def law_audit(law: np.ndarray) -> dict:
    law = np.asarray(law, dtype=float)
    finite = bool(np.all(np.isfinite(law)))
    minimum = float(np.min(law)) if law.size else float("nan")
    total = float(np.sum(law)) if law.size else 0.
    return {
        "shape": list(law.shape),
        "shape_ok": bool(law.shape == (Q,)),
        "finite": finite,
        "minimum": minimum,
        "nonnegative": bool(finite and minimum >= -TOL),
        "sum": total,
        "normalized": bool(finite and abs(total - 1.) < TOL),
    }


def column_audit(matrix: np.ndarray) -> dict:
    # Rows are the literal indexed branch columns phi_e; storage is
    # (exponent, work-label), so normalization is along axis 1.
    norms = np.sum(np.abs(matrix) ** 2, axis=1)
    finite = bool(np.all(np.isfinite(norms)))
    return {
        "shape": list(matrix.shape),
        "finite": finite,
        "minimum_norm": float(np.min(norms)),
        "maximum_norm": float(np.max(norms)),
        "max_norm_error": float(np.max(np.abs(norms - 1.))),
        "normalized": bool(finite and np.max(np.abs(norms - 1.)) < TOL),
    }


def ensure_budget(counters: dict, bounds: dict, additions: dict, label: str) -> None:
    """Check cumulative ER/dual work before invoking a loop or call."""
    for key, amount in additions.items():
        if counters[key] + int(amount) > bounds[key]:
            raise MemoryError(f"{label}: {key} preflight cap")


def preflight() -> dict:
    branch_entries = Q * PERIOD
    components = {
        "direct_cycle_dual_column_arrays": 3 * branch_entries * 16,
        "FFT_complex_and_absolute_temporaries": 2 * branch_entries * 16,
        "retained_direct_cycle_dual_geometric_laws": 4 * len(INSERTIONS) * Q * 8,
        "geometric_complex_accumulator": Q * 16,
        "repeated_work_blocks_and_setup": 4 * PERIOD * PERIOD * 16,
        "component_tuples_and_scalar_lists": 3 * len(INSERTIONS) * PERIOD * Q * 16,
        "numeric_safety_reserve": 512 * 1024,
    }
    payload = guard_bytes(sum(components.values()), "dual-cover numeric preflight")
    claimed_bound = H * R * math.ceil(L / PERIOD)
    safe_bound = H * L
    sum_k = sum(L >> insertion for insertion in INSERTIONS)
    max_cycle_period = max(er.generic_phase_period(v) for v in INSERTIONS)
    wrong_k = L >> 3
    operation_bounds = {
        "direct_column_storage_assignments": len(INSERTIONS) * Q,
        "direct_phase_queries": len(INSERTIONS) * Q * 2 * PERIOD,
        "direct_block_products": len(INSERTIONS) * Q * (PERIOD // BLOCK) * BLOCK**2,
        "direct_shift_ops": len(INSERTIONS) * Q * 3,
        "row_local_terms": len(INSERTIONS) * PERIOD * H * BLOCK**2,
        "early_phase_terms": len(INSERTIONS) * PERIOD * H * BLOCK**2 * max_cycle_period,
        "early_phase_queries": len(INSERTIONS) * PERIOD * H * BLOCK**2 * max_cycle_period,
        "late_phase_queries": len(INSERTIONS) * PERIOD * H,
        "modular_pow_queries": (len(INSERTIONS) * Q * 2 * PERIOD
                                + len(INSERTIONS) * PERIOD * H
                                  * (BLOCK**2 * max_cycle_period + 1)),
        "expansion_terms": len(INSERTIONS) * PERIOD * Q,
        "dual_local_terms": H * PERIOD * BLOCK**2 * sum_k + H * wrong_k * BLOCK**2,
        "dual_phase_queries": (H * PERIOD * (sum_k * BLOCK**2 + len(INSERTIONS))
                               + H * (wrong_k * BLOCK**2 + 1)),
        "dual_modular_power_queries": (H * PERIOD * (sum_k * BLOCK**2 + len(INSERTIONS))
                                       + H * (wrong_k * BLOCK**2 + 1)),
        "dual_expansion_terms": len(INSERTIONS) * PERIOD * Q,
        "dual_wrong_expansion_terms": Q,
        "cycle_matrix_assignments": len(INSERTIONS) * PERIOD * Q,
        "dual_matrix_assignments": len(INSERTIONS) * PERIOD * Q,
        "fft_calls": len(INSERTIONS) * 3,
        "fft_input_entries": len(INSERTIONS) * 3 * branch_entries,
        "fft_output_entries": len(INSERTIONS) * 3 * Q,
        "geometric_component_terms": len(INSERTIONS) * PERIOD * Q * claimed_bound,
        "geometric_sum_queries": len(INSERTIONS) * PERIOD * Q * claimed_bound,
        "geometric_phase_queries": len(INSERTIONS) * PERIOD * Q * claimed_bound,
    }
    for key, value in operation_bounds.items():
        if value > TERM_CAP and key not in {"fft_calls", "fft_input_entries",
                                            "fft_output_entries"}:
            raise MemoryError(f"{key} operation preflight exceeds cap")
    if operation_bounds["fft_calls"] > len(INSERTIONS) * 3:
        raise MemoryError("FFT call preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": payload,
        "payload_components": components,
        "claimed_component_bound_per_row": claimed_bound,
        "safe_distinct_exponent_bound_per_row": safe_bound,
        "claimed_root_scalar_bound": len(INSERTIONS) * PERIOD * Q * claimed_bound,
        "safe_root_scalar_bound": len(INSERTIONS) * PERIOD * Q * safe_bound,
        "operation_bounds": operation_bounds,
        "sum_dual_K": sum_k,
        "max_cycle_period": max_cycle_period,
        "max_numeric_payload_bytes": BYTE_CAP,
        "max_scalar_terms": TERM_CAP,
    }


def law_from_matrix(matrix: np.ndarray, counters: dict, bounds: dict,
                    label: str) -> np.ndarray:
    add_count(counters, "fft_calls", 1, bounds["fft_calls"], label=f"{label} FFT call")
    add_count(counters, "fft_input_entries", matrix.size,
              bounds["fft_input_entries"], label=f"{label} FFT input")
    # The branch matrix contains normalized phi_e columns.  The uniform
    # exponent amplitude and unitary QFT together give FFT(matrix)/Q.
    add_count(counters, "fft_output_entries", Q,
              bounds["fft_output_entries"], label=f"{label} FFT output")
    transform = np.fft.fft(matrix, axis=0) / Q
    result = np.sum(np.abs(transform) ** 2, axis=1)
    return result


def geometric_law(component_rows: list[list[tuple]], counters: dict,
                  bounds: dict) -> np.ndarray:
    """Evaluate the dual nonwrapping components without expanding labels."""
    result = np.zeros(Q, dtype=float)
    for components in component_rows:
        for output in range(Q):
            total = 0j
            for start, stride, count, gamma in components:
                add_count(counters, "geometric_component_terms", 1,
                          bounds["geometric_component_terms"], label="geometric term")
                add_count(counters, "geometric_phase_queries", 1,
                          bounds["geometric_phase_queries"], label="geometric phase")
                add_count(counters, "geometric_sum_queries", 1,
                          bounds["geometric_sum_queries"], label="geometric sum")
                total += gamma * unit_phase(-output * int(start), Q) * geometric_sum(
                    int(count), -output * int(stride), Q)
            result[output] += abs(total) ** 2 / (Q * Q)
    return result


def tv(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.sum(np.abs(first - second)) / 2.)


def main() -> None:
    exp = Experiment("earlier_phase_dual", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "cycle and dual rows reconstruct literal amplitudes")
    exp.predict("P2", "both covers agree on complete Fourier laws")
    exp.predict("P3", "the stated dual component bound is explicitly audited")
    exp.must_fail("C1", "omitting -A*k from the dual phase changes an interior witness")
    started = time.perf_counter()
    counters = {
        "column_local_terms": 0, "row_local_terms": 0,
        "early_phase_terms": 0, "late_phase_queries": 0,
        "early_phase_queries": 0, "direct_phase_queries": 0,
        "formula_phase_queries": 0, "modular_pow_queries": 0,
        "direct_block_products": 0, "direct_shift_ops": 0,
        "expansion_terms": 0,
        "dual_local_terms": 0, "dual_phase_queries": 0,
        "dual_modular_power_queries": 0, "dual_expansion_terms": 0,
        "dual_wrong_expansion_terms": 0,
        "direct_column_storage_assignments": 0,
        "cycle_matrix_assignments": 0, "dual_matrix_assignments": 0,
        "fft_calls": 0, "fft_input_entries": 0, "fft_output_entries": 0,
        "geometric_component_terms": 0, "geometric_sum_queries": 0,
        "geometric_phase_queries": 0,
    }
    report = {"status": "FAIL", "rows": [], "predictions": {
        "covers": "cycle and dual both exact",
        "count_comparison": "measured, no dominance assumed",
    }}
    p1 = p2 = p3 = c1 = False
    try:
        report["fixture"] = {
            "N": N, "base": BASE, "period": PERIOD, "block_size": BLOCK,
            "width": WIDTH, "split": SPLIT, "Q": Q, "L": L, "H": H,
            "insertions": list(INSERTIONS), "R": R,
        }
        report["preflight"] = preflight()
        bounds = report["preflight"]["operation_bounds"]
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        if not er.valid_unitary(W0) or not er.valid_unitary(W1):
            raise AssertionError("invalid frozen work blocks")
        all_rows = []
        direct_laws, cycle_laws, dual_laws = {}, {}, {}
        geometric_laws = {}
        law_audits = {}
        column_audits = {}
        count_rows = {}
        max_cycle_column_error = max_dual_column_error = 0.
        max_direct_norm_error = max_cycle_column_norm_error = max_dual_column_norm_error = 0.
        max_cycle_error = max_dual_error = 0.
        max_cycle_norm_error = max_dual_norm_error = 0.
        wrong_witness = None
        for insertion in INSERTIONS:
            guard_shape((Q, PERIOD), label=f"v={insertion} direct columns")
            guard_shape((Q, PERIOD), label=f"v={insertion} dual columns")
            direct = np.zeros((Q, PERIOD), dtype=np.complex128)
            dual = np.zeros((Q, PERIOD), dtype=np.complex128)
            cycle = np.zeros((Q, PERIOD), dtype=np.complex128)
            for exponent in range(Q):
                ensure_budget(counters, bounds, {
                    "direct_phase_queries": 2 * PERIOD,
                    "modular_pow_queries": 2 * PERIOD,
                    "direct_block_products": (PERIOD // BLOCK) * BLOCK**2,
                    "direct_shift_ops": 3,
                    "direct_column_storage_assignments": 1,
                }, "ER direct column")
                counters["direct_column_storage_assignments"] += 1
                direct[exponent] = er.direct_column(
                    exponent, insertion, W0, W1, counters)
            direct_column_report = column_audit(direct)
            column_norm_error = direct_column_report["max_norm_error"]
            row_errors = []
            cycle_errors = []
            dual_errors = []
            cycle_norm_errors = []
            dual_norm_errors = []
            cycle_counts, dual_counts = [], []
            dual_component_rows = []
            disjoint_ok = True
            for work in range(PERIOD):
                direct_row = direct[:, work]
                cycle_period = er.generic_phase_period(insertion)
                ensure_budget(counters, bounds, {
                    "row_local_terms": H * BLOCK**2,
                    "early_phase_terms": H * BLOCK**2 * cycle_period,
                    "early_phase_queries": H * BLOCK**2 * cycle_period,
                    "late_phase_queries": H,
                    "modular_pow_queries": H * (BLOCK**2 * cycle_period + 1),
                    "expansion_terms": Q,
                }, "ER refined row")
                cycle_components = er.refined_row(work, insertion, W0, W1, counters)
                cycle_expanded = er.expand_row(cycle_components, counters)
                A, K = 1 << insertion, L >> insertion
                ensure_budget(counters, bounds, {
                    "dual_local_terms": H * K * BLOCK**2,
                    "dual_phase_queries": H * (K * BLOCK**2 + 1),
                    "dual_modular_power_queries": H * (K * BLOCK**2 + 1),
                    "dual_expansion_terms": Q,
                }, "dual row")
                dual_components = dual_row(work, insertion, W0, W1, counters)
                dual_expanded = expand_components(dual_components, counters)
                dual_component_rows.append(dual_components)
                cycle_support = set(cycle_expanded)
                dual_support = set(dual_expanded)
                disjoint_ok &= (len(cycle_support) == sum(
                    int(c[2]) for c in cycle_components))
                disjoint_ok &= (len(dual_support) == sum(
                    int(c[2]) for c in dual_components))
                ensure_budget(counters, bounds, {
                    "cycle_matrix_assignments": Q,
                    "dual_matrix_assignments": Q,
                }, "row matrix storage")
                for exponent in range(Q):
                    counters["cycle_matrix_assignments"] += 1
                    counters["dual_matrix_assignments"] += 1
                    cycle[exponent, work] = cycle_expanded.get(exponent, 0j)
                    dual[exponent, work] = dual_expanded.get(exponent, 0j)
                row_errors.append(max(abs(cycle_expanded.get(e, 0j)
                                         - direct_row[e]) for e in range(Q)))
                dual_errors.append(max(abs(dual_expanded.get(e, 0j)
                                          - direct_row[e]) for e in range(Q)))
                cycle_norm = sum(int(c[2]) * abs(c[3])**2
                                 for c in cycle_components)
                dual_norm = sum(int(c[2]) * abs(c[3])**2
                                for c in dual_components)
                direct_norm = float(np.sum(abs(direct_row)**2))
                cycle_norm_errors.append(abs(cycle_norm - direct_norm))
                dual_norm_errors.append(abs(dual_norm - direct_norm))
                cycle_counts.append(len(cycle_components))
                dual_counts.append(len(dual_components))
                if insertion == 3 and work == 0 and wrong_witness is None:
                    wrong_k = L >> 3
                    ensure_budget(counters, bounds, {
                        "dual_local_terms": H * wrong_k * BLOCK**2,
                        "dual_phase_queries": H * (wrong_k * BLOCK**2 + 1),
                        "dual_modular_power_queries": H * (wrong_k * BLOCK**2 + 1),
                        "dual_wrong_expansion_terms": Q,
                    }, "wrong dual phase control")
                    wrong_components = dual_row(work, insertion, W0, W1,
                                                counters, wrong_shift=True)
                    wrong_expanded = expand_components(
                        wrong_components, counters, key="dual_wrong_expansion_terms")
                    for exponent in range(Q):
                        if abs(wrong_expanded.get(exponent, 0j)
                               - direct_row[exponent]) > 1e-8:
                            wrong_witness = {
                                "insertion": insertion, "work": work,
                                "exponent": exponent,
                                "high": int(exponent // L),
                                "k": int((exponent % L) // (1 << insertion)),
                                "a": int((exponent % L) % (1 << insertion)),
                                "correct": [float(dual_expanded.get(exponent, 0j).real),
                                             float(dual_expanded.get(exponent, 0j).imag)],
                                "wrong": [float(wrong_expanded.get(exponent, 0j).real),
                                           float(wrong_expanded.get(exponent, 0j).imag)],
                                "direct": [float(direct_row[exponent].real),
                                            float(direct_row[exponent].imag)],
                                "absolute_error": float(abs(
                                    wrong_expanded.get(exponent, 0j)
                                    - direct_row[exponent])),
                            }
                            break
            cycle_component_report = column_audit(cycle)
            dual_column_report = column_audit(dual)
            direct_law = law_from_matrix(direct, counters, bounds, "direct")
            cycle_law = law_from_matrix(cycle, counters, bounds, "cycle")
            dual_law = law_from_matrix(dual, counters, bounds, "dual")
            geometric_terms = Q * sum(len(row) for row in dual_component_rows)
            ensure_budget(counters, bounds, {
                "geometric_component_terms": geometric_terms,
                "geometric_sum_queries": geometric_terms,
                "geometric_phase_queries": geometric_terms,
            }, "geometric dual law")
            geometric_law_value = geometric_law(dual_component_rows, counters, bounds)
            direct_laws[str(insertion)] = direct_law.tolist()
            cycle_laws[str(insertion)] = cycle_law.tolist()
            dual_laws[str(insertion)] = dual_law.tolist()
            geometric_laws[str(insertion)] = geometric_law_value.tolist()
            law_audits[str(insertion)] = {
                "direct": law_audit(direct_law),
                "cycle": law_audit(cycle_law),
                "dual": law_audit(dual_law),
                "geometric_dual": law_audit(geometric_law_value),
            }
            column_audits[str(insertion)] = {
                "direct": direct_column_report,
                "cycle": cycle_component_report,
                "dual": dual_column_report,
            }
            max_cycle_column_error = max(max_cycle_column_error,
                                         float(np.max(abs(cycle - direct))))
            max_dual_column_error = max(max_dual_column_error,
                                        float(np.max(abs(dual - direct))))
            max_direct_norm_error = max(max_direct_norm_error,
                                        direct_column_report["max_norm_error"])
            max_cycle_column_norm_error = max(max_cycle_column_norm_error,
                                              cycle_component_report["max_norm_error"])
            max_dual_column_norm_error = max(max_dual_column_norm_error,
                                             dual_column_report["max_norm_error"])
            max_cycle_error = max(max_cycle_error, max(row_errors))
            max_dual_error = max(max_dual_error, max(dual_errors))
            max_cycle_norm_error = max(max_cycle_norm_error, max(cycle_norm_errors))
            max_dual_norm_error = max(max_dual_norm_error, max(dual_norm_errors))
            claimed = report["preflight"]["claimed_component_bound_per_row"]
            count_rows[str(insertion)] = {
                "cycle_counts": cycle_counts,
                "dual_counts": dual_counts,
                "cycle_max_components": max(cycle_counts),
                "dual_max_components": max(dual_counts),
                "claimed_bound_per_row": claimed,
                "cycle_within_claimed_bound": max(cycle_counts) <= claimed,
                "dual_within_claimed_bound": max(dual_counts) <= claimed,
                "cycle_structural_bound": H * min(L, R * er.effective_phase_period(insertion)),
                "dual_structural_bound": H * K * min(A, R),
                "dual_support_disjoint": disjoint_ok,
                "cycle_law_tv_vs_direct": tv(cycle_law, direct_law),
                "dual_law_tv_vs_direct": tv(dual_law, direct_law),
                "geometric_law_tv_vs_direct": tv(geometric_law_value, direct_law),
            }
            all_rows.append({
                "insertion": insertion,
                "max_cycle_row_error": max(row_errors),
                "max_dual_row_error": max(dual_errors),
                "max_cycle_norm_error": max(cycle_norm_errors),
                "max_dual_norm_error": max(dual_norm_errors),
                "cycle_component_max": max(cycle_counts),
                "dual_component_max": max(dual_counts),
                "cycle_law_tv": tv(cycle_law, direct_law),
                "dual_law_tv": tv(dual_law, direct_law),
                "geometric_law_tv": tv(geometric_law_value, direct_law),
            })
            del direct, cycle, dual
        claimed_holds = all(
            row["cycle_within_claimed_bound"] and row["dual_within_claimed_bound"]
            and row["cycle_max_components"] <= row["cycle_structural_bound"]
            and row["dual_max_components"] <= row["dual_structural_bound"]
            for row in count_rows.values())
        p1 = bool(max_cycle_error < TOL and max_dual_error < TOL
                  and max_cycle_norm_error < TOL and max_dual_norm_error < TOL
                  and max_cycle_column_error < TOL and max_dual_column_error < TOL
                  and max_direct_norm_error < TOL
                  and max_cycle_column_norm_error < TOL
                  and max_dual_column_norm_error < TOL
                  and all(audit["finite"] and audit["normalized"]
                          for per_insertion in column_audits.values()
                          for audit in per_insertion.values())
                  and all(row["dual_support_disjoint"] for row in count_rows.values()))
        p2 = bool(all(row["cycle_law_tv"] < TOL
                      and row["dual_law_tv"] < TOL
                      and row["geometric_law_tv"] < TOL
                      and all(law_audits[str(row["insertion"])][name]["shape_ok"]
                              and law_audits[str(row["insertion"])][name]["finite"]
                              and law_audits[str(row["insertion"])][name]["normalized"]
                              and law_audits[str(row["insertion"])][name]["nonnegative"]
                              for name in ("direct", "cycle", "dual", "geometric_dual"))
                      for row in all_rows))
        operation_keys = tuple(bounds)
        operation_bound_checks = {
            key: {
                "actual": int(counters[key]),
                "bound": int(bounds[key]),
                "ok": bool(counters[key] <= bounds[key]),
            }
            for key in operation_keys
        }
        cost_ok = all(item["ok"] for item in operation_bound_checks.values())
        exact_count_checks = {
            "fft_calls": counters["fft_calls"] == bounds["fft_calls"],
            "fft_input_entries": counters["fft_input_entries"] == bounds["fft_input_entries"],
            "direct_column_storage_assignments": (
                counters["direct_column_storage_assignments"]
                == bounds["direct_column_storage_assignments"]),
            "cycle_matrix_assignments": counters["cycle_matrix_assignments"] == bounds["cycle_matrix_assignments"],
            "dual_matrix_assignments": counters["dual_matrix_assignments"] == bounds["dual_matrix_assignments"],
            "fft_output_entries": counters["fft_output_entries"] == 3 * len(INSERTIONS) * Q,
            "row_local_terms": counters["row_local_terms"] == bounds["row_local_terms"],
            "dual_local_terms": counters["dual_local_terms"] == bounds["dual_local_terms"],
            "direct_block_products": counters["direct_block_products"] == bounds["direct_block_products"],
            "modular_powers": counters["modular_pow_queries"] == (
                counters["early_phase_queries"] + counters["late_phase_queries"]
                + counters["direct_phase_queries"]),
            "geometric_terms": counters["geometric_component_terms"] == (
                Q * sum(sum(row["dual_counts"]) for row in count_rows.values())),
            "geometric_subcounts": counters["geometric_component_terms"] ==
                counters["geometric_sum_queries"] == counters["geometric_phase_queries"],
        }
        p3 = bool(claimed_holds and cost_ok
                  and all(exact_count_checks.values()))
        c1 = bool(wrong_witness is not None)
        p1 = bool(p1 and cost_ok)
        actual_counter_sum_with_overlaps = int(sum(counters.values()))
        report.update({
            "rows": all_rows,
            "component_counts": count_rows,
            "direct_full_laws": direct_laws,
            "cycle_full_laws": cycle_laws,
            "dual_full_laws": dual_laws,
            "geometric_dual_full_laws": geometric_laws,
            "law_audits": law_audits,
            "column_audits": column_audits,
            "max_cycle_column_error": max_cycle_column_error,
            "max_dual_column_error": max_dual_column_error,
            "max_direct_column_norm_error": max_direct_norm_error,
            "max_cycle_column_norm_error": max_cycle_column_norm_error,
            "max_dual_column_norm_error": max_dual_column_norm_error,
            "max_cycle_row_error": max_cycle_error,
            "max_dual_row_error": max_dual_error,
            "max_cycle_norm_error": max_cycle_norm_error,
            "max_dual_norm_error": max_dual_norm_error,
            "claimed_component_bound_holds": claimed_holds,
            "wrong_phase_witness": wrong_witness,
            "counters": counters,
            "operation_bound_checks": operation_bound_checks,
            "exact_count_checks": exact_count_checks,
            "actual_counter_sum_with_overlaps": actual_counter_sum_with_overlaps,
            "cost_within_cap": cost_ok,
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "cycle and dual rows reconstruct literal amplitudes")
        exp.check("P2", p2, "cycle/dual full Fourier laws agree")
        exp.check("P3", p3, "claimed per-row component bound holds")
        exp.fail_check("C1", c1, "wrong -A*k phase changes the fixed interior witness")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before dual-cover audit")
        exp.check("P2", False, "exception before dual-cover audit")
        exp.check("P3", False, "exception before dual-cover audit")
        exp.fail_check("C1", False, "exception before wrong-phase control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": report.get("fixture"),
        "reference": "existing ER direct_column/refined_row plus independent dual rows",
        "full_output": "FFT of literal indexed columns and both covers",
        "no_generic_propagator": True,
        "max_numeric_payload_bytes": BYTE_CAP,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
