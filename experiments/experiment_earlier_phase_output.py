"""Complete-output comparison for the earlier-phase C79 schedule.

The fixture is the frozen N=61, a=2, r=60, b=3, t=8, split s=7 schedule
with W0=work_block(pi/4), W1=work_block(-pi/10), late phase g, and a second
identical phase inserted after v=0,...,7 low controls.  Existing ER direct
columns and refined rows are used.  This experiment compares their complete
256-output laws, without constructing a production sampler.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The tiny direct-column FFT and refined geometric-component laws agree
      as complete output distributions for every insertion position.
  P2  Every target, omitted-earlier-phase, and component-only law is finite,
      nonnegative, and normalized; the target first output bit is uniform.
  P3  The report records omission and component-only TV at thresholds 1e-3
      and 1e-2 without predicting which classification they take.
  C1  Dropping coherent Fourier cross-component terms changes at least one
      complete output law.

This is a bounded float law check, not a numerical certificate, timing claim,
sampler implementation, or hardness result.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.fourier_sampling import geometric_sum, unit_phase
from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_earlier_phase_cycles import (
    direct_column, refined_row, generic_phase_period, effective_phase_period,
)


N, BASE, PERIOD, BLOCK, WIDTH, SPLIT = 61, 2, 60, 3, 8, 7
Q, L = 1 << WIDTH, 1 << SPLIT
MAX_BYTES = 16 << 20
MAX_ROOT_TERMS = 5_000_000
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"earlier_phase_output_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def valid_law(values):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (Q,) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right):
    return float(np.sum(np.abs(np.asarray(left)-np.asarray(right))) / 2)


def preflight():
    # Direct columns, FFT output/temporary, and comparison magnitudes are
    # retained only for one insertion at a time; output laws are retained for
    # all eight positions. Python report containers are count-bounded below.
    matrix = Q * PERIOD * np.dtype(np.complex128).itemsize
    numeric = 5 * matrix + 32 * PERIOD * np.dtype(np.complex128).itemsize
    # jsonable temporarily retains both the old and converted scalar lists.
    law_vectors = 2 * 8 * 4 * Q * np.dtype(np.float64).itemsize
    # Both formula and omission component tuples, work blocks, FFT magnitude
    # temporaries, and converted law-list scalars can coexist at a boundary.
    component_entry_bytes = np.dtype(np.complex128).itemsize + 3 * 8
    component_scalars = 2 * PERIOD * 30 * component_entry_bytes
    work_blocks = 2 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize
    fft_vectors = 16 * Q * np.dtype(np.float64).itemsize
    guard_bytes(numeric + law_vectors + component_scalars + work_blocks
                + fft_vectors,
                "aggregate complete-law payload")
    # The omitted-earlier-phase law is independent of insertion position and
    # is evaluated once; the refined target is evaluated at all eight v.
    root_terms = 9 * PERIOD * Q * 30
    direct_blocks = 8 * Q * (PERIOD // BLOCK) * BLOCK * BLOCK
    direct_phase = 8 * Q * 2 * PERIOD
    fft_entries = 8 * Q * PERIOD
    if root_terms > MAX_ROOT_TERMS:
        raise MemoryError("geometric root-component preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": (numeric + law_vectors
                                           + component_scalars + work_blocks
                                           + fft_vectors),
        "payload_components": {"direct_fft_comparison_arrays": 5*matrix,
                                "work_vector_and_local_scalar_reserve": 32*PERIOD*16,
                                "retained_output_laws": law_vectors,
                                "component_tuple_numeric_reserve": component_scalars,
                                "work_block_reserve": work_blocks,
                                "fft_vector_temporaries": fft_vectors},
        "geometric_root_term_upper_bound": root_terms,
        "direct_block_product_preflight": direct_blocks,
        "direct_phase_query_upper_bound": direct_phase,
        "row_local_term_upper_bound": 9 * PERIOD * 2 * BLOCK * BLOCK,
        "early_phase_term_upper_bound": 9 * PERIOD * 2 * BLOCK * BLOCK * 4,
        "late_phase_query_upper_bound": 9 * PERIOD * 2,
        "fft_calls": SPLIT + 1,
        "fft_input_entries": fft_entries,
        "law_list_entries": 8 * 4 * Q,
        "component_tuple_entries": 2 * PERIOD * 30,
        "category_caps": {
            "direct_block_products": direct_blocks,
            "direct_phase_queries": direct_phase,
            "direct_shift_ops": 8 * Q * 3,
            "row_local_terms": 9 * PERIOD * 2 * BLOCK * BLOCK,
            "early_phase_terms": 9 * PERIOD * 2 * BLOCK * BLOCK * 4,
            "early_phase_queries": 9 * PERIOD * 2 * BLOCK * BLOCK * 4,
            "late_phase_queries": 9 * PERIOD * 2,
            "formula_phase_queries": 0,
            "modular_pow_queries": (direct_phase
                                     + 9 * PERIOD * 2 * BLOCK * BLOCK * 4
                                     + 9 * PERIOD * 2),
            "geometric_component_terms": root_terms,
            "fft_calls": SPLIT + 1,
            "fft_input_entries": fft_entries,
            "expansion_terms": 0,
        },
        "max_root_terms": MAX_ROOT_TERMS,
        "max_numeric_payload_bytes": MAX_BYTES,
    }


def geometric_laws(components_by_work, counters, reference_rows=None,
                   component_cap=MAX_ROOT_TERMS):
    target = np.zeros(Q, dtype=float)
    proposal = np.zeros(Q, dtype=float)
    row_norm_errors = []
    for work, components in enumerate(components_by_work):
        norm = sum(int(count) * abs(gamma)**2
                   for _start, _stride, count, gamma in components)
        if not np.isfinite(norm) or norm <= 0:
            raise ArithmeticError("invalid refined row norm")
        for y in range(Q):
            coherent = 0j
            incoherent = 0.
            for start, stride, count, gamma in components:
                if counters["geometric_component_terms"] >= component_cap:
                    raise MemoryError("geometric root term cap before phase")
                counters["geometric_component_terms"] += 1
                term = (gamma * unit_phase(-y * int(start), Q)
                        * geometric_sum(int(count), -y * int(stride), Q))
                coherent += term
                incoherent += abs(term)**2
            target[y] += abs(coherent)**2 / (Q*Q)
            proposal[y] += incoherent / (Q*Q)
        if reference_rows is None:
            row_norm_errors.append(float(norm))
        else:
            expected = float(np.sum(np.abs(reference_rows[:, work])**2))
            row_norm_errors.append(abs(norm - expected))
    return target, proposal, max(row_norm_errors)


def ensure_call_budget(counters, caps, increments, label):
    """Reject a call before entering it if its conservative work can overflow."""
    for key, amount in increments.items():
        if counters.get(key, 0) + int(amount) > caps[key]:
            raise MemoryError(
                f"{label} would exceed cumulative {key} cap: "
                f"{counters.get(key, 0)}+{amount}>{caps[key]}"
            )


def main():
    exp = Experiment("earlier_phase_output", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "direct-column FFT and refined geometric laws agree")
    exp.predict("P2", "complete laws normalize and target first bit is uniform")
    exp.predict("P3", "omission/control TVs are recorded without preclassifying them")
    exp.must_fail("C1", "component-only proposal law drops coherent cross terms")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        category_caps = report["preflight"]["category_caps"]
        W0, W1 = work_block(math.pi/4), work_block(-math.pi/10)
        all_laws = {}
        insertion_rows = {}
        counters = {"direct_block_products": 0, "direct_phase_queries": 0,
                    "direct_shift_ops": 0, "geometric_component_terms": 0,
                    "row_local_terms": 0, "formula_phase_queries": 0,
                    "early_phase_queries": 0, "late_phase_queries": 0,
                    "early_phase_terms": 0, "expansion_terms": 0,
                    "modular_pow_queries": 0, "fft_calls": 0,
                    "fft_input_entries": 0}
        max_fft_geom_error = 0.
        max_fft_geom_tv = 0.
        max_norm_error = 0.
        max_column_norm_error = 0.
        omission_tvs = {}
        proposal_tvs = {}
        first_bit_errors = []
        direct_phase_counters = {}
        # Omission of the early phase gives the same C78 law for every v.
        # Build it once at v=0 and retain its complete output vector.
        omission_rows = []
        for work in range(PERIOD):
            ensure_call_budget(
                counters, category_caps,
                {"row_local_terms": 2 * BLOCK * BLOCK,
                 "early_phase_terms": 2 * BLOCK * BLOCK,
                 "late_phase_queries": 2,
                 "modular_pow_queries": 2},
                "omitted-earlier row")
            omission_rows.append(refined_row(work, 0, W0, W1, counters,
                                              omit_early=True))
        omit_law, _, omission_component_norm_max = geometric_laws(
            omission_rows, counters,
            component_cap=category_caps["geometric_component_terms"])
        omission_law = omit_law.copy()
        for insertion in range(SPLIT + 1):
            direct_matrix = np.zeros((Q, PERIOD), dtype=np.complex128)
            formula_rows = []
            for exponent in range(Q):
                ensure_call_budget(
                    counters, category_caps,
                    {"direct_block_products": (PERIOD // BLOCK) * BLOCK * BLOCK,
                     "direct_phase_queries": 2 * PERIOD,
                     "direct_shift_ops": 3,
                     "modular_pow_queries": 2 * PERIOD},
                    "direct column")
                before = {"direct_phase_queries": counters["direct_phase_queries"],
                          "modular_pow_queries": counters["modular_pow_queries"],
                          "direct_block_products": counters["direct_block_products"],
                          "direct_shift_ops": counters["direct_shift_ops"]}
                direct_matrix[exponent] = direct_column(
                    exponent, insertion, W0, W1, counters)
                for key, value in before.items():
                    direct_phase_counters[key] = direct_phase_counters.get(key, 0) + (
                        counters[key] - value)
            ensure_call_budget(counters, category_caps,
                               {"fft_calls": 1, "fft_input_entries": Q * PERIOD},
                               "literal FFT")
            counters["fft_calls"] += 1
            counters["fft_input_entries"] += Q * PERIOD
            fft_amplitudes = np.fft.fft(direct_matrix, axis=0)
            fft_law = np.sum(np.abs(fft_amplitudes)**2, axis=1) / (Q*Q)
            column_norm_error = max(
                abs(float(np.sum(np.abs(direct_matrix[e])**2)) - 1.)
                for e in range(Q))
            max_column_norm_error = max(max_column_norm_error,
                                        column_norm_error)
            for work in range(PERIOD):
                cycle_classes = min(effective_phase_period(insertion),
                                    (L + PERIOD - 1) // PERIOD)
                early_term_bound = 2 * BLOCK * BLOCK * cycle_classes
                ensure_call_budget(
                    counters, category_caps,
                    {"row_local_terms": 2 * BLOCK * BLOCK,
                     "early_phase_terms": early_term_bound,
                     "early_phase_queries": early_term_bound,
                     "late_phase_queries": 2,
                     "modular_pow_queries": early_term_bound + 2},
                    "refined row")
                refined = refined_row(work, insertion, W0, W1, counters)
                formula_rows.append(refined)
                direct_norm = float(np.sum(np.abs(direct_matrix[:, work])**2))
                refined_norm = sum(int(n)*abs(gamma)**2
                                   for _s, _d, n, gamma in refined)
                max_norm_error = max(max_norm_error, abs(direct_norm-refined_norm))
            geom_law, proposal_law, row_norm_error = geometric_laws(
                formula_rows, counters, direct_matrix,
                category_caps["geometric_component_terms"])
            proposal_tvs[str(insertion)] = tv(geom_law, proposal_law)
            omission_tvs[str(insertion)] = tv(fft_law, omission_law)
            max_fft_geom_error = max(max_fft_geom_error,
                                     float(np.max(np.abs(fft_law-geom_law))))
            max_fft_geom_tv = max(max_fft_geom_tv, tv(fft_law, geom_law))
            first_bit_errors.append(max(abs(float(fft_law[z::2].sum())-.5)
                                        for z in range(2)))
            all_laws[str(insertion)] = {
                "fft_literal": fft_law.tolist(),
                "geometric_refined": geom_law.tolist(),
                "component_only_proposal": proposal_law.tolist(),
                "omitted_early": omit_law.tolist(),
            }
            insertion_rows[str(insertion)] = {
                "generic_P": generic_phase_period(insertion),
                "effective_P": effective_phase_period(insertion),
                "fft_geom_max_error": float(np.max(np.abs(fft_law-geom_law))),
                "fft_geom_tv": tv(fft_law, geom_law),
                "omission_tv": omission_tvs[str(insertion)],
                "component_only_tv": proposal_tvs[str(insertion)],
                "fft_mass": float(fft_law.sum()),
                "geom_mass": float(geom_law.sum()),
                "proposal_mass": float(proposal_law.sum()),
                "omission_mass": float(omission_law.sum()),
                "row_norm_error": row_norm_error,
                "direct_column_norm_error": float(column_norm_error),
                "first_bit_error": first_bit_errors[-1],
            }
            insertion_rows[str(insertion)]["omission_thresholds"] = {
                "tv_le_1e-3": bool(omission_tvs[str(insertion)] <= 1e-3),
                "tv_le_1e-2": bool(omission_tvs[str(insertion)] <= 1e-2),
            }
            insertion_rows[str(insertion)]["component_thresholds"] = {
                "tv_le_1e-3": bool(proposal_tvs[str(insertion)] <= 1e-3),
                "tv_le_1e-2": bool(proposal_tvs[str(insertion)] <= 1e-2),
            }
        all_valid = all(valid_law(all_laws[str(v)]["fft_literal"])
                        and valid_law(all_laws[str(v)]["geometric_refined"])
                        and valid_law(all_laws[str(v)]["component_only_proposal"])
                        and valid_law(all_laws[str(v)]["omitted_early"])
                        for v in range(SPLIT + 1))
        category_reconciliation = {
            key: {"actual": int(counters.get(key, 0)),
                  "cap": int(cap),
                  "within": bool(counters.get(key, 0) <= cap)}
            for key, cap in category_caps.items()
        }
        c1 = max(proposal_tvs.values()) > 1e-8
        p1 = max_fft_geom_error < TOL and max_fft_geom_tv < TOL
        p2 = (all_valid and max(first_bit_errors) < TOL
              and max_norm_error < TOL and max_column_norm_error < TOL
              and counters["fft_calls"] == report["preflight"]["fft_calls"]
              and counters["fft_input_entries"] == report["preflight"]["fft_input_entries"]
              and counters["direct_block_products"] == report["preflight"]["direct_block_product_preflight"]
              and counters["row_local_terms"] == report["preflight"]["row_local_term_upper_bound"])
        threshold_match = all(
            insertion_rows[str(v)]["omission_thresholds"] == {
                "tv_le_1e-3": bool(omission_tvs[str(v)] <= 1e-3),
                "tv_le_1e-2": bool(omission_tvs[str(v)] <= 1e-2),
            }
            and insertion_rows[str(v)]["component_thresholds"] == {
                "tv_le_1e-3": bool(proposal_tvs[str(v)] <= 1e-3),
                "tv_le_1e-2": bool(proposal_tvs[str(v)] <= 1e-2),
            }
            for v in range(SPLIT + 1)
        )
        p3 = (all(0 <= omission_tvs[str(v)] <= 1+TOL
                  and 0 <= proposal_tvs[str(v)] <= 1+TOL
                  for v in range(SPLIT + 1))
              and threshold_match
              and all(item["within"] for item in category_reconciliation.values()))
        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L},
            "laws": all_laws,
            "insertion_rows": insertion_rows,
            "omission_tvs": omission_tvs,
            "component_only_tvs": proposal_tvs,
            "category_reconciliation": category_reconciliation,
            "threshold_classifications_match_measured_tvs": threshold_match,
            "threshold_classifications": {
                "omission": {key: insertion_rows[key]["omission_thresholds"]
                              for key in insertion_rows},
                "component_only": {key: insertion_rows[key]["component_thresholds"]
                                    for key in insertion_rows},
            },
            "omission_component_norm_max": float(omission_component_norm_max),
            "max_fft_geom_error": max_fft_geom_error,
            "max_fft_geom_tv": max_fft_geom_tv,
            "max_norm_error": max_norm_error,
            "max_direct_column_norm_error": max_column_norm_error,
            "max_first_bit_error": max(first_bit_errors),
            "counters": counters,
            "direct_phase_counter_delta": direct_phase_counters,
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "tv_thresholds_recorded": [1e-3, 1e-2],
        })
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "complete FFT and geometric laws agree")
        exp.check("P2", p2, "all laws normalize and first bit is uniform")
        exp.check("P3", p3, "TV baselines are finite and recorded")
        exp.fail_check("C1", c1, "component-only proposal loses cross-component coherence")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before complete-law audit")
        exp.fail_check("C1", False, "exception before coherence control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    def jsonable(value):
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {key: jsonable(item) for key, item in value.items()}
        if isinstance(value, list):
            return [jsonable(item) for item in value]
        if isinstance(value, tuple):
            return [jsonable(item) for item in value]
        return value

    report = jsonable(report)
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                    "t": WIDTH, "s": SPLIT},
        "reference": "ER direct columns plus refined geometric components",
        "complete_probability_vectors_retained": True,
        "no_sampler": True, "no_timing_claim": True,
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_geometric_root_terms": MAX_ROOT_TERMS,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
