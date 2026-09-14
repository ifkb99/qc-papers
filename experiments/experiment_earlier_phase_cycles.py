"""Coefficient audit for one earlier pointwise phase insertion.

The frozen schedule has N=61, a=2, r=60, b=3, t=8, split s=7,
W0=work_block(pi/4), W1=work_block(-pi/10), and the late phase
g(j)=exp(2*pi*i*(2**j mod 61)/61).  A second identical pointwise phase is
inserted after v low controls, v=0,...,7, before the unchanged late mixer and
phase.  This file tests coefficients and amplitudes only; it does not claim a
full-output sampler.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The cyclic-refined column formula agrees with an independent direct
      vector-shift/diagonal/block-mixer schedule for all 256 columns at every
      insertion position.
  P2  Refined row components, expanded amplitudes, and norm sums agree with
      that same direct schedule for all 60 work labels.  The period bound is
      an upper bound, not a prediction of tightness.
  P3  v=0,1,2 have P=1; the endpoint v=7 is constant in n despite its loose
      period bound.  Interior positions are tested for actual variation.
  C1  Omitting the earlier phase from the row coefficients fails at an
      interior insertion position.

The direct schedule is a tiny independent reference using cyclic vector
shifts, diagonal pointwise phases, and b-by-b block multiplication.  The
refined formula never reads an orbit/physical-label table.  This is a bounded
float coefficient test, not a production sampler, full Fourier-law test,
timing claim, or general multi-defect result.
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
from lab.fourier_sampling import unit_phase
from experiments.experiment_clean_orbit_output import work_block


N, BASE, PERIOD, BLOCK, WIDTH, SPLIT = 61, 2, 60, 3, 8, 7
Q, L = 1 << WIDTH, 1 << SPLIT
MAX_BYTES = 16 << 20
MAX_TERMS = 1_000_000
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"earlier_phase_cycles_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def valid_unitary(block):
    return (isinstance(block, np.ndarray) and block.shape == (BLOCK, BLOCK)
            and np.all(np.isfinite(block))
            and np.allclose(block.conj().T @ block, np.eye(BLOCK), rtol=0, atol=1e-12))


def phase(index, counters, key="formula_phase_queries"):
    if counters[key] >= MAX_TERMS:
        raise MemoryError(f"{key} cap before modular power")
    counters[key] += 1
    counters["modular_pow_queries"] += 1
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def preflight():
    # Retained/replaced matrices, a previous row's view and comparison temps.
    # Numeric scalar lists get a separate conservative reserve. Python object
    # headers/container bytes are not claimed as process RSS.
    matrix_bytes = 6 * Q * PERIOD * np.dtype(np.complex128).itemsize
    vectors = 8 * PERIOD * np.dtype(np.complex128).itemsize
    setup = 4 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize
    scalar_lists = 8 * PERIOD * 30 * 16 + 128*1024
    guard_bytes(matrix_bytes + vectors + setup + scalar_lists, "aggregate cycle-audit payload")
    local_column_terms = 8 * Q * BLOCK * BLOCK
    local_row_terms = 3 * 8 * PERIOD * 2 * BLOCK * BLOCK
    direct_matrix_terms = 8 * Q * ((PERIOD // BLOCK) * BLOCK * BLOCK)
    expansion_terms = 2 * 8 * PERIOD * 2 * (2 * BLOCK - 1) * 3
    early_phase_terms = 2 * 8 * PERIOD * 2 * BLOCK * BLOCK * 3
    formula_phase_terms = (8 * Q * 2 * BLOCK * BLOCK
                           + 8 * PERIOD * 2 * 2
                           + early_phase_terms // 2)
    direct_phase_terms = 8 * Q * 2 * PERIOD
    total = (local_column_terms + local_row_terms + direct_matrix_terms
             + expansion_terms + early_phase_terms + formula_phase_terms
             + direct_phase_terms)
    if total > MAX_TERMS:
        raise MemoryError("earlier-phase scalar-work preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": matrix_bytes + vectors + setup + scalar_lists,
        "payload_components": {"direct_formula_and_comparison_temporaries": matrix_bytes,
                                "working_vectors": vectors,
                                "block_setup_temporaries": setup,
                                "numeric_scalar_lists": scalar_lists},
        "column_local_term_preflight": local_column_terms,
        "row_local_term_preflight": local_row_terms,
        "direct_matrix_term_preflight": direct_matrix_terms,
        "expansion_term_preflight": expansion_terms,
        "early_phase_term_preflight": early_phase_terms,
        "formula_phase_query_upper_bound": formula_phase_terms,
        "direct_phase_query_upper_bound": direct_phase_terms,
        "total_scalar_work_preflight": total,
        "max_scalar_work": MAX_TERMS,
        "max_numeric_payload_bytes": MAX_BYTES,
    }


def generic_phase_period(v):
    return (1 << v) // math.gcd(PERIOD, 1 << v)


def effective_phase_period(v):
    return 1 if v == SPLIT else generic_phase_period(v)


def direct_column(exponent, insertion, W0, W1, counters, *, split=SPLIT,
                  early_insertions=None):
    """Existing tiny literal schedule; optional positions share the G1 phase.

    The single-insertion default and its operation order are unchanged.
    Multiple ascending positions (including duplicates) are for bounded
    independent reference checks, not a compressed sampler.
    """
    if (not isinstance(split, (int, np.integer)) or isinstance(split, (bool, np.bool_))
            or not 0 <= insertion <= split <= 63):
        raise ValueError("require 0<=insertion<=split<=63 for literal columns")
    positions = (insertion,) if early_insertions is None else tuple(early_insertions)
    if (any(not isinstance(v, (int, np.integer)) or isinstance(v, (bool, np.bool_))
            or not 0 <= v <= split for v in positions)
            or tuple(sorted(positions)) != positions):
        raise ValueError("literal phase positions must be ascending integers within split")
    length = 1 << int(split)
    h, l = divmod(int(exponent), length)
    vector = np.zeros(PERIOD, dtype=np.complex128)
    vector[:BLOCK] = W0[:, 0]
    previous_prefix = 0
    for position in positions:
        prefix = l % (1 << int(position))
        counters["direct_shift_ops"] += 1
        vector = np.roll(vector, prefix - previous_prefix)
        for index in range(PERIOD):
            if vector[index] != 0:
                vector[index] *= phase(index, counters, "direct_phase_queries")
        previous_prefix = prefix
    remainder = l - previous_prefix
    counters["direct_shift_ops"] += 1
    vector = np.roll(vector, remainder)
    for cell in range(PERIOD // BLOCK):
        if counters["direct_block_products"] + BLOCK * BLOCK > MAX_TERMS:
            raise MemoryError("direct block-product cap before multiplication")
        counters["direct_block_products"] += BLOCK * BLOCK
        start = cell * BLOCK
        vector[start:start+BLOCK] = W1 @ vector[start:start+BLOCK]
    for index in range(PERIOD):
        if vector[index] != 0:
            vector[index] *= phase(index, counters, "direct_phase_queries")
    counters["direct_shift_ops"] += 1
    return np.roll(vector, length * h)


def refined_row(work, insertion, W0, W1, counters, omit_early=False):
    """Return (start,stride,count,gamma) cyclic-refined components."""
    components = []
    cycle_bound = effective_phase_period(insertion)
    for high in range(2):
        v = (int(work) - L * high) % PERIOD
        cell, p = divmod(v, BLOCK)
        pairs = {}
        for u in range(BLOCK):
            for q in range(BLOCK):
                if counters["row_local_terms"] >= MAX_TERMS:
                    raise MemoryError("row local-term cap before pair")
                counters["row_local_terms"] += 1
                rho = (cell * BLOCK + q - u) % PERIOD
                pairs.setdefault(rho, []).append((u, q, W1[p, q] * W0[u, 0]))
        late = phase(v, counters, "late_phase_queries")
        for rho, local_pairs in sorted(pairs.items()):
            old_count = max(0, 1 + (L - 1 - rho) // PERIOD)
            for cycle in range(min(cycle_bound, old_count)):
                count = 1 + (old_count - 1 - cycle) // cycle_bound
                if not count:
                    continue
                gamma = 0j
                for u, _q, amplitude in local_pairs:
                    if counters["early_phase_terms"] >= MAX_TERMS:
                        raise MemoryError("early-phase term cap before phase")
                    counters["early_phase_terms"] += 1
                    early = 1.+0j if omit_early else phase(
                        u + ((rho + cycle * PERIOD) % (1 << insertion)),
                        counters, "early_phase_queries")
                    gamma += amplitude * early
                gamma *= late
                if gamma != 0:
                    components.append((L * high + rho + cycle * PERIOD,
                                       PERIOD * cycle_bound, count, gamma))
    return components


def expand_row(components, counters):
    result = {}
    for start, stride, count, gamma in components:
        for n in range(int(count)):
            if counters["expansion_terms"] >= MAX_TERMS:
                raise MemoryError("row expansion cap before append")
            counters["expansion_terms"] += 1
            exponent = int(start) + int(stride) * n
            if exponent in result:
                raise AssertionError("duplicate refined exponent")
            result[exponent] = gamma
    return result


def cancellation_revival_control():
    """Tiny algebraic guard: a phase can revive a canceled candidate residue."""
    local = (1.+0j, -1.+0j)
    without_phase = sum(local)
    with_phase = local[0] * (1.+0j) + local[1] * (1j)
    return abs(without_phase) == 0. and abs(with_phase) > 1e-8


def main():
    exp = Experiment("earlier_phase_cycles", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "refined columns match the independent direct schedule")
    exp.predict("P2", "refined rows and norms match direct amplitudes")
    exp.predict("P3", "period bounds and endpoint/low-v simplifications hold")
    exp.must_fail("C1", "omitting the earlier phase changes an interior row")
    exp.must_fail("C2", "discarding canceled candidate residues misses phase revival")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = c2 = False
    try:
        if (N, BASE, PERIOD, BLOCK, WIDTH, SPLIT) != (61, 2, 60, 3, 8, 7):
            raise AssertionError("frozen earlier-phase fixture changed")
        report["preflight"] = preflight()
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        if not valid_unitary(W0) or not valid_unitary(W1):
            raise ValueError("invalid frozen work blocks")
        counters = {"column_local_terms": 0, "row_local_terms": 0, "early_phase_terms": 0,
                    "late_phase_queries": 0, "early_phase_queries": 0,
                    "direct_phase_queries": 0, "formula_phase_queries": 0,
                    "modular_pow_queries": 0, "direct_block_products": 0,
                    "direct_shift_ops": 0, "expansion_terms": 0}
        insertion_rows = {}
        column_errors, column_norm_errors, row_errors, norm_errors = [], [], [], []
        unrefined_errors = []
        endpoint_counts = {}
        stride_values = {}
        variation_witnesses = {}
        for insertion in range(SPLIT + 1):
            component_counts = []
            direct_matrix = np.zeros((Q, PERIOD), dtype=np.complex128)
            formula_matrix = np.zeros((Q, PERIOD), dtype=np.complex128)
            for exponent in range(Q):
                direct_matrix[exponent] = direct_column(
                    exponent, insertion, W0, W1, counters)
                h, l = divmod(exponent, L)
                prefix = l % (1 << insertion)
                for u in range(BLOCK):
                    source = (u + l) % PERIOD
                    cell, q = divmod(source, BLOCK)
                    early = phase(u + prefix, counters, "formula_phase_queries")
                    for p in range(BLOCK):
                        if counters["column_local_terms"] >= MAX_TERMS:
                            raise MemoryError("column local-term cap before product")
                        counters["column_local_terms"] += 1
                        before = cell * BLOCK + p
                        formula_matrix[exponent, (before + L * h) % PERIOD] += (
                            W0[u, 0] * early * W1[p, q]
                            * phase(before, counters, "formula_phase_queries"))
            column_errors.append(float(np.max(abs(direct_matrix-formula_matrix))))
            column_norm_errors.append(float(max(
                max(abs(np.sum(abs(direct_matrix[e])**2)-1.) for e in range(Q)),
                max(abs(np.sum(abs(formula_matrix[e])**2)-1.) for e in range(Q)))))
            for work in range(PERIOD):
                row_counters = counters
                refined = refined_row(work, insertion, W0, W1, row_counters)
                component_counts.append(len(refined))
                expanded = expand_row(refined, row_counters)
                omitted = expand_row(refined_row(work, insertion, W0, W1,
                                                 row_counters, omit_early=True),
                                     row_counters)
                direct_row = direct_matrix[:, work]
                row_errors.append(float(max(
                    [abs(expanded.get(e, 0j) - direct_row[e]) for e in range(Q)])))
                unrefined_errors.append(float(max(
                    [abs(omitted.get(e, 0j) - direct_row[e]) for e in range(Q)])))
                norm = sum(int(count) * abs(gamma)**2
                           for _start, _stride, count, gamma in refined)
                norm_errors.append(abs(norm - float(np.sum(abs(direct_row)**2))))
                counts = [int(count) for _start, _stride, count, _gamma in refined]
                endpoint_counts[str(insertion)] = endpoint_counts.get(str(insertion), []) + counts
                stride_values.setdefault(str(insertion), set()).update(
                    int(stride) for _start, stride, _count, _gamma in refined)
                if insertion in range(3, 7):
                    by_base = {}
                    for start, _stride, _count, gamma in refined:
                        high0, low0 = divmod(int(start), L)
                        base_key = (high0, low0 % PERIOD)
                        by_base.setdefault(base_key, []).append((start,gamma))
                    for (high0,rho0), values in by_base.items():
                        for start_a, a in values:
                            for start_b, b in values:
                                if (str(insertion) not in variation_witnesses
                                        and abs(a-b) > 1e-8):
                                    variation_witnesses[str(insertion)] = {
                                        "work": work, "high": high0, "rho": rho0,
                                        "exponents": [int(start_a),int(start_b)],
                                        "coefficients": [[float(a.real),float(a.imag)],
                                                         [float(b.real),float(b.imag)]],
                                        "absolute_difference": float(abs(a-b))}
            insertion_rows[str(insertion)] = {
                "bound_P": generic_phase_period(insertion),
                "column_max_error": column_errors[-1],
                "column_norm_error": column_norm_errors[-1],
                "row_max_error": max(row_errors[-PERIOD:]),
                "max_norm_error": max(norm_errors[-PERIOD:]),
                "max_omitted_phase_error": max(unrefined_errors[-PERIOD:]),
                "max_component_count": max(component_counts),
                "generic_P": generic_phase_period(insertion),
                "effective_P": effective_phase_period(insertion),
                "component_bound": 2 * generic_phase_period(insertion) * min(PERIOD, 2*BLOCK-1),
                "stride_values": sorted(stride_values[str(insertion)]),
                "has_multiterm_progression": any(
                    value > 1 for value in endpoint_counts[str(insertion)]),
            }
        interior_fail = max(
            max(unrefined_errors[insertion * PERIOD:(insertion + 1) * PERIOD])
            for insertion in range(3, 7))
        pre = report["preflight"]
        formula_phase_total = (counters["formula_phase_queries"]
                               + counters["early_phase_queries"]
                               + counters["late_phase_queries"])
        actual_scalar = (counters["column_local_terms"]
                         + counters["row_local_terms"]
                         + counters["direct_block_products"]
                         + counters["expansion_terms"]
                         + counters["early_phase_terms"]
                         + formula_phase_total
                         + counters["direct_phase_queries"])
        cost_ok = bool(
            counters["column_local_terms"] <= pre["column_local_term_preflight"]
            and counters["row_local_terms"] <= pre["row_local_term_preflight"]
            and counters["direct_block_products"] <= pre["direct_matrix_term_preflight"]
            and counters["expansion_terms"] <= pre["expansion_term_preflight"]
            and counters["early_phase_terms"] <= pre["early_phase_term_preflight"]
            and formula_phase_total <= pre["formula_phase_query_upper_bound"]
            and counters["direct_phase_queries"] <= pre["direct_phase_query_upper_bound"]
            and actual_scalar <= pre["total_scalar_work_preflight"]
            and counters["row_local_terms"] == 2*8*PERIOD*2*BLOCK**2
            and counters["column_local_terms"] == 8*Q*BLOCK**2
            and counters["direct_block_products"] == 8*Q*PERIOD*BLOCK
            and counters["direct_shift_ops"] == 3*8*Q
            and counters["modular_pow_queries"] == formula_phase_total+counters["direct_phase_queries"])
        c2 = cancellation_revival_control()
        p1 = bool(max(column_errors) < TOL and max(column_norm_errors) < TOL
                  and cost_ok)
        p2 = bool(max(row_errors) < TOL and max(norm_errors) < TOL)
        p3 = bool(all(insertion_rows[str(v)]["generic_P"] == 1 for v in range(3))
                  and all(insertion_rows[str(v)]["max_component_count"]
                          <= insertion_rows[str(v)]["component_bound"]
                          for v in range(SPLIT + 1))
                  and insertion_rows["7"]["effective_P"] == 1
                  and insertion_rows["7"]["stride_values"] == [PERIOD])
        c1 = bool(interior_fail > 1e-8)
        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L},
            "insertion_rows": insertion_rows,
            "max_column_error": max(column_errors),
            "max_column_norm_error": max(column_norm_errors),
            "max_row_error": max(row_errors),
            "max_norm_error": max(norm_errors),
            "max_omitted_phase_interior_error": interior_fail,
            "endpoint_progression_counts": endpoint_counts["7"],
            "variation_witnesses": variation_witnesses,
            "counters": counters,
            "actual_scalar_work": actual_scalar,
            "cost_reconciliation_pass": cost_ok,
            "cancellation_revival_control": c2,
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1, "C2": c2}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "all refined columns match direct schedule")
        exp.check("P2", p2, "all refined rows and norms match direct amplitudes")
        exp.check("P3", p3, "period and endpoint/low-v predictions hold")
        exp.fail_check("C1", c1, "omitted earlier phase fails at an interior insertion")
        exp.fail_check("C2", c2, "candidate cancellation can be revived by a phase")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before cycle audit")
        exp.fail_check("C1", False, "exception before omitted-phase control")
        exp.fail_check("C2", False, "exception before cancellation-revival control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                    "t": WIDTH, "s": SPLIT},
        "direct_reference": "vector shifts, diagonal phases, and repeated b-block W",
        "no_full_fourier_laws": True, "no_sampler": True,
        "max_numeric_payload_bytes": MAX_BYTES,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
