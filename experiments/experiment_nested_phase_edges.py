"""Small edge-only amplitude audit for TODO40 nested phase cuts.

This supersedes the earlier broad r60 pilot in this file. It uses only the
existing indexed literal-column reference and experiment-local ``cut_row``
prototype; no new
statevector or vector-shift propagator is implemented here. Old reports remain
historical superseded artifacts.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every frozen edge fixture, cut, and work row agrees in amplitudes, norm,
      and (when norm is positive) the complete normalized Fourier law.
  P2  Empty, initial, duplicate, endpoint, non-divisor-period, and r=b alias
      schedules agree with their existing literal/folded references.
  P3  The r6 Hadamard phase fixture revives the known rho=0 candidate at
      exponent 0 with coefficient 1/2+i/2; removing it changes that amplitude.
  C1  Removing that actual cut-row candidate must fail, rather than merely
      comparing two reference constructions.
  C2  In the preserved near-canceled r3 fixture, small absolute amplitude
      error does NOT imply an accurate conditional law: normalization
      magnifies roundoff. This is not a claim of appreciable joint-law error.

This is a bounded complex128 diagnostic, not a numerical certificate,
sampler, timing claim, or generic multi-phase propagator.
"""
from __future__ import annotations

import math
import platform
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.fourier_sampling import unit_phase
from experiments import experiment_earlier_phase_sampler_edges as edge_ref

MAX_BYTES = 4 << 20
MAX_TERMS = 1_000_000
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_edges_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def safe(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [safe(v) for v in value]
    return value


def reserve(ledger, key, amount, cap=MAX_TERMS):
    amount = int(amount)
    if amount < 0 or ledger.get(key, 0) + amount > cap:
        raise MemoryError(f"{key} cap before operation")
    ledger[key] = ledger.get(key, 0) + amount


def parity(index):
    index = int(index)
    if not 0 <= index < 6:
        raise AssertionError("r6 parity oracle received out-of-range index")
    return -1. + 0j if index & 1 else 1. + 0j


def f2(index):
    index = int(index)
    if not 0 <= index < 6:
        raise AssertionError("r6 f2 oracle received out-of-range index")
    return (1j) ** (index & 1)


def f3(index):
    index = int(index)
    if not 0 <= index < 3:
        raise AssertionError("r3 oracle received out-of-range index")
    return unit_phase(index, 3)


def f9(index):
    # This non-divisor-period assertion catches accidental global-60 reduction.
    index = int(index)
    if not 0 <= index < 9:
        raise AssertionError("r9 oracle received out-of-range index")
    return unit_phase(index, 9)


def identity(_index):
    return 1. + 0j


def h2():
    return np.asarray([[1., 1.], [1., -1.]], dtype=np.complex128) / math.sqrt(2.)


def f3_matrix():
    return np.asarray([[unit_phase(row * col, 3) for col in range(3)]
                       for row in range(3)], dtype=np.complex128) / math.sqrt(3.)


def preflight():
    fixtures = {"A_r6_revival": (6, 2, 2, 2),
                "B_r6_empty": (6, 2, 3, 2),
                "C_r6_initial": (6, 2, 3, 2),
                "D_r9_duplicate": (9, 3, 4, 3),
                "E_r_equals_b": (3, 3, 3, 2),
                # Reserve a full tiny family for one retained negative row.
                "F_r3_roundoff_control": (3, 3, 3, 2)}
    reference_columns = sum((1 << width) * period
                             for period, _b, width, _s in fixtures.values())
    helper_calls = reference_columns
    group_bound = coeff_bound = phase_call_bound = phase_product_bound = expansion_bound = 0
    for period, block, width, split in fixtures.values():
        q, h = 1 << width, 1 << (width - split)
        for cut in range(split + 1):
            k = (1 << split) // (1 << cut)
            p = 1 << max(0, cut - 1)
            f = min(p, ((1 << cut) + period - 1) // period)
            coef = h * k * block * block * f
            group_bound += period * h * k * block * block
            coeff_bound += period * coef
            phase_call_bound += period * (h + 4 * coef)
            phase_product_bound += period * (5 * coef)
            expansion_bound += period * h * (1 << split)
    reference_block = sum((1 << width) * (period // block) * block * block
                          for period, block, width, _s in fixtures.values())
    reference_phase = sum((1 << width) * (period + block)
                          for period, block, width, _s in fixtures.values())
    fft_calls = 2 * helper_calls
    fft_entries = sum(2 * (1 << width) * (1 << width) * period
                      for period, _b, width, _s in fixtures.values())
    comparison_visits = sum(4 * (1 << width) * (1 << width) * period
                            for period, _b, width, _s in fixtures.values())
    norm_visits = sum((4 * (1 << width) * (1 << width) + 2 * (1 << width)) * period
                      for period, _b, width, _s in fixtures.values())
    # Four folded-matrix products, six matrix-entry constructions and twelve
    # setup phase evaluations, followed by all frozen unitarity checks.
    setup_terms = (2 * 2**3 + 2 * 3**3 + 39 + 12
                   + sum(4 * block**3 for _period, block, _width, _s in fixtures.values()))
    reference_bounds = {
        "column_calls": sum(1 << width for _r, _b, width, _s in fixtures.values()),
        "early_phase_queries": sum((1 << width) * block
                                    for period, block, width, _s in fixtures.values()),
        "column_local_terms": reference_block,
        "late_phase_queries": sum((1 << width) * period
                                   for period, _b, width, _s in fixtures.values()),
        "shift_operations": sum(1 << width for _r, _b, width, _s in fixtures.values()),
    }
    control_reserve = 32 + 64 + 64 + 8 * 4
    reserved = (reference_columns + reference_block + reference_phase
                + group_bound + coeff_bound + phase_call_bound
                + phase_product_bound + expansion_bound
                + helper_calls + fft_calls + fft_entries + comparison_visits
                + norm_visits + setup_terms + control_reserve)
    if reserved > MAX_TERMS:
        raise MemoryError(f"edge preflight {reserved}>{MAX_TERMS}")
    payload = (3 * max((1 << width) * period * 16
                       for period, _b, width, _s in fixtures.values())
               + 256 * 1024)
    if payload > MAX_BYTES:
        raise MemoryError(f"edge payload {payload}>{MAX_BYTES}")
    return {"fixtures": fixtures, "reference_columns": reference_columns,
            "helper_calls": helper_calls,
            "reserved_reference_block_terms": reference_block,
            "reserved_reference_phase_terms": reference_phase,
            "reserved_group_pair_visits": group_bound,
            "reserved_coefficient_pair_visits": coeff_bound,
            "reserved_phase_call_terms": phase_call_bound,
            "reserved_phase_product_terms": phase_product_bound,
            "reserved_expansion_terms": expansion_bound,
            "reserved_fft_calls": fft_calls,
            "reserved_fft_entries": fft_entries,
            "reserved_comparison_visits": comparison_visits,
            "reserved_norm_visits": norm_visits,
            "reserved_setup_terms": setup_terms,
            "reserved_control_terms": control_reserve,
            "reference_bounds": reference_bounds,
            "reserved_total_terms": reserved,
            "planned_numeric_payload_bytes": payload,
            "max_numeric_payload_bytes": MAX_BYTES,
            "max_total_terms": MAX_TERMS}


def unpack_row(row):
    if not isinstance(row, dict) or "components" not in row:
        raise TypeError("cut_row must return a components dict")
    components = []
    for item in row["components"]:
        if len(item) != 3:
            raise ValueError("component is not (start,count,coefficient)")
        start, count, coefficient = item
        components.append((int(start), int(count), complex(coefficient)))
    return components, int(row["stride"]), float(row["norm"])


def expand(components, stride, q, actual):
    values = np.zeros(int(q), dtype=np.complex128)
    for start, count, coefficient in components:
        actual["expansion_terms"] += int(count)
        for n in range(count):
            exponent = start + stride * n
            if not 0 <= exponent < q or values[exponent] != 0:
                raise AssertionError("overlapping/out-of-range edge support")
            values[exponent] = coefficient
    return values


def law(amplitudes):
    values = np.abs(np.asarray(amplitudes)) ** 2
    total = float(values.sum())
    if not math.isfinite(total) or np.any(~np.isfinite(values)):
        raise AssertionError("nonfinite Fourier law")
    if total == 0.0:
        return None
    return values / total


def run_fixture(name, *, period, block, width, split, initial, middle,
                phase, early_phases, reference_initial, reference_middle,
                reference_early_split, reference_early_phase, cuts, cut_row,
                reservations, actual, details):
    q = 1 << width
    reserve(reservations, "reference_columns", q * period)
    reserve(reservations, "reference_block_terms", q * (period // block) * block * block)
    reserve(reservations, "reference_phase_terms", q * (period + block))
    reference = edge_ref.literal_columns(
        period, block, width, split, reference_initial, reference_middle,
        identity, reference_early_split, reference_early_phase)
    if reference.shape != (q, period):
        raise AssertionError(f"{name} literal shape mismatch")
    reserve(reservations, "norm_visits", 2 * q * period)
    actual["norm_visits"] += q * period
    column_norms = np.sum(np.abs(reference) ** 2, axis=1)
    if not np.all(np.isfinite(column_norms)) or np.max(np.abs(column_norms - 1.)) > TOL:
        raise AssertionError(f"{name} unnormalized literal columns")
    actual["norm_visits"] += q * period
    work_mass = float(np.sum(np.abs(reference) ** 2)) / q
    if abs(work_mass - 1.) > TOL:
        raise AssertionError(f"{name} work marginal does not normalize")
    max_amp = max_norm = max_tv = 0.
    for cut in cuts:
        for work in range(period):
            local = defaultdict(int)
            reserve(reservations, "cut_row_calls", 1)
            q, h = 1 << width, 1 << (width - split)
            k = (1 << split) // (1 << cut)
            left = [int(v) for v, _oracle in early_phases if int(v) < cut]
            p = ((1 << max(left)) // math.gcd(period, 1 << max(left))) if left else 1
            f = min(p, ((1 << cut) + period - 1) // period)
            candidates = h * k * block * block * f
            reserve(reservations, "reserved_group_pair_visits", h * k * block * block)
            reserve(reservations, "reserved_coefficient_pair_visits", candidates)
            reserve(reservations, "reserved_phase_call_terms",
                    h + len(early_phases) * candidates)
            reserve(reservations, "reserved_phase_product_terms",
                    (len(early_phases) + 1) * candidates + candidates)
            reserve(reservations, "reserved_expansion_terms", h * (1 << split))
            reserve(reservations, "fft_calls", 2)
            reserve(reservations, "fft_entries", 2 * q)
            reserve(reservations, "comparison_visits", 4 * q)
            reserve(reservations, "norm_visits", 4 * q)
            try:
                row = cut_row(work, period=period, block=block, width=width,
                              split=split, initial=initial, middle=middle,
                              phase=phase, early_phases=early_phases, cut=cut,
                              counters=local, max_terms=MAX_TERMS)
            finally:
                for key, value in local.items():
                    actual[key] += int(value)
            components, stride, reported_norm = unpack_row(row)
            expanded = expand(components, stride, q, actual)
            expected = reference[:, work]
            actual["comparison_visits"] += 4 * q
            actual["norm_visits"] += 2 * q
            amp = float(np.max(np.abs(expanded - expected)))
            norm = float(np.sum(np.abs(expanded) ** 2))
            expected_norm = float(np.sum(np.abs(expected) ** 2))
            norm_error = max(abs(norm - expected_norm),
                             abs(norm - reported_norm))
            actual["fft_calls"] += 2
            actual["fft_entries"] += 2 * q
            actual["norm_visits"] += 2 * q
            observed_law = law(np.fft.fft(expanded) / q)
            reference_law = law(np.fft.fft(expected) / q)
            if observed_law is None or reference_law is None:
                if not (observed_law is None and reference_law is None
                        and norm == 0. and expected_norm == 0.):
                    raise AssertionError(f"zero-law mismatch {name} cut={cut} work={work}")
                tv = 0.
            else:
                tv = float(np.sum(np.abs(observed_law - reference_law)) / 2.)
            max_amp = max(max_amp, amp)
            max_norm = max(max_norm, norm_error)
            max_tv = max(max_tv, tv)
            if amp > TOL or norm_error > TOL or tv > TOL:
                raise AssertionError(f"edge mismatch {name} cut={cut} work={work}")
    details[name] = {"cuts": list(cuts), "works": period,
                     "work_marginal_total": work_mass,
                     "max_literal_column_norm_error": float(np.max(np.abs(column_norms - 1.))),
                     "max_amplitude_error": max_amp,
                     "max_norm_error": max_norm, "max_conditional_tv": max_tv}
    return max_amp, max_norm, max_tv


def main():
    exp = Experiment("nested_phase_edges", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "all small edge rows match literal amplitudes and laws")
    exp.predict("P2", "edge schedules use their supplied period")
    exp.predict("P3", "r6 phase revives the canceled rho=0 candidate")
    exp.must_fail("C1", "removing that actual candidate changes amplitude")
    exp.must_fail("C2", "tiny absolute error need not give a stable near-zero conditional law")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = c2 = False
    reservations = defaultdict(int)
    actual = defaultdict(int)
    reference_actual = defaultdict(int)
    old_reference_ledger = edge_ref.REFERENCE_LEDGER
    old_reference_bounds = edge_ref.REFERENCE_BOUNDS
    old_setup_ledger = edge_ref.SETUP_LEDGER
    try:
        report["preflight"] = preflight()
        reference_bounds = report["preflight"]["reference_bounds"]
        edge_ref.REFERENCE_LEDGER = reference_actual
        edge_ref.REFERENCE_BOUNDS = reference_bounds
        edge_ref.SETUP_LEDGER = {"constructor_calls": 0,
                                 "unitary_check_product_terms": 0}
        reserve(reservations, "setup_terms",
                report["preflight"]["reserved_setup_terms"])
        from experiments.experiment_nested_phase_amplitudes import cut_row
        actual["setup_terms"] += 39 + 12 + 2 * 2**3 + 2 * 3**3
        H = h2()
        folded_h = H @ np.diag([1., 1j])
        w03 = f3_matrix()
        # Column-phased DFT keeps this complex fixture away from accidental
        # near-zero Fourier rows, while remaining exactly unitary.
        w13 = w03 @ np.diag([1., 1j, -1j])
        initial_c = np.diag([1., 1j]) @ H
        initial_e = np.diag([f3(j) for j in range(3)]) @ w03
        duplicate9 = lambda index: f9(index) * f9(index)
        fixtures = {
            "A_r6_revival": dict(period=6, block=2, width=2, split=2,
                initial=H, middle=H, phase=identity,
                early_phases=((1, parity), (2, f2)), reference_initial=H,
                reference_middle=folded_h, reference_early_split=1,
                reference_early_phase=parity, cuts=range(3)),
            "B_r6_empty": dict(period=6, block=2, width=3, split=2,
                initial=H, middle=H, phase=identity, early_phases=(),
                reference_initial=H, reference_middle=H,
                reference_early_split=1, reference_early_phase=identity,
                cuts=range(3)),
            "C_r6_initial": dict(period=6, block=2, width=3, split=2,
                initial=H, middle=H, phase=identity,
                early_phases=((0, f2), (1, parity)),
                reference_initial=initial_c, reference_middle=H,
                reference_early_split=1, reference_early_phase=parity,
                cuts=range(3)),
            "D_r9_duplicate": dict(period=9, block=3, width=4, split=3,
                initial=w03, middle=w13, phase=identity,
                early_phases=((1, f9), (1, f9)),
                reference_initial=w03, reference_middle=w13,
                reference_early_split=1, reference_early_phase=duplicate9,
                cuts=range(4)),
            "E_r_equals_b": dict(period=3, block=3, width=3, split=2,
                initial=w03, middle=w13, phase=identity,
                early_phases=((0, f3), (1, f3)),
                reference_initial=initial_e, reference_middle=w13,
                reference_early_split=1, reference_early_phase=f3,
                cuts=range(3)),
        }
        # Charge and perform the fixed matrix/setup checks before any row call.
        for fixture in (fixtures.values()):
            block = int(fixture["block"])
            for matrix in (fixture["initial"], fixture["middle"],
                           fixture["reference_initial"], fixture["reference_middle"]):
                if matrix.shape != (block, block) or not np.all(np.isfinite(matrix)):
                    raise AssertionError("nonfinite/malformed edge block")
                actual["setup_terms"] += block ** 3
                if not np.allclose(matrix.conj().T @ matrix,
                                   np.eye(block), atol=1e-12, rtol=0):
                    raise AssertionError("edge block failed unitary checksum")
        details = {}
        maxima = [run_fixture(name, cut_row=cut_row,
                               reservations=reservations, actual=actual,
                               details=details, **fixture)
                   for name, fixture in fixtures.items()]
        max_amp = max(value[0] for value in maxima)
        max_norm = max(value[1] for value in maxima)
        max_tv = max(value[2] for value in maxima)
        local = defaultdict(int)
        reserve(reservations, "cut_row_calls", 1)
        reserve(reservations, "reserved_group_pair_visits", 2 * 1 * 2 * 2)
        reserve(reservations, "reserved_coefficient_pair_visits", 2 * 1 * 1 * 2 * 2)
        reserve(reservations, "reserved_phase_call_terms", 2 + 2 * 4)
        reserve(reservations, "reserved_phase_product_terms", 3 * 4 + 4)
        reserve(reservations, "reserved_expansion_terms", 2 * 4)
        try:
            actual_row = cut_row(1, period=6, block=2, width=2, split=2,
                                 initial=H, middle=H, phase=identity,
                                 early_phases=((1, parity), (2, f2)), cut=2,
                                 counters=local, max_terms=MAX_TERMS)
        finally:
            for key, value in local.items():
                actual[f"helper_{key}"] += int(value)
        components, stride, _ = unpack_row(actual_row)
        actual_values = expand(components, stride, 4, actual)
        revived = actual_values[0]
        if abs(revived - (0.5 + 0.5j)) > TOL:
            raise AssertionError(f"r6 revival coefficient changed: {revived}")
        removed = [item for item in components if item[0] == 0]
        wrong_values = expand([item for item in components if item[0] != 0],
                              stride, 4, actual)
        wrong_error = abs(wrong_values[0] - revived)
        c1 = bool(removed and wrong_error > TOL)
        # Reproduce the superseded fixture as a NEGATIVE precision control,
        # not as a silently omitted positive case. Its tiny total work mass
        # makes this conditional discrepancy negligible in the joint law.
        actual["setup_terms"] += 27 + 9
        old_middle = np.diag([1., 1j, -1j]) @ w03
        reserve(reservations, "reference_columns", 8 * 3)
        reserve(reservations, "reference_block_terms", 8 * 9)
        reserve(reservations, "reference_phase_terms", 8 * (3 + 3))
        roundoff_reference = edge_ref.literal_columns(
            3, 3, 3, 2, initial_e, old_middle, identity, 1, f3)[:, 0]
        reserve(reservations, "cut_row_calls", 1)
        for key, amount in {"reserved_group_pair_visits": 72,
                            "reserved_coefficient_pair_visits": 72,
                            "reserved_phase_call_terms": 146,
                            "reserved_phase_product_terms": 224,
                            "reserved_expansion_terms": 8,
                            "fft_calls": 2, "fft_entries": 16,
                            "comparison_visits": 32, "norm_visits": 32}.items():
            reserve(reservations, key, amount)
        local = defaultdict(int)
        try:
            roundoff_row = cut_row(
                0, period=3, block=3, width=3, split=2, initial=w03,
                middle=old_middle, phase=identity, early_phases=((0,f3),(1,f3)),
                cut=0, counters=local, max_terms=MAX_TERMS)
        finally:
            for key, value in local.items():
                actual[key] += int(value)
        comp, stride, _ = unpack_row(roundoff_row)
        roundoff_values = expand(comp, stride, 8, actual)
        actual["comparison_visits"] += 32
        actual["norm_visits"] += 32
        actual["fft_calls"] += 2
        actual["fft_entries"] += 16
        roundoff = {
            "amplitude_error": float(np.max(np.abs(roundoff_values-roundoff_reference))),
            "formula_norm": float(np.sum(np.abs(roundoff_values)**2)),
            "reference_norm": float(np.sum(np.abs(roundoff_reference)**2)),
            "conditional_tv": float(np.sum(np.abs(
                law(np.fft.fft(roundoff_values))-law(np.fft.fft(roundoff_reference))))/2.),
        }
        c2 = (roundoff["amplitude_error"] < TOL
              and 0 < roundoff["formula_norm"] < 1e-25
              and 0 < roundoff["reference_norm"] < 1e-25
              and roundoff["conditional_tv"] > 1e-3)
        p1 = max_amp <= TOL and max_norm <= TOL and max_tv <= TOL
        p2 = p1
        actual_total = (sum(int(value) for value in actual.values())
                        + sum(int(value) for value in reference_actual.values()))
        helper_actual = defaultdict(int)
        for key, value in actual.items():
            helper_actual[key.removeprefix("helper_")] += int(value)
        budget_checks = {
            "setup_terms": actual.get("setup_terms", 0)
                          <= reservations["setup_terms"],
            "group_pair_visits": helper_actual["group_pair_visits"]
                                 <= reservations["reserved_group_pair_visits"],
            "coefficient_pair_visits": helper_actual["coefficient_pair_visits"]
                                       <= reservations["reserved_coefficient_pair_visits"],
            "phase_terms": (helper_actual["phase_calls"]
                             <= reservations["reserved_phase_call_terms"]
                             and helper_actual["phase_products"]
                             <= reservations["reserved_phase_product_terms"]),
            "expansion_terms": actual.get("expansion_terms", 0)
                              <= reservations["reserved_expansion_terms"],
            "fft_calls": actual.get("fft_calls", 0) <= reservations["fft_calls"],
            "fft_entries": actual.get("fft_entries", 0) <= reservations["fft_entries"],
            "comparison_visits": actual.get("comparison_visits", 0)
                                 <= reservations["comparison_visits"],
            "norm_visits": actual.get("norm_visits", 0)
                           <= reservations["norm_visits"],
        }
        budget_checks["reference"] = all(
            int(value) <= int(reference_bounds[key])
            for key, value in reference_actual.items())
        # Reconcile reservations themselves against the workload preflight;
        # actual<=reservation alone would accept an inflated reservation.
        pre = report["preflight"]
        reservation_names = {
            "reference_columns": "reference_columns",
            "reference_block_terms": "reserved_reference_block_terms",
            "reference_phase_terms": "reserved_reference_phase_terms",
            "cut_row_calls": "helper_calls",
            "setup_terms": "reserved_setup_terms",
            "fft_calls": "reserved_fft_calls", "fft_entries": "reserved_fft_entries",
            "comparison_visits": "reserved_comparison_visits",
            "norm_visits": "reserved_norm_visits",
        }
        for key in reservations:
            limit = pre[reservation_names.get(key, key)]
            # One additional control row was reserved separately before work.
            if key.startswith("reserved_") or key == "cut_row_calls":
                limit += pre["reserved_control_terms"]
            budget_checks[f"reservation_{key}"] = reservations[key] <= limit
        budget_checks["reservation_total"] = sum(reservations.values()) <= pre["reserved_total_terms"]
        budget_checks["actual_total"] = actual_total <= pre["reserved_total_terms"]
        revival_ok = abs(revived - (0.5 + 0.5j)) <= TOL
        p3 = revival_ok and actual_total <= MAX_TERMS and all(budget_checks.values())
        report.update({"fixtures": details,
            "metrics": {"max_amplitude_error": max_amp,
                        "max_norm_error": max_norm,
                        "max_conditional_tv": max_tv,
                        "revived_coefficient": safe(revived),
                        "wrong_removed_candidate_error": float(wrong_error)},
            "wrong_control": {"removed_components": safe(removed),
                              "correct_e0": safe(revived),
                              "wrong_e0": safe(wrong_values[0]),
                              "error": float(wrong_error)},
            "near_zero_roundoff_control": roundoff,
            "budget_reconciliation": budget_checks,
            "counters": safe(dict(actual)),
            "reference_counters": safe(dict(reference_actual)),
            "reservations": safe(dict(reservations)),
            "actual_total_terms": actual_total,
            "counter_semantics": "named entry/product/query ledger units, not FLOPs or timing",
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1, "C2": c2},
            "status": "PASS" if p1 and p2 and p3 and c1 and c2 else "FAIL",
            "python_version": sys.version, "numpy_version": np.__version__,
            "platform": platform.platform()})
        exp.check("P1", p1, "small edge amplitudes/norms/laws")
        exp.check("P2", p2, "period and schedule edge identities")
        exp.check("P3", p3, "bounded actual edge workload")
        exp.fail_check("C1", c1, "removed actual candidate changes e0")
        exp.fail_check("C2", c2, "near-canceled row has small amplitude error but unstable conditional law")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        report["counters_partial"] = safe(dict(actual))
        report["reference_counters_partial"] = safe(dict(reference_actual))
        report["reservations_partial"] = safe(dict(reservations))
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before small-edge audit")
        exp.fail_check("C1", False, "exception before canceled-candidate control")
        exp.fail_check("C2", False, "exception before roundoff control")
    finally:
        edge_ref.REFERENCE_LEDGER = old_reference_ledger
        edge_ref.REFERENCE_BOUNDS = old_reference_bounds
        edge_ref.SETUP_LEDGER = old_setup_ledger
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[safe(report)], metadata={
        "scope": "small edge-only nested phase amplitude audit",
        "reference": "existing edge_ref.literal_columns and experiment-local cut_row prototype",
        "no_new_propagator": True, "max_numeric_payload_bytes": MAX_BYTES,
        "max_total_terms": MAX_TERMS})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
