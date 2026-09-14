"""Bounded draw-level validation of the C78 late-work progression sampler.

This is the first measurement of the opt-in ``LateWorkProgressions`` helper.
The frozen fixtures are N=61, a=2, r=60, b=3, W0=work_block(pi/4),
W1=work_block(-pi/10), and g(j)=exp(2*pi*i*2**j/61), with t=6 and t=8
and split=t-1.  The t=6 forced joint law is compared with the independent
pre-final dense UP oracle; the t=8 checks use the previously written literal
progression formula because no t=8 dense circuit oracle is retained.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  All forced t=6 joint probabilities agree with the dense pre-final
      reference, and both fixtures' columns/rows agree with the literal
      progression formulas.
  P2  Forced conditional/proposal/accepted laws normalize; acceptance is at
      most one and the exact mean number of proposals is the component count.
  P3  Bounded seeded draws exercise both work selection and rejection; their
      returned ranges and counters reconcile with each draw's reported path.
  C1  Replacing coherent row probabilities by the incoherent proposal law
      changes the joint output law.
  C2  Redrawing the work label after each rejected output changes the law.

This is a finite float smoke test, not a statistical distribution certificate,
timing benchmark, certified arbitrary-width sampler, or general propagator.
The sampled work label is kept fixed during retries by the intended sampler.
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
from experiments import experiment_uniform_prefix_statevec as up
from experiments.experiment_work_first_progressions import (
    row_components as formula_row_components,
    sparse_column as formula_sparse_column,
    progression_laws as formula_progression_laws,
)


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
T_VALUES = (6, 8)
MAX_BYTES = 16 << 20
MAX_TERMS = 1_000_000
MAX_SAMPLE_PROPOSALS = 2 * 512 * 256
TOL = 3e-10
SAVED_UP = Path("out/uniform_prefix_output_20260911T101723682629Z.json")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"work_first_sampler_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def valid_law(values, size):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (int(size,),) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right):
    return float(np.sum(np.abs(np.asarray(left) - np.asarray(right))) / 2)


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


def phase(index):
    return unit_phase(pow(BASE, int(index), N), N)


def phase_formula_counters():
    return {"column_local_terms": 0, "row_local_terms": 0,
            "progression_expansion_terms": 0,
            "fourier_component_terms": 0, "direct_root_terms": 0,
            "modular_pow_queries": 0, "label_pow_queries": 0}


def compare_components(actual, expected):
    actual = list(actual)
    expected = list(expected)
    if len(actual) != len(expected):
        return float("inf")
    error = 0.
    for lhs, rhs in zip(actual, expected):
        error = max(error, abs(int(lhs[0]) - int(rhs[0])))
        error = max(error, abs(int(lhs[1]) - int(rhs[1])))
        error = max(error, abs(complex(lhs[2]) - complex(rhs[2])))
    return float(error)


def add_counters(total, counters):
    for key, value in counters.items():
        if isinstance(value, (int, np.integer)):
            total[key] = total.get(key, 0) + int(value)


def preflight():
    """Guard numeric payload and scalar work before any fixture/oracle call."""
    work_dim = 1 << ((N - 1).bit_length())
    dense_branches = 2 * 64 * work_dim * np.dtype(np.complex128).itemsize
    dense_laws = PERIOD * 64 * np.dtype(np.float64).itemsize
    dense_setup = 3 * PERIOD * PERIOD * np.dtype(np.complex128).itemsize
    dense_temporaries = 4 * PERIOD * PERIOD * np.dtype(np.complex128).itemsize
    fixture_matrices = len(T_VALUES) * 2 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize
    identity_matrix = BLOCK * BLOCK * np.dtype(np.complex128).itemsize
    law_vectors = sum((2 * (1 << width) + PERIOD) * np.dtype(np.float64).itemsize
                      for width in T_VALUES)
    target_tables = sum(PERIOD * (1 << width) * np.dtype(np.float64).itemsize
                        for width in T_VALUES)
    scalar_reserve = 64 * (64 + 256) * 16
    total_payload = (dense_branches + dense_laws + dense_setup
                     + dense_temporaries + fixture_matrices
                     + identity_matrix + law_vectors + target_tables
                     + scalar_reserve)
    guard_bytes(total_payload, "aggregate sampler numeric payload")
    dense_fourier_terms = PERIOD * 64 * 64
    forced_local_terms = sum(PERIOD * (1 << width) * 2 * BLOCK * BLOCK
                             for width in T_VALUES)
    forced_phase_terms = sum(PERIOD * (1 << width) * (1 << (width-(width-1)))
                             for width in T_VALUES)
    formula_local_terms = sum((1 << width) * BLOCK * BLOCK
                              + PERIOD * 2 * BLOCK * BLOCK
                              for width in T_VALUES) + PERIOD*2*BLOCK*BLOCK
    formula_fourier_terms = PERIOD * 64 * 10
    api_diagnostic_local_terms = sum((1 << width) * BLOCK * BLOCK
                                     + PERIOD * 2 * BLOCK * BLOCK
                                     for width in T_VALUES)
    api_diagnostic_phase_terms = sum((1 << width) * BLOCK * BLOCK
                                     + PERIOD * 2 * 2
                                     for width in T_VALUES)
    api_forced_fourier_terms = PERIOD * sum(1 << width for width in T_VALUES) * 10
    sample_cap = MAX_SAMPLE_PROPOSALS
    total_scalar_work = (dense_fourier_terms + forced_local_terms
                         + forced_phase_terms + formula_local_terms
                         + formula_fourier_terms + api_diagnostic_local_terms
                         + api_diagnostic_phase_terms + api_forced_fourier_terms)
    if total_scalar_work > MAX_TERMS:
        raise MemoryError("sampler scalar-work preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": total_payload,
        "payload_components": {
            "dense_branches_pre_and_final": dense_branches,
            "dense_joint_law": dense_laws,
            "dense_setup_matrices": dense_setup,
            "dense_construction_temporaries": dense_temporaries,
            "fixture_matrices": fixture_matrices,
            "identity_control_matrix": identity_matrix,
            "law_vectors": law_vectors,
            "target_tables": target_tables,
            "scalar_reserve": scalar_reserve},
        "dense_fourier_term_preflight": dense_fourier_terms,
        "forced_local_term_preflight": forced_local_terms,
        "forced_phase_term_preflight": forced_phase_terms,
        "formula_local_term_preflight": formula_local_terms,
        "formula_fourier_term_preflight": formula_fourier_terms,
        "api_diagnostic_local_term_preflight": api_diagnostic_local_terms,
        "api_diagnostic_phase_term_upper_bound": api_diagnostic_phase_terms,
        "api_forced_fourier_term_upper_bound": api_forced_fourier_terms,
        "sample_proposal_cap_preflight": sample_cap,
        "total_scalar_work_preflight": total_scalar_work,
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_scalar_work": MAX_TERMS,
    }


def dense_joint_reference():
    """Independent t=6 joint law from the existing retained pre-final oracle."""
    orbit = up.orbit_labels()
    if len(orbit) != PERIOD or len(set(orbit)) != PERIOD:
        raise AssertionError("invalid frozen orbit")
    work_dim = 1 << ((N - 1).bit_length())
    guard_bytes(64 * work_dim * 16, "dense t=6 oracle")
    work = {"branch_columns": 0, "work_matvec_terms": 0,
            "shift_entry_updates": 0, "phase_entries": 0,
            "setup_repeated_matrices": 0, "setup_phase_values": 0,
            "qft_gate_entry_updates": 0, "qft_gates": 0,
            "dense_fourier_terms": 0}
    if 64 * 2 * PERIOD * PERIOD > 500_000:
        raise MemoryError("dense oracle work cap")
    _, pre = up.branch_vectors(np.asarray(orbit, dtype=np.int64),
                               identity_last=False, retain_pre_final=True,
                               work=work)
    pre = np.asarray(pre)
    law = np.zeros((PERIOD, 64), dtype=float)
    work_probability = np.zeros(PERIOD, dtype=float)
    for j, physical in enumerate(orbit):
        row = pre[:, int(physical)]
        z = float(np.sum(np.abs(row) ** 2))
        work_probability[j] = z / 64.
        for y in range(64):
            amp = 0j
            for e in range(64):
                if work["dense_fourier_terms"] >= MAX_TERMS:
                    raise MemoryError("dense Fourier term cap before phase")
                work["dense_fourier_terms"] += 1
                amp += row[e] * unit_phase(-y * e, 64)
            law[j, y] = abs(amp) ** 2 / (64. * 64.)
    return orbit, law, work_probability, work


def fixture(width, W0, W1):
    from lab.work_first import LateWorkProgressions
    return LateWorkProgressions(PERIOD, BLOCK, width, width - 1, W0, W1,
                                phase, max_local_terms=MAX_TERMS,
                                max_components=4096,
                                max_payload_bytes=MAX_BYTES)


def main():
    exp = Experiment("work_first_sampler", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "forced laws and column/row formulas agree with the frozen references")
    exp.predict("P2", "all conditional laws normalize and rejection has mean m")
    exp.predict("P3", "bounded seeded draws exercise the intended fixed-work sampler")
    exp.must_fail("C1", "dephasing each conditional row changes the coherent joint law")
    exp.must_fail("C2", "redrawing work after rejection changes the output law")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = c2 = False
    try:
        if (N, BASE, PERIOD, BLOCK) != (61, 2, 60, 3):
            raise AssertionError("frozen fixture changed")
        report["preflight"] = preflight()
        # All dense/reference arrays are guarded before construction.
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        guard_bytes(2 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize,
                    "work-block matrices")
        orbit, dense_t6, dense_work, oracle_work = dense_joint_reference()
        fixtures = {}
        rows = {}
        forced_errors = []
        column_errors = []
        row_errors = []
        formula_term_totals = {"column_local_terms": 0, "row_local_terms": 0,
                               "fourier_component_terms": 0}
        formula_counters = phase_formula_counters()
        wrong_redraw = np.zeros(64, dtype=float)
        target_tables = {}
        for width in T_VALUES:
            Q = 1 << width
            c = fixture(width, W0, W1)
            fixtures[str(width)] = c
            row_data = []
            api_counter_totals = {}
            coherent_joint = np.zeros(Q, dtype=float)
            incoherent_joint = np.zeros(Q, dtype=float)
            work_mass = np.zeros(PERIOD, dtype=float)
            guard_bytes(PERIOD * Q * np.dtype(np.float64).itemsize,
                        f"width-{width} target conditional table")
            target_table = np.zeros((PERIOD, Q), dtype=float)
            target_tables[str(width)] = target_table
            forced_row_checks = []
            for e in range(Q):
                actual = c.column(e)
                expected = formula_sparse_column(e, W0, W1,
                                                 1 << (width - 1), formula_counters,
                                                 MAX_TERMS)
                got = {int(k): complex(v) for k, v in actual["amplitudes"].items()}
                add_counters(api_counter_totals, actual.get("counters", {}))
                column_errors.append(max(abs(got.get(k, 0j) - v)
                                        for k, v in expected.items()))
                column_errors.append(max((abs(got.get(k, 0j)) for k in set(got)-set(expected)),
                                         default=0.))
            for j in range(PERIOD):
                actual = c.row(j)
                expected = formula_row_components(j, W0, W1,
                                                  1 << (width - 1), formula_counters,
                                                  MAX_TERMS)
                got_components = actual["components"]
                add_counters(api_counter_totals, actual.get("counters", {}))
                row_errors.append(compare_components(got_components, expected))
                norm = float(actual["norm"])
                if not np.isfinite(norm) or norm <= 0:
                    raise ArithmeticError("invalid row norm")
                row_data.append({"work": j, "component_count": len(got_components),
                                 "norm": norm,
                                 "components": json_safe(got_components),
                                 "counters": json_safe(actual.get("counters", {}))})
                work_mass[j] = norm / Q
                conditional_values = []
                proposal_values = []
                accepted_values = []
                acceptance_values = []
                for y in range(Q):
                    forced = c.forced_joint(j, y)
                    add_counters(api_counter_totals, forced.get("counters", {}))
                    conditional_values.append(float(forced["conditional_probability"] or 0.))
                    target_table[j, y] = conditional_values[-1]
                    proposal_values.append(float(forced["proposal_probability"]))
                    acceptance_values.append(float(forced["acceptance"]))
                    accepted_values.append(float(forced["proposal_probability"])
                                           * float(forced["acceptance"]))
                    if width == 6:
                        expected_joint = float(dense_t6[j, y])
                        forced_errors.append(abs(float(forced["joint_probability"])
                                                 - expected_joint))
                        wrong_redraw[y] += (float(forced["work_probability"])
                                            * float(forced["proposal_probability"])
                                            * float(forced["acceptance"]))
                    coherent_joint[y] += float(forced["joint_probability"])
                    incoherent_joint[y] += (float(forced["work_probability"])
                                            * float(forced["proposal_probability"]))
                forced_row_checks.append({
                    "work": j,
                    "component_count": len(got_components),
                    "conditional_mass": float(sum(conditional_values)),
                    "proposal_mass": float(sum(proposal_values)),
                    "accepted_mass": float(sum(accepted_values)),
                    "accepted_normalized_law_error": float(max(
                        abs(len(got_components) * accepted_values[index]
                            - conditional_values[index])
                        for index in range(Q))),
                    "max_acceptance": max(acceptance_values),
                    "mean_attempts": float(1./sum(accepted_values))
                    if sum(accepted_values) > 0 else float("inf"),
                    "conditional_min": min(conditional_values),
                    "proposal_min": min(proposal_values),
                })
            rows[str(width)] = {"row_data": row_data,
                                "coherent_joint": coherent_joint.tolist(),
                                "incoherent_joint": incoherent_joint.tolist(),
                                "work_mass": work_mass.tolist(),
                                "coherent_mass": float(coherent_joint.sum()),
                                "incoherent_mass": float(incoherent_joint.sum()),
                                "max_forced_error": max(forced_errors) if width == 6 else None,
                                "max_column_error": max(column_errors),
                                "max_row_error": max(row_errors),
                                "forced_row_checks": forced_row_checks,
                                "api_counter_totals": api_counter_totals,
                                "stats": json_safe(c.stats())}

        # Exact forced-law controls from t=6: proposal is the incoherent
        # Born mixture, while API joint_probability is the coherent target.
        primary = rows["6"]
        target = np.asarray(primary["coherent_joint"])
        incoherent = np.asarray(primary["incoherent_joint"])
        # The loop above added work_probability*proposal_probability, but the
        # intended joint target is already accumulated from forced_joint.
        p_work = np.asarray(primary["work_mass"])
        dense_work_tv = tv(dense_work, p_work)
        dense_work_max_error = float(np.max(np.abs(dense_work - p_work)))
        wrong_mass = float(wrong_redraw.sum())
        if wrong_mass > 0:
            wrong_redraw /= wrong_mass
        # Recompute coherent/incoherent laws from the literal row formula for
        # an independent check and to avoid treating the API as its own oracle.
        literal_coherent = np.zeros(64, dtype=float)
        literal_incoherent = np.zeros(64, dtype=float)
        for j in range(PERIOD):
            components = formula_row_components(j, W0, W1,
                                                32, formula_counters, MAX_TERMS)
            z = sum(int(n) * abs(gamma) ** 2 for _, n, gamma in components)
            coherent, proposal, _, _ = formula_progression_laws(
                components, z, 64, formula_counters, MAX_TERMS)
            literal_coherent += (z / 64.) * coherent
            literal_incoherent += (z / 64.) * proposal
        rows["6"]["literal_coherent_tv"] = tv(target, literal_coherent)
        rows["6"]["literal_incoherent_tv"] = tv(incoherent, literal_incoherent)
        rows["6"]["wrong_redraw_mass"] = wrong_mass
        rows["6"]["wrong_redraw_tv"] = tv(target, wrong_redraw)

        all_forced_rows = [row for width in T_VALUES
                           for row in rows[str(width)]["forced_row_checks"]]
        forced_rows_valid = all(
            np.isfinite(row["conditional_mass"])
            and np.isfinite(row["proposal_mass"])
            and np.isfinite(row["accepted_mass"])
            and np.isfinite(row["mean_attempts"])
            and row["conditional_min"] >= -TOL
            and row["proposal_min"] >= -TOL
            and abs(row["conditional_mass"] - 1.) < TOL
            and abs(row["proposal_mass"] - 1.) < TOL
            and abs(row["accepted_mass"] - 1./row["component_count"]) < TOL
            and abs(row["mean_attempts"] - row["component_count"]) < TOL
            and row["accepted_normalized_law_error"] < TOL
            for row in all_forced_rows)
        report["forced_row_checks_valid"] = forced_rows_valid
        report["max_forced_accepted_mass_error"] = max(
            abs(row["accepted_mass"] - 1./row["component_count"])
            for row in all_forced_rows)
        report["max_forced_mean_attempt_error"] = max(
            abs(row["mean_attempts"] - row["component_count"])
            for row in all_forced_rows)
        report["max_forced_accepted_normalized_law_error"] = max(
            row["accepted_normalized_law_error"] for row in all_forced_rows)
        report["max_forced_acceptance"] = max(
            row["max_acceptance"] for row in all_forced_rows)

        samples = {}
        actual_sample_proposals = 0
        for width, seed in ((6, 62406), (8, 62408)):
            Q = 1 << width
            draws = []
            sample_counter_totals = {}
            counter_checks = []
            positive_target_draws = 0
            c = fixtures[str(width)]
            for offset in range(512):
                if actual_sample_proposals+256 > MAX_SAMPLE_PROPOSALS:
                    raise MemoryError("sample-proposal reserve before bounded draw")
                result = c.sample(np.random.default_rng(seed + offset), max_attempts=256)
                add_counters(sample_counter_totals, result.get("counters", {}))
                actual_sample_proposals += int(result["counters"]["progression_proposals"])
                counters = result["counters"]
                attempts = int(result["attempts"])
                expected_marginal_queries = 2 * (width - 2) * attempts
                m = int(result["component_count"])
                checks = {
                    "column_local_terms": counters.get("column_local_terms") == BLOCK*BLOCK,
                    "row_local_terms": counters.get("row_local_terms") == 2*BLOCK*BLOCK,
                    "phase_queries": (0 < counters.get("phase_queries", 0)
                                      <= BLOCK*BLOCK + 2),
                    "work_draws": counters.get("work_draws") == 1,
                    "component_draws": counters.get("component_draws") == attempts,
                    "progression_proposals": counters.get("progression_proposals") == attempts,
                    "acceptance_draws": counters.get("acceptance_draws") == attempts,
                    "fourier_component_terms": counters.get("fourier_component_terms") == m*attempts,
                    "marginal_queries": counters.get("marginal_queries") == expected_marginal_queries,
                    "ranges": (0 <= int(result["work"]) < PERIOD
                               and 0 <= int(result["output"]) < Q
                               and 1 <= attempts <= 256),
                }
                counter_checks.append(all(checks.values()))
                positive_target_draws += int(
                    target_tables[str(width)][int(result["work"]), int(result["output"])] > 0)
                draws.append({"work": int(result["work"]),
                              "output": int(result["output"]),
                              "attempts": int(result["attempts"]),
                              "component_count": int(result["component_count"])})
            samples[str(width)] = {
                "count": len(draws),
                "work_distinct": len({d["work"] for d in draws}),
                "output_distinct": len({d["output"] for d in draws}),
                "min_attempts": min(d["attempts"] for d in draws),
                "max_attempts": max(d["attempts"] for d in draws),
                "mean_attempts": float(np.mean([d["attempts"] for d in draws])),
                "component_counts": sorted({d["component_count"] for d in draws}),
                "counter_totals": sample_counter_totals,
                "all_counter_checks": all(counter_checks),
                "positive_target_draws": positive_target_draws,
                "stats": json_safe(c.stats()),
            }

        invalid = {}
        try:
            fixtures["6"].sample(np.random.default_rng(1), max_attempts=0)
            invalid["zero_attempt_cap_raises"] = False
        except Exception:
            invalid["zero_attempt_cap_raises"] = True
        identity = None
        try:
            from lab.work_first import LateWorkProgressions
            identity = LateWorkProgressions(60, 3, 4, 3, np.eye(3), np.eye(3),
                                            lambda _: 1.+0j,
                                            max_local_terms=MAX_TERMS,
                                            max_components=4096,
                                            max_payload_bytes=MAX_BYTES)
            zero = identity.forced_joint(30, 0)
            invalid["zero_row_returns_zero"] = bool(
                zero["conditional_probability"] is None
                and zero["proposal_probability"] == 0.
                and zero["acceptance"] == 0.
                and zero["joint_probability"] == 0.)
        except Exception:
            invalid["zero_row_returns_zero"] = False

        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "widths": list(T_VALUES), "split": "width-1"},
            "rows": rows, "samples": samples, "invalid_input_controls": invalid,
            "forced_t6_max_error": max(forced_errors),
            "dense_work_tv": dense_work_tv,
            "dense_work_max_error": dense_work_max_error,
            "max_column_error": max(column_errors),
            "max_row_error": max(row_errors),
            "formula_term_totals": dict(formula_counters),
            "actual_sample_proposals": actual_sample_proposals,
            "dense_oracle_work": json_safe(oracle_work),
            "python_version": sys.version, "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        # Reconcile actual instrumented diagnostics against the same named
        # preflight categories. Dense matvec/reference setup and RNG proposals
        # are separate budgets, not hidden inside this scalar-term sum.
        diagnostic_totals = {}
        for width in T_VALUES:
            add_counters(diagnostic_totals, rows[str(width)]["api_counter_totals"])
        actual_scalar_terms = (oracle_work["dense_fourier_terms"]
            + sum(diagnostic_totals.get(k,0) for k in
                  ("column_local_terms","row_local_terms","phase_queries","fourier_component_terms"))
            + sum(formula_counters[k] for k in
                  ("column_local_terms","row_local_terms","fourier_component_terms")))
        report["actual_diagnostic_scalar_terms"] = actual_scalar_terms
        report["diagnostic_cost_reconciliation"] = (
            actual_scalar_terms <= report["preflight"]["total_scalar_work_preflight"]
            and oracle_work["dense_fourier_terms"] == PERIOD*64*64
            and diagnostic_totals["row_local_terms"] == PERIOD*(64+256+2)*18)
        p1 = bool(report["forced_t6_max_error"] < TOL
                   and report["max_column_error"] < TOL
                   and report["max_row_error"] < TOL
                   and report["dense_work_max_error"] < TOL
                   and rows["6"]["literal_coherent_tv"] < TOL
                   and valid_law(np.asarray(rows["6"]["coherent_joint"]), 64)
                   and valid_law(np.asarray(rows["8"]["coherent_joint"]), 256))
        p2 = bool(valid_law(target, 64)
                   and report["diagnostic_cost_reconciliation"]
                   and valid_law(incoherent, 64)
                   and rows["6"]["literal_incoherent_tv"] < TOL
                   and rows["6"]["wrong_redraw_mass"] > 0
                   and report["forced_row_checks_valid"]
                   and report["invalid_input_controls"]["zero_attempt_cap_raises"]
                   and report["invalid_input_controls"]["zero_row_returns_zero"])
        p3 = bool(report["max_forced_acceptance"] <= 1+TOL
                   and all(samples[str(w)]["all_counter_checks"] for w in T_VALUES)
                   and report["actual_sample_proposals"] <= MAX_SAMPLE_PROPOSALS
                   and all(samples[str(w)]["positive_target_draws"] == 512 for w in T_VALUES)
                   and all(samples[str(w)]["work_distinct"] > 1 for w in T_VALUES)
                   and all(samples[str(w)]["output_distinct"] > 1 for w in T_VALUES)
                   and all(samples[str(w)]["max_attempts"] > 1 for w in T_VALUES))
        c1 = bool(tv(target, incoherent) > 1e-8)
        c2 = bool(report["rows"]["6"]["wrong_redraw_tv"] > 1e-8)
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1, "C2": c2}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "forced laws and literal formula references agree")
        exp.check("P2", p2, "conditional/proposal laws and input gates are valid")
        exp.check("P3", p3, "bounded draws exercise work and rejection")
        exp.fail_check("C1", c1, "coherence cannot be deleted")
        exp.fail_check("C2", c2, "work cannot be redrawn after rejection")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before sampler audit")
        for name in ("C1", "C2"):
            exp.fail_check(name, False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                    "widths": list(T_VALUES)},
        "reference": "dense t6 pre-final UP plus literal progression formulas",
        "draws_are_smoke_only": True, "no_timing_claim": True,
        "max_numeric_payload_bytes": MAX_BYTES,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
