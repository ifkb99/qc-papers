"""Bounded work-first sparse-row Fourier audit for the C77 pilot.

This uses the frozen N=61, a=2, r=60, b=3, Q=64 fixture.  The existing
closed boundary-row helper supplies one dense diagnostic column per exponent;
the actual work-first calculation only enumerates the exponent labels allowed
by the R=4 pre-terminal support cone.  A common terminal W is stripped by
using the retained pre-final rows, which is valid because work is traced.

For a fixed work label j, let a_e=phi_e(j), Z_j=sum_e |a_e|^2, and
e=c+Lh with L=Q/H.  Summing the inverse-QFT probabilities over the unused
suffix gives

  p(z|j) = 1/(H Z_j) sum_c |sum_{e mod L=c}
             a_e exp(-2*pi*i*z*floor(e/L)/H)|^2.

The global exp(-2*pi*i*z*c/Q) cancels inside each modulus.  The complete
output law is the H=Q case, mixed with p(j)=Z_j/Q.  This is an exact finite
row identity, not a scalable column oracle or a new propagator.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The sparse-row complete law agrees with the saved independent UP law.
  P2  Every sparse-row conditional and every grouped prefix law normalizes;
      the work marginal is nonuniform or uniform, whichever the frozen rows
      actually determine.
  C1  Replacing each coherent row by its |a_e|^2 mixture makes the output
      uniform and must disagree with the nonuniform coherent law.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from experiments.experiment_uniform_prefix_statevec import (
    BASE, BLOCK, N, PERIOD, Q, WIDTH, RADIUS_BEFORE_FINAL,
    branch_vectors, orbit_labels,
)


MAX_BYTES = 16 * 1024 * 1024
MAX_PREFIX_TERMS = 100_000
MAX_ROW_QUERY_TERMS = 100_000
TOL = 3e-10
SAVED_UP_REPORT = Path("out/uniform_prefix_output_20260911T101723682629Z.json")


def report_path() -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"work_first_rows_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(payload: int, label: str) -> int:
    payload = int(payload)
    if payload < 0 or payload > MAX_BYTES:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def guard_shape(shape, dtype=np.complex128, label="array") -> int:
    return guard_bytes(math.prod(int(x) for x in shape)
                       * np.dtype(dtype).itemsize, label)


def valid_law(values, size=Q) -> bool:
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (size,) and np.all(np.isfinite(values))
                and np.all(values >= -TOL)
                and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right) -> float:
    return float(np.sum(np.abs(np.asarray(left) - np.asarray(right))) / 2)


def json_safe(value):
    """Convert numpy scalar metadata without changing measured values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def candidate_exponents(work_label: int, radius: int = RADIUS_BEFORE_FINAL):
    """Enumerate e by modular progressions, not by scanning all Q per row."""
    values = set()
    for delta in range(-int(radius), int(radius) + 1):
        residue = (int(work_label) - delta) % PERIOD
        for exponent in range(residue, Q, PERIOD):
            values.add(exponent)
    return tuple(sorted(values))


def sparse_row_prefix(row, candidates, bits: int, counts, kind):
    H, L = 1 << bits, Q >> bits
    cap = MAX_PREFIX_TERMS if kind == "prefix_terms" else MAX_ROW_QUERY_TERMS
    if counts[kind] + H*len(candidates) > cap:
        raise MemoryError("cumulative sparse-row term cap BEFORE calculation")
    groups = {}
    for exponent in candidates:
        groups.setdefault(int(exponent % L), []).append(int(exponent))
    weights = np.abs(row[list(candidates)]) ** 2
    z_norm = float(weights.sum())
    if not np.isfinite(z_norm) or z_norm <= 0:
        raise ArithmeticError("zero or invalid work-row mass")
    result = np.zeros(H, dtype=float)
    terms = 0
    for z in range(H):
        for exponents in groups.values():
            total = 0j
            for exponent in exponents:
                if counts[kind] >= cap:
                    raise MemoryError("sparse Fourier term cap")
                total += row[exponent] * np.exp(
                    -2j * np.pi * z * (exponent // L) / H)
                counts[kind] += 1
                terms += 1
            result[z] += abs(total) ** 2
        result[z] /= H * z_norm
    return result, terms, z_norm


def main():
    exp = Experiment("work_first_rows", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "sparse-row grouped Fourier law matches saved full UP law")
    exp.predict("P2", "row conditionals, prefixes, and work marginal normalize")
    exp.must_fail("C1", "incoherent Born mixture loses the nonuniform full law")
    started = time.perf_counter()
    report = {"status": "FAIL", "rows": []}
    p1 = p2 = c1 = False
    try:
        if (N, BASE, PERIOD, BLOCK, WIDTH, Q) != (61, 2, 60, 3, 6, 64):
            raise AssertionError("work-first fixture changed")
        orbit = orbit_labels()
        if len(orbit) != PERIOD or len(set(orbit)) != PERIOD:
            raise AssertionError("orbit order/labels invalid")
        max_candidates = (2 * RADIUS_BEFORE_FINAL + 1) * math.ceil(Q / PERIOD)
        if PERIOD*max_candidates > 2_000:
            raise MemoryError("candidate enumeration cap before constructing rows")
        candidates_by_work = {
            int(physical): candidate_exponents(orbit_index)
            for orbit_index, physical in enumerate(orbit)
        }
        candidate_entries = sum(len(values) for values in candidates_by_work.values())
        expected_prefix_terms = candidate_entries * sum(1 << d for d in (1, 2, 3))
        expected_row_terms = candidate_entries * Q
        expected_mixture_terms = candidate_entries * Q
        # Dense closed-row oracle plus retained work laws and all prefix outputs.
        work_dim = 1 << ((N - 1).bit_length())
        branch_setup = 3 * guard_shape((work_dim, work_dim), label="three branch setup matrices")
        branch_retained = 2 * guard_shape((Q, work_dim), label="retained final/pre-final rows")
        branch_temporaries = 4 * guard_shape((work_dim, work_dim), label="branch construction temporaries")
        retained_laws = (guard_shape((Q,), np.float64, "reference output law")
                         + guard_shape((work_dim,), np.float64, "work marginal")
                         + 4 * guard_shape((Q,), np.float64, "sparse output laws"))
        planned = branch_setup + branch_retained + branch_temporaries + retained_laws
        guard_bytes(planned, "aggregate work-first payload before arrays")
        prefix_work_upper = len(orbit) * max_candidates * sum(1 << d for d in (1, 2, 3))
        row_work_upper = len(orbit) * max_candidates * Q
        if (expected_prefix_terms > MAX_PREFIX_TERMS
                or expected_row_terms > MAX_ROW_QUERY_TERMS
                or expected_mixture_terms > MAX_ROW_QUERY_TERMS):
            raise MemoryError("sparse-row operation preflight exceeds cap")
        report.update({"fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                                    "Q": Q, "t": WIDTH},
                       "support_radius": RADIUS_BEFORE_FINAL,
                       "planned_payload_bytes": planned,
                       "payload_components": {"branch_setup": branch_setup,
                                              "branch_retained": branch_retained,
                                              "branch_temporaries": branch_temporaries,
                                              "retained_laws": retained_laws},
                       "max_candidate_entries_per_work": max_candidates,
                       "candidate_entry_upper_bound": len(orbit) * max_candidates,
                       "candidate_entry_preflight": candidate_entries,
                       "prefix_term_preflight": expected_prefix_terms,
                       "row_query_term_preflight": expected_row_terms,
                       "mixture_term_preflight": expected_mixture_terms,
                       "prefix_term_upper_bound": prefix_work_upper,
                       "row_query_term_upper_bound": row_work_upper})

        work = {"branch_columns": 0, "work_matvec_terms": 0,
                "shift_entry_updates": 0, "phase_entries": 0,
                "setup_repeated_matrices": 0, "setup_phase_values": 0,
                "qft_gate_entry_updates": 0, "qft_gates": 0}
        expected_branch_work = {"branch_columns": Q,
                                "work_matvec_terms": Q * 2 * PERIOD**2,
                                "shift_entry_updates": Q * 2 * PERIOD,
                                "phase_entries": Q * PERIOD,
                                "setup_repeated_matrices": 3,
                                "setup_phase_values": PERIOD,
                                "qft_gate_entry_updates": 0,
                                "qft_gates": 0}
        if Q > 64 or expected_branch_work["work_matvec_terms"] > 500_000:
            raise MemoryError("dense oracle work cap before calling the helper")
        _, pre_final = branch_vectors(np.asarray(orbit, dtype=np.int64),
                                      identity_last=False,
                                      retain_pre_final=True, work=work)
        pre_final = np.asarray(pre_final)
        candidate_entries = sum(len(v) for v in candidates_by_work.values())
        outside_max = 0.
        outside_nonzero = 0
        for j in orbit:
            allowed = set(candidates_by_work[int(j)])
            outside = [abs(pre_final[e, int(j)]) for e in range(Q)
                       if e not in allowed]
            if outside:
                outside_max = max(outside_max, float(max(outside)))
                outside_nonzero += sum(value > 2e-12 for value in outside)
        report["candidate_entries_actual"] = candidate_entries
        report["outside_support_max_abs"] = outside_max
        report["outside_support_nonzero_count"] = outside_nonzero

        work_marginal = np.zeros(1 << ((N - 1).bit_length()), dtype=float)
        candidate_marginal = np.zeros_like(work_marginal)
        coherent_law = np.zeros(Q, dtype=float)
        incoherent_law = np.zeros(Q, dtype=float)
        uniform_j_law = np.zeros(Q, dtype=float)
        row_prefixes = {str(d): np.zeros(1 << d, dtype=float) for d in (1, 2, 3)}
        term_counts = {"prefix_terms": 0, "row_terms": 0, "mixture_terms": 0}
        row_count = 0
        max_row_norm_error = 0.
        max_prefix_norm_error = 0.
        for j in orbit:
            j = int(j)
            allowed = candidates_by_work[j]
            row = pre_final[:, j]
            z_dense = float(np.sum(np.abs(row) ** 2))
            z_sparse = float(np.sum(np.abs(row[list(allowed)]) ** 2))
            work_marginal[j] = z_dense / Q
            candidate_marginal[j] = z_sparse / Q
            if z_sparse <= 0:
                raise ArithmeticError("work row has zero sparse mass")
            full_row, full_terms, z_used = sparse_row_prefix(
                row, allowed, WIDTH, term_counts, "row_terms")
            max_row_norm_error = max(max_row_norm_error, abs(float(full_row.sum())-1.))
            if not valid_law(full_row):
                raise ArithmeticError("invalid conditional full-row law")
            coherent_law += (z_sparse / Q) * full_row
            row_mixture = np.zeros(Q, dtype=float)
            for exponent in allowed:
                if term_counts["mixture_terms"]+Q > MAX_ROW_QUERY_TERMS:
                    raise MemoryError("Born-mixture term cap BEFORE vector update")
                row_mixture += (abs(row[exponent]) ** 2 / z_sparse) * (1. / Q)
                term_counts["mixture_terms"] += Q
            incoherent_law += (z_sparse / Q) * row_mixture
            uniform_j_law += full_row / len(orbit)
            for d in (1, 2, 3):
                prefix, terms, _ = sparse_row_prefix(
                    row, allowed, d, term_counts, "prefix_terms")
                if not valid_law(prefix, 1 << d):
                    raise ArithmeticError("invalid conditional sparse prefix")
                max_prefix_norm_error = max(max_prefix_norm_error,
                                            abs(float(prefix.sum())-1.))
                row_prefixes[str(d)] += (z_sparse / Q) * prefix
            row_count += 1
        if float(work_marginal.sum()) <= 0:
            raise ArithmeticError("empty work marginal")
        report["work_marginal"] = work_marginal.tolist()
        report["candidate_work_marginal"] = candidate_marginal.tolist()
        report["work_marginal_mass"] = float(work_marginal.sum())
        report["candidate_work_marginal_mass"] = float(candidate_marginal.sum())
        report["work_marginal_uniform_tv"] = tv(
            work_marginal[orbit], np.full(len(orbit), 1 / len(orbit)))
        report["candidate_marginal_error"] = tv(work_marginal, candidate_marginal)
        report["coherent_law"] = coherent_law.tolist()
        report["uniform_j_law"] = uniform_j_law.tolist()
        report["coherent_law_mass"] = float(coherent_law.sum())
        report["uniform_j_law_mass"] = float(uniform_j_law.sum())
        report["row_count"] = row_count
        report["max_row_norm_error"] = max_row_norm_error
        report["max_prefix_norm_error"] = max_prefix_norm_error
        report["prefix_laws"] = {key: value.tolist()
                                  for key, value in row_prefixes.items()}
        report["prefix_law_masses"] = {key: float(value.sum())
                                       for key, value in row_prefixes.items()}

        mixture_law = incoherent_law
        report["incoherent_mixture_law"] = mixture_law.tolist()
        report["incoherent_mixture_tv"] = tv(coherent_law, mixture_law)
        report["incoherent_mixture_mass"] = float(mixture_law.sum())
        report["uniform_j_tv_vs_coherent"] = tv(uniform_j_law, coherent_law)
        report["uniform_j_classification"] = (
            "different" if report["uniform_j_tv_vs_coherent"] > 1e-10
            else "invisible")

        saved = json.loads(SAVED_UP_REPORT.read_text())
        saved_law = np.asarray(saved["rows"][0]["full_output"], dtype=float)
        if saved_law.shape != (Q,):
            raise ValueError("saved UP law has wrong shape")
        report["saved_reference_report"] = str(SAVED_UP_REPORT)
        report["saved_reference_valid"] = valid_law(saved_law)
        report["saved_law"] = saved_law.tolist()
        report["sparse_vs_saved_tv"] = tv(coherent_law, saved_law)
        report["uniform64_vs_saved_tv"] = tv(mixture_law, saved_law)
        report["prefix_vs_saved_tv"] = {
            str(d): tv(row_prefixes[str(d)],
                       np.array([saved_law[z::1 << d].sum() for z in range(1 << d)]))
            for d in (1, 2, 3)}
        report["work_counts"] = {key: int(value) for key, value in work.items()}
        report["sparse_row_term_counts"] = dict(term_counts)
        report["expected_sparse_row_term_counts"] = {
            "prefix_terms": expected_prefix_terms,
            "row_terms": expected_row_terms,
            "mixture_terms": expected_mixture_terms,
        }
        report["actual_candidate_entries"] = candidate_entries
        report["actual_prefix_terms"] = term_counts["prefix_terms"]
        report["actual_row_query_terms"] = term_counts["row_terms"]
        report["actual_mixture_terms"] = term_counts["mixture_terms"]
        column_norms = np.sum(np.abs(pre_final) ** 2, axis=1)
        report["column_norm_max_error"] = float(np.max(np.abs(column_norms-1.)))

        p1 = (report["saved_reference_valid"] and valid_law(coherent_law)
              and report["candidate_marginal_error"] < TOL
              and report["outside_support_nonzero_count"] == 0
              and report["sparse_vs_saved_tv"] < 3e-10
              and report["coherent_law_mass"] > 1-TOL)
        p2 = (row_count == len(orbit) and max_row_norm_error < TOL
              and max_prefix_norm_error < TOL
              and report["column_norm_max_error"] < TOL
              and report["sparse_row_term_counts"] == report["expected_sparse_row_term_counts"]
              and work == expected_branch_work
              and all(abs(value-1.) < TOL for value in report["prefix_law_masses"].values())
              and all(valid_law(row_prefixes[str(d)], 1 << d) for d in (1, 2, 3))
              and valid_law(mixture_law)
              and abs(report["incoherent_mixture_mass"]-1.) < TOL
              and abs(report["work_marginal_mass"]-1.) < TOL
              and abs(report["candidate_work_marginal_mass"]-1.) < TOL
              and all(value < TOL for value in report["prefix_vs_saved_tv"].values()))
        c1 = report["incoherent_mixture_tv"] > 1e-4
        report["checks"] = {"P1": p1, "P2": p2, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and c1 else "FAIL"
        exp.check("P1", p1, "coherent sparse-row law matches saved UP law")
        exp.check("P2", p2, "row/prefix/work laws normalize")
        exp.fail_check("C1", c1, "Born-mixture control differs from coherent law")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK, "Q": Q},
        "support_radius": RADIUS_BEFORE_FINAL,
        "max_dense_bytes": MAX_BYTES,
        "max_prefix_terms": MAX_PREFIX_TERMS,
        "max_row_query_terms": MAX_ROW_QUERY_TERMS,
        "oracle": "existing closed branch_vectors; dense finite diagnostic only",
        "no_new_propagator": True, "no_scalable_column_oracle_claim": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
