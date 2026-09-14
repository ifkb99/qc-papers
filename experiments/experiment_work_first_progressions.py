"""Closed sparse-column/progression audit for the frozen late-mixer pilot.

The pre-final work column has the literal schedule

    U^(32h) D_g W_1^rep U^l W_0|0>,   e=l+32h,

with N=61, a=2, r=60, b=3, t=6, W_0=Rx/work_block(pi/4),
W_1=work_block(-pi/10), and g(j)=exp(2*pi*i*(2^j mod 61)/61).  The common final
W(pi/11) is stripped.  This experiment evaluates C78's column and row
progression formulas using only 3-by-3 local terms, then compares them with
the existing dense pre-final branch matrix as a charged tiny oracle.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  All 64 literal sparse columns and all 60 progression-row expansions
      match the existing pre-final branch oracle, including row norms.
  P2  Coherent progression Fourier laws, component proposals and accepted
      laws normalize; their mixture matches both frozen full-law references.
  P3  Every component envelope is <=1 and each row's exact mean proposal
      attempts equals its component count m.  Component counts stay <=10.
      Follow-up derived after the first primary run: the separately frozen
      width-eight, penultimate-control-relative insertion has genuine two-
      and three-term progressions. Their geometric laws agree with direct
      finite-root sums, and the component proposal retains mean m attempts.
  C1  Omitting coherent cross-component interference changes the full law, and
      the initial work-label marginal is not uniform.

No Q- or r-sized array is created by the specialized column/row helpers.
The dense branch matrix is retained only as the explicitly capped reference;
there are no random samples or timing sweep. The secondary row-only check
varies width under the explicitly fixed relative-insertion rule; it is not
an independently compiled physical-circuit check.
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
from experiments import experiment_uniform_prefix_statevec as up


N, BASE, PERIOD, BLOCK, WIDTH = 61, 2, 60, 3, 6
Q, LATE_H, L = 1 << WIDTH, 2, 1 << 5
SECONDARY_WIDTH, SECONDARY_Q, SECONDARY_L = 8, 1 << 8, 1 << 7
MAX_BYTES = 16 << 20
MAX_TERMS = 1_000_000
TOL = 3e-10
SAVED_UP = Path("out/uniform_prefix_output_20260911T101723682629Z.json")
SAVED_WORK_FIRST = Path("out/work_first_rows_20260911T103811510460Z.json")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"work_first_progressions_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def guard_shape(shape, dtype=np.complex128, label="array"):
    return guard_bytes(math.prod(int(x) for x in shape)
                       * np.dtype(dtype).itemsize, label)


def g_phase(orbit_index, counters):
    """Pointwise phase g(j)=exp(2*pi*i*a^j/N), without an orbit table."""
    if counters["modular_pow_queries"] >= 4_000:
        raise MemoryError("pointwise modular-phase query cap")
    counters["modular_pow_queries"] += 1
    return unit_phase(pow(BASE, int(orbit_index), N), N)


def sparse_column(exponent, W0, W1, late_L, counters, term_cap):
    """Return orbit-index -> amplitude terms, without an r-sized array."""
    l, h = int(exponent) % late_L, int(exponent) // late_L
    terms = {}
    for u in range(BLOCK):
        source = (u + l) % PERIOD
        cell, q = divmod(source, BLOCK)
        for p in range(BLOCK):
            if counters["column_local_terms"] >= term_cap:
                raise MemoryError("column local-term cap before multiplication")
            counters["column_local_terms"] += 1
            amp = W0[u, 0] * W1[p, q]
            if amp == 0:
                continue
            before_d = BLOCK * cell + p
            amp *= g_phase(before_d, counters)
            output = (before_d + late_L * h) % PERIOD
            terms[output] = terms.get(output, 0j) + amp
    return terms


def row_components(work_index, W0, W1, late_L, counters, term_cap):
    """Return (start,count,gamma) progression terms, with no Q-sized row."""
    result = []
    for h in range(LATE_H):
        v = (int(work_index) - late_L * h) % PERIOD
        cell, p = divmod(v, BLOCK)
        by_rho = {}
        for u in range(BLOCK):
            for q in range(BLOCK):
                if counters["row_local_terms"] >= term_cap:
                    raise MemoryError("row local-term cap before multiplication")
                counters["row_local_terms"] += 1
                rho = (BLOCK * cell + q - u) % PERIOD
                by_rho[rho] = by_rho.get(rho, 0j) + W1[p, q] * W0[u, 0]
        gv = g_phase(v, counters)
        for rho, local in sorted(by_rho.items()):
            gamma = gv * local
            count = max(0, 1 + (late_L - 1 - rho) // PERIOD)
            if count and gamma != 0:
                result.append((int(late_L) * h + rho, count, gamma))
    if len(result) > LATE_H * min(PERIOD, 2 * BLOCK - 1):
        raise AssertionError("progression component count exceeded C78 bound")
    return result


def expand_components(components, counters, term_cap):
    """Reference-only expansion for comparisons; helper remains progression-only."""
    expanded = {}
    if counters["progression_expansion_terms"]+sum(int(n) for _,n,_ in components) > term_cap:
        raise MemoryError("progression expansion cap BEFORE row construction")
    for start, count, gamma in components:
        for n in range(int(count)):
            if counters["progression_expansion_terms"] >= term_cap:
                raise MemoryError("progression expansion cap before append")
            counters["progression_expansion_terms"] += 1
            exponent = int(start) + n * PERIOD
            if exponent in expanded:
                raise AssertionError("duplicate exponent in progression expansion")
            expanded[exponent] = gamma
    return expanded


def progression_laws(components, z_norm, q_size, counters, term_cap):
    """Return coherent, incoherent proposal and accepted laws for one row."""
    m = len(components)
    if m <= 0 or z_norm <= 0:
        raise ArithmeticError("zero progression row")
    guard_bytes(4 * int(q_size) * np.dtype(np.float64).itemsize,
                "progression law arrays")
    proposal = np.zeros(q_size, dtype=float)
    coherent = np.zeros(q_size, dtype=float)
    for y in range(q_size):
        values = []
        for start, count, gamma in components:
            if counters["fourier_component_terms"] >= term_cap:
                raise MemoryError("Fourier component cap before evaluation")
            counters["fourier_component_terms"] += 1
            first = unit_phase(-y * int(start), q_size)
            values.append(gamma * first
                          * geometric_sum(count, -y * PERIOD, q_size))
        coherent[y] = abs(sum(values)) ** 2 / (q_size * z_norm)
        proposal[y] = sum(abs(value) ** 2 for value in values) / (q_size * z_norm)
    acceptance = np.zeros(q_size, dtype=float)
    for y in range(q_size):
        if proposal[y] > 0:
            acceptance[y] = coherent[y] / (m * proposal[y])
        else:
            if abs(coherent[y]) > 3e-12:
                raise ArithmeticError("zero proposal denominator with nonzero coherent mass")
            acceptance[y] = 0.
    accepted = proposal * acceptance
    return coherent, proposal, accepted, acceptance


def direct_row_law(expanded, z_norm, q_size, counters, term_cap):
    """Reference-only finite-root sum over an already expanded secondary row."""
    guard_bytes(int(q_size) * np.dtype(np.float64).itemsize,
                "direct finite-root law")
    result = np.zeros(q_size, dtype=float)
    for y in range(q_size):
        total = 0j
        for exponent, amplitude in expanded.items():
            if counters["direct_root_terms"] >= term_cap:
                raise MemoryError("direct finite-root term cap before summation")
            counters["direct_root_terms"] += 1
            total += amplitude * unit_phase(-y * int(exponent), q_size)
        result[y] = abs(total) ** 2 / (q_size * z_norm)
    return result


def preflight():
    # Aggregate numerical payload, including reference-only secondary columns.
    # Python dictionaries/report objects are separately bounded in count, not
    # claimed as byte-exact heap/RSS measurements.
    work_dim = 1 << ((N - 1).bit_length())
    branch = guard_shape((Q, work_dim), label="dense pre-final branches")
    pre_final = branch
    matrices = 3 * guard_shape((PERIOD, PERIOD), label="three repeated setup matrices")
    temporaries = 4 * guard_shape((PERIOD, PERIOD), label="branch setup temporaries")
    laws = 5 * guard_shape((Q,), np.float64, label="reference and law arrays")
    rows = PERIOD * 10 * 3 * 16
    secondary_columns = SECONDARY_Q * BLOCK**2 * (16+8)
    secondary_laws = 8 * SECONDARY_Q * 8
    scalar_reserve = 64 * Q * 16
    total = (branch + pre_final + matrices + temporaries + laws + rows
             + secondary_columns + secondary_laws + scalar_reserve)
    guard_bytes(total, "aggregate progression audit payload")
    column_terms = Q * BLOCK * BLOCK
    row_terms = PERIOD * LATE_H * BLOCK * BLOCK
    formula_terms = PERIOD * Q * LATE_H * min(PERIOD, 2 * BLOCK - 1)
    secondary_column_terms = SECONDARY_Q * BLOCK * BLOCK
    secondary_row_terms = PERIOD * LATE_H * BLOCK * BLOCK
    secondary_components = PERIOD * LATE_H * min(PERIOD, 2 * BLOCK - 1)
    secondary_formula_terms = secondary_components * SECONDARY_Q
    secondary_direct_terms = (PERIOD * LATE_H * min(PERIOD, 2 * BLOCK - 1)
                               * math.ceil(SECONDARY_L / PERIOD) * SECONDARY_Q)
    expansion_upper = PERIOD*LATE_H*min(PERIOD,2*BLOCK-1)*(1+math.ceil(SECONDARY_L/PERIOD))
    total_terms = (column_terms + row_terms + formula_terms + expansion_upper
                   + secondary_column_terms + secondary_row_terms
                   + secondary_formula_terms + secondary_direct_terms)
    if total_terms > MAX_TERMS:
        raise MemoryError("specialized progression term preflight exceeds cap")
    return {
        "planned_payload_bytes": total,
        "payload_components": {"branch": branch, "pre_final": pre_final,
                                "matrices": matrices, "temporaries": temporaries,
                                "laws": laws, "row_component_records": rows,
                                "secondary_column_numeric_entries":secondary_columns,
                                "secondary_law_arrays":secondary_laws,
                                "other_numeric_temporaries":scalar_reserve},
        "column_local_terms": column_terms,
        "row_local_terms": row_terms,
        "fourier_component_term_upper_bound": formula_terms,
        "secondary_column_local_terms": secondary_column_terms,
        "secondary_row_local_terms": secondary_row_terms,
        "secondary_fourier_component_term_upper_bound": secondary_formula_terms,
        "secondary_direct_root_term_upper_bound": secondary_direct_terms,
        "progression_expansion_upper_bound": expansion_upper,
        "total_specialized_term_upper_bound": total_terms,
        "max_terms": MAX_TERMS,
    }


def valid_law(values, expected_size=Q):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (int(expected_size),) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right):
    return float(np.sum(abs(np.asarray(left) - np.asarray(right))) / 2)


def load_saved(path, key):
    payload = json.loads(path.read_text())
    row = payload["rows"][0]
    return np.asarray(row[key], dtype=float)


def secondary_width8_audit(W0, W1, counters):
    """Row-only L=128 control; all sparse columns are a CHARGED reference."""
    columns = [sparse_column(e, W0, W1, SECONDARY_L, counters, MAX_TERMS)
               for e in range(SECONDARY_Q)]
    row_errors, fourier_errors = [], []
    coherent_mass_errors, direct_mass_errors = [], []
    coherent_mins, direct_mins = [], []
    component_counts, progression_counts = [], []
    accepted_mass_errors, acceptance_maxima, accepted_law_errors = [], [], []
    secondary_column_norm_error = max(abs(sum(abs(a)**2 for a in col.values())-1)
                                      for col in columns)
    for j in range(PERIOD):
        components = row_components(j, W0, W1, SECONDARY_L,
                                    counters, MAX_TERMS)
        expanded = expand_components(components, counters, MAX_TERMS)
        component_counts.append(len(components))
        progression_counts.extend(int(count) for _, count, _ in components)
        row_errors.append(max(abs(expanded.get(e, 0j)
                                  - columns[e].get(j, 0j))
                              for e in range(SECONDARY_Q)))
        z_norm = sum(int(count) * abs(gamma)**2
                     for _, count, gamma in components)
        coherent, proposal, accepted, acceptance = progression_laws(
            components, z_norm, SECONDARY_Q, counters, MAX_TERMS)
        if not valid_law(proposal,SECONDARY_Q) or not valid_law(coherent,SECONDARY_Q):
            raise ArithmeticError("invalid secondary proposal/coherent law")
        accepted_mass_errors.append(abs(float(accepted.sum())-1/len(components)))
        accepted_law_errors.append(float(np.max(abs(len(components)*accepted-coherent))))
        acceptance_maxima.append(float(acceptance.max()))
        direct = direct_row_law(expanded, z_norm, SECONDARY_Q,
                                counters, MAX_TERMS)
        coherent_mass = float(coherent.sum())
        direct_mass = float(direct.sum())
        if (not np.all(np.isfinite(coherent)) or not np.all(np.isfinite(direct))
                or float(np.min(coherent)) < -TOL
                or float(np.min(direct)) < -TOL):
            raise AssertionError("secondary law has nonfinite/negative mass")
        fourier_errors.append(float(np.max(abs(coherent-direct))))
        coherent_mass_errors.append(abs(coherent_mass - 1.0))
        direct_mass_errors.append(abs(direct_mass - 1.0))
        coherent_mins.append(float(np.min(coherent)))
        direct_mins.append(float(np.min(direct)))
    return {
        "column_row_max_error": float(max(row_errors)),
        "column_norm_max_error":float(secondary_column_norm_error),
        "geometric_vs_direct_max_error": max(fourier_errors),
        "coherent_mass_max_error": max(coherent_mass_errors),
        "direct_mass_max_error": max(direct_mass_errors),
        "coherent_min": min(coherent_mins),
        "direct_min": min(direct_mins),
        "max_component_count": max(component_counts),
        "component_count_histogram": {
            str(value): component_counts.count(value)
            for value in sorted(set(component_counts))},
        "progression_count_values": sorted(set(progression_counts)),
        "all_rows_have_long_progressions": all(value > 1 for value in progression_counts),
        "has_two_and_three_term_progressions": {2, 3} <= set(progression_counts),
        "max_acceptance":max(acceptance_maxima),
        "max_accepted_mass_error":max(accepted_mass_errors),
        "max_accepted_law_error":max(accepted_law_errors),
        "columns_generated": len(columns), "rows_checked": PERIOD,
    }


def main():
    exp = Experiment("work_first_progressions", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "literal sparse columns and row progressions match the charged pre-final oracle")
    exp.predict("P2", "coherent progression mixture matches both saved full-law references")
    exp.predict("P3", "progression proposals have envelope <=1 and exact mean attempts m")
    exp.must_fail("C1", "dropping coherent cross-component interference and assuming uniform work rows fails")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        orbit_list = up.orbit_labels()
        if len(orbit_list) != PERIOD or len(set(orbit_list)) != PERIOD:
            raise AssertionError("invalid frozen orbit")
        orbit = np.asarray(orbit_list, dtype=np.int64)
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        # Existing literal branch formula is the charged dense reference oracle.
        work = {"orbit": orbit, "shift_entry_updates": 0,
                "work_matvec_terms": 0, "branch_columns": 0, "phase_entries": 0,
                "setup_repeated_matrices": 0, "setup_phase_values": 0,
                "qft_gate_entry_updates": 0, "qft_gates": 0}
        if Q > 64 or 2*Q*PERIOD**2 > 500_000:
            raise MemoryError("dense reference column-work cap before helper call")
        branches, oracle_pre = up.branch_vectors(
            orbit, identity_last=False, retain_pre_final=True, work=work)
        oracle_pre = np.asarray(oracle_pre)
        orbit_set = set(int(value) for value in orbit_list)
        column_errors = []
        column_norm_errors = []
        column_components = []
        specialized_counts = {"column_local_terms": 0,
                              "row_local_terms": 0,
                              "progression_expansion_terms": 0,
                              "fourier_component_terms": 0,
                              "direct_root_terms": 0,
                              "modular_pow_queries": 0,
                              "label_pow_queries": 0}
        for e in range(Q):
            terms = sparse_column(e, W0, W1, L, specialized_counts, MAX_TERMS)
            column_components.append(len(terms))
            err = 0.
            for index, physical in enumerate(orbit_list):
                err = max(err, abs(terms.get(index, 0j)
                               - oracle_pre[e, int(physical)]))
            outside = [abs(oracle_pre[e, j]) for j in range(oracle_pre.shape[1])
                       if j not in orbit_set]
            if outside:
                err = max(err, max(outside))
            column_errors.append(float(err))
            column_norm_errors.append(abs(sum(abs(value)**2 for value in terms.values())-1.))

        row_records = []
        guard_bytes(Q * np.dtype(np.float64).itemsize,
                    "coherent joint law")
        coherent_joint = np.zeros(Q, dtype=float)
        guard_bytes(Q * np.dtype(np.float64).itemsize,
                    "incoherent joint law")
        incoherent_joint = np.zeros(Q, dtype=float)
        guard_bytes(PERIOD * np.dtype(np.float64).itemsize,
                    "work marginal")
        work_marginal = np.zeros(PERIOD, dtype=float)
        guard_bytes(Q * np.dtype(np.float64).itemsize,
                    "uniform-work output law")
        uniform_work_joint = np.zeros(Q, dtype=float)
        max_row_error = max_norm_error = max_envelope = 0.
        max_component_count = 0
        accepted_mass_errors = []
        coherent_control_tvs = []
        duplicate_count = 0
        for j in range(PERIOD):
            specialized_counts["label_pow_queries"] += 1
            physical = pow(BASE, j, N)
            components = row_components(j, W0, W1, L, specialized_counts, MAX_TERMS)
            expanded = expand_components(components, specialized_counts, MAX_TERMS)
            if len(expanded) != sum(int(count) for _, count, _ in components):
                duplicate_count += 1
            dense_row = oracle_pre[:, physical]
            max_row_error = max(max_row_error, max(
                [abs(expanded.get(e, 0j) - dense_row[e]) for e in range(Q)]))
            z_formula = sum(int(count) * abs(gamma)**2
                            for _, count, gamma in components)
            z_dense = float(np.sum(abs(dense_row)**2))
            max_norm_error = max(max_norm_error, abs(z_formula-z_dense))
            work_marginal[j] = z_formula / Q
            coherent, proposal, accepted, acceptance = progression_laws(
                components, z_formula, Q, specialized_counts, MAX_TERMS)
            coherent_joint += work_marginal[j] * coherent
            incoherent_joint += work_marginal[j] * proposal
            uniform_work_joint += coherent / PERIOD
            accepted_mass = float(accepted.sum())
            accepted_mass_errors.append(abs(accepted_mass - 1/len(components)))
            accepted_law_error = float(np.max(abs(
                accepted * len(components) - coherent)))
            max_envelope = max(max_envelope, float(acceptance.max()))
            coherent_control_tvs.append(tv(coherent, proposal))
            max_component_count = max(max_component_count, len(components))
            row_records.append({
                "work_index": j, "physical_label": physical,
                "component_count": len(components), "Z_dense": z_dense,
                "Z_formula": z_formula, "max_row_amplitude_error": float(
                    max(abs(expanded.get(e, 0j)-dense_row[e]) for e in range(Q))),
                "components": [{"start": int(start), "count": int(count),
                                "gamma": [float(gamma.real), float(gamma.imag)]}
                               for start, count, gamma in components],
                "proposal_mass": float(proposal.sum()),
                "coherent_mass": float(coherent.sum()),
                "proposal_min": float(np.min(proposal)),
                "coherent_min": float(np.min(coherent)),
                "accepted_mass": accepted_mass,
                "accepted_normalized_mass": float(accepted.sum() * len(components)),
                "accepted_normalized_law_error": accepted_law_error,
                "mean_attempts": float(1 / accepted_mass),
                "max_acceptance": float(acceptance.max()),
                "coherent_vs_incoherent_tv": coherent_control_tvs[-1],
            })

        secondary = secondary_width8_audit(W0, W1, specialized_counts)
        primary_components = sum(row["component_count"] for row in row_records)
        primary_expanded = sum(c["count"] for row in row_records for c in row["components"])
        secondary_components = sum(int(m)*n for m,n in secondary["component_count_histogram"].items())
        expected_work = {"shift_entry_updates":2*Q*PERIOD,
            "work_matvec_terms":2*Q*PERIOD**2,"branch_columns":Q,
            "phase_entries":Q*PERIOD,"setup_repeated_matrices":3,
            "setup_phase_values":PERIOD,"qft_gate_entry_updates":0,"qft_gates":0}
        work_actual = {key:value for key,value in work.items() if key!="orbit"}
        cost_ok = (work_actual == expected_work
            and specialized_counts["column_local_terms"] == (Q+SECONDARY_Q)*BLOCK**2
            and specialized_counts["row_local_terms"] == 2*PERIOD*LATE_H*BLOCK**2
            and specialized_counts["fourier_component_terms"] == Q*primary_components+SECONDARY_Q*secondary_components
            and specialized_counts["direct_root_terms"] == SECONDARY_Q*(specialized_counts["progression_expansion_terms"]-primary_expanded)
            and specialized_counts["label_pow_queries"] == PERIOD
            and specialized_counts["modular_pow_queries"] <= (Q+SECONDARY_Q)*BLOCK**2+2*PERIOD*LATE_H
            and sum(specialized_counts[k] for k in ("column_local_terms","row_local_terms",
                "fourier_component_terms","direct_root_terms","progression_expansion_terms"))
                <= report["preflight"]["total_specialized_term_upper_bound"])
        saved_up = load_saved(SAVED_UP, "full_output")
        saved_work = load_saved(SAVED_WORK_FIRST, "coherent_law")
        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "L": L, "H": LATE_H,
                        "stripped_final_angle": math.pi / 11},
            "sparse_column_errors": column_errors,
            "column_component_counts": column_components,
            "row_records": row_records,
            "max_column_error": max(column_errors),
            "max_row_component_error": max_row_error,
            "max_row_norm_error": max_norm_error,
            "max_component_count": max_component_count,
            "duplicate_exponent_count": duplicate_count,
            "work_marginal": work_marginal.tolist(),
            "work_marginal_mass": float(work_marginal.sum()),
            "work_marginal_uniform_tv": tv(work_marginal,
                                             np.full(PERIOD, 1 / PERIOD)),
            "uniform_work_output_law": uniform_work_joint.tolist(),
            "uniform_work_output_mass": float(uniform_work_joint.sum()),
            "uniform_work_output_tv": tv(uniform_work_joint, coherent_joint),
            "coherent_joint_law": coherent_joint.tolist(),
            "incoherent_joint_law": incoherent_joint.tolist(),
            "coherent_joint_mass": float(coherent_joint.sum()),
            "incoherent_joint_mass": float(incoherent_joint.sum()),
            "incoherent_vs_coherent_tv": tv(incoherent_joint, coherent_joint),
            "coherent_vs_saved_up_tv": tv(coherent_joint, saved_up),
            "coherent_vs_saved_work_first_tv": tv(coherent_joint, saved_work),
            "saved_up_report": str(SAVED_UP),
            "saved_work_first_report": str(SAVED_WORK_FIRST),
            "max_acceptance_envelope": max_envelope,
            "max_accepted_mass_error": max(accepted_mass_errors),
            "max_accepted_normalized_law_error": max(
                row["accepted_normalized_law_error"] for row in row_records),
            "max_mean_attempt_error": max(abs(row["mean_attempts"]
                                                - row["component_count"])
                                           for row in row_records),
            "max_coherent_incoherent_row_tv": max(coherent_control_tvs),
            "max_column_norm_error": max(column_norm_errors),
            "secondary_width8": secondary,
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "oracle_work_counts": {k: (v if isinstance(v, int) else v)
                                    for k, v in work.items() if k != "orbit"},
            "oracle_branch_norm_max_error": float(np.max(
                abs(np.sum(abs(branches)**2, axis=1)-1))),
            "specialized_term_counts": specialized_counts,
            "cost_reconciliation_pass":cost_ok,
        })
        p1 = bool(report["max_column_error"] < TOL
                   and report["max_column_norm_error"] < TOL
                   and report["max_row_component_error"] < TOL
                   and report["max_row_norm_error"] < TOL
                   and report["duplicate_exponent_count"] == 0
                   and report["max_component_count"] <= 10
                   and cost_ok
                   and valid_law(work_marginal,PERIOD)
                   and report["oracle_branch_norm_max_error"] < TOL)
        row_laws_valid = all(
            np.isfinite(row["proposal_mass"])
            and np.isfinite(row["coherent_mass"])
            and row["proposal_min"] >= -TOL
            and row["coherent_min"] >= -TOL
            and row["proposal_mass"] >= -TOL
            and row["coherent_mass"] >= -TOL
            and abs(row["proposal_mass"] - 1.0) < TOL
            and abs(row["coherent_mass"] - 1.0) < TOL
            for row in row_records)
        p2 = bool(valid_law(coherent_joint) and valid_law(incoherent_joint)
                   and valid_law(uniform_work_joint)
                   and row_laws_valid
                   and report["coherent_joint_mass"] > 1-TOL
                   and report["coherent_vs_saved_up_tv"] < TOL
                   and report["coherent_vs_saved_work_first_tv"] < TOL
                   and report["incoherent_joint_mass"] > 1-TOL)
        p3 = bool(report["max_acceptance_envelope"] <= 1+TOL
                   and report["max_accepted_mass_error"] < TOL
                   and report["max_accepted_normalized_law_error"] < TOL
                   and report["max_mean_attempt_error"] < TOL
                   and all(row["proposal_mass"] > 1-TOL for row in row_records)
                   and all(row["accepted_mass"] > 0 for row in row_records)
                   and report["secondary_width8"]["column_row_max_error"] < TOL
                   and report["secondary_width8"]["column_norm_max_error"] < TOL
                   and report["secondary_width8"]["geometric_vs_direct_max_error"] < TOL
                   and report["secondary_width8"]["coherent_mass_max_error"] < TOL
                   and report["secondary_width8"]["direct_mass_max_error"] < TOL
                   and report["secondary_width8"]["all_rows_have_long_progressions"]
                   and report["secondary_width8"]["has_two_and_three_term_progressions"])
        p3 = bool(p3 and secondary["max_component_count"]<=10
                  and secondary["max_acceptance"]<=1+TOL
                  and secondary["max_accepted_mass_error"]<TOL
                  and secondary["max_accepted_law_error"]<TOL)
        c1 = bool(report["incoherent_vs_coherent_tv"] > 1e-8
                   and report["work_marginal_uniform_tv"] > 1e-8
                   and report["uniform_work_output_tv"] > 1e-8
                   and report["max_coherent_incoherent_row_tv"] > 1e-8)
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "columns, progression rows, norms, and component caps match oracle")
        exp.check("P2", p2, "coherent accepted mixture matches frozen full laws")
        exp.check("P3", p3, "proposal envelope and exact mean attempts hold for all rows")
        exp.fail_check("C1", c1, "incoherent cross-component law and uniform-work shortcut fail")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before progression comparison")
        exp.check("P2", False, "exception before full-law mixture")
        exp.check("P3", False, "exception before envelope audit")
        exp.fail_check("C1", False, "exception before negative control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK, "t": WIDTH},
        "reference": "existing uniform_prefix_statevec branch_vectors pre-final matrix",
        "saved_references": [str(SAVED_UP), str(SAVED_WORK_FIRST)],
        "specialized_helpers_no_q_or_r_arrays": True,
        "dense_reference_only": True, "max_dense_bytes": MAX_BYTES,
        "no_random_samples": True, "no_wide_sweep": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
