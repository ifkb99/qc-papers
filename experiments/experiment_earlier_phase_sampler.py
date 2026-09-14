"""Bounded sampler audit for the opt-in one-earlier-phase progression helper.

The frozen fixture is the C79/ER N=61, a=2, r=60, b=3, t=8, s=7 schedule
with W0=work_block(pi/4), W1=work_block(-pi/10), late phase g, and a second
identical g inserted after v=0,...,7 low controls.  The independent target is
assembled from ER's literal direct columns and an FFT; no generic propagator
or orbit table is constructed here.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  EarlierPhaseProgressions(auto) agrees with the independent complete
      joint work/output laws for all eight insertion positions.  Explicit
      cycle/dual modes agree on selected public rows.
  P2  Cached conditional/proposal/accepted laws normalize, row norms and
      columns agree, and the exact accepted conditional law is recovered after
      multiplying proposal*acceptance by its component count.
  P3  The genuine seeded sampler returns in-range outputs with one work draw,
      bounded attempts, and reconciled counters; these samples are not a
      histogram accuracy certificate.
  C1  Omitting coherent acceptance (using the component-only proposal law)
      changes the complete output law.

The exact v=0 folded-W0 C78 baseline and archived omission/feedback-aware
baselines are reported with fixed source paths.  This is a bounded float
validation, not a numerical certificate, timing benchmark, or production
sampler integration.
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
from lab.fourier_sampling import unit_phase
from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_earlier_phase_cycles import direct_column


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT, Q, L, H = 8, 7, 256, 128, 2
INSERTIONS = tuple(range(SPLIT + 1))
MAX_BYTES = 16 << 20
MAX_ROOT_TERMS = 5_000_000
MAX_SAMPLE_ATTEMPTS = len(INSERTIONS) * 32 * 512
MAX_SAMPLE_FOURIER_TERMS = MAX_SAMPLE_ATTEMPTS * 30
TOL = 3e-10

# These are frozen evidence inputs, not dynamically selected "latest" reports.
OUTPUT_REFERENCE_REPORT = Path("out/earlier_phase_output_20260911T114719371212Z.json")
FEEDBACK_REFERENCE_REPORT = Path("out/earlier_phase_feedback_20260911T115103487134Z.json")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"earlier_phase_sampler_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


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


def tv(left, right):
    return float(np.sum(np.abs(np.asarray(left) - np.asarray(right))) / 2.)


def valid_law(values, shape):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == tuple(shape) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) < TOL)


def cover_bounds(insertion, cover):
    """Pure integer conservative bounds matching EarlierPhaseProgressions."""
    A = 1 << int(insertion)
    P = A // math.gcd(PERIOD, A)
    effective = 1 if insertion == SPLIT else P
    K = L // A
    residue_count = min(PERIOD, 2 * BLOCK - 1)
    occupied = H * residue_count * ((L + PERIOD - 1) // PERIOD)
    if cover == "cycle":
        components = min(H * min(L, residue_count * effective), occupied, Q)
        row_pairs = H * BLOCK * BLOCK * (1 + min(effective,
                                                  (L + PERIOD - 1) // PERIOD))
    elif cover == "dual":
        components = min(H * K * min(A, residue_count), occupied, Q)
        row_pairs = H * K * BLOCK * BLOCK
    else:
        raise ValueError("cover bound requires cycle or dual")
    return {"components": components, "row_pairs": row_pairs,
            "phase_queries": row_pairs + H + 2,
            "P": P, "effective_P": effective, "K": K}


def phase(index):
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def preflight():
    """Guard all retained numeric objects and named work before construction."""
    if not OUTPUT_REFERENCE_REPORT.exists() or not FEEDBACK_REFERENCE_REPORT.exists():
        raise FileNotFoundError("frozen archived output/feedback reports missing")
    archive_file_bytes = (OUTPUT_REFERENCE_REPORT.stat().st_size
                          + FEEDBACK_REFERENCE_REPORT.stat().st_size)
    # json.loads creates Python containers larger than the source bytes. This
    # is a conservative parsed-source reserve, not an RSS claim.
    archived_source_reserve = 8 * archive_file_bytes
    matrix = Q * PERIOD * np.dtype(np.complex128).itemsize
    direct_and_helper_matrices = 6 * matrix
    # Direct, auto-target, and auto-proposal arrays are all retained.
    retained_joint_laws = (3 * len(INSERTIONS) * PERIOD * Q
                           * np.dtype(np.float64).itemsize)
    converted_joint_law_lists = retained_joint_laws
    retained_output_laws = 4 * len(INSERTIONS) * Q * np.dtype(np.float64).itemsize
    fft_vectors = 4 * Q * np.dtype(np.float64).itemsize
    work_blocks = 4 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize
    # Auto rows and one C78 baseline can coexist while their component tuples
    # are converted to the retained report form.
    tuple_bytes = 2 * PERIOD * 30 * (16 + 3 * 8)
    scalar_lists = 512 * 64 * 16
    source_baseline_vectors = 2 * len(INSERTIONS) * Q * np.dtype(np.float64).itemsize
    payload = (direct_and_helper_matrices + retained_joint_laws
               + converted_joint_law_lists
               + retained_output_laws + fft_vectors + work_blocks
               + tuple_bytes + scalar_lists + source_baseline_vectors
               + archived_source_reserve)
    guard_bytes(payload, "aggregate sampler numeric payload")

    direct_blocks = len(INSERTIONS) * Q * (PERIOD // BLOCK) * BLOCK * BLOCK
    direct_phase = len(INSERTIONS) * Q * 2 * PERIOD
    direct_shifts = len(INSERTIONS) * Q * 3
    direct_fft_entries = len(INSERTIONS) * Q * PERIOD
    helper_column_local = len(INSERTIONS) * Q * BLOCK * BLOCK
    helper_column_phase = len(INSERTIONS) * Q * (BLOCK + PERIOD)
    cover_rows = {mode: [cover_bounds(v, mode) for v in INSERTIONS]
                  for mode in ("cycle", "dual")}
    # The production auto mode chooses the cover by row-pair work, while the
    # diagnostic below evaluates both explicit covers independently.
    auto_cover_names = [
        "cycle" if cover_rows["cycle"][v]["row_pairs"]
        <= cover_rows["dual"][v]["row_pairs"] else "dual"
        for v in INSERTIONS]
    auto_rows = [cover_rows[auto_cover_names[v]][v]["row_pairs"]
                 for v in INSERTIONS]
    auto_row_local = sum(PERIOD * value for value in auto_rows)
    diagnostic_row_local = sum(PERIOD * cover_rows[mode][v]["row_pairs"]
                               for mode in ("cycle", "dual")
                               for v in INSERTIONS)
    forced_row_local = sum(3 * 4 * auto_rows[v] for v in INSERTIONS)
    helper_row_local = auto_row_local + diagnostic_row_local + forced_row_local
    auto_row_phase = sum(PERIOD * cover_rows[auto_cover_names[v]][v][
        "phase_queries"] for v in INSERTIONS)
    diagnostic_row_phase = sum(PERIOD * cover_rows[mode][v]["phase_queries"]
                               for mode in ("cycle", "dual")
                               for v in INSERTIONS)
    forced_row_phase = sum(3 * 4 * cover_rows[auto_cover_names[v]][v][
        "phase_queries"] for v in INSERTIONS)
    helper_row_phase = auto_row_phase + diagnostic_row_phase + forced_row_phase
    helper_phase = helper_column_phase + helper_row_phase
    # Full Fourier evaluation is done only for auto; cycle/dual are selected
    # Full Fourier evaluation is done only for auto; cycle/dual are selected
    # public-row diagnostics. These are exact structural component bounds.
    auto_fourier_terms = sum(
        PERIOD * Q * cover_rows[auto_cover_names[v]][v]["components"]
        for v in range(len(INSERTIONS)))
    diagnostic_fourier_terms = sum(
        PERIOD * 4 * (cover_rows["cycle"][v]["components"]
                      + cover_rows["dual"][v]["components"])
        for v in INSERTIONS)
    forced_fourier_terms = sum(3 * 4 * cover_rows[auto_cover_names[v]][v][
        "components"] for v in INSERTIONS)
    baseline_fourier_terms = PERIOD * Q * (H * min(PERIOD, 2 * BLOCK - 1))
    public_forced_calls = len(INSERTIONS) * 3 * 4
    setup_helper_constructions = (len(INSERTIONS) * 3 + 1) * 2 * BLOCK**3
    setup_folded_phase_queries = BLOCK
    sample_calls = len(INSERTIONS) * 32
    sample_column_local = sample_calls * BLOCK * BLOCK
    sample_row_local = sample_calls * max(auto_rows)
    sample_phase_queries = sample_calls * (BLOCK + PERIOD
                                           + max(cover_rows["cycle"][v]["phase_queries"]
                                                 for v in INSERTIONS))
    total = (direct_blocks + direct_phase + direct_shifts
             + helper_column_local + helper_row_local + helper_phase
             + auto_fourier_terms + diagnostic_fourier_terms
             + forced_fourier_terms + baseline_fourier_terms
             + public_forced_calls + MAX_SAMPLE_ATTEMPTS
             + MAX_SAMPLE_FOURIER_TERMS + sample_column_local
             + sample_row_local + sample_phase_queries
             + 2 * WIDTH * MAX_SAMPLE_ATTEMPTS
             + setup_helper_constructions + setup_folded_phase_queries)
    if auto_fourier_terms > MAX_ROOT_TERMS or baseline_fourier_terms > MAX_ROOT_TERMS:
        raise MemoryError("root Fourier component preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": payload,
        "payload_components": {
            "direct_and_helper_live_matrices": direct_and_helper_matrices,
            "retained_joint_law_arrays": retained_joint_laws,
            "converted_joint_law_lists": converted_joint_law_lists,
            "retained_output_laws": retained_output_laws,
            "fft_abs_square_vectors": fft_vectors,
            "work_blocks_and_identity_temps": work_blocks,
            "live_component_tuple_reserve": tuple_bytes,
            "converted_scalar_lists": scalar_lists,
            "archived_baseline_vectors": source_baseline_vectors,
        },
        "work_bounds": {
            "direct_block_products": direct_blocks,
            "direct_phase_queries": direct_phase,
            "direct_shift_ops": direct_shifts,
            "direct_fft_calls": len(INSERTIONS),
            "direct_fft_entries": direct_fft_entries,
            "helper_column_local_terms": helper_column_local,
            "helper_column_phase_queries": helper_column_phase,
            "helper_row_local_terms": helper_row_local,
            "helper_row_phase_queries": helper_row_phase,
            "helper_phase_queries": helper_phase,
            "auto_fourier_component_terms": auto_fourier_terms,
            "diagnostic_fourier_component_terms": diagnostic_fourier_terms,
            "forced_fourier_component_terms": forced_fourier_terms,
            "baseline_fourier_component_terms": baseline_fourier_terms,
            "public_forced_calls": public_forced_calls,
            "setup_helper_constructions": setup_helper_constructions,
            "setup_folded_phase_queries": setup_folded_phase_queries,
            "sample_attempts": MAX_SAMPLE_ATTEMPTS,
            "sample_fourier_component_terms": MAX_SAMPLE_FOURIER_TERMS,
            "sample_column_local_terms": sample_column_local,
            "sample_row_local_terms": sample_row_local,
            "sample_phase_queries": sample_phase_queries,
            "sample_marginal_queries": 2 * WIDTH * MAX_SAMPLE_ATTEMPTS,
            "total_scalar_work": total,
            "archived_source_file_bytes": archive_file_bytes,
            "archived_parsed_source_reserve": archived_source_reserve,
            "auto_cover_names": auto_cover_names,
        },
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_root_terms": MAX_ROOT_TERMS,
    }


def add_counts(total, counters):
    for key, value in (counters or {}).items():
        if isinstance(value, (int, np.integer)):
            total[key] = total.get(key, 0) + int(value)


def ensure_counter(total, increments, caps, label):
    for key, amount in increments.items():
        if total.get(key, 0) + int(amount) > caps[key]:
            raise MemoryError(
                f"{label} exceeds {key} cap: "
                f"{total.get(key, 0)}+{amount}>{caps[key]}"
            )


def direct_joint_law(insertion, W0, W1, counters, caps):
    matrix = np.zeros((Q, PERIOD), dtype=np.complex128)
    for exponent in range(Q):
        ensure_counter(counters, {
            "direct_block_products": (PERIOD // BLOCK) * BLOCK * BLOCK,
            "direct_phase_queries": 2 * PERIOD,
            "direct_shift_ops": 3,
            "modular_pow_queries": 2 * PERIOD,
        }, caps, "ER direct column")
        matrix[exponent] = direct_column(exponent, insertion, W0, W1, counters)
    ensure_counter(counters, {"direct_fft_calls": 1,
                              "direct_fft_entries": Q * PERIOD}, caps,
                   "ER FFT")
    counters["direct_fft_calls"] = counters.get("direct_fft_calls", 0) + 1
    counters["direct_fft_entries"] = counters.get("direct_fft_entries", 0) + Q * PERIOD
    transform = np.fft.fft(matrix, axis=0)
    # FFT is indexed (output, work); public sampler laws are (work, output).
    joint = np.abs(transform) ** 2 / (Q * Q)
    return matrix, joint.T


def cached_sampler_law(sampler, selected_outputs=None, *, ledger=None,
                       category_caps=None, row_bound=None, label="helper"):
    """Use one cached row per work and the helper's Fourier oracle thereafter."""
    counters = defaultdict(int)
    ledger = defaultdict(int) if ledger is None else ledger
    joint = np.zeros((PERIOD, Q), dtype=float)
    proposal = np.zeros((PERIOD, Q), dtype=float)
    accepted = np.zeros((PERIOD, Q), dtype=float)
    work_mass = np.zeros(PERIOD, dtype=float)
    row_summaries = []
    row_checks = []
    selected = []
    outputs = tuple(range(Q)) if selected_outputs is None else tuple(selected_outputs)
    for work in range(PERIOD):
        if category_caps is not None:
            ensure_counter(ledger, row_bound, category_caps,
                           f"{label} row")
        row = sampler.row(work)
        add_counts(counters, row.get("counters", {}))
        add_counts(ledger, row.get("counters", {}))
        norm = float(row["norm"])
        components = tuple(row["components"])
        m = len(components)
        work_mass[work] = norm / Q
        local = defaultdict(int)
        for output in outputs:
            if category_caps is not None:
                ensure_counter(ledger, {
                    "fourier_component_terms": m,
                }, category_caps, f"{label} Fourier")
            before_terms = local.get("fourier_component_terms", 0)
            value = sampler._fourier(row, output, local)
            actual_terms = local.get("fourier_component_terms", 0) - before_terms
            if actual_terms != m:
                raise AssertionError(
                    f"{label} Fourier counter mismatch: expected {m}, "
                    f"got {actual_terms}")
            ledger["fourier_component_terms"] += m
            conditional = float(value.get("conditional_probability") or 0.)
            prop = float(value.get("proposal_probability") or 0.)
            acc = float(value.get("acceptance") or 0.)
            joint[work, output] = work_mass[work] * conditional
            proposal[work, output] = work_mass[work] * prop
            accepted[work, output] = work_mass[work] * prop * acc
            if selected_outputs is not None:
                selected.append({"work": work, "output": output,
                                 "conditional": conditional,
                                 "proposal": prop, "acceptance": acc,
                                 "joint": joint[work, output]})
        add_counts(counters, local)
        row_summaries.append({
            "work": work,
            "norm": norm,
            "component_count": m,
            "stride": int(row.get("stride", 0)),
            "component_tuples": json_safe(components),
        })
        if selected_outputs is None:
            if work_mass[work] == 0:
                row_checks.append({"work": work, "zero_mass": True,
                                   "component_count": m})
            else:
                row_checks.append({
                    "work": work,
                    "zero_mass": False,
                    "component_count": m,
                    "conditional_mass": float(joint[work].sum()
                                               / work_mass[work]),
                    "proposal_mass": float(proposal[work].sum()
                                             / work_mass[work]),
                    "accepted_mass": float(accepted[work].sum()
                                             / work_mass[work]),
                    "accepted_law_error": float(np.max(
                        np.abs(m * accepted[work] - joint[work]))),
                })
    if selected_outputs is None:
        return {"joint": joint, "proposal": proposal, "accepted": accepted,
                "work_mass": work_mass, "rows": row_summaries,
                "row_checks": row_checks, "counters": dict(counters)}
    return {"selected": selected, "row_count": PERIOD,
            "selected_strides": sorted({row["stride"] for row in row_summaries}),
            "counters": dict(counters)}


def baseline_law(sampler, *, ledger=None, category_caps=None):
    cached = cached_sampler_law(
        sampler, ledger=ledger, category_caps=category_caps,
        row_bound={"row_local_terms": H * BLOCK * BLOCK,
                   "phase_queries": H + 2,
                   "early_phase_queries": H + 2,
                   "late_phase_queries": H + 2},
        label="folded C78 baseline")
    return cached


def load_archived_baselines():
    if not OUTPUT_REFERENCE_REPORT.exists() or not FEEDBACK_REFERENCE_REPORT.exists():
        raise FileNotFoundError("frozen archived output/feedback reports missing")
    output = json.loads(OUTPUT_REFERENCE_REPORT.read_text())
    feedback = json.loads(FEEDBACK_REFERENCE_REPORT.read_text())
    output_row = output["rows"][0]
    feedback_row = feedback["rows"][0]
    omitted = np.asarray([
        output_row["laws"][str(v)]["omitted_early"] for v in INSERTIONS],
        dtype=float)
    dephased = np.asarray(feedback_row["dephased_laws"], dtype=float)
    if omitted.shape != (len(INSERTIONS), Q) or dephased.shape != (len(INSERTIONS), Q):
        raise ValueError("archived baseline shape changed")
    if not all(valid_law(row, (Q,)) for row in omitted):
        raise ValueError("archived omission law is not normalized")
    if not all(valid_law(row, (Q,)) for row in dephased):
        raise ValueError("archived feedback law is not normalized")
    return omitted, dephased


def main():
    exp = Experiment("earlier_phase_sampler", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "auto complete joint laws agree with independent ER FFT")
    exp.predict("P2", "cached conditional/proposal/accepted laws normalize")
    exp.predict("P3", "seeded draws and cumulative counters reconcile")
    exp.must_fail("C1", "omitting coherent acceptance changes the joint law")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        work_bounds = report["preflight"]["work_bounds"]
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        if not np.allclose(W0.conj().T @ W0, np.eye(BLOCK), atol=1e-12):
            raise ValueError("invalid W0")
        if not np.allclose(W1.conj().T @ W1, np.eye(BLOCK), atol=1e-12):
            raise ValueError("invalid W1")
        omitted, dephased = load_archived_baselines()
        # The new API is imported only after preflight and archived-source
        # validation, so an unavailable API cannot create science arrays.
        from lab.work_first import EarlierPhaseProgressions, LateWorkProgressions

        direct_laws = {}
        helper_laws, helper_proposals = {}, {}
        mode_errors = {mode: {} for mode in ("auto", "cycle", "dual")}
        row_norm_errors, column_norm_errors = {}, {}
        auto_row_checks = {}
        helper_counters = {}
        helper_stats = {}
        helper_column_counters = defaultdict(int)
        helper_column_errors = {}
        direct_counters = defaultdict(int)
        cover_ledger = defaultdict(int)
        cover_caps = {
            "row_local_terms": work_bounds["helper_row_local_terms"],
            "phase_queries": work_bounds["helper_row_phase_queries"],
            "early_phase_queries": work_bounds["helper_row_phase_queries"],
            "late_phase_queries": work_bounds["helper_row_phase_queries"],
            "fourier_component_terms": (
                work_bounds["auto_fourier_component_terms"]
                + work_bounds["diagnostic_fourier_component_terms"]
                + work_bounds["forced_fourier_component_terms"]),
            "public_forced_calls": work_bounds["public_forced_calls"],
        }
        setup_ledger = defaultdict(int)
        setup_caps = {
            "helper_constructions": work_bounds["setup_helper_constructions"],
            "folded_phase_queries": work_bounds["setup_folded_phase_queries"],
        }
        sample_counter_totals = defaultdict(int)
        sample_ledger = defaultdict(int)
        sample_attempts_total = 0
        public_forced_call_count = 0
        direct_caps = {
            "direct_block_products": work_bounds["direct_block_products"],
            "direct_phase_queries": work_bounds["direct_phase_queries"],
            "direct_shift_ops": work_bounds["direct_shift_ops"],
            "direct_fft_calls": work_bounds["direct_fft_calls"],
            "direct_fft_entries": work_bounds["direct_fft_entries"],
            "modular_pow_queries": work_bounds["direct_phase_queries"],
        }
        selected_checks = {}
        forced_public_checks = {}
        sample_rows = {}
        output_tvs = {"omission": {}, "feedback": {}, "proposal": {}}

        for insertion in INSERTIONS:
            matrix, direct_joint = direct_joint_law(
                insertion, W0, W1, direct_counters, direct_caps)
            direct_laws[str(insertion)] = direct_joint
            column_norm_errors[str(insertion)] = float(
                np.max(np.abs(np.sum(np.abs(matrix) ** 2, axis=1) - 1.)))

            ensure_counter(setup_ledger, {"helper_constructions": 2 * BLOCK**3},
                           setup_caps, "primary helper construction")
            setup_ledger["helper_constructions"] += 2 * BLOCK**3
            primary = EarlierPhaseProgressions(
                PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, phase,
                early_split=insertion, early_phase=phase, cover="auto",
                max_local_terms=1_000_000, max_components=4096,
                max_payload_bytes=MAX_BYTES)
            helper_stats[str(insertion)] = json_safe(primary.stats())
            column_error = 0.
            for exponent in range(Q):
                ensure_counter(helper_column_counters, {
                    "column_local_terms": BLOCK * BLOCK,
                    "early_phase_queries": BLOCK,
                    "late_phase_queries": PERIOD,
                    "phase_queries": BLOCK + PERIOD,
                }, {
                    "column_local_terms": work_bounds["helper_column_local_terms"],
                    "early_phase_queries": work_bounds["helper_column_phase_queries"],
                    "late_phase_queries": work_bounds["helper_column_phase_queries"],
                    "phase_queries": work_bounds["helper_column_phase_queries"],
                }, "helper column")
                observed = primary.column(exponent)
                add_counts(helper_column_counters, observed.get("counters", {}))
                actual = {int(k): complex(v)
                          for k, v in observed["amplitudes"].items()}
                expected = matrix[exponent]
                column_error = max(column_error, max(
                    [abs(actual.get(index, 0j) - expected[index])
                     for index in range(PERIOD)]
                    + [abs(value) for index, value in actual.items()
                       if not 0 <= index < PERIOD]))
            helper_column_errors[str(insertion)] = float(column_error)
            selected_cover = primary.cover
            primary_bound = cover_bounds(insertion, selected_cover)
            cached = cached_sampler_law(
                primary, ledger=cover_ledger, category_caps=cover_caps,
                row_bound={"row_local_terms": primary_bound["row_pairs"],
                           "phase_queries": primary_bound["phase_queries"],
                           "early_phase_queries": primary_bound["phase_queries"],
                           "late_phase_queries": primary_bound["phase_queries"]},
                label=f"auto v={insertion}")
            helper_laws[str(insertion)] = cached["joint"]
            helper_proposals[str(insertion)] = cached["proposal"]
            helper_counters[str(insertion)] = cached["counters"]
            auto_row_checks[str(insertion)] = cached["row_checks"]
            forced_public = []
            for work in (0, 1, PERIOD - 1):
                for output in (0, 1, Q // 2, Q - 1):
                    ensure_counter(cover_ledger, {
                        "row_local_terms": primary_bound["row_pairs"],
                        "phase_queries": primary_bound["phase_queries"],
                        "fourier_component_terms": primary_bound["components"],
                        "public_forced_calls": 1,
                    }, cover_caps, "public forced_joint")
                    forced = primary.forced_joint(work, output)
                    add_counts(cover_ledger, forced.get("counters", {}))
                    cover_ledger["public_forced_calls"] += 1
                    expected = direct_joint[work, output]
                    forced_public.append({
                        "work": work, "output": output,
                        "joint": float(forced["joint_probability"]),
                        "error": abs(float(forced["joint_probability"])
                                      - float(expected)),
                        "component_count": int(forced["component_count"]),
                    })
                    public_forced_call_count += 1
            forced_public_checks[str(insertion)] = forced_public
            row_norm_errors[str(insertion)] = float(np.max(np.abs(
                np.asarray([r["norm"] for r in cached["rows"]])
                - np.sum(np.abs(matrix) ** 2, axis=0))))
            selected_checks[str(insertion)] = {}
            for mode in ("cycle", "dual"):
                ensure_counter(setup_ledger, {"helper_constructions": 2 * BLOCK**3},
                               setup_caps, f"{mode} helper construction")
                setup_ledger["helper_constructions"] += 2 * BLOCK**3
                diagnostic = EarlierPhaseProgressions(
                    PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, phase,
                    early_split=insertion, early_phase=phase, cover=mode,
                    max_local_terms=1_000_000, max_components=4096,
                    max_payload_bytes=MAX_BYTES)
                diagnostic_bound = cover_bounds(insertion, mode)
                selected = cached_sampler_law(
                    diagnostic, selected_outputs=(0, 1, Q // 2, Q-1),
                    ledger=cover_ledger, category_caps=cover_caps,
                    row_bound={"row_local_terms": diagnostic_bound["row_pairs"],
                               "phase_queries": diagnostic_bound["phase_queries"],
                               "early_phase_queries": diagnostic_bound["phase_queries"],
                               "late_phase_queries": diagnostic_bound["phase_queries"]},
                    label=f"{mode} v={insertion}")
                selected["stats"] = json_safe(diagnostic.stats())
                selected_checks[str(insertion)][mode] = selected
                # Compare the selected public values with the independent law.
                mode_errors[mode][str(insertion)] = max(
                    abs(item["joint"] - direct_joint[item["work"], item["output"]])
                    for item in selected["selected"])
            auto_error = float(np.max(np.abs(cached["joint"] - direct_joint)))
            mode_errors["auto"][str(insertion)] = auto_error
            output_tvs["omission"][str(insertion)] = tv(
                direct_joint.sum(axis=0), omitted[insertion])
            output_tvs["feedback"][str(insertion)] = tv(
                direct_joint.sum(axis=0), dephased[insertion])
            output_tvs["proposal"][str(insertion)] = tv(
                direct_joint.sum(axis=0), cached["proposal"].sum(axis=0))

            rng = np.random.default_rng(624 + insertion)
            draws = []
            for _ in range(32):
                if sample_attempts_total + 512 > MAX_SAMPLE_ATTEMPTS:
                    raise MemoryError("sample attempt cap before sampler call")
                primary_bound = cover_bounds(insertion, primary.cover)
                ensure_counter(sample_ledger, {
                    "attempts": 512,
                    "fourier_component_terms": 512 * primary_bound["components"],
                    "row_local_terms": primary_bound["row_pairs"],
                    "phase_queries": primary_bound["phase_queries"] + BLOCK + PERIOD,
                    "column_local_terms": BLOCK * BLOCK,
                }, {
                    "attempts": MAX_SAMPLE_ATTEMPTS,
                    "fourier_component_terms": MAX_SAMPLE_FOURIER_TERMS,
                    "row_local_terms": work_bounds["sample_row_local_terms"],
                    "phase_queries": work_bounds["sample_phase_queries"],
                    "column_local_terms": work_bounds["sample_column_local_terms"],
                }, f"sample v={insertion}")
                ensure_counter(sample_ledger, {
                    "marginal_queries": 2 * WIDTH * 512,
                }, {
                    "marginal_queries": work_bounds["sample_marginal_queries"],
                }, f"sample v={insertion} marginal queries")
                sample = primary.sample(rng, max_attempts=512)
                counters = sample.get("counters", {})
                add_counts(sample_counter_totals, counters)
                sample_attempts_total += int(sample["attempts"])
                sample_ledger["attempts"] += int(sample["attempts"])
                add_counts(sample_ledger, counters)
                draws.append({"work": int(sample["work"]),
                              "output": int(sample["output"]),
                              "attempts": int(sample["attempts"]),
                              "component_count": int(sample["component_count"]),
                              "counters": json_safe(counters)})
            sample_rows[str(insertion)] = draws

        # Exact v=0 C78 baseline: fold the phase g(0..2) into the initial block.
        ensure_counter(setup_ledger, {"folded_phase_queries": BLOCK},
                       setup_caps, "folded W0 phase queries")
        setup_ledger["folded_phase_queries"] += BLOCK
        folded = np.diag(np.asarray([phase(u) for u in range(BLOCK)])) @ W0
        ensure_counter(setup_ledger, {"helper_constructions": 2 * BLOCK**3},
                       setup_caps, "baseline helper construction")
        setup_ledger["helper_constructions"] += 2 * BLOCK**3
        baseline = LateWorkProgressions(
            PERIOD, BLOCK, WIDTH, SPLIT, folded, W1, phase,
            max_local_terms=1_000_000, max_components=4096,
            max_payload_bytes=MAX_BYTES)
        baseline_ledger = defaultdict(int)
        baseline_cached = baseline_law(
            baseline, ledger=baseline_ledger,
            category_caps={
                "row_local_terms": PERIOD * H * BLOCK * BLOCK,
                "phase_queries": PERIOD * (H + 2),
                "early_phase_queries": PERIOD * (H + 2),
                "late_phase_queries": PERIOD * (H + 2),
                "fourier_component_terms": work_bounds["baseline_fourier_component_terms"],
            })
        folded_error = float(np.max(np.abs(
            baseline_cached["joint"] - direct_laws["0"])))
        folded_baseline_valid = valid_law(baseline_cached["joint"], (PERIOD, Q))

        no_acceptance_tvs = {
            str(v): tv(direct_laws[str(v)], helper_proposals[str(v)])
            for v in INSERTIONS
        }
        # Equal-mixture TV is also a normalized quantity; concatenated rows
        # would spuriously multiply TV by the number of insertion positions.
        no_acceptance_tv = float(tv(
            sum((direct_laws[str(v)] for v in INSERTIONS), np.zeros((PERIOD, Q)))
            / len(INSERTIONS),
            sum((helper_proposals[str(v)] for v in INSERTIONS),
                np.zeros((PERIOD, Q))) / len(INSERTIONS)))
        c1 = bool(no_acceptance_tv > 1e-8)

        valid_targets = all(valid_law(direct_laws[str(v)], (PERIOD, Q))
                            and valid_law(helper_laws[str(v)], (PERIOD, Q))
                            for v in INSERTIONS)
        p1 = bool(valid_targets
                  and max(mode_errors["auto"].values()) < TOL
                  and max(mode_errors["cycle"].values()) < TOL
                  and max(mode_errors["dual"].values()) < TOL
                  and max(helper_column_errors.values()) < TOL
                  and max(item["error"] for rows in forced_public_checks.values()
                          for item in rows) < TOL
                  and folded_error < TOL
                  and max(column_norm_errors.values()) < TOL
                  and max(row_norm_errors.values()) < TOL)
        all_row_checks = [row for rows in auto_row_checks.values()
                          for row in rows]
        auto_counter_totals = defaultdict(int)
        for values in helper_counters.values():
            add_counts(auto_counter_totals, values)
        counter_reconciliation = {
            "direct": {
                "direct_block_products": direct_counters.get("direct_block_products", 0)
                    <= work_bounds["direct_block_products"],
                "direct_phase_queries": direct_counters.get("direct_phase_queries", 0)
                    <= work_bounds["direct_phase_queries"],
                "direct_fft_entries": direct_counters.get("direct_fft_entries", 0)
                    <= work_bounds["direct_fft_entries"],
            },
            "helper_columns": {
                "column_local_terms": helper_column_counters.get("column_local_terms", 0)
                    <= work_bounds["helper_column_local_terms"],
                "phase_queries": helper_column_counters.get("phase_queries", 0)
                    <= work_bounds["helper_column_phase_queries"],
                "early_phase_queries": helper_column_counters.get("early_phase_queries", 0)
                    <= work_bounds["helper_column_phase_queries"],
                "late_phase_queries": helper_column_counters.get("late_phase_queries", 0)
                    <= work_bounds["helper_column_phase_queries"],
            },
            "helper_auto_rows": {
                "row_local_terms": auto_counter_totals.get("row_local_terms", 0)
                    <= work_bounds["helper_row_local_terms"],
                "phase_queries": auto_counter_totals.get("phase_queries", 0)
                    <= work_bounds["helper_row_phase_queries"],
                "fourier_component_terms": auto_counter_totals.get(
                    "fourier_component_terms", 0)
                    <= work_bounds["auto_fourier_component_terms"],
            },
            "all_cover_calls": {
                "row_local_terms": cover_ledger.get("row_local_terms", 0)
                    <= work_bounds["helper_row_local_terms"],
                "phase_queries": cover_ledger.get("phase_queries", 0)
                    <= work_bounds["helper_row_phase_queries"],
                "fourier_component_terms": cover_ledger.get(
                    "fourier_component_terms", 0)
                    <= (work_bounds["auto_fourier_component_terms"]
                        + work_bounds["diagnostic_fourier_component_terms"]
                        + work_bounds["forced_fourier_component_terms"]),
                "public_forced_calls": cover_ledger.get(
                    "public_forced_calls", 0)
                    == work_bounds["public_forced_calls"],
            },
            "setup": {
                "helper_constructions": setup_ledger.get(
                    "helper_constructions", 0)
                    == work_bounds["setup_helper_constructions"],
                "folded_phase_queries": setup_ledger.get(
                    "folded_phase_queries", 0)
                    == work_bounds["setup_folded_phase_queries"],
            },
            "folded_baseline": {
                "row_local_terms": baseline_ledger.get("row_local_terms", 0)
                    <= PERIOD * H * BLOCK * BLOCK,
                "phase_queries": baseline_ledger.get("phase_queries", 0)
                    <= PERIOD * (H + 2),
                "fourier_component_terms": baseline_ledger.get(
                    "fourier_component_terms", 0)
                    <= work_bounds["baseline_fourier_component_terms"],
            },
            "samples": {
                "attempts": sample_attempts_total <= work_bounds["sample_attempts"],
                "fourier_component_terms": sample_counter_totals.get(
                    "fourier_component_terms", 0)
                    <= work_bounds["sample_fourier_component_terms"],
                "column_local_terms": sample_ledger.get("column_local_terms", 0)
                    <= work_bounds["sample_column_local_terms"],
                "row_local_terms": sample_ledger.get("row_local_terms", 0)
                    <= work_bounds["sample_row_local_terms"],
                "phase_queries": sample_ledger.get("phase_queries", 0)
                    <= work_bounds["sample_phase_queries"],
                "marginal_queries": sample_ledger.get("marginal_queries", 0)
                    <= work_bounds["sample_marginal_queries"],
            },
        }
        counters_reconciled = all(
            value for category in counter_reconciliation.values()
            for value in category.values())
        row_checks_valid = all(
            item.get("zero_mass", False)
            or (abs(item["conditional_mass"] - 1.) < TOL
                and abs(item["proposal_mass"] - 1.) < TOL
                and abs(item["accepted_mass"] - 1./item["component_count"]) < TOL
                and item["accepted_law_error"] < TOL)
            for item in all_row_checks)
        accepted_rows_finite = all(
            item.get("zero_mass", False)
            or math.isfinite(item["accepted_mass"])
            for item in all_row_checks)
        p2 = bool(all(np.all(np.isfinite(helper_proposals[str(v)]))
                      and np.min(helper_proposals[str(v)]) >= -TOL
                      for v in INSERTIONS)
                  and all(np.all(np.isfinite(helper_laws[str(v)]))
                          and np.min(helper_laws[str(v)]) >= -TOL
                          for v in INSERTIONS)
                  and row_checks_valid
                  and accepted_rows_finite
                  and folded_baseline_valid)
        sample_flat = [item for rows in sample_rows.values() for item in rows]
        p3 = bool(sample_flat and all(
            0 <= item["work"] < PERIOD and 0 <= item["output"] < Q
            and 1 <= item["attempts"] <= 512
            and item["component_count"] >= 1
            and item["counters"].get("work_draws", 0) == 1
            and item["counters"].get("component_draws", 0) == item["attempts"]
            and item["counters"].get("acceptance_draws", 0) == item["attempts"]
            for item in sample_flat))
        p3 = bool(p3 and counters_reconciled)
        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L, "H": H,
                        "insertions": list(INSERTIONS)},
            "direct_joint_laws": {k: v.tolist() for k, v in direct_laws.items()},
            "auto_joint_laws": {k: v.tolist() for k, v in helper_laws.items()},
            "auto_proposal_laws": {k: v.tolist() for k, v in helper_proposals.items()},
            "mode_max_joint_errors": mode_errors,
            "folded_v0_max_joint_error": folded_error,
            "folded_baseline_valid": folded_baseline_valid,
            "column_norm_errors": column_norm_errors,
            "helper_column_errors": helper_column_errors,
            "row_norm_errors": row_norm_errors,
            "output_tvs": output_tvs,
            "no_acceptance_joint_tv": no_acceptance_tv,
            "no_acceptance_joint_tvs_by_insertion": no_acceptance_tvs,
            "selected_cover_checks": selected_checks,
            "forced_public_checks": forced_public_checks,
            "sample_rows": sample_rows,
            "sample_count": len(sample_flat),
            "sample_attempts_total": sample_attempts_total,
            "sample_counter_totals": dict(sample_counter_totals),
            "counter_reconciliation": counter_reconciliation,
            "counters_reconciled": counters_reconciled,
            "public_forced_call_count": public_forced_call_count,
            "row_check_count": len(all_row_checks),
            "row_checks_valid": row_checks_valid,
            "accepted_rows_finite": accepted_rows_finite,
            "auto_row_checks": auto_row_checks,
            "archived_baselines": {
                "output_reference_report": str(OUTPUT_REFERENCE_REPORT),
                "feedback_reference_report": str(FEEDBACK_REFERENCE_REPORT),
                "omission_laws_normalized": True,
                "feedback_dephased_laws_normalized": True,
            },
            "direct_counters": dict(direct_counters),
            "helper_column_counters": dict(helper_column_counters),
            "helper_row_fourier_counters": helper_counters,
            "helper_stats": helper_stats,
            "cover_ledger": dict(cover_ledger),
            "setup_ledger": dict(setup_ledger),
            "setup_caps": setup_caps,
            "baseline_ledger": dict(baseline_ledger),
            "sample_ledger": dict(sample_ledger),
            "folded_baseline_stats": json_safe(baseline.stats()),
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "sample_failure_policy": (
                "A capped sampler failure aborts this audit; no counters from "
                "the failed call are credited."),
            "numeric_payload_is_not_rss": True,
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
        })
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "auto/cycle/dual laws and folded baseline")
        exp.check("P2", p2, "normalized conditional/proposal/accepted laws")
        exp.check("P3", p3, "seeded sampler ranges and counters")
        exp.fail_check("C1", c1, "no-acceptance proposal differs from target")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before sampler audit")
        exp.fail_check("C1", False, "exception before no-acceptance control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report = json_safe(report)
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "r": PERIOD, "b": BLOCK, "t": WIDTH, "s": SPLIT},
        "reference": "ER direct columns/FFT plus opt-in EarlierPhaseProgressions",
        "complete_joint_laws_retained": True,
        "archived_baselines_are_provenance_labeled": True,
        "no_timing_claim": True,
        "max_numeric_payload_bytes": MAX_BYTES,
        "numeric_payload_is_not_rss": True,
        "sample_failure_policy": "failed capped calls abort and are not credited",
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
