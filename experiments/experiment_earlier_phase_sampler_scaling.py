"""Bounded sampler comparison for the C79 earlier-phase covers.

This is the frozen untruncated follow-up from TODO39: N=61, a=2, r=60,
b=3, W0=work_block(pi/4), W1=work_block(-pi/10), g(j)=exp(2*pi*i*2**j/61),
s=7..12, t=s+1, and early insertion v=s-1.  H=K=2 throughout.  The
opt-in EarlierPhaseProgressions helper is compared in auto/cycle/dual modes
at work labels 0,1,59.  No Q-by-r orbit/state array is constructed.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Auto, cycle, and dual rows produce the same complete conditional Fourier
      laws as literal indexed columns; auto chooses the smaller safe cover.
  P2  The dual description remains within its fixed component bound 20 for
      every width, and the v=0 early phase agrees with the exact folded-W0
      C78 baseline.
  P3  Bounded auto samples return valid outputs while retaining work across
      rejection; counts and finite attempt caps remain explicit diagnostics.
  C1  A constant/omitted early phase at the interior insertion changes at
      least one selected conditional law.

Literal references call the existing ER direct_column once per exponent in
the union of the proven support residues for work labels 0,1,59, then extract
all three work entries.  Complete laws use small Q-by-3 row matrices and FFT,
not a full Q-by-r array.  Selected outputs additionally use the helper
geometric progression formula, including row strides and gcd lifts.  The
experiment is a float diagnostic, not a timing or numerical certificate and
does not claim superiority over the width-linear full-r sequential baseline.
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
from lab.work_first import EarlierPhaseProgressions, LateWorkProgressions


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTHS = tuple(range(7, 13))  # split s, not register width t
WORKS = (0, 1, 59)
COVERS = ("auto", "cycle", "dual")
SAMPLE_COUNT = 32
SAMPLE_ATTEMPT_CAP = 512
BYTE_CAP = 16 * 1024 * 1024
TERM_CAP = 5_000_000
TOL = 3e-9


def report_path() -> Path:
    path = Path("out") / (
        "earlier_phase_sampler_scaling_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%fZ")
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


def phase(index: int) -> complex:
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def charged_call(worker, method, counters, *args, **kwargs):
    """Reserve before invocation; retain actual counters even on exhaustion.

    Instrumentation only, not another propagator. Early/late phase subcounts
    are never added again to the inclusive phase-query total.
    """
    attempts = kwargs.get("max_attempts", 0) if method == "sample" else 0
    increments = {
        "helper_row_local_terms": worker.local_term_bound-worker.b**2,
        "helper_column_local_terms": worker.b**2 if attempts else 0,
        "helper_phase_queries": worker.stats()["phase_queries_per_draw_bound"],
        "sample_fourier_terms": attempts*worker.component_bound,
        "sample_marginal_queries": 2*worker.width*attempts,
    }
    for key, amount in increments.items():
        if counters[key]+amount > counters[key+"_cap"]:
            raise MemoryError(f"{key} cumulative reserve before {method}")
    captured = []
    original = worker._new_counters
    def capture():
        values = original()
        captured.append(values)
        return values
    worker._new_counters = capture
    try:
        return getattr(worker, method)(*args, **kwargs)
    finally:
        worker._new_counters = original
        mapping = {"row_local_terms": "helper_row_local_terms",
                   "column_local_terms": "helper_column_local_terms",
                   "phase_queries": "helper_phase_queries",
                   "fourier_component_terms": "sample_fourier_terms",
                   "marginal_queries": "sample_marginal_queries"}
        for values in captured:
            for source, target in mapping.items():
                amount = values[source]
                if amount > increments[target]:
                    raise AssertionError(f"actual {source} exceeded reserved bound")
                counters[target] += amount
            counters["helper_early_phase_queries"] += values.get("early_phase_queries", 0)
            counters["helper_late_phase_queries"] += values.get("late_phase_queries", values["phase_queries"])


def valid_law(law: np.ndarray) -> bool:
    return (law.shape and law.ndim == 1 and np.all(np.isfinite(law))
            and np.min(law) >= -TOL and abs(float(np.sum(law)) - 1.) < TOL)


def tv(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.sum(np.abs(first - second)) / 2.)


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def support_exponents(width: int) -> list[int]:
    residues = sorted({(work - delta) % PERIOD
                       for work in WORKS for delta in range(-2, 5)})
    return [rho + PERIOD * n for rho in residues
            for n in range(((1 << width) - 1 - rho) // PERIOD + 1)]


def preflight() -> dict:
    sum_q = sum(1 << (s + 1) for s in WIDTHS)
    max_q = 1 << (max(WIDTHS) + 1)
    direct_calls = 9 * sum(math.ceil((1 << (s + 1)) / PERIOD) for s in WIDTHS)
    direct_block = direct_calls * (PERIOD // BLOCK) * BLOCK**2
    direct_phase = direct_calls * 2 * PERIOD
    base_helper_rows = len(WIDTHS) * len(COVERS) * len(WORKS)
    # Include the folded-v=0 pair and the three-row wrong-phase control in
    # the allocation/operation budget: these are actual helper evaluations.
    extra_helper_rows = 2 * len(WIDTHS) * len(WORKS) + len(WORKS)
    helper_rows = base_helper_rows + extra_helper_rows
    # The auto/dual component bound is 20; cycle's safe bound is derived
    # before construction for the frozen v=s-1 schedule.
    component_bounds = {}
    for s in WIDTHS:
        L = 1 << s
        A = 1 << (s - 1)
        P = A // math.gcd(PERIOD, A)
        cycle_classes = min(P, math.ceil(L / PERIOD))
        cycle = 2 * min(L, 5 * cycle_classes)
        component_bounds[str(s)] = {"auto": 20, "cycle": cycle, "dual": 20}
    max_components = max(v[c] for v in component_bounds.values() for c in COVERS)
    law_fft_entries = helper_rows * max_q
    # A row's progressions are disjoint, so its expanded assignments are
    # bounded by Q regardless of how many components describe it.  Using
    # component_count*Q here would charge the same exponent repeatedly and
    # reject the frozen run for a purely bookkeeping overestimate.
    helper_expansions = helper_rows * 2*5*math.ceil((max_q//2)/PERIOD)
    selected_roots = sum(
        4 * component_bounds[str(s)][cover] * len(WORKS)
        for s in WIDTHS for cover in COVERS
    ) + extra_helper_rows * 4 * 20
    components = {
        "reference_q_by_three_rows": max_q * len(WORKS) * 16,
        "helper_expanded_row": max_q * 16,
        "helper_law_and_fft_temporaries": 4 * max_q * 16,
        "helper_copies_and_check_temporaries": 32 * BLOCK * BLOCK * 16,
        "retained_law_vectors": 4 * sum_q * len(WORKS) * 8,
        "scalar_and_component_reserve": 2 * 1024 * 1024,
    }
    payload = guard_bytes(sum(components.values()), "sampler scaling preflight")
    operation_bounds = {
        "reference_direct_calls": direct_calls,
        "reference_block_products": direct_block,
        "reference_phase_queries": direct_phase,
        "reference_shift_ops": 3 * direct_calls,
        "helper_rows": helper_rows,
        "helper_expansion_assignments": helper_expansions,
        "helper_fft_calls": helper_rows,
        "helper_fft_input_entries": law_fft_entries,
        "reference_fft_calls": len(WIDTHS),
        "reference_fft_input_entries": sum_q * len(WORKS),
        "selected_geometric_root_queries": selected_roots,
        "sample_calls": len(WIDTHS) * SAMPLE_COUNT,
        "sample_attempts": len(WIDTHS) * SAMPLE_COUNT * SAMPLE_ATTEMPT_CAP,
        "max_components": max_components,
        "helper_row_local_terms": helper_rows*18*(1+math.ceil((max_q//2)/PERIOD))
                                  + len(WIDTHS)*SAMPLE_COUNT*36,
        "helper_column_local_terms": len(WIDTHS)*SAMPLE_COUNT*BLOCK**2,
        "helper_phase_queries": helper_rows*(18*math.ceil((max_q//2)/PERIOD)+14)
                                + len(WIDTHS)*SAMPLE_COUNT*26,
        "sample_fourier_terms": len(WIDTHS)*SAMPLE_COUNT*SAMPLE_ATTEMPT_CAP*20,
        "sample_marginal_queries": len(WIDTHS)*SAMPLE_COUNT*SAMPLE_ATTEMPT_CAP*2*(max(WIDTHS)+1),
        "unitary_check_product_terms": (len(WIDTHS)*5+2)*2*BLOCK**3,
        "fold_setup_phase_queries": len(WIDTHS)*BLOCK,
    }
    for key, value in operation_bounds.items():
        if value > TERM_CAP and key not in {"sample_attempts", "max_components"}:
            raise MemoryError(f"{key} preflight exceeds scalar cap")
    return {
        "planned_numeric_payload_bytes": payload,
        "payload_components": components,
        "operation_bounds": operation_bounds,
        "component_bounds_by_width": component_bounds,
        "support_union_residue_count": 9,
        "support_union_call_bound": direct_calls,
        "max_numeric_payload_bytes": BYTE_CAP,
        "max_scalar_terms": TERM_CAP,
    }


def reference_laws(width: int, W0: np.ndarray, W1: np.ndarray,
                   counters: dict) -> tuple[dict, dict, dict]:
    Q = 1 << width
    support = support_exponents(width)
    guard_shape((Q, len(WORKS)), label=f"reference width {width} Q-by-3 rows")
    rows = np.zeros((Q, len(WORKS)), dtype=np.complex128)
    for exponent in support:
        if counters["reference_direct_calls"] + 1 > counters["reference_call_cap"]:
            raise MemoryError("reference direct-column cap before call")
        if counters["reference_block_products"] + (PERIOD // BLOCK) * BLOCK**2 > counters["reference_block_cap"]:
            raise MemoryError("reference block-product cap before call")
        counters["reference_direct_calls"] += 1
        counters["reference_block_products"] += (PERIOD // BLOCK) * BLOCK**2
        counters["reference_shift_ops"] += 3
        # Pass one cumulative ER counter dictionary.  Resetting it for every
        # column would make the existing one-million-term guard meaningless.
        er_cap = int(er.MAX_TERMS)
        if (counters["direct_phase_queries"] + 2 * PERIOD > er_cap
                or counters["direct_block_products"]
                + (PERIOD // BLOCK) * BLOCK**2 > er_cap
                or counters["direct_shift_ops"] + 3 > er_cap):
            raise MemoryError("ER direct-column cumulative cap before call")
        vector = er.direct_column(exponent, width - 2, W0, W1,
                                  counters, split=width - 1)
        for pos, work in enumerate(WORKS):
            rows[exponent, pos] = vector[work]
    norms = np.sum(np.abs(rows) ** 2, axis=0)
    add_count(counters, "reference_fft_calls", 1,
              counters["reference_fft_call_cap"], label="reference FFT")
    add_count(counters, "reference_fft_input_entries", Q * len(WORKS),
              counters["reference_fft_input_cap"], label="reference FFT input")
    raw = np.fft.fft(rows, axis=0)
    laws = {str(work): np.abs(raw[:, pos]) ** 2 / (Q * norms[pos])
            for pos, work in enumerate(WORKS)}
    return laws, {str(work): float(norms[pos] / Q)
                  for pos, work in enumerate(WORKS)}, rows


def expanded_law(row: dict, width: int, counters: dict, label: str,
                 reference: np.ndarray | None = None) -> tuple[np.ndarray, float | None]:
    Q = 1 << width
    guard_shape((Q,), label=f"{label} expanded row")
    expanded = np.zeros(Q, dtype=np.complex128)
    for start, count, gamma in row["components"]:
        for n in range(int(count)):
            add_count(counters, "helper_expansion_assignments", 1,
                      counters["helper_expansion_cap"], label="helper expansion")
            exponent = int(start) + int(row["stride"]) * n
            if not 0 <= exponent < Q or expanded[exponent] != 0:
                raise AssertionError("invalid or colliding helper progression")
            expanded[exponent] = gamma
    add_count(counters, "helper_fft_calls", 1, counters["helper_fft_call_cap"],
              label="helper FFT")
    add_count(counters, "helper_fft_input_entries", Q,
              counters["helper_fft_input_cap"], label="helper FFT input")
    amplitude_error = (None if reference is None else
                       float(np.max(np.abs(expanded - reference))))
    raw = np.fft.fft(expanded)
    return np.abs(raw) ** 2 / (Q * float(row["norm"])), amplitude_error


def selected_geometric(row: dict, width: int, output: int,
                       counters: dict) -> float:
    Q = 1 << width
    total = 0j
    for start, count, gamma in row["components"]:
        add_count(counters, "selected_geometric_root_queries", 1,
                  counters["selected_root_cap"], label="selected geometric root")
        total += gamma * unit_phase(-output * int(start), Q) * geometric_sum(
            int(count), -output * int(row["stride"]), Q)
    return float(abs(total) ** 2 / (Q * float(row["norm"])))


def main() -> None:
    exp = Experiment("earlier_phase_sampler_scaling", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "all covers match literal conditional Fourier laws")
    exp.predict("P2", "dual bound and folded v=0 baseline hold")
    exp.predict("P3", "bounded auto samples exercise fixed-work rejection")
    exp.must_fail("C1", "constant early phase differs at the interior insertion")
    started = time.perf_counter()
    counters = {
        "reference_direct_calls": 0, "reference_block_products": 0,
        "reference_phase_queries": 0, "reference_modular_pow_queries": 0,
        "reference_shift_ops": 0, "reference_fft_calls": 0,
        "reference_fft_input_entries": 0, "helper_rows": 0,
        "direct_phase_queries": 0, "modular_pow_queries": 0,
        "direct_block_products": 0, "direct_shift_ops": 0,
        "helper_row_local_terms": 0, "helper_early_phase_queries": 0,
        "helper_column_local_terms": 0, "helper_phase_queries": 0,
        "sample_fourier_terms": 0, "sample_marginal_queries": 0,
        "unitary_check_product_terms": 0, "fold_setup_phase_queries": 0,
        "helper_late_phase_queries": 0, "helper_expansion_assignments": 0,
        "helper_fft_calls": 0, "helper_fft_input_entries": 0,
        "selected_geometric_root_queries": 0, "sample_calls": 0,
        "sample_attempts": 0, "sample_rejection_steps": 0,
    }
    report = {"status": "FAIL", "rows": [], "predictions": {
        "covers": "auto/cycle/dual laws agree",
        "width_schedule": "s=7..12, t=s+1, early v=s-1, H=K=2",
    }}
    p1 = p2 = p3 = c1 = False
    try:
        report["fixture"] = {"N": N, "base": BASE, "period": PERIOD,
                              "block_size": BLOCK, "splits": list(WIDTHS),
                              "works": list(WORKS), "covers": list(COVERS)}
        report["preflight"] = preflight()
        bounds = report["preflight"]["operation_bounds"]
        for key in ("helper_row_local_terms", "helper_column_local_terms",
                    "helper_phase_queries", "sample_fourier_terms", "sample_marginal_queries"):
            counters[key+"_cap"] = bounds[key]
        counters.update({
            "reference_call_cap": bounds["reference_direct_calls"],
            "reference_block_cap": bounds["reference_block_products"],
            "reference_fft_call_cap": bounds["reference_fft_calls"],
            "reference_fft_input_cap": bounds["reference_fft_input_entries"],
            "helper_expansion_cap": bounds["helper_expansion_assignments"],
            "helper_fft_call_cap": bounds["helper_fft_calls"],
            "helper_fft_input_cap": bounds["helper_fft_input_entries"],
            "helper_row_cap": bounds["helper_rows"],
            "selected_root_cap": bounds["selected_geometric_root_queries"],
        })
        W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
        add_count(counters,"unitary_check_product_terms",2*BLOCK**3,
                  bounds["unitary_check_product_terms"])
        if not np.allclose(W0.conj().T @ W0, np.eye(BLOCK), atol=1e-12, rtol=0):
            raise AssertionError("W0 is not unitary")
        if not np.allclose(W1.conj().T @ W1, np.eye(BLOCK), atol=1e-12, rtol=0):
            raise AssertionError("W1 is not unitary")
        rows_report = {}
        reference_law_report = {}
        helper_law_report = {}
        folded_report = {}
        auto_rows = {}
        for s in WIDTHS:
            # WIDTHS stores split s; the exponent register has t=s+1 bits
            # and the extra phase is inserted at v=s-1.
            width = s + 1
            row_key = str(s)
            Q = 1 << width
            reference_laws_by_work, work_mass, _reference_rows = reference_laws(
                width, W0, W1, counters)
            reference_law_report[row_key] = {
                str(work): law.tolist() for work, law in reference_laws_by_work.items()}
            workers = {}
            stats_by_cover = {}
            for cover in COVERS:
                add_count(counters,"unitary_check_product_terms",2*BLOCK**3,
                          bounds["unitary_check_product_terms"])
                worker = EarlierPhaseProgressions(
                    PERIOD, BLOCK, width, s, W0, W1, phase,
                    early_split=s - 1, early_phase=phase, cover=cover,
                    max_local_terms=1_000_000, max_components=4096,
                    max_payload_bytes=BYTE_CAP)
                workers[cover] = worker
                stats_by_cover[cover] = json_safe(worker.stats())
            width_rows = {}
            helper_law_report[row_key] = {}
            auto_rows[row_key] = {}
            for cover, worker in workers.items():
                helper_law_report[row_key][cover] = {}
                width_rows[cover] = {}
                for pos, work in enumerate(WORKS):
                    add_count(counters, "helper_rows", 1, counters["helper_row_cap"],
                              label="helper row")
                    row = charged_call(worker,"row",counters,work)
                    law, amplitude_error = expanded_law(
                        row, width, counters,
                        f"width={width},cover={cover},work={work}",
                        reference=_reference_rows[:, pos])
                    selected = {}
                    for output in (0, 1, Q // 3, Q - 1):
                        selected[str(output)] = selected_geometric(
                            row, width, output, counters)
                    target = reference_laws_by_work[str(work)]
                    law_error = tv(law, target)
                    selected_error = max(abs(selected[str(output)] - target[output])
                                         for output in map(int, selected))
                    width_rows[cover][str(work)] = {
                        "cover": row.get("cover", cover),
                        "stride": int(row["stride"]),
                        "component_count": len(row["components"]),
                        "longest_progression": max(count for _,count,_ in row["components"]),
                        "occupied_exponents": sum(count for _,count,_ in row["components"]),
                        "construction_counters": dict(row["counters"]),
                        "norm": float(row["norm"]),
                        "work_probability": float(row["norm"] / Q),
                        "work_probability_error": abs(
                            float(row["norm"] / Q) - work_mass[str(work)]),
                        "amplitude_max_error": amplitude_error,
                        "conditional_tv_vs_reference": law_error,
                        "selected_geometric_max_error": float(selected_error),
                        "selected_geometric_outputs": selected,
                        "law_normalized": bool(valid_law(law)),
                        "law": law.tolist(),
                    }
                    if cover == "auto":
                        auto_rows[row_key][str(work)] = row
                # A fresh row audit is used only for data; no Q-by-r table is kept.
            # Exact v=0 folded-W0 C78 baseline for this width.
            add_count(counters,"fold_setup_phase_queries",BLOCK,bounds["fold_setup_phase_queries"])
            folded_initial = np.array([phase(u) for u in range(BLOCK)])[:,None] * W0
            add_count(counters,"unitary_check_product_terms",4*BLOCK**3,
                      bounds["unitary_check_product_terms"])
            folded = LateWorkProgressions(
                PERIOD, BLOCK, width, s, folded_initial, W1, phase,
                max_local_terms=1_000_000, max_components=4096,
                max_payload_bytes=BYTE_CAP)
            folded_rows = {}
            v0_worker = EarlierPhaseProgressions(
                PERIOD, BLOCK, width, s, W0, W1, phase,
                early_split=0, early_phase=phase, cover="auto",
                max_local_terms=1_000_000, max_components=4096,
                max_payload_bytes=BYTE_CAP)
            for work in WORKS:
                add_count(counters, "helper_rows", 1, counters["helper_row_cap"],
                          label="folded helper row")
                folded_row = charged_call(folded,"row",counters,work)
                add_count(counters, "helper_rows", 1, counters["helper_row_cap"],
                          label="v0 helper row")
                early_row = charged_call(v0_worker,"row",counters,work)
                folded_law, _ = expanded_law(folded_row, width, counters, "folded v0")
                early_law, _ = expanded_law(early_row, width, counters, "early v0")
                folded_rows[str(work)] = {
                    "tv_early_vs_folded": tv(early_law, folded_law),
                    "tv_folded_vs_reference": tv(folded_law, reference_laws_by_work[str(work)]),
                    "stride_early": int(early_row["stride"]),
                    "stride_folded": int(folded_row["stride"]),
                }
            folded_report[row_key] = folded_rows
            rows_report[row_key] = {
                "reference_work_mass": work_mass,
                "stats_by_cover": stats_by_cover,
                "covers": width_rows,
                "auto_selected_covers": {
                    cover: stats_by_cover[cover]["cover"] for cover in COVERS},
            }
            # Bounded actual auto samples, preserving work across retries.
            sample_rows = []
            for index in range(SAMPLE_COUNT):
                add_count(counters,"sample_calls",1,bounds["sample_calls"])
                remaining = len(WIDTHS) * SAMPLE_COUNT * SAMPLE_ATTEMPT_CAP - counters["sample_attempts"]
                attempt_cap = min(SAMPLE_ATTEMPT_CAP, remaining)
                if attempt_cap <= 0:
                    raise MemoryError("sample attempt cap before call")
                result = charged_call(workers["auto"],"sample",counters,
                    np.random.default_rng(91_000 + width * 100 + index),
                    max_attempts=attempt_cap)
                attempts = int(result["attempts"])
                counters["sample_attempts"] += attempts
                counters["sample_rejection_steps"] += max(0, attempts - 1)
                sample_counters = result["counters"]
                sample_ok = bool(
                    sample_counters.get("work_draws") == 1
                    and sample_counters.get("component_draws") == attempts
                    and sample_counters.get("progression_proposals") == attempts
                    and sample_counters.get("acceptance_draws") == attempts
                    and sample_counters["fourier_component_terms"] == attempts*result["component_count"]
                    and sample_counters["phase_queries"] == sample_counters["early_phase_queries"]+sample_counters["late_phase_queries"]
                    and 0 <= result["output"] < Q and 0 <= result["work"] < PERIOD
                    # Sampling draws an arbitrary work label from the full
                    # period; only 0,1,59 have retained reference rows.  The
                    # structural bound, rather than a missing selected-row
                    # lookup, is the valid independent check here.
                    and 0 < int(result["component_count"])
                    <= int(stats_by_cover["auto"]["component_bound"]))
                sample_rows.append({
                    "seed": 91_000 + width * 100 + index,
                    "work": int(result["work"]), "output": int(result["output"]),
                    "attempts": attempts, "attempt_cap": attempt_cap,
                    "component_count": int(result["component_count"]),
                    "sample_counter_ok": sample_ok,
                    "counters": json_safe(sample_counters),
                })
            rows_report[row_key]["samples"] = sample_rows
        add_count(counters,"unitary_check_product_terms",2*BLOCK**3,
                  bounds["unitary_check_product_terms"])
        wrong = EarlierPhaseProgressions(
            PERIOD, BLOCK, WIDTHS[0] + 1, WIDTHS[0], W0, W1, phase,
            early_split=WIDTHS[0] - 1, early_phase=lambda _index: 1.+0j,
            cover="auto", max_local_terms=1_000_000,
            max_components=4096, max_payload_bytes=BYTE_CAP)
        wrong_laws = []
        for work in WORKS:
            add_count(counters, "helper_rows", 1, counters["helper_row_cap"],
                      label="wrong-phase helper row")
            wrong_row = charged_call(wrong,"row",counters,work)
            wrong_law, _ = expanded_law(wrong_row, WIDTHS[0] + 1, counters, "wrong phase")
            wrong_laws.append(tv(wrong_law,
                                 np.asarray(reference_law_report[str(WIDTHS[0])][str(work)])))
        wrong_phase_tv = max(wrong_laws)
        p1 = all(
            rows_report[str(s)]["covers"][cover][str(work)]["conditional_tv_vs_reference"] < TOL
            and rows_report[str(s)]["covers"][cover][str(work)]["selected_geometric_max_error"] < TOL
            and rows_report[str(s)]["covers"][cover][str(work)]["work_probability_error"] < TOL
            and rows_report[str(s)]["covers"][cover][str(work)]["amplitude_max_error"] < TOL
            and rows_report[str(s)]["covers"][cover][str(work)]["law_normalized"]
            for s in WIDTHS for cover in COVERS for work in WORKS)
        p2 = all(
            rows_report[str(s)]["covers"][cover][str(work)]["component_count"]
                <= 20
            and rows_report[str(s)]["covers"][cover][str(work)]["stride"] > 0
            and folded_report[str(s)][str(work)]["tv_early_vs_folded"] < TOL
            for s in WIDTHS for cover in ("auto", "dual") for work in WORKS)
        p2 = p2 and all(rows_report[str(s)]["stats_by_cover"]["auto"]["cover"] == "dual"
                        for s in WIDTHS)
        p2 = p2 and all(rows_report[str(WIDTHS[-1])]["covers"]["dual"][str(work)]["longest_progression"] > 1
                        for work in WORKS)
        p3 = (all(row["sample_counter_ok"] for s in WIDTHS
                  for row in rows_report[str(s)]["samples"])
              and counters["sample_calls"] == len(WIDTHS) * SAMPLE_COUNT
              and counters["sample_attempts"] <= bounds["sample_attempts"]
              and all(counters[key] <= bound for key,bound in bounds.items()
                      if key in counters)
              and counters["helper_phase_queries"] == counters["helper_early_phase_queries"]+counters["helper_late_phase_queries"])
        c1 = wrong_phase_tv > 1e-6
        # Keep the report's reference-labelled totals synchronized with the
        # cumulative counters consumed by er.direct_column.
        counters["reference_phase_queries"] = counters["direct_phase_queries"]
        counters["reference_modular_pow_queries"] = counters["modular_pow_queries"]
        report.update({
            "rows": rows_report,
            "reference_laws": reference_law_report,
            "folded_v0_baseline": folded_report,
            "wrong_phase_tvs": wrong_laws,
            "wrong_phase_max_tv": wrong_phase_tv,
            "counters": counters,
            "cost_bounds": bounds,
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "all complete conditional laws and selected roots agree")
        exp.check("P2", p2, "dual bound and folded v=0 baseline hold")
        exp.check("P3", p3, "bounded auto samples exercise fixed-work rejection")
        exp.fail_check("C1", c1, "constant early phase changes the interior law")
    except Exception as exc:
        report["counters"] = counters
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before sampler scaling audit")
        exp.check("P2", False, "exception before sampler scaling audit")
        exp.check("P3", False, "exception before sampler scaling audit")
        exp.fail_check("C1", False, "exception before wrong phase control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": report.get("fixture"),
        "reference": "ER direct columns on nine support residues, Q-by-3 rows and FFT",
        "helper": "EarlierPhaseProgressions auto/cycle/dual",
        "no_Q_by_r_array": True,
        "no_timing_claim": True,
        "max_numeric_payload_bytes": BYTE_CAP,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
