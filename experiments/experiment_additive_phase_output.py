"""Independent full-output audit for an additive physical phase insertion.

Fixed indexed fixture: N=13, a=2, b=3, r=12, exponent width 3.  This
reuses the already audited work-block, repeated-block and shift-matrix
constructors, but assembles its own branch schedule and calls the existing
sequential_path oracle.  The released physical builder is imported only for
the compiled comparison.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Physical output agrees with the independent full-r reference for
      interior k=0,1,2 and the end-insertion k=1 row; all laws normalize.
  P2  k=0 interior and k=1 end insertion are null controls, matching the
      no-kick background output. Ideal order finding is a separate baseline,
      not predicted to equal the circuit with fine-work mixers.
  C1  Omitting the nonzero interior kick changes the output law.
  C2  Replacing the coherent initial sector state by coarse dephasing changes
      at least one nonzero interior-kick output law.

The mathematical references are evaluated in floating point with explicit
finite orbit-sized arrays, capped at 16 MiB. They record sequential_path
call counts and t*d^3 work units. This is not a new
propagator, sampler, synthesis claim, or runtime result.
"""
from __future__ import annotations

import importlib
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.semiclassical import sequential_path
from experiments.experiment_clean_orbit_output import (
    BASE, BLOCK, INITIAL_W_ANGLE, N, PERIOD, WIDTH,
    repeated_block, shift_matrix, work_block, planned_reference_bytes,
    reference_path, REFERENCE_WORK,
)
from experiments.experiment_route_regrouping import sector_basis


MAX_DENSE_BYTES = 16 * 1024 * 1024
MAX_PHYSICAL_UPDATES = 100_000_000
K_VALUES = (0, 1, 2)

def guard_bytes(payload, label):
    payload = int(payload)
    if payload < 0 or payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def guard(shape, dtype, label):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    guard_bytes(payload, label)
    return payload


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"additive_phase_output_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def additive_phase(k):
    guard((PERIOD, PERIOD), np.complex128, "additive phase matrix")
    return np.diag([np.exp(2j * math.pi * k * pow(BASE, j, N) / N)
                    for j in range(PERIOD)]).astype(complex)


def full_pairs(k, location="interior"):
    if location not in ("interior", "end"):
        raise ValueError("location must be interior or end")
    posts = {}
    for i in range(WIDTH):
        posts[i] = repeated_block(work_block(((-1)**i) * math.pi / (5+i)))
    if location == "interior":
        posts[1] = additive_phase(k) @ posts[1]
    else:
        posts[2] = additive_phase(k) @ posts[2]
    pairs = []
    for i in range(WIDTH):
        shift = shift_matrix(1 << i)
        pairs.append((posts[i], posts[i] @ shift))
    initial = repeated_block(work_block(INITIAL_W_ANGLE))
    return initial, pairs


def full_marginal(k, location="interior"):
    initial, pairs = full_pairs(k, location)
    initial_state = np.zeros(PERIOD, dtype=complex)
    initial_state[0] = 1.
    initial_state = initial @ initial_state
    result = np.zeros(1 << WIDTH, dtype=float)
    guard((1 << WIDTH,), np.float64, "full additive output law")
    for y in range(1 << WIDTH):
        result[y] = reference_path(pairs, initial_state, y)
    return result


def coarse_dephased_marginal(k):
    """Mixture of initial alpha sectors with identical post-kick branches."""
    initial, pairs = full_pairs(k, "interior")
    M = PERIOD // BLOCK
    result = np.zeros(1 << WIDTH, dtype=float)
    work_initial = initial[:BLOCK, :BLOCK] @ np.eye(BLOCK)[:, 0]
    for alpha in range(M):
        basis = sector_basis(PERIOD, BLOCK, alpha)
        state = basis @ work_initial
        for y in range(1 << WIDTH):
            result[y] += reference_path(pairs, state, y) / M
    return result


def ideal_orderfinding_marginal():
    identity = np.eye(PERIOD, dtype=complex)
    pairs = []
    for i in range(WIDTH):
        shift = shift_matrix(1 << i)
        pairs.append((identity, shift))
    state = np.zeros(PERIOD, dtype=complex)
    state[0] = 1.
    return np.array([reference_path(pairs, state, y)
                     for y in range(1 << WIDTH)])


def valid_law(values):
    values = np.asarray(values, dtype=float)
    return (values.shape == (1 << WIDTH,)
            and np.all(np.isfinite(values))
            and np.all(values >= -3e-12)
            and abs(float(values.sum()) - 1.) < 3e-10)


def field(result, name):
    return result[name] if isinstance(result, dict) else getattr(result, name)


def builder_output(k, location, budget):
    module = importlib.import_module("experiments.experiment_additive_phase")
    return module.physical_output(k, location=location, budget=budget)


def main():
    exp = Experiment("additive_phase_output", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "compiled additive-phase rows match the independent full-r law")
    exp.predict("P2", "zero additive phase and end insertion preserve the background law")
    exp.must_fail("C1", "omitting the nonzero interior kick changes output")
    exp.must_fail("C2", "coarse initial sector dephasing changes an interior-kick law")
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "physical_rows": 0}
    p1 = p2 = c1 = c2 = False
    try:
        report["reference_payload_bytes"] = planned_reference_bytes()
        report["reference_work_scope"] = "sequential_path calls and width*dimension^3 units"
        report["physical_update_cap"] = MAX_PHYSICAL_UPDATES
        REFERENCE_WORK.clear()
        REFERENCE_WORK.update(sequential_calls=0, width_dimension_cubed_units=0)
        physical_budget = [0]
        ideal = ideal_orderfinding_marginal()
        report["ideal_orderfinding"] = ideal.tolist()
        background = full_marginal(0, "interior")
        report["background_no_kick"] = background.tolist()
        references = {0: background}
        for k in K_VALUES:
            reference = full_marginal(k, "interior")
            references[k] = reference
            built = builder_output(k, "interior", physical_budget)
            probabilities = np.asarray(field(built, "probabilities"), dtype=float)
            report["rows"].append({
                "k": k, "location": "interior",
                "reference_probabilities": reference.tolist(),
                "built_probabilities": probabilities.tolist(),
                "max_abs_error": float(np.max(np.abs(reference-probabilities))),
                "reference_valid": bool(valid_law(reference)),
                "built_valid": bool(valid_law(probabilities)),
                "leakage": float(field(built, "leakage")),
                "norm": float(field(built, "norm")),
                "gate_count": int(field(built, "gate_count")),
                "gate_entry_updates": int(field(built, "gate_entry_updates")),
                "qubits": int(field(built, "qubits")),
                "payload_bound_bytes": int(field(built, "payload_bound_bytes")),
                "phase_rotation_count": int(field(built, "phase_rotation_count")),
                "interior_vs_background_tv": float(np.sum(np.abs(reference-background))/2),
                "interior_vs_ideal_tv": float(np.sum(np.abs(reference-ideal))/2),
            })
        end_reference = full_marginal(1, "end")
        end_built = builder_output(1, "end", physical_budget)
        end_probabilities = np.asarray(field(end_built, "probabilities"), dtype=float)
        report["rows"].append({
            "k": 1, "location": "end",
            "reference_probabilities": end_reference.tolist(),
            "built_probabilities": end_probabilities.tolist(),
            "max_abs_error": float(np.max(np.abs(end_reference-end_probabilities))),
            "reference_valid": bool(valid_law(end_reference)),
            "built_valid": bool(valid_law(end_probabilities)),
            "leakage": float(field(end_built, "leakage")),
            "norm": float(field(end_built, "norm")),
            "gate_count": int(field(end_built, "gate_count")),
            "gate_entry_updates": int(field(end_built, "gate_entry_updates")),
            "qubits": int(field(end_built, "qubits")),
            "payload_bound_bytes": int(field(end_built, "payload_bound_bytes")),
            "phase_rotation_count": int(field(end_built, "phase_rotation_count")),
        })
        report["physical_rows"] = len(report["rows"])
        dephased = {k: coarse_dephased_marginal(k) for k in (0, 1, 2)}
        report["dephased"] = {str(k): law.tolist() for k, law in dephased.items()}
        report["dephased_valid"] = {str(k): bool(valid_law(law))
                                    for k, law in dephased.items()}
        report["coarse_dephased_tv"] = {
            str(k): float(np.sum(np.abs(references[k]-dephased[k]))/2)
            for k in (1, 2)}
        report["end_null_tv"] = float(np.sum(np.abs(end_reference-background))/2)
        report["interior_k0_background_tv"] = float(
            np.sum(np.abs(references[0]-background))/2)
        report["dephased_k0_background_tv"] = float(
            np.sum(np.abs(dephased[0]-background))/2)
        report["ideal_background_tv"] = float(np.sum(np.abs(ideal-background))/2)
        interior_probabilities = {
            row["k"]: np.asarray(row["built_probabilities"], dtype=float)
            for row in report["rows"] if row["location"] == "interior"}
        end_probabilities = np.asarray(
            next(row for row in report["rows"] if row["location"] == "end")["built_probabilities"],
            dtype=float)
        report["compiled_k0_end_tv"] = float(
            np.sum(np.abs(interior_probabilities[0]-end_probabilities))/2)
        report["compiled_nonzero_vs_k0_tv"] = {
            str(k): float(np.sum(np.abs(interior_probabilities[k]
                                      - interior_probabilities[0]))/2)
            for k in (1, 2)}
        # Snapshot reference counters only after the final reference call.
        report["sequential_reference_work"] = dict(REFERENCE_WORK)
        report["physical_budget_final"] = physical_budget[0]
        p1 = (len(report["rows"]) == 4
              and all(row["max_abs_error"] < 3e-10
                      and row["reference_valid"] and row["built_valid"]
                      and row["leakage"] < 3e-10
                      and row["qubits"] == 10 and row["phase_rotation_count"] == 5
                      and row["payload_bound_bytes"] <= MAX_DENSE_BYTES
                      and abs(row["norm"]-1.) < 3e-10
                      for row in report["rows"]))
        p2 = (report["interior_k0_background_tv"] < 3e-10
              and report["end_null_tv"] < 3e-10
              and report["dephased_k0_background_tv"] < 3e-10
              and all(report["dephased_valid"].values())
              and valid_law(ideal))
        c1 = (any(value > 1e-8 for value in report["compiled_nonzero_vs_k0_tv"].values())
              and report["compiled_k0_end_tv"] < 3e-10)
        c2 = any(value > 1e-8 for value in report["coarse_dephased_tv"].values())
        total_updates = sum(int(row["gate_entry_updates"]) for row in report["rows"])
        report["physical_gate_entry_updates_sum"] = total_updates
        p1 = (p1 and total_updates <= MAX_PHYSICAL_UPDATES
              and total_updates == physical_budget[0])
        report["status"] = "PASS" if p1 and p2 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "compiled additive-phase rows match independent full-r laws")
        exp.check("P2", p2, "k=0 and end insertion are null controls")
        exp.fail_check("C1", c1, "omitting interior additive kick changes output")
        exp.fail_check("C2", c2, "coarse initial dephasing changes interior kick law")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2"):
            exp.check(name, False, "exception before completion")
        for name in ("C1", "C2"):
            exp.fail_check(name, False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "b": BLOCK, "r": PERIOD, "width": WIDTH},
        "k_values": K_VALUES, "locations": ("interior", "end"),
        "max_dense_bytes": MAX_DENSE_BYTES,
        "max_physical_gate_entry_updates": MAX_PHYSICAL_UPDATES,
        "reference": "independent full-r matrices plus existing sequential_path",
        "no_new_propagator": True, "no_sampler": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
