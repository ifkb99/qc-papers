"""Independent full-output audit of the fixed N=13 clean-orbit fixture.

The released tiny-fixture physical builder is checked against a separately assembled
12-dimensional matrix/sequential_path reference.  A second reference groups
the q=0 reflection sectors {alpha,-alpha} using the existing sector-basis
helpers.  This is a finite circuit check only: no generic propagator, new
physical synthesis, or scaling claim is made.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The independent full-r sequential_path marginal agrees with the
      builder's physical output for theta=0, pi/7, and pi/2, with no leakage
      and unit mass.
  P2  The independently projected static q=0 groups agree with the full-r
      marginal and have no cross-group branch leakage.
  C1  Omitting the second reflection changes at least one nonzero-theta output
      law; the zero-theta row is retained but is not used as the control.

All dense reference allocations are preflight-capped at 16 MiB.  The fixed
fixture has r=12 and width=3; the report records a conservative payload,
actual sequential calls and t*d^3 work units, not timing, FLOPs or RSS.
The physical builder is an exponential tiny-fixture compiler only.
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
from experiments.experiment_route_regrouping import reflection_orbit, sector_basis


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, BLOCK, PERIOD, WIDTH = 13, 2, 3, 12, 3
THETAS = (0.0, math.pi / 7, math.pi / 2)
FIRST_REFLECTION_ANGLE = math.pi / 3
INITIAL_W_ANGLE = math.pi / 4

REFERENCE_WORK = {}


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
    path = Path("out") / f"clean_orbit_output_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def rx_block(theta, p, q):
    guard((BLOCK, BLOCK), np.complex128, "work Rx block")
    result = np.eye(BLOCK, dtype=complex)
    c, s = math.cos(theta / 2), -1j * math.sin(theta / 2)
    result[p, p] = result[q, q] = c
    result[p, q] = result[q, p] = s
    return result


def work_block(theta):
    """W(theta)=diag(1,e^(i theta/3),1) Rx12(theta/2) Rx01(theta)."""
    diagonal = np.diag((1., np.exp(1j * theta / 3), 1.)).astype(complex)
    return diagonal @ rx_block(theta / 2, 1, 2) @ rx_block(theta, 0, 1)


def repeated_block(block):
    guard((PERIOD, PERIOD), np.complex128, "repeated full-r work block")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for m in range(PERIOD // BLOCK):
        start = m * BLOCK
        result[start:start + BLOCK, start:start + BLOCK] = block
    return result


def shift_matrix(power):
    guard((PERIOD, PERIOD), np.complex128, "controlled modular shift")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for source in range(PERIOD):
        result[(source + power) % PERIOD, source] = 1.
    return result


def cell_reflection():
    guard((PERIOD, PERIOD), np.complex128, "cell reflection")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    M = PERIOD // BLOCK
    for m in range(M):
        for p in range(BLOCK):
            source = BLOCK*m + p
            target = BLOCK*((-m) % M) + p
            result[target, source] = 1.
    return result


def reflection_gate(angle):
    identity = np.eye(PERIOD, dtype=complex)
    R = cell_reflection()
    return math.cos(angle / 2) * identity - 1j * math.sin(angle / 2) * R


def planned_reference_bytes():
    # This deliberately overestimates retained pairs, projections and the
    # internal sequential_path QR/effect work before any reference allocation.
    total = np.dtype(np.complex128).itemsize * (
        80 * PERIOD * PERIOD + 12 * (1 << WIDTH))
    guard_bytes(total, "combined clean-orbit reference payload")
    return total


def reference_path(pairs, initial, output):
    # These instrument actual calls; t*d^3 is an explicitly named work unit,
    # not a claimed exact FLOP count for QR or native arithmetic.
    dimension = len(initial)
    REFERENCE_WORK["sequential_calls"] += 1
    REFERENCE_WORK["width_dimension_cubed_units"] += len(pairs)*dimension**3
    return sequential_path(pairs, initial, output=output)["conditional_path_probability"]


def full_pairs(theta, omit_second=False):
    posts = {}
    initial = repeated_block(work_block(INITIAL_W_ANGLE))
    for i in range(WIDTH):
        w = repeated_block(work_block(((-1)**i) * math.pi / (5 + i)))
        post = w
        if i == 0:
            post = reflection_gate(FIRST_REFLECTION_ANGLE) @ post
        if i == 1 and not omit_second:
            post = reflection_gate(theta) @ post
        posts[i] = post
    pairs = []
    for i in range(WIDTH):
        shift = shift_matrix(1 << i)
        pairs.append((posts[i], posts[i] @ shift))
    return initial, pairs


def full_r_marginal(theta, omit_second=False):
    initial, pairs = full_pairs(theta, omit_second=omit_second)
    initial_state = np.zeros(PERIOD, dtype=complex)
    initial_state[0] = 1.
    initial_state = initial @ initial_state
    result = np.zeros(1 << WIDTH, dtype=float)
    guard((1 << WIDTH,), np.float64, "full-r output marginal")
    for y in range(1 << WIDTH):
        result[y] = reference_path(pairs, initial_state, y)
    return result


def grouped_marginal(theta):
    initial, pairs = full_pairs(theta)
    M = PERIOD // BLOCK
    groups = []
    unseen = set(range(M))
    while unseen:
        alpha = min(unseen)
        group = reflection_orbit(M, (0,), alpha)
        groups.append(group)
        unseen.difference_update(group)
    result = np.zeros(1 << WIDTH, dtype=float)
    leakage = 0.
    max_dimension = 0
    initial_work = initial[:BLOCK, :BLOCK] @ np.eye(BLOCK)[:, 0]
    identity = np.eye(PERIOD, dtype=complex)
    for group in groups:
        g = len(group)
        dimension = g * BLOCK
        max_dimension = max(max_dimension, dimension)
        if max_dimension > 6:
            raise AssertionError("q=0 static group exceeds expected dimension 6")
        guard((PERIOD, dimension), np.complex128, "static group projection")
        basis = np.column_stack([sector_basis(PERIOD, BLOCK, alpha)
                                 for alpha in group])
        projector = basis @ basis.conj().T
        for B0, B1 in pairs:
            leakage = max(leakage,
                          float(np.linalg.norm((identity-projector) @ B0 @ basis)),
                          float(np.linalg.norm((identity-projector) @ B1 @ basis)))
        local_pairs = [(basis.conj().T @ B0 @ basis,
                        basis.conj().T @ B1 @ basis)
                       for B0, B1 in pairs]
        local_initial = np.zeros(dimension, dtype=complex)
        for index in range(g):
            local_initial[index*BLOCK:(index+1)*BLOCK] = initial_work / math.sqrt(g)
        local = np.array([
            reference_path(local_pairs, local_initial, y)
            for y in range(1 << WIDTH)])
        result += (g / M) * local
    return result, leakage, [list(group) for group in groups], max_dimension


def builder_output(theta, omit_second=False):
    module = importlib.import_module("experiments.experiment_clean_orbit_gates")
    fixture = module.Fixture(N=N, a=BASE, b=BLOCK, r=PERIOD)
    return module.physical_output(fixture, WIDTH, theta,
                                  omit_second=omit_second)


def field(result, name):
    if isinstance(result, dict):
        return result[name]
    return getattr(result, name)


def valid_probability_law(values):
    values = np.asarray(values, dtype=float)
    if values.shape != (1 << WIDTH,) or not np.all(np.isfinite(values)):
        return False
    return bool(np.all(values >= -3e-12)
                and abs(float(values.sum()) - 1.) < 3e-10)


def main():
    exp = Experiment("clean_orbit_output", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "independent full-r sequential_path matches compiled physical output")
    exp.predict("P2", "q=0 static sector groups match and have zero cross-group leakage")
    exp.must_fail("C1", "omitting the second nonzero-angle reflection changes the output law")
    REFERENCE_WORK.clear()
    REFERENCE_WORK.update(sequential_calls=0, width_dimension_cubed_units=0)
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "reference_payload_bytes": 0}
    p1 = p2 = c1 = False
    try:
        report["reference_payload_bytes"] = planned_reference_bytes()
        grouped_rows = []
        for theta in THETAS:
            reference = full_r_marginal(theta)
            grouped, leakage, groups, max_dimension = grouped_marginal(theta)
            built = builder_output(theta)
            probabilities = np.asarray(field(built, "probabilities"), dtype=float)
            mass = float(probabilities.sum())
            finite = bool(np.all(np.isfinite(probabilities)))
            compiled_error = float(np.max(np.abs(probabilities-reference)))
            grouped_error = float(np.max(np.abs(grouped-reference)))
            grouped_rows.append({"theta": theta,
                                 "compiled_max_abs_error": compiled_error,
                                 "grouped_max_abs_error": grouped_error,
                                 "grouped_leakage": leakage,
                                 "groups": groups, "max_group_dimension": max_dimension,
                                 "reference_probabilities": reference.tolist(),
                                 "grouped_probabilities": grouped.tolist(),
                                 "built_probabilities": probabilities.tolist(),
                                 "mass": mass, "finite": finite,
                                 "reference_valid": valid_probability_law(reference),
                                 "grouped_valid": valid_probability_law(grouped),
                                 "built_valid": valid_probability_law(probabilities),
                                 "builder_leakage": float(field(built, "leakage")),
                                 "builder_norm": float(field(built, "norm")),
                                 "builder_qubits": int(field(built, "qubits")),
                                 "builder_costs": {key:field(built,key) for key in
                                     ("payload_bound_bytes","state_entries","gate_entry_updates")},
                                 "gate_count": int(field(built, "gate_count"))})
            report["rows"].append(grouped_rows[-1])
        omit_differences = []
        for theta in THETAS[1:]:
            reference = full_r_marginal(theta)
            omitted = full_r_marginal(theta, omit_second=True)
            omit_differences.append(float(np.sum(np.abs(reference-omitted))/2))
        report["omitted_second_tv"] = omit_differences
        omitted_theta = THETAS[1]
        omitted_reference = full_r_marginal(omitted_theta, omit_second=True)
        omitted_built = builder_output(omitted_theta, omit_second=True)
        omitted_probabilities = np.asarray(field(omitted_built, "probabilities"), dtype=float)
        report["compiled_omitted"] = {
            "theta": omitted_theta,
            "reference_probabilities": omitted_reference.tolist(),
            "built_probabilities": omitted_probabilities.tolist(),
            "max_abs_error": float(np.max(np.abs(omitted_reference-omitted_probabilities))),
            "valid": valid_probability_law(omitted_probabilities),
            "leakage": float(field(omitted_built, "leakage")),
            "norm": float(field(omitted_built, "norm")),
            "qubits": int(field(omitted_built, "qubits")),
            "builder_costs": {key:field(omitted_built,key) for key in
                ("payload_bound_bytes","state_entries","gate_entry_updates")},
            "gate_count": int(field(omitted_built, "gate_count")),
        }
        report["compiled_omitted_tv"] = float(np.sum(abs(
            np.asarray(report["rows"][1]["built_probabilities"])-omitted_probabilities))/2)
        report["total_physical_gate_entry_updates"] = sum(
            row["builder_costs"]["gate_entry_updates"] for row in report["rows"]
        ) + report["compiled_omitted"]["builder_costs"]["gate_entry_updates"]
        p1 = all(row["compiled_max_abs_error"] < 3e-10
                 and row["reference_valid"] and row["grouped_valid"]
                 and row["built_valid"] and row["finite"]
                 and abs(row["builder_norm"] - 1.) < 3e-10
                 and row["builder_leakage"] < 3e-10
                 and row["builder_qubits"] == 10
                 for row in report["rows"])
        p2 = all(row["grouped_max_abs_error"] < 3e-10
                 and row["grouped_leakage"] < 3e-10
                 and row["max_group_dimension"] <= 6
                 for row in report["rows"])
        c1 = (report["compiled_omitted_tv"] > 1e-8
              and report["compiled_omitted"]["max_abs_error"] < 3e-10
              and report["compiled_omitted"]["valid"]
              and valid_probability_law(omitted_reference)
              and abs(report["compiled_omitted"]["norm"]-1) < 3e-10
              and report["compiled_omitted"]["qubits"] == 10
              and report["compiled_omitted"]["leakage"] < 3e-10)
        report["status"] = "PASS" if p1 and p2 and c1 else "FAIL"
        exp.check("P1", p1, "compiled physical output matches independent full-r reference")
        exp.check("P2", p2, "static q=0 groups match with no cross-group leakage")
        exp.fail_check("C1", c1, "omitting the second nonzero-theta reflection changes output")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["reference_work"] = dict(REFERENCE_WORK)
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "b": BLOCK, "r": PERIOD, "width": WIDTH},
        "theta_values": THETAS,
        "max_dense_bytes": MAX_DENSE_BYTES,
        "payload_scope": "conservative preflight upper bound including retained matrices, projections and sequential_path internals; not RSS",
        "scalar_work_scope": "actual sequential-call and t*d^3 units, not FLOPs/bit runtime; physical gate-entry counts include all synthesized gates",
        "reference_construction": "independent full-r matrices plus existing sequential_path",
        "no_generic_propagator": True,
        "no_new_gate_synthesis_in_reference": True,
        "builder_role": "released capped tiny-fixture gate builder; not implemented here",
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
