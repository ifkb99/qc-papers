"""Independent full-output reference for the TODO36 uniform-prefix pilot.

Frozen indexed fixture: N=61, a=2, r=60, b=3, and six controls.  The
reference assembles 60-by-60 repeated work blocks and controlled shifts, then
uses the existing ``sequential_path`` instrument for every complete output.
Low-output prefixes are obtained only by grouping those same 64 probabilities;
no prefix propagator or physical truth-table compiler is introduced.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The first output bit is uniform.  The proposed support certificate has
      R=6 and the shift 32 has circular distance 28 modulo 60, exceeding 2R.
  P2  The d=2 and d=3 grouped laws are recorded without a forced outcome:
      their separation test reaches the boundary or is otherwise uncertified.
      All law, helper, branch-unitarity and work-budget predicates hold.
  P3  Main's terminal-W cancellation refinement, derived independently before
      main saw output and after the agent's initial run: R=4 certifies d=2
      on the SAME circuit. Retain the original R=6 evidence as conservative.
  C1  Replacing only the final controlled shift (32) by the identity destroys
      the first-bit uniformity and yields a deterministic parity bit.
  C2  Follow-up independence diagnostic: the observed uniform first two bits
      do not factor from the remaining output; replacing the full law by its
      product of marginals fails. Its marginals are computed from this tiny
      enumerated law, not presumed cheaply available to a sampler.

The last W is kept in its original post-control position.  Its unitary
normalization is checked, and its cancellation after the final work trace is
recorded structurally; no gate is commuted in the reference construction.
All dense allocations and the two 64-output reference runs are preflighted.
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
from lab.semiclassical import sequential_path
from experiments.experiment_clean_orbit_output import work_block


MAX_DENSE_BYTES = 16 * 1024 * 1024
MAX_REFERENCE_CALLS = 128
N = 61
BASE = 2
BLOCK = 3
PERIOD = 60
WIDTH = 6
PREFIX_BITS = (1, 2, 3)
INITIAL_ANGLE = math.pi / 4
POST_ANGLE_4 = -math.pi / 10
POST_ANGLE_5 = math.pi / 11
RADIUS = 6
SETUP = {}


def _report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"uniform_prefix_output_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {value} bytes exceeds 16 MiB")
    return value


def _guard(shape, dtype, label):
    return _guard_bytes(math.prod(int(x) for x in shape)
                        * np.dtype(dtype).itemsize, label)


def _exact_order(base, modulus):
    value = 1
    for order in range(1, modulus):
        value = (value * base) % modulus
        if value == 1:
            return order
    return None


def _repeated_block(block):
    SETUP["repeated_blocks"] = SETUP.get("repeated_blocks", 0)+1
    _guard((PERIOD, PERIOD), np.complex128, "repeated 60x60 work block")
    if block.shape != (BLOCK, BLOCK):
        raise ValueError("work block has the wrong dimension")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for m in range(PERIOD // BLOCK):
        lo = m * BLOCK
        result[lo:lo + BLOCK, lo:lo + BLOCK] = block
    return result


def _shift(power):
    SETUP["shift_matrices"] = SETUP.get("shift_matrices", 0)+1
    _guard((PERIOD, PERIOD), np.complex128, "60x60 controlled shift")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for source in range(PERIOD):
        result[(source + power) % PERIOD, source] = 1.
    return result


def _additive_phase(k=1):
    SETUP["modular_phase_values"] = SETUP.get("modular_phase_values", 0)+PERIOD
    _guard((PERIOD, PERIOD), np.complex128, "G1 phase")
    return np.diag([np.exp(2j * math.pi * k * pow(BASE, j, N) / N)
                    for j in range(PERIOD)]).astype(complex)


def _planned_payload_bytes():
    # sequential_path's documented internal estimate plus retained pairs,
    # initial state, outputs, and a conservative matrix-construction reserve.
    internal = 16 * (5 * WIDTH + 16) * PERIOD * PERIOD
    pairs = 4 * WIDTH * 16 * PERIOD * PERIOD  # both fixtures remain live
    matrices = 8 * 16 * PERIOD * PERIOD
    vectors = 4 * 16 * PERIOD + 2 * (1 << WIDTH) * 8
    total = internal + pairs + matrices + vectors
    return _guard_bytes(total, "aggregate uniform-prefix reference payload")


def _unitarity_error(matrix):
    SETUP["explicit_unitarity_checks"] = SETUP.get("explicit_unitarity_checks", 0)+1
    identity = np.eye(matrix.shape[0], dtype=complex)
    return float(np.linalg.norm(matrix.conj().T @ matrix - identity, ord="fro"))


def _valid_law(values):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (1 << WIDTH,)
                and np.all(np.isfinite(values))
                and np.all(values >= -3e-12)
                and abs(float(values.sum()) - 1.) < 3e-10)


def _group_prefix(values, bits):
    values = np.asarray(values, dtype=float)
    groups = np.zeros(1 << bits, dtype=float)
    for y, probability in enumerate(values):
        groups[y % (1 << bits)] += probability
    return groups


def _direct_collision_count(period, width, bits, radius):
    histories = 1 << bits
    step = (1 << (width - bits)) % period
    return sum(1 for q in range(1, histories)
               if min((step*q) % period, (-step*q) % period) <= 2*radius)


def _pairs(*, final_shift=True):
    initial = _repeated_block(work_block(INITIAL_ANGLE))
    identity = np.eye(PERIOD, dtype=complex)
    post = []
    for index in range(WIDTH):
        if index == 4:
            block = _repeated_block(work_block(POST_ANGLE_4))
            block = _additive_phase(1) @ block
            SETUP["construction_matrix_products"] = SETUP.get("construction_matrix_products", 0)+1
        elif index == 5:
            block = _repeated_block(work_block(POST_ANGLE_5))
        else:
            block = identity
        post.append(block)
    pairs = []
    for index, after in enumerate(post):
        shift = _shift(1 << index) if (index != 5 or final_shift) else identity
        pairs.append((after, after @ shift))
        SETUP["construction_matrix_products"] = SETUP.get("construction_matrix_products", 0)+1
    return initial, pairs, post


def _reference_law(initial, pairs, counters):
    state = np.zeros(PERIOD, dtype=complex)
    state[0] = 1.
    state = initial @ state
    SETUP["initial_matvec_calls"] = SETUP.get("initial_matvec_calls", 0)+1
    _guard((1 << WIDTH,), np.float64, "64-output law")
    result = np.zeros(1 << WIDTH, dtype=float)
    for output in range(1 << WIDTH):
        if counters["sequential_calls"] >= MAX_REFERENCE_CALLS:
            raise MemoryError("reference call cap before sequential_path")
        counters["sequential_calls"] += 1
        counters["width_dimension_cubed_units"] += WIDTH * PERIOD**3
        result[output] = sequential_path(pairs, state, output=output,
                                        max_payload_bytes=MAX_DENSE_BYTES)[
            "conditional_path_probability"]
    return result


def _helper_certificate_rows(radius=RADIUS):
    # This is only a cross-check of the independently counted near collisions;
    # the output reference above does not depend on the helper.
    from lab.periodic import uniform_prefix_certificate
    rows = {}
    for bits in PREFIX_BITS:
        direct = _direct_collision_count(PERIOD, WIDTH, bits, radius)
        helper = uniform_prefix_certificate(PERIOD, WIDTH, bits, radius)
        rows[str(bits)] = {
            "direct_near_collision_count": direct,
            "helper_near_collision_count": int(helper["near_collision_count"]),
            "helper_certified_uniform": bool(helper["certified_uniform"]),
            "step_mod_period": int(helper["step_mod_period"]),
            "match": (direct == int(helper["near_collision_count"])
                      and helper["certified_uniform"] == (direct == 0)),
        }
    return rows


def main():
    exp = Experiment("uniform_prefix_output", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "d=1 grouped output is uniform from the separated-support pilot")
    exp.predict("P2", "d=2 and d=3 are measured but not assigned a forced outcome")
    exp.predict("P3", "terminal-W-cancelled radius four certifies the observed d=2 uniformity")
    exp.must_fail("C1", "identity final shift destroys the d=1 uniformity")
    exp.must_fail("C2", "uniform prefix permits discarding prefix/suffix correlations")
    started = time.perf_counter()
    counters = {"sequential_calls": 0, "width_dimension_cubed_units": 0}
    report = {
        "status": "FAIL", "fixture": {"N": N, "a": BASE, "b": BLOCK,
                                         "r": PERIOD, "width": WIDTH},
        "rows": [], "predictions": {
            "d1": "uniform",
            "d2": "open_boundary",
            "d3": "open_boundary",
        },
    }
    p1 = p2 = p3 = c1 = c2 = False
    try:
        SETUP.clear()
        order = _exact_order(BASE, N)
        report["multiplicative_order"] = order
        if order != PERIOD:
            raise AssertionError(f"expected order {PERIOD}, got {order}")
        # All aggregate guards occur before any retained reference arrays.
        report["planned_payload_bytes"] = _planned_payload_bytes()
        report["reference_call_cap"] = MAX_REFERENCE_CALLS
        report["work_units_cap"] = 128 * WIDTH * PERIOD**3
        report["radius"] = RADIUS
        report["helper_certificate"] = _helper_certificate_rows()
        report["terminal_cancelled_certificate"] = _helper_certificate_rows(4)

        initial, pairs, posts = _pairs(final_shift=True)
        initial_control, control_pairs, control_posts = _pairs(final_shift=False)
        report["initial_unitarity_error"] = _unitarity_error(initial)
        report["control_initial_unitarity_error"] = _unitarity_error(initial_control)
        report["post_unitarity_errors"] = [_unitarity_error(matrix) for matrix in posts]
        report["control_post_unitarity_errors"] = [_unitarity_error(matrix)
                                                   for matrix in control_posts]
        branch_errors = []
        control_branch_errors = []
        for pair in pairs:
            branch_errors.extend(_unitarity_error(matrix) for matrix in pair)
        for pair in control_pairs:
            control_branch_errors.extend(_unitarity_error(matrix) for matrix in pair)
        report["branch_unitarity_max_error"] = max(branch_errors)
        report["control_branch_unitarity_max_error"] = max(control_branch_errors)
        report["final_work_gate"] = {
            "kept_after_control_index": 5,
            "post_is_last_work_gate": True,
            "post_unitarity_error": _unitarity_error(posts[5]),
            "trace_cancellation_identity_error": _unitarity_error(posts[5]),
            "commuted": False,
        }

        full = _reference_law(initial, pairs, counters)
        control = _reference_law(initial_control, control_pairs, counters)
        report["full_output"] = full.tolist()
        report["control_output_final_shift_identity"] = control.tolist()
        report["full_valid"] = _valid_law(full)
        report["control_valid"] = _valid_law(control)
        report["full_mass"] = float(full.sum())
        report["control_mass"] = float(control.sum())
        suffix_size = 1 << (WIDTH-2)
        joint = full.reshape(suffix_size, 4)
        prefix_mass, suffix_mass = joint.sum(axis=0), joint.sum(axis=1)
        product = np.outer(suffix_mass, prefix_mass).reshape(-1)
        conditionals = joint/prefix_mass
        report["independence_diagnostic"] = {
            "prefix_mass": prefix_mass.tolist(), "suffix_mass": suffix_mass.tolist(),
            "product_marginals_tv": float(np.sum(abs(full-product))/2),
            "uniform_full_tv": float(np.sum(abs(full-1/(1 << WIDTH)))/2),
            "conditional_suffix_laws": conditionals.T.tolist(),
            "max_pairwise_conditional_tv": max(float(np.sum(abs(conditionals[:,z]-conditionals[:,zp]))/2)
                                                for z in range(4) for zp in range(z)),
            "scope": "full-law diagnostic; no oracle-free sampler for the suffix supplied",
        }
        report["grouped"] = {}
        report["control_grouped"] = {}
        for bits in PREFIX_BITS:
            grouped = _group_prefix(full, bits)
            control_grouped = _group_prefix(control, bits)
            uniform = np.full(1 << bits, 1. / (1 << bits))
            report["grouped"][str(bits)] = {
                "law": grouped.tolist(),
                "uniform_tv": float(np.sum(np.abs(grouped-uniform))/2),
                "classification": ("uniform" if np.sum(np.abs(grouped-uniform))/2 < 3e-10
                                    else "nonuniform"),
                "valid": bool(np.all(np.isfinite(grouped))
                              and np.all(grouped >= -3e-12)
                              and abs(float(grouped.sum())-1.) < 3e-10),
            }
            report["control_grouped"][str(bits)] = {
                "law": control_grouped.tolist(),
                "uniform_tv": float(np.sum(np.abs(control_grouped-uniform))/2),
                "max_probability": float(control_grouped.max()),
                "valid": bool(np.all(np.isfinite(control_grouped))
                              and np.all(control_grouped >= -3e-12)
                              and abs(float(control_grouped.sum())-1.) < 3e-10),
            }
        report["reference_work"] = dict(counters)
        report["setup_work"] = dict(SETUP)
        report["final_reference_calls"] = counters["sequential_calls"]
        report["final_work_units"] = counters["width_dimension_cubed_units"]
        report["expected_work_units"] = 128 * WIDTH * PERIOD**3
        report["helper_all_matches"] = all(row["match"]
                                            for row in report["helper_certificate"].values())
        d1 = report["grouped"]["1"]
        p1 = (report["full_valid"] and report["control_valid"]
              and d1["classification"] == "uniform"
              and d1["uniform_tv"] < 3e-10)
        p2 = (all(report["grouped"][str(bits)]["valid"] for bits in PREFIX_BITS)
              and all(report["control_grouped"][str(bits)]["valid"] for bits in PREFIX_BITS)
              and report["helper_all_matches"]
              and report["branch_unitarity_max_error"] < 3e-12
              and report["control_branch_unitarity_max_error"] < 3e-12
              and report["initial_unitarity_error"] < 3e-12
              and report["control_initial_unitarity_error"] < 3e-12
              and max(report["post_unitarity_errors"]+report["control_post_unitarity_errors"]) < 3e-12
              and counters["sequential_calls"] == MAX_REFERENCE_CALLS
              and counters["width_dimension_cubed_units"] == 128*WIDTH*PERIOD**3
              and SETUP == {"repeated_blocks": 6, "shift_matrices": 11,
                            "modular_phase_values": 120, "construction_matrix_products": 14,
                            "explicit_unitarity_checks": 40, "initial_matvec_calls": 2})
        p3 = (all(row["match"] for row in report["terminal_cancelled_certificate"].values())
              and report["terminal_cancelled_certificate"]["2"]["helper_certified_uniform"]
              and not report["helper_certificate"]["2"]["helper_certified_uniform"]
              and report["grouped"]["2"]["uniform_tv"] < 3e-10)
        c1 = (report["control_grouped"]["1"]["law"][0] > 1-3e-10
              and report["control_grouped"]["1"]["law"][1] < 3e-10
              and report["control_grouped"]["1"]["uniform_tv"] > 0.49)
        c2 = (report["independence_diagnostic"]["product_marginals_tv"] > 1e-6
              and report["independence_diagnostic"]["uniform_full_tv"] > 1e-6
              and report["grouped"]["2"]["uniform_tv"] < 3e-10
              and conditionals.shape == (suffix_size,4) and np.all(np.isfinite(conditionals))
              and np.min(conditionals) >= 0
              and np.max(abs(conditionals.sum(axis=0)-1.)) < 3e-10
              and abs(float(product.sum())-1.) < 3e-10)
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1, "C2": c2}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "d=1 grouped law is uniform")
        exp.check("P2", p2, "all grouped laws and independent certificate cross-checks valid")
        exp.check("P3", p3, "terminal-W cancellation sharpens the same circuit's support certificate")
        exp.fail_check("C1", c1, "final-shift identity control is deterministic/nonuniform")
        exp.fail_check("C2", c2, "uniform prefix still correlates with remaining output")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
        exp.fail_check("C2", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["reference_call_cap_ok"] = counters["sequential_calls"] <= MAX_REFERENCE_CALLS
    report["work_units_cap_ok"] = counters["width_dimension_cubed_units"] <= 128 * WIDTH * PERIOD**3
    report_path = _report_path()
    ok = exp.finish(report_path=report_path, rows=[report], metadata={
        "fixture": report["fixture"], "prefix_bits": PREFIX_BITS,
        "max_dense_bytes": MAX_DENSE_BYTES,
        "max_reference_calls": MAX_REFERENCE_CALLS,
        "reference": "independent repeated 60x60 matrices + existing sequential_path",
        "no_new_propagator": True, "no_truth_table_compiler": True,
        "final_work_gate_not_commuted": True,
    })
    print(f"report: {report_path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
