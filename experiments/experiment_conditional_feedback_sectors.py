"""Bounded test of a high-history-dephased conditional suffix candidate.

This keeps the known inverse-QFT feedback after conditioning on the first two
output bits.  On the frozen N=61, a=2, r=60, b=3, t=6 fixture, the candidate
uses the existing 3-dimensional coarse-sector shifts for the four remaining
controls and averages their exact sequential_path laws over all 20 sectors.
An independently assembled 60-dimensional version of the same candidate is
also evaluated.  The complete six-control target and a target with G_1 omitted
are evaluated with the existing repeated-block reference construction.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The complete reference, all four candidate conditional laws, and their
      joint laws normalize; the first-two-bit target masses are uniform.
  P2  Sector averaging and the independent full-r candidate agree.  Candidate
      versus target TV is recorded and classified at 1e-3 and 1e-2, without
      predicting either classification.
  C1  A normalized one-low-input Fourier-basis control exposes the feedback
      sign: omitting the nonzero feedback phase gives a different shifted
      finite-Fourier law.

The experiment uses only the existing ``PeriodicOrbitCircuit`` shift builder,
``sequential_path``, and the existing 60-dimensional fixture constructors.
All 1472 forced calls, their dimension-cubed work units, and retained numeric
payload are preflighted below 16 MiB.  This is a bounded comparator, not a
sampler or an exact-target claim.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import sequential_path
from experiments import experiment_uniform_prefix_output as up


MAX_DENSE_BYTES = 16 << 20
MAX_CALLS = 2_000
N, BASE, BLOCK, PERIOD, WIDTH = 61, 2, 3, 60, 6
PREFIX_BITS, H, LOW_WIDTH = 2, 4, 4
Q = 1 << WIDTH
SECTORS = PERIOD // BLOCK
INITIAL_ANGLE = math.pi / 4
POST_ANGLE_4 = -math.pi / 10
POST_ANGLE_5 = math.pi / 11
SYNTH_CALLS = 3 * (1 << LOW_WIDTH)
EXPECTED_CALLS = 64 + 64 + (H * SECTORS * (1 << LOW_WIDTH)) + 64 + SYNTH_CALLS
MAX_WORK_UNITS = 230_000_000
UP_REPORT_PATH = Path("out/uniform_prefix_output_20260911T101723682629Z.json")
SETUP_WORK = {}


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"conditional_feedback_sectors_{stamp}.json"


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def guard_matrix(shape, label):
    return guard_bytes(math.prod(int(x) for x in shape) * 16, label)


def preflight():
    # Retain target pairs, omission pairs, and full-r candidate pairs together;
    # sector matrices are one 3x3 pair at a time.  Include two simultaneous
    # sequential workspaces and all result vectors.
    target_retained = (1 + 2 * WIDTH) * PERIOD * PERIOD * 16
    omission_retained = 2 * WIDTH * PERIOD * PERIOD * 16
    candidate_retained = (1 + 2 * LOW_WIDTH) * PERIOD * PERIOD * 16
    max_target_workspace = 16 * (5 * WIDTH + 16) * PERIOD * PERIOD
    max_candidate_workspace = 16 * (5 * LOW_WIDTH + 16) * PERIOD * PERIOD
    outputs = (2 * Q * 8) + (5 * H * (1 << LOW_WIDTH) * 8)
    sectors_one = 2 * LOW_WIDTH * BLOCK * BLOCK * 16
    matrix_construction_reserve = 16 * (4 * PERIOD * PERIOD + 2 * BLOCK * BLOCK)
    total = (target_retained + omission_retained + candidate_retained
             + max_target_workspace + max_candidate_workspace + outputs
             + sectors_one + matrix_construction_reserve)
    guard_bytes(total, "aggregate reference/candidate numerical payload")
    if EXPECTED_CALLS > MAX_CALLS:
        raise MemoryError("forced sequential call preflight exceeds cap")
    target_units = Q * WIDTH * PERIOD**3
    candidate_units = H * (1 << LOW_WIDTH) * LOW_WIDTH * PERIOD**3
    sector_units = H * SECTORS * (1 << LOW_WIDTH) * LOW_WIDTH * BLOCK**3
    omission_units = Q * WIDTH * PERIOD**3
    synthetic_units = SYNTH_CALLS * LOW_WIDTH
    if target_units+candidate_units+sector_units+omission_units+synthetic_units > MAX_WORK_UNITS:
        raise MemoryError("forced work-unit preflight exceeds cap")
    return {
        "planned_payload_bytes": total,
        "payload_components": {
            "target_retained": target_retained,
            "omission_retained": omission_retained,
            "candidate_full_retained": candidate_retained,
            "max_target_workspace": max_target_workspace,
            "max_candidate_workspace": max_candidate_workspace,
            "outputs": outputs,
            "one_sector_pair": sectors_one,
            "matrix_construction_reserve": matrix_construction_reserve,
        },
        "planned_calls": EXPECTED_CALLS,
        "planned_work_units": (target_units + candidate_units + sector_units
                                + omission_units + synthetic_units),
        "work_components": {
            "target": target_units, "candidate_full": candidate_units,
            "candidate_sectors": sector_units, "omission": omission_units,
            "synthetic": synthetic_units,
        },
    }


def state_from_initial(initial, dim, label):
    guard_matrix((dim, dim), f"{label} initial matrix")
    guard_bytes(dim * 16, f"{label} initial vector")
    SETUP_WORK["initial_matvec_calls"] = SETUP_WORK.get("initial_matvec_calls", 0) + 1
    state = np.zeros(dim, dtype=complex)
    state[0] = 1.
    state = initial @ state
    if abs(float(np.vdot(state, state).real) - 1.) > 3e-12:
        raise ArithmeticError(f"{label} initial state is not normalized")
    return state


def forced_law(initial, pairs, width, dim, counters, label):
    """Enumerate a complete finite-bit law through existing sequential_path."""
    if len(pairs) != width:
        raise ValueError(f"{label} pair width mismatch")
    guard_bytes((1 << width) * 8, f"{label} output law")
    result = np.zeros(1 << width, dtype=float)
    for output in range(1 << width):
        next_work = counters["work_units"] + width * dim**3
        if counters["calls"] >= MAX_CALLS:
            raise MemoryError("sequential call cap before forced call")
        if next_work > MAX_WORK_UNITS:
            raise MemoryError("sequential work-unit cap before forced call")
        counters["calls"] += 1
        counters["work_units"] = next_work
        result[output] = sequential_path(
            pairs, initial, output=output, max_payload_bytes=MAX_DENSE_BYTES
        )["conditional_path_probability"]
    if (not np.all(np.isfinite(result)) or np.min(result) < -3e-12
            or abs(float(result.sum()) - 1.) > 3e-10):
        raise ArithmeticError(f"{label} law failed normalization")
    return result


def phase_pair(identity, shift, z, bit, sign=-1):
    guard_matrix(identity.shape, "phase-pair matrix")
    if (SETUP_WORK.get("phase_pair_calls",0) >= 400
            or SETUP_WORK.get("phase_pair_matrix_entries",0)+identity.size > 65_000):
        raise MemoryError("phase-pair setup cap before multiplication")
    SETUP_WORK["phase_pair_calls"] = SETUP_WORK.get("phase_pair_calls", 0) + 1
    SETUP_WORK["phase_pair_matrix_entries"] = (
        SETUP_WORK.get("phase_pair_matrix_entries", 0) + identity.shape[0]**2)
    phase = np.exp(sign * 2j * np.pi * z * (1 << bit) / Q)
    return identity, phase * shift


def candidate_sector_laws(counters):
    block = up.work_block(INITIAL_ANGLE)
    circuit = PeriodicOrbitCircuit(PERIOD, BLOCK, LOW_WIDTH, {0: block})
    identity = np.eye(BLOCK, dtype=complex)
    initial = np.asarray(block[:, 0], dtype=complex)
    guard_bytes(H * (1 << LOW_WIDTH) * 8, "sector candidate laws")
    laws = np.zeros((H, 1 << LOW_WIDTH), dtype=float)
    SETUP_WORK["sector_initial_vector_extractions"] = 1
    sector_rows = []
    for z in range(H):
        guard_bytes((1 << LOW_WIDTH) * 8, f"sector z={z} mixture")
        mixture = np.zeros(1 << LOW_WIDTH, dtype=float)
        for alpha in range(SECTORS):
            pairs = []
            for bit in range(LOW_WIDTH):
                SETUP_WORK["sector_shift_calls"] = SETUP_WORK.get("sector_shift_calls", 0) + 1
                pairs.append(phase_pair(identity, circuit.shift(1 << bit, alpha), z, bit))
            law = forced_law(initial, pairs, LOW_WIDTH, BLOCK, counters,
                             f"sector z={z} alpha={alpha}")
            mixture += law / SECTORS
        laws[z] = mixture
        sector_rows.append({"z": z, "mass": float(mixture.sum()),
                            "min": float(mixture.min())})
    return laws, sector_rows, circuit.stats()


def candidate_full_laws(counters):
    block = up.work_block(INITIAL_ANGLE)
    initial_matrix = up._repeated_block(block)
    initial = state_from_initial(initial_matrix, PERIOD, "full candidate")
    identity = np.eye(PERIOD, dtype=complex)
    shifts = [up._shift(1 << bit) for bit in range(LOW_WIDTH)]
    guard_bytes(H * (1 << LOW_WIDTH) * 8, "full candidate laws")
    laws = np.zeros((H, 1 << LOW_WIDTH), dtype=float)
    for z in range(H):
        pairs = [phase_pair(identity, shifts[bit], z, bit)
                 for bit in range(LOW_WIDTH)]
        laws[z] = forced_law(initial, pairs, LOW_WIDTH, PERIOD, counters,
                             f"full candidate z={z}")
    return laws


def target_and_omission(counters):
    # Existing audited complete fixture; no new propagator is introduced.
    initial_matrix, target_pairs, _ = up._pairs(final_shift=True)
    target_initial = state_from_initial(initial_matrix, PERIOD, "target")
    target = forced_law(target_initial, target_pairs, WIDTH, PERIOD, counters,
                        "complete target")

    # Remove only G_1 from the existing i=4 branches by its exact inverse;
    # all other schedule matrices remain unchanged.
    g_inverse = up._additive_phase(-1)
    omission_pairs = list(target_pairs)
    b0, b1 = omission_pairs[4]
    SETUP_WORK["omission_matrix_products"] = SETUP_WORK.get("omission_matrix_products", 0) + 2
    omission_pairs[4] = (g_inverse @ b0, g_inverse @ b1)
    omission = forced_law(target_initial, omission_pairs, WIDTH, PERIOD, counters,
                          "G1 omission baseline")
    return target, omission, {
        "target_initial_norm2": float(np.vdot(target_initial, target_initial).real),
        "target_pair_count": len(target_pairs),
        "omission_pair_count": len(omission_pairs),
    }


def synthetic_feedback_control(counters):
    """Flat low-input Fourier state, z=1, with the same Q=64 feedback."""
    q = 1 << LOW_WIDTH
    z = 1
    identity = np.eye(1, dtype=complex)
    initial = np.ones(1, dtype=complex)
    correct_pairs = [phase_pair(identity, identity, z, bit)
                     for bit in range(LOW_WIDTH)]
    omitted_pairs = [(identity, identity) for _ in range(LOW_WIDTH)]
    wrong_pairs = [phase_pair(identity, identity, z, bit, sign=1)
                   for bit in range(LOW_WIDTH)]
    correct = forced_law(initial, correct_pairs, LOW_WIDTH, 1, counters,
                         "synthetic feedback")
    omitted = forced_law(initial, omitted_pairs, LOW_WIDTH, 1, counters,
                         "synthetic no-feedback")
    wrong_sign = forced_law(initial, wrong_pairs, LOW_WIDTH, 1, counters,
                            "synthetic wrong-sign")
    guard_bytes(3 * q * 8, "synthetic geometric laws")
    expected_correct = np.zeros(q, dtype=float)
    expected_omitted = np.zeros(q, dtype=float)
    expected_wrong = np.zeros(q, dtype=float)
    SETUP_WORK["synthetic_formula_terms"] = 0
    for y in range(q):
        for offset, destination in ((z,expected_correct),(0,expected_omitted),
                                    (-z,expected_wrong)):
            total=0j
            for l in range(q):
                if SETUP_WORK["synthetic_formula_terms"] >= 768:
                    raise MemoryError("synthetic geometric-term cap")
                total += np.exp(-2j*np.pi*(4*y+offset)*l/Q)
                SETUP_WORK["synthetic_formula_terms"] += 1
            destination[y] = abs(total/q)**2
    return {
        "q": q, "Q": Q, "z": z, "boundary": "flat low-input Fourier state",
        "correct_law": correct.tolist(), "omitted_law": omitted.tolist(),
        "wrong_sign_law": wrong_sign.tolist(),
        "expected_correct_law": expected_correct.tolist(),
        "expected_omitted_law": expected_omitted.tolist(),
        "expected_wrong_sign_law": expected_wrong.tolist(),
        "correct_formula_max_abs_error": float(np.max(abs(correct-expected_correct))),
        "omitted_formula_max_abs_error": float(np.max(abs(omitted-expected_omitted))),
        "wrong_formula_max_abs_error": float(np.max(abs(wrong_sign-expected_wrong))),
        "correct_argmax": int(np.argmax(correct)),
        "omitted_argmax": int(np.argmax(omitted)),
        "wrong_sign_argmax": int(np.argmax(wrong_sign)),
        "correct_mass": float(correct.sum()),
        "omitted_mass": float(omitted.sum()),
        "wrong_sign_mass": float(wrong_sign.sum()),
        "omission_tv": float(np.sum(abs(correct-omitted))/2),
        "wrong_sign_tv": float(np.sum(abs(correct-wrong_sign))/2),
        "formula_terms": 3 * q * q,
    }


def source_report_crosscheck(target):
    if not UP_REPORT_PATH.exists():
        return {"available": False}
    payload = json.load(open(UP_REPORT_PATH))
    source = payload["rows"][0]["full_output"]
    guard_bytes(len(source) * 8, "frozen source full law")
    source = np.asarray(source, dtype=float)
    return {"available": True, "path": str(UP_REPORT_PATH),
            "max_abs_difference": float(np.max(abs(source-target))),
            "source_mass": float(source.sum())}


def main():
    exp = Experiment("conditional_feedback_sectors", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "reference, candidate, and conditioned laws normalize with uniform target prefix")
    exp.predict("P2", "sector and full-r candidate laws agree; target TV is reported, not preselected")
    exp.must_fail("C1", "feedback omission/sign control differs for a nonzero Fourier label")
    SETUP_WORK.clear()
    up.SETUP.clear()
    started = time.perf_counter()
    report = {"status": "FAIL", "fixture": {
        "N": N, "a": BASE, "b": BLOCK, "r": PERIOD, "width": WIDTH,
        "prefix_bits": PREFIX_BITS, "low_width": LOW_WIDTH,
        "initial_angle": INITIAL_ANGLE, "post_angle_4": POST_ANGLE_4,
        "post_angle_5": POST_ANGLE_5,
    }}
    p1 = p2 = c1 = False
    counters = {"calls": 0, "work_units": 0}
    try:
        report["preflight"] = preflight()
        if report["preflight"]["planned_calls"] != EXPECTED_CALLS:
            raise AssertionError("call preflight inconsistency")
        if report["preflight"]["planned_work_units"] > MAX_WORK_UNITS:
            raise MemoryError("planned work units exceed bounded experiment cap")
        target, omission, setup = target_and_omission(counters)
        report["target_law"] = target.tolist()
        report["omission_law"] = omission.tolist()
        report["target_mass"] = float(target.sum())
        report["omission_mass"] = float(omission.sum())
        report["target_prefix_masses"] = [
            float(target[z::H].sum()) for z in range(H)]
        report["omission_vs_target_tv"] = float(np.sum(abs(omission-target))/2)
        report["target_setup"] = setup
        report["source_report_crosscheck"] = source_report_crosscheck(target)

        candidate_sector, sector_rows, sector_stats = candidate_sector_laws(counters)
        candidate_full = candidate_full_laws(counters)
        report["candidate_sector_laws"] = candidate_sector.tolist()
        report["candidate_full_laws"] = candidate_full.tolist()
        report["candidate_sector_rows"] = sector_rows
        report["sector_stats"] = sector_stats

        guard_bytes(H * (1 << LOW_WIDTH) * 8, "target conditional laws")
        target_cond = np.zeros((H, 1 << LOW_WIDTH), dtype=float)
        for z in range(H):
            target_cond[z] = H * target[z::H]
        guard_bytes(H * (1 << LOW_WIDTH) * 8, "candidate joint law")
        candidate_joint = candidate_sector / H
        # Candidate rows are indexed (z,w), corresponding to physical output
        # y=z+H*w; reshape(target) would incorrectly use contiguous y blocks.
        guard_bytes(H * (1 << LOW_WIDTH) * 8, "target joint law")
        target_joint = target_cond / H
        guard_bytes(H * (1 << LOW_WIDTH) * 8, "candidate full joint law")
        candidate_full_joint = candidate_full / H
        report["target_conditional_laws"] = target_cond.tolist()
        report["candidate_joint_law_zw"] = candidate_joint.tolist()
        report["candidate_full_joint_law_zw"] = candidate_full_joint.tolist()
        # Physical output ordering is y=z+H*w, i.e. transpose the (z,w)
        # table before flattening.  Keep the row-major z,w tables above too.
        report["candidate_joint_law"] = (candidate_joint.T.reshape(-1)).tolist()
        report["candidate_full_joint_law"] = (candidate_full_joint.T.reshape(-1)).tolist()
        report["candidate_tv_by_z"] = [
            float(np.sum(abs(candidate_sector[z]-target_cond[z]))/2)
            for z in range(H)]
        report["candidate_sector_full_max_abs"] = float(
            np.max(abs(candidate_sector-candidate_full)))
        report["candidate_joint_tv"] = float(
            np.sum(abs(candidate_joint-target_joint))/2)
        report["candidate_full_joint_tv"] = float(
            np.sum(abs(candidate_full_joint-target_joint))/2)
        report["candidate_tv_classification"] = {
            "z0_to_1e-3": report["candidate_tv_by_z"][0] <= 1e-3,
            "all_to_1e-3": all(v <= 1e-3 for v in report["candidate_tv_by_z"]),
            "z0_to_1e-2": report["candidate_tv_by_z"][0] <= 1e-2,
            "all_to_1e-2": all(v <= 1e-2 for v in report["candidate_tv_by_z"]),
        }
        synthetic = synthetic_feedback_control(counters)
        report["synthetic_feedback_control"] = synthetic
        report["counters"] = dict(counters)
        report["setup_work"] = dict(SETUP_WORK)
        report["up_setup_work"] = dict(up.SETUP)
        report["expected_calls"] = EXPECTED_CALLS
        report["expected_work_units"] = report["preflight"]["planned_work_units"]
        report["call_count_match"] = counters["calls"] == EXPECTED_CALLS
        report["work_units_match"] = counters["work_units"] == report["expected_work_units"]
        report["setup_counts_match"] = (
            SETUP_WORK == {"initial_matvec_calls":2,"omission_matrix_products":2,
                "sector_initial_vector_extractions":1,"sector_shift_calls":320,
                "phase_pair_calls":344,"phase_pair_matrix_entries":60488,
                "synthetic_formula_terms":768}
            and up.SETUP == {"repeated_blocks":4,"modular_phase_values":120,
                "construction_matrix_products":7,"shift_matrices":10})

        p1 = bool(abs(report["target_mass"]-1.) < 3e-10
              and abs(report["omission_mass"]-1.) < 3e-10
              and max(abs(x-.25) for x in report["target_prefix_masses"]) < 3e-10
              and all(abs(row["mass"]-1.) < 3e-10 for row in sector_rows)
              and np.all(np.isfinite(candidate_full))
              and np.min(candidate_full) >= -3e-12)
        p2 = bool(report["candidate_sector_full_max_abs"] < 3e-10
              and len(report["candidate_joint_law"]) == Q
              and abs(sum(report["candidate_joint_law"])-1) < 3e-10
              and report["setup_counts_match"]
              and report["call_count_match"] and report["work_units_match"]
              and report["source_report_crosscheck"].get("available", False)
              and report["source_report_crosscheck"]["max_abs_difference"] < 3e-10)
        c1 = bool(synthetic["correct_formula_max_abs_error"] < 3e-12
              and synthetic["omitted_formula_max_abs_error"] < 3e-12
              and synthetic["wrong_formula_max_abs_error"] < 3e-12
              and abs(synthetic["correct_mass"]-1) < 3e-12
              and abs(synthetic["omitted_mass"]-1) < 3e-12
              and abs(synthetic["wrong_sign_mass"]-1) < 3e-12
              and abs(synthetic["omitted_law"][0]-1) < 3e-12
              and synthetic["omission_tv"] > 1e-3
              and synthetic["wrong_sign_tv"] > 1e-3
              and synthetic["correct_argmax"] == 0
              and synthetic["omitted_argmax"] == 0
              and synthetic["wrong_sign_argmax"] == 0)
        report["checks"] = {"P1": p1, "P2": p2, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and c1 else "FAIL"
        exp.check("P1", p1, "target, omission, and candidate laws normalize; target prefix is uniform")
        exp.check("P2", p2, "sector candidate equals independent full-r candidate and work counters match")
        exp.fail_check("C1", c1, "feedback omission and sign controls differ in the finite Fourier law")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before candidate completion")
        exp.check("P2", False, "exception before candidate completion")
        exp.fail_check("C1", False, "exception before sign control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": report["fixture"],
        "reference": "existing uniform_prefix_output repeated-block matrices + sequential_path",
        "candidate": "PeriodicOrbitCircuit sector shifts and independent full-r sequential_path",
        "feedback": "exp(-2*pi*i*z*2**i/64) on candidate branch one",
        "max_dense_bytes": MAX_DENSE_BYTES, "max_calls": MAX_CALLS,
        "max_work_units": MAX_WORK_UNITS,
        "frozen_source_report": str(UP_REPORT_PATH),
        "no_new_propagator": True, "no_sampler_claim": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
