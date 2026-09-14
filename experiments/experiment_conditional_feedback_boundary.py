"""Bounded conditional-boundary test for Fourier feedback (TODO 37).

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  On the frozen N=61,a=2,r=60,b=3,t=6 boundary schedule, the exact
      eta_z construction and the complete six-bit inverse-QFT law agree for
      k=0 and k=1 after conditioning on y mod 4.  The eta sign is the one
      derived in TODO 37: exp(-2*pi*i*z*h/4) exp(-2*pi*i*z*l/64).
  P2  The high-history-dephased comparator retains the low-input feedback
      exp(-2*pi*i*z*l/64), and its complete conditional-law TV against the
      target is recorded at 1e-3 and 1e-2 rather than assumed to pass.
  P3  The full k=0 run is the explicit omission-of-G baseline; its conditional
      distance from k=1 is recorded without assuming that it must be large.
  P4  The synthetic normalized flat low-input state has the derived geometric
      Fourier law at z=1, including nonzero w != 0 mass.

C1  The feedback law must differ from the no-feedback z=1 delta-at-w=0 law;
    the no-feedback law itself must retain probability one at w=0.  This
    checks the convention and sign, not a physical no-feedback claim.

The physical rows use the literal closed branch formula

  phi_(l+16h) = W_final U^(32 h1) G_k W_pre U^(16 h0) U^l W_init|0>,
  h=h0+2*h1,

followed by the existing Circuit/statevec inverse QFT.  The eta and candidate
rows use the existing four-qubit QFT on supplied boundary vectors.  This is a
bounded boundary/QFT reference, not a generic propagator or sampler.
"""
from __future__ import annotations

import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from circuits import Circuit
from experiments.experiment_clean_orbit_output import work_block
from lab import Experiment
import statevec


N = 61
BASE = 2
PERIOD = 60
BLOCK = 3
WIDTH = 6
LOW_BITS = 4
LOW_SIZE = 1 << LOW_BITS
HIGH_SIZE = 1 << (WIDTH - LOW_BITS)
WORK_BITS = (N - 1).bit_length()
WORK_DIM = 1 << WORK_BITS
Q = 1 << WIDTH
BYTE_CAP = 16 * 1024 * 1024
MAX_QFT_UPDATES = 2_000_000
MAX_BRANCH_COLUMNS = 2 * Q
MAX_WORK_TERMS = 1_000_000
MAX_TRANSFORM_TERMS = 100_000
TOL = 3e-10


def report_path() -> Path:
    path = Path("out") / (
        "conditional_feedback_boundary_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        + ".json"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(payload: int, label: str) -> int:
    payload = int(payload)
    if payload < 0 or payload > BYTE_CAP:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def guard_shape(shape, dtype=np.complex128, label="array") -> int:
    return guard_bytes(math.prod(int(x) for x in shape)
                       * np.dtype(dtype).itemsize, label)


def orbit_labels() -> list[int]:
    values = []
    x = 1
    for _ in range(PERIOD):
        values.append(x)
        x = (x * BASE) % N
    if x != 1 or len(set(values)) != PERIOD:
        raise ValueError("N=61, a=2 does not have order 60")
    return values


def repeated_block(block: np.ndarray, counters: dict | None = None) -> np.ndarray:
    if block.shape != (BLOCK, BLOCK):
        raise ValueError("invalid work block")
    guard_shape((PERIOD, PERIOD), label="repeated work block")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    if counters is not None:
        counters["setup_repeated_matrix_entries"] += PERIOD * PERIOD
    for m in range(PERIOD // BLOCK):
        start = m * BLOCK
        result[start:start + BLOCK, start:start + BLOCK] = block
        if counters is not None:
            counters["setup_repeated_block_writes"] += BLOCK * BLOCK
    return result


def planned_qft(width: int, work_bits: int) -> tuple[Circuit, int, int]:
    if not 0 <= width <= 12 or work_bits < 0:
        raise ValueError("invalid QFT dimensions")
    qc = Circuit(width + work_bits)
    qc.qft(list(range(work_bits, work_bits + width)), inverse=True)
    entries = 1 << (width + work_bits)
    return qc, len(qc.gates), entries


def qft_law(state_matrix: np.ndarray, *, width: int, work_bits: int,
            counters: dict, label: str) -> np.ndarray:
    """Run an existing inverse QFT on a supplied exponent/work boundary."""
    expected_shape = (1 << width, WORK_DIM)
    if state_matrix.shape != expected_shape:
        raise ValueError(f"{label}: expected {expected_shape}, got {state_matrix.shape}")
    qc, gate_count, entries = planned_qft(width, work_bits)
    updates = gate_count * entries
    if counters["qft_calls"] + 1 > counters["qft_call_cap"]:
        raise MemoryError("QFT call cap before state allocation")
    if counters["qft_gate_entry_updates"] + updates > MAX_QFT_UPDATES:
        raise MemoryError(f"{label}: QFT update cap before state allocation")
    guard_shape((entries,), label=f"{label} state")
    state = np.asarray(state_matrix, dtype=complex).reshape(-1).copy()
    result = statevec.run(qc, psi=state)
    counters["qft_calls"] += 1
    counters["qft_gate_entry_updates"] += updates
    counters.setdefault("qft_rows", []).append({
        "label": label, "width": width, "work_bits": work_bits,
        "gates": gate_count, "state_entries": entries,
        "gate_entry_updates": updates,
    })
    matrix = result.reshape(1 << width, WORK_DIM)
    return np.sum(abs(matrix) ** 2, axis=1)


def build_branches(orbit: np.ndarray, k: int, matrices: tuple[np.ndarray, ...],
                   counters: dict) -> np.ndarray:
    """Evaluate the frozen branch formula with rolls, not a new propagator."""
    init, pre, final = matrices
    guard_shape((Q, WORK_DIM), label=f"k={k} branch matrix")
    branches = np.zeros((Q, WORK_DIM), dtype=complex)
    phases = np.exp(2j * math.pi * k * orbit / N)
    counters["setup_phase_values"] += PERIOD
    if k == 0:
        counters["k0_phase_identity_max_error"] = float(np.max(abs(phases - 1.)))
    initial = init[:, 0]
    for e in range(Q):
        if counters["branch_columns"] + 1 > MAX_BRANCH_COLUMNS:
            raise MemoryError("branch-column cap before evaluation")
        if counters["work_matvec_terms"] + 2 * PERIOD**2 > MAX_WORK_TERMS:
            raise MemoryError("work matvec cap before evaluation")
        l = e % LOW_SIZE
        h0 = (e // LOW_SIZE) % 2
        h1 = e // (2 * LOW_SIZE)
        shifted = np.roll(initial, l + LOW_SIZE * h0)
        vector = pre @ shifted
        vector = phases * vector
        vector = np.roll(vector, 2 * LOW_SIZE * h1)
        vector = final @ vector
        branches[e, orbit] = vector
        counters["branch_columns"] += 1
        counters["work_matvec_terms"] += 2 * PERIOD**2
        counters["shift_entry_updates"] += 2 * PERIOD
        counters["phase_entries"] += PERIOD
    return branches


def eta_boundary(branches: np.ndarray, z: int, counters: dict) -> np.ndarray:
    if branches.shape != (Q, WORK_DIM) or not 0 <= z < HIGH_SIZE:
        raise ValueError("invalid eta inputs")
    guard_shape((LOW_SIZE, WORK_DIM), label="eta boundary")
    eta = np.zeros((LOW_SIZE, WORK_DIM), dtype=complex)
    norm = math.sqrt(Q)
    for l in range(LOW_SIZE):
        feedback = np.exp(-2j * math.pi * z * l / Q)
        if counters["eta_feedback_phase_terms"] + 1 > MAX_TRANSFORM_TERMS:
            raise MemoryError("eta feedback phase cap before loop action")
        counters["eta_feedback_phase_terms"] += 1
        for h in range(HIGH_SIZE):
            if counters["eta_transform_terms"] + WORK_DIM > MAX_TRANSFORM_TERMS:
                raise MemoryError("eta transform cap before loop action")
            eta[l] += (feedback
                       * np.exp(-2j * math.pi * z * h / HIGH_SIZE)
                       * branches[l + LOW_SIZE * h]) / norm
            counters["eta_transform_terms"] += WORK_DIM
            counters["eta_high_phase_terms"] += 1
    return eta


def candidate_boundary(initial: np.ndarray, orbit: np.ndarray, z: int,
                       counters: dict) -> np.ndarray:
    guard_shape((LOW_SIZE, WORK_DIM), label="dephased candidate boundary")
    candidate = np.zeros((LOW_SIZE, WORK_DIM), dtype=complex)
    norm = math.sqrt(LOW_SIZE)
    for l in range(LOW_SIZE):
        if counters["candidate_orbit_entry_terms"] + PERIOD > MAX_TRANSFORM_TERMS:
            raise MemoryError("candidate orbit transform cap before loop action")
        if counters["candidate_roll_entry_updates"] + PERIOD > MAX_TRANSFORM_TERMS:
            raise MemoryError("candidate roll cap before loop action")
        # The closed formula is indexed in orbit order; embed its 60 labels
        # into the 64-state work register exactly as the physical branches do.
        candidate[l, orbit] = (np.exp(-2j * math.pi * z * l / Q)
                               * np.roll(initial, l)) / norm
        counters["candidate_orbit_entry_terms"] += PERIOD
        counters["candidate_roll_entry_updates"] += PERIOD
        counters["candidate_feedback_phase_terms"] += 1
    return candidate


def synthetic_boundary(z: int, feedback: bool, counters: dict | None = None) -> np.ndarray:
    guard_shape((LOW_SIZE, WORK_DIM), label="synthetic boundary")
    result = np.zeros((LOW_SIZE, WORK_DIM), dtype=complex)
    for l in range(LOW_SIZE):
        if counters is not None:
            if counters["synthetic_phase_terms"] + 1 > MAX_TRANSFORM_TERMS:
                raise MemoryError("synthetic phase cap before loop action")
            counters["synthetic_phase_terms"] += 1
        phase = np.exp(-2j * math.pi * z * l / Q) if feedback else 1.
        result[l, 0] = phase / math.sqrt(LOW_SIZE)
    return result


def tv(a: np.ndarray, b: np.ndarray) -> float:
    if a.shape != b.shape:
        raise ValueError("TV shape mismatch")
    return float(np.sum(np.abs(a - b)) / 2)


def geometric_synthetic(z: int, counters: dict | None = None) -> np.ndarray:
    """Derived law for the feedback state, before running Circuit/statevec."""
    out = np.zeros(LOW_SIZE, dtype=float)
    for w in range(LOW_SIZE):
        terms = []
        for l in range(LOW_SIZE):
            if counters is not None:
                if counters["geometric_phase_terms"] + 1 > MAX_TRANSFORM_TERMS:
                    raise MemoryError("geometric phase cap before loop action")
                counters["geometric_phase_terms"] += 1
            terms.append(np.exp(-2j * math.pi * l * (z / Q + w / LOW_SIZE)))
        amp = sum(terms) / LOW_SIZE
        out[w] = abs(amp) ** 2
    return out


def law_audit(law: np.ndarray, expected_size: int) -> dict:
    arr = np.asarray(law)
    finite = bool(np.all(np.isfinite(arr)))
    minimum = float(np.min(arr)) if arr.size else float("nan")
    total = float(np.sum(arr)) if arr.size else 0.0
    return {
        "shape": list(arr.shape),
        "expected_shape": [expected_size],
        "shape_ok": bool(arr.shape == (expected_size,)),
        "finite": finite,
        "minimum": minimum,
        "nonnegative": bool(finite and minimum >= -TOL),
        "sum": total,
        "normalized": bool(finite and abs(total - 1.) < TOL),
    }


def interleave_conditionals(conditionals: dict[int, np.ndarray]) -> np.ndarray:
    result = np.empty(Q, dtype=float)
    for z in range(HIGH_SIZE):
        # The conditional law is p(w|z); the complete law has p(z)=1/H.
        result[z::HIGH_SIZE] = conditionals[z] / HIGH_SIZE
    return result


def main() -> None:
    exp = Experiment("conditional_feedback_boundary", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "eta_z conditional law agrees with the complete QFT law for k=0,1")
    exp.predict("P2", "feedback-aware high-history dephasing is compared at both fixed TV thresholds")
    exp.predict("P3", "the explicit k=0 omission-of-G baseline is reported against k=1")
    exp.predict("P4", "synthetic feedback law agrees with its derived geometric Fourier formula")
    exp.must_fail("C1", "dropping feedback on synthetic z=1 changes the geometric law")

    counters = {
        "branch_columns": 0, "work_matvec_terms": 0,
        "shift_entry_updates": 0, "phase_entries": 0,
        "setup_repeated_matrix_entries": 0,
        "setup_repeated_block_writes": 0, "setup_phase_values": 0,
        "eta_transform_terms": 0, "eta_feedback_phase_terms": 0,
        "eta_high_phase_terms": 0,
        "candidate_orbit_entry_terms": 0,
        "candidate_roll_entry_updates": 0,
        "candidate_feedback_phase_terms": 0,
        "synthetic_phase_terms": 0, "geometric_phase_terms": 0,
        "qft_calls": 0, "qft_gate_entry_updates": 0,
        "qft_call_cap": 20, "qft_rows": [],
        "k0_phase_identity_max_error": 0.,
    }
    report = {"status": "FAIL", "rows": [], "counters": counters}
    p1 = p2 = p3 = p4 = c1 = False
    try:
        # Retained branches for k=0,1; two full-QFT buffers; low-QFT
        # conditional/candidate work; setup matrices and conservative
        # temporary space.  This is checked before any numerical allocation.
        full_entries = 1 << (WIDTH + WORK_BITS)
        low_entries = 1 << (LOW_BITS + WORK_BITS)
        preflight_components = {
            "retained_branch_matrices": 2 * Q * WORK_DIM * 16,
            "full_state_division": full_entries * 16,
            "full_qft_input_copy": full_entries * 16,
            # Conservative reserve for statevec.run's input/output,
            # apply_rotation/apply_pauli temporaries and index arrays.
            "full_statevec_temporaries": 12 * full_entries * 16,
            "low_boundary": LOW_SIZE * WORK_DIM * 16,
            "low_qft_input_copy": low_entries * 16,
            "low_statevec_temporaries": 12 * low_entries * 16,
            "repeated_matrix_allocations": 3 * PERIOD * PERIOD * 16,
            "matrix_construction_temporaries": 5 * PERIOD * PERIOD * 16,
            "law_outputs_and_reports": 32 * Q * 16,
        }
        planned = sum(preflight_components.values())
        guard_bytes(planned, "aggregate numerical preflight")
        if 2 * Q * 2 * PERIOD**2 > MAX_WORK_TERMS:
            raise MemoryError("branch work preflight")
        _, full_gates, planned_full_entries = planned_qft(WIDTH, WORK_BITS)
        _, low_gates, planned_low_entries = planned_qft(LOW_BITS, WORK_BITS)
        if (planned_full_entries != full_entries
                or planned_low_entries != low_entries):
            raise AssertionError("QFT preflight dimensions changed")
        expected_updates = 2 * full_gates * full_entries + 14 * low_gates * low_entries
        if expected_updates > MAX_QFT_UPDATES:
            raise MemoryError("QFT update preflight")
        expected_eta_terms = 2 * HIGH_SIZE * LOW_SIZE * HIGH_SIZE * WORK_DIM
        expected_candidate_terms = HIGH_SIZE * LOW_SIZE * PERIOD
        expected_candidate_rolls = HIGH_SIZE * LOW_SIZE * PERIOD
        expected_synthetic_phases = 2 * LOW_SIZE
        expected_geometric_terms = LOW_SIZE * LOW_SIZE
        if (expected_eta_terms > MAX_TRANSFORM_TERMS
                or expected_candidate_terms > MAX_TRANSFORM_TERMS
                or expected_candidate_rolls > MAX_TRANSFORM_TERMS):
            raise MemoryError("transform preflight")

        orbit_list = orbit_labels()
        orbit = np.asarray(orbit_list, dtype=np.int64)
        init = repeated_block(work_block(math.pi / 4), counters)
        pre = repeated_block(work_block(-math.pi / 10), counters)
        final = repeated_block(work_block(math.pi / 11), counters)
        matrices = (init, pre, final)
        branches_by_k = {}
        for k in (0, 1):
            branches_by_k[k] = build_branches(orbit, k, matrices, counters)
        branch_norm_errors = {
            str(k): float(np.max(np.abs(np.sum(abs(branches_by_k[k]) ** 2, axis=1) - 1.)))
            for k in (0, 1)
        }

        full_laws = {}
        eta_laws = {0: {}, 1: {}}
        candidate_laws = {}
        target_conditionals = {0: {}, 1: {}}
        initial = init[:, 0]
        for k in (0, 1):
            # The complete pre-QFT state has a uniform 1/sqrt(Q) amplitude
            # over its Q exponent branches; eta/candidate boundaries are
            # already normalized by their own formulas.
            full_state = branches_by_k[k] / math.sqrt(Q)
            full_laws[k] = qft_law(full_state, width=WIDTH, work_bits=WORK_BITS,
                                   counters=counters, label=f"full k={k}")
            for z in range(HIGH_SIZE):
                target_mass = float(np.sum(full_laws[k][z::HIGH_SIZE]))
                if not target_mass > 0:
                    raise ArithmeticError(
                        f"exact-zero conditional prefix k={k}, z={z} is out of scope")
                target_conditionals[k][z] = (
                    full_laws[k][z::HIGH_SIZE] / target_mass
                )
                eta = eta_boundary(branches_by_k[k], z, counters)
                eta_laws[k][z] = qft_law(eta, width=LOW_BITS, work_bits=WORK_BITS,
                                         counters=counters, label=f"eta k={k} z={z}")

        for z in range(HIGH_SIZE):
            candidate = candidate_boundary(initial, orbit, z, counters)
            candidate_laws[z] = qft_law(candidate, width=LOW_BITS,
                                        work_bits=WORK_BITS, counters=counters,
                                        label=f"candidate z={z}")

        synthetic_feedback = synthetic_boundary(1, True, counters)
        synthetic_no_feedback = synthetic_boundary(1, False, counters)
        synthetic_feedback_law = qft_law(
            synthetic_feedback, width=LOW_BITS, work_bits=WORK_BITS,
            counters=counters, label="synthetic feedback z=1")
        synthetic_no_feedback_law = qft_law(
            synthetic_no_feedback, width=LOW_BITS, work_bits=WORK_BITS,
            counters=counters, label="synthetic no-feedback z=1")
        synthetic_expected = geometric_synthetic(1, counters)

        eta_errors = []
        target_masses = {str(k): {} for k in (0, 1)}
        comparison_rows = []
        for k in (0, 1):
            for z in range(HIGH_SIZE):
                mass = float(np.sum(full_laws[k][z::HIGH_SIZE]))
                target_masses[str(k)][str(z)] = mass
                if not mass > 0:
                    raise ArithmeticError(
                        f"exact-zero conditional prefix k={k}, z={z} is out of scope")
                eta_error = tv(eta_laws[k][z], target_conditionals[k][z])
                eta_errors.append(eta_error)
                candidate_error = tv(candidate_laws[z], target_conditionals[k][z])
                comparison_rows.append({
                    "k": k, "z": z, "target_mass": mass,
                    "eta_target_tv": eta_error,
                    "candidate_target_tv": candidate_error,
                    "candidate_meets_1e-3": candidate_error <= 1e-3,
                    "candidate_meets_1e-2": candidate_error <= 1e-2,
                    "candidate_vs_omit_g_tv": tv(candidate_laws[z], target_conditionals[0][z]),
                })
        omit_rows = []
        for z in range(HIGH_SIZE):
            omit_rows.append({
                "z": z,
                "omit_g_vs_k1_tv": tv(target_conditionals[0][z], target_conditionals[1][z]),
            })

        # The complete law is interleaved as y=z+H*w.  With the measured
        # prefix masses equal to 1/H, complete-law TV must equal the average
        # conditional TV; this checks the comparison rather than merely
        # recording finite scalar distances.
        interleaving_errors = {}
        candidate_full_comparisons = {}
        for k in (0, 1):
            target_mix = interleave_conditionals(target_conditionals[k])
            interleaving_errors[str(k)] = tv(full_laws[k], target_mix)
            candidate_mix = interleave_conditionals(candidate_laws)
            conditional_tvs = [tv(candidate_laws[z], target_conditionals[k][z])
                               for z in range(HIGH_SIZE)]
            candidate_full_comparisons[str(k)] = {
                "full_law_tv": tv(full_laws[k], candidate_mix),
                "average_conditional_tv": float(np.mean(conditional_tvs)),
                "tv_identity_error": abs(tv(full_laws[k], candidate_mix)
                                          - float(np.mean(conditional_tvs))),
            }
        omit_full_tv = tv(full_laws[0], full_laws[1])
        omit_average_conditional_tv = float(np.mean(
            [row["omit_g_vs_k1_tv"] for row in omit_rows]))
        omit_tv_identity_error = abs(omit_full_tv - omit_average_conditional_tv)

        full_audits = {str(k): law_audit(full_laws[k], Q) for k in (0, 1)}
        target_audits = {
            str(k): {str(z): law_audit(target_conditionals[k][z], LOW_SIZE)
                     for z in range(HIGH_SIZE)}
            for k in (0, 1)
        }
        eta_audits = {
            str(k): {str(z): law_audit(eta_laws[k][z], LOW_SIZE)
                     for z in range(HIGH_SIZE)}
            for k in (0, 1)
        }
        candidate_audits = {
            str(z): law_audit(candidate_laws[z], LOW_SIZE)
            for z in range(HIGH_SIZE)
        }
        expected_costs = {
            "branch_columns": 2 * Q,
            "work_matvec_terms": 2 * Q * 2 * PERIOD**2,
            "shift_entry_updates": 2 * Q * 2 * PERIOD,
            "phase_entries": 2 * Q * PERIOD,
            "setup_repeated_matrix_entries": 3 * PERIOD**2,
            "setup_repeated_block_writes": 3 * (PERIOD // BLOCK) * BLOCK**2,
            "setup_phase_values": 2 * PERIOD,
            "eta_transform_terms": 2 * HIGH_SIZE * LOW_SIZE * HIGH_SIZE * WORK_DIM,
            "eta_feedback_phase_terms": 2 * HIGH_SIZE * LOW_SIZE,
            "eta_high_phase_terms": 2 * HIGH_SIZE * LOW_SIZE * HIGH_SIZE,
            "candidate_orbit_entry_terms": HIGH_SIZE * LOW_SIZE * PERIOD,
            "candidate_roll_entry_updates": HIGH_SIZE * LOW_SIZE * PERIOD,
            "candidate_feedback_phase_terms": HIGH_SIZE * LOW_SIZE,
            "synthetic_phase_terms": 2 * LOW_SIZE,
            "geometric_phase_terms": LOW_SIZE**2,
            "qft_calls": 16,
            "qft_gate_entry_updates": expected_updates,
        }
        cost_reconciliation = {
            key: {"actual": counters[key], "expected": value,
                  "ok": counters[key] == value}
            for key, value in expected_costs.items()
        }
        cost_ok = all(item["ok"] for item in cost_reconciliation.values())
        all_full_valid = all(item["shape_ok"] and item["finite"]
                             and item["nonnegative"] and item["normalized"]
                             for item in full_audits.values())
        all_target_valid = all(item["shape_ok"] and item["finite"]
                               and item["nonnegative"] and item["normalized"]
                               for group in target_audits.values()
                               for item in group.values())
        all_eta_valid = all(item["shape_ok"] and item["finite"]
                            and item["nonnegative"] and item["normalized"]
                            for group in eta_audits.values()
                            for item in group.values())
        all_candidate_valid = all(item["shape_ok"] and item["finite"]
                                  and item["nonnegative"] and item["normalized"]
                                  for item in candidate_audits.values())

        p1 = (max(eta_errors) < TOL
              and all(abs(m - 0.25) < TOL
                      for values in target_masses.values() for m in values.values())
              and all_full_valid and all_target_valid and all_eta_valid
              and all(error < TOL for error in interleaving_errors.values()))
        p2 = (len(comparison_rows) == 8
              and all_candidate_valid
              and all(np.isfinite(row["candidate_target_tv"]) for row in comparison_rows)
              and all("candidate_meets_1e-3" in row and "candidate_meets_1e-2" in row
                      for row in comparison_rows)
              and all(item["tv_identity_error"] < TOL
                      for item in candidate_full_comparisons.values()))
        p3 = (len(omit_rows) == HIGH_SIZE
              and all(np.isfinite(row["omit_g_vs_k1_tv"]) for row in omit_rows)
              and max(branch_norm_errors.values()) < TOL
              and counters["k0_phase_identity_max_error"] < TOL
              and omit_tv_identity_error < TOL and cost_ok
              and counters["qft_calls"] == 16
              and len(counters["qft_rows"]) == 16
              and counters["qft_gate_entry_updates"] == expected_updates)
        synthetic_feedback_audit = law_audit(synthetic_feedback_law, LOW_SIZE)
        synthetic_no_feedback_audit = law_audit(synthetic_no_feedback_law, LOW_SIZE)
        p4 = (tv(synthetic_feedback_law, synthetic_expected) < TOL
              and synthetic_feedback_audit["shape_ok"]
              and synthetic_no_feedback_audit["shape_ok"]
              and synthetic_feedback_audit["normalized"]
              and synthetic_no_feedback_audit["normalized"]
              and synthetic_feedback_audit["nonnegative"]
              and synthetic_no_feedback_audit["nonnegative"])
        feedback_tv = tv(synthetic_feedback_law, synthetic_no_feedback_law)
        c1 = (feedback_tv > 1e-3
              and synthetic_no_feedback_law[0] > 1. - TOL
              and float(np.max(synthetic_no_feedback_law[1:])) < TOL)

        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "Q": Q, "L": LOW_SIZE, "H": HIGH_SIZE,
                        "work_bits": WORK_BITS, "work_dim": WORK_DIM},
            "preflight_payload_bytes": planned,
            "preflight_components": preflight_components,
            "expected_qft_gate_entry_updates": expected_updates,
            "full_laws": {str(k): full_laws[k].tolist() for k in (0, 1)},
            "target_conditionals": {
                str(k): {str(z): target_conditionals[k][z].tolist()
                         for z in range(HIGH_SIZE)}
                for k in (0, 1)
            },
            "eta_conditionals": {
                str(k): {str(z): eta_laws[k][z].tolist()
                         for z in range(HIGH_SIZE)}
                for k in (0, 1)
            },
            "candidate_conditionals": {
                str(z): candidate_laws[z].tolist() for z in range(HIGH_SIZE)
            },
            "law_audits": {
                "full": full_audits, "target": target_audits,
                "eta": eta_audits, "candidate": candidate_audits,
                "synthetic_feedback": synthetic_feedback_audit,
                "synthetic_no_feedback": synthetic_no_feedback_audit,
            },
            "target_masses": target_masses,
            "interleaving_errors": interleaving_errors,
            "candidate_full_comparisons": candidate_full_comparisons,
            "omit_full_tv": omit_full_tv,
            "omit_average_conditional_tv": omit_average_conditional_tv,
            "omit_tv_identity_error": omit_tv_identity_error,
            "branch_norm_errors": branch_norm_errors,
            "cost_reconciliation": cost_reconciliation,
            "eta_errors_max": max(eta_errors),
            "comparisons": comparison_rows,
            "omit_g_rows": omit_rows,
            "synthetic": {
                "z": 1, "feedback_law": synthetic_feedback_law.tolist(),
                "no_feedback_law": synthetic_no_feedback_law.tolist(),
                "derived_geometric_law": synthetic_expected.tolist(),
                "feedback_vs_no_feedback_tv": feedback_tv,
                "feedback_vs_derived_tv": tv(synthetic_feedback_law, synthetic_expected),
            },
            "counters": counters,
        })
        report["status"] = "PASS" if p1 and p2 and p3 and p4 and c1 else "FAIL"
        exp.check("P1", p1, "eta_z agrees with complete QFT conditionals")
        exp.check("P2", p2, "feedback-aware dephased candidate has complete TV rows")
        exp.check("P3", p3, "k=0 omission-of-G baseline is explicitly reported")
        exp.check("P4", p4, "synthetic feedback matches derived geometric law")
        exp.fail_check("C1", c1, "synthetic no-feedback law differs from z=1 feedback law")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.check("P4", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")

    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK, "t": WIDTH},
        "payload_cap_bytes": BYTE_CAP,
        "qft_update_cap": MAX_QFT_UPDATES,
        "scope": "closed boundary formula plus existing Circuit/statevec QFT",
        "float_diagnostic": True,
        "no_generic_propagator": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
