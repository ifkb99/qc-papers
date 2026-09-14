"""Bounded QFT-only test of a finite-displacement uniform output prefix.

DERIVED BEFORE MEASURING (TODO 36).  Freeze N=61, a=2, r=60, b=3 and six
control bits.  The boundary work vectors are built directly from the closed
schedule

    phi_e = W_final U^(32 h) G_1 W_pre U^(l+16 z) W_init |j=0>,
    e = l + 16 z + 32 h.

Only the known orbit-label shifts, repeated 3-by-3 blocks, one diagonal phase,
and the existing inverse QFT are used.  This is not a physical arithmetic
compiler or a generic propagator.

P1: the d=1 output prefix is uniform and the two high-history work states have
    zero cross-overlap, as predicted by the R=6 support separation.
P2: d=2 and d=3 are recorded without a uniformity prediction; their actual
    prefix laws match the independently contracted high-history Gram matrices.
P3: a pre-terminal R=4 cone certifies d=2; the terminal common W preserves
    every measured Gram matrix. Main derived this refinement independently
    before reading the initial output; the original R=6 record is retained.
C1: replacing the final controlled U^32 shift by identity violates separation
    and must produce deterministic even parity.

The first d output bits mean the integer low-d-bit prefix y mod 2^d under the
repository's swapped/standard inverse-QFT convention.  Gates are QFT-only
    references from a supplied boundary state; no full arithmetic circuit is
    implied.  Every retained state, branch matrix and QFT update count is
    preflight-capped.
"""
from __future__ import annotations

import json
import math
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
WORK_DIM = 1 << ((N - 1).bit_length())  # 64 physical work labels
Q = 1 << WIDTH
MAX_BYTES = 16 * 1024 * 1024
MAX_QFT_UPDATES = 5_000_000
TOL = 3e-10
SUPPORT_TOL = 2e-12
RADIUS = 3 * (BLOCK - 1)
RADIUS_BEFORE_FINAL = 2 * (BLOCK - 1)


def report_path() -> Path:
    path = Path("out") / (
        "uniform_prefix_statevec_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        + ".json"
    )
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


def orbit_labels() -> list[int]:
    values = []
    x = 1
    for _ in range(PERIOD):
        values.append(x)
        x = (x * BASE) % N
    if x != 1 or len(set(values)) != PERIOD:
        raise ValueError("N=61, a=2 does not have the supplied order 60")
    return values


def branch_vectors(orbit: np.ndarray, *, identity_last: bool,
                   retain_pre_final: bool, work: dict):
    """Evaluate the ONE displayed closed schedule, not a gate-action engine.

    Dense repeated matrices and NumPy circular index shifts implement the
    finite formula literally. No new parser or reusable propagator is added.
    Pre-terminal columns are retained in the SAME pass, not re-propagated.
    """
    guard_shape((Q, WORK_DIM), np.complex128, "branch matrix")
    init = np.kron(np.eye(PERIOD//BLOCK), work_block(math.pi/4))
    pre = np.kron(np.eye(PERIOD//BLOCK), work_block(-math.pi/10))
    final = np.kron(np.eye(PERIOD//BLOCK), work_block(math.pi/11))
    work["setup_repeated_matrices"] += 3
    phases = np.exp(2j * math.pi * orbit / N)
    work["setup_phase_values"] += PERIOD
    branches = np.zeros((Q, WORK_DIM), dtype=complex)
    pre_final = np.zeros_like(branches) if retain_pre_final else None
    initial_column = init[:, 0]
    for e in range(Q):
        if work["branch_columns"] >= 128 or work["work_matvec_terms"]+2*PERIOD**2 > 1_000_000:
            raise MemoryError("closed branch formula budget before evaluation")
        l = e & 15
        z = (e >> 4) & 1
        h = (e >> 5) & 1
        shifted_low = np.roll(initial_column, l+16*z)
        vector = np.roll(phases*(pre @ shifted_low), 0 if identity_last else 32*h)
        if pre_final is not None:
            pre_final[e, orbit] = vector
        branches[e, orbit] = final @ vector
        work["branch_columns"] += 1
        work["work_matvec_terms"] += 2*PERIOD**2
        work["shift_entry_updates"] += 2*PERIOD
        work["phase_entries"] += PERIOD
    return branches, pre_final


def full_qft_law(branches: np.ndarray, *, work: dict) -> tuple[np.ndarray, dict]:
    """Embed normalized branch columns and run only the inverse QFT."""
    guard_shape((1 << (WIDTH + (N - 1).bit_length()),), np.complex128,
                "joint pre-QFT state")
    state = np.zeros(Q * WORK_DIM, dtype=complex)
    orbit = work["orbit"]
    for e in range(Q):
        state[e * WORK_DIM + orbit] = branches[e, orbit] / math.sqrt(Q)
    qc = Circuit(WIDTH + (N - 1).bit_length())
    qc.qft(list(range((N - 1).bit_length(), WIDTH + (N - 1).bit_length())),
           inverse=True)
    updates = len(qc.gates) * state.size
    if work["qft_gate_entry_updates"] + updates > MAX_QFT_UPDATES:
        raise MemoryError("QFT gate-entry update cap")
    work["qft_gate_entry_updates"] += updates
    work["qft_gates"] = len(qc.gates)
    result = statevec.run(qc, psi=state)
    matrix = result.reshape(Q, WORK_DIM)
    law = np.sum(abs(matrix) ** 2, axis=1)
    return law, {"gate_count": len(qc.gates), "state_entries": state.size,
                 "gate_entry_updates": updates,
                 "norm": float(np.vdot(result, result).real)}


def prefix_law(law: np.ndarray, d: int) -> np.ndarray:
    H = 1 << d
    return np.array([float(np.sum(law[z::H])) for z in range(H)])


def circular_distance(a: int, b: int) -> int:
    delta = (int(a) - int(b)) % PERIOD
    return min(delta, PERIOD - delta)


def support_audit(branches: np.ndarray, orbit: np.ndarray,
                  radius: int) -> dict:
    violations = []
    max_distance = 0
    support_counts = []
    for e, vector in enumerate(branches):
        active = np.flatnonzero(np.abs(vector[orbit]) > SUPPORT_TOL)
        distances = [circular_distance(int(j), e % PERIOD) for j in active]
        if distances:
            max_distance = max(max_distance, max(distances))
        if any(distance > radius for distance in distances):
            violations.append(e)
        support_counts.append(int(len(active)))
    return {"radius": radius, "max_observed_distance": max_distance,
            "support_counts": support_counts, "violating_branches": violations}


def overlap_audit(branches: np.ndarray, d: int, audit_work: dict) -> tuple[np.ndarray, dict]:
    L = 1 << (WIDTH - d)
    H = 1 << d
    matrices = []
    max_offdiag = 0.
    terms = 0
    for l in range(L):
        if audit_work["gram_terms"]+H*H*WORK_DIM > 150_000:
            raise MemoryError("Gram product terms before allocation/product")
        cols = branches[[l + L * h for h in range(H)]]
        gram = cols.conj() @ cols.T
        terms += H*H*WORK_DIM
        audit_work["gram_terms"] += H*H*WORK_DIM
        matrices.append(gram)
        offdiag = gram - np.diag(np.diag(gram))
        max_offdiag = max(max_offdiag, float(np.max(np.abs(offdiag))))
    return np.asarray(matrices), {"L": L, "H": H,
                                  "matrix_product_terms": terms,
                                  "max_offdiag_abs": max_offdiag,
                                  "max_diag_error": float(np.max(
                                      abs(np.diagonal(np.asarray(matrices), axis1=1, axis2=2)-1)))}


def separation_summary(radius: int) -> dict:
    result = {}
    for d in (1, 2, 3):
        L, H = 1 << (WIDTH - d), 1 << d
        distances = [circular_distance(L * q, 0) for q in range(1, H)]
        result[str(d)] = {"L": L, "H": H,
                          "min_center_distance": min(distances) if distances else None,
                          "strictly_separated": bool(distances and min(distances) > 2*radius)}
    return result


def complex_json(matrix: np.ndarray):
    return [[[float(value.real), float(value.imag)] for value in row]
            for row in matrix]


def counter_view(work: dict) -> dict:
    return {key: (int(value) if isinstance(value, (int, np.integer)) else value)
            for key, value in work.items() if key != "orbit"}


def main():
    exp = Experiment("uniform_prefix_statevec", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "d=1 prefix uniformity follows from separated high histories")
    exp.predict("P2", "prefix laws match high-history Gram contractions; all counters reconcile")
    exp.predict("P3", "terminal common W preserves overlaps; R4 certifies d2")
    exp.must_fail("C1", "identity final shift destroys separation and gives even parity")
    report = {"status": "PASS", "rows": []}
    p1 = p2 = p3 = c1 = False
    try:
        # Simultaneously retained branch/Gram arrays, QFT working buffers,
        # matrix setup and construction temporaries; not Python-object RSS.
        planned = 16*(32*Q*WORK_DIM+16*PERIOD**2)+8*(4*Q*WORK_DIM+4*Q)
        guard_bytes(planned, "aggregate numerical preflight BEFORE arrays")
        expected_qft_updates = 2*90*Q*WORK_DIM
        if expected_qft_updates > MAX_QFT_UPDATES:
            raise MemoryError("two QFT runs exceed preflight updates")
        expected_gram_terms = 2*Q*WORK_DIM*sum(1<<d for d in (1,2,3))
        expected_probability_terms = Q*sum((1<<d)**2+(1<<d) for d in (1,2,3))
        if expected_gram_terms > 150_000 or expected_probability_terms > 10_000:
            raise MemoryError("Gram audit preflight")
        audit_work={"gram_terms":0, "probability_terms":0}
        orbit_list = orbit_labels()
        orbit = np.asarray(orbit_list, dtype=np.int64)
        work = {"orbit": orbit, "shift_entry_updates": 0,
                "work_matvec_terms": 0, "branch_columns": 0, "phase_entries": 0,
                "setup_repeated_matrices": 0, "setup_phase_values": 0,
                "qft_gate_entry_updates": 0, "qft_gates": 0}
        branches, pre_final_branches = branch_vectors(
            orbit, identity_last=False, retain_pre_final=True, work=work)
        norms = np.sum(abs(branches) ** 2, axis=1)
        support = support_audit(branches, orbit, RADIUS)
        support_before_final = support_audit(pre_final_branches, orbit,
                                             RADIUS_BEFORE_FINAL)
        law, qft_cost = full_qft_law(branches, work=work)
        prefix_rows = []
        overlaps = {}
        gram_terms = 0
        for d in (1, 2, 3):
            prefix = prefix_law(law, d)
            grams, info = overlap_audit(branches, d, audit_work)
            before_grams, before_info = overlap_audit(pre_final_branches, d, audit_work)
            gram_terms += info["matrix_product_terms"]+before_info["matrix_product_terms"]
            H,L=1<<d,1<<(WIDTH-d)
            formula=[]
            for value in range(H):
                v=np.exp(-2j*np.pi*value*np.arange(H)/H)
                mass=0.
                for g in grams:
                    if audit_work["probability_terms"]+H*H+H > 10_000:
                        raise MemoryError("Gram probability product cap")
                    mass += float(np.vdot(v, g@v).real)
                    audit_work["probability_terms"] += H*H+H
                formula.append(mass/(L*H*H))
            prefix_rows.append({"d": d, "prefix_probabilities": prefix.tolist(),
                                "max_uniform_error": float(np.max(abs(prefix-1/(1<<d)))),
                                "overlap": info,
                                "closed_gram_probabilities": formula,
                                "gram_formula_error": float(np.max(abs(prefix-np.array(formula)))),
                                "terminal_overlap_error": float(np.max(abs(grams-before_grams))),
                                "overlap_matrices": [complex_json(g) for g in grams]})
            overlaps[str(d)] = info
        baseline_counts=counter_view(work)
        control_branches, _ = branch_vectors(orbit, identity_last=True,
                                             retain_pre_final=False, work=work)
        control_law, control_qft_cost = full_qft_law(control_branches,
                                                     work=work)
        report.update({"fixture": {"N": N, "a": BASE, "r": PERIOD,
                                    "b": BLOCK, "t": WIDTH,
                                    "work_qubits": (N-1).bit_length(),
                                    "total_qubits": WIDTH+(N-1).bit_length()},
                       "predicted_radius": RADIUS,
                       "predicted_pre_terminal_radius": RADIUS_BEFORE_FINAL,
                       "separation_R6": separation_summary(RADIUS),
                       "separation_R4_before_final": separation_summary(RADIUS_BEFORE_FINAL),
                       "orbit_labels": orbit_list,
                       "branch_norm_max_error": float(np.max(abs(norms-1))),
                       "support": support,
                       "support_before_final": support_before_final,
                       "full_output_law": law.tolist(),
                       "prefix_rows": prefix_rows,
                       "control_identity_last_law": control_law.tolist(),
                       "control_even_probability": float(np.sum(control_law[::2])),
                       "qft_cost": qft_cost,
                       "control_qft_cost": control_qft_cost,
                       "work_counts": counter_view(work),
                       "baseline_work_counts": baseline_counts,
                       "gram_matrix_product_terms": gram_terms,
                       "gram_probability_terms": audit_work["probability_terms"],
                       "uniform_full_output_tv": float(np.sum(abs(law-1/Q))/2),
                       "payload_bound_bytes": planned})
        d1 = prefix_rows[0]
        p1 = (report["branch_norm_max_error"] < TOL
              and not support["violating_branches"]
              and not support_before_final["violating_branches"]
              and d1["max_uniform_error"] < TOL
              and d1["overlap"]["max_offdiag_abs"] < TOL
              and d1["overlap"]["max_diag_error"] < TOL
              and abs(float(law.sum())-1.) < TOL
              and qft_cost["norm"] - 1. < TOL
              and qft_cost["norm"] - 1. > -TOL)
        p2 = (all(len(row["prefix_probabilities"]) == 1<<row["d"]
                  and abs(sum(row["prefix_probabilities"])-1.) < TOL
                  and min(row["prefix_probabilities"]) >= -TOL
                  and row["gram_formula_error"] < TOL
                  and row["overlap"]["max_diag_error"] < TOL
                  and len(row["overlap_matrices"]) == (1 << (WIDTH-row["d"]))
                  for row in prefix_rows)
              and all(np.isfinite(law)) and len(law) == Q and min(law) >= -TOL
              and work["branch_columns"] == 128
              and work["work_matvec_terms"] == 128*2*PERIOD**2
              and work["shift_entry_updates"] == 128*2*PERIOD
              and work["phase_entries"] == 128*PERIOD
              and work["setup_repeated_matrices"] == 6 and work["setup_phase_values"] == 120
              and work["qft_gate_entry_updates"] == expected_qft_updates
              and gram_terms == audit_work["gram_terms"] == expected_gram_terms
              and audit_work["probability_terms"] == expected_probability_terms)
        p3 = (all(row["terminal_overlap_error"] < TOL for row in prefix_rows)
              and prefix_rows[1]["max_uniform_error"] < TOL
              and prefix_rows[1]["overlap"]["max_offdiag_abs"] < TOL
              and report["separation_R4_before_final"]["2"]["strictly_separated"]
              and not report["separation_R6"]["2"]["strictly_separated"])
        c1 = (float(np.sum(control_law[::2])) > 1.-TOL
              and abs(float(control_law.sum())-1.) < TOL)
        exp.check("P1", p1, "d=1 uniform prefix and separated high histories")
        exp.check("P2", p2, "all prefix Gram identities and explicit cost counters agree")
        exp.check("P3", p3, "common terminal W cancels and sharper support certifies d2")
        exp.fail_check("C1", c1, "identity last shift gives deterministic even parity")
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK, "t": WIDTH},
        "support_radius": RADIUS,
        "max_dense_bytes": MAX_BYTES,
        "max_qft_gate_entry_updates": MAX_QFT_UPDATES,
        "payload_scope": "closed boundary branches plus QFT state/output; no arithmetic compiler",
        "reference_scope": "existing Circuit/statevec inverse QFT only",
        "float_diagnostic": True,
        "no_sampler": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
