"""Selected full-output reference for an earlier pointwise phase (TODO 39).

The frozen indexed fixture is N=61, a=2, r=60, b=3, t=8, s=7, with
W0=work_block(pi/4), W1=work_block(-pi/10), and
g(j)=exp(2*pi*i*(2**j mod 61)/61).  The second copy of g is inserted after
v=0, 2, 3, or 7 low ascending controls.  This is a selected full-output
reference, not the C78 sampler and not a physical arithmetic compiler.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The four complete output laws from Circuit/statevec agree with an
      independent NumPy FFT of the same indexed branch columns.  The first
      output bit is uniform for every selected insertion.
  P2  The complete laws are normalized and the literal branch columns have
      unit norm; the raw 256-entry laws are retained for every insertion.
  C1  Omitting the exponent inverse-QFT, and C2 decohering exponent histories,
      must fail to reproduce the coherent output law whenever its nonuniform
      structure is visible.  These controls test the output observable rather
      than only coefficient variation.

The branch formula is assembled literally with existing 60-by-60 repeated
work matrices, cyclic vector rolls, and pointwise modular phases.  Work labels
are indexed basis labels 0..59 embedded in a 64-state work register; they are
not physical residues 2**j mod 61.  Circuit/statevec supplies the only QFT
reference.  Numeric guards include retained branch matrices, full joint-state
temporaries, repeated matrices, QFT gate-entry work, and modular phase terms.
This is a float64 diagnostic without a finite-bit certificate or timing claim.
"""
from __future__ import annotations

import json
import math
import os
import platform
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from circuits import Circuit
from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_uniform_prefix_output import _repeated_block as repeated_block
from lab import Experiment
from lab.fourier_sampling import unit_phase
import statevec


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT = 8, 7
Q, L, HIGH = 1 << WIDTH, 1 << SPLIT, 1 << (WIDTH - SPLIT)
WORK_BITS, WORK_DIM = 6, 1 << 6
INSERTIONS = (0, 2, 3, 7)
BYTE_CAP = 16 * 1024 * 1024
QFT_UPDATE_CAP = 12_000_000
# The four selected rows perform one dense 60-by-60 W1 matvec per branch
# column: 4*256*3600 plus the initial W0 matvec.  Keep this separate from the
# QFT gate-entry cap; it is a bounded arithmetic counter, not a FLOP claim.
TERM_CAP = 5_000_000
TOL = 3e-10


def report_path() -> Path:
    path = Path("out") / (
        "earlier_phase_qft_"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
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


def phase_at(index: int, counters: dict, key: str) -> complex:
    index = int(index)
    if not 0 <= index < PERIOD:
        raise ValueError("phase index outside indexed orbit")
    add_count(counters, key, 1, label=key)
    add_count(counters, "modular_power_queries", 1, label="modular power")
    return unit_phase(pow(BASE, index, N), N)


def tv(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.sum(np.abs(first - second)) / 2.)


def law_audit(law: np.ndarray) -> dict:
    law = np.asarray(law, dtype=float)
    minimum = float(np.min(law))
    total = float(np.sum(law))
    return {
        "shape": list(law.shape),
        "shape_ok": bool(law.shape == (Q,)),
        "finite": bool(np.all(np.isfinite(law))),
        "minimum": minimum,
        "nonnegative": bool(np.all(np.isfinite(law)) and minimum >= -TOL),
        "sum": total,
        "normalized": bool(np.all(np.isfinite(law)) and abs(total - 1.) < TOL),
    }


def preflight() -> dict:
    """Reserve all retained numeric buffers and named work before allocation."""
    full_entries = Q * WORK_DIM
    branch_entries = Q * PERIOD
    components = {
        # Four selected branch matrices are budgeted even though processing is
        # serial; this keeps the reserve conservative if rows are retained.
        "selected_branch_matrices": len(INSERTIONS) * branch_entries * 16,
        "full_joint_state_temporaries_12": 12 * full_entries * 16,
        "fft_and_qft_working_states": 2 * full_entries * 16,
        "repeated_60x60_buffers": 4 * PERIOD * PERIOD * 16,
        "work_vector_temporaries": 16 * PERIOD * 16,
        "raw_laws_and_controls": 16 * Q * 8,
        "scalar_and_container_reserve": 512 * 1024,
    }
    payload = guard_bytes(sum(components.values()), "aggregate numerical preflight")
    qft = Circuit(WORK_BITS + WIDTH)
    qft.qft(list(range(WORK_BITS, WORK_BITS + WIDTH)), inverse=True)
    gate_count = len(qft.gates)
    qft_updates = len(INSERTIONS) * gate_count * full_entries
    if qft_updates > QFT_UPDATE_CAP:
        raise MemoryError("QFT gate-entry update preflight exceeds cap")
    branch_columns = len(INSERTIONS) * Q
    matvec_terms = PERIOD * PERIOD + branch_columns * PERIOD * PERIOD
    phase_terms = len(INSERTIONS) * Q * 2 * PERIOD
    roll_terms = len(INSERTIONS) * Q * 3 * PERIOD
    if matvec_terms > TERM_CAP or phase_terms > TERM_CAP or roll_terms > TERM_CAP:
        raise MemoryError("literal branch work preflight exceeds cap")
    return {
        "planned_numeric_payload_bytes": payload,
        "payload_components": components,
        "qft_gate_count": gate_count,
        "full_joint_entries": full_entries,
        "branch_entries": branch_entries,
        "qft_gate_entry_updates": qft_updates,
        "qft_update_cap": QFT_UPDATE_CAP,
        "branch_column_count": branch_columns,
        "dense_matvec_term_preflight": matvec_terms,
        "phase_query_preflight": phase_terms,
        "modular_power_preflight": phase_terms,
        "vector_roll_entry_preflight": roll_terms,
        "fft_call_preflight": len(INSERTIONS),
        "fft_input_entry_preflight": len(INSERTIONS) * full_entries,
        "embedding_entry_preflight": branch_columns * PERIOD,
        "max_scalar_terms": TERM_CAP,
        "max_numeric_payload_bytes": BYTE_CAP,
    }


def literal_column(exponent: int, insertion: int, initial: np.ndarray,
                   middle: np.ndarray, counters: dict) -> np.ndarray:
    """Build one indexed work column using the frozen schedule literally."""
    high, low = divmod(int(exponent), L)
    prefix = low % (1 << insertion)
    remainder = low - prefix
    add_count(counters, "vector_roll_entry_updates", PERIOD, label="prefix roll")
    vector = np.roll(initial, prefix)
    for index in range(PERIOD):
        vector[index] *= phase_at(index, counters, "early_phase_queries")
    add_count(counters, "vector_roll_entry_updates", PERIOD, label="remainder roll")
    vector = np.roll(vector, remainder)
    add_count(counters, "dense_matvec_terms", PERIOD * PERIOD,
              label="repeated W1 matvec")
    vector = middle @ vector
    for index in range(PERIOD):
        vector[index] *= phase_at(index, counters, "late_phase_queries")
    add_count(counters, "vector_roll_entry_updates", PERIOD, label="late roll")
    vector = np.roll(vector, L * high)
    return vector


def main() -> None:
    exp = Experiment("earlier_phase_qft", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "statevec inverse-QFT and independent FFT laws agree")
    exp.predict("P2", "all selected full laws and columns are valid; first bit is uniform")
    exp.must_fail("C1", "omitting the exponent inverse-QFT changes the output law")
    exp.must_fail("C2", "decohering exponent histories changes coherent output")
    started = time.perf_counter()
    counters = {
        "branch_columns": 0, "dense_matvec_terms": 0,
        "vector_roll_entry_updates": 0, "early_phase_queries": 0,
        "late_phase_queries": 0, "modular_power_queries": 0,
        "qft_calls": 0, "qft_gate_entry_updates": 0,
        "fft_calls": 0, "fft_transform_entries": 0, "state_embedding_entries": 0,
    }
    report = {"status": "FAIL", "rows": [], "predictions": {
        "first_output_bit": "uniform for v=0,2,3,7",
        "full_laws": "statevec and independent FFT retained",
    }}
    p1 = p2 = c1 = c2 = False
    try:
        report["fixture"] = {
            "N": N, "base": BASE, "period": PERIOD, "block_size": BLOCK,
            "width": WIDTH, "split": SPLIT, "Q": Q, "L": L, "H": HIGH,
            "work_dimension": WORK_DIM, "insertions": list(INSERTIONS),
            "work_basis": "indexed labels 0..59 embedded in 64 labels",
        }
        report["preflight"] = preflight()
        W0 = work_block(math.pi / 4)
        W1 = work_block(-math.pi / 10)
        if (W0.shape != (BLOCK, BLOCK) or W1.shape != (BLOCK, BLOCK)
                or not np.allclose(W0.conj().T @ W0, np.eye(BLOCK), atol=1e-12, rtol=0)
                or not np.allclose(W1.conj().T @ W1, np.eye(BLOCK), atol=1e-12, rtol=0)):
            raise AssertionError("frozen work blocks are not unitary")
        guard_shape((PERIOD,), label="initial work vector")
        basis0 = np.zeros(PERIOD, dtype=complex)
        basis0[0] = 1.
        add_count(counters, "dense_matvec_terms", PERIOD * PERIOD,
                  label="repeated W0 matvec")
        initial_matrix = repeated_block(W0)
        middle_matrix = repeated_block(W1)
        initial = initial_matrix @ basis0
        qft = Circuit(WORK_BITS + WIDTH)
        qft.qft(list(range(WORK_BITS, WORK_BITS + WIDTH)), inverse=True)
        gate_count = len(qft.gates)
        full_entries = Q * WORK_DIM
        full_laws = {}
        fft_laws = {}
        no_qft_laws = {}
        dephased_laws = {}
        row_summaries = []
        for insertion in INSERTIONS:
            guard_shape((Q, PERIOD), label=f"v={insertion} branch matrix")
            branches = np.zeros((Q, PERIOD), dtype=complex)
            for exponent in range(Q):
                add_count(counters, "branch_columns", 1, len(INSERTIONS) * Q,
                          "branch columns")
                branches[exponent] = literal_column(
                    exponent, insertion, initial, middle_matrix, counters)
            column_norm_error = float(np.max(
                np.abs(np.sum(np.abs(branches) ** 2, axis=1) - 1.)))
            guard_shape((Q, WORK_DIM), label=f"v={insertion} full joint state")
            joint = np.zeros((Q, WORK_DIM), dtype=complex)
            add_count(counters, "state_embedding_entries", Q * PERIOD,
                      report["preflight"]["embedding_entry_preflight"],
                      "indexed work embedding")
            joint[:, :PERIOD] = branches / math.sqrt(Q)
            no_qft = np.sum(np.abs(joint) ** 2, axis=1)
            # Explicit history dephasing has a uniform output law: each fixed
            # exponent basis state Fourier-transforms to a uniform exponent
            # distribution after tracing work.
            dephased = np.full(Q, 1. / Q, dtype=float)
            add_count(counters, "fft_calls", 1, len(INSERTIONS))
            add_count(counters, "fft_transform_entries", full_entries,
                      report["preflight"]["fft_input_entry_preflight"],
                      "independent FFT")
            fft_state = np.fft.fft(joint, axis=0) / math.sqrt(Q)
            fft_law = np.sum(np.abs(fft_state) ** 2, axis=1)
            qft_updates = gate_count * full_entries
            if counters["qft_gate_entry_updates"] + qft_updates > QFT_UPDATE_CAP:
                raise MemoryError("QFT update cap before statevec call")
            counters["qft_gate_entry_updates"] += qft_updates
            counters["qft_calls"] += 1
            transformed = statevec.run(qft, psi=joint.reshape(-1))
            qft_law = np.sum(np.abs(transformed.reshape(Q, WORK_DIM)) ** 2, axis=1)
            full_laws[str(insertion)] = qft_law.tolist()
            fft_laws[str(insertion)] = fft_law.tolist()
            no_qft_laws[str(insertion)] = no_qft.tolist()
            dephased_laws[str(insertion)] = dephased.tolist()
            first = np.array([qft_law[0::2].sum(), qft_law[1::2].sum()])
            row_summaries.append({
                "insertion": insertion,
                "branch_column_norm_error": column_norm_error,
                "joint_norm": float(np.sum(np.abs(joint) ** 2)),
                "statevec_law_audit": law_audit(qft_law),
                "fft_law_audit": law_audit(fft_law),
                "no_qft_law_audit": law_audit(no_qft),
                "dephased_law_audit": law_audit(dephased),
                "statevec_vs_fft_tv": tv(qft_law, fft_law),
                "statevec_vs_fft_max_error": float(np.max(np.abs(qft_law - fft_law))),
                "first_bit_law": first.tolist(),
                "first_bit_tv": float(np.sum(np.abs(first - .5)) / 2),
                "first_bit_max_error": float(np.max(np.abs(first - .5))),
                "no_qft_tv": tv(qft_law, no_qft),
                "dephased_tv": tv(qft_law, dephased),
            })
            del branches, joint, fft_state, transformed
        report["rows"] = row_summaries
        report["full_laws"] = full_laws
        report["independent_fft_laws"] = fft_laws
        report["no_qft_control_laws"] = no_qft_laws
        report["dephased_history_control_laws"] = dephased_laws
        report["counters"] = dict(counters)
        report["work_blocks"] = {
            "W0_unitarity_error": float(np.linalg.norm(W0.conj().T @ W0 - np.eye(BLOCK))),
            "W1_unitarity_error": float(np.linalg.norm(W1.conj().T @ W1 - np.eye(BLOCK))),
            "repeated_matrix_payload_bytes": int(initial_matrix.nbytes + middle_matrix.nbytes),
        }
        p1 = (all(row["statevec_vs_fft_tv"] < 3e-10 for row in row_summaries)
              and all(row["statevec_vs_fft_max_error"] < 3e-10 for row in row_summaries))
        p2 = (all(row["branch_column_norm_error"] < 3e-10 for row in row_summaries)
              and all(row[audit][key]
                      for row in row_summaries
                      for audit in ("statevec_law_audit", "fft_law_audit",
                                    "no_qft_law_audit", "dephased_law_audit")
                      for key in ("shape_ok", "finite", "normalized", "nonnegative"))
              and all(row["first_bit_max_error"] < TOL for row in row_summaries)
              and counters["qft_calls"] == len(INSERTIONS)
              and counters["qft_gate_entry_updates"] == len(INSERTIONS) * gate_count * full_entries
              and counters["fft_calls"] == len(INSERTIONS)
              and counters["fft_transform_entries"] == report["preflight"]["fft_input_entry_preflight"]
              and counters["state_embedding_entries"] == report["preflight"]["embedding_entry_preflight"]
              and counters["dense_matvec_terms"] == report["preflight"]["dense_matvec_term_preflight"]
              and counters["modular_power_queries"] == report["preflight"]["modular_power_preflight"]
              and counters["vector_roll_entry_updates"] == report["preflight"]["vector_roll_entry_preflight"])
        c1 = max(row["no_qft_tv"] for row in row_summaries) > 1e-6
        c2 = max(row["dephased_tv"] for row in row_summaries) > 1e-6
        report["checks"] = {"P1": p1, "P2": p2, "C1": c1, "C2": c2}
        report["status"] = "PASS" if p1 and p2 and c1 and c2 else "FAIL"
        exp.check("P1", p1, "statevec and independent FFT laws agree")
        exp.check("P2", p2, "columns, laws, first-bit prediction, and QFT work are valid")
        exp.fail_check("C1", c1, "omitting exponent QFT changes the output law")
        exp.fail_check("C2", c2, "history dephasing changes the output law")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
        exp.fail_check("C2", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["environment"] = {
        "python": sys.version,
        "numpy": np.__version__,
        "platform": platform.platform(),
        "openblas_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
    }
    report["caps"] = {
        "numeric_payload_bytes": BYTE_CAP,
        "qft_gate_entry_updates": QFT_UPDATE_CAP,
        "scalar_terms": TERM_CAP,
        "actual_qft_gate_entry_updates": counters["qft_gate_entry_updates"],
        "actual_scalar_terms": counters["dense_matvec_terms"]
            + counters["vector_roll_entry_updates"]
            + counters["early_phase_queries"]
            + counters["late_phase_queries"]
            + counters["fft_transform_entries"],
    }
    path = report_path()
    # Retain the complete raw report, including every 256-entry law, controls,
    # counters, preflight and environment.  Passing only report["rows"] would
    # silently discard the full-output evidence while still yielding PASS.
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": report.get("fixture"),
        "insertions": list(INSERTIONS),
        "reference": "indexed literal columns + existing Circuit/statevec inverse QFT",
        "independent_crosscheck": "NumPy FFT of same indexed columns",
        "index_basis_not_physical_residue": True,
        "no_new_propagator": True,
        "no_physical_arithmetic_compiler": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
