"""Bounded diagnostic for arbitrary early-exponent unitary transfer.

For the frozen N=13,a=2,b=3,r=12,width=3 background, form the four
pre-G work vectors after the first two controls, normalized by sqrt(4).  For
G_k inserted at that split, compare B=G_k A with A times the best right
unitary.  The optimized ROOT fidelity is not an output-TV optimum.  Main's
follow-up embeds each candidate at the split and runs the existing physical
Circuit/statevec suffix, comparing the target with the independent full-r law.
It does not implement another propagator or a sampler.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  SVD residual, diagonal-only residual, omission overlap, commutator and
      positivity checks obey their exact finite-dimensional identities.
  P2  k=0 is an exact null transfer; actual k=1,2 are recorded as open
      diagnostic cases rather than promoted to a theorem.
  P3  Physical target outputs match the independent full-r law; every
      candidate law normalizes and its TV is bounded by joint trace distance.
      No output-error ordering between candidates is predicted.
  P4  The independent closed Fourier formula agrees with the physical suffix.
      The final parity marginal depends only on work rho and is invariant
      under every early-exponent unitary. A target/background parity gap is
      a lower bound on output TV for this whole restricted replacement class;
      whether the fixed kick produces a nonzero gap is an open diagnostic.
  C1  An embedded nontrivial DFT4/2 control fails diagonal-only transfer while
      succeeding under a general exponent-side unitary.

The initial unnormalized synthetic-control failure is retained in its JSON.
Dense arrays and reference matrices are preflight-capped jointly at 16 MiB.
Matrix products and physical gate-entry updates are instrumented separately;
SVD/eigen dimensions are recorded, not called native FLOPs or bit complexity.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from circuits import Circuit
from experiments.experiment_clean_orbit_gates import (
    Fixture, append, checked_run, controlled_multiplier, mixer,
)
from experiments.experiment_additive_phase_output import (
    BASE, BLOCK, N, PERIOD, WIDTH, additive_phase, full_pairs, full_marginal,
    planned_reference_bytes, REFERENCE_WORK,
)


MAX_DENSE_BYTES = 16 * 1024 * 1024
L = 4
MAX_MATRIX_TERMS = 200_000
MAX_PHYSICAL_UPDATES = 100_000_000
WORK = {"matrix_multiply_calls": 0, "matrix_product_terms": 0,
        "svd_dimensions": [], "eigen_dimensions": [], "parity_inner_product_entries": 0}


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
    path = Path("out") / f"additive_phase_transfer_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def mm(left, right):
    """Count actual matrix-product scalar terms, not all native arithmetic."""
    if left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[0]:
        raise ValueError("matrix product shape")
    terms = left.shape[0]*left.shape[1]*right.shape[1]
    if WORK["matrix_product_terms"] + terms > MAX_MATRIX_TERMS:
        raise MemoryError("matrix-product term cap, before multiplication")
    guard((left.shape[0], right.shape[1]), complex, "matrix product")
    WORK["matrix_multiply_calls"] += 1
    WORK["matrix_product_terms"] += terms
    return left @ right


def preflight():
    # Includes returned metrics/candidates, matrix references and conservative
    # LAPACK workspace; separately add physical checked_run and reference work.
    matrix_bytes = 16*(80*PERIOD**2 + 40*PERIOD*L + 40*L**2)
    physical_bytes = 16*64*(1 << 10)
    reference_bytes = planned_reference_bytes()
    total = matrix_bytes + physical_bytes + reference_bytes
    guard_bytes(total, "aggregate numerical payload before scientific arrays")
    return dict(matrix_bytes=matrix_bytes, physical_bytes=physical_bytes,
                reference_bytes=reference_bytes, total_numeric_bytes=total)


def pre_g_matrix():
    """A's four columns, from the actual first-two-control branch matrices."""
    initial, pairs = full_pairs(0, "interior")
    state = np.zeros((PERIOD, 1), dtype=complex)
    state[0, 0] = 1.
    state = mm(initial, state)
    guard((PERIOD, L), np.complex128, "pre-G column matrix")
    A = np.zeros((PERIOD, L), dtype=complex)
    for bits in range(L):
        vector = state.copy()
        vector = mm(pairs[0][(bits >> 0) & 1], vector)
        vector = mm(pairs[1][(bits >> 1) & 1], vector)
        A[:, bits] = vector[:, 0] / math.sqrt(L)
    return A


def transfer_metrics(A, G, label):
    guard((PERIOD, L), np.complex128, f"{label} B matrix")
    B = mm(G, A)
    C = mm(A.conj().T, B)
    WORK["svd_dimensions"].append(L)
    U, singular, Vh = np.linalg.svd(C)
    V = mm(U, Vh)
    E = V.T
    F = float(np.sum(singular))
    Fdiag = float(np.sum(np.abs(np.diag(C))))
    Fomit = float(abs(np.trace(C)))
    general = mm(A, V)
    # Independently apply E to the exponent-row convention to check transpose.
    exponent_form = mm(E, A.T).T
    residual = float(np.linalg.norm(B - general, "fro")**2)
    exponent_residual = float(np.linalg.norm(B - exponent_form, "fro")**2)
    rho = mm(A, A.conj().T)
    commutator = float(np.linalg.norm(mm(mm(G, rho), G.conj().T) - rho, "fro"))
    WORK["eigen_dimensions"].append(PERIOD)
    rho_eigs, rho_basis = np.linalg.eigh((rho + rho.conj().T) / 2)
    sqrt_rho = mm(rho_basis * np.sqrt(np.maximum(rho_eigs, 0)), rho_basis.conj().T)
    fidelity_matrix = mm(mm(sqrt_rho, G), sqrt_rho)
    WORK["svd_dimensions"].append(PERIOD)
    independent_F = float(np.linalg.svd(fidelity_matrix, compute_uv=False).sum())
    general_unitarity = float(np.linalg.norm(mm(V.conj().T, V) - np.eye(L), "fro"))
    diagonal = np.diag(C)
    phases = np.ones(L, dtype=complex)
    nonzero = abs(diagonal) > 0
    phases[nonzero] = diagonal[nonzero] / abs(diagonal[nonzero])
    diagonal_state = A * phases
    diag_residual = float(np.linalg.norm(B - diagonal_state, "fro")**2)
    omit_phase = np.trace(C)/abs(np.trace(C)) if abs(np.trace(C)) > 0 else 1.
    omit_residual = float(np.linalg.norm(B - omit_phase*A, "fro")**2)
    metrics = {
        "label": label, "singular_values": singular.tolist(), "F": F,
        "Fdiag": Fdiag, "Fomit": Fomit,
        "independent_root_fidelity": independent_F,
        "independent_fidelity_error": abs(F-independent_F),
        "residual_squared": residual,
        "residual_identity_error": abs(residual-(2.-2.*F)),
        "exponent_residual_squared": exponent_residual,
        "exponent_convention_error": float(np.max(abs(general-exponent_form))),
        "diagonal_residual_squared": diag_residual,
        "diagonal_identity_error": abs(diag_residual-(2.-2.*Fdiag)),
        "omission_identity_error": abs(omit_residual-(2.-2.*Fomit)),
        "best_joint_trace_bound": math.sqrt(max(0., 1.-F*F)),
        "commutator_frobenius": commutator,
        "rho_trace": float(np.trace(rho).real),
        "rho_min_eigenvalue": float(np.min(rho_eigs)),
        "rho_hermitian_error": float(np.linalg.norm(rho-rho.conj().T, "fro")),
        "general_unitarity_error": general_unitarity,
        "G_unitarity_error": float(np.linalg.norm(mm(G.conj().T, G)-np.eye(PERIOD), "fro")),
        "input_norm_squared": float(np.linalg.norm(A, "fro")**2),
        "target_norm_squared": float(np.linalg.norm(B, "fro")**2),
        "candidate_norms_squared": [float(np.linalg.norm(v, "fro")**2)
                                    for v in (general, diagonal_state, A)],
        "right_unitary": [[[float(z.real), float(z.imag)] for z in row] for row in V],
    }
    return metrics, {"target": B, "general": general,
                     "diagonal": diagonal_state, "omission": A}


def synthetic_control():
    guard((PERIOD, L), np.complex128, "synthetic DFT matrix")
    A = np.zeros((PERIOD, L), dtype=complex)
    for row in range(L):
        for col in range(L):
            # Normalized Fourier UNITARY, then divide by sqrt(L) once more
            # to obtain a normalized purification (Frobenius norm one).
            A[row, col] = np.exp(2j * math.pi * row * col / L) / L
    G = additive_phase(1)
    metrics, _ = transfer_metrics(A, G, "synthetic_normalized_Fourier_purification")
    return metrics


def physical_suffix():
    f = Fixture(13, 2, 3, 12)
    qc = Circuit(10)
    controlled_multiplier(qc, f, 4, 9)
    append(qc, mixer(f, math.pi/7))
    qc.qft([7, 8, 9], inverse=True)
    expected = 12*len(qc.gates)*(1 << qc.n)
    if expected > MAX_PHYSICAL_UPDATES:
        raise MemoryError("all twelve suffix propagations exceed original update cap")
    return qc, f, expected


def suffix_law(A, qc, f, budget):
    state = np.zeros(1 << qc.n, dtype=complex)
    for high in (0, 1):
        for low in range(L):
            for j in range(PERIOD):
                state[((low+L*high) << 7) | pow(BASE, j, N)] = A[j, low]/math.sqrt(2)
    result, costs = checked_run(qc, state, budget)
    law = np.sum(abs(result.reshape(1 << WIDTH, 1 << 7))**2, axis=1)
    idx = np.arange(result.size)
    dirty = ((idx >> f.nw) & ((1 << (f.np+1))-1)) != 0
    leakage = float(np.sum(abs(result[dirty])**2))
    if (not np.all(np.isfinite(law)) or np.min(law) < 0
            or abs(float(law.sum())-1.) > 3e-10 or leakage > 3e-10):
        raise ArithmeticError("invalid suffix law or dirty ancilla")
    return dict(probabilities=law.tolist(), leakage=leakage, **costs)


def closed_output_formula(A):
    """Finite boundary formula, not another gate/state propagation engine.

    The terminal common W cancels in the norm without commuting with U^4.
    Sum_{even y} p_y=(||A||_F^2+Re Tr(A^dag U^4 A))/2.
    """
    shifted = np.roll(A, 4, axis=0)
    probabilities = []
    for y in range(1 << WIDTH):
        phases = np.exp(-2j*np.pi*y*np.arange(L)/(1 << WIDTH)).reshape(L, 1)
        amplitude = mm(A+(-1)**y*shifted, phases)/4
        probabilities.append(float(np.vdot(amplitude, amplitude).real))
    if WORK["parity_inner_product_entries"] + A.size > 12*PERIOD*L:
        raise MemoryError("parity inner product pre-call cap")
    WORK["parity_inner_product_entries"] += A.size
    even_probability = float((np.linalg.norm(A, "fro")**2+np.vdot(A, shifted).real)/2)
    return np.array(probabilities), even_probability


def main():
    exp = Experiment("additive_phase_transfer", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "explicit residual identities, independent root fidelity and matrix checks hold")
    exp.predict("P2", "k=0 is an exact null transfer")
    exp.predict("P3", "physical target/reference laws and candidate TV bounds agree")
    exp.predict("P4", "closed Fourier law and work-only parity marginal agree; all early-unitary candidates preserve parity")
    exp.must_fail("C1", "synthetic general transfer succeeds but diagonal-only transfer fails")
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "output_comparison": "existing physical suffix"}
    p1 = p2 = p3 = p4 = c1 = False
    try:
        report["preflight"] = preflight()
        WORK.update(matrix_multiply_calls=0, matrix_product_terms=0,
                    svd_dimensions=[], eigen_dimensions=[], parity_inner_product_entries=0)
        REFERENCE_WORK.update(sequential_calls=0, width_dimension_cubed_units=0)
        expected_matrix_terms = 9*PERIOD**2 + 4*(
            2*PERIOD**2*L + 3*PERIOD*L**2 + 2*L**3 + 6*PERIOD**3)
        expected_matrix_terms += 12*(1 << WIDTH)*PERIOD*L
        expected_matrix_calls = 9 + 4*13 + 12*(1 << WIDTH)
        if expected_matrix_terms > MAX_MATRIX_TERMS:
            raise MemoryError("preflight matrix product budget")
        qc, fixture, expected_physical_updates = physical_suffix()
        physical_budget = [0]
        A = pre_g_matrix()
        report["input"] = {"N": N, "a": BASE, "b": BLOCK, "r": PERIOD,
                            "width": WIDTH, "columns": L,
                            "A_frobenius_squared": float(np.linalg.norm(A, "fro")**2)}
        for k in (0, 1, 2):
            metrics, candidates = transfer_metrics(A, additive_phase(k), f"actual_k{k}")
            metrics["classification"] = (
                "numerically_exact" if metrics["residual_squared"] < 3e-12
                else "nonexact_diagnostic")
            metrics["physical"] = {name: suffix_law(state, qc, fixture, physical_budget)
                                   for name, state in candidates.items()}
            for name, state in candidates.items():
                law, even = closed_output_formula(state)
                cell = metrics["physical"][name]
                cell["closed_formula_probabilities"] = law.tolist()
                cell["closed_formula_error"] = float(np.max(abs(law-np.asarray(cell["probabilities"]))))
                cell["even_probability_formula"] = even
                cell["even_probability_measured"] = sum(cell["probabilities"][::2])
                cell["even_formula_error"] = abs(even-cell["even_probability_measured"])
            metrics["early_exponent_only_tv_lower_bound"] = abs(
                metrics["physical"]["target"]["even_probability_formula"]
                - metrics["physical"]["omission"]["even_probability_formula"])
            reference = full_marginal(k)
            target = np.asarray(metrics["physical"]["target"]["probabilities"])
            metrics["target_reference_error"] = float(np.max(abs(target-reference)))
            metrics["reference_probabilities"] = reference.tolist()
            for name in ("general", "diagonal", "omission"):
                row = metrics["physical"][name]
                row["tv_from_target"] = float(np.sum(abs(np.asarray(row["probabilities"])-target))/2)
                overlap = float(abs(np.vdot(candidates["target"], candidates[name])))
                row["joint_trace_bound"] = math.sqrt(max(0., 1.-overlap**2))
                row["joint_overlap"] = overlap
            metrics["state_optimum_is_output_best_among_tested"] = bool(
                metrics["physical"]["general"]["tv_from_target"] <= min(
                    metrics["physical"][name]["tv_from_target"]
                    for name in ("diagonal", "omission")) + 3e-10)
            report["rows"].append(metrics)
        synthetic = synthetic_control()
        report["synthetic"] = synthetic
        actual0 = report["rows"][0]
        p1 = (all(row["F"] + 1e-12 >= row["Fdiag"] >= row["Fomit"] - 1e-12
                   for row in report["rows"] + [synthetic])
              and all(row["rho_min_eigenvalue"] >= -3e-12
                      and abs(row["rho_trace"]-1.) < 3e-12
                      and row["rho_hermitian_error"] < 3e-12
                      and row["general_unitarity_error"] < 3e-12
                      and row["G_unitarity_error"] < 3e-12
                      and 0 <= row["F"] <= 1.+3e-12
                      and abs(row["input_norm_squared"]-1.) < 3e-12
                      and abs(row["target_norm_squared"]-1.) < 3e-12
                      and all(abs(v-1.) < 3e-12 for v in row["candidate_norms_squared"])
                      and row["residual_identity_error"] < 3e-12
                      and row["diagonal_identity_error"] < 3e-12
                      and row["omission_identity_error"] < 3e-12
                      and row["exponent_convention_error"] < 3e-12
                      and row["independent_fidelity_error"] < 2e-7
                      for row in report["rows"] + [synthetic]))
        p2 = (actual0["residual_squared"] < 3e-12
              and actual0["diagonal_residual_squared"] < 3e-12
              and actual0["commutator_frobenius"] < 3e-12
              and max(actual0["physical"][name]["tv_from_target"]
                      for name in ("general", "diagonal", "omission")) < 3e-10)
        c1 = (synthetic["residual_squared"] < 3e-12
              and synthetic["commutator_frobenius"] < 3e-12
              and abs(synthetic["F"]-1.) < 3e-12
              and synthetic["diagonal_residual_squared"] > 1e-8)
        p3 = all(row["target_reference_error"] < 3e-10
                 and len(row["reference_probabilities"]) == (1 << WIDTH)
                 and abs(sum(row["reference_probabilities"])-1.) < 3e-10
                 and all(row["physical"][name]["tv_from_target"]
                         <= row["physical"][name]["joint_trace_bound"]+3e-10
                         for name in ("general", "diagonal", "omission"))
                 for row in report["rows"])
        p4 = all(all(cell["closed_formula_error"] < 3e-10
                     and cell["even_formula_error"] < 3e-10
                     and abs(sum(cell["closed_formula_probabilities"])-1.) < 3e-12
                     for cell in row["physical"].values())
                 and all(abs(row["physical"][name]["even_probability_formula"]
                             -row["physical"]["omission"]["even_probability_formula"]) < 3e-12
                         and row["physical"][name]["tv_from_target"]+3e-10
                             >= row["early_exponent_only_tv_lower_bound"]
                         for name in ("general", "diagonal", "omission"))
                 for row in report["rows"])
        report["matrix_work"] = dict(WORK)
        report["reference_work"] = dict(REFERENCE_WORK)
        report["physical_updates"] = physical_budget[0]
        report["expected_physical_updates"] = expected_physical_updates
        report["expected_matrix_terms"] = expected_matrix_terms
        p1 = (p1 and WORK["matrix_product_terms"] == expected_matrix_terms
              and WORK["matrix_multiply_calls"] == expected_matrix_calls
              and WORK["parity_inner_product_entries"] == 12*PERIOD*L
              and len(WORK["svd_dimensions"]) == 8 and len(WORK["eigen_dimensions"]) == 4)
        p3 = (p3 and physical_budget[0] == expected_physical_updates
              and REFERENCE_WORK["sequential_calls"] == 24)
        report["status"] = "PASS" if p1 and p2 and p3 and p4 and c1 else "FAIL"
        exp.check("P1", p1, "SVD, positivity, unitarity and transfer inequalities hold")
        exp.check("P2", p2, "k=0 exact joint-state and output null")
        exp.check("P3", p3, "complete physical laws, joint TV bounds and work counters")
        exp.check("P4", p4, "closed formula and early-register-only parity invariance")
        exp.fail_check("C1", c1, "synthetic general transfer succeeds while diagonal-only fails")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.check("P4", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "b": BLOCK, "r": PERIOD, "width": WIDTH},
        "max_dense_bytes": MAX_DENSE_BYTES,
        "output_comparison": "normalized candidate inserted into existing physical suffix; independent full-r target",
        "work_scope": "instrumented matrix products, SVD/eigen dimensions, physical gate-entry updates; not native FLOPs/bit cost/RSS",
        "root_fidelity_convention": "F is unsquared; output-TV optimum NOT implied",
        "limits": {"matrix_product_terms": MAX_MATRIX_TERMS,
                   "physical_updates": MAX_PHYSICAL_UPDATES},
        "precision": "float diagnostic; independent square-root route uses 2e-7 tolerance near rank deficiency",
        "no_production_sampler": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
