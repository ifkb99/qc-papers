"""Bounded Gauss-kernel and additive physical-phase gate audit.

For the frozen primitive-root prime fixtures, the normalized multiplicative
Fourier coefficients are computed directly and compared with the known Gauss
norm identity.  The forthcoming ``additive_gate(N,k,qubits)`` builder is then
checked column-by-column, including its global phase, and under composition.
The b=3 calculation is kept separate: its subgroup kernel is derived as a
sector projection, not inferred from the b=1 flat-magnitude formula.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The b=1 direct kernels obey the k=0 delta control and the known k!=0
      Gauss norms, while physical builder columns match the additive phase.
  P2  Physical G_1 G_-1 is identity and G_1 G_1 is G_2, including global
      phases; the b=3 projected kernel agrees with its direct formula.
  P3  Every requested direct-sum and retained-array budget remains bounded.
  C1  Erasing complex Fourier-kernel phases breaks inverse/composition.

This is known finite-field mathematics and a tiny gate audit, not a new
Gauss-sum theorem, output sampler, or hardness claim.
"""
from __future__ import annotations

import cmath
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment


MAX_TERMS = 50_000
MAX_BYTES = 16 << 20
EXPECTED_TERMS = 1_496  # b1:872 + b3:192 + composition:432
COMBINED_NUMERIC_BYTES = 16 * (20 * 32**2 + 20 * 16**2)
TOL = 3e-11

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "Gauss norms and physical additive-gate columns agree")
exp.predict("P2", "physical composition and independent b=3 subgroup projection agree")
exp.predict("P3", "direct-term and retained-array budgets stay bounded")
exp.must_fail("C1", "phase-erased Fourier kernels preserve inverse/composition")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"gauss_kernels_{stamp}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def orbit(N, a):
    r = N - 1
    values = []
    x = 1
    for _ in range(r):
        values.append(x)
        x = x * a % N
    if x != 1 or len(set(values)) != r:
        raise ValueError("fixture is not a primitive-root full orbit")
    return values


def gauss_kernel(N, a, k, *, terms):
    values = orbit(N, a)
    r = len(values)
    if terms[0] + r * r > MAX_TERMS:
        raise MemoryError("Gauss direct-term budget exceeded")
    terms[0] += r * r
    coeff = np.zeros(r, complex)
    for delta in range(r):
        coeff[delta] = sum(
            cmath.exp(2j * math.pi * k * values[j] / N)
            * cmath.exp(2j * math.pi * delta * j / r)
            for j in range(r)) / r
    return coeff


def kernel_matrix(coeff):
    r = len(coeff)
    guard((r, r), label="sector kernel matrix")
    # Fourier coefficients act by circular convolution in the sector basis.
    return np.array([[coeff[(beta - alpha) % r]
                      for alpha in range(r)] for beta in range(r)], complex)


def physical_matrix(N, k):
    """Run the supplied additive_gate on every computational basis column."""
    from circuits import Circuit
    import statevec
    from experiments.experiment_additive_phase import additive_gate
    n = max(1, (N - 1).bit_length())
    guard((1 << n, 1 << n), label="physical gate matrix")
    supplied = additive_gate(N, k, list(range(n)))
    global_phase = 1.0 + 0j
    if isinstance(supplied, tuple):
        supplied, global_phase = supplied
    global_phase *= complex(getattr(supplied, "global_phase", 1.0 + 0j))
    if not isinstance(supplied, Circuit):
        raise TypeError("additive_gate must return Circuit (or (Circuit, global_phase))")
    columns = []
    for x in range(1 << n):
        columns.append(statevec.run(supplied, statevec.basis(n, x)))
    matrix = global_phase * np.column_stack(columns)
    if not np.all(np.isfinite(matrix)):
        raise ArithmeticError("nonfinite physical additive-gate matrix")
    return matrix, n, len(supplied.gates), len(supplied.gates) * (1 << n)


def physical_expected(N, k):
    n = max(1, (N - 1).bit_length())
    return np.diag([cmath.exp(2j * math.pi * k * x / N)
                    for x in range(1 << n)])


def b1_audit(terms):
    rows = []
    for N, a in ((7, 3), (13, 2), (17, 3)):
        r = N - 1
        for k in (0, 1):
            coeff = gauss_kernel(N, a, k, terms=terms)
            probs = np.abs(coeff) ** 2
            if k == 0:
                expected = np.array([1.0] + [0.0] * (r - 1))
            else:
                expected = np.array([1 / r**2] + [N / r**2] * (r - 1))
            rows.append({"N": N, "a": a, "k": k,
                         "norm": float(probs.sum()),
                         "max_norm_error": float(abs(probs.sum() - 1)),
                         "max_gauss_probability_error": float(np.max(np.abs(probs - expected))),
                         "zero_probability": float(probs[0]),
                         "nonzero_min": float(np.min(probs[1:])) if r > 1 else 0.0,
                         "nonzero_max": float(np.max(probs[1:])) if r > 1 else 0.0})
            matrix, n, gate_count, column_gate_terms = physical_matrix(N, k)
            expected_matrix = physical_expected(N, k)
            rows[-1].update({"physical_qubits": n,
                             "physical_max_abs_error": float(np.max(np.abs(matrix - expected_matrix))),
                             "physical_unitarity_error": float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(1 << n)))),
                             "physical_gate_entries": gate_count,
                             "physical_column_gate_terms": column_gate_terms,
                             "physical_gate_entry_updates": column_gate_terms * (1 << n),
                             "kernel": coeff.tolist()})
    return rows


def b3_audit(terms):
    # For |alpha,p> = M^-1/2 sum_m exp(-2pi i alpha m/M)|bm+p>,
    # <beta,p|G_k|alpha,p> equals the displayed H(beta,alpha;p).
    from experiments.experiment_coherent_route_sampling import sector_basis
    N, a, b, r, M, k = 13, 2, 3, 12, 4, 1
    values = [pow(a, j, N) for j in range(r)]
    guard((r, r), label="b3 physical diagonal")
    physical = np.diag([cmath.exp(2j * math.pi * k * x / N) for x in values])
    direct = np.zeros((M, M, b), complex)
    max_projection_error = 0.0
    max_off_p = 0.0
    direct_terms = M * M * b * M
    if terms[0] + direct_terms > MAX_TERMS:
        raise MemoryError("b3 subgroup direct-term budget exceeded")
    terms[0] += direct_terms
    for beta in range(M):
        Bbeta = sector_basis(r, b, beta)
        for alpha in range(M):
            Balpha = sector_basis(r, b, alpha)
            projected = Bbeta.conj().T @ physical @ Balpha
            direct[beta, alpha, :] = np.diag(projected)
            max_off_p = max(max_off_p, max(abs(projected[p, q])
                                           for p in range(b) for q in range(b) if p != q))
            for p in range(b):
                direct_formula = sum(
                    cmath.exp(2j * math.pi * k * values[b*m+p] / N)
                    * cmath.exp(2j * math.pi * (beta-alpha) * m / M)
                    for m in range(M)) / M
                max_projection_error = max(max_projection_error,
                                           abs(direct[beta, alpha, p] - direct_formula))
                if abs(direct[beta, alpha, p] - direct_formula) > TOL:
                    raise AssertionError("b3 subgroup kernel formula mismatch")
    if not np.all(np.isfinite(direct)):
        raise ArithmeticError("nonfinite b3 subgroup kernel")
    return {"N": N, "a": a, "b": b, "r": r, "M": M, "k": k,
            "direct_terms": direct_terms, "max_projection_formula_error": float(max_projection_error),
            "off_p_max_abs": float(max_off_p),
            "kernel_norm_sum": float(np.sum(np.abs(direct) ** 2)),
            "kernel": direct.tolist()}


def composition_and_control(terms):
    # Sector-kernel composition is checked independently of the physical
    # matrix columns; the physical columns were checked in b1_audit.
    N, a = 13, 2
    h1 = gauss_kernel(N, a, 1, terms=terms)
    hm = gauss_kernel(N, a, -1, terms=terms)
    h2 = gauss_kernel(N, a, 2, terms=terms)
    K1, Km, K2 = map(kernel_matrix, (h1, hm, h2))
    I = np.eye(N - 1)
    erased = kernel_matrix(np.abs(h1))
    erased_inverse_error = float(np.max(np.abs(erased @ erased - I)))
    physical_one, _, gate1, colterms1 = physical_matrix(N, 1)
    physical_minus_one, _, gateminus, coltermsminus = physical_matrix(N, -1)
    physical_two, _, gate2, colterms2 = physical_matrix(N, 2)
    physical_I = np.eye(1 << max(1, (N - 1).bit_length()))
    return {"G1Gminus1_error": float(np.max(np.abs(K1 @ Km - I))),
            "G1G1G2_error": float(np.max(np.abs(K1 @ K1 - K2))),
            "physical_G1Gminus1_error": float(np.max(np.abs(physical_one @ physical_minus_one - physical_I))),
            "physical_G1G1G2_error": float(np.max(np.abs(physical_one @ physical_one - physical_two))),
            "phase_erased_inverse_error": erased_inverse_error,
            "phase_erased_control_fails": erased_inverse_error > 1e-6,
            "h1": h1.tolist(), "hminus1": hm.tolist(), "h2": h2.tolist(),
            "physical_gate_entries": {"G1": gate1, "Gminus1": gateminus, "G2": gate2},
            "physical_column_gate_terms": {"G1": colterms1, "Gminus1": coltermsminus, "G2": colterms2},
            "physical_gate_entry_updates": (colterms1+coltermsminus+colterms2)*(1 << (N-1).bit_length()),
            "terms": terms[0]}


def main():
    report = {"status": "PASS", "b1": [], "b3": {}, "composition": {}}
    p1 = p2 = p3 = c1 = False
    try:
        if EXPECTED_TERMS > MAX_TERMS:
            raise MemoryError("upfront direct-term budget exceeds cap")
        planned_physical_updates = sum(2*((N-1).bit_length()+1)*(1 << (2*(N-1).bit_length()))
                                       for N in (7,13,17)) + 3*5*16*16
        if planned_physical_updates > 100_000:
            raise MemoryError("upfront physical gate-entry cap")
        terms = [0]
        b1 = b1_audit(terms)
        report["b1"] = b1
        report["b3"] = b3_audit(terms)
        report["composition"] = composition_and_control(terms)
        report["physical_gate_entry_updates"] = sum(row["physical_gate_entry_updates"] for row in b1) + report["composition"]["physical_gate_entry_updates"]
        p1 = (all(row["max_gauss_probability_error"] < TOL
                  and row["max_norm_error"] < TOL
                  and row["physical_max_abs_error"] < TOL
                  for row in b1)
              and all(row["physical_unitarity_error"] < TOL for row in b1))
        p2 = (report["composition"]["G1Gminus1_error"] < TOL
              and report["composition"]["G1G1G2_error"] < TOL
              and report["composition"]["physical_G1Gminus1_error"] < TOL
              and report["composition"]["physical_G1G1G2_error"] < TOL
              and report["b3"]["max_projection_formula_error"] < TOL
              and report["b3"]["off_p_max_abs"] < TOL
              and abs(report["b3"]["kernel_norm_sum"] - 12.0) < TOL)
        p3 = (terms[0] == EXPECTED_TERMS
              and terms[0] <= MAX_TERMS
              and COMBINED_NUMERIC_BYTES <= MAX_BYTES
              and report["physical_gate_entry_updates"] == planned_physical_updates
              and all(row["physical_gate_entries"] == row["physical_qubits"]+1 for row in b1))
        c1 = report["composition"]["phase_erased_control_fails"]
        report["status"] = "PASS" if all((p1, p2, p3, c1)) else "FAIL"
        exp.check("P1", p1, "b=1 Gauss norms and physical columns")
        exp.check("P2", p2, "composition and independently projected b=3 kernel")
        exp.check("P3", p3, "summand/array budgets")
        exp.fail_check("C1", c1, f"phase-erased inverse error={report['composition']['phase_erased_inverse_error']:.6g}")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "max_direct_terms": MAX_TERMS, "array_payload_cap_bytes": MAX_BYTES,
        "expected_direct_terms": EXPECTED_TERMS,
        "combined_numeric_payload_bound_bytes": COMBINED_NUMERIC_BYTES,
        "fixtures": [(7, 3), (13, 2), (17, 3)],
        "b3_fixture": (13, 2, 3, 12), "physical_builder": "experiment_additive_phase.additive_gate",
        "no_generic_propagator": True, "known_gauss_identity": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
