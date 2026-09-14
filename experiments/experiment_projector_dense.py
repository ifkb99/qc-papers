"""Bounded dense test of the single-Toffoli projector criterion.

Question: for a pure stabilizer state psi and T(a,b;t)=I-2 Pi, with
Pi=(I-Za)(I-Zb)(I-Xt)/8, does the output remain a stabilizer state exactly
when p=<psi|Pi|psi> is different from 1/8?

The inputs are prepared by actual Clifford circuits through the repository's
dense Circuit implementation, including S phases.  The Toffoli under test is
then applied independently by basis-index permutation.  Output stabilizer
status is detected by enumerating every Hermitian Pauli expectation and
counting values whose magnitude is one; the detector does not use p.

PREDICTIONS, WRITTEN BEFORE MEASURING.
P1  Every tested input has p in {0,1/8,1/4,1/2,1}; output stabilizer status
    equals (p != 1/8), including complex Clifford states and n=3..5.
P2  Explicit p=0 and p=1 rows remain stabilizer; p=1 is a global sign under
    the Toffoli when the target is |-> and both controls are |1>.
P3  Multiplying an input by a global -1 changes every output amplitude by -1,
    while preserving p and the independently measured stabilizer count.
C1  A naive rule that every Toffoli preserves every stabilizer input must fail
    on a p=1/8 row (the |++0> input is the smallest such control).

Resource bound: 32x32 dense matrices at most, 33 rows. The detector uses one
Pauli matrix at a time; the projector uses three, identity and temporaries.
Integrated from the independent board draft, with corrected detector gap and
explicit relative-complex coverage. Run from research:
uv run python -m experiments.experiment_projector_dense
"""
from __future__ import annotations

import platform
import sys
from pathlib import Path

import numpy as np

from circuits import Circuit
from lab import Experiment
from pauli import to_matrix


OUT = Path(__file__).resolve().parents[1] / "out"
REPORT = OUT / "projector_dense_main.json"
TOL = 2.0e-9

exp = Experiment("projector_dense", doc=__doc__)
exp.predict("P1", "p has only 0, 1/8, 1/4, 1/2, 1 and stabilizer iff p != 1/8")
exp.predict("P2", "explicit p=0 and p=1 outputs stabilize; p=1 is a sign flip")
exp.predict("P3", "global input minus transports as an output minus exactly")
exp.predict("P4", "corpus includes phases not removable by a global scalar, and detector nonunit expectations are separated from its tolerance")
exp.must_fail("C1", "all-Toffoli-preserve rule fails on p=1/8")


def prepared(n: int, ops: tuple[tuple, ...]) -> np.ndarray:
    """Prepare |0...0> with named repository Clifford gates."""
    qc = Circuit(n)
    for op in ops:
        kind, *args = op
        getattr(qc, kind)(*args)
    zero = np.zeros(1 << n, dtype=complex)
    zero[0] = 1.0
    psi = qc.apply(zero)
    if not np.isclose(np.vdot(psi, psi).real, 1.0, atol=2e-12):
        raise AssertionError("Clifford preparation lost norm")
    return psi


def direct_toffoli(psi: np.ndarray, a: int, b: int, t: int) -> np.ndarray:
    """Apply CCX by direct basis-index permutation, independently of Circuit."""
    n = (psi.size.bit_length() - 1)
    out = np.empty_like(psi)
    am, bm, tm = 1 << a, 1 << b, 1 << t
    for i, amplitude in enumerate(psi):
        j = i ^ tm if (i & am) and (i & bm) else i
        out[j] = amplitude
    if not np.isclose(np.vdot(out, out).real, 1.0, atol=2e-12):
        raise AssertionError("direct Toffoli lost norm")
    return out


def projector_expectation(psi: np.ndarray, n: int, a: int, b: int, t: int) -> float:
    """Evaluate <Pi> from the displayed Pauli projector formula."""
    ident = np.eye(1 << n, dtype=complex)
    za = to_matrix((0, 1 << a), n)
    zb = to_matrix((0, 1 << b), n)
    xt = to_matrix((1 << t, 0), n)
    pi = (ident - za) @ (ident - zb) @ (ident - xt) / 8.0
    value = np.vdot(psi, pi @ psi)
    if abs(value.imag) > 2e-10:
        raise AssertionError(f"projector expectation not real: {value}")
    return float(value.real)


def projector_overlap_check(psi: np.ndarray, n: int, a: int, b: int, t: int) -> float:
    """Second p calculation using control projection and target |-> overlap."""
    control = np.zeros_like(psi)
    mask = (1 << a) | (1 << b)
    for i, value in enumerate(psi):
        if (i & mask) == mask:
            control[i] = value
    anti_x = control.copy()
    for i in range(1 << n):
        anti_x[i] = control[i] - control[i ^ (1 << t)]
    return float(np.vdot(anti_x, anti_x).real / 4.0)


def stabilizer_expectation_count(psi: np.ndarray, n: int) -> tuple[int, float]:
    """Enumerate all 4^n Hermitian Pauli expectations, one matrix at a time."""
    count = 0
    nearest_nonunit = 2.0
    for x in range(1 << n):
        for z in range(1 << n):
            value = np.vdot(psi, to_matrix((x, z), n) @ psi)
            if abs(value.imag) > 3e-8:
                raise AssertionError(f"Hermitian expectation has imaginary part {value}")
            distance = abs(abs(value.real) - 1.0)
            if distance > TOL:
                nearest_nonunit = min(nearest_nonunit, distance)
            if distance <= TOL:
                count += 1
    return count, nearest_nonunit


# Fixed corpus (33 cases): 3, 4, and 5 qubits, with product, entangled, and
# complex S-phase stabilizers.  Extra qubits are spectators or are entangled
# through Clifford CNOTs; no random draw is used, so the rows are reproducible.
CASES: list[tuple[str, int, tuple[tuple, ...]]] = [
    ("n3_zero_p0", 3, ()),
    ("n3_controls11_target0_p1/2", 3, (("x", 0), ("x", 1))),
    ("n3_controls11_targetminus_p1", 3, (("x", 0), ("x", 1), ("h", 2), ("z", 2))),
    ("n3_controls11_targetplus_p0", 3, (("x", 0), ("x", 1), ("h", 2))),
    ("n3_equal_controls_target0_p1/4", 3, (("h", 0), ("cnot", 0, 1))),
    ("n3_pp0_p1/8", 3, (("h", 0), ("h", 1))),
    ("n3_complex_product_p1/8", 3, (("h", 0), ("s", 0), ("h", 1), ("sdg", 1), ("s", 2))),
    ("n3_complex_entangled_p1/8", 3, (("h", 0), ("s", 0), ("cnot", 0, 1), ("h", 2), ("s", 2))),
    ("n3_ghz_targetminus_p1/2", 3, (("h", 0), ("cnot", 0, 1), ("cnot", 0, 2), ("z", 2))),
    ("n3_bell_spectator_p1/8", 3, (("h", 0), ("cnot", 0, 1), ("h", 2))),
    ("n3_y_target_p1/4", 3, (("x", 0), ("x", 1), ("h", 2), ("s", 2))),
    ("n3_signed_pp0_p1/8", 3, (("h", 0), ("z", 0), ("h", 1), ("x", 2))),
    ("n3_signed_equal_p1/4", 3, (("h", 0), ("z", 0), ("cnot", 0, 1))),
    ("n4_zero_p0", 4, ()),
    ("n4_pp0_spectator_p1/8", 4, (("h", 0), ("h", 1), ("s", 3))),
    ("n4_complex_pp0_p1/8", 4, (("h", 0), ("s", 0), ("h", 1), ("sdg", 1), ("s", 2), ("h", 3))),
    ("n4_equal_controls_target0_p1/4", 4, (("h", 0), ("cnot", 0, 1), ("s", 3))),
    ("n4_controls11_targetminus_p1", 4, (("x", 0), ("x", 1), ("h", 2), ("z", 2), ("h", 3), ("s", 3))),
    ("n4_controls11_target0_p1/2", 4, (("x", 0), ("x", 1), ("x", 3))),
    ("n4_ghz_target0_p1/4", 4, (("h", 0), ("cnot", 0, 1), ("cnot", 0, 3))),
    ("n4_entangled_complex_p1/8", 4, (("h", 0), ("cnot", 0, 1), ("s", 1), ("h", 2), ("s", 2), ("h", 3))),
    ("n4_signed_pp0_p1/8", 4, (("h", 0), ("z", 0), ("h", 1), ("z", 1))),
    ("n5_zero_p0", 5, ()),
    ("n5_pp0_spectators_p1/8", 5, (("h", 0), ("h", 1), ("s", 3), ("h", 4))),
    ("n5_complex_product_p1/8", 5, (("h", 0), ("s", 0), ("h", 1), ("sdg", 1), ("s", 2), ("h", 3), ("s", 4))),
    ("n5_equal_controls_target0_p1/4", 5, (("h", 0), ("cnot", 0, 1), ("s", 3), ("h", 4))),
    ("n5_controls11_targetminus_p1", 5, (("x", 0), ("x", 1), ("h", 2), ("z", 2), ("h", 3), ("h", 4))),
    ("n5_controls11_target0_p1/2", 5, (("x", 0), ("x", 1), ("x", 3), ("x", 4))),
    ("n5_ghz_targetminus_p1/2", 5, (("h", 0), ("cnot", 0, 1), ("cnot", 0, 3), ("cnot", 3, 4), ("z", 2))),
    ("n5_complex_entangled_p1/8", 5, (("h", 0), ("s", 0), ("cnot", 0, 1), ("h", 2), ("s", 2), ("cnot", 2, 3), ("h", 4))),
    ("n5_y_target_p1/4", 5, (("x", 0), ("x", 1), ("h", 2), ("s", 2), ("h", 3), ("s", 4))),
    ("n5_signed_pp0_p1/8", 5, (("h", 0), ("z", 0), ("h", 1), ("z", 1), ("s", 3), ("z", 4))),
    ("n5_signed_equal_p1/4", 5, (("h", 0), ("z", 0), ("cnot", 0, 1), ("h", 4), ("s", 4))),
]


rows = []
for name, n, ops in CASES:
    if n > 5:
        raise AssertionError("corpus exceeds n<=5 bound")
    psi = prepared(n, ops)
    pivot = psi[np.argmax(np.abs(psi))]
    relative_complex = bool(np.max(np.abs((psi * np.conj(pivot / abs(pivot))).imag)) > 1e-8)
    psi_minus = -psi
    p = projector_expectation(psi, n, 0, 1, 2)
    p_overlap = projector_overlap_check(psi, n, 0, 1, 2)
    if abs(p - p_overlap) > 5e-10:
        raise AssertionError(f"independent projector calculations disagree: {p} vs {p_overlap}")
    out = direct_toffoli(psi, 0, 1, 2)
    out_minus = direct_toffoli(psi_minus, 0, 1, 2)
    count, nearest_nonunit = stabilizer_expectation_count(out, n)
    count_minus, _ = stabilizer_expectation_count(out_minus, n)
    stabilizer = count == (1 << n)
    p_class = min((0.0, 0.125, 0.25, 0.5, 1.0), key=lambda q: abs(p - q))
    rows.append({
        "name": name,
        "relative_complex": relative_complex,
        "n": n,
        "p": p,
        "p_class": p_class,
        "p_error": abs(p - p_class),
        "pauli_unit_expectation_count": count,
        "pauli_unit_expectation_count_global_minus": count_minus,
        "expected_stabilizer": p_class != 0.125,
        "detected_stabilizer": stabilizer,
        "nearest_nonunit_distance": nearest_nonunit,
        "global_minus_output_error": float(np.max(np.abs(out_minus + out))),
        "p1_toffoli_global_sign_error": float(np.max(np.abs(out + psi))) if p_class == 1.0 else None,
        "input_norm_error": abs(float(np.vdot(psi, psi).real) - 1.0),
    })

expected_p = {0.0, 0.125, 0.25, 0.5, 1.0}
p_classes_ok = all(row["p_error"] <= 5e-9 and row["p_class"] in expected_p for row in rows)
criterion_ok = all(row["detected_stabilizer"] == row["expected_stabilizer"] for row in rows)
sign_ok = all(row["global_minus_output_error"] <= 5e-11 and
              row["pauli_unit_expectation_count_global_minus"] == row["pauli_unit_expectation_count"]
              for row in rows)
p0_rows = [row for row in rows if row["p_class"] == 0.0]
p1_rows = [row for row in rows if row["p_class"] == 1.0]
p0_p1_ok = bool(p0_rows and p1_rows and all(r["detected_stabilizer"] for r in p0_rows + p1_rows)
                 and max(r["p1_toffoli_global_sign_error"] for r in p1_rows) <= 5e-10)
naive_failed = any(row["p_class"] == 0.125 and not row["detected_stabilizer"] for row in rows)

exp.section("P1 criterion and p classes")
exp.check("P1", p_classes_ok and criterion_ok,
          f"rows={len(rows)}, p_classes={sorted(set(row['p_class'] for row in rows))}, "
          f"criterion={criterion_ok}")
exp.section("P2 explicit p=0, p=1")
exp.check("P2", p0_p1_ok,
          f"p0_rows={len(p0_rows)}, p1_rows={len(p1_rows)}, all_stabilizer={p0_p1_ok}")
exp.section("P3 global minus transport")
exp.check("P3", sign_ok,
          f"max global-minus output error={max(row['global_minus_output_error'] for row in rows):.3g}")
exp.section("P4 independent detector coverage")
exp.check("P4", any(row["relative_complex"] for row in rows) and min(row["nearest_nonunit_distance"] for row in rows) > 100*TOL,
          f"relative-complex inputs={sum(row['relative_complex'] for row in rows)}, nonunit distance={min(row['nearest_nonunit_distance'] for row in rows)}")
exp.section("C1 naive all-Toffoli-preserve control")
exp.fail_check("C1", naive_failed,
               "p=1/8 nonstabilizer row observed" if naive_failed else "all rows stabilized")

metadata = {
    "n_cases": len(rows),
    "n_range": [min(row["n"] for row in rows), max(row["n"] for row in rows)],
    "max_dense_matrix_dimension": 32,
    "detector_max_pauli_matrices_live": 1,
    "projector_workspace": "three Pauli matrices, identity and bounded matrix-product temporaries",
    "detector": "all Hermitian Pauli expectations, count |expectation|=1",
    "tolerance": TOL,
    "python": sys.version,
    "numpy": np.__version__,
    "platform": platform.platform(),
}
exp.finish(report_path=REPORT, rows=rows, metadata=metadata)
