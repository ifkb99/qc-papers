"""Audit two coherent orbit-character reflections with an actual circuit.

DERIVED BEFORE MEASUREMENT.  Use the indexed orbit r=8=b*M with b=2,M=4,
three work bits (p, m_LSB, m_MSB), and ascending controlled addition mod 8.
The cell reflection R maps m to -m mod 4 via CNOT(m_LSB,m_MSB).  For q=0,1,
J_q=D_q R is Hermitian and unitary.  q=0 has a four-term commuting Pauli
expansion.  q=1 has two orthogonal projector groups, each with two commuting
Pauli terms.  Both exponentials can therefore be compiled exactly as products
of existing Circuit.rot gates, including the identity/global phase term.

The q=1 expansion contains noncommuting Pauli terms globally, but splits into
two orthogonal projector groups: J_1=P(c=0) Z_t-P(c=1) Y_t.  Each group's
two Pauli terms commute, so its exponential is compiled exactly in four
existing Circuit.rot gates.  No generic gate compiler or propagator is added.

P1: exact small matrices, grouped Pauli expansions, and the controlled-add
    truth table on every clean-ancilla work/control basis state agree before
    any full-circuit comparison.
P2: the existing Circuit/statevec circuit agrees with an independent full-r
    sequential_path marginal for two distinct coherent reflections, while one
    angle is varied and the original ascending order is retained.
P3: theta=0 recovers the no-reflection route and every compiled variant keeps
    the ancilla clean; the two distinct q labels remain observable with the
    repeated p-qubit mixers retained.

C1: omitting the second reflection must change a selected output probability.
C2: moving the second reflection before its controlled arithmetic insertion
    must change a selected output probability.

This is a bounded abstract indexed-orbit validation, not a physical locality
claim, order discovery, or generic quantum simulator.  Dense allocation is
guarded at 16 MiB before construction; t<=5 and the statevector has only the
three work bits, one increment ancilla and exponent bits.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_coherent_route_circuit
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.semiclassical import sequential_path
from pauli import commutes, to_matrix
import statevec


R, B, M = 8, 2, 4
MAX_T = 5
MAX_BYTES = 16 * 1024 * 1024
WORK_BITS = 3
ANCILLA = 3
EXP_OFFSET = 4


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def validate_width(width: int) -> int:
    if not isinstance(width, (int, np.integer)) or not 1 <= int(width) <= MAX_T:
        raise ValueError(f"width must satisfy 1 <= width <= {MAX_T}")
    return int(width)


def report_path(prefix: str = "coherent_route_circuit") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    return root / f"{prefix}_{stamp}.json"


def cell_reflection_matrix() -> np.ndarray:
    rcell = np.zeros((R, R), dtype=complex)
    for j in range(R):
        p = j & 1
        m = j >> 1
        rcell[p | (((-m) % M) << 1), j] = 1.0
    return rcell


def character_matrix(q: int) -> np.ndarray:
    phases = np.exp(2j * np.pi * q * (np.arange(R) >> 1) / M)
    return np.diag(phases)


def j_matrix(q: int) -> np.ndarray:
    return character_matrix(q) @ cell_reflection_matrix()


def pauli_groups_for_q(q: int) -> list[list[tuple[tuple[int, int], float]]]:
    """Pauli groups for J_q on local cell bits.

    Local bit 0 is m_LSB and local bit 1 is m_MSB.  Terms are returned with
    coefficients +/-1/2, so a coefficient sign changes Circuit.rot's angle.
    Groups are internally commuting.  For q=1 the groups are orthogonal
    projector blocks, so their generated exponentials commute even though
    cross-group Pauli representatives need not commute.
    """
    if q == 0:
        return [[((0, 0), 0.5), ((0, 1), 0.5),
                 ((2, 0), 0.5), ((2, 1), -0.5)]]
    if q == 1:
        # c=m_LSB, t=m_MSB: P(c=0) Z_t and -P(c=1) Y_t.
        return [[((0, 2), 0.5), ((0, 3), 0.5)],
                [((2, 2), -0.5), ((2, 3), 0.5)]]
    if q == 2:
        return [[((0, 0), 0.5), ((0, 1), 0.5),
                 ((2, 0), -0.5), ((2, 1), 0.5)]]
    raise ValueError("only q=0, q=1, and q=2 have verified expansions")


def pauli_terms_for_q(q: int) -> list[tuple[tuple[int, int], float]]:
    """Flatten q's grouped expansion for exact local matrix checks."""
    return [term for group in pauli_groups_for_q(q) for term in group]


def expand_local_terms(q: int) -> list[tuple[tuple[int, int], float]]:
    groups = pauli_groups_for_q(q)
    for group in groups:
        if not all(commutes(a, b) for a, _ in group for b, _ in group):
            raise AssertionError(f"q={q} expansion group must commute")
    local = np.zeros((4, 4), dtype=complex)
    for group in groups:
        for sigma, coeff in group:
            local += coeff * to_matrix(sigma, 2)
    # Build the local D_q R matrix independently of the full orbit p factor.
    rcell = np.zeros((4, 4), dtype=complex)
    for m in range(4):
        rcell[(-m) % 4, m] = 1.0
    d = np.diag(np.exp(2j * np.pi * q * np.arange(4) / 4))
    if np.max(np.abs(local - d @ rcell)) > 2e-12:
        raise AssertionError(f"Pauli expansion mismatch for q={q}")
    return pauli_terms_for_q(q)


def append_coherent_reflection(qc: Circuit, q: int, theta: float) -> None:
    """Append exp(-i theta J_q/2) on cell work bits 1,2.

    The q=1 groups are emitted consecutively.  Terms commute within each
    group; the corresponding projector blocks are orthogonal, so this product
    is the exact exponential despite cross-group Pauli noncommutation.
    """
    for group in pauli_groups_for_q(q):
        for (x, z), coeff in group:
            # Local bit 0 -> global work q=1, local bit 1 -> global work q=2.
            gx = ((x & 1) << 1) | ((x & 2) << 1)
            gz = ((z & 1) << 1) | ((z & 2) << 1)
            qc.rot((gx, gz), coeff * theta)


def coherent_gate_matrix(q: int, theta: float) -> np.ndarray:
    j = j_matrix(q)
    return np.cos(theta / 2) * np.eye(R) - 1j * np.sin(theta / 2) * j


def append_controlled_add_mod8(qc: Circuit, control: int, power_index: int) -> None:
    """Controlled U^(2^power_index), U is increment modulo eight.

    The temporary work ancilla is q=3 and is returned to |0>.  These are only
    existing X/CNOT/Toffoli primitives, not a new state propagator.
    """
    x0, x1, x2, anc = 0, 1, 2, ANCILLA
    if power_index == 0:
        qc.toffoli(control, x0, anc)
        qc.toffoli(anc, x1, x2)
        qc.toffoli(control, x0, anc)
        qc.toffoli(control, x0, x1)
        qc.cnot(control, x0)
    elif power_index == 1:
        qc.toffoli(control, x1, x2)
        qc.cnot(control, x1)
    elif power_index == 2:
        qc.cnot(control, x2)
    # U^8 is identity modulo 8, so power_index >=3 has no gate.


def controlled_add_truth_table(power_index: int) -> tuple[float, float]:
    """Check every clean-ancilla |control,work> input for one controlled power."""
    if not isinstance(power_index, (int, np.integer)) or not 0 <= int(power_index) < MAX_T:
        raise ValueError("truth-table power outside bounded fixture")
    control = 4
    qc = Circuit(5)
    append_controlled_add_mod8(qc, control, int(power_index))
    guard((1 << qc.n,), np.complex128, "controlled-add truth-table statevector")
    max_probability_error = 0.0
    max_physical_sum_error = 0.0
    for c in (0, 1):
        for label in range(R):
            initial = np.zeros(1 << qc.n, dtype=complex)
            initial[label | (c << control)] = 1.0
            # Run the existing Circuit/statevec engine from this basis input.
            output = statevec.run(qc, psi=initial)
            target = (label + c * (1 << int(power_index))) % R
            expected_index = target | (c << control)
            expected_probability = float(abs(output[expected_index]) ** 2)
            max_probability_error = max(max_probability_error,
                                        abs(expected_probability - 1.0))
            max_physical_sum_error = max(max_physical_sum_error,
                                         abs(float(np.vdot(output, output).real) - 1.0))
    return max_probability_error, max_physical_sum_error


def mixer_angle(power_index: int) -> float:
    """Fixed repeated p-block mixer angles for this bounded fixture."""
    return (np.pi / 7, -np.pi / 5, np.pi / 9, -np.pi / 11,
            np.pi / 13)[power_index]


def append_p_mixer(qc: Circuit, power_index: int) -> None:
    """Append the same-branch p-qubit W_i after arithmetic insertion i."""
    if power_index % 2 == 0:
        qc.rx(0, mixer_angle(power_index))
    else:
        qc.rz(0, mixer_angle(power_index))


def circuit_output(theta0: float, theta1: float, *, move_second: bool = False,
                   omit_second: bool = False, width: int = 4) -> tuple[np.ndarray, float]:
    width = validate_width(width)
    n = EXP_OFFSET + width
    guard((1 << n,), np.complex128, "coherent-route statevector")
    qc = Circuit(n)
    exponent = list(range(EXP_OFFSET, EXP_OFFSET + width))
    for q in exponent:
        qc.h(q)
    for i, control in enumerate(exponent):
        append_controlled_add_mod8(qc, control, i)
        append_p_mixer(qc, i)
        if i + 1 == 1:
            append_coherent_reflection(qc, 0, theta0)
        if i + 1 == 2 and not omit_second and not move_second:
            append_coherent_reflection(qc, 1, theta1)
        if i + 1 == 1 and not omit_second and move_second:
            # Deliberately wrong schedule: J_1 is moved before its normal
            # controlled U^2 insertion.  It follows J_0 at this insertion.
            append_coherent_reflection(qc, 1, theta1)
    qc.qft(exponent, inverse=True)
    psi = statevec.run(qc)
    probabilities = np.zeros(1 << width, dtype=float)
    for index, amplitude in enumerate(psi):
        probabilities[index >> EXP_OFFSET] += abs(amplitude) ** 2
    ancilla_leakage = 0.0
    for index, amplitude in enumerate(psi):
        if (index >> ANCILLA) & 1:
            ancilla_leakage += abs(amplitude) ** 2
    return probabilities, ancilla_leakage


def shift_matrix(power: int) -> np.ndarray:
    guard((R, R), np.complex128, "orbit shift")
    result = np.zeros((R, R), dtype=complex)
    for j in range(R):
        result[(j + power) % R, j] = 1.0
    return result


def full_reference(theta0: float, theta1: float, *, move_second: bool = False,
                   omit_second: bool = False, width: int = 4) -> np.ndarray:
    width = validate_width(width)
    v0 = coherent_gate_matrix(0, theta0)
    v1 = coherent_gate_matrix(1, theta1)
    identity = np.eye(R, dtype=complex)
    initial = np.zeros(R, dtype=complex)
    initial[0] = 1.0
    pairs = []
    for i in range(width):
        mixer = Circuit(WORK_BITS)
        append_p_mixer(mixer, i)
        gate = mixer.to_unitary()
        if i + 1 == 1:
            gate = v0 @ gate
        if i + 1 == 2 and not omit_second and not move_second:
            gate = v1 @ gate
        if i + 1 == 1 and not omit_second and move_second:
            gate = v1 @ gate
        pairs.append((gate, gate @ shift_matrix(1 << i)))
    guard((1 << width,), np.float64, "coherent-route marginal")
    return np.array([
        sequential_path(pairs, initial, output=y)["conditional_path_probability"]
        for y in range(1 << width)
    ])


def production_joint_marginal(theta0: float, theta1: float, *, width: int = 4) -> np.ndarray:
    """Marginalize the production coherent-route joint sector/output API.

    The production helper takes b-by-b repeated mixers at insertion s.  The
    circuit fixture's W_i is a one-qubit p block, so each block is built by
    the existing one-qubit Circuit and supplied at insertion i+1.
    """
    width = validate_width(width)
    blocks = {}
    for i in range(width):
        mixer = Circuit(1)
        if i % 2 == 0:
            mixer.rx(0, mixer_angle(i))
        else:
            mixer.rz(0, mixer_angle(i))
        blocks[i + 1] = mixer.to_unitary()
    routed = CoherentReflectionCircuit(
        R, B, width, blocks, {1: (0, theta0), 2: (1, theta1)})
    return np.array([
        sum(routed.joint_probability(sector, y) for sector in range(M))
        for y in range(1 << width)
    ])


def main() -> None:
    exp = Experiment("coherent_route_circuit", doc=__doc__)
    exp.predict("P1", "q=0 and q=1 grouped expansions reproduce exact Hermitian reflections")
    exp.predict("P2", "Circuit/statevec matches full-r sequential_path under two coherent reflections")
    exp.predict("P3", "zero angle and clean ancilla controls behave as derived")
    exp.must_fail("C1", "omitting J_1 changes a visible output probability")
    exp.must_fail("C2", "moving J_1 to the wrong insertion changes a visible output probability")

    started = time.perf_counter()
    rows: list[dict] = []
    exp.section("P1 exact reflection matrices and commuting expansions")
    p1_max = 0.0
    for q in (0, 1):
        terms = expand_local_terms(q)
        j = j_matrix(q)
        matrix_error = float(np.max(np.abs(j - j.conj().T)))
        involution_error = float(np.max(np.abs(j @ j - np.eye(R))))
        for theta in (0.0, np.pi / 7, -np.pi / 3):
            qc = Circuit(WORK_BITS)
            append_coherent_reflection(qc, q, theta)
            compiled = qc.to_unitary()
            expected = coherent_gate_matrix(q, theta)
            error = float(np.max(np.abs(compiled - expected)))
            p1_max = max(p1_max, error, matrix_error, involution_error)
            exp.check("P1", matrix_error < 3e-11 and involution_error < 3e-11
                      and error < 3e-10,
                      f"q={q},theta/pi={theta/np.pi:.6g}: "
                      f"Herm={matrix_error:.2e},J2={involution_error:.2e},gate={error:.2e}")
            rows.append(dict(series="reflection_matrix", q=q,
                             theta=float(theta), theta_over_pi=float(theta / np.pi),
                             term_count=len(terms), hermitian_error=matrix_error,
                             involution_error=involution_error, gate_error=error))
    q1_groups = pauli_groups_for_q(1)
    q1_cross_commuting = all(commutes(a[0], b[0])
                             for a in q1_groups[0] for b in q1_groups[1])
    q1_group_matrices = [sum((coeff * to_matrix(sigma, 2)
                              for sigma, coeff in group), np.zeros((4, 4), complex))
                         for group in q1_groups]
    q1_ab_error = float(np.max(np.abs(q1_group_matrices[0] @ q1_group_matrices[1])))
    q1_commutator_error = float(np.max(np.abs(
        q1_group_matrices[0] @ q1_group_matrices[1]
        - q1_group_matrices[1] @ q1_group_matrices[0])))
    exp.check("P1", q1_ab_error < 3e-11 and q1_commutator_error < 3e-11,
              f"q=1 orthogonal groups: AB={q1_ab_error:.2e}, "
              f"[A,B]={q1_commutator_error:.2e}, "
              f"cross_terms_commuting={q1_cross_commuting}")
    exp.log("q=1 grouped compilation: 4 Pauli terms, "
            f"within_groups=True, AB={q1_ab_error:.2e}, "
            f"[A,B]={q1_commutator_error:.2e}, "
            f"cross_terms_commuting={q1_cross_commuting}; "
            "orthogonal projector groups are emitted separately")

    exp.section("P1 controlled-addition truth table")
    max_add_probability_error = 0.0
    max_add_physical_sum_error = 0.0
    for power_index in range(MAX_T):
        probability_error, physical_sum_error = controlled_add_truth_table(power_index)
        max_add_probability_error = max(max_add_probability_error, probability_error)
        max_add_physical_sum_error = max(max_add_physical_sum_error, physical_sum_error)
        exp.check("P1", probability_error < 3e-11 and physical_sum_error < 3e-11,
                  f"controlled U^{1 << power_index}: basis_probability_error="
                  f"{probability_error:.2e}, physical_sum_error={physical_sum_error:.2e}")
        rows.append(dict(series="controlled_add_truth_table", power_index=power_index,
                         power=1 << power_index,
                         max_expected_basis_probability_error=probability_error,
                         max_physical_sum_error=physical_sum_error,
                         inputs="control=0,1; work=0..7; ancilla=0"))

    exp.section("P2 original-order Circuit/statevec versus full-r reference")
    theta0 = np.pi / 5
    max_circuit_error = 0.0
    for theta1 in (-np.pi / 3, -np.pi / 7, 0.0, np.pi / 8, np.pi / 3):
        reference = full_reference(theta0, theta1)
        circuit, leakage = circuit_output(theta0, theta1)
        error = float(np.max(np.abs(reference - circuit)))
        production = production_joint_marginal(theta0, theta1)
        production_error = float(np.max(np.abs(reference - production)))
        # The marginal already includes BOTH ancilla values.
        physical_sum = float(circuit.sum())
        max_circuit_error = max(max_circuit_error, error)
        exp.check("P2", error < 3e-9 and production_error < 3e-9
                  and leakage < 3e-10
                  and abs(physical_sum - 1.0) < 3e-10
                  and abs(reference.sum() - 1) < 3e-10,
                  f"theta1/pi={theta1/np.pi:.6g}: circuit={error:.2e}, "
                  f"production={production_error:.2e}, ancilla leakage={leakage:.2e}, "
                  f"physical_sum={physical_sum:.12g}, reference_sum={reference.sum():.12g}")
        rows.append(dict(series="two_reflections", theta0=float(theta0),
                         theta1=float(theta1), theta1_over_pi=float(theta1 / np.pi),
                         max_probability_error=error, production_probability_error=production_error,
                         ancilla_leakage=leakage, physical_sum=physical_sum,
                         reference=reference.tolist(), circuit=circuit.tolist()))

    exp.section("P3 zero-angle and distinct-q visibility")
    zero_ref = full_reference(0.0, 0.0)
    zero_circuit, zero_leakage = circuit_output(0.0, 0.0)
    zero_physical_sum = float(zero_circuit.sum())
    ideal_circuit, _ = circuit_output(0.0, 0.0, omit_second=True)
    exp.check("P3", np.max(np.abs(zero_ref - zero_circuit)) < 3e-9
              and zero_leakage < 3e-10 and abs(zero_physical_sum - 1.0) < 3e-10,
              f"zero-angle circuit/reference={np.max(np.abs(zero_ref-zero_circuit)):.2e}, "
              f"leakage={zero_leakage:.2e}, physical_sum={zero_physical_sum:.12g}")
    q_visibility_ref = full_reference(theta0, np.pi / 3)
    q0_only_ref = full_reference(theta0, np.pi / 3, omit_second=True)
    visibility = float(np.max(np.abs(q_visibility_ref - q0_only_ref)))
    exp.check("P3", visibility > 1e-5,
              f"distinct q=0/q=1 reflection visibility={visibility:.6g}")
    rows.append(dict(series="controls_positive", zero_leakage=zero_leakage,
                     q_visibility=visibility,
                     zero_physical_sum=zero_physical_sum,
                     zero_error=float(np.max(np.abs(zero_ref-zero_circuit))),
                     ideal_without_second=ideal_circuit.tolist()))

    exp.section("must-fail controls")
    selected_theta1 = np.pi / 3
    reference = full_reference(theta0, selected_theta1)
    omitted = full_reference(theta0, selected_theta1, omit_second=True)
    moved = full_reference(theta0, selected_theta1, move_second=True)
    omitted_error = float(np.max(np.abs(reference - omitted)))
    moved_error = float(np.max(np.abs(reference - moved)))
    exp.fail_check("C1", omitted_error > 1e-5,
                   f"omitted J_1 probability error={omitted_error:.6g}")
    exp.fail_check("C2", moved_error > 1e-5,
                   f"moved J_1 schedule probability error={moved_error:.6g}")
    rows.append(dict(series="must_fail_controls", omitted_error=omitted_error,
                     moved_error=moved_error))

    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(
                   orbit_period=R, block_size=B, cell_count=M, max_t=MAX_T,
                   max_dense_bytes=MAX_BYTES, max_reflection_gate_error=p1_max,
                   max_circuit_probability_error=max_circuit_error,
                   max_controlled_add_probability_error=max_add_probability_error,
                   max_controlled_add_physical_sum_error=max_add_physical_sum_error,
                   q1_group_product_error=q1_ab_error,
                   q1_group_commutator_error=q1_commutator_error,
                   q1_compilation="orthogonal projector groups; each group has commuting Pauli terms",
                   references="existing Circuit/statevec vs existing full-r sequential_path",
                   assumption="abstract indexed r=8 orbit; controlled mod-8 addition, repeated p-block mixers, and coherent J_q gate synthesis are supplied constructions",
                   elapsed_seconds=time.perf_counter() - started,
               ))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("coherent_route_circuit_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
