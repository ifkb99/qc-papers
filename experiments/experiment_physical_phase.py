"""Physical work-qubit Rz phase amid repeated orbit mixers.

Fix N=7, a=3, orbit [1,3,2,6,4,5], b=3 and t=5.  The background is the
validated repeated schedule W01@s1, W12@s3, W01@s4, with Rx(pi/2) blocks.
At s=2 insert the physical work-qubit Rz(theta) on work bit 1.  Its orbit
phase is

    z[j] = (-1)**bit1(3**j mod 7) = [+,-,-,-,+,+],
    z[j+3] = -z[j].

Thus theta=pi gives K=-i Z, which permutes the two U^3 coarse sectors even
though [K,U^3] != 0; initial coarse dephasing must be restored at that angle.
Intermediate angles can couple sectors and are tested without assuming a
uniform failure.  The physical circuit and the existing full-r sequential_path
contraction are the two references; no generic propagator is introduced.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  Existing Circuit/statevec output agrees with full-r sequential_path for
      theta=0, +/-pi/8, +/-pi/4, +/-pi/2, +/-pi.
  P2  Wrong initial coarse dephasing using the SAME full-r branches fails at
      at least one intermediate angle, but is restored at theta=pi by the
      derived sector swap (and also holds at theta=0).
  P3  A physical phase inserted at the end is output-invisible, and a
      genuinely repeated diagonal orbit phase agrees with the periodic-sector
      helper and its coarse-dephased reference.
  P4  Replacing the physical phase by Rz(theta) on exponent bit 1 after s=2
      while leaving earlier background mixers unchanged fails either in output
      law or full state vector at an intermediate angle.

  C1  Initial coarse dephasing must fail at one intermediate angle.
  C2  Naive exponent-only phase transfer must fail for the background circuit.

This is a finite physical validation, not a hardness, locality, novelty or
generic simulation claim.  State allocations are preflight-capped at 16 MiB;
no amplitude cutoff is used.  Rz phases on invalid work labels are allowed,
but must not leak out of their basis states.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_physical_phase
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
from experiments import experiment_symmetry_kick as kick
from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import sequential_path
from modexp import ModExp
import statevec


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, WIDTH, PERIOD, BLOCK = 7, 3, 5, 6, 3
ORBIT = np.array([1, 3, 2, 6, 4, 5], dtype=np.int64)
BACKGROUND_THETA = np.pi / 2
ANGLES = np.pi * np.array([0., 1/8, -1/8, 1/4, -1/4, 1/2, -1/2, 1., -1.])


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def report_path(prefix: str = "physical_phase") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def complex_payload(matrix: np.ndarray) -> dict:
    return dict(real=np.asarray(matrix).real.tolist(), imag=np.asarray(matrix).imag.tolist())


def orbit_basis(index: int) -> np.ndarray:
    result = np.zeros(PERIOD, dtype=complex)
    result[index] = 1.
    return result


def work_phase_orbit(theta: float) -> np.ndarray:
    bits = np.array([(int(label) >> 1) & 1 for label in ORBIT])
    z = 1 - 2 * bits
    return np.diag(np.exp(-0.5j * theta * z))


def background_full_blocks() -> dict[int, np.ndarray]:
    return {1: kick.repeated_orbit_matrix(kick.block_rotation("W01", BACKGROUND_THETA)),
            3: kick.repeated_orbit_matrix(kick.block_rotation("W12", BACKGROUND_THETA)),
            4: kick.repeated_orbit_matrix(kick.block_rotation("W01", BACKGROUND_THETA))}


def background_sector_blocks() -> dict[int, np.ndarray]:
    return {1: kick.block_rotation("W01", BACKGROUND_THETA),
            3: kick.block_rotation("W12", BACKGROUND_THETA),
            4: kick.block_rotation("W01", BACKGROUND_THETA)}


def full_r_branches(theta: float, *, phase_insertion: int = 2,
                    phase_matrix: np.ndarray | None = None,
                    width: int = WIDTH) -> list[tuple[np.ndarray, np.ndarray]]:
    if phase_insertion not in (2, width):
        raise ValueError("phase insertion must be the intended s=2 or end position")
    blocks = background_full_blocks()
    # An end-only work phase is traced out and therefore omitted from the
    # branch matrices.  At s=2 it is composed as the explicit post-control
    # operation; no coincident background insertion is overwritten silently.
    if phase_insertion < width:
        if phase_insertion in blocks:
            raise ValueError("phase insertion collides with a background defect")
        blocks[phase_insertion] = work_phase_orbit(theta) if phase_matrix is None else phase_matrix
    identity = np.eye(PERIOD, dtype=complex)
    pairs = []
    for i in range(width):
        w = blocks.get(i + 1, identity)
        shift = np.zeros((PERIOD, PERIOD), dtype=complex)
        for j in range(PERIOD):
            shift[(j + (1 << i)) % PERIOD, j] = 1.
        pairs.append((w, w @ shift))
    return pairs


def sequential_marginal(theta: float, *, initial_label: int = 0,
                        phase_insertion: int = 2,
                        phase_matrix: np.ndarray | None = None,
                        width: int = WIDTH) -> np.ndarray:
    pairs = full_r_branches(theta, phase_insertion=phase_insertion,
                            phase_matrix=phase_matrix, width=width)
    result = np.zeros(1 << width, dtype=float)
    initial = orbit_basis(initial_label)
    for y in range(1 << width):
        result[y] = sequential_path(pairs, initial, output=y)["conditional_path_probability"]
    return result


def physical_circuit(theta: float, *, end_phase: bool = False,
                     exponent_transfer: bool = False, width: int = WIDTH) -> tuple[Circuit, ModExp]:
    me = ModExp(N, BASE, n_exp=width)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == 1:
            qc.extend(kick.background_work_gate("W01", BACKGROUND_THETA, me.x))
        elif i + 1 == 2:
            if not end_phase:
                if exponent_transfer:
                    qc.rz(me.exp[1], theta)
                else:
                    qc.rz(me.x[1], theta)
        elif i + 1 == 3:
            qc.extend(kick.background_work_gate("W12", BACKGROUND_THETA, me.x))
        elif i + 1 == 4:
            qc.extend(kick.background_work_gate("W01", BACKGROUND_THETA, me.x))
    if end_phase:
        qc.rz(me.x[1], theta)
    qc.qft(me.exp, inverse=True)
    return qc, me


def physical_output(theta: float, *, end_phase: bool = False,
                    exponent_transfer: bool = False,
                    width: int = WIDTH) -> tuple[np.ndarray, np.ndarray, dict]:
    qc, me = physical_circuit(theta, end_phase=end_phase,
                              exponent_transfer=exponent_transfer, width=width)
    payload = (1 << me.n_qubits) * np.dtype(np.complex128).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"physical state allocation {payload} bytes exceeds 16 MiB")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    psi = statevec.run(qc, psi)
    probabilities = np.sum(np.abs(psi.reshape(1 << width, -1)) ** 2, axis=1)
    return probabilities, psi, dict(qubits=me.n_qubits, vector_payload_bytes=psi.nbytes,
                                    norm_error=abs(float(np.vdot(psi, psi).real) - 1),
                                    gates=len(qc.gates), rotations=qc.n_nonclifford())


def periodic_diag_reference(theta: float) -> tuple[np.ndarray, np.ndarray]:
    block = np.diag(np.exp(-0.5j * theta * np.array([1., -1., 1.])))
    blocks = background_sector_blocks()
    blocks[2] = block
    sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, WIDTH, blocks)
    result = np.zeros(1 << WIDTH, dtype=float)
    for alpha in range(sampler.sectors):
        for y in range(1 << WIDTH):
            result[y] += sampler.forced_joint(alpha, y)["joint_latent_output_probability"]
    full = np.zeros(1 << WIDTH, dtype=float)
    branch_pairs = full_r_branches(theta, phase_matrix=np.kron(np.eye(2), block))
    for y in range(1 << WIDTH):
        full[y] = sequential_path(branch_pairs, orbit_basis(0), output=y)["conditional_path_probability"]
    return result, full


# ---------------------------------------------------------------------------
# Independent capped N=13 extension.  The two-level rotations are synthesized
# only for this fixed four-work-bit fixture: a Gray CNOT path maps each pair to
# one target bit, then a finite Pauli expansion implements the required
# multi-controlled Rx and unmaps it.  This is not a generic gate compiler.
# ---------------------------------------------------------------------------
N13, BASE13, PERIOD13 = 13, 2, 12
ORBIT13 = np.array([1, 2, 4, 8, 3, 6, 12, 11, 9, 5, 10, 7], dtype=np.int64)
EXTENSION_GATE_BUDGET = 30_000


def controlled_rx_fixture(qc: Circuit, target: int, controls: list[int],
                          values: list[int], theta: float) -> None:
    if len(controls) != len(values) or target in controls:
        raise ValueError("invalid fixed controlled-Rx fixture")
    for mask in range(1 << len(controls)):
        zmask = 0
        sign = 1
        for j, control in enumerate(controls):
            if mask & (1 << j):
                zmask |= 1 << control
                sign *= -1 if values[j] else 1
        qc.rot((1 << target, zmask), sign * theta / (1 << len(controls)))


def two_level_rotation_fixture(left: int, right: int, theta: float,
                               work_qubits: list[int]) -> Circuit:
    """Rx(theta) on |left>,|right> only, via a fixed Gray-path expansion."""
    if not (0 <= left < 1 << len(work_qubits) and 0 <= right < 1 << len(work_qubits)):
        raise ValueError("two-level labels outside work register")
    diff = left ^ right
    differing = [bit for bit in range(len(work_qubits)) if diff & (1 << bit)]
    if not differing:
        raise ValueError("two-level endpoints must differ")
    target = differing[0]
    qc = Circuit(max(work_qubits) + 1)
    # CNOT target -> k removes every additional differing bit from the pair.
    for bit in differing[1:]:
        qc.cnot(work_qubits[target], work_qubits[bit])
    transformed = left
    for bit in differing[1:]:
        if (transformed >> target) & 1:
            transformed ^= 1 << bit
    controls = [bit for bit in range(len(work_qubits)) if bit != target]
    values = [(transformed >> bit) & 1 for bit in controls]
    controlled_rx_fixture(qc, work_qubits[target], [work_qubits[b] for b in controls],
                          values, theta)
    for bit in reversed(differing[1:]):
        qc.cnot(work_qubits[target], work_qubits[bit])
    return qc


def n13_repeated_block_gate(kind: str, theta: float,
                            work_qubits: list[int] | None = None) -> Circuit:
    if work_qubits is None:
        work_qubits = [0, 1, 2, 3]
    if len(work_qubits) != 4 or kind not in ("W01", "W12"):
        raise ValueError("N=13 fixture requires four work qubits")
    pairs = ((1, 2), (8, 3), (12, 11), (5, 10)) if kind == "W01" else \
            ((2, 4), (3, 6), (11, 9), (10, 7))
    qc = Circuit(max(work_qubits) + 1)
    for left, right in pairs:
        qc.extend(two_level_rotation_fixture(left, right, theta, work_qubits))
    return qc


def n13_physical_circuit(theta: float) -> tuple[Circuit, ModExp]:
    me = ModExp(N13, BASE13, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE13, 1 << i, N13)))
        if i + 1 == 1:
            qc.extend(n13_repeated_block_gate("W01", BACKGROUND_THETA, me.x))
        elif i + 1 == 2:
            qc.rz(me.x[1], theta)
        elif i + 1 == 3:
            qc.extend(n13_repeated_block_gate("W12", BACKGROUND_THETA, me.x))
        elif i + 1 == 4:
            qc.extend(n13_repeated_block_gate("W01", BACKGROUND_THETA, me.x))
    qc.qft(me.exp, inverse=True)
    return qc, me


def n13_physical_output(theta: float) -> tuple[np.ndarray, np.ndarray, dict]:
    qc, me = n13_physical_circuit(theta)
    payload = (1 << me.n_qubits) * np.dtype(np.complex128).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"N=13 physical state allocation {payload} bytes exceeds 16 MiB")
    if len(qc.gates) > EXTENSION_GATE_BUDGET:
        raise MemoryError(f"N=13 fixture gate count {len(qc.gates)} exceeds budget")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    psi = statevec.run(qc, psi)
    return (np.sum(np.abs(psi.reshape(1 << WIDTH, -1)) ** 2, axis=1), psi,
            dict(qubits=me.n_qubits, vector_payload_bytes=psi.nbytes,
                 gates=len(qc.gates), rotations=qc.n_nonclifford(),
                 norm_error=abs(float(np.vdot(psi, psi).real) - 1)))


def generic_full_r_marginal(period: int, orbit: np.ndarray, blocks: dict[int, np.ndarray],
                            width: int = WIDTH) -> np.ndarray:
    if not 1 <= period <= 18 or not 4 <= width <= 7 or orbit.shape != (period,):
        raise ValueError("full-r extension exceeds finite fixture caps")
    guard((period, period), np.complex128, "extension full-r matrix")
    guard((1 << width,), np.float64, "extension output table")
    pairs = []
    identity = np.eye(period, dtype=complex)
    for i in range(width):
        w = blocks.get(i + 1, identity)
        shift = np.zeros((period, period), dtype=complex)
        for j in range(period):
            shift[(j + (1 << i)) % period, j] = 1.
        pairs.append((w, w @ shift))
    initial = np.zeros(period, dtype=complex)
    initial[0] = 1.
    result = np.zeros(1 << width, dtype=float)
    for y in range(1 << width):
        result[y] = sequential_path(pairs, initial, output=y)["conditional_path_probability"]
    return result


def n13_grouped_reference(theta: float) -> tuple[np.ndarray, np.ndarray]:
    w01_3 = kick.block_rotation("W01", BACKGROUND_THETA)
    w12_3 = kick.block_rotation("W12", BACKGROUND_THETA)
    w01_6 = np.kron(np.eye(2), w01_3)
    w12_6 = np.kron(np.eye(2), w12_3)
    z = 1 - 2 * ((ORBIT13 >> 1) & 1)
    phase_6 = np.diag(np.exp(-0.5j * theta * z[:6]))
    grouped = PeriodicOrbitCircuit(
        PERIOD13, 6, WIDTH, {1: w01_6, 2: phase_6, 3: w12_6, 4: w01_6})
    grouped_output = np.zeros(1 << WIDTH, dtype=float)
    for alpha in range(grouped.sectors):
        for y in range(1 << WIDTH):
            grouped_output[y] += grouped.forced_joint(alpha, y)["joint_latent_output_probability"]
    full_blocks = {
        1: np.kron(np.eye(4), w01_3),
        2: np.diag(np.exp(-0.5j * theta * z)),
        3: np.kron(np.eye(4), w12_3),
        4: np.kron(np.eye(4), w01_3),
    }
    return grouped_output, generic_full_r_marginal(PERIOD13, ORBIT13, full_blocks)


def validate_n13_work_gates(exp: Experiment) -> dict:
    rows = []
    for kind in ("W01", "W12"):
        local = n13_repeated_block_gate(kind, BACKGROUND_THETA)
        actual = local.to_unitary()
        expected = np.eye(16, dtype=complex)
        pairs = ((1, 2), (8, 3), (12, 11), (5, 10)) if kind == "W01" else \
                ((2, 4), (3, 6), (11, 9), (10, 7))
        c, s = np.cos(BACKGROUND_THETA / 2), np.sin(BACKGROUND_THETA / 2)
        for left, right in pairs:
            expected[left, left] = expected[right, right] = c
            expected[left, right] = expected[right, left] = -1j * s
        matrix_error = float(np.max(np.abs(actual - expected)))
        statevec_error, invalid_error = 0., 0.
        orbit_set = set(int(x) for x in ORBIT13)
        for label in range(16):
            basis = statevec.basis(4, label)
            evolved = statevec.run(local, basis)
            statevec_error = max(statevec_error,
                                 float(np.max(np.abs(evolved - actual @ basis))))
            if label not in orbit_set:
                invalid_error = max(invalid_error,
                                    float(np.max(np.abs(evolved - basis))))
        exp.check("P5", matrix_error < 3e-12 and statevec_error < 3e-12
                  and invalid_error < 3e-12,
                  f"N=13 {kind}: matrix={matrix_error:.2e}, statevec={statevec_error:.2e}, "
                  f"invalid={invalid_error:.2e}")
        rows.append(dict(kind=kind, matrix_error=matrix_error,
                         statevec_error=statevec_error, invalid_error=invalid_error,
                         gates=len(local.gates)))
    rz = Circuit(4).rz(1, np.pi / 4)
    rz_unitary = rz.to_unitary()
    rz_expected = np.diag(np.exp(-.5j*(np.pi/4)*(1-2*((np.arange(16) >> 1) & 1))))
    rz_matrix_error = float(np.max(np.abs(rz_unitary-rz_expected)))
    rz_leakage = 0.
    for label in range(16):
        evolved = statevec.run(rz, statevec.basis(4, label))
        rz_leakage = max(rz_leakage,
                         float(np.sum(np.abs(evolved[np.arange(16) != label]) ** 2)))
    exp.check("P5", rz_leakage < 3e-12 and rz_matrix_error < 3e-12,
              f"N=13 Rz work bit1 leakage={rz_leakage:.2e}, matrix={rz_matrix_error:.2e}")
    return dict(rows=rows, rz_leakage=rz_leakage, rz_matrix_error=rz_matrix_error,
                orbit=ORBIT13.tolist())


def preflight_dimensions(exp: Experiment) -> dict:
    rows = []
    for width in range(4, 8):
        me = ModExp(N, BASE, n_exp=width)
        payload = (1 << me.n_qubits) * np.dtype(np.complex128).itemsize
        rows.append(dict(width=width, qubits=me.n_qubits, state_payload_bytes=payload))
    ok = all(row["state_payload_bytes"] <= MAX_DENSE_BYTES for row in rows)
    ok = ok and PERIOD <= 18 and PERIOD13 <= 18 and ORBIT.size == PERIOD
    exp.check("P5", ok, f"preflight widths4..7={rows}, periods=6,12 within r<=18")
    return dict(rows=rows, max_period=PERIOD13, dense_budget_bytes=MAX_DENSE_BYTES)


def validate_rz_matrix(exp: Experiment) -> dict:
    local = Circuit(3).rz(1, float(np.pi / 2))
    actual = local.to_unitary()
    expected = np.diag(np.exp(-0.5j * (np.pi / 2) * (1 - 2*((np.arange(8) >> 1) & 1))))
    unitary_error = float(np.max(np.abs(actual.conj().T @ actual - np.eye(8))))
    matrix_error = float(np.max(np.abs(actual - expected)))
    statevec_error, leakage = 0., 0.
    for label in range(8):
        basis = statevec.basis(3, label)
        evolved = statevec.run(local, basis)
        statevec_error = max(statevec_error,
                             float(np.max(np.abs(evolved - actual @ basis))))
        leakage = max(leakage, float(np.sum(np.abs(evolved[np.arange(8) != label]) ** 2)))
    exp.check("P4", unitary_error < 2e-14 and matrix_error < 2e-14
              and statevec_error < 2e-14 and leakage < 2e-14,
              f"Rz work bit1: unitary={unitary_error:.2e}, matrix={matrix_error:.2e}, "
              f"statevec={statevec_error:.2e}, leakage={leakage:.2e}")
    return dict(unitary_error=unitary_error, matrix_error=matrix_error,
                statevec_error=statevec_error, leakage=leakage,
                actual=complex_payload(actual), expected=complex_payload(expected))


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(np.asarray(a) - np.asarray(b))) / 2)


def main() -> None:
    exp = Experiment("physical_phase", doc=__doc__)
    exp.predict("P1", "physical statevec agrees with full-r sequential_path at every angle")
    exp.predict("P2", "coarse dephasing fails intermediate but is restored at theta=pi")
    exp.predict("P3", "end phase and repeated diagonal orbit phase are valid controls")
    exp.predict("P4", "physical Rz has no work-label leakage and matches its matrix")
    exp.predict("P5", "capped N=13 physical fixture agrees with full-r and grouped-b6 references")
    exp.must_fail("C1", "wrong initial coarse dephasing at an intermediate angle")
    exp.must_fail("C2", "naive exponent-only phase transfer with earlier mixers unchanged")
    start = time.perf_counter()
    rows = []
    preflight = preflight_dimensions(exp)
    setup = validate_rz_matrix(exp)
    max_dephasing_tv = 0.
    max_transfer_tv = max_transfer_state = 0.
    base_output, base_state, base_info = physical_output(0.)
    end_output, _, end_info = physical_output(np.pi / 2, end_phase=True)
    end_error = float(np.max(np.abs(end_output - base_output)))
    exp.check("P3", end_error < 3e-10,
              f"end phase output change={end_error:.2e}, norm={end_info['norm_error']:.2e}")
    for theta in ANGLES:
        theta = float(theta)
        physical, state, physical_info = physical_output(theta)
        sequential = sequential_marginal(theta)
        seq_error = float(np.max(np.abs(physical - sequential)))
        exp.check("P1", seq_error < 3e-10 and physical_info["norm_error"] < 2e-12,
                  f"theta/pi={theta/np.pi:g}: statevec/sequential={seq_error:.2e}, "
                  f"norm={physical_info['norm_error']:.2e}")
        wrong = (sequential_marginal(theta, initial_label=0)
                 + sequential_marginal(theta, initial_label=3)) / 2
        dephasing_tv = tv(wrong, physical)
        max_dephasing_tv = max(max_dephasing_tv, dephasing_tv)
        if abs(abs(theta) - np.pi) < 1e-14:
            exp.check("P2", dephasing_tv < 3e-10,
                      f"theta/pi={theta/np.pi:g} sector-swap dephasing TV={dephasing_tv:.2e}")
        elif abs(theta) < 1e-14:
            exp.check("P2", dephasing_tv < 3e-10,
                      f"theta=0 dephasing TV={dephasing_tv:.2e}")
        transfer, transfer_state, transfer_info = physical_output(theta, exponent_transfer=True)
        transfer_tv = tv(transfer, physical)
        transfer_state_error = float(np.linalg.norm(transfer_state - state))
        max_transfer_tv = max(max_transfer_tv, transfer_tv)
        max_transfer_state = max(max_transfer_state, transfer_state_error)
        if abs(theta) > 1e-14:
            exp.check("P4", transfer_tv > 1e-7 or transfer_state_error > 1e-7,
                      f"theta/pi={theta/np.pi:g}: naive transfer TV={transfer_tv:.3e}, "
                      f"state L2={transfer_state_error:.3e}")
        periodic, periodic_full = periodic_diag_reference(theta)
        periodic_error = float(np.max(np.abs(periodic - periodic_full)))
        periodic_block = np.kron(
            np.eye(2), np.diag(np.exp(-0.5j * theta * np.array([1., -1., 1.]))))
        periodic_dephased = (sequential_marginal(theta, phase_matrix=periodic_block,
                                                 initial_label=0)
                             + sequential_marginal(theta, phase_matrix=periodic_block,
                                                   initial_label=3)) / 2
        dephased_error = float(np.max(np.abs(periodic_full - periodic_dephased)))
        exp.check("P3", periodic_error < 3e-10 and dephased_error < 3e-10,
                  f"theta/pi={theta/np.pi:g}: periodic helper/full-r={periodic_error:.2e}, "
                  f"dephased={dephased_error:.2e}")
        rows.append(dict(theta=theta, theta_over_pi=theta/np.pi,
                         physical=physical.tolist(), sequential=sequential.tolist(),
                         wrong_coarse=wrong.tolist(), transfer=transfer.tolist(),
                         statevec_sequential_error=seq_error, dephasing_tv=dephasing_tv,
                         transfer_tv=transfer_tv, transfer_state_l2=transfer_state_error,
                         periodic_diag_error=periodic_error,
                         periodic_diag_dephased_error=dephased_error,
                         physical_info=physical_info, transfer_info=transfer_info))
    exp.fail_check("C1", max_dephasing_tv > 1e-5,
                   f"max intermediate coarse-dephasing TV={max_dephasing_tv:.6g}")
    exp.fail_check("C2", max_transfer_tv > 1e-5 or max_transfer_state > 1e-5,
                   f"naive transfer max TV={max_transfer_tv:.6g}, state L2={max_transfer_state:.6g}")
    n13_setup = validate_n13_work_gates(exp)
    n13_physical, n13_state, n13_info = n13_physical_output(np.pi / 4)
    n13_grouped, n13_full = n13_grouped_reference(np.pi / 4)
    n13_full_error = float(np.max(np.abs(n13_physical - n13_full)))
    n13_grouped_error = float(np.max(np.abs(n13_grouped - n13_full)))
    exp.check("P5", n13_full_error < 3e-9 and n13_grouped_error < 3e-9
              and n13_info["norm_error"] < 2e-12,
              f"N=13 theta=pi/4: physical/full-r={n13_full_error:.2e}, "
              f"grouped-b6/full-r={n13_grouped_error:.2e}, norm={n13_info['norm_error']:.2e}")
    rows.append(dict(n13_theta_over_pi=.25, n13_physical=n13_physical.tolist(),
                     n13_full_r=n13_full.tolist(), n13_grouped_b6=n13_grouped.tolist(),
                     n13_physical_full_r_error=n13_full_error,
                     n13_grouped_full_r_error=n13_grouped_error,
                     n13_info=n13_info, n13_setup=n13_setup))
    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(N=N, base=BASE, orbit=ORBIT.tolist(), period=PERIOD,
                             block_size=BLOCK, width=WIDTH, insertion=2,
                             background="W01@s1,W12@s3,W01@s4",
                             angles=[float(x) for x in ANGLES],
                             derived_z=[1, -1, -1, -1, 1, 1],
                             dense_budget_bytes=MAX_DENSE_BYTES,
                             preflight=preflight,
                             extension="N=13,a=2,r=12,b3 physical W gates, grouped b6 reference",
                             references="existing Circuit/statevec and sequential_path",
                             no_amplitude_cutoff=True, elapsed_seconds=elapsed,
                             numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("physical_phase_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
