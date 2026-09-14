"""One localized symmetry-breaking kick amid periodic orbit mixers.

TODO 21's smallest physical test.  The known N=7, a=3 orbit is
[1,3,2,6,4,5], with b=3 and t=5.  Fixed repeated background blocks are
W01 after s=1, W12 after s=3, and W01 after s=4, all Rx(pi/2).  A single
localized Rx(theta) on orbit positions (0,1), i.e. physical labels (1,3) in
only the first block, is inserted after s=2.  This kick does not commute with
U^3 and breaks the shared coarse-sector promise.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The original ascending-order Circuit/statevec output agrees at every
      kick angle with the existing full-r sequential_path contraction using
      the same supplied branch matrices.
  P2  At theta=0 the periodic coarse-sector sampler agrees with the physical
      and full-r references; a repeated kick in every block restores that
      agreement at nonzero theta.
  P3  The intentionally wrong initial coarse dephasing
      (|0><0|+|3><3|)/2, propagated by the SAME full-r branches, differs from
      the physical output for at least one nonzero angle.
  P4  The localized four-term physical kick is the intended 8x8 unitary on
      every work input and leaves invalid labels 0 and 7 unchanged.

  C1  Wrong initial coarse dephasing must fail at one tested nonzero angle.
  C2  Treating the localized kick as a repeated block must fail at one
      tested nonzero angle (the restored repeated-kick case is separately
      checked as a positive control).

This is an indexed known-order validation, not a generic propagator, locality
claim, hardness result, or sampler for the broken-symmetry family.  Dense
physical states and full-r sequential factors are capped before allocation;
no amplitude cutoff is used.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_symmetry_kick
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
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import sequential_path
from modexp import ModExp
import statevec


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, WIDTH, PERIOD, BLOCK = 7, 3, 5, 6, 3
ORBIT = np.array([1, 3, 2, 6, 4, 5], dtype=np.int64)
BACKGROUND_THETA = np.pi / 2
KICK_ANGLES = np.pi * np.array([0., 1/8, 1/4, 1/2, 1.])


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def guard_width(width: int) -> None:
    if not isinstance(width, (int, np.integer)) or not 4 <= width <= 7:
        raise ValueError("fixed background requires 4<=width<=7 in this experiment")


def report_path(prefix: str = "symmetry_kick") -> Path:
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


# Small reusable fixtures for the main/bounds agent.  They intentionally live
# here rather than changing a production helper or another experiment file.
def orbit_fixture() -> np.ndarray:
    return ORBIT.copy()


def block_rotation(kind: str, theta: float) -> np.ndarray:
    if kind not in ("W01", "W12"):
        raise ValueError("unknown block kind")
    target, other = ((0, 1) if kind == "W01" else (1, 2))
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    result = np.eye(BLOCK, dtype=complex)
    result[target, target] = result[other, other] = c
    result[target, other] = result[other, target] = -1j * s
    return result


def localized_kick_block(theta: float) -> np.ndarray:
    """The 3x3 block acting on p=0,1,2 in the first orbit block."""
    return block_rotation("W01", theta)


def repeated_orbit_matrix(block: np.ndarray) -> np.ndarray:
    if np.asarray(block).shape != (BLOCK, BLOCK):
        raise ValueError("block must be 3x3")
    return np.kron(np.eye(PERIOD // BLOCK, dtype=complex), np.asarray(block, dtype=complex))


def localized_orbit_matrix(theta: float) -> np.ndarray:
    result = np.eye(PERIOD, dtype=complex)
    result[:BLOCK, :BLOCK] = localized_kick_block(theta)
    return result


def background_orbit_blocks() -> dict[int, np.ndarray]:
    return {1: repeated_orbit_matrix(block_rotation("W01", BACKGROUND_THETA)),
            3: repeated_orbit_matrix(block_rotation("W12", BACKGROUND_THETA)),
            4: repeated_orbit_matrix(block_rotation("W01", BACKGROUND_THETA))}


def background_sector_blocks() -> dict[int, np.ndarray]:
    return {1: block_rotation("W01", BACKGROUND_THETA),
            3: block_rotation("W12", BACKGROUND_THETA),
            4: block_rotation("W01", BACKGROUND_THETA)}


def background_work_gate(kind: str, theta: float, work_qubits: list[int]) -> Circuit:
    if kind == "W01":
        target = work_qubits[1]
        spectators = (1 << work_qubits[0]) | (1 << work_qubits[2])
    elif kind == "W12":
        target = work_qubits[0]
        spectators = (1 << work_qubits[1]) | (1 << work_qubits[2])
    else:
        raise ValueError("unknown background block")
    qc = Circuit(max(work_qubits) + 1)
    qc.rot((1 << target, 0), theta / 2)
    qc.rot((1 << target, spectators), -theta / 2)
    return qc


def localized_kick_gate(theta: float, work_qubits: list[int]) -> Circuit:
    """Physical X1(I-Z0)(I+Z2)/4 kick, active only on labels 1 and 3."""
    x1 = 1 << work_qubits[1]
    z0, z2 = 1 << work_qubits[0], 1 << work_qubits[2]
    qc = Circuit(max(work_qubits) + 1)
    for zmask, sign in ((0, +1), (z2, +1), (z0, -1), (z0 | z2, -1)):
        qc.rot((x1, zmask), sign * theta / 4)
    return qc


def _append_work_gate(qc: Circuit, kind: str, theta: float, work_qubits: list[int],
                      *, localized: bool = False) -> None:
    gate = (localized_kick_gate(theta, work_qubits) if localized
            else background_work_gate(kind, theta, work_qubits))
    qc.extend(gate)


def physical_circuit(theta: float, *, repeated_kick: bool = False,
                     width: int = WIDTH) -> tuple[Circuit, ModExp]:
    guard_width(width)
    me = ModExp(N, BASE, n_exp=width)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == 1:
            _append_work_gate(qc, "W01", BACKGROUND_THETA, me.x)
        elif i + 1 == 2:
            _append_work_gate(qc, "W01", theta, me.x, localized=not repeated_kick)
        elif i + 1 == 3:
            _append_work_gate(qc, "W12", BACKGROUND_THETA, me.x)
        elif i + 1 == 4:
            _append_work_gate(qc, "W01", BACKGROUND_THETA, me.x)
    qc.qft(me.exp, inverse=True)
    return qc, me


def physical_output(theta: float, *, repeated_kick: bool = False,
                    width: int = WIDTH) -> tuple[np.ndarray, dict]:
    qc, me = physical_circuit(theta, repeated_kick=repeated_kick, width=width)
    payload = (1 << me.n_qubits) * np.dtype(np.complex128).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"physical state allocation {payload} bytes exceeds 16 MiB")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    psi = statevec.run(qc, psi)
    probabilities = np.sum(np.abs(psi.reshape(1 << width, -1)) ** 2, axis=1)
    return probabilities, dict(qubits=me.n_qubits, vector_payload_bytes=psi.nbytes,
                               norm_error=abs(float(np.vdot(psi, psi).real) - 1),
                               gates=len(qc.gates), rotations=qc.n_nonclifford())


def full_orbit_branches(theta: float, *, repeated_kick: bool = False,
                        width: int = WIDTH) -> list[tuple[np.ndarray, np.ndarray]]:
    """Existing sequential_path branch matrices in original ascending order."""
    guard_width(width)
    blocks = background_orbit_blocks()
    blocks[2] = (repeated_orbit_matrix(localized_kick_block(theta))
                 if repeated_kick else localized_orbit_matrix(theta))
    identity = np.eye(PERIOD, dtype=complex)
    pairs = []
    for i in range(width):
        w = blocks.get(i + 1, identity)
        shift = np.zeros((PERIOD, PERIOD), dtype=complex)
        for j in range(PERIOD):
            shift[(j + (1 << i)) % PERIOD, j] = 1.
        pairs.append((w, w @ shift))
    return pairs


def orbit_basis(index: int) -> np.ndarray:
    if not 0 <= index < PERIOD:
        raise ValueError("orbit basis index outside supplied orbit")
    result = np.zeros(PERIOD, dtype=complex)
    result[index] = 1.
    return result


def sequential_marginal(theta: float, *, initial_label: int = 0,
                        repeated_kick: bool = False, width: int = WIDTH) -> np.ndarray:
    guard_width(width)
    guard((1 << width,), np.float64, "full output enumeration")
    pairs = full_orbit_branches(theta, repeated_kick=repeated_kick, width=width)
    initial = orbit_basis(initial_label)
    result = np.zeros(1 << width, dtype=float)
    for y in range(1 << width):
        result[y] = sequential_path(pairs, initial, output=y)["conditional_path_probability"]
    return result


def wrong_coarse_dephased(theta: float, *, repeated_kick: bool = False,
                          width: int = WIDTH) -> np.ndarray:
    # Same supplied full-r branches; only the initial density is replaced by
    # (|0><0|+|3><3|)/2.  Each pure component is contracted independently.
    return (sequential_marginal(theta, initial_label=0, repeated_kick=repeated_kick,
                                width=width)
            + sequential_marginal(theta, initial_label=3, repeated_kick=repeated_kick,
                                  width=width)) / 2


def periodic_marginal(theta: float, *, repeated_kick: bool = False,
                      width: int = WIDTH) -> np.ndarray:
    guard_width(width)
    guard((1 << width,), np.float64, "periodic output enumeration")
    blocks = background_sector_blocks()
    kick = (localized_kick_block(theta) if repeated_kick
            else None)
    if kick is not None:
        blocks[2] = kick
    elif 2 in blocks:
        del blocks[2]
    sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, width, blocks)
    result = np.zeros(1 << width, dtype=float)
    for alpha in range(sampler.sectors):
        for y in range(1 << width):
            result[y] += sampler.forced_joint(alpha, y)["joint_latent_output_probability"]
    return result


def validate_kick_matrix(exp: Experiment) -> dict:
    local = localized_kick_gate(BACKGROUND_THETA, [0, 1, 2])
    guard((8, 8), np.complex128, "localized kick work matrix")
    actual = local.to_unitary()
    expected = np.eye(8, dtype=complex)
    c, s = np.cos(BACKGROUND_THETA / 2), np.sin(BACKGROUND_THETA / 2)
    expected[1, 1] = expected[3, 3] = c
    expected[1, 3] = expected[3, 1] = -1j * s
    unitary_error = float(np.max(np.abs(actual.conj().T @ actual - np.eye(8))))
    matrix_error = float(np.max(np.abs(actual - expected)))
    statevec_error, invalid_error = 0., 0.
    for label in range(8):
        basis = statevec.basis(3, label)
        evolved = statevec.run(local, basis)
        statevec_error = max(statevec_error,
                             float(np.max(np.abs(evolved - actual @ basis))))
        if label in (0, 7):
            invalid_error = max(invalid_error, float(np.max(np.abs(evolved - basis))))
    exp.check("P4", unitary_error < 2e-14 and matrix_error < 2e-14
              and statevec_error < 2e-14 and invalid_error < 2e-14,
              f"kick: unitary={unitary_error:.2e}, matrix={matrix_error:.2e}, "
              f"statevec={statevec_error:.2e}, invalid={invalid_error:.2e}")
    return dict(unitary_error=unitary_error, matrix_error=matrix_error,
                statevec_error=statevec_error, invalid_label_error=invalid_error,
                actual=complex_payload(actual), expected=complex_payload(expected))


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(np.asarray(a) - np.asarray(b))) / 2)


def main() -> None:
    exp = Experiment("symmetry_kick", doc=__doc__)
    exp.predict("P1", "physical statevec agrees with full-r sequential_path for every angle")
    exp.predict("P2", "zero kick and repeated kick recover valid periodic references")
    exp.predict("P3", "wrong initial coarse dephasing loses observable interference somewhere")
    exp.predict("P4", "localized kick unitary matches all work inputs, including invalid labels")
    exp.must_fail("C1", "wrong initial coarse dephasing")
    exp.must_fail("C2", "treating localized kick as repeated in every block")
    start = time.perf_counter()
    rows = []
    setup = validate_kick_matrix(exp)
    max_wrong_tv = 0.
    max_repeated_mismatch = 0.
    for theta in KICK_ANGLES:
        theta = float(theta)
        physical, physical_info = physical_output(theta)
        sequential = sequential_marginal(theta)
        seq_error = float(np.max(np.abs(physical - sequential)))
        wrong = wrong_coarse_dephased(theta)
        wrong_tv = tv(wrong, physical)
        max_wrong_tv = max(max_wrong_tv, wrong_tv)
        exp.check("P1", seq_error < 3e-10 and physical_info["norm_error"] < 2e-12,
                  f"theta/pi={theta/np.pi:g}: statevec/sequential={seq_error:.2e}, "
                  f"norm={physical_info['norm_error']:.2e}")
        if theta == 0:
            periodic = periodic_marginal(theta)
            periodic_error = float(np.max(np.abs(periodic - physical)))
            exp.check("P2", periodic_error < 3e-10,
                      f"theta=0 periodic/physical={periodic_error:.2e}")
            exp.check("P2", wrong_tv < 3e-10,
                      f"theta=0 initial coarse dephasing is valid, TV={wrong_tv:.2e}")
        else:
            periodic_error = None
        repeated_physical, repeated_info = physical_output(theta, repeated_kick=True)
        repeated_seq = sequential_marginal(theta, repeated_kick=True)
        repeated_periodic = periodic_marginal(theta, repeated_kick=True)
        repeated_wrong = wrong_coarse_dephased(theta, repeated_kick=True)
        repeated_error = float(np.max(np.abs(repeated_physical - repeated_periodic)))
        repeated_seq_error = float(np.max(np.abs(repeated_physical - repeated_seq)))
        repeated_dephasing_error = float(np.max(np.abs(repeated_physical - repeated_wrong)))
        max_repeated_mismatch = max(max_repeated_mismatch, tv(repeated_physical, physical))
        exp.check("P2", repeated_error < 3e-10 and repeated_seq_error < 3e-10
                  and repeated_dephasing_error < 3e-10,
                  f"theta/pi={theta/np.pi:g}: repeated periodic={repeated_error:.2e}, "
                  f"sequential={repeated_seq_error:.2e}, "
                  f"coarse-dephasing={repeated_dephasing_error:.2e}")
        rows.append(dict(theta=theta, theta_over_pi=theta/np.pi,
                         physical=physical.tolist(), sequential=sequential.tolist(),
                         wrong_coarse=wrong.tolist(), repeated_physical=repeated_physical.tolist(),
                         repeated_periodic=repeated_periodic.tolist(),
                         repeated_wrong_coarse=repeated_wrong.tolist(),
                         statevec_sequential_error=seq_error, wrong_coarse_tv=wrong_tv,
                         periodic_error=periodic_error, repeated_error=repeated_error,
                         repeated_sequential_error=repeated_seq_error,
                         repeated_dephasing_error=repeated_dephasing_error,
                         physical_info=physical_info, repeated_info=repeated_info))

    exp.check("P3", max_wrong_tv > 1e-5,
              f"max wrong-coarse TV from physical={max_wrong_tv:.6g}")
    exp.fail_check("C1", max_wrong_tv > 1e-5,
                   f"max wrong-coarse TV from physical={max_wrong_tv:.6g}")
    exp.fail_check("C2", max_repeated_mismatch > 1e-5,
                   f"max repeated-kick vs localized TV={max_repeated_mismatch:.6g}")
    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(N=N, base=BASE, orbit=ORBIT.tolist(), period=PERIOD,
                             block_size=BLOCK, width=WIDTH,
                             background_theta=float(BACKGROUND_THETA),
                             kick_angles=[float(x) for x in KICK_ANGLES],
                             insertion_schedule="W01@s1, localized kick@s2, W12@s3, W01@s4",
                             dense_budget_bytes=MAX_DENSE_BYTES, physical_setup=setup,
                             full_r_branches="existing sequential_path; no generic propagator",
                             no_amplitude_cutoff=True, elapsed_seconds=elapsed,
                             numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("symmetry_kick_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
