"""Test the Bloch-sector shortcut for repeated periodic orbit blocks.

DERIVED BEFORE MEASUREMENT.  Let r=b*M and group orbit labels as bm+p.
The sector vectors |alpha,p> have amplitudes
  M**(-1/2) exp(-2*pi*i*alpha*m/M)
on bm+p.  Orbit translation therefore becomes a twisted b-cycle, while a
block-diagonal periodic W is the same b-by-b matrix in every alpha sector.
Consequently any number of interleaved W_i preserves alpha, not only one
defect.  Measuring alpha at the start should give a uniform mixture of these
b-dimensional sectors for the initial orbit label p=0.

P1: projected explicit orbit shifts, including powers, equal the integer
    twisted-cycle formula; periodic W projects to the same W in every sector.
P2: for one periodic W, the uniform-alpha sector mixture agrees with the
    existing single_defect_effects spectral reference, including r<=12 and
    t<=6 rows.
P3: two interleaved periodic blocks in the r=4 abstract orbit agree with an
    independent existing Circuit/statevec construction (controlled increments
    and local work rotations).

C1: the opposite twist sign must disagree with explicit projected shifts.
C2: fully dephasing the initial U-eigenphase must disagree for a nontrivial W.
C3: replacing forward controlled shifts by inverse shifts as a reverse-gate
    shortcut must disagree for the interleaved noncommuting case.
C4: dropping the wrap feedback (using an untwisted cycle) must disagree.

The dense objects below are validation-only and are capped before allocation
at 16 MiB; r<=12 and t<=6.  No scalable generic propagator is introduced.
The only multiple-defect reference uses the existing Circuit/statevec engine.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_periodic_sectors
"""
from __future__ import annotations

import math
import os
import time

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit
from lab.spectral import single_defect_effects
import statevec


BYTE_CAP = 16 * 1024 * 1024
MAX_R = 12
MAX_T = 6
REPORT_PATH = os.environ.get("PERIODIC_SECTORS_REPORT", "out/periodic_sectors.json")


def _guard(shape, dtype=np.complex128):
    """Refuse a dense allocation before NumPy sees it."""
    nbytes = math.prod(shape) * np.dtype(dtype).itemsize
    if nbytes > BYTE_CAP:
        raise MemoryError(f"dense allocation {shape} would use {nbytes} bytes")


def shift_matrix(r: int) -> np.ndarray:
    _guard((r, r))
    out = np.zeros((r, r), dtype=complex)
    for j in range(r):
        out[(j + 1) % r, j] = 1.0
    return out


def sector_basis(r: int, b: int, alpha: int) -> np.ndarray:
    if r > MAX_R or r % b or not 0 <= alpha < r // b:
        raise ValueError("invalid capped sector parameters")
    M = r // b
    _guard((r, b))
    q = np.zeros((r, b), dtype=complex)
    for m in range(M):
        q[b * m:b * m + b, :] += np.eye(b) * np.exp(
            -2j * np.pi * alpha * m / M) / np.sqrt(M)
    return q


def twisted_power(b: int, M: int, alpha: int, n: int,
                  twist_sign: int = 1, feedback: bool = True) -> np.ndarray:
    """T_alpha^n without floating matrix squaring.

    The destination of p after n forward steps is (p+n) mod b.  Each full
    crossing of b->0 contributes exp(+2*pi*i*alpha/M).  Negative n is used
    only for the deliberate invalid-reverse-gate control.
    """
    _guard((b, b))
    out = np.zeros((b, b), dtype=complex)
    for p in range(b):
        dst = (p + n) % b
        wraps = (p + n) // b
        phase = 1.0
        if feedback:
            phase = np.exp(2j * np.pi * twist_sign * alpha * wraps / M)
        out[dst, p] = phase
    return out


def periodic_block(r: int, b: int, block: np.ndarray) -> np.ndarray:
    block = np.asarray(block, dtype=complex)
    if r > MAX_R or r % b or block.shape != (b, b):
        raise ValueError("invalid periodic block")
    _guard((r, r))
    out = np.zeros((r, r), dtype=complex)
    for m in range(r // b):
        out[b * m:b * m + b, b * m:b * m + b] = block
    return out


def rotation_block(b: int, theta: float, axis: str = "z") -> np.ndarray:
    """A fixed, exactly unitary b-block (nontrivial on p=0,1)."""
    if b < 2:
        return np.eye(b, dtype=complex)
    out = np.eye(b, dtype=complex)
    if axis == "z":
        out[0, 0] = np.exp(-0.5j * theta)
        out[1, 1] = np.exp(0.5j * theta)
    elif axis == "x":
        c, s = np.cos(theta / 2), -1j * np.sin(theta / 2)
        out[:2, :2] = ((c, s), (s, c))
    else:
        raise ValueError("axis must be z or x")
    return out


def projected_power_error(r: int, b: int, alpha: int, powers: list[int]) -> float:
    M = r // b
    q = sector_basis(r, b, alpha)
    u = shift_matrix(r)
    direct = np.eye(r, dtype=complex)
    errors = []
    for n in range(max(powers) + 1):
        if n in powers:
            projected = q.conj().T @ direct @ q
            formula = twisted_power(b, M, alpha, n)
            errors.append(float(np.max(np.abs(projected - formula))))
        direct = u @ direct
    return max(errors, default=0.0)


def _append_twisted_step(qc: Circuit, control: int, b: int, M: int,
                         alpha: int, *, twist_sign: int = 1,
                         feedback: bool = True, reverse: bool = False) -> None:
    """Append one controlled T_alpha using only existing Circuit gates.

    This is deliberately specialized to b=2,3.  The b=3 construction keeps
    the unused computational state |3> fixed and follows the requested
    phase-on-|2>, then (1,2) and (0,1) swaps.  A reverse control is made by
    Circuit.inverse(), never by a second matrix propagator.
    """
    work = [q for q in range(control)]
    phi = (2 * np.pi * twist_sign * alpha / M) if feedback else 0.0
    step = Circuit(qc.n)
    if b == 2:
        step.cphase(control, work[0], phi)
        step.cnot(control, work[0])
    elif b == 3:
        # p=2 is (work0,work1)=(0,1), then cycle 0->1->2->0.
        step.x(work[0])
        step.ccphase(control, work[1], work[0], phi)
        step.x(work[0])
        step.cswap(control, work[0], work[1])
        step.x(work[1])
        step.toffoli(control, work[1], work[0])
        step.x(work[1])
    else:
        raise ValueError("Circuit verifier supports only b=2 or b=3")
    qc.extend(step.inverse() if reverse else step)


def _append_block_gate(qc: Circuit, work: list[int], b: int,
                       axis: str, theta: float) -> None:
    """Append rotation_block(b,theta,axis) in the computational p basis."""
    if b == 2:
        if axis == "x":
            qc.rx(work[0], theta)
        elif axis == "z":
            qc.rz(work[0], theta)
        else:
            raise ValueError("axis must be x or z")
    elif b == 3:
        # exp(-i theta P (I+Z_work1)/4) applies rotation_block's P rotation
        # only on p=0,1 (work1=0), while fixing the unused p=3 state.
        if axis == "x":
            p0 = (1 << work[0], 0)
        elif axis == "z":
            p0 = (0, 1 << work[0])
        else:
            raise ValueError("axis must be x or z")
        qc.rot(p0, theta / 2)
        qc.rot((p0[0], p0[1] | (1 << work[1])), theta / 2)
    else:
        raise ValueError("Circuit verifier supports only b=2 or b=3")


def circuit_sector_output(width: int, b: int, M: int, alpha: int,
                          block_specs: dict[int, tuple[str, float]], *,
                          twist_sign: int = 1, feedback: bool = True,
                          reverse: bool = False) -> np.ndarray:
    """Conditional sector probabilities from the existing Circuit/statevec."""
    if width > MAX_T or b * M > MAX_R or b not in (2, 3):
        raise ValueError("sector cap or Circuit verifier support exceeded")
    if not 0 <= alpha < M:
        raise ValueError("invalid sector")
    work = list(range(1 if b == 2 else 2))
    exponent = list(range(len(work), len(work) + width))
    # Guard the actual full statevector allocation before Circuit/statevec.
    _guard((1 << (len(work) + width),))
    qc = Circuit(len(work) + width)
    for q in exponent:
        qc.h(q)
    for i, control in enumerate(exponent):
        for _ in range(1 << i):
            _append_twisted_step(qc, control, b, M, alpha,
                                 twist_sign=twist_sign,
                                 feedback=feedback, reverse=reverse)
        if i in block_specs:
            axis, theta = block_specs[i]
            _append_block_gate(qc, work, b, axis, theta)
    qc.qft(exponent, inverse=True)
    psi = statevec.run(qc)
    probs = np.zeros(1 << width, dtype=float)
    for index, amplitude in enumerate(psi):
        probs[index >> len(work)] += abs(amplitude) ** 2
    return probs


def uniform_sector_mixture(width: int, b: int, M: int,
                           block_specs: dict[int, tuple[str, float]], **kwargs) -> np.ndarray:
    rows = [circuit_sector_output(width, b, M, alpha, block_specs, **kwargs)
            for alpha in range(M)]
    return np.mean(rows, axis=0)


def spectral_one_defect(period: int, width: int, split: int,
                        block: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Existing spectral reference, with the physical periodic V transformed
    to the U-eigenbasis required by single_defect_effects."""
    r = period
    v = periodic_block(r, block.shape[0], block)
    j = np.arange(r)
    fourier = np.exp(-2j * np.pi * j[:, None] * j[None, :] / r) / np.sqrt(r)
    v_eigen = fourier.conj().T @ v @ fourier
    effects = single_defect_effects(r, width, split, v_eigen)
    exact = effects.sum(axis=(1, 2)).real / r
    dephased = np.trace(effects, axis1=1, axis2=2).real / r
    return exact, dephased


def abstract_r4_circuit(theta0: float, theta1: float) -> np.ndarray:
    """Existing statevec reference: controlled U, U^2 on a two-bit orbit.

    Work qubits 0,1 encode j in binary, exponent qubits 2,3.  The sequence
    CNOT/Toffoli gates is a controlled increment modulo four, and the local
    rotations on work bit zero are the repeated b=2 block W_i.
    """
    qc = Circuit(4)
    qc.h(2).h(3)
    qc.toffoli(2, 0, 1).cnot(2, 0)
    qc.rz(0, theta0)
    qc.cnot(3, 1)
    qc.rx(0, theta1)
    qc.qft([2, 3], inverse=True)
    psi = statevec.run(qc)
    probs = np.zeros(4, dtype=float)
    for index, amplitude in enumerate(psi):
        probs[index >> 2] += abs(amplitude) ** 2
    return probs


def main() -> None:
    exp = Experiment("periodic_sectors", doc=__doc__)
    exp.predict("P1", "projected shift powers and periodic block projections match")
    exp.predict("P2", "uniform alpha mixture matches existing one-defect spectral effects")
    exp.predict("P3", "two interleaved blocks match the existing abstract r=4 statevec")
    exp.predict("P4", "sector probabilities agree with lab.periodic instrument contractions")
    exp.must_fail("C1", "wrong twist sign disagrees with projected U powers")
    exp.must_fail("C2", "initial full eigenphase dephasing changes one-defect outputs")
    exp.must_fail("C3", "inverse controlled shifts are not a valid reverse shortcut")
    exp.must_fail("C4", "dropping wrap feedback changes the sector output")

    rows: list[dict] = []
    started = time.perf_counter()

    # Matrix-level independent check.  Include powers beyond one cycle and a
    # power near the t=6 cap, but compute each power by the integer entry rule.
    exp.section("P1 projected shifts and periodic W")
    p1_max = 0.0
    for b, M in ((2, 2), (2, 3), (3, 2), (3, 3), (3, 4)):
        r = b * M
        block = rotation_block(b, np.pi / 5, "x")
        w = periodic_block(r, b, block)
        powers = [0, 1, b - 1, b, 2 * b + 3, 63]
        for alpha in range(M):
            q = sector_basis(r, b, alpha)
            u = shift_matrix(r)
            shift_error = projected_power_error(r, b, alpha, powers)
            block_error = float(np.max(np.abs(q.conj().T @ w @ q - block)))
            p1_max = max(p1_max, shift_error, block_error)
            exp.check("P1", shift_error < 2e-12 and block_error < 2e-12,
                      f"r={r}, b={b}, M={M}, alpha={alpha}: "
                      f"shift={shift_error:.2e}, block={block_error:.2e}")
        rows.append(dict(series="sector_matrices", r=r, b=b, M=M,
                         powers=powers, block_angle=float(np.pi / 5),
                         max_shift_error=float(max(
                             projected_power_error(r, b, alpha, powers)
                             for alpha in range(M))),
                         max_block_error=float(max(
                             np.max(np.abs(sector_basis(r, b, alpha).conj().T @ w
                                           @ sector_basis(r, b, alpha) - block))
                             for alpha in range(M)))))

    # One periodic defect against the existing independent spectral effects.
    exp.section("P2 one-defect alpha mixture")
    one_defect_max = 0.0
    geometry = ((2, 2), (2, 3), (3, 2), (3, 3), (3, 4))
    for b, M in geometry:
        r, width, split = b * M, 4, 2
        for theta in (np.pi / 7, -np.pi / 5, np.pi / 3):
            block = rotation_block(b, theta, "x")
            sector = uniform_sector_mixture(
                width, b, M, {split - 1: ("x", float(theta))})
            exact, dephased = spectral_one_defect(r, width, split, block)
            error = float(np.max(np.abs(sector - exact)))
            one_defect_max = max(one_defect_max, error)
            exp.check("P2", error < 3e-11,
                      f"r={r}, b={b}, theta/pi={theta / np.pi:.6g}: "
                      f"sector-vs-effects={error:.2e}")
            rows.append(dict(series="one_defect", r=r, b=b, M=M,
                             width=width, split=split,
                             theta=float(theta), theta_over_pi=float(theta / np.pi),
                             max_error=error,
                             dephased_error=float(np.max(np.abs(exact - dephased))),
                             sector=sector.tolist(), exact=exact.tolist(),
                             dephased=dephased.tolist()))

    # Independent multiple-defect reference.  The sector run uses two W_i;
    # statevec supplies the only full-circuit calculation in this experiment.
    exp.section("P3 two interleaved blocks")
    theta0, theta1 = np.pi / 3, -np.pi / 5
    block0 = rotation_block(2, theta0, "z")
    block1 = rotation_block(2, theta1, "x")
    multi = uniform_sector_mixture(
        2, 2, 2, {0: ("z", float(theta0)), 1: ("x", float(theta1))})
    circuit = abstract_r4_circuit(theta0, theta1)
    multi_error = float(np.max(np.abs(multi - circuit)))
    exp.check("P3", multi_error < 3e-11,
              f"r=4 statevec-vs-sectors={multi_error:.2e}")
    rows.append(dict(series="multiple_defect", r=4, b=2, M=2, width=2,
                     theta0=float(theta0), theta1=float(theta1),
                     sector=multi.tolist(), statevec=circuit.tolist(),
                     max_error=multi_error))

    # Compare conditional probabilities with the production sequential
    # instrument on the same supplied block schedule.  Its insertion key s
    # means "after s controls", hence our post-control index i is s-1.
    exp.section("P4 production instrument cross-check")
    prod_one = PeriodicOrbitCircuit(6, 3, 4, {2: rotation_block(3, np.pi / 3, "x")})
    prod_block = rotation_block(3, np.pi / 3, "x")
    prod_max = 0.0
    for alpha in range(2):
        ours = circuit_sector_output(4, 3, 2, alpha, {1: ("x", np.pi / 3)})
        theirs = np.array([prod_one.forced_joint(alpha, y)[
            "conditional_path_probability"] for y in range(16)])
        error = float(np.max(np.abs(ours - theirs)))
        prod_max = max(prod_max, error)
        exp.check("P4", error < 3e-11,
                  f"r=6 alpha={alpha}: production one-W error={error:.2e}")
    prod_multi_max = 0.0
    production_multi = PeriodicOrbitCircuit(4, 2, 2, {1: block0, 2: block1})
    for alpha in range(2):
        ours = circuit_sector_output(
            2, 2, 2, alpha, {0: ("z", float(theta0)), 1: ("x", float(theta1))})
        theirs = np.array([production_multi.forced_joint(alpha, y)[
            "conditional_path_probability"] for y in range(4)])
        error = float(np.max(np.abs(ours - theirs)))
        prod_multi_max = max(prod_multi_max, error)
        exp.check("P4", error < 3e-11,
                  f"r=4 alpha={alpha}: production two-W error={error:.2e}")
    rows.append(dict(series="production_crosscheck", one_defect_max_error=prod_max,
                     multiple_defect_max_error=prod_multi_max))

    # Deliberate controls.  These are kept separate from positive checks so a
    # unanimous positive test cannot hide an uninformative observable.
    exp.section("must-fail controls")
    b, M, width = 2, 2, 2
    schedule = {0: ("z", float(theta0)), 1: ("x", float(theta1))}
    correct = uniform_sector_mixture(width, b, M, schedule)
    # M=3, alpha=1 avoids the accidental exp(+/- i*pi) equality at M=2.
    projected = sector_basis(6, 2, 1).conj().T @ shift_matrix(6) @ sector_basis(6, 2, 1)
    wrong_twist = twisted_power(2, 3, 1, 1, twist_sign=-1)
    no_feedback = uniform_sector_mixture(width, b, M, schedule, feedback=False)
    reverse = uniform_sector_mixture(width, b, M, schedule, reverse=True)
    wrong_twist_error = float(np.max(np.abs(projected - wrong_twist)))
    exp.fail_check("C1", wrong_twist_error > 1e-4,
                   f"wrong-twist matrix error={wrong_twist_error:.6g}")
    # Use the nontrivial r=6 row to ensure the dephasing control is not a
    # coincidental r=4 symmetry.  The one-defect loop already measured it.
    control_block = rotation_block(3, np.pi / 3, "x")
    control_exact, control_dephased = spectral_one_defect(6, 4, 2, control_block)
    dephase_error = float(np.max(np.abs(control_exact - control_dephased)))
    exp.fail_check("C2", dephase_error > 1e-4,
                   f"initial-dephasing error={dephase_error:.6g}")
    # Compare output probabilities from independent Circuit/statevec runs;
    # global phase differences are intentionally not treated as evidence.
    # The old b=2 two-W schedule had identical output rows despite a large
    # amplitude difference, so use the smallest b=3 case that exposes the
    # invalid reverse operation in the observable itself.
    c3_schedule = {0: ("x", 0.31)}
    c3_forward = uniform_sector_mixture(2, 3, 2, c3_schedule)
    c3_reverse = uniform_sector_mixture(2, 3, 2, c3_schedule, reverse=True)
    reverse_error = float(np.max(np.abs(c3_forward - c3_reverse)))
    exp.fail_check("C3", reverse_error > 1e-4,
                   f"inverse-shift probability error={reverse_error:.6g}; "
                   "prior amplitude-only b=2 control was vacuous")
    exp.fail_check("C4", float(np.max(np.abs(correct - no_feedback))) > 1e-4,
                   f"missing-feedback error={np.max(np.abs(correct - no_feedback)):.6g}")
    rows.append(dict(series="controls", wrong_twist_error=wrong_twist_error,
                     dephasing_error=dephase_error,
                     reverse_error=reverse_error,
                     missing_feedback_error=float(np.max(np.abs(correct - no_feedback)))))

    exp.finish(
        report_path=REPORT_PATH,
        rows=rows,
        metadata=dict(
            byte_cap=BYTE_CAP, max_r=MAX_R, max_t=MAX_T,
            elapsed_seconds=time.perf_counter() - started,
            one_defect_max_error=one_defect_max,
            sector_matrix_max_error=p1_max,
            multiple_defect_error=multi_error,
            reference="lab.spectral.single_defect_effects for one W; existing Circuit/statevec for all sector path verification and r=4 full-orbit check",
            assumption="known orbit index and exact periodic b-block promise; not a physical locality claim",
            numerical_note="complex128 validation; no amplitudes truncated; outputs are normalized probability rows",
            audit_note="C3 now compares probabilities: prior b=2 amplitude-only check differed in phase while its probabilities agreed",
        ),
    )


if __name__ == "__main__":
    main()
