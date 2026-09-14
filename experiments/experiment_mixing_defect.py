"""A genuine orbit-preserving work-label mixing defect through the inverse QFT.

DERIVED BEFORE MEASUREMENT.  Fix N=7, a=3, t=4, split=2 and retain the
ascending-power order U, U^2; V(theta); U^4, U^8; inverse QFT.  The physical
work operation is the two-level rotation mixing |1> and |5>, identity on the
other six three-bit work labels.  Its generator is

    X_work2 (I-Z_work0) (I+Z_work1) / 4,

expanded into four commuting Pauli rotations with signed angle theta/4.

P1: the four-term Circuit implementation has the stated full 8x8 work
    unitary, preserves the six-state clean modular-multiplication orbit, and
    the spectral effects agree with both the existing full state-vector
    circuit and an independent clean-basis Fourier calculation.
P2: the effects form a complete POVM; their off-diagonal response reconstructs
    the difference between the coherent clean input and initial spectral
    dephasing, with the response recorded rather than interpreted as memory.
P3: the zero-angle row is ideal, while nonzero mixing produces a measurable
    output change and is not reproduced by the old diagonal phase-offset
    scalar model.
P4: extended-precision spectral arithmetic agrees with the float64 result and
    with an extended-state gate-level replay.

C1: the old ideal latent-eigenphase scalar model must fail on the sweep.
C2: the initial-eigenphase-dephased scalar model must fail on the sweep.
C3: moving V to the end must fail (the end operation is traced out and returns
    the ideal output, unlike the original time ordering).
C4: the old diagonal work-Rz phase-offset scalar model must fail; this mixing
    operation is not silently a product exponent phase.

Only theta varies in the main sweep: -pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi.
This is a perturbed-circuit validation experiment, not ideal Shor and not a
speed benchmark.  Dense arrays are capped at tiny dimensions and all setup,
reference and leakage costs are recorded.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_mixing_defect
"""
from __future__ import annotations

import time

import numpy as np

from circuits import Circuit
from experiments.experiment_spectral_defect import (
    orbit_by_return, orbit_fourier, probability_rows, scalar_mixture,
)
from lab import Experiment
from lab.spectral import coherence_response, single_defect_effects
from modexp import ModExp
import statevec


N, BASE, WIDTH, SPLIT = 7, 3, 4, 2
THETAS = np.pi * np.array([-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0])


def mixing_terms(x):
    """Four commuting Pauli strings for X_2(I-Z_0)(I+Z_1)/4."""
    x0, x1, x2 = x
    xm = 1 << x2
    z0, z1 = 1 << x0, 1 << x1
    return ((xm, 0, +1), (xm, z1, +1), (xm, z0, -1),
            (xm, z0 | z1, -1))


def mixing_circuit(theta, x):
    qc = Circuit(max(x) + 1)
    for xmask, zmask, sign in mixing_terms(x):
        qc.rot((xmask, zmask), sign * theta / 4)
    return qc


def mixing_work_matrix(theta, dtype=np.float64):
    """Independent 8x8 computational-basis matrix of the stated defect."""
    real = np.dtype(dtype).type
    c, s = np.cos(real(theta) / 2), np.sin(real(theta) / 2)
    V = np.eye(8, dtype=np.result_type(dtype, 1j))
    # The active projector is work bit 0 = 1, bit 1 = 0, so X_2 swaps 1/5.
    V[1, 1] = V[5, 5] = c
    V[1, 5] = V[5, 1] = -1j * s
    return V


def mixing_orbit_matrix(orbit, theta, dtype=np.float64):
    lookup = {int(label): i for i, label in enumerate(orbit)}
    V8 = mixing_work_matrix(theta, dtype)
    return V8[np.ix_(orbit, orbit)], lookup


def coherent_reference(theta, *, split=SPLIT, dtype=np.complex128):
    """Existing gate-level Fourier arithmetic, in original time order."""
    me = ModExp(N, BASE, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == split:
            qc.extend(mixing_circuit(theta, me.x))
    qc.qft(me.exp, inverse=True)
    psi = np.zeros(1 << me.n_qubits, dtype=dtype)
    psi[0] = 1
    psi = statevec.run(qc, psi)
    probabilities = np.sum(np.abs(psi.reshape(1 << WIDTH, -1)) ** 2, axis=1)
    clean_mask = sum(1 << q for q in me.x + me.exp)
    ids = np.arange(psi.size)
    leakage = np.sum(np.abs(psi[(ids & ~clean_mask) != 0]) ** 2)
    invalid = ((ids >> me.x[0]) & ((1 << me.n) - 1)) >= N
    invalid_leakage = np.sum(np.abs(psi[invalid]) ** 2)
    return probabilities, dict(
        qubits=qc.n, rotations=len(qc.gates), vector_payload_bytes=psi.nbytes,
        scratch_leakage=float(leakage), invalid_work_leakage=float(invalid_leakage),
    )


def clean_fourier_reference(theta, *, split=SPLIT, initial_label=1,
                            dtype=np.complex128):
    """Independent computational-basis clean-orbit formula plus inverse QFT."""
    Q, L = 1 << WIDTH, 1 << split
    orbit = orbit_by_return(N, BASE)
    V, lookup = mixing_orbit_matrix(orbit, theta, np.float64)
    r = len(orbit)
    joint = np.zeros((Q, r), dtype=dtype)
    for e in range(Q):
        l, h = e % L, e // L
        # V acts after U^l, then U^(Lh) shifts the orbit label.
        source = lookup[initial_label]
        for j in range(r):
            final = (j + L * h) % r
            joint[e, final] += V[j, (source + l) % r] / np.sqrt(Q)
    return np.sum(np.abs(np.fft.fft(joint, axis=0) / np.sqrt(Q)) ** 2, axis=1), dict(
        clean_joint_payload_bytes=joint.nbytes,
    )


def mixing_eigenbasis(orbit, theta, dtype=np.float64):
    V, _ = mixing_orbit_matrix(orbit, theta, dtype)
    F = orbit_fourier(len(orbit), dtype)
    return F.conj().T @ V @ F


def expected_full_unitary(theta, n, x):
    """Full-register expected matrix, identity outside the 1/5 work pair."""
    U = np.eye(1 << n, dtype=complex)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    for spectator in range(1 << n):
        if ((spectator >> x[0]) & 7) != 0:
            continue
        # spectator has all work bits clear; enumerate the non-work part.
        base = spectator
        i1 = base | (1 << x[0])
        i5 = base | (1 << x[0]) | (1 << x[2])
        U[i1, i1] = U[i5, i5] = c
        U[i1, i5] = U[i5, i1] = -1j * s
    return U


def main():
    exp = Experiment("mixing_defect", doc=__doc__)
    exp.predict("P1", "four-term unitary, orbit preservation, spectral/statevec/clean-FFT agreement")
    exp.predict("P2", "complete effects and recorded off-diagonal response reconstruct dephasing")
    exp.predict("P3", "zero is ideal; nonzero mixing differs and diagonal offset is not exact")
    exp.predict("P4", "extended precision agrees with float64 and extended gate replay")
    exp.must_fail("C1", "old ideal latent-eigenphase scalar model")
    exp.must_fail("C2", "initial physical-work eigenphase dephasing")
    exp.must_fail("C3", "forbidden move of defect to the end")
    exp.must_fail("C4", "old diagonal work-Rz phase-offset scalar model")

    orbit = orbit_by_return(N, BASE)
    period = len(orbit)
    me = ModExp(N, BASE, n_exp=WIDTH)
    ideal = scalar_mixture(period)
    rows = []
    max_ideal_error = max_dephasing_error = 0.0
    max_diagonal_error = max_end_error = 0.0

    # Full-register decomposition/unitarity and clean-orbit checks precede the sweep.
    theta_check = np.pi / 2
    start = time.perf_counter()
    physical = mixing_circuit(theta_check, me.x)
    dense_budget_bytes = 16 * 1024 * 1024
    full_unitary_shape = (1 << physical.n, 1 << physical.n)
    planned_full_unitary_payload_bytes = (full_unitary_shape[0] * full_unitary_shape[1]
                                          * np.dtype(np.complex128).itemsize)
    if planned_full_unitary_payload_bytes > dense_budget_bytes:
        raise MemoryError("full defect unitary exceeds the 16 MiB experiment budget")
    actual_full = physical.to_unitary()
    full_unitary_payload_bytes = actual_full.nbytes
    expected_full = expected_full_unitary(theta_check, physical.n, me.x)
    full_unitary_error = float(np.max(np.abs(actual_full - expected_full)))
    work_actual = actual_full[np.ix_(
        [label << me.x[0] for label in range(8)],
        [label << me.x[0] for label in range(8)])]
    work_expected = mixing_work_matrix(theta_check)
    work_unitary_error = float(np.max(np.abs(work_actual - work_expected)))
    orbit_complement = [i for i in range(1 << me.n) if i not in set(orbit)]
    orbit_preservation_error = float(np.max(np.abs(
        mixing_work_matrix(theta_check)[np.ix_(orbit_complement, orbit)])))
    unitary_seconds = time.perf_counter() - start
    exp.check("P1", full_unitary_error < 1e-12 and work_unitary_error < 1e-12
              and orbit_preservation_error < 1e-12,
              f"full={full_unitary_error:.2e}, work={work_unitary_error:.2e}, "
              f"orbit leakage={orbit_preservation_error:.2e}")

    setup = dict(N=N, a=BASE, width=WIDTH, split=SPLIT, period=period,
                 orbit=orbit.tolist(), orbit_labels=[int(x) for x in orbit],
                 exponent_order="ascending powers 1,2,4,8",
                 insertion="after powers 1,2 and before powers 4,8",
                 physical_generator="X_work2*(I-Z_work0)*(I+Z_work1)/4",
                 pauli_terms=[dict(xmask=xm, zmask=zm, signed_angle_factor=sign)
                              for xm, zm, sign in mixing_terms(me.x)],
                 theta_values=[float(x) for x in THETAS],
                 unitary_check_theta=float(theta_check),
                 full_unitary_error=full_unitary_error,
                 work_unitary_error=work_unitary_error,
                 orbit_preservation_error=orbit_preservation_error,
                 unitary_seconds=unitary_seconds,
                 dense_array_budget_bytes=dense_budget_bytes,
                 max_dense_payload_bytes=full_unitary_payload_bytes,
                 max_dense_payload_shape=list(full_unitary_shape),
                 max_dense_payload_dtype="complex128",
                 task="perturbed clean-input exponent-output distribution",
                 cost="dense effects and full state vector are validation only; no speed claim")

    for theta in THETAS:
        theta = float(theta)
        start = time.perf_counter()
        eigen_defect = mixing_eigenbasis(orbit, theta)
        effects = single_defect_effects(period, WIDTH, SPLIT, eigen_defect)
        effect_seconds = time.perf_counter() - start
        exact, dephased, contributions = probability_rows(effects)
        start = time.perf_counter()
        gate, gate_info = coherent_reference(theta)
        gate_seconds = time.perf_counter() - start
        clean_fft, fft_info = clean_fourier_reference(theta)
        spectral_gate_error = float(np.max(np.abs(exact - gate)))
        spectral_fft_error = float(np.max(np.abs(exact - clean_fft)))
        exp.check("P1", spectral_gate_error < 1e-9 and spectral_fft_error < 1e-11
                  and max(gate_info["scratch_leakage"], gate_info["invalid_work_leakage"]) < 1e-18,
                  f"theta/pi={theta/np.pi:g}: gate={spectral_gate_error:.2e}, "
                  f"clean FFT={spectral_fft_error:.2e}, leakage="
                  f"{gate_info['scratch_leakage']:.2e}")
        completeness = float(np.max(np.abs(effects.sum(axis=0) - np.eye(period))))
        minimum_eigenvalue = float(np.linalg.eigvalsh(effects).min())
        reconstruction = float(np.max(np.abs(exact - dephased - contributions.sum(axis=1))))
        response = coherence_response(effects)
        exp.check("P2", completeness < 1e-10 and minimum_eigenvalue > -1e-12
                  and reconstruction < 1e-12,
                  f"theta/pi={theta/np.pi:g}: complete={completeness:.2e}, "
                  f"min eig={minimum_eigenvalue:.2e}, pairs={reconstruction:.2e}")
        diagonal_offset = scalar_mixture(period, theta)
        dephasing_error = float(np.max(np.abs(exact - dephased)))
        ideal_error = float(np.max(np.abs(exact - ideal)))
        diagonal_error = float(np.max(np.abs(exact - diagonal_offset)))
        max_dephasing_error = max(max_dephasing_error, dephasing_error)
        max_ideal_error = max(max_ideal_error, ideal_error)
        max_diagonal_error = max(max_diagonal_error, diagonal_error)
        if abs(theta) < 1e-15:
            exp.check("P3", ideal_error < 1e-11 and dephasing_error < 1e-11,
                      "zero mixing recovers ideal and dephased baselines")
        else:
            exp.check("P3", ideal_error > 1e-6 and diagonal_error > 1e-6,
                      f"theta/pi={theta/np.pi:g}: ideal error={ideal_error:.2e}, "
                      f"diagonal-offset error={diagonal_error:.2e}")
        rows.append(dict(
            series="strength", theta=theta, theta_over_pi=theta / np.pi,
            probabilities=exact.tolist(), dephased=dephased.tolist(), ideal=ideal.tolist(),
            diagonal_offset=diagonal_offset.tolist(), gate_reference=gate.tolist(),
            clean_fft_reference=clean_fft.tolist(), pair_contributions=contributions.tolist(),
            response=response, spectral_gate_error=spectral_gate_error,
            spectral_fft_error=spectral_fft_error, completeness_error=completeness,
            minimum_eigenvalue=minimum_eigenvalue, dephasing_error=dephasing_error,
            dephasing_tv=float(np.abs(exact - dephased).sum() / 2),
            ideal_error=ideal_error, diagonal_offset_error=diagonal_error,
            effect_payload_bytes=effects.nbytes, clean_fft_payload_bytes=fft_info["clean_joint_payload_bytes"],
            effect_seconds=effect_seconds, gate_seconds=gate_seconds,
            **gate_info))
        exp.log(f"theta/pi={theta/np.pi:g}: dephasing max={dephasing_error:.6g}, "
                f"ideal max={ideal_error:.6g}, diagonal max={diagonal_error:.6g}, "
                f"detectable pairs={response['detectable_pairs']}, rank={response['rank']}")

    # The invalid end-slid circuit is ideal because the final work register is traced.
    end, end_info = coherent_reference(np.pi / 2, split=WIDTH)
    middle = next(row for row in rows if abs(row["theta_over_pi"] - 0.5) < 1e-15)
    max_end_error = float(np.max(np.abs(middle["probabilities"] - end)))
    exp.fail_check("C1", max_ideal_error > 1e-3,
                   f"largest ideal scalar error={max_ideal_error:.6g}")
    exp.fail_check("C2", max_dephasing_error > 1e-3,
                   f"largest initial-dephasing error={max_dephasing_error:.6g}")
    exp.fail_check("C3", max_end_error > 1e-3 and np.max(np.abs(end - ideal)) < 1e-9,
                   f"middle-vs-end error={max_end_error:.6g}, end-vs-ideal="
                   f"{np.max(np.abs(end - ideal)):.2e}")
    exp.fail_check("C4", max_diagonal_error > 1e-3,
                   f"largest old diagonal-offset error={max_diagonal_error:.6g}")

    # Precision bug check at theta=pi/2: repeat the independent spectral route
    # in extended precision and replay the existing circuit with an extended state.
    theta_ld = np.arccos(np.longdouble(-1)) / 2
    ext_effects = single_defect_effects(
        period, WIDTH, SPLIT, mixing_eigenbasis(orbit, theta_ld, np.longdouble),
        dtype=np.longdouble)
    ext_probs, _, _ = probability_rows(ext_effects)
    ext_gate, _ = coherent_reference(theta_ld, dtype=np.clongdouble)
    precision_row = next(row for row in rows if abs(row["theta_over_pi"] - 0.5) < 1e-15)
    spectral_delta = float(np.max(np.abs(ext_probs - precision_row["probabilities"])))
    extended_gate_error = float(np.max(np.abs(ext_probs - ext_gate)))
    exp.check("P4", spectral_delta < 1e-12 and extended_gate_error < 1e-9,
              f"extended spectral delta={spectral_delta:.2e}, gate={extended_gate_error:.2e}")
    rows.append(dict(series="precision", float_mantissa=np.finfo(float).nmant,
                     extended_mantissa=np.finfo(np.longdouble).nmant,
                     spectral_delta=spectral_delta, extended_gate_error=extended_gate_error))

    exp.finish(report_path="out/mixing_defect.json", rows=rows,
               metadata={**setup, "numpy": np.__version__,
                         "end_control_info": end_info,
                         "precision": "float64/complex128 plus longdouble/clongdouble replay"})


if __name__ == "__main__":
    main()
