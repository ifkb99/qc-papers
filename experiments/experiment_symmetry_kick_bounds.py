"""Actual-state distance bounds for one localized kick amid periodic mixers.

DERIVED BEFORE MEASUREMENT.  Use the fixed physical fixture N=7, a=3,
orbit [1,3,2,6,4,5], b=3, t=5, with repeated W01 after s=1, W12 after
s=3, and W01 after s=4.  Insert one local K=Rx(theta) mixing orbit indices
0,1 (physical labels 1,3) after s=2.  The unperturbed periodic sampler is
the cheap baseline.

Let rho_pre be the ACTUAL work state just before K, including the earlier
W01 mixer, and P_S the support projector onto physical labels 1,3.  The
candidate distance audit is

    F = Tr(P_S rho_pre)
    delta^2 = Tr((K-I)^dagger(K-I) rho_pre)
            = 4 sin(theta/4)^2 F
    g = Tr(K rho_pre)
      = 1 + (cos(theta/2)-1)F - i sin(theta/2)Tr(X_S rho_pre).

The output TV distance is bounded by both sqrt(delta^2) and the pure-state
trace distance sqrt(1-|g|^2), capped at one.  F is measured from the actual
pre-kick Circuit/statevec state; it is not replaced by bare residue counts.
The same existing Circuit/statevec construction computes the kicked output.

P1: theta=0 reproduces the unperturbed physical output and the cheap
    PeriodicOrbitCircuit sampler; all kicked rows remain normalized.
P2: actual F, overlap g and delta^2 identities agree with direct statevector
    measurements, and the kicked output TV obeys both distance bounds.
P3: the angle sweep theta=0,pi/8,pi/4,pi/2,pi records whether the cheap
    unperturbed sampler is an informative approximation.

C1: omitting the angle factor and asserting delta^2=4F must fail on a
    nonzero-angle row.  This is an identifiable false distance identity,
    not a claim that the correct bound is always tight.

This is a tiny physical validation, not a generic propagator or an
asymptotic claim.  Dense state allocations are guarded at 16 MiB.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_symmetry_kick_bounds
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit
from modexp import ModExp
import statevec

from experiments.experiment_periodic_sampler import (
    BLOCK, N, BASE, WIDTH, block_rotation, physical_output,
    sampler_marginal, work_mixer,
)


MAX_DENSE_BYTES = 16 * 1024 * 1024
PERIOD = 6
THETAS = (0.0, np.pi / 8, np.pi / 4, np.pi / 2, np.pi)
BACKGROUND = {1: "W01", 3: "W12", 4: "W01"}
SUPPORT_LABELS = (1, 3)


def guard(shape, dtype, label: str) -> int:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")
    return payload


def report_path(prefix: str = "symmetry_kick_bounds") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def local_kick(work: list[int], theta: float) -> Circuit:
    """K=Rx(theta) on orbit labels 0,1, i.e. physical work labels 1,3."""
    if len(work) != 3:
        raise ValueError("the physical fixture has three work qubits")
    target, require_one, require_zero = work[1], work[0], work[2]
    # X_target*(I-Z_require_one)*(I+Z_require_zero)/4.  Circuit.rot(theta)
    # implements exp(-i theta P/2), so each Pauli term gets angle theta/4.
    qc = Circuit(max(work) + 1)
    qc.rot((1 << target, 0), theta / 4)
    qc.rot((1 << target, 1 << require_zero), theta / 4)
    qc.rot((1 << target, 1 << require_one), -theta / 4)
    qc.rot((1 << target, (1 << require_one) | (1 << require_zero)), -theta / 4)
    return qc


def build_prefix_state() -> tuple[np.ndarray, ModExp]:
    """State after U^2 and W01@s1, immediately before the local kick."""
    me = ModExp(N, BASE, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp[:2]):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == 1:
            qc.extend(work_mixer("W01", np.pi / 2, me.x))
    guard((1 << me.n_qubits,), np.complex128, "pre-kick state")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    return statevec.run(qc, psi), me


def build_full(theta: float) -> tuple[np.ndarray, ModExp]:
    """Full original-time-order circuit with the local kick after s=2."""
    me = ModExp(N, BASE, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == 1:
            qc.extend(work_mixer("W01", np.pi / 2, me.x))
        elif i + 1 == 2:
            qc.extend(local_kick(me.x, theta))
        elif i + 1 == 3:
            qc.extend(work_mixer("W12", np.pi / 2, me.x))
        elif i + 1 == 4:
            qc.extend(work_mixer("W01", np.pi / 2, me.x))
    qc.qft(me.exp, inverse=True)
    guard((1 << me.n_qubits,), np.complex128, "full state")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    return statevec.run(qc, psi), me


def exponent_marginal(psi: np.ndarray, me: ModExp) -> np.ndarray:
    guard((1 << WIDTH,), np.float64, "exponent marginal")
    out = np.zeros(1 << WIDTH, dtype=float)
    for index, amplitude in enumerate(psi):
        out[statevec.read_register(index, me.exp)] += abs(amplitude) ** 2
    return out


def pre_support_observables(psi: np.ndarray, me: ModExp) -> tuple[float, complex, float]:
    """F, Tr(X_S rho_pre), and pre-state norm from the actual statevector."""
    f = 0.0
    x_expectation = 0j
    flip = 1 << me.x[1]
    for index, amplitude in enumerate(psi):
        label = statevec.read_register(index, me.x)
        if label in SUPPORT_LABELS:
            f += float(abs(amplitude) ** 2)
            x_expectation += np.conj(psi[index ^ flip]) * amplitude
    return f, complex(x_expectation), float(np.vdot(psi, psi).real)


def kick_direct_observables(psi: np.ndarray, me: ModExp, theta: float):
    kicked = statevec.run(local_kick(me.x, theta), psi.copy())
    overlap = complex(np.vdot(psi, kicked))
    delta2 = float(np.vdot(kicked - psi, kicked - psi).real)
    return kicked, overlap, delta2, float(np.vdot(kicked, kicked).real)


def tv(left: np.ndarray, right: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(left - right)))


def main() -> None:
    exp = Experiment("symmetry_kick_bounds", doc=__doc__)
    exp.predict("P1", "zero kick matches the periodic sampler baseline and all rows normalize")
    exp.predict("P2", "actual F/g/delta identities hold and output TV obeys both bounds")
    exp.predict("P3", "angle sweep records the informativeness of the unperturbed approximation")
    exp.must_fail("C1", "the angle-free delta²=4F identity is false")

    # Existing cheap periodic sampler baseline, independently checked against
    # its existing physical Circuit/statevec reference.
    blocks = {1: block_rotation("W01", np.pi / 2),
              3: block_rotation("W12", np.pi / 2),
              4: block_rotation("W01", np.pi / 2)}
    sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, WIDTH, blocks)
    cheap, sampler_mass_error = sampler_marginal(sampler)
    unperturbed, unperturbed_info = physical_output(BACKGROUND, WIDTH)
    baseline_error = float(np.max(np.abs(cheap - unperturbed)))
    exp.check("P1", baseline_error < 3e-10 and sampler_mass_error < 3e-11
              and unperturbed_info["norm_error"] < 3e-12,
              f"cheap baseline vs Circuit={baseline_error:.2e}, "
              f"mass={sampler_mass_error:.2e}, norm={unperturbed_info['norm_error']:.2e}")

    pre, me = build_prefix_state()
    F, x_expectation, pre_norm = pre_support_observables(pre, me)
    rows = []
    max_identity_error = max_overlap_error = max_delta_error = 0.0
    max_bound_violation = 0.0
    for theta in THETAS:
        theta = float(theta)
        kicked_prefix, overlap_direct, delta_direct, kick_norm = kick_direct_observables(pre, me, theta)
        overlap_formula = (1. + (np.cos(theta / 2.) - 1.) * F
                           - 1j * np.sin(theta / 2.) * x_expectation)
        delta_formula = 4. * np.sin(theta / 4.) ** 2 * F
        vector_bound = min(1., math.sqrt(max(0., delta_formula)))
        overlap_bound = min(1., math.sqrt(max(0., 1. - abs(overlap_formula) ** 2)))
        output_state, out_me = build_full(theta)
        output = exponent_marginal(output_state, out_me)
        output_tv = tv(output, cheap)
        bound = min(vector_bound, overlap_bound)
        violation = output_tv - bound
        max_bound_violation = max(max_bound_violation, violation)
        identity_error = abs(delta_direct - delta_formula)
        overlap_error = abs(overlap_direct - overlap_formula)
        delta_error = abs(float(np.vdot(kicked_prefix - pre, kicked_prefix - pre).real)
                          - delta_formula)
        max_identity_error = max(max_identity_error, identity_error)
        max_overlap_error = max(max_overlap_error, overlap_error)
        max_delta_error = max(max_delta_error, delta_error)
        exp.check("P2", identity_error < 3e-11 and overlap_error < 3e-11
                  and delta_error < 3e-11 and abs(pre_norm - 1.) < 3e-12
                  and abs(kick_norm - 1.) < 3e-12
                  and abs(output_state.dot(output_state.conj()).real - 1.) < 3e-12
                  and output_tv <= bound + 3e-10,
                  f"theta/pi={theta / np.pi:.6g}: F={F:.9g}, TV={output_tv:.9g}, "
                  f"vector={vector_bound:.9g}, overlap={overlap_bound:.9g}, "
                  f"gerr={overlap_error:.2e}, derr={delta_error:.2e}")
        rows.append(dict(theta=theta, theta_over_pi=theta / np.pi, F=F,
                         x_expectation_real=x_expectation.real,
                         x_expectation_imag=x_expectation.imag,
                         pre_norm=pre_norm, kick_norm=kick_norm,
                         overlap_direct_real=overlap_direct.real,
                         overlap_direct_imag=overlap_direct.imag,
                         overlap_formula_real=overlap_formula.real,
                         overlap_formula_imag=overlap_formula.imag,
                         delta2_direct=delta_direct, delta2_formula=delta_formula,
                         vector_bound=vector_bound, overlap_bound=overlap_bound,
                         output_tv_from_cheap=output_tv,
                         bound_violation=violation, output=output.tolist()))
    cover = sampler.localized_kick_bound(2, (0, 1))
    cover_fraction = cover["hit_probability_numerator"]/cover["hit_probability_denominator"]
    exp.check("P2", F <= cover_fraction + 3e-12,
              f"actual F={F:.9g} <= co-moving support cover={cover_fraction:.9g}")
    exp.check("P3", rows[0]["output_tv_from_cheap"] < 3e-10
              and rows[-1]["output_tv_from_cheap"] > 1e-6,
              f"theta=0 TV={rows[0]['output_tv_from_cheap']:.3e}; "
              f"theta=pi TV={rows[-1]['output_tv_from_cheap']:.6g}; F={F:.6g}")

    # Deliberately omit the angle factor.  At theta=pi/4 this predicts a
    # wrong squared state distance 4F instead of 4 sin²(theta/4)F.
    wrong_angle_error = abs(rows[2]["delta2_direct"] - 4. * F)
    exp.fail_check("C1", wrong_angle_error > 1e-5,
                   f"theta=pi/4 omitted-angle delta² error={wrong_angle_error:.9g}")

    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(
        N=N, base=BASE, period=PERIOD, block_size=BLOCK, width=WIDTH,
        background=BACKGROUND, support_labels=SUPPORT_LABELS,
        kick="Rx(theta) on physical labels 1,3 after s=2",
        theta_values=[float(x) for x in THETAS], actual_pre_kick_F=F,
        pre_kick_X_expectation_real=x_expectation.real,
        pre_kick_X_expectation_imag=x_expectation.imag,
        baseline_cheap_vs_physical_error=baseline_error,
        sampler_mass_error=sampler_mass_error,
        max_overlap_identity_error=max_overlap_error,
        max_delta_identity_error=max_delta_error,
        max_output_bound_violation=max_bound_violation,
        omitted_angle_control_error=wrong_angle_error,
        production_lightcone_cover=cover,
        dense_budget_bytes=MAX_DENSE_BYTES,
        reference="existing Circuit/statevec plus PeriodicOrbitCircuit baseline",
        no_new_propagator=True))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("symmetry_kick_bounds_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
