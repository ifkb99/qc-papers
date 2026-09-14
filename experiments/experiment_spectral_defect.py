"""Does one noncommuting work rotation actually defeat cheap output sampling?

DERIVED BEFORE MEASUREMENT. Fix N=7,a=3,t=4, split=2: ascending powers
1,2; V(theta); powers 4,8; inverse QFT. V=Rz(x0) preserves the clean orbit.
With Q=2^t,L=2^split, output Kraus matrices in the U eigenbasis factor as
K_y=diag(B_y) V diag(A_y), the two finite Fourier filters in lab.spectral.
Thus p(y)=sum(E_y)/r; initial spectral dephasing predicts trace(E_y)/r.

P1: spectral effects agree with the UNCHANGED existing statevec engine on
Fourier-compiled arithmetic, and a computational-basis FFT reference.
P2: effects are a complete POVM; paired off-diagonals reconstruct dephasing
error; a different clean input tests complex coordinates beyond the |1> input.
P3: theta=0 and the nontrivial commuting exp(-i theta X0X1X2/2) recover ideal
statistics (X0X1X2 = U^3 on this orbit). Precision sweep is a separate check.
P4: the surprising stronger baseline: low exponents l=0,1,2,3 reach 1,3,2,6,
whose low work bit is 1-l_1. Hence the work Rz is EXACTLY an exponent-bit-1
phase -theta on this input. A phase-adjusted scalar eigenphase sampler should
still succeed, even if the original physical-work dephasing model fails.
C1: discarding initial eigenphase coherences must fail somewhere at split=2.
C2: keeping the old ideal sampler unchanged must fail somewhere.
C3: moving the defect to the end (an invalid reordering) must fail somewhere.
C4: noncommutation alone must NOT imply changed outputs: the separate split=1
control applies only a global phase on reachable |1>,|3> before the defect.

Only theta varies in the main sweep. Controls are separate fixed-instance
comparisons. This is a PERTURBED circuit, not ideal Shor, and not a speedup
benchmark. Order/orbit/filter setup is explicit, dense, and paid for; no claim
that counting response rank bounds sampling memory. No new propagator.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_spectral_defect
"""
from __future__ import annotations
import time
import numpy as np
from circuits import Circuit
from lab import Experiment
from lab.reachable import compiled_pairs, action_stats
from lab.semiclassical import eigenphase_path
from lab.spectral import single_defect_effects, coherence_response
from modexp import ModExp
import statevec


N, BASE, WIDTH, SPLIT = 7, 3, 4, 2


def orbit_by_return(N, a):
    orbit, value = [1], a % N
    while value != 1:
        orbit.append(value)
        value = value*a % N
    return np.array(orbit)


def orbit_fourier(period, dtype=np.float64):
    real = np.dtype(dtype).type
    pi = np.arccos(real(-1))
    j = np.arange(period, dtype=dtype)
    # Column k has eigenvalue exp(+2*pi*i*k/r) for the forward orbit shift.
    return np.exp(-2j*pi*j[:, None]*j[None, :]/period) / np.sqrt(real(period))


def defect_matrix(orbit, theta, kind, dtype=np.float64):
    real = np.dtype(dtype).type
    theta = real(theta)
    r = len(orbit)
    if kind == "z":
        generator = np.diag(1 - 2*(orbit & 1)).astype(dtype)
    elif kind == "xxx":
        lookup = {int(x): j for j, x in enumerate(orbit)}
        generator = np.zeros((r, r), dtype=dtype)
        generator[[lookup[int(x) ^ 7] for x in orbit], np.arange(r)] = 1
    else:
        raise ValueError("unknown defect")
    F = orbit_fourier(r, dtype)
    V = np.cos(theta/2)*np.eye(r) - 1j*np.sin(theta/2)*generator
    return F.conj().T @ V @ F


def coherent_reference(theta, kind="z", split=SPLIT, dtype=np.complex128):
    """Existing gate-level Fourier arithmetic, in its original time ordering."""
    me = ModExp(N, BASE, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i+1 == split:
            if kind == "z":
                qc.rz(me.x[0], theta)
            else:
                qc.rot((sum(1 << qx for qx in me.x), 0), theta)
    qc.qft(me.exp, inverse=True)
    psi = np.zeros(1 << me.n_qubits, dtype=dtype)
    psi[0] = 1
    psi = statevec.run(qc, psi)
    probabilities = np.sum(np.abs(psi.reshape(1 << WIDTH, -1))**2, axis=1)
    clean_mask = sum(1 << q for q in me.x + me.exp)
    ids = np.arange(psi.size)
    leakage = np.sum(np.abs(psi[(ids & ~clean_mask) != 0])**2)
    invalid = ((ids >> me.x[0]) & ((1 << me.n)-1)) >= N
    invalid_leakage = np.sum(np.abs(psi[invalid])**2)
    return probabilities, dict(qubits=qc.n, rotations=len(qc.gates),
                               vector_payload_bytes=psi.nbytes,
                               scratch_leakage=float(leakage),
                               invalid_work_leakage=float(invalid_leakage))


def phase_fft_reference(theta, *, split=SPLIT, initial_label=1):
    """Independent clean basis calculation: Rz changes phases, not labels."""
    Q, L = 1 << WIDTH, 1 << split
    joint = np.zeros((Q, N), dtype=complex)
    for e in range(Q):
        early = initial_label*pow(BASE, e % L, N) % N
        final = initial_label*pow(BASE, e, N) % N
        joint[e, final] = np.exp(-.5j*theta*(1-2*(early & 1))) / np.sqrt(Q)
    return np.sum(np.abs(np.fft.fft(joint, axis=0)/np.sqrt(Q))**2, axis=1)


def scalar_mixture(period, theta=0.):
    phases = np.zeros(WIDTH)
    phases[1] = -theta
    return np.array([sum(eigenphase_path(period, WIDTH, k, output=y,
                                         input_phases=phases)
                         ["conditional_path_probability"] for k in range(period))/period
                     for y in range(1 << WIDTH)])


def probability_rows(effects):
    r = effects.shape[-1]
    exact = effects.sum(axis=(1, 2)).real / r
    dephased = np.trace(effects, axis1=1, axis2=2).real / r
    pairs = np.triu_indices(r, 1)
    contributions = 2*effects[:, pairs[0], pairs[1]].real / r
    return exact, dephased, contributions


def main():
    exp = Experiment("spectral_defect", doc=__doc__)
    exp.predict("P1", "spectral effects match coherent compiled and clean FFT references")
    exp.predict("P2", "complete positive effects and pair contributions explain dephasing error")
    exp.predict("P3", "zero/commuting controls recover ideal; precision agrees")
    exp.predict("P4", "proved exponent-phase replacement restores scalar sampling")
    exp.must_fail("C1", "initial physical-work eigenphase dephasing is not sufficient")
    exp.must_fail("C2", "unmodified ideal scalar sampler is not sufficient")
    exp.must_fail("C3", "sliding the noncommuting defect past later arithmetic is invalid")
    exp.must_fail("C4", "noncommutation alone does not force observable output change")
    rows = []
    start = time.perf_counter()
    orbit = orbit_by_return(N, BASE)
    orbit_seconds = time.perf_counter()-start
    period = len(orbit)
    start = time.perf_counter()
    me, pairs = compiled_pairs(N, BASE, WIDTH)
    compile_seconds = time.perf_counter()-start
    start = time.perf_counter()
    for i, pair in enumerate(pairs):
        ids = np.arange(N) << me.x[0]
        exp.check("P1", np.array_equal(pair[0](ids), ids) and
                  np.array_equal(pair[1](ids),
                                 (np.arange(N)*pow(BASE, 1 << i, N) % N) << me.x[0]),
                  f"existing Toffoli block {i}: all clean inputs, both branches")
    validation_seconds = time.perf_counter()-start
    # Check the input-specific exponent-phase replacement BEFORE using it.
    low_bits = np.array([pow(BASE, l, N) & 1 for l in range(1 << SPLIT)])
    exp.check("P4", np.array_equal(low_bits, 1-((np.arange(4) >> 1) & 1)),
              f"reached work LSB={low_bits.tolist()} equals 1 minus exponent bit 1")
    ideal = scalar_mixture(period)
    max_dephasing_error, max_ideal_error = 0., 0.
    for theta in np.pi*np.array([-1., -.5, -.25, -.125, 0., .125, .25, .5, 1.]):
        start = time.perf_counter()
        V = defect_matrix(orbit, theta, "z")
        effects = single_defect_effects(period, WIDTH, SPLIT, V)
        effect_seconds = time.perf_counter()-start
        exact, dephased, contributions = probability_rows(effects)
        start = time.perf_counter()
        reference, ref_info = coherent_reference(theta)
        reference_seconds = time.perf_counter()-start
        fft = phase_fft_reference(theta)
        error = float(np.max(np.abs(exact-reference)))
        fft_error = float(np.max(np.abs(exact-fft)))
        exp.check("P1", max(error, fft_error) < 1e-9
                  and max(ref_info["scratch_leakage"], ref_info["invalid_work_leakage"]) < 1e-18,
                  f"theta/pi={theta/np.pi:g}: gate error={error:.2e}, FFT={fft_error:.2e}")
        completeness = float(np.max(np.abs(effects.sum(axis=0)-np.eye(period))))
        minimum_eigenvalue = float(np.linalg.eigvalsh(effects).min())
        reconstruction = float(np.max(np.abs(exact-dephased-contributions.sum(axis=1))))
        exp.check("P2", completeness < 1e-10 and minimum_eigenvalue > -1e-12
                  and reconstruction < 1e-12,
                  f"complete={completeness:.2e}, min eig={minimum_eigenvalue:.2e}, pairs={reconstruction:.2e}")
        # A different clean input prevents the all-ones spectral rho from
        # hiding an orientation error in E_y or its imaginary coordinates.
        F = orbit_fourier(period)
        psi2 = F.conj().T[:, 1]  # work |a>, not |1>
        probe = np.einsum("i,yij,j->y", psi2.conj(), effects, psi2).real
        exp.check("P2", np.max(np.abs(probe-phase_fft_reference(theta, initial_label=BASE))) < 1e-11,
                  "different orbit input checks complex effect orientation")
        start = time.perf_counter()
        adjusted = scalar_mixture(period, theta)
        scalar_seconds = time.perf_counter()-start
        exp.check("P4", np.max(np.abs(adjusted-reference)) < 1e-9,
                  f"phase-adjusted scalar mixture error={np.max(np.abs(adjusted-reference)):.2e}")
        # The nontrivial commuting control is run at every same theta.
        commuting = single_defect_effects(period, WIDTH, SPLIT,
                                         defect_matrix(orbit, theta, "xxx"))
        pcomm, dcomm, _ = probability_rows(commuting)
        exp.check("P3", np.max(np.abs(pcomm-ideal)) < 1e-11
                  and np.max(np.abs(dcomm-ideal)) < 1e-11,
                  f"theta/pi={theta/np.pi:g}: commuting orbit rotation recovers ideal")
        if theta == 0:
            exp.check("P3", np.max(np.abs(exact-ideal)) < 1e-11
                      and np.max(np.abs(dephased-ideal)) < 1e-11,
                      "zero perturbation recovers original scalar baseline")
        dephasing_error = float(np.max(np.abs(exact-dephased)))
        ideal_error = float(np.max(np.abs(exact-ideal)))
        max_dephasing_error = max(max_dephasing_error, dephasing_error)
        max_ideal_error = max(max_ideal_error, ideal_error)
        response = coherence_response(effects)
        rows.append(dict(series="strength", theta=float(theta), theta_over_pi=float(theta/np.pi),
                         probabilities=exact.tolist(), dephased=dephased.tolist(),
                         adjusted_scalar=adjusted.tolist(), reference=reference.tolist(),
                         gate_error=error, fft_error=fft_error, dephasing_error=dephasing_error,
                         dephasing_tv=float(np.abs(exact-dephased).sum()/2), ideal_error=ideal_error,
                         pair_contributions=contributions.tolist(), response=response,
                         completeness_error=completeness, minimum_eigenvalue=minimum_eigenvalue,
                         effect_payload_bytes=effects.nbytes, effect_seconds=effect_seconds,
                         reference_seconds=reference_seconds, scalar_full_mixture_seconds=scalar_seconds,
                         **ref_info))
        exp.log(f"dephasing max error={dephasing_error:.6g}, TV={rows[-1]['dephasing_tv']:.6g}, "
                f"detectable pairs={response['detectable_pairs']}, response rank={response['rank']}")

    exp.fail_check("C1", max_dephasing_error > 1e-3, f"largest dephasing error={max_dephasing_error:.6g}")
    exp.fail_check("C2", max_ideal_error > 1e-3, f"largest ideal-model error={max_ideal_error:.6g}")
    # Moving a work-only unitary to the end cannot change exponent statistics.
    end, _ = coherent_reference(np.pi/2, split=WIDTH)
    middle = next(row for row in rows if row["theta_over_pi"] == .5)
    exp.fail_check("C3", np.max(np.abs(end-middle["probabilities"])) > 1e-3
                   and np.max(np.abs(end-ideal)) < 1e-9,
                   f"wrong-time-order error={np.max(np.abs(end-middle['probabilities'])):.6g}")
    # Separate fixed-location control: noncommuting globally on the orbit,
    # yet scalar on the actual early reachable states. No main-sweep confound.
    V = defect_matrix(orbit, np.pi/2, "z")
    U = np.diag(np.exp(2j*np.pi*np.arange(period)/period))
    commutator = float(np.max(np.abs(U @ V - V @ U)))
    early, _ = coherent_reference(np.pi/2, split=1)
    early_effects = single_defect_effects(period, WIDTH, 1, V)
    early_exact, early_dephased, _ = probability_rows(early_effects)
    exp.fail_check("C4", commutator > .1 and np.max(np.abs(early-ideal)) < 1e-9
                   and np.max(np.abs(early-early_exact)) < 1e-9,
                   f"noncommutator={commutator:.6g}, unchanged-output error={np.max(np.abs(early-ideal)):.2e}")
    rows.append(dict(series="location_control", split=1, theta_over_pi=.5,
                     commutator_max=commutator, ideal_error=float(np.max(np.abs(early-ideal))),
                     dephasing_error=float(np.max(np.abs(early_exact-early_dephased)))))
    comm_gate, _ = coherent_reference(np.pi/2, kind="xxx")
    exp.check("P3", np.max(np.abs(comm_gate-ideal)) < 1e-9,
              "nontrivial commuting rotation checked by existing compiled statevec")
    # Extended precision spectral arithmetic; same existing coherent engine,
    # with its original float64 gate constants, also runs with extended state.
    theta_ld = np.arccos(np.longdouble(-1))/2
    extended = single_defect_effects(period, WIDTH, SPLIT,
                                    defect_matrix(orbit, theta_ld, "z", np.longdouble),
                                    dtype=np.longdouble)
    extended_probs, _, _ = probability_rows(extended)
    extended_gate, _ = coherent_reference(theta_ld, dtype=np.clongdouble)
    precision_error = float(np.max(np.abs(extended_probs-middle["probabilities"])))
    extended_gate_error = float(np.max(np.abs(extended_probs-extended_gate)))
    exp.check("P3", precision_error < 1e-12 and extended_gate_error < 1e-9,
              f"extended spectral delta={precision_error:.2e}, extended-state gate error={extended_gate_error:.2e}")
    rows.append(dict(series="precision", float_mantissa=np.finfo(float).nmant,
                     extended_mantissa=np.finfo(np.longdouble).nmant,
                     spectral_delta=precision_error, extended_gate_error=extended_gate_error))
    exp.finish(report_path="out/spectral_defect.json", rows=rows,
               metadata=dict(N=N, a=BASE, width=WIDTH, split=SPLIT, period=period,
                             orbit=orbit.tolist(), orbit_setup_seconds=orbit_seconds,
                             compile_seconds=compile_seconds, validation_seconds=validation_seconds,
                             compiled_validation=action_stats(pairs), numpy=np.__version__,
                             task="perturbed clean-input exponent-output distribution",
                             precision="float64/complex128; extended-state check retains original gate constants",
                             cost="dense Q*r*r effects and full distributions are validation only; no speed benchmark"))


if __name__ == "__main__":
    main()
