"""Finite-rank Fourier formulas for a localized orbit defect.

DERIVED BEFORE MEASUREMENT.  In an orbit basis of period r, let V be the
identity outside a known set S of d indices.  For the ascending-power circuit
U^0,...,U^(L-1), V, U^L,...,U^(Q-L), write

    h_k(l) = exp(2*pi*i*k*l/r) + sum_p delta_p(k) 1[l mod r=p],
    delta_p(k) = sum_q (V[q,p]-1[q=p]) exp(2*pi*i*k*q/r).

The conditional row norm is the residue-count expression

    S_k = L + sum_p n_p (|exp(2*pi*i*k*p/r)+delta_p(k)|^2-1).

The output amplitude is a late geometric sum multiplying one plane-wave
geometric sum plus d truncated progression sums.  This experiment checks that
formula against tiny direct sums, lab.spectral, SparseOrbitPrefix, and the
existing original-order full coherent mixer for N=7,a=3.

P1: geometric plane/progression amplitudes equal direct finite sums and the
    localized norm formula agrees with the direct row norm, including L<r,
    L>r, and a zero-weight phase.
P2: formula marginals agree with lab.spectral and SparseOrbitPrefix for the
    physical r=6 mixer and seeded complex localized 2x2/3x3 blocks.
P3: the representative original ascending-order coherent gate circuit agrees
    with the formula at every insertion s=1,...,7 for r=6,t=8,theta=pi/2.
P4: the independently computed component-rejection envelope has weighted mean
    at most 5m, with zero-weight rows omitted; no pointwise C_k<=5m claim is
    made.

C1: deleting coherent cross terms must fail.
C2: using floor(L/r) for every progression length must fail.
C3: omitting inverse-QFT feedback from the early factor must fail.

Only the insertion index varies in the main r=6 sweep. Other rows are fixed
controls for periods, angles, and localized seeded blocks. This is a bounded
formula/reference test, not a general propagator or a sampler claim.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_localized_formulas
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
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects
from modexp import ModExp
import statevec


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, WIDTH = 7, 3, 8
PHYSICAL_ORBIT = np.array([1, 3, 2, 6, 4, 5], dtype=np.int64)


def guard(shape, dtype, label: str) -> None:
    """Reject a dense allocation before allocating it."""
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def alloc(shape, dtype, label: str) -> np.ndarray:
    guard(shape, dtype, label)
    return np.empty(shape, dtype=dtype)


def report_path(prefix: str = "localized_formulas") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def orbit_fourier(period: int) -> np.ndarray:
    guard((period, period), np.complex128, "orbit Fourier")
    j = np.arange(period, dtype=float)
    return np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)


def localized_unitary(period: int, support: tuple[int, ...], seed: int) -> np.ndarray:
    """Seeded complex unitary block, identity on the complement."""
    if not support or len(set(support)) != len(support) or min(support) < 0 or max(support) >= period:
        raise ValueError("invalid localized support")
    rng = np.random.default_rng(seed)
    guard((len(support), len(support)), np.complex128, "localized QR block")
    z = rng.normal(size=(len(support), len(support))) + 1j * rng.normal(size=(len(support), len(support)))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    q *= np.where(np.abs(diagonal) > 0, np.conj(diagonal) / np.abs(diagonal), 1.)
    guard((period, period), np.complex128, "localized defect")
    result = np.eye(period, dtype=complex)
    result[np.ix_(support, support)] = q
    return result


def zero_weight_unitary(period: int = 3) -> np.ndarray:
    """A 2x2 Hadamard block gives S_0=0 when L=1."""
    guard((period, period), np.complex128, "zero-weight defect")
    result = np.eye(period, dtype=complex)
    guard((2, 2), np.complex128, "zero-weight block")
    block = np.array([[1., 1.], [-1., 1.]], dtype=complex) / np.sqrt(2.)
    result[np.ix_((0, 1), (0, 1))] = block
    return result


def physical_work_matrix(theta: float) -> np.ndarray:
    """The existing N=7,a=3 mixer, restricted to orbit indices."""
    guard((6, 6), np.complex128, "physical defect")
    result = np.eye(6, dtype=complex)
    i, j = 0, 5                         # labels 1 and 5 in PHYSICAL_ORBIT
    result[i, i] = result[j, j] = np.cos(theta / 2)
    result[i, j] = result[j, i] = -1j * np.sin(theta / 2)
    return result


def geom_sum(length: int, first_phase: float, step_phase: float) -> complex:
    """sum_j exp(2*pi*i*(first_phase+step_phase*j)), stably at ratio 1."""
    if length <= 0:
        return 0j
    first = np.exp(2j * np.pi * first_phase)
    ratio = np.exp(2j * np.pi * step_phase)
    if abs(ratio - 1) < 2e-14:
        return complex(length) * first
    return first * (1 - ratio ** length) / (1 - ratio)


def deltas(v_orbit: np.ndarray, period: int) -> tuple[tuple[int, ...], np.ndarray]:
    support = tuple(int(p) for p in range(period)
                    if np.max(np.abs(v_orbit[:, p] -
                                     np.array([1. if q == p else 0. for q in range(period)], dtype=complex))) != 0)
    if not support:
        return support, np.zeros((period, 0), dtype=complex)
    delta = alloc((period, len(support)), np.complex128, "localized deltas")
    for k in range(period):
        for ip, p in enumerate(support):
            delta[k, ip] = sum((v_orbit[q, p] - (1. if q == p else 0.))
                               * np.exp(2j * np.pi * k * q / period) for q in range(period))
    return support, delta


def norm_rows(v_orbit: np.ndarray, period: int, length: int):
    support, delta = deltas(v_orbit, period)
    counts = {p: sum(1 for l in range(length) if l % period == p) for p in support}
    rows = alloc((period, length), np.complex128, "conditional rows")
    norms = alloc((period,), np.float64, "row norms")
    predicted = alloc((period,), np.float64, "residue-count norms")
    for k in range(period):
        for l in range(length):
            value = np.exp(2j * np.pi * k * l / period)
            if l % period in support:
                value += delta[k, support.index(l % period)]
            rows[k, l] = value
        norms[k] = float(np.vdot(rows[k], rows[k]).real)
        predicted[k] = length + sum(counts[p] * (
            abs(np.exp(2j * np.pi * k * p / period) + delta[k, ip]) ** 2 - 1)
            for ip, p in enumerate(support))
    return rows, norms, predicted, support, delta, counts


def formula_joint(v_orbit: np.ndarray, period: int, width: int, split: int,
                  *, feedback: bool = True, wrong_progression: bool = False,
                  deleted_interference: bool = False) -> np.ndarray:
    """Return joint p(k,y) using geometric plane/progression terms."""
    qsize, length = 1 << width, 1 << split
    height = qsize // length
    support, delta = deltas(v_orbit, period)
    guard((period, qsize), np.float64, "formula joint")
    result = alloc((period, qsize), np.float64, "formula joint")
    for k in range(period):
        for y in range(qsize):
            frequency = k / period - y / qsize
            late = geom_sum(height, 0., length * frequency)
            plane = geom_sum(length, 0., frequency)
            pieces = [plane]
            for ip, p in enumerate(support):
                count = sum(1 for l in range(length) if l % period == p)
                if wrong_progression:
                    count = length // period
                progression = geom_sum(count, -y * p / qsize, -y * period / qsize)
                pieces.append(delta[k, ip] * progression)
            if deleted_interference:
                magnitude = abs(pieces[0]) ** 2 + sum(abs(x) ** 2 for x in pieces[1:])
            else:
                magnitude = abs(sum(pieces)) ** 2
            # Missing feedback means the early inverse-QFT kernel uses z/L,
            # while q remains the late scalar output.  This branch is kept
            # explicit so the control cannot accidentally test the same sum.
            if not feedback:
                z = y // height
                q = y % height
                frequency_early = -z / length
                plane = geom_sum(length, 0., k / period + frequency_early)
                pieces = [plane]
                for ip, p in enumerate(support):
                    count = sum(1 for l in range(length) if l % period == p)
                    progression = geom_sum(count, frequency_early * p,
                                           frequency_early * period)
                    pieces.append(delta[k, ip] * progression)
                magnitude = (abs(pieces[0]) ** 2 + sum(abs(x) ** 2 for x in pieces[1:])
                             if deleted_interference else abs(sum(pieces)) ** 2)
                late = geom_sum(height, 0., length * (k / period - q / height / length))
            result[k, y] = abs(late) ** 2 * magnitude / (qsize * qsize * period)
    return result


def direct_joint(v_orbit: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    """Independent tiny finite sums, used only as the initial formula check."""
    qsize, length = 1 << width, 1 << split
    height = qsize // length
    # Independent of the plane-plus-delta formula: multiply actual V columns
    # by final eigenbras. Both initial and final exponent normalizations cost
    # sqrt(Q), hence the probability denominator Q^2*r.
    guard((period, length), np.complex128, "direct conditional row")
    rows = np.sqrt(period) * orbit_fourier(period).conj().T @ v_orbit[:, np.arange(length) % period]
    result = alloc((period, qsize), np.float64, "direct joint")
    for k in range(period):
        for y in range(qsize):
            amplitude = sum(rows[k, l] * np.exp(2j * np.pi * k * length * h / period)
                            * np.exp(-2j * np.pi * y * (l + length * h) / qsize)
                            for l in range(length) for h in range(height))
            result[k, y] = abs(amplitude) ** 2 / (qsize * qsize * period)
    return result


def spectral_marginal(v_orbit: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    fourier = orbit_fourier(period)
    v_eigen = fourier.conj().T @ v_orbit @ fourier
    effects = single_defect_effects(period, width, split, v_eigen,
                                    max_entries=1_000_000)
    return effects.sum(axis=(1, 2)).real / period


def sparse_oracle(v_orbit: np.ndarray, period: int, length: int):
    def column(l):
        values = v_orbit[:, l % period]
        keep = np.flatnonzero(values != 0)
        return keep.astype(np.int64), values[keep]
    return column


def prefix_marginal(v_orbit: np.ndarray, period: int, width: int, split: int) -> tuple[np.ndarray, dict]:
    _, _, _, support, _, _ = norm_rows(v_orbit, period, 1 << split)
    max_terms = max(1, len(support))
    prefix = SparseOrbitPrefix(period, split,
                               sparse_oracle(v_orbit, period, 1 << split), max_terms)
    guard((1 << width,), np.float64, "prefix marginal")
    result = np.zeros(1 << width, dtype=float)
    phase_error = 0.
    direct = direct_joint(v_orbit, period, width, split)
    for k in range(period):
        expected = float(np.sum(direct[k]))
        phase_error = max(phase_error, abs(prefix.phase_probability(k) - expected))
        for y in range(1 << width):
            result[y] += prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
    return result, dict(phase_error=phase_error, stats=prefix.stats())


def bound_audit(v_orbit: np.ndarray, period: int, length: int):
    rows, norms, _, support, delta, counts = norm_rows(v_orbit, period, length)
    m = len(support) + 1
    tvals = np.array([length + sum(counts[p] * abs(delta[k, ip]) ** 2
                                  for ip, p in enumerate(support)) for k in range(period)])
    cvals = np.zeros(period)
    nonzero = norms > 0
    cvals[nonzero] = m * tvals[nonzero] / norms[nonzero]
    weights = norms / (length * period)
    weighted_mean = float(np.sum(weights[nonzero] * cvals[nonzero]))
    parseval_mean = float(np.sum(tvals) / (length * period))
    return dict(m=m, T=tvals.tolist(), S=norms.tolist(), weights=weights.tolist(),
                C=cvals.tolist(), zero_weight_indices=np.flatnonzero(~nonzero).tolist(),
                weighted_mean=weighted_mean, parseval_T_mean=parseval_mean,
                bound=5 * m, max_C=float(np.max(cvals)),
                independent_identity_error=abs(weighted_mean - m * float(np.sum(tvals[nonzero])) / (length * period)))


def mixing_terms(x):
    x0, x1, x2 = x
    xm = 1 << x2
    z0, z1 = 1 << x0, 1 << x1
    return ((xm, 0, +1), (xm, z1, +1), (xm, z0, -1), (xm, z0 | z1, -1))


def mixing_circuit(theta: float, x: list[int]) -> Circuit:
    qc = Circuit(max(x) + 1)
    for xmask, zmask, sign in mixing_terms(x):
        qc.rot((xmask, zmask), sign * theta / 4)
    return qc


def coherent_reference(theta: float, split: int) -> tuple[np.ndarray, dict]:
    """Existing gate-level arithmetic in its original ascending order."""
    me = ModExp(N, BASE, n_exp=WIDTH)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    for i, q in enumerate(me.exp):
        qc.extend(me.u_a(q, pow(BASE, 1 << i, N)))
        if i + 1 == split:
            qc.extend(mixing_circuit(theta, me.x))
    qc.qft(me.exp, inverse=True)
    shape = (1 << me.n_qubits,)
    guard(shape, np.complex128, "coherent state vector")
    psi = np.zeros(shape[0], dtype=complex)
    psi[0] = 1.
    psi = statevec.run(qc, psi)
    guard((1 << WIDTH, 1 << me.n), np.float64, "coherent probability view")
    probabilities = np.sum(np.abs(psi.reshape(1 << WIDTH, -1)) ** 2, axis=1)
    return probabilities, dict(qubits=qc.n, rotations=len(qc.gates), vector_payload_bytes=psi.nbytes,
                               leakage=float(abs(np.vdot(psi, psi).real - 1.)))


def evaluate_case(name: str, v_orbit: np.ndarray, period: int, width: int, split: int,
                  exp: Experiment, *, physical_theta: float | None = None):
    length = 1 << split
    direct = direct_joint(v_orbit, period, width, split)
    formula = formula_joint(v_orbit, period, width, split)
    rows, norms, predicted_norms, support, delta, counts = norm_rows(v_orbit, period, length)
    spectral = spectral_marginal(v_orbit, period, width, split)
    prefix, prefix_info = prefix_marginal(v_orbit, period, width, split)
    direct_marginal = direct.sum(axis=0)
    formula_marginal = formula.sum(axis=0)
    formula_error = float(np.max(np.abs(formula - direct)))
    norm_error = float(np.max(np.abs(norms - predicted_norms)))
    spectral_error = float(np.max(np.abs(formula_marginal - spectral)))
    prefix_error = float(np.max(np.abs(formula_marginal - prefix)))
    direct_total_error = float(abs(direct.sum() - 1.))
    exp.check("P1", formula_error < 2e-10 and norm_error < 2e-11 and direct_total_error < 2e-11,
              f"{name}: geometric/direct={formula_error:.2e}, norm={norm_error:.2e}, total={direct_total_error:.2e}")
    exp.check("P2", spectral_error < 2e-10 and prefix_error < 2e-10 and prefix_info["phase_error"] < 2e-10,
              f"{name}: spectral={spectral_error:.2e}, prefix={prefix_error:.2e}, phase={prefix_info['phase_error']:.2e}")
    bound = bound_audit(v_orbit, period, length)
    bound_ok = bound["weighted_mean"] <= bound["bound"] + 2e-10 and bound["independent_identity_error"] < 2e-10
    exp.check("P4", bound_ok,
              f"{name}: weighted C={bound['weighted_mean']:.6g} <= {bound['bound']}, "
              f"zero={bound['zero_weight_indices']}, identity={bound['independent_identity_error']:.2e}")
    if name == "zero_weight_d2":
        exp.check("P4", 0 in bound["zero_weight_indices"], "constructed phase k=0 really has zero weight")
    gate_error = None
    gate_info = None
    if physical_theta is not None:
        gate, gate_info = coherent_reference(physical_theta, split)
        gate_error = float(np.max(np.abs(gate - formula_marginal)))
        exp.check("P3", gate_error < 2e-8 and gate_info["leakage"] < 2e-12,
                  f"{name}: original-order gate={gate_error:.2e}, leakage={gate_info['leakage']:.2e}")
    return dict(name=name, period=period, width=width, split=split, L=length,
                support=list(support), counts={str(k): v for k, v in counts.items()},
                formula_error=formula_error, norm_error=norm_error,
                spectral_error=spectral_error, prefix_error=prefix_error,
                prefix_phase_error=prefix_info["phase_error"], direct_total_error=direct_total_error,
                formula_marginal=formula_marginal.tolist(), spectral=spectral.tolist(),
                bound=bound, gate_error=gate_error, gate_info=gate_info,
                max_dense_bytes=MAX_DENSE_BYTES)


def main() -> None:
    exp = Experiment("localized_formulas", doc=__doc__)
    exp.predict("P1", "geometric plane/progression sums and residue-count norms match direct finite sums")
    exp.predict("P2", "localized formula matches lab.spectral and SparseOrbitPrefix")
    exp.predict("P3", "physical r=6 original-order coherent circuit agrees at every insertion")
    exp.predict("P4", "weighted component-rejection mean obeys 5m via Parseval, including zero weights")
    exp.must_fail("C1", "deleting coherent cross terms")
    exp.must_fail("C2", "using an incorrect floor progression length")
    exp.must_fail("C3", "omitting inverse-QFT early feedback")

    start = time.perf_counter()
    rows = []
    max_deleted = max_wrong_progression = max_no_feedback = 0.

    # Main one-parameter series: only insertion s changes.
    for split in range(1, 8):
        v = physical_work_matrix(np.pi / 2)
        rows.append(evaluate_case(f"physical_r6_t8_s{split}", v, 6, WIDTH, split, exp,
                                  physical_theta=np.pi / 2))

    # Fixed controls for L<r and L>r, angle signs/zero, and seeded localized blocks.
    controls = [
        ("abstract_r3_Llt", localized_unitary(3, (0, 1), 3103), 3, 5, 1, None),
        ("abstract_r3_Lgt", localized_unitary(3, (0, 1), 3103), 3, 5, 3, None),
        ("abstract_r5_Llt", localized_unitary(5, (0, 1), 3505), 5, 5, 1, None),
        ("abstract_r5_Lgt", localized_unitary(5, (0, 1), 3505), 5, 5, 3, None),
        ("angle_minus_pi", physical_work_matrix(-np.pi), 6, 5, 3, None),
        ("angle_minus_half_pi", physical_work_matrix(-np.pi / 2), 6, 5, 3, None),
        ("angle_zero", physical_work_matrix(0.), 6, 5, 3, None),
        ("angle_plus_pi", physical_work_matrix(np.pi), 6, 5, 3, None),
        ("seeded_complex_d2", localized_unitary(5, (1, 4), 4202), 5, 5, 3, None),
        ("seeded_complex_d3", localized_unitary(5, (0, 2, 4), 4303), 5, 5, 3, None),
        ("zero_weight_d2", zero_weight_unitary(3), 3, 3, 0, None),
    ]
    for name, v, period, width, split, physical_theta in controls:
        rows.append(evaluate_case(name, v, period, width, split, exp,
                                  physical_theta=physical_theta))

    # Controls are evaluated only after all positive formula checks and only
    # on the main nontrivial physical row, so they cannot be vacuous.
    control_v = physical_work_matrix(np.pi / 2)
    control_period, control_width, control_split = 6, WIDTH, 3
    correct = formula_joint(control_v, control_period, control_width, control_split)
    deleted = formula_joint(control_v, control_period, control_width, control_split,
                            deleted_interference=True)
    wrong = formula_joint(control_v, control_period, control_width, control_split,
                          wrong_progression=True)
    no_feedback = formula_joint(control_v, control_period, control_width, control_split,
                                feedback=False)
    correct_m = correct.sum(axis=0)
    max_deleted = float(np.max(np.abs(deleted.sum(axis=0) - correct_m)))
    max_wrong_progression = float(np.max(np.abs(wrong.sum(axis=0) - correct_m)))
    max_no_feedback = float(np.max(np.abs(no_feedback.sum(axis=0) - correct_m)))
    exp.fail_check("C1", max_deleted > 1e-5, f"deleted-interference max error={max_deleted:.6g}")
    exp.fail_check("C2", max_wrong_progression > 1e-5,
                   f"wrong progression-length max error={max_wrong_progression:.6g}")
    exp.fail_check("C3", max_no_feedback > 1e-5,
                   f"missing-feedback max error={max_no_feedback:.6g}")

    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(main_series="r=6,t=8,theta=pi/2,s=1..7", controls=len(controls),
                             max_deleted_interference_error=max_deleted,
                             max_wrong_progression_error=max_wrong_progression,
                             max_missing_feedback_error=max_no_feedback,
                             dense_budget_bytes=MAX_DENSE_BYTES,
                             bound_statement="weighted mean of C_k is <=5m; no pointwise bound asserted",
                             setup="known order/orbit indices; formula and references are finite validation routes",
                             elapsed_seconds=elapsed, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        # Never overwrite an earlier failure artifact.  This is useful when a
        # later rerun is attempted after a failed prediction/control.
        failure = report_path("localized_formulas_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
