"""Few-harmonic Fourier formulas for a repeated orbit block mixer.

DERIVED BEFORE MEASUREMENT.  Let r be divisible by b and apply the same
unitary W to every consecutive orbit block bm,...,bm+b-1 after L=2**s
ascending-power controls.  For p=l mod b, put

    A_p(k) = sum_q W[q,p] exp(2*pi*i*k*q/r)
    f_p(k) = exp(-2*pi*i*k*p/r) A_p(k)
    h_k(l) = exp(2*pi*i*k*l/r) f_{l mod b}(k).

The period-b multiplier has the positive-sign Fourier expansion
f_p(k)=sum_u c_u(k) exp(2*pi*i*u*p/b).  Therefore the conditional row is a
sum of b plane waves.  This test compares the resulting full joint p(k,y)
to the existing SparseOrbitPrefix implementation and its marginal to
lab.spectral.single_defect_effects.  It does not implement a generic
propagator or use a state-vector reference.

P1: with one seeded complex b=3 block held fixed, r=6,t=8 and s=0..8,
    the formula agrees with SparseOrbitPrefix for every joint (k,y), and
    the output marginal agrees with lab.spectral.
P2: c_u uses the stated positive-sign DFT and sum_u |c_u|^2=1; the phase
    norm is S_k=sum_p n_p |A_p(k)|^2, including incomplete early periods.
P3: the b-component envelope has T_k=L and its mean rejection cost over the
    actually sampled phase weights is at most b.  Pointwise C_k<=b is not
    asserted; any observed pointwise excess is recorded.
P4: identity W, r=9 with the same W, and s=t controls preserve the same
    formulas and normalization without changing the seeded block.

C1: deleting harmonic interference must fail.
C2: replacing final phase weights by uniform weights must fail.
C3: omitting inverse-QFT feedback from the early factors must fail.
C4: applying the unmodified b=3 formula to r=7 with an incomplete identity
    block must fail, guarding the b|r promise.

All dense allocations are preflight-capped with math.prod: r<=18, t<=8 and
the planned simultaneous payload is <=16 MiB.  The known orbit indices and
the block matrix are supplied inputs; no order-discovery or locality claim
is made.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_periodic_formulas
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects


MAX_DENSE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_DENSE_BYTES = 16 * 1024 * 1024
MAX_R = 18
MAX_T = 8
B = 3
SEED = 2030910


def guard(shape, dtype, label: str) -> int:
    """Return payload and reject before a dense allocation."""
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")
    return payload


def report_path(prefix: str = "periodic_formulas") -> Path:
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


def seeded_block(size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    guard((size, size), np.complex128, "seeded block")
    z = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    q *= np.where(np.abs(diagonal) > 0, np.conj(diagonal) / np.abs(diagonal), 1.)
    return q


def repeated_block(period: int, block_size: int, block: np.ndarray) -> np.ndarray:
    if period < 1 or period > MAX_R or block_size < 1 or period % block_size:
        raise ValueError("period must be a positive multiple of block size")
    if block.shape != (block_size, block_size):
        raise ValueError("block shape mismatch")
    if np.max(np.abs(block.conj().T @ block - np.eye(block_size))) > 2e-12:
        raise ValueError("block must be unitary")
    guard((period, period), np.complex128, "repeated block")
    result = np.eye(period, dtype=complex)
    for start in range(0, period, block_size):
        result[start:start + block_size, start:start + block_size] = block
    return result


def incomplete_block(period: int, block_size: int, block: np.ndarray) -> np.ndarray:
    """The deliberate b∤r control: complete blocks, identity remainder."""
    if period < 1 or period > MAX_R or block_size < 1 or block_size > period:
        raise ValueError("invalid incomplete-block dimensions")
    result = np.eye(period, dtype=complex)
    for start in range(0, period - block_size + 1, block_size):
        result[start:start + block_size, start:start + block_size] = block
    return result


def geom(length: int, frequency: float) -> complex:
    if length <= 0:
        return 0j
    ratio = np.exp(2j * np.pi * frequency)
    if abs(ratio - 1) < 2e-14:
        return complex(length)
    return (1 - ratio ** length) / (1 - ratio)


def periodic_coefficients(period: int, block_size: int, block: np.ndarray,
                          k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return A_p, f_p and positive-sign DFT c_u for one k."""
    p = np.arange(block_size)
    q = np.arange(block_size)
    phase = np.exp(2j * np.pi * k * q / period)
    a = block.T @ phase
    f = np.exp(-2j * np.pi * k * p / period) * a
    c = np.array([np.sum(f * np.exp(-2j * np.pi * u * p / block_size))
                  / block_size for u in range(block_size)])
    return a, f, c


def formula_joint(period: int, block_size: int, block: np.ndarray,
                  width: int, split: int, *, deleted_interference: bool = False,
                  uniform_phase_weights: bool = False,
                  feedback: bool = True) -> np.ndarray:
    """Candidate p(k,y), with no generic circuit propagation."""
    if period % block_size:
        # This function intentionally permits the invalid control so the
        # must-fail check measures the consequence of violating b|r.
        pass
    if not 1 <= period <= MAX_R or not 1 <= block_size <= MAX_R or not 0 <= split <= width <= MAX_T:
        raise ValueError("invalid width or split")
    Q, L, H = 1 << width, 1 << split, 1 << (width - split)
    planned = (guard((period, Q), np.float64, "joint formula")
               + guard((period, block_size), np.complex128, "coefficient table"))
    if planned > MAX_TOTAL_DENSE_BYTES:
        raise MemoryError("planned dense formula payload exceeds 16 MiB")
    output = np.zeros((period, Q), dtype=float)
    rows = np.empty((period, Q), dtype=float)
    phase_norm = np.empty(period, dtype=float)
    for k in range(period):
        a, _, c = periodic_coefficients(period, block_size, block, k)
        counts = np.array([sum(1 for l in range(L) if l % block_size == p)
                           for p in range(block_size)])
        phase_norm[k] = float(np.sum(counts * np.abs(a) ** 2))
        for y in range(Q):
            z, q = divmod(y, H)
            late = geom(H, L * (k / period - y / Q))
            early = np.empty(block_size, dtype=complex)
            for u in range(block_size):
                early_frequency = k / period + u / block_size - y / Q
                if not feedback:
                    # Deliberately omit q*l/Q feedback while retaining the
                    # late scalar phase; z is the low-register output.
                    early_frequency = k / period + u / block_size - z / L
                early[u] = geom(L, early_frequency)
            pieces = c * early
            if deleted_interference:
                magnitude = float(np.sum(np.abs(pieces) ** 2))
            else:
                magnitude = float(abs(np.sum(pieces)) ** 2)
            rows[k, y] = abs(late) ** 2 * magnitude / (Q * Q * period)
    if uniform_phase_weights:
        # Preserve each conditional row, but replace w_k by 1/r.
        weights = rows.sum(axis=1)
        for k in range(period):
            if weights[k] > 0:
                output[k] = rows[k] / weights[k] / period
        return output
    return rows


def prefix_joint(v_orbit: np.ndarray, period: int, block_size: int,
                 width: int, split: int) -> tuple[np.ndarray, SparseOrbitPrefix]:
    def column(l):
        values = v_orbit[:, l % period]
        keep = np.flatnonzero(values != 0)
        return keep.astype(np.int64), values[keep]

    prefix = SparseOrbitPrefix(period, split, column, block_size)
    Q = 1 << width
    guard((period, Q), np.float64, "prefix joint")
    result = np.zeros((period, Q), dtype=float)
    for k in range(period):
        for y in range(Q):
            result[k, y] = prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
    return result, prefix


def spectral_marginal(v_orbit: np.ndarray, period: int, width: int,
                      split: int) -> np.ndarray:
    fourier = orbit_fourier(period)
    v_eigen = fourier.conj().T @ v_orbit @ fourier
    effects = single_defect_effects(period, width, split, v_eigen,
                                    max_entries=1_000_000)
    return effects.sum(axis=(1, 2)).real / period


def evaluate_case(name: str, period: int, block: np.ndarray, v_orbit: np.ndarray,
                  width: int, split: int, exp: Experiment) -> dict:
    formula = formula_joint(period, B, block, width, split)
    reference, prefix = prefix_joint(v_orbit, period, B, width, split)
    spectral = spectral_marginal(v_orbit, period, width, split)
    marginal = formula.sum(axis=0)
    ref_marginal = reference.sum(axis=0)
    formula_error = float(np.max(np.abs(formula - reference)))
    spectral_error = float(np.max(np.abs(marginal - spectral)))
    prefix_marginal_error = float(np.max(np.abs(marginal - ref_marginal)))
    total_error = float(abs(formula.sum() - 1.))

    norm_errors = []
    c_norm_errors = []
    weights = []
    for k in range(period):
        a, f, c = periodic_coefficients(period, B, block, k)
        counts = np.array([sum(1 for l in range(1 << split) if l % B == p)
                           for p in range(B)])
        direct_s = float(np.sum(counts * np.abs(a) ** 2))
        # Prefix gives the same final-eigenphase norm independently.
        formula_s = float(reference[k].sum() * (1 << split) * period)
        norm_errors.append(abs(direct_s - formula_s))
        c_norm_errors.append(abs(float(np.sum(np.abs(c) ** 2)) - 1.))
        weights.append(float(reference[k].sum()))
    weights = np.asarray(weights)
    max_norm_error = float(max(norm_errors))
    max_c_norm_error = float(max(c_norm_errors))
    L = 1 << split
    nonzero = weights > 0
    cvals = np.zeros(period)
    envelope_error = 0.
    # Actually normalize the Fourier-component proposal and integrate its
    # acceptance, rather than obtaining the mean solely from w*(b/(r*w)).
    H, Q = 1 << (width-split), 1 << width
    late_q = min(1, H-1)
    for k in np.flatnonzero(nonzero):
        _, _, c = periodic_coefficients(period, B, block, int(k))
        pieces = np.array([[c[u]*geom(L, k/period+u/B-(late_q+H*z)/Q)/L
                            for u in range(B)] for z in range(L)])
        diagonal = np.sum(np.abs(pieces)**2, axis=1)
        proposal = diagonal/float(np.sum(np.abs(c)**2))
        acceptance = np.zeros(L)
        nz = diagonal > 0
        acceptance[nz] = np.abs(pieces[nz].sum(axis=1))**2/(B*diagonal[nz])
        accepted = proposal*acceptance
        mass = float(accepted.sum())
        cvals[k] = 1/mass
        row = prefix.phase_row(int(k))
        expected = np.abs(np.fft.fft(row*np.exp(-2j*np.pi*late_q*np.arange(L)/Q)))**2
        expected /= L*float(np.vdot(row, row).real)
        envelope_error = max(envelope_error, abs(proposal.sum()-1),
                             float(np.max(np.abs(accepted/mass-expected))))
    mean_cost = float(np.sum(weights[nonzero] * cvals[nonzero]))
    max_cost = float(np.max(cvals))

    exp.check("P1", formula_error < 3e-10 and spectral_error < 3e-10
              and prefix_marginal_error < 3e-10 and total_error < 3e-10,
              f"{name}: joint={formula_error:.2e}, spectral={spectral_error:.2e}, "
              f"prefix marginal={prefix_marginal_error:.2e}, total={total_error:.2e}")
    exp.check("P2", max_norm_error < 3e-10 and max_c_norm_error < 3e-12,
              f"{name}: S={max_norm_error:.2e}, sum|c|²={max_c_norm_error:.2e}")
    exp.check("P3", mean_cost <= B + 3e-10 and envelope_error < 3e-10,
              f"{name}: weighted rejection mean={mean_cost:.9g} <= {B}, "
              f"pointwise max={max_cost:.9g}, proposal/rejection={envelope_error:.2e}, "
              f"zero phases={np.flatnonzero(~nonzero).tolist()}")
    return dict(name=name, period=period, width=width, split=split, L=L,
                joint_error=formula_error, spectral_error=spectral_error,
                prefix_marginal_error=prefix_marginal_error, total_error=total_error,
                max_phase_norm_error=max_norm_error, max_c_norm_error=max_c_norm_error,
                phase_weights=weights.tolist(), weighted_rejection_mean=mean_cost,
                pointwise_rejection_max=max_cost,
                envelope_error=float(envelope_error),
                pointwise_exceeds_b=bool(max_cost > B + 1e-8),
                zero_weight_indices=np.flatnonzero(~nonzero).tolist(),
                prefix_stats=prefix.stats())


def main() -> None:
    exp = Experiment("periodic_formulas", doc=__doc__)
    exp.predict("P1", "few-harmonic p(k,y) agrees with SparseOrbitPrefix and lab.spectral")
    exp.predict("P2", "positive-sign coefficients are normalized and S_k uses residue counts")
    exp.predict("P3", "averaged early rejection cost is at most b, without a pointwise claim")
    exp.predict("P4", "identity, r=9 and s=t controls preserve the formula")
    exp.must_fail("C1", "deleted harmonic interference changes the output")
    exp.must_fail("C2", "uniform final phase weights are not the measured weights")
    exp.must_fail("C3", "omitting inverse-QFT feedback changes the output")
    exp.must_fail("C4", "the b|r formula is invalid for an incomplete identity block")

    start = time.perf_counter()
    block = seeded_block(B, SEED)
    identity = np.eye(B, dtype=complex)
    rows = []
    # Main fixed r=6,t=8,W series; only insertion split varies by one.
    v6 = repeated_block(6, B, block)
    for split in range(0, 9):
        rows.append(evaluate_case(f"seeded_r6_t8_s{split}", 6, block, v6, 8, split, exp))

    # Fixed block and small controls: no fresh random matrix is drawn.
    rows.append(evaluate_case("seeded_r9_t8_s4", 9, block,
                               repeated_block(9, B, block), 8, 4, exp))
    rows.append(evaluate_case("identity_r6_t8_s4", 6, identity,
                               repeated_block(6, B, identity), 8, 4, exp))
    rows.append(evaluate_case("seeded_r6_t8_st", 6, block, v6, 8, 8, exp))
    control_rows = [rows[-3], rows[-2], rows[-1]]
    exp.check("P4", all(row["joint_error"] < 3e-10 and
                         row["spectral_error"] < 3e-10 for row in control_rows),
              "identity/r=9/s=t controls retain formula and spectral agreement")

    # Controls are evaluated against the nontrivial main row and retain full
    # coherent sums until the deliberately deleted-interference branch.
    control = formula_joint(6, B, block, 8, 4)
    deleted = formula_joint(6, B, block, 8, 4, deleted_interference=True)
    uniform = formula_joint(6, B, block, 8, 4, uniform_phase_weights=True)
    no_feedback = formula_joint(6, B, block, 8, 4, feedback=False)
    correct_marginal = control.sum(axis=0)
    deleted_error = float(np.max(np.abs(deleted.sum(axis=0) - correct_marginal)))
    uniform_error = float(np.max(np.abs(uniform.sum(axis=0) - correct_marginal)))
    no_feedback_error = float(np.max(np.abs(no_feedback.sum(axis=0) - correct_marginal)))
    exp.fail_check("C1", deleted_error > 1e-5,
                   f"deleted-interference max error={deleted_error:.9g}")
    exp.fail_check("C2", uniform_error > 1e-5,
                   f"uniform-phase max error={uniform_error:.9g}")
    exp.fail_check("C3", no_feedback_error > 1e-5,
                   f"missing-feedback max error={no_feedback_error:.9g}")

    # Deliberate b∤r control.  The physical matrix has identity on the final
    # incomplete block, while the formula still assumes the same W there.
    bad_period = 7
    bad_v = incomplete_block(bad_period, B, block)
    bad_exact, _ = prefix_joint(bad_v, bad_period, B, 8, 4)
    bad_formula = formula_joint(bad_period, B, block, 8, 4)
    bad_error = float(np.max(np.abs(bad_formula.sum(axis=0) - bad_exact.sum(axis=0)))
                      )
    exp.fail_check("C4", bad_error > 1e-5,
                   f"r=7 incomplete identity block max marginal error={bad_error:.9g}")

    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(
        seeded_block_size=B, seeded_block_seed=SEED, main_series="r=6,t=8,s=0..8",
        fixed_control="r=9,t=8,s=4", dense_budget_bytes=MAX_DENSE_BYTES,
        max_total_dense_bytes=MAX_TOTAL_DENSE_BYTES, max_r=MAX_R, max_t=MAX_T,
        deleted_interference_error=deleted_error,
        uniform_phase_error=uniform_error,
        missing_feedback_error=no_feedback_error,
        nondividing_period_error=bad_error,
        phase_cost_statement="mean over actual w_k <= b; pointwise C_k may exceed b",
        known_order_and_orbit_indices=True, generic_propagator=False,
        elapsed_seconds=elapsed, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("periodic_formulas_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
