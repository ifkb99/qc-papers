"""Bounded scaling audit for the sparse orbit-prefix sampler.

DERIVED BEFORE MEASUREMENT.  Fix split s=2 (L=4) and use the known-index
two-level mixer

  V|0> = cos(theta/2)|0> - i sin(theta/2)|r-1>,
  V|r-1> = cos(theta/2)|r-1> - i sin(theta/2)|0>,

with all other orbit-index columns unchanged, at theta=pi/2.  The oracle is
given orbit indices, not computational labels.  The small series uses r=6
and widths 2..8,16,32,63; the abstract series uses width 32 and
r=6,30,308,1000000007.  The latter is a cyclic-orbit oracle with known period
and support indices, not compiled modular arithmetic or a physical-qubit
simulation.

P1: for r=6 and widths <=8, forced-joint enumeration from SparseOrbitPrefix
    agrees with lab.spectral.single_defect_effects; width 4 also agrees with
    an independent computational-basis control/work Fourier calculation.
P2: every one of 16 fixed-seed samples per row reports the same joint latent
    probability as forced_joint(k,y), and wide samples agree with an
    independent integer-reduced geometric tail plus direct length-4 DFT.
P3: setup and sample accounting shows only the four early amplitudes and
    sparse oracle calls are retained; all measured probabilities are finite.
P4: the API rejects a computational label outside the supplied orbit-index
    range, catching labels/indices confusion.
C1: replacing the weighted final-k draw by uniform k must fail on the small
    width-4 row (a nonuniform final weight check precedes this control).
C2: forgetting the early feedback phase exp(-2*pi*i*q*l/Q) must fail on that
    same row.

The sparse rejection envelope is an expected-at-most-D statement in exact
arithmetic, not a deterministic bound on 16 random attempts.  Order discovery,
finding orbit indices, and oracle setup remain external costs; no speedup or
physical-memory claim is made.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_prefix_scaling
"""
from __future__ import annotations

import time

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects
import statevec


SPLIT = 2
WIDTHS = (2, 3, 4, 5, 6, 7, 8, 16, 32, 63)
ABSTRACT_PERIODS = (6, 30, 308, 1_000_000_007)
ABSTRACT_WIDTH = 32
THETA = np.pi / 2
SAMPLES = 16
MAX_REFERENCE_ENTRIES = 1_000_000


def mixer_matrix(period: int) -> np.ndarray:
    if period != 6:
        raise ValueError("dense mixer reference is capped to the r=6 validation case")
    c, s = np.cos(THETA / 2), np.sin(THETA / 2)
    matrix = np.eye(period, dtype=complex)
    matrix[0, 0] = matrix[-1, -1] = c
    matrix[-1, 0] = matrix[0, -1] = -1j * s
    return matrix


def mixer_oracle(period: int):
    """Return known orbit-index amplitudes for V|a^(l mod r)>."""
    c, s = np.cos(THETA / 2), np.sin(THETA / 2)

    def column(l):
        j = l % period
        if j == 0:
            return np.array([0, period - 1], dtype=np.int64), np.array([c, -1j * s])
        if j == period - 1:
            return np.array([period - 1, 0], dtype=np.int64), np.array([c, -1j * s])
        return np.array([j], dtype=np.int64), np.array([1. + 0j])

    return column


def orbit_fourier(period: int) -> np.ndarray:
    j = np.arange(period)
    return np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)


def spectral_outputs(period: int, width: int, v_orbit: np.ndarray) -> np.ndarray:
    if period != 6 or width > 8 or (1 << width) * period * period > MAX_REFERENCE_ENTRIES:
        raise ValueError("spectral reference allocation exceeds the explicit r=6,width<=8 cap")
    fourier = orbit_fourier(period)
    v_eigen = fourier.conj().T @ v_orbit @ fourier
    effects = single_defect_effects(period, width, SPLIT, v_eigen)
    return effects.sum(axis=(1, 2)).real / period


def clean_fourier_reference(period: int, width: int, v_orbit: np.ndarray) -> np.ndarray:
    """Independent control/work basis replay followed by control-only FFT."""
    if period != 6 or width > 8 or (1 << width) * period > MAX_REFERENCE_ENTRIES:
        raise ValueError("clean Fourier reference allocation exceeds the explicit r=6,width<=8 cap")
    qsize, length, height = 1 << width, 1 << SPLIT, 1 << (width - SPLIT)
    joint = np.zeros((qsize, period), dtype=complex)
    for e in range(qsize):
        l, h = e % length, e // length
        for j in range(period):
            joint[e, (j + length * h) % period] += v_orbit[j, l % period] / np.sqrt(qsize)
    transformed = np.fft.fft(joint, axis=0) / np.sqrt(qsize)
    return np.sum(np.abs(transformed) ** 2, axis=1)


def early_row(period: int, k: int, length: int):
    """Independent h_k(l), with exact integer phase reduction."""
    c, s = np.cos(THETA / 2), np.sin(THETA / 2)
    row = np.empty(length, dtype=complex)
    for l in range(length):
        j = l % period
        if j == 0:
            ids, amps = (0, period - 1), (c, -1j * s)
        elif j == period - 1:
            ids, amps = (period - 1, 0), (c, -1j * s)
        else:
            ids, amps = (j,), (1. + 0j,)
        row[l] = sum(a * np.exp(2j * np.pi * (((k * q) % period) / period))
                     for q, a in zip(ids, amps))
    return row


def geometric_tail_probability(period: int, width: int, k: int, q: int, length: int) -> float:
    """Finite H-term geometric norm using integer-reduced phase numerator."""
    height = 1 << (width - SPLIT)
    denominator = period * height
    numerator = (k * length * height - q * period) % denominator
    numerator = min(numerator, denominator - numerator)
    if numerator == 0:
        return 1.
    # Reduce the numerator of sin(pi*d/r) independently to avoid huge angles.
    dmod = numerator % (2 * period)
    dmod = min(dmod, 2 * period - dmod)
    top = np.sin(np.pi * (dmod / period))
    bottom = np.sin(np.pi * (numerator / denominator))
    return float((top / (height * bottom)) ** 2)


def independent_joint(period: int, width: int, k: int, y: int, *, feedback: bool = True) -> float:
    """Independent p(k,y): integer-reduced late sum, direct 4-point DFT."""
    length, height = 1 << SPLIT, 1 << (width - SPLIT)
    q, z = y % height, y // height
    row = early_row(period, k, length)
    norm2 = float(np.vdot(row, row).real)
    phase_probability = norm2 / (length * period)
    early_amplitude = 0j
    for l in range(length):
        phase = 1.
        if feedback:
            phase = np.exp(-2j * np.pi * (((q * l) % (length * height)) / (length * height)))
        early_amplitude += row[l] / np.sqrt(norm2) * phase * np.exp(-2j * np.pi * l * z / length)
    early_probability = float(abs(early_amplitude / np.sqrt(length)) ** 2)
    return phase_probability * geometric_tail_probability(period, width, k, q, length) * early_probability


def independent_marginal(period: int, width: int, *, feedback: bool = True) -> np.ndarray:
    if period != 6 or width > 8 or (1 << width) > MAX_REFERENCE_ENTRIES:
        raise ValueError("independent marginal allocation exceeds the explicit r=6,width<=8 cap")
    return np.array([sum(independent_joint(period, width, k, y, feedback=feedback)
                         for k in range(period)) for y in range(1 << width)])


def one_row(exp: Experiment, name: str, period: int, width: int, *, full_check: bool,
            reference_width4: bool = False):
    length = 1 << SPLIT
    # Abstract large-period rows are oracle-only; never allocate a dense
    # period-by-period matrix for them.
    if full_check and (period != 6 or width > 8):
        raise ValueError("full reference requested outside the explicit r=6,width<=8 cap")
    v_orbit = mixer_matrix(period) if full_check else None
    max_terms = 2
    started = time.perf_counter()
    prefix = SparseOrbitPrefix(period, SPLIT, mixer_oracle(period), max_terms)
    setup_seconds = time.perf_counter() - started
    setup_stats = prefix.stats().copy()
    rng = np.random.default_rng(0x510000 + period + width)
    samples = []
    max_forced_error = 0.
    max_independent_abs_error = 0.
    max_independent_rel_error = 0.
    positive_reference_count = 0
    zero_reference_count = 0
    minimum_positive_reference = None
    for sample_id in range(SAMPLES):
        before = prefix.stats().copy()
        result = prefix.sample(width, rng)
        after = prefix.stats().copy()
        forced = prefix.forced_joint(width, result["eigenphase"], result["output"])
        forced_error = abs(result["joint_latent_output_probability"]
                           - forced["joint_latent_output_probability"])
        independent = independent_joint(period, width, result["eigenphase"], result["output"])
        independent_abs_error = abs(forced["joint_latent_output_probability"] - independent)
        independent_positive = bool(independent > 0.)
        independent_rel_error = (independent_abs_error / independent
                                 if independent_positive else None)
        if independent_rel_error is not None:
            max_independent_rel_error = max(max_independent_rel_error, independent_rel_error)
            positive_reference_count += 1
            minimum_positive_reference = (independent if minimum_positive_reference is None
                                          else min(minimum_positive_reference, independent))
        else:
            zero_reference_count += 1
        max_forced_error = max(max_forced_error, forced_error)
        max_independent_abs_error = max(max_independent_abs_error, independent_abs_error)
        samples.append(dict(sample_id=sample_id, eigenphase=result["eigenphase"], output=result["output"],
                             rejection_attempts=result["rejection_attempts"],
                             joint_probability=result["joint_latent_output_probability"],
                             forced_joint_probability=forced["joint_latent_output_probability"],
                             independent_joint_probability=independent,
                             forced_error=float(forced_error), independent_abs_error=float(independent_abs_error),
                             independent_rel_error=None if independent_rel_error is None else float(independent_rel_error),
                             independent_reference_positive=independent_positive,
                             oracle_column_calls=after["column_calls"] - before["column_calls"],
                             oracle_phase_rows=after["phase_rows"] - before["phase_rows"],
                             early_amplitudes=result["early_amplitudes"],
                             phase_row_payload_bytes=after["phase_row_payload_bytes"]))
    finite = all(np.isfinite([s["joint_probability"], s["forced_joint_probability"],
                              s["independent_joint_probability"]]).all() for s in samples)
    forced_ok = max_forced_error < 2e-12
    independent_ok = max_independent_abs_error < 2e-10
    exp.check("P2", forced_ok and independent_ok,
              f"{name}: max forced error={max_forced_error:.2e}, independent abs={max_independent_abs_error:.2e}, "
              f"rel={max_independent_rel_error:.2e}")
    exp.check("P3", finite and all(s["early_amplitudes"] == length for s in samples),
              f"{name}: finite samples, early payload={length}, setup columns={setup_stats['column_calls']}")

    exact = None
    full_error = None
    clean_error = None
    uniform_error = None
    no_feedback_error = None
    weights = None
    if full_check:
        exact = spectral_outputs(period, width, v_orbit)
        forced_marginal = np.zeros(1 << width)
        conditional = np.zeros((period, 1 << width))
        weights = np.array([prefix.phase_probability(k) for k in range(period)])
        for k in range(period):
            for y in range(1 << width):
                joint = prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
                forced_marginal[y] += joint
                conditional[k, y] = joint / weights[k]
        full_error = float(np.max(np.abs(forced_marginal - exact)))
        exp.check("P1", full_error < 2e-10,
                  f"{name}: forced marginal vs spectral max error={full_error:.2e}")
        if reference_width4:
            clean_error = float(np.max(np.abs(clean_fourier_reference(period, width, v_orbit) - exact)))
            exp.check("P1", clean_error < 2e-10,
                      f"{name}: independent clean Fourier vs spectral max error={clean_error:.2e}")
        if width == 4:
            uniform = np.mean(conditional, axis=0)
            uniform_error = float(np.max(np.abs(uniform - exact)))
            no_feedback = independent_marginal(period, width, feedback=False)
            no_feedback_error = float(np.max(np.abs(no_feedback - exact)))
            exp.fail_check("C1", np.max(weights) - np.min(weights) > 1e-10 and uniform_error > 1e-8,
                           f"weight range={weights.min():.6g}..{weights.max():.6g}, uniform error={uniform_error:.6g}")
            exp.fail_check("C2", no_feedback_error > 1e-8,
                           f"feedback-erased marginal error={no_feedback_error:.6g}")

    return dict(series=name, period=period, width=width, split=SPLIT, full_enumeration=full_check,
                setup_seconds=setup_seconds, setup_stats=setup_stats, samples=samples,
                max_forced_error=float(max_forced_error), max_independent_abs_error=float(max_independent_abs_error),
                max_independent_rel_error=float(max_independent_rel_error),
                positive_reference_count=positive_reference_count,
                zero_reference_count=zero_reference_count,
                minimum_positive_reference=None if minimum_positive_reference is None
                else float(minimum_positive_reference),
                full_marginal_error=full_error, clean_reference_error=clean_error,
                uniform_k_error=uniform_error, no_feedback_error=no_feedback_error,
                weight_range=None if weights is None else [float(weights.min()), float(weights.max())],
                oracle_contract="known orbit indices; no physical labels or discrete logs",
                note="abstract cyclic-orbit oracle; no compiled arithmetic or physical-qubit claim")


def main() -> None:
    exp = Experiment("prefix_scaling", doc=__doc__)
    exp.predict("P1", "small r=6 forced marginals match spectral and clean Fourier references")
    exp.predict("P2", "fixed-seed samples match forced joints and independent split formula")
    exp.predict("P3", "sample accounting remains finite and stores only the early payload")
    exp.predict("P4", "orbit-index validation rejects an out-of-range computational label")
    exp.must_fail("C1", "uniform final-k control fails on the nonuniform small row")
    exp.must_fail("C2", "feedback-erased early control changes the small-row marginal")

    rows = []
    started = time.perf_counter()
    # Invalid physical-label-as-orbit-index oracle: label 6 is outside r=6's
    # index range [0,5], and must be rejected before any sampling.
    try:
        SparseOrbitPrefix(6, SPLIT, lambda l: (np.array([6]), np.array([1.+0j])), 1)
        rejected = False
    except ValueError:
        rejected = True
    exp.check("P4", rejected, "out-of-range computational label rejected by orbit-index validator")

    for width in WIDTHS:
        rows.append(one_row(exp, f"r6_width{width}", 6, width,
                            full_check=width <= 8, reference_width4=width == 4))
    for period in ABSTRACT_PERIODS:
        rows.append(one_row(exp, f"abstract_r{period}_width{ABSTRACT_WIDTH}", period,
                            ABSTRACT_WIDTH, full_check=False))

    exp.finish(report_path="out/prefix_scaling.json", rows=rows,
               metadata=dict(split=SPLIT, length=1 << SPLIT, theta_over_pi=.5,
                             small_period=6, small_widths=list(WIDTHS),
                             abstract_periods=list(ABSTRACT_PERIODS), abstract_width=ABSTRACT_WIDTH,
                             samples_per_row=SAMPLES, dense_enumeration_cap="only r=6,width<=8",
                             independent_wide_reference="integer-reduced finite geometric H-tail plus direct length-4 DFT",
                             oracle="known cyclic orbit indices; setup/order discovery/discrete logs external",
                             no_physical_qubit_or_speedup_claim=True,
                             max_runtime_seconds=10, elapsed_seconds=time.perf_counter() - started,
                             numpy=np.__version__))


if __name__ == "__main__":
    main()
