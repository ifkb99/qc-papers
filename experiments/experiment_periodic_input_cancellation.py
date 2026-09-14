"""Challenge the L-divides-r input-subspace cancellation lemma.

DERIVED BEFORE MEASUREMENT.  Use the ascending-power order with one work
unitary V after s early controls, whose range has L=2**s labels.  The candidate
lemma is:

    if L divides the known orbit period r, and V preserves the span of
    orbit labels 0,...,L-1 unitarily, then tracing work leaves the ideal
    exponent output unchanged.  No promise is made about V on its orthogonal
    complement.

The reason to test this separately is that a repeated b=3 block can look
invisible at one split for accidental divisibility, while being visible at a
nearby split.  This experiment uses only the existing lab.spectral effects
reference: it constructs finite orbit matrices and compares their output
marginals with the identity defect.  It does not add a circuit propagator.

P1: seeded genuinely complex block-diagonal V, with arbitrary seeded unitary
    action on the unreachable complement, is exactly ideal whenever L|r.
P2: the result survives L=1 and L=r endpoint rows, with normalized output
    probabilities and no dense allocation beyond r<=18,t<=7.
P3: the initially failed leakage control (0<->L at r=2L) is ideal because
    its rotation commutes with U^L; retain this post-audit prediction explicitly.

C1: L not dividing r with a nontrivial first-span block must be detectable.
C2: mixing a reached label with an unreached label must be detectable even
    when L|r.

The controls are deliberately narrow: a particular non-dividing or leakage
matrix could still be accidentally invisible.  Such a control is recorded,
not promoted to a universal converse.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_periodic_input_cancellation
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.spectral import single_defect_effects


MAX_DENSE_BYTES = 16 * 1024 * 1024
MAX_R = 18
MAX_T = 7


def guard(shape, dtype, label: str) -> int:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")
    return payload


def report_path(prefix: str = "periodic_input_cancellation") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def seeded_unitary(size: int, seed: int) -> np.ndarray:
    if not 1 <= size <= MAX_R:
        raise ValueError("seeded unitary dimension outside experiment cap")
    guard((size, size), np.complex128, "seeded unitary")
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    q *= np.where(np.abs(diagonal) > 0, np.conj(diagonal) / np.abs(diagonal), 1.)
    return q


def embed_preserving(period: int, reached: int, first: np.ndarray,
                     complement: np.ndarray) -> np.ndarray:
    if period < 1 or period > MAX_R or not 1 <= reached <= period:
        raise ValueError("invalid period/reached span")
    if first.shape != (reached, reached) or complement.shape != (period-reached, period-reached):
        raise ValueError("block shapes do not match")
    guard((period, period), np.complex128, "embedded defect")
    result = np.zeros((period, period), dtype=complex)
    result[:reached, :reached] = first
    if period > reached:
        result[reached:, reached:] = complement
    return result


def leakage_rotation(period: int, reached: int, theta: float = np.pi / 3,
                     *, offset: int = 1) -> np.ndarray:
    """A unitary mixing one reached label with one unreached label."""
    # Choosing reached+1 rather than reached itself avoids the special
    # translation-by-L pairing when the first unreached label is exactly one
    # block away.  The caller still supplies the reached-span size.
    j = reached + offset
    if not 1 <= period <= MAX_R or not 1 <= reached <= j < period:
        raise ValueError("leakage control needs an unreached label")
    guard((period, period), np.complex128, "leakage rotation")
    result = np.eye(period, dtype=complex)
    c, s = np.cos(theta / 2), -1j * np.sin(theta / 2)
    result[0, 0] = result[j, j] = c
    result[0, j] = result[j, 0] = s
    return result


def orbit_fourier(period: int) -> np.ndarray:
    if not 1 <= period <= MAX_R:
        raise ValueError("orbit Fourier outside experiment cap")
    guard((period, period), np.complex128, "orbit Fourier")
    j = np.arange(period, dtype=float)
    return np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)


def marginal(defect_orbit: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    """Existing spectral-effects route for one intervening orbit defect."""
    if not 1 <= period <= MAX_R or not 0 <= width <= MAX_T or not 0 <= split <= width:
        raise ValueError("experiment cap or split violation")
    if not isinstance(defect_orbit, np.ndarray) or defect_orbit.shape != (period, period):
        raise ValueError("defect must have the capped orbit shape")
    guard((1 << width, period, period), np.complex128, "spectral effects")
    guard((1 << width,), np.float64, "output marginal")
    fourier = orbit_fourier(period)
    defect_eigen = fourier.conj().T @ defect_orbit @ fourier
    effects = single_defect_effects(period, width, split, defect_eigen,
                                    max_entries=1_000_000)
    result = effects.sum(axis=(1, 2)).real / period
    return result


def tv_distance(left: np.ndarray, right: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(left - right)))


def check_case(name: str, period: int, split: int, defect: np.ndarray,
               ideal: np.ndarray, exp: Experiment, *, expected_ideal: bool) -> dict:
    width = MAX_T
    actual = marginal(defect, period, width, split)
    max_error = float(np.max(np.abs(actual - ideal)))
    tv = tv_distance(actual, ideal)
    total_error = float(abs(actual.sum() - 1.))
    unitary_error = float(np.max(np.abs(defect.conj().T @ defect - np.eye(period))))
    if expected_ideal:
        exp.check("P1", max_error < 3e-11,
                  f"{name}: max ideal-output error={max_error:.3e}, TV={tv:.3e}")
        exp.check("P2", total_error < 3e-11 and unitary_error < 3e-12,
                  f"{name}: total={total_error:.3e}, unitarity={unitary_error:.3e}")
    return dict(name=name, period=period, width=width, split=split, L=1 << split,
                max_error=max_error, tv=tv, total_error=total_error,
                unitary_error=unitary_error, expected_ideal=expected_ideal)


def main() -> None:
    exp = Experiment("periodic_input_cancellation", doc=__doc__)
    exp.predict("P1", "L|r and preservation of the reached span imply ideal exponent output")
    exp.predict("P2", "endpoint rows and arbitrary complement action remain normalized and ideal")
    exp.predict("P3", "the old 0-to-L leakage control commutes with U^L and remains ideal")
    exp.must_fail("C1", "a non-dividing period control is visibly non-ideal")
    exp.must_fail("C2", "mixing a reached and unreached label is visibly non-ideal")

    rows = []
    # Block diagonal W on the reached span and an independent genuinely
    # complex Z on its complement.  Only the first block is ever reached
    # before V; Z tests the lemma's explicit no-assumption-outside-S scope.
    for period, split, seed_w, seed_z in ((6, 1, 6101, 6102),
                                           (8, 2, 6201, 6202),
                                           (8, 3, 6301, 6302),
                                           (12, 2, 6401, 6402),
                                           (8, 0, 6501, 6502),
                                           (8, 3, 6601, 6602)):
        reached = 1 << split
        first = seeded_unitary(reached, seed_w)
        complement = seeded_unitary(period-reached, seed_z) if period > reached else np.empty((0, 0), complex)
        defect = embed_preserving(period, reached, first, complement)
        ideal = marginal(np.eye(period, dtype=complex), period, MAX_T, split)
        rows.append(check_case(f"preserving_r{period}_s{split}", period, split,
                               defect, ideal, exp, expected_ideal=True))

    # Endpoint L=r with a full genuinely complex W is a useful stress row:
    # there is no complement at all, yet the statement still predicts ideal.
    period, split = 8, 3
    full = seeded_unitary(period, 6708)
    ideal = marginal(np.eye(period, dtype=complex), period, MAX_T, split)
    rows.append(check_case("preserving_full_r8_s3", period, split, full,
                           ideal, exp, expected_ideal=True))

    # C1: L does not divide r.  The span-preserving construction is still a
    # valid unitary, but the shifted copies of the reached span overlap in a
    # different residue pattern.  This particular seeded row must deviate.
    period, split = 6, 2  # L=4 does not divide 6
    reached = 1 << split
    nondividing = embed_preserving(period, reached, seeded_unitary(reached, 6804),
                                   seeded_unitary(period-reached, 6805))
    ideal = marginal(np.eye(period, dtype=complex), period, MAX_T, split)
    actual = marginal(nondividing, period, MAX_T, split)
    nondividing_error = float(np.max(np.abs(actual - ideal)))
    nondividing_tv = tv_distance(actual, ideal)
    exp.fail_check("C1", nondividing_error > 1e-7,
                   f"r={period}, L={1 << split}: max error={nondividing_error:.9g}, "
                   f"TV={nondividing_tv:.9g}")
    rows.append(dict(name="nondividing_r6_s2", period=period, split=split,
                     L=1 << split, max_error=nondividing_error,
                     tv=nondividing_tv, expected_ideal=False))

    # C2: L|r, but V rotates a reached label with an unreached one.  This
    # violates exactly the span-preservation hypothesis while retaining the
    # arithmetic divisibility condition.
    period, split = 8, 2  # L=4 divides 8; labels 0 and 5 are mixed
    commuting_leak = leakage_rotation(period, 1 << split, offset=0)
    late_shift = np.roll(np.eye(period, dtype=complex), 1 << split, axis=0)
    old_control = marginal(commuting_leak, period, MAX_T, split)
    ideal = marginal(np.eye(period, dtype=complex), period, MAX_T, split)
    old_control_error = float(np.max(np.abs(old_control-ideal)))
    commutator_error = float(np.max(np.abs(commuting_leak @ late_shift-late_shift @ commuting_leak)))
    exp.check("P3", old_control_error < 3e-11 and commutator_error < 3e-12,
              f"old 0-to-L leakage: ideal error={old_control_error:.3e}, commutator={commutator_error:.3e}")
    rows.append(dict(name="commuting_leakage_r8_s2", period=period, split=split,
                     L=1 << split, max_error=old_control_error,
                     commutator_error=commutator_error, expected_ideal=True,
                     violates_span_preservation=True))
    leaking = leakage_rotation(period, 1 << split)
    ideal = marginal(np.eye(period, dtype=complex), period, MAX_T, split)
    actual = marginal(leaking, period, MAX_T, split)
    leakage_error = float(np.max(np.abs(actual - ideal)))
    leakage_tv = tv_distance(actual, ideal)
    exp.fail_check("C2", leakage_error > 1e-7,
                   f"r={period}, L={1 << split}: max error={leakage_error:.9g}, "
                   f"TV={leakage_tv:.9g}")
    rows.append(dict(name="leakage_r8_s2", period=period, split=split,
                     L=1 << split, max_error=leakage_error,
                     tv=leakage_tv, expected_ideal=False))

    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(
        width=MAX_T, max_r=MAX_R, dense_budget_bytes=MAX_DENSE_BYTES,
        reference="lab.spectral.single_defect_effects", generic_propagator=False,
        lemma_scope="one intervening work unitary, known orbit, L divides r, V preserves first L span",
        nondividing_error=nondividing_error, nondividing_tv=nondividing_tv,
        leakage_error=leakage_error, leakage_tv=leakage_tv))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("periodic_input_cancellation_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
