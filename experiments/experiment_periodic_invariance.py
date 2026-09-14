"""Repeated orbit-block mixers: 2-adic invisibility versus odd-period visibility.

TODO 20 candidate, derived before measurement.  On an indexed orbit of period
r, let V repeat the same b-by-b unitary W on every consecutive block.  When
b|r it commutes with the orbit translation U^b.  If b|L=2^s, every arithmetic
power after the insertion is a multiple of b, so V can be moved to the end of
the ascending-power circuit and cannot change exponent output probabilities.
For odd b>1, b never divides a power of two; the same move is not licensed.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  For b=2,r=6 and b=4,r=12, the output equals the ideal distribution at
      and after s=log2(b), while the row immediately below threshold is
      non-ideal.  The same fixed seeded W is used across each threshold.
  P2  For fixed seeded b=3 W, r=6 and r=9 at s=3 retain a nonzero TV
      distance from ideal; this is the full-orbit effect, not a finite defect.
  P3  The repeated-block row decomposition h_k(l)=exp(2*pi*i*k*l/r)f[l mod
      b](k), and its b-periodic Fourier expansion, agrees with the streamed
      prefix reference without amplitude truncation.
  P4  Float64 spectral output agrees with an independent long-double spectral
      evaluation and with the prefix reference on representative rows.

  C1  Moving an odd-b mixer to the end (using ideal output) must fail at a
      visible row.
  C2  A b-by-b construction with b not dividing r must fail the U^b
      commutation check; its incomplete-block row is not admitted as a
      periodic promise.
  C3  The nontrivial result must disappear for identity W.

This is an indexed cyclic-orbit experiment.  It is not a physical local-gate
claim, order discovery, or a second propagator.  References are the existing
original-time-order spectral effects and SparseOrbitPrefix sampler.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_periodic_invariance
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
from lab.semiclassical import order_finding_probability
from lab.spectral import single_defect_effects


MAX_DENSE_BYTES = 16 * 1024 * 1024
WIDTH = 8


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def alloc(shape, dtype, label: str) -> np.ndarray:
    guard(shape, dtype, label)
    return np.empty(shape, dtype=dtype)


def report_path(prefix: str = "periodic_invariance") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def seeded_unitary(b: int, seed: int) -> np.ndarray:
    """One fixed complex unitary, reused across every row for this b."""
    rng = np.random.default_rng(seed)
    guard((b, b), np.complex128, "seeded W")
    z = rng.normal(size=(b, b)) + 1j * rng.normal(size=(b, b))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    q *= np.where(np.abs(diagonal) > 0, np.conj(diagonal) / np.abs(diagonal), 1.)
    return q


def periodic_block_unitary(period: int, b: int, w: np.ndarray) -> np.ndarray:
    if not 1 <= period <= 18 or not 1 <= b <= 18 or period % b:
        raise ValueError("periodic construction requires b dividing r")
    w = np.asarray(w, dtype=complex)
    if w.shape != (b, b) or not np.all(np.isfinite(w)):
        raise ValueError("W shape/finite check failed")
    if np.max(np.abs(w.conj().T @ w - np.eye(b))) > 2e-12:
        raise ValueError("W must be unitary")
    blocks = period // b
    guard((period, period), np.complex128, "periodic orbit unitary")
    return np.kron(np.eye(blocks, dtype=complex), w)


def incomplete_block_unitary(period: int, b: int, w: np.ndarray) -> np.ndarray:
    """A deliberately invalid floor-block construction for must-fail control."""
    if period % b == 0:
        raise ValueError("control requires an incomplete final block")
    guard((period, period), np.complex128, "incomplete orbit unitary")
    result = np.eye(period, dtype=complex)
    for start in range(0, period - b + 1, b):
        result[start:start + b, start:start + b] = w
    return result


def translation(period: int, step: int) -> np.ndarray:
    guard((period, period), np.complex128, "orbit translation")
    result = np.zeros((period, period), dtype=complex)
    for j in range(period):
        result[(j + step) % period, j] = 1.
    return result


def commutator_error(v: np.ndarray, period: int, step: int) -> float:
    u = translation(period, step)
    return float(np.max(np.abs(v @ u - u @ v)))


def orbit_fourier(period: int, dtype=np.float64) -> np.ndarray:
    guard((period, period), np.result_type(dtype, 1j), "orbit Fourier")
    real = np.dtype(dtype).type
    j = np.arange(period, dtype=dtype)
    pi = np.arccos(real(-1))
    return np.exp(-2j * pi * j[:, None] * j[None, :] / period) / np.sqrt(real(period))


def spectral_output(v_orbit: np.ndarray, period: int, width: int, split: int,
                    *, dtype=np.float64) -> np.ndarray:
    """Existing spectral route, with a selectable precision for the audit."""
    qsize = 1 << width
    guard((qsize, period, period), np.result_type(dtype, 1j), "spectral effects")
    f = orbit_fourier(period, dtype=dtype)
    ve = f.conj().T @ np.asarray(v_orbit, dtype=np.result_type(dtype, 1j)) @ f
    effects = single_defect_effects(period, width, split, ve, dtype=dtype,
                                    max_entries=1_000_000)
    return effects.sum(axis=(1, 2)).real / period


def prefix_output(v_orbit: np.ndarray, period: int, width: int, split: int):
    """Existing streamed sparse-column prefix reference, no second propagator."""
    length = 1 << split

    def column(l):
        values = v_orbit[:, l % period]
        keep = np.flatnonzero(values != 0)
        return keep.astype(np.int64), values[keep]

    prefix = SparseOrbitPrefix(period, split, column, max_terms=min(period, 64))
    guard((1 << width,), np.float64, "prefix output")
    result = np.zeros(1 << width, dtype=float)
    for k in range(period):
        for y in range(1 << width):
            result[y] += prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
    return result, prefix.stats()


def ideal_output(period: int, width: int) -> np.ndarray:
    return np.array([order_finding_probability(y, width, period)
                     for y in range(1 << width)], dtype=float)


def row_harmonic_error(v_orbit: np.ndarray, period: int, b: int, split: int,
                       k: int) -> tuple[float, float]:
    """Check the derived period-b multiplier and its b Fourier modes."""
    length = 1 << split
    # Prefix's column convention is the orbit basis: h_k(l)=<phi_k|V|l>.
    def h_at(l: int) -> complex:
        return sum(v_orbit[q, l % period] *
                   np.exp(2j * np.pi * k * q / period)
                   for q in range(period))

    h = np.array([h_at(l) for l in range(length)], complex)
    base = np.exp(2j * np.pi * k * np.arange(length) / period)
    base_b = np.exp(2j * np.pi * k * np.arange(b) / period)
    f = np.array([h_at(p) / base_b[p] for p in range(b)], complex)
    c = np.array([np.sum(f * np.exp(-2j * np.pi * u * np.arange(b) / b)) / b
                  for u in range(b)], complex)
    expanded = base * np.array([sum(c[u] * np.exp(2j * np.pi * u * l / b)
                                     for u in range(b)) for l in range(length)])
    return float(np.max(np.abs(h - base * f[np.arange(length) % b]))), float(np.max(np.abs(h - expanded)))


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(np.asarray(a) - np.asarray(b))))


def complex_payload(matrix: np.ndarray) -> dict:
    """JSON-safe report payload without discarding any amplitudes."""
    return dict(real=np.asarray(matrix).real.tolist(), imag=np.asarray(matrix).imag.tolist())


def evaluate_case(name: str, period: int, width: int, split: int, block: np.ndarray,
                  b: int, exp: Experiment, *, expect_ideal: bool,
                  fixed_precision_audit: bool = False) -> dict:
    v = periodic_block_unitary(period, b, block)
    ideal = ideal_output(period, width)
    spectral = spectral_output(v, period, width, split)
    prefix, stats = prefix_output(v, period, width, split)
    spectral_prefix_error = float(np.max(np.abs(spectral - prefix)))
    total_error = max(abs(float(spectral.sum()) - 1.), abs(float(prefix.sum()) - 1.))
    tv_ideal = tv(spectral, ideal)
    comm = commutator_error(v, period, b)
    harmonic = [row_harmonic_error(v, period, b, split, k) for k in range(period)]
    harmonic_error = float(max(max(x) for x in harmonic))
    exp.check("P3", harmonic_error < 3e-11 and comm < 3e-11,
              f"{name}: comm(U^{b})={comm:.2e}, harmonic={harmonic_error:.2e}")
    exp.check("P4", spectral_prefix_error < 3e-10 and total_error < 3e-10,
              f"{name}: spectral/prefix={spectral_prefix_error:.2e}, total={total_error:.2e}")
    if expect_ideal:
        exp.check("P1", tv_ideal < 3e-9,
                  f"{name}: TV(defect,ideal)={tv_ideal:.3e} (invisible expected)")
    else:
        exp.check("P2", tv_ideal > 1e-5,
                  f"{name}: TV(defect,ideal)={tv_ideal:.6g} (visible expected)")
    precision_error = None
    if fixed_precision_audit:
        extended = spectral_output(v, period, width, split, dtype=np.longdouble)
        precision_error = float(np.max(np.abs(spectral - extended)))
        exp.check("P4", precision_error < 2e-9,
                  f"{name}: float64/longdouble={precision_error:.2e}")
    return dict(name=name, period=period, width=width, split=split, L=1 << split,
                b=b, commutator_error=comm, tv_from_ideal=tv_ideal,
                spectral_prefix_error=spectral_prefix_error, total_error=total_error,
                harmonic_error=harmonic_error, precision_error=precision_error,
                spectral=spectral.tolist(), prefix=prefix.tolist(), ideal=ideal.tolist(),
                prefix_stats=stats, W=complex_payload(block))


def main() -> None:
    exp = Experiment("periodic_invariance", doc=__doc__)
    exp.predict("P1", "power-of-two b becomes output-invisible at and beyond s=log2(b)")
    exp.predict("P2", "fixed odd-b full-orbit mixers retain nonzero TV from ideal")
    exp.predict("P3", "period-b row/harmonic identity and U^b commutation hold when b|r")
    exp.predict("P4", "spectral and prefix references agree, including precision audit")
    exp.must_fail("C1", "moving an odd-b mixer to the end must fail at visible s")
    exp.must_fail("C2", "incomplete b-block with b not dividing r must fail commutation")
    exp.must_fail("C3", "identity W removes the nontrivial periodic effect")

    start = time.perf_counter()
    rows = []
    w2 = seeded_unitary(2, 20260220)
    w4 = seeded_unitary(4, 20260240)
    w3 = seeded_unitary(3, 20260230)

    # One-parameter threshold rows: each W is fixed while only s changes.
    for b, r, threshold, w in ((2, 6, 1, w2), (4, 12, 2, w4)):
        for split in (threshold - 1, threshold, threshold + 1):
            rows.append(evaluate_case(f"power2_b{b}_r{r}_s{split}", r, WIDTH, split,
                                      w, b, exp, expect_ideal=(split >= threshold),
                                      fixed_precision_audit=(split == threshold)))

    # Odd b: same W for both periods and a fixed s.  These are abstract indexed
    # orbits; the order series below is deliberately descriptive, not monotonic.
    for r in (6, 9):
        rows.append(evaluate_case(f"odd_b3_r{r}_s3", r, WIDTH, 3, w3, 3, exp,
                                  expect_ideal=False, fixed_precision_audit=(r == 6)))
    for r in range(3, 19, 3):
        rows.append(evaluate_case(f"odd_b3_order_series_r{r}_s3", r, WIDTH, 3,
                                  w3, 3, exp, expect_ideal=False))

    # Controls are measured only after the positive rows and use nontrivial W.
    odd_v = periodic_block_unitary(6, 3, w3)
    odd_spectral = spectral_output(odd_v, 6, WIDTH, 3)
    odd_ideal = ideal_output(6, WIDTH)
    moved_end_error = tv(odd_spectral, odd_ideal)
    exp.fail_check("C1", moved_end_error > 1e-5,
                   f"odd b=3 moved-to-end/ideal TV={moved_end_error:.6g}")

    incomplete = incomplete_block_unitary(7, 3, w3)
    incomplete_comm = commutator_error(incomplete, 7, 3)
    exp.fail_check("C2", incomplete_comm > 1e-6,
                   f"incomplete r=7,b=3 commutator max={incomplete_comm:.6g}")

    identity = np.eye(3, dtype=complex)
    identity_v = periodic_block_unitary(6, 3, identity)
    identity_spectral = spectral_output(identity_v, 6, WIDTH, 3)
    identity_tv = tv(identity_spectral, ideal_output(6, WIDTH))
    exp.fail_check("C3", identity_tv < 3e-10,
                   f"identity W TV from ideal={identity_tv:.3e} (effect vanished)")

    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(width=WIDTH, threshold_rows="s=log2(b)-1,log2(b),log2(b)+1",
                             fixed_seeds={"b2": 20260220, "b3": 20260230, "b4": 20260240},
                             max_dense_bytes=MAX_DENSE_BYTES,
                             moved_end_control_tv=moved_end_error,
                             incomplete_commutator=incomplete_comm,
                             identity_control_tv=identity_tv,
                             order_series="r=3,6,...,18; b=3,s=3; no monotonicity asserted",
                             references="existing lab.spectral and lab.prefix; no second propagator",
                             elapsed_seconds=elapsed, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("periodic_invariance_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
