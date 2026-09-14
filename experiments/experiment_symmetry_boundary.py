"""Fixed-output low-rank boundary formula for one localized symmetry break.

This is a formula test, not a sampler.  Keep the C56 periodic background
W01@s1, W12@s3, W01@s4 on the r=6, b=3 orbit, and insert one localized kick
K=I+E Delta E^dagger after split=2, where E selects orbit labels 0 and 1.
For a complete fixed output y, Q=2^t and M=r/b, the coarse-sector factors
are evaluated independently:

  A_i(alpha,y) = W_(i+1) [I + exp(-2*pi*i*y*2^i/Q) T_alpha^(2^i)] / 2
  P_alpha = A_(split-1)...A_0
  S_alpha = A_(t-1)...A_split
  c(y) = sum_alpha H_alpha^dagger P_alpha |0>/sqrt(M)
  a_gamma = S_gamma [P_gamma |0>/sqrt(M) + H_gamma Delta c(y)]
  p(y) = sum_gamma ||a_gamma||^2.

The formula is derived from the pre-kick sector amplitudes and the rank-d
update K-I=E Delta E^dagger.  It streams the M sectors and keeps no M-sized
array.  Its cost still sums over M=r/b and it is not an exact sampler.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The streamed boundary formula agrees for every y with both the existing
      original-time-order statevec output and full-r sequential_path reference
      at theta=0, pi/4, pi/2, pi.
  P2  The formula and each reference normalize to one without amplitude
      truncation, including cancellation diagnostics.
  P3  At theta=0 the correction vanishes; replacing K by a repeated block is
      not required for this finite-support formula.

  C1  Deleting baseline/correction cross terms must fail on a nonzero angle.
  C2  Replacing cross-sector c(y) by the same-sector contribution must fail.

All arrays are capped before allocation (t<=7,r<=18,16 MiB).  This records a
finite fixed-output identity only; low rank does not make the cost independent
of r and no compressed sampler or generic propagator is introduced.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_symmetry_boundary
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
from lab.fourier_sampling import unit_phase
from experiments import experiment_symmetry_kick as kick


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, WIDTH, PERIOD, BLOCK, SPLIT = 7, 3, 5, 6, 3, 2
M = PERIOD // BLOCK
THETAS = np.pi * np.array([0., .25, .5, 1.])


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def alloc(shape, dtype, label: str) -> np.ndarray:
    guard(shape, dtype, label)
    return np.empty(shape, dtype=dtype)


def report_path(prefix: str = "symmetry_boundary") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def sector_shift(power: int, alpha: int) -> np.ndarray:
    """T_alpha^power in the b-dimensional coarse sector."""
    result = alloc((BLOCK, BLOCK), np.complex128, "sector shift")
    result.fill(0.)
    for p in range(BLOCK):
        wraps, target = divmod(p + power, BLOCK)
        result[target, p] = unit_phase(alpha * wraps, M)
    return result


def sector_embedding(alpha: int) -> np.ndarray:
    """H_alpha=F_alpha^dagger E for E={orbit labels 0,1}."""
    result = alloc((BLOCK, 2), np.complex128, "boundary embedding")
    result.fill(0.)
    for j in (0, 1):
        m, p = divmod(j, BLOCK)
        result[p, j] = unit_phase(alpha * m, M) / np.sqrt(M)
    return result


def background_sector_blocks() -> dict[int, np.ndarray]:
    return kick.background_sector_blocks()


def boundary_formula(theta: float, *, interference: bool = True,
                     cross_sector: bool = True) -> tuple[np.ndarray, dict]:
    """Stream all sectors for every output, retaining only small factors."""
    if not 1 <= PERIOD <= 18 or not 0 <= SPLIT <= WIDTH <= 7 or PERIOD != BLOCK*M:
        raise ValueError("bounded formula fixture dimensions")
    qsize = 1 << WIDTH
    guard((qsize,), np.float64, "validation-only output enumeration")
    identity = np.eye(BLOCK, dtype=complex)
    blocks = background_sector_blocks()
    delta = kick.localized_kick_block(theta)[:2, :2] - np.eye(2, dtype=complex)
    v = np.array([1., 0., 0.], dtype=complex)
    probabilities = np.zeros(qsize, dtype=float)
    max_norm = max_cancellation = 0.
    max_correction = max_baseline = 0.
    for y in range(qsize):
        # First pass: c(y)=sum_alpha H_alpha^dagger P_alpha v/sqrt(M).
        c = np.zeros(2, dtype=complex)
        for alpha in range(M):
            p_factor = identity.copy()
            for i in range(SPLIT):
                w = blocks.get(i + 1, identity)
                phase = unit_phase(-y * (1 << i), qsize)
                a = w @ (identity + phase * sector_shift(1 << i, alpha)) / 2
                p_factor = a @ p_factor
            h = sector_embedding(alpha)
            baseline = p_factor @ v / np.sqrt(M)
            c += h.conj().T @ baseline
        if not cross_sector:
            # Deliberately wrong control: c_gamma is replaced later by the
            # same sector's H_gamma^dagger baseline.
            c = None
        # Second pass RECOMPUTES each early factor. Keeping the first-pass
        # sector factors in a list would take O(M) storage, not streaming.
        for alpha in range(M):
            p_factor = identity.copy()
            for i in range(SPLIT):
                w = blocks.get(i + 1, identity)
                phase = unit_phase(-y * (1 << i), qsize)
                a = w @ (identity + phase * sector_shift(1 << i, alpha)) / 2
                p_factor = a @ p_factor
            h = sector_embedding(alpha)
            baseline = p_factor @ v / np.sqrt(M)
            s_factor = identity.copy()
            for i in range(SPLIT, WIDTH):
                w = blocks.get(i + 1, identity)
                phase = unit_phase(-y * (1 << i), qsize)
                a = w @ (identity + phase * sector_shift(1 << i, alpha)) / 2
                s_factor = a @ s_factor
            if cross_sector:
                c_used = c
            else:
                c_used = h.conj().T @ baseline
            correction = h @ delta @ c_used
            max_baseline = max(max_baseline, float(np.linalg.norm(s_factor @ baseline)))
            max_correction = max(max_correction, float(np.linalg.norm(s_factor @ correction)))
            if interference:
                amplitude = s_factor @ (baseline + correction)
            else:
                # Return the incoherent baseline-plus-correction control law.
                amplitude = s_factor @ baseline
                probabilities[y] += float(np.vdot(amplitude, amplitude).real)
                probabilities[y] += float(np.vdot(s_factor @ correction,
                                                   s_factor @ correction).real)
                continue
            sector_probability = float(np.vdot(amplitude, amplitude).real)
            probabilities[y] += sector_probability
            max_norm = max(max_norm, sector_probability)
        max_cancellation = max(max_cancellation,
                               float(np.linalg.norm(delta @ (c if c is not None
                                                             else np.zeros(2)))))
    return probabilities, dict(max_probability=float(np.max(probabilities)),
                               min_probability=float(np.min(probabilities)),
                               total=float(probabilities.sum()),
                               normalization_error=abs(float(probabilities.sum()) - 1),
                               max_sector_probability=max_norm,
                               max_boundary_correction_norm=max_cancellation,
                               max_baseline_amplitude=max_baseline,
                               max_correction_amplitude=max_correction,
                               delta_norm=float(np.linalg.norm(delta)),
                               sectors=M, support_dimension=2,
                               streamed_sector_storage="one sector at a time")


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(np.asarray(a) - np.asarray(b))) / 2)


def main() -> None:
    exp = Experiment("symmetry_boundary", doc=__doc__)
    exp.predict("P1", "boundary formula matches physical and full-r references for every output")
    exp.predict("P2", "formula/reference rows normalize and retain cancellation diagnostics")
    exp.predict("P3", "zero kick has vanishing correction while remaining a valid fixed-output row")
    exp.must_fail("C1", "delete baseline/correction interference")
    exp.must_fail("C2", "replace cross-sector c(y) with same-sector contribution")
    start = time.perf_counter()
    rows = []
    max_deleted_tv = max_same_sector_tv = 0.
    for theta in THETAS:
        theta = float(theta)
        formula, info = boundary_formula(theta)
        physical, physical_info = kick.physical_output(theta)
        sequential = kick.sequential_marginal(theta)
        formula_physical_error = float(np.max(np.abs(formula - physical)))
        formula_sequential_error = float(np.max(np.abs(formula - sequential)))
        reference_norm_error = max(abs(float(physical.sum()) - 1),
                                   abs(float(sequential.sum()) - 1))
        exp.check("P1", formula_physical_error < 3e-10
                  and formula_sequential_error < 3e-10,
                  f"theta/pi={theta/np.pi:g}: formula/physical={formula_physical_error:.2e}, "
                  f"formula/sequential={formula_sequential_error:.2e}")
        exp.check("P2", info["normalization_error"] < 3e-10
                  and reference_norm_error < 3e-10,
                  f"theta/pi={theta/np.pi:g}: formula norm={info['normalization_error']:.2e}, "
                  f"reference norm={reference_norm_error:.2e}")
        if theta == 0:
            exp.check("P3", info["max_correction_amplitude"] < 3e-14,
                      f"theta=0 correction={info['max_correction_amplitude']:.2e}")
        deleted, deleted_info = boundary_formula(theta, interference=False)
        same_sector, same_info = boundary_formula(theta, cross_sector=False)
        # These wrong weights are not generally normalized laws. Keep their
        # raw half-L1 discrepancy, then compare normalized laws for real TV.
        deleted_raw_half_l1 = tv(deleted, physical)
        same_raw_half_l1 = tv(same_sector, physical)
        deleted_tv = tv(deleted/deleted.sum(), physical/physical.sum())
        same_sector_tv = tv(same_sector/same_sector.sum(), physical/physical.sum())
        max_deleted_tv = max(max_deleted_tv, deleted_tv)
        max_same_sector_tv = max(max_same_sector_tv, same_sector_tv)
        rows.append(dict(theta=theta, theta_over_pi=theta/np.pi,
                         formula=formula.tolist(), physical=physical.tolist(),
                         sequential=sequential.tolist(), deleted_interference=deleted.tolist(),
                         same_sector=same_sector.tolist(),
                         formula_physical_error=formula_physical_error,
                         formula_sequential_error=formula_sequential_error,
                         deleted_interference_tv=deleted_tv,
                         same_sector_tv=same_sector_tv,
                         deleted_interference_raw_half_l1=deleted_raw_half_l1,
                         same_sector_raw_half_l1=same_raw_half_l1,
                         formula_info=info, deleted_info=deleted_info,
                         same_sector_info=same_info, physical_info=physical_info))
    exp.fail_check("C1", max_deleted_tv > 1e-5,
                   f"max deleted-interference TV from physical={max_deleted_tv:.6g}")
    exp.fail_check("C2", max_same_sector_tv > 1e-5,
                   f"max same-sector-c TV from physical={max_same_sector_tv:.6g}")
    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(N=N, base=BASE, period=PERIOD, block_size=BLOCK,
                             width=WIDTH, split=SPLIT, support_labels=[0, 1],
                             kick_angles=[float(x) for x in THETAS],
                             dense_budget_bytes=MAX_DENSE_BYTES,
                             cost_scope="fixed complete-output formula; O(M*t*b^3), streamed M",
                             not_sampler=True, no_amplitude_cutoff=True,
                             elapsed_seconds=elapsed, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("symmetry_boundary_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
