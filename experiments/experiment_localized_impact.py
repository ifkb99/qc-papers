"""Impact check: does finite basis support become negligible at large order?

AUDIT: intentionally exits nonzero for the refuted monotone-dilution intuition
P2. The first agent run left P2 unresolved; its second compared only endpoints
(the last was a binary resonance). Main replaced that inadequate predicate by
an adjacent-row check. The original reports are preserved; this is a recorded
negative result, not a failing implementation regression.

PREDICTIONS WRITTEN BEFORE MEASUREMENT (TODO 19 follow-up).  For a uniform
pre-defect prefix of length L and a defect differing from identity only on
support S, let

    F = sum_{p in S} n_p/L,
    n_p = #{0 <= l < L : l mod r = p}.

The clean and defective pre-QFT states differ only on those F branches.  Since
||(V-I)|p>|| <= 2 for a unitary V, their state-vector distance is at most
2 sqrt(F), and measurement total variation distance obeys

    TV(output_defect, output_ideal) <= min(1, 2 sqrt(F)).

This bounded abstract-order sweep fixes t=8,s=5 and uses the same endpoint
two-level mixer at r=3,...,16.  It compares the existing spectral reference
and SparseOrbitPrefix; it does not compile or propagate a large-r circuit.

P1: spectral and SparseOrbitPrefix marginals agree, and every interior-row TV
    obeys the independently calculated min(1,2 sqrt(F)) bound.
P2 (exploratory intuition, NOT derived): as r grows at fixed L, the observed
    ideal-sampler TV is nonincreasing. Tested on every adjacent pair; resonant
    exact cancellations refute this intuition. F itself need not tend to zero
    when L is fixed: the endpoint p=0 always has positive branch weight.
P3: zero angle is exactly ideal, while placing the same unitary at endpoint
    insertion s=t is exactly invisible after tracing out work.

C1: replacing every interior defect output by the ideal sampler must fail on a
    nonzero-angle, non-endpoint row.

The order and endpoint indices are supplied inputs.  Dense references are
guarded before allocation at 16 MiB and are limited to t=8.  This tests the
impact of the strongest simple ideal-sampler approximation, not a runtime or
order-discovery comparison.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_localized_impact
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects


MAX_DENSE_BYTES = 16 * 1024 * 1024
WIDTH = 8
SPLIT = 5


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def mixer(theta: float) -> np.ndarray:
    guard((2, 2), np.complex128, "endpoint mixer")
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex)


def embedded(period: int, indices: tuple[int, ...], block: np.ndarray) -> np.ndarray:
    guard((period, period), np.complex128, "abstract defect")
    result = np.eye(period, dtype=complex)
    if indices:
        result[np.ix_(indices, indices)] = block
    return result


def spectral_marginal(v: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    guard((period, period), np.complex128, "orbit Fourier")
    j = np.arange(period, dtype=float)
    fourier = np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)
    # single_defect_effects has its own entry cap; this outer guard records
    # the same payload discipline for the input and expected output shape.
    guard((1 << width, period, period), np.complex128, "spectral effects")
    effects = single_defect_effects(period, width, split,
                                    fourier.conj().T @ v @ fourier,
                                    max_entries=1_000_000)
    return effects.sum(axis=(1, 2)).real / period


def prefix_marginal(v: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    length = 1 << split

    def column(l):
        values = v[:, l % period]
        ids = np.flatnonzero(values != 0).astype(np.int64)
        return ids, values[ids]

    prefix = SparseOrbitPrefix(period, split, column,
                               max(1, int(np.max(np.count_nonzero(v != 0, axis=0)))))
    guard((1 << width,), np.float64, "prefix marginal")
    marginal = np.zeros(1 << width, dtype=float)
    for k in range(period):
        for y in range(1 << width):
            marginal[y] += prefix.forced_joint(width, k, y)["joint_latent_output_probability"]
    return marginal


def hit_fraction(period: int, indices: tuple[int, ...], split: int) -> tuple[float, tuple[int, ...]]:
    length = 1 << split
    counts = tuple(0 if p >= length else 1 + (length - 1 - p) // period for p in indices)
    return float(sum(counts) / length), counts


def tv_distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.abs(a - b).sum() / 2)


def unique_path(prefix: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def main() -> None:
    exp = Experiment("localized_impact", doc=__doc__)
    exp.predict("P1", "spectral and SparseOrbitPrefix marginals agree and TV obeys min(1,2sqrt(F))")
    exp.predict("P2", "exploratory monotone-dilution intuition: ideal-sampler TV is nonincreasing across adjacent orders")
    exp.predict("P3", "zero-angle and endpoint-insertion controls are exactly ideal")
    exp.must_fail("C1", "the ideal sampler is not exact for every interior nonzero defect")

    support_by_order = lambda r: (0, r - 1)
    theta = np.pi / 2
    rows = []
    max_reference_error = 0.
    max_bound_violation = 0.
    max_ideal_approximation = 0.
    interior_tvs = []

    # One parameter varies here: only the abstract period r changes.
    for period in range(3, 17):
        indices = support_by_order(period)
        defect = embedded(period, indices, mixer(theta))
        actual = spectral_marginal(defect, period, WIDTH, SPLIT)
        prefix = prefix_marginal(defect, period, WIDTH, SPLIT)
        ideal = spectral_marginal(np.eye(period, dtype=complex), period, WIDTH, SPLIT)
        reference_error = max(float(np.max(np.abs(actual - prefix))),
                              float(abs(actual.sum() - 1.)))
        tv = tv_distance(actual, ideal)
        fraction, counts = hit_fraction(period, indices, SPLIT)
        bound = min(1., 2. * math.sqrt(fraction))
        violation = max(0., tv - bound)
        max_reference_error = max(max_reference_error, reference_error)
        max_bound_violation = max(max_bound_violation, violation)
        max_ideal_approximation = max(max_ideal_approximation, tv)
        interior_tvs.append(tv)
        exp.check("P1", reference_error < 2e-10 and violation < 2e-10,
                  f"r={period}: ref={reference_error:.2e}, F={fraction:.6g}, TV={tv:.6g}, bound={bound:.6g}")
        rows.append(dict(kind="period_sweep", period=period, width=WIDTH, split=SPLIT,
                         support=list(indices), counts=list(counts), hit_fraction=fraction,
                         tv_defect_vs_ideal=tv, tv_bound=bound, bound_violation=violation,
                         reference_error=reference_error, ideal_sampler_approximation_tv=tv))

    # Fixed controls: angle zero and endpoint placement s=t.  These are not
    # folded into the one-parameter r sweep.
    # r=6 has a visible INTERIOR effect in the sweep; r=8 was already invisible
    # there and would make endpoint invariance a less discriminating control.
    period = 6
    indices = support_by_order(period)
    zero = spectral_marginal(embedded(period, indices, mixer(0.)), period, WIDTH, SPLIT)
    ideal = spectral_marginal(np.eye(period, dtype=complex), period, WIDTH, SPLIT)
    zero_tv = tv_distance(zero, ideal)
    exp.check("P3", zero_tv < 2e-11, f"zero angle r={period},s=5: TV={zero_tv:.2e}")
    rows.append(dict(kind="zero_angle", period=period, width=WIDTH, split=SPLIT,
                     tv_defect_vs_ideal=zero_tv))

    endpoint = spectral_marginal(embedded(period, indices, mixer(theta)), period, WIDTH, WIDTH)
    ideal_endpoint = spectral_marginal(np.eye(period, dtype=complex), period, WIDTH, WIDTH)
    endpoint_tv = tv_distance(endpoint, ideal_endpoint)
    exp.check("P3", endpoint_tv < 2e-10, f"endpoint insertion r={period},s=t=8: TV={endpoint_tv:.2e}")
    rows.append(dict(kind="endpoint_insertion", period=period, width=WIDTH, split=WIDTH,
                     tv_defect_vs_ideal=endpoint_tv))

    increases = [(r, r+1, float(b-a)) for r, (a, b) in
                 enumerate(zip(interior_tvs, interior_tvs[1:]), start=3) if b > a+1e-12]
    exp.check("P2", not increases, f"adjacent TV increases contradict monotone dilution: {increases}")

    # The must-fail control is deliberately tested only after the positive
    # rows establish nonzero interior discrepancies.
    max_interior = max(interior_tvs)
    exp.fail_check("C1", max_interior > 1e-7,
                   f"largest interior ideal-sampler TV={max_interior:.6g}")

    # A compact impact summary makes the simple-baseline question explicit.
    small_thresholds = {str(threshold): [row["period"] for row in rows
                                         if row.get("kind") == "period_sweep"
                                         and row["tv_defect_vs_ideal"] <= threshold]
                        for threshold in (1e-1, 1e-2, 1e-3)}
    path = unique_path("localized_impact")
    exp.finish(report_path=path, rows=rows,
               metadata=dict(width=WIDTH, split=SPLIT, theta_over_pi=.5,
                             support="endpoint indices (0,r-1), d=2",
                             max_reference_error=max_reference_error,
                             max_bound_violation=max_bound_violation,
                             max_interior_tv=max_interior,
                             tv_r3=interior_tvs[0], tv_r16=interior_tvs[-1],
                             tv_is_monotone_nonincreasing=all(a >= b - 1e-12
                                                              for a, b in zip(interior_tvs, interior_tvs[1:])),
                             threshold_crossings=small_thresholds,
                             interpretation="ideal output is an approximation; no order discovery or runtime claim",
                             dense_budget_bytes=MAX_DENSE_BYTES, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = unique_path("localized_impact_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc()), indent=2) + "\n")
        raise
