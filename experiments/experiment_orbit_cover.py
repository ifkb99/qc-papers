"""Exact interval arithmetic for the co-moving localized-kick certificate.

PREDICTIONS WRITTEN BEFORE MEASUREMENT.
P1: merging circular radius-R neighborhoods counts exactly the brute-force
    eligible early exponent labels, including wrapping and incomplete cycles.
P2: the helper handles width 63 and large supplied r with O(d) intervals;
    zero support gives zero error, full coverage gives only the trivial bound.
P3: the background method counts only mixers before (or at) the specified
    kick, and its rational squared-TV bound remains separate from float
    sampling error. A seeded wide background draw needs no new state arrays.
C1: dropping the wrapped part of an interval gives a wrong cover count.
C2: replacing a truncated traversal by uniform residues gives a wrong count.

The brute force is VALIDATION ONLY, r<=18,width<=7,radius<=19. It stores
no state vector or propagator. Production is exact integer interval arithmetic,
not certification of complex128 sampling or an exact broken-symmetry sampler.
Run from research/: uv run --no-project --python 3.12 --with 'numpy<2.5'
python -m experiments.experiment_orbit_cover
"""
from __future__ import annotations

from datetime import datetime, timezone
import numpy as np
from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit, lightcone_cover


def brute_count(period, width, support, radius):
    if not 1 <= period <= 18 or not 0 <= width <= 7 or not 0 <= radius <= 19:
        raise ValueError("brute-force validation cap")
    return sum(any(min((l-p) % period, (p-l) % period) <= radius for p in support)
               for l in range(1 << width))


def main():
    exp = Experiment("orbit_cover", doc=__doc__)
    exp.predict("P1", "exact merged-interval counts agree with capped brute force")
    exp.predict("P2", "integer edge cases and large supplied period avoid tables")
    exp.predict("P3", "background insertion convention and wide sample remain bounded")
    exp.must_fail("C1", "dropping the wrapped tail changes the true count")
    exp.must_fail("C2", "uniform residue mass loses the incomplete traversal")
    cases, max_error = 0, 0
    for period in range(1, 19):
        supports = {(), (0,), (period-1,), tuple(sorted({0, period-1})),
                    tuple(range(period))}
        for width in range(8):
            for radius in sorted({0, 1, 2, period//2, period, 19}):
                for support in sorted(supports):
                    got = lightcone_cover(period, width, support, radius)
                    want = brute_count(period, width, support, radius)
                    max_error = max(max_error, abs(got["covered_prefix_labels"]-want))
                    cases += 1
    exp.check("P1", max_error == 0, f"{cases} exact integer count cases; max error={max_error}")
    empty = lightcone_cover(3_000_000_021, 63, (), 100)
    whole = lightcone_cover(3_000_000_021, 63, (0,), 3_000_000_021)
    exp.check("P2", empty["covered_prefix_labels"] == 0
              and whole["covered_prefix_labels"] == 1 << 63,
              "empty and whole-orbit support at width 63 have exact expected counts")
    X = np.array([[0., 1., 0.], [1., 0., 0.], [0., 0., 0.]])
    W = np.eye(3, dtype=complex)
    W[:2, :2] = np.cos(.7/2)*np.eye(2)-1j*np.sin(.7/2)*X[:2, :2]
    circuit = PeriodicOrbitCircuit(3_000_000_021, 3, 63, {0: W, 16: W, 48: W})
    cover = circuit.localized_kick_bound(32, (0, 1))
    sample = circuit.sample(np.random.default_rng(5721))
    exp.check("P3", cover["radius"] == 4 and cover["pre_kick_mixer_count"] == 2
              and cover["covered_prefix_labels"] == 16
              and cover["tv_squared_numerator"] == 64
              and cover["tv_squared_denominator"] == 1 << 32
              and len(cover["intervals"]) == 2
              and not cover["numerical_sampler_error_included"]
              and sample["matrix_dimension"] == 3,
              "two preceding blocks: exact mathematical TV bound 2^-13; width-63 background draw")
    wrap = lightcone_cover(9, 4, (0,), 1)
    wrong_wrap_count = sum(l % 9 < 2 for l in range(16))
    exp.fail_check("C1", wrap["covered_prefix_labels"] != wrong_wrap_count,
                   f"true count={wrap['covered_prefix_labels']}; dropped wrap={wrong_wrap_count}")
    incomplete = lightcone_cover(9, 2, (0,), 0)
    exp.fail_check("C2", incomplete["covered_prefix_labels"]*9 != 4,
                   "true hit mass=1/4; incorrectly uniform mass=1/9")
    path = "out/orbit_cover_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")+".json"
    exp.finish(report_path=path, metadata=dict(count_cases=cases,
               count_max_integer_error=max_error, wide_cover=cover, wide_sample=sample,
               known_order_and_index=True, broken_symmetry_sample=False))
    print(f"report: {path}")


if __name__ == "__main__":
    main()
