"""Does a certified rejection comparator remain useful after charging proposals?

PREDICTIONS BEFORE MEASUREMENT (TEMPLATE.py workflow).
P1: Exact-zero and increasingly small proposal masses terminate under the
    ABSOLUTE interval-width rule, and downward acceptance cannot over-accept.
    For r6,t2 with final K(q1,theta), gamma0,y1, the identity history vanishes;
    P=sin(theta/2)^2*n and D=(|cos|+|sin|)|sin|*n, n>0. Hence theta->0
    creates a small denominator but not a need for relative accuracy.
P2: Both verified samplers obey the same requested joint-output TV promise
    on the supplied-wide fixed fixture, with all attempts/refinements charged.
H-time: Rejection pays retries but contracts one history per proposal prefix.
    Either route can win: record matched alternating-order repeated timings,
    without predicting a speedup or inferring worst-case runtime from seeds.
C1: Flooring tiny positive acceptance to zero is NOT an exact local sampler;
    the nonzero discrepancy must be charged by accepted mass, not concealed.
C2: An uncharged finite proposal cap cannot inherit the uncapped certificate.

No new propagator, no dense arrays or full-r reference in this experiment.
The tiny FULL-LAW comparator lives in experiment_verified_rejection. Here
fixed-size scalar objects and at most twelve timing rows are retained; native
backend memory is not measured or inferred from scalar counts. Timing starts
before sampler construction and covers planning, history setup, all failed
attempts, acceptance and refinement. PRNG bits are only a diagnostic driver.
"""
from __future__ import annotations

import json
import time
import traceback
from fractions import Fraction
from random import Random
from statistics import median

from flint import ctx
from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_norm_sampling import wide
from experiments.experiment_coherent_route_sampling import stamp


class CountedBits:
    def __init__(self, seed):
        self.driver, self.bits, self.calls = Random(seed), 0, 0

    def getrandbits(self, bits):
        self.bits += bits
        self.calls += 1
        return self.driver.getrandbits(bits)


def simple(value):
    return {k: str(v) if isinstance(v, Fraction) else v for k, v in value.items()}


def main():
    exp = Experiment("rejection_comparison", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "zero/near-zero absolute interval decisions terminate and under-accept")
    exp.predict("P2", "both wide samplers carry the same requested TV contract and charge all work")
    exp.predict("H-time", "report matched runtime distribution; no predicted speed ordering")
    exp.must_fail("C1", "zeroing tiny positive acceptance is exact")
    exp.must_fail("C2", "an uncharged finite attempt cap inherits the uncapped certificate")
    rows, underaccepted = [], []
    tolerance, bits = Fraction(1, 10**9), 30
    for power in (None, 4, 16, 64, 256):
        theta = 0 if power is None else Fraction(1, 1 << power)
        sampler = VerifiedRejectionSampler(VerifiedReflectionCircuit(6, 2, {}, {2: (1, theta)}))
        decision = sampler.acceptance(0, 1, interval_tolerance=tolerance, random_bits=bits)
        P, D = sampler.acceptance_enclosure(0, 1, working_precision=768)
        with ctx.workprec(768):
            pl, pu = binary_fraction(P.lower()), binary_fraction(P.upper())
            dl, du = binary_fraction(D.lower()), binary_fraction(D.upper())
        a = decision.probability
        # Strong exact rational sufficient conditions, not float magnitudes.
        lower_deficit = pl-du*a
        upper_deficit = pu-dl*a
        bound = decision.width_sum+du*Fraction(1, 1 << bits)
        zero = power is None and P.is_zero() and D.is_zero() and a == 0
        positive = power is not None and pl > 0 and dl > 0
        exp.check("P1", decision.width_sum <= tolerance and (zero or positive)
                  and lower_deficit >= 0 and upper_deficit <= bound,
                  f"angle pi/{'infinity' if power is None else '2^'+str(power)}: "
                  f"P={decision.working_precision}, retries={decision.refinements}, a={a}")
        if positive and a == 0:
            underaccepted.append(power)
        rows.append(dict(series="near_zero", angle=str(theta), power=power,
            acceptance=str(a), numerator_high_bounds=[str(pl),str(pu)],
            denominator_high_bounds=[str(dl),str(du)], width_sum=str(decision.width_sum),
            upper_deficit=str(upper_deficit), deficit_bound=str(bound),
            working_precision=decision.working_precision, refinements=decision.refinements))
    exp.fail_check("C1", bool(underaccepted),
                   f"strictly positive ideal acceptance rounded to zero at powers {underaccepted}; charged, not exact")

    timings = {"prefix": [], "rejection": []}
    records = {name: [] for name in timings}
    target = Fraction(1, 1_000_000)
    for repeat in range(2):
        for seed in (624, 625, 626):
            names = ("prefix", "rejection") if (repeat+seed)%2 == 0 else ("rejection", "prefix")
            for name in names:
                rng = CountedBits(seed)
                started = time.perf_counter()
                circuit = wide("rectangular")
                result = (circuit.sample(rng, target_tv=target) if name == "prefix" else
                          VerifiedRejectionSampler(circuit).sample(rng, target_tv=target))
                elapsed = time.perf_counter()-started
                timings[name].append(elapsed)
                records[name].append(result)
                rows.append(dict(series="wide", method=name, repeat=repeat, seed=seed,
                    elapsed_seconds=elapsed, random_bits_requested=rng.bits,
                    random_word_calls=rng.calls, result=simple(result)))
    all_results = [r for values in records.values() for r in values]
    exp.check("P2", all(r["total_tv_upper_bound"] <= r["target_tv"] == target
              and not r["orbit_or_output_tables"] and r["precision_cap"] is None
              for r in all_results)
              and all(r["acceptance_enclosure_evaluations"] >= r["attempts"] >= 1
                      and r["component_history_count"] == 1 and r["attempt_cap"] is None
                      for r in records["rejection"]),
              f"rejection attempts: {[r['attempts'] for r in records['rejection']]}")
    medians = {name: median(values) for name, values in timings.items()}
    exp.check("H-time", all(len(v) == 6 and all(t > 0 for t in v) for v in timings.values()),
              f"median seconds={medians}; rejection/prefix={medians['rejection']/medians['prefix']:.3g}")
    # Exact two-cell model: target(1/2,1/2), q(1/4,3/4), C2, acceptance(1,1/3).
    # For a one-attempt cap with fallback cell0, law=(3/4,1/4), TV=1/4.
    capped_law = (Fraction(3,4), Fraction(1,4))
    cap_tv = sum(abs(v-Fraction(1,2)) for v in capped_law)/2
    exp.fail_check("C2", cap_tv == Fraction(1,4) > target,
                   f"exact ideal acceptance plus uncharged fallback has TV={cap_tv}")
    path = stamp("rejection_comparison")
    ok = exp.finish(report_path=path, rows=rows, metadata=dict(
        optional_backend="python-flint==0.9.0", medians_seconds=medians,
        rejection_over_prefix_median=medians["rejection"]/medians["prefix"],
        times_seconds=timings, same_requested_joint_output_tv=str(target),
        component_proposal="C61 prefix sampler, NOT a verified O(t) finite-work instrument",
        native_memory_measured=False, seeded_rng_proves_unbiasedness=False,
        no_runtime_or_novelty_theorem=True))
    print(json.dumps(dict(report=str(path),ok=ok)))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        path = stamp("rejection_comparison_failure")
        path.write_text(json.dumps(dict(error=repr(error), traceback=traceback.format_exc()), indent=2)+"\n")
        raise
