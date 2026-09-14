"""Do norm-aware enclosures preserve the full certified law at lower cost?

PREDICTIONS BEFORE THIS COMPARISON (after a disclosed seed624 smoke test).
P1: Both modes' complete exact tiny transition laws normalize and obey the
    same rational TV budgets as requested accuracy varies, fixed exact gates.
P2: On a fixed wide fixture and seed set, norm mode uses fewer refinement
    retries/working bits than rectangular mode without changing the input or
    accuracy contract. This is measured, not a uniform precision theorem.
H-time: extra norm bookkeeping may outweigh saved precision/retries. No speed
    claim is predicted; report matched alternating-order repeated timings.
C1: Merely discarding local-error radii cannot certify midpoint arithmetic;
    an independent norm-prefix probe owns that control. Here a claimed
    'doubling endpoint equals minimum sufficient precision' must FAIL on the
    frozen query identified by experiment_precision_growth.

Uses EXISTING exact tiny Markov enumeration, capped at r10,t4, and outward
final-ball TV, not another state propagator. The wide runs have no frontier,
trace or orbit/output table. Timing includes local ball/radius computation
and refinement; Python/Arb native memory is not inferred from scalar counts.
"""
from __future__ import annotations

import json
import time
import traceback
from fractions import Fraction
from random import Random
from statistics import median

from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction
from experiments.experiment_verified_sampling import fixtures, exact_transition_law, enclosed_tv
from experiments.experiment_coherent_route_sampling import stamp


def wide(mode):
    return VerifiedReflectionCircuit((1 << 61)-2, 63,
        {1: ("x", Fraction(1, 7)), 31: ("z", Fraction(1, 5)), 63: ("x", Fraction(1, 7))},
        {2: (0, Fraction(1, 5)), 32: (1, Fraction(1, 5)), 63: (2, Fraction(1, 5))},
        enclosure_mode=mode)


def main():
    exp = Experiment("norm_sampling", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "both modes' exact finite-bit laws meet identical requested TV budgets")
    exp.predict("P2", "norm mode reduces working bits/retries on the frozen wide seed set")
    exp.predict("H-time", "report matched timing without assuming fewer bits means faster execution")
    exp.must_fail("C1", "the adaptive doubling endpoint is the minimum sufficient mantissa")
    start = time.perf_counter()
    tiny, _ = fixtures()
    circuits = {mode: VerifiedReflectionCircuit(10, 4, dict(tiny.backgrounds),
                dict(tiny.reflections), enclosure_mode=mode) for mode in ("rectangular", "norm")}
    rows = []
    for target in (Fraction(1, 1000), Fraction(1, 1_000_000), Fraction(1, 10**12)):
        laws = {}
        for mode, circuit in circuits.items():
            plan = circuit.accuracy_plan(target)
            law, stats = exact_transition_law(circuit, plan["absolute_accuracy_bits"], plan["random_bits"])
            error, _ = enclosed_tv(circuit, law)
            laws[mode] = law
            rows.append(dict(series="tiny", mode=mode, target=str(target), mass=str(sum(law.values())),
                outward_tv_upper=str(error), outward_tv_upper_float=float(error),
                planned_bound=str(plan["total_tv_upper_bound"]), **stats))
            exp.check("P1", sum(law.values()) == 1 and error <= plan["total_tv_upper_bound"] <= target,
                      f"{mode} target={target}, outward TV<={float(error):.6g}")
        rows.append(dict(series="exact_law_comparison", target=str(target), identical=laws["rectangular"] == laws["norm"]))
    wides = {mode: wide(mode) for mode in circuits}
    timings = {mode: [] for mode in circuits}
    samples = {mode: [] for mode in circuits}
    for repeat in range(2):
        for seed in (624, 625, 626):
            order = ("rectangular", "norm") if (repeat+seed) % 2 == 0 else ("norm", "rectangular")
            for mode in order:
                started = time.perf_counter()
                result = wides[mode].sample(Random(seed), target_tv=Fraction(1, 1_000_000))
                elapsed = time.perf_counter()-started
                timings[mode].append(elapsed)
                samples[mode].append(result)
                rows.append(dict(series="wide", repeat=repeat, seed=seed, mode=mode,
                    elapsed_seconds=elapsed, result={k: str(v) if isinstance(v,Fraction) else v for k,v in result.items()}))
    exp.check("P2", all(n["max_working_precision"] < r["max_working_precision"]
                         and n["refinement_retries"] < r["refinement_retries"]
                         and n["total_tv_upper_bound"] == r["total_tv_upper_bound"]
                         for n,r in zip(samples["norm"],samples["rectangular"])),
              f"max P: rectangular={max(s['max_working_precision'] for s in samples['rectangular'])}, "
              f"norm={max(s['max_working_precision'] for s in samples['norm'])}")
    medians = {mode: median(values) for mode,values in timings.items()}
    exp.check("H-time", all(len(values) == 6 and all(t > 0 for t in values) for values in timings.values()),
              f"median seconds={medians}; norm/rectangular={medians['norm']/medians['rectangular']:.3g}")
    label = (1038651310644072125, 2246884356896903187, 63)
    kw = dict(boundary="reflection", measured=45, output=16994140961528)
    result = wides["rectangular"].prefix_dyadic(*label, accuracy_bits=61, **kw)
    ball = wides["rectangular"].prefix_enclosure(*label, working_precision=78, **kw)
    radius = max(binary_fraction(z.rad().upper()) for v in ball for z in (v.real,v.imag))
    exp.fail_check("C1", radius <= Fraction(1, 1 << 62) and result.working_precision > 78,
                   f"78 bits already meet radius tolerance; adaptive endpoint={result.working_precision}")
    path = stamp("norm_sampling")
    ok = exp.finish(report_path=path, rows=rows, metadata=dict(
        optional_backend="python-flint==0.9.0", times_seconds=timings, medians_seconds=medians,
        norm_over_rectangular_median=medians["norm"]/medians["rectangular"],
        runtime_theorem=False, seed_set_worst_case_proof=False, elapsed_seconds=time.perf_counter()-start,
        scope="same exact-input and TV promise; enclosure strategy comparison, not rejection benchmark"))
    print(json.dumps(dict(report=str(path),ok=ok)))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        path = stamp("norm_sampling_failure")
        path.write_text(json.dumps(dict(error=repr(error),traceback=traceback.format_exc()),indent=2)+"\n")
        raise
