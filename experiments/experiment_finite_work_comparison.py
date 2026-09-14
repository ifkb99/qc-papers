"""Can verified linear-depth and b=2 scalar proposals replace prefix proposals?

Predictions before this comparison; a finite-work seed624 smoke timing was
already observed, not a speed prediction. P1: both new proposals inside C63
rejection meet the same full-law contract, including W0 Rx(pi/4) sensitivity.
P2: a valid loose enclosure injected AFTER one bit forces rebuild/replay
without changing that bit or querying RNG again. P3: actual wide runs charge
all attempts, forward rebuilds, scalar evaluations and working precision.
H-time: measure three explicit algorithms, no predicted speed ordering.
C1: omitting backward replay preserves the unnormalized prefix-mass oracle
must FAIL (even if a resulting common scalar sometimes cancels in ratios).

Frozen tiny r10,t4 reference guards are inherited from the existing complete
accepted-law evaluator. Only requested accuracy varies within either fixture;
the W0 fixture is a separately declared sensitivity control, not a width
sweep. Wide timing uses the original fixed fixture and three seeds repeated
twice, rotating execution order. No full orbit/output tables or native RSS
claims. Seeded bits are diagnostics, not a source-of-randomness proof.
"""
from __future__ import annotations

import json
import time
import traceback
from fractions import Fraction
from statistics import median

import numpy as np
from flint import acb,acb_mat,arb,ctx
from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_finite_work import VerifiedFiniteWork
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_verified_sampling import fixtures
from experiments.experiment_verified_rejection import run_target,json_safe
from experiments.experiment_rejection_comparison import CountedBits
from experiments.experiment_norm_sampling import wide
from experiments.experiment_coherent_route_sampling import rx,stamp


def main():
    exp = Experiment("finite_work_comparison",doc=__doc__,exit_on_fail=False)
    exp.predict("P1","both verified proposals yield normalized accepted laws within all error/success budgets")
    exp.predict("P2","refinement replays a previously selected bit without redrawing")
    exp.predict("P3","wide executions charge actual work without orbit/output tables")
    exp.predict("H-time","matched three-method timing, without a predicted ordering")
    exp.must_fail("C1","omitting replay preserves unnormalized prefix weights")
    rows = []
    original,float_original = fixtures()
    augmented = VerifiedReflectionCircuit(10,4,{0:("x",Fraction(1,4)),**dict(original.backgrounds)},dict(original.reflections))
    augmented_float = CoherentReflectionCircuit(10,2,4,
        {0:rx(np.pi/4),**float_original.background.defects},dict(float_original.reflections))
    for name,c,f in (("original",original,float_original),("with_W0",augmented,augmented_float)):
        for target in (Fraction(1,1000),Fraction(1,1_000_000)):
            for method in ("finite_work","scalar"):
                row = run_target(c,f,target,component_method=method)
                row.update(series="tiny",fixture=name)
                rows.append(row)
                exp.check("P1", Fraction(row["proposal_mass"]) == Fraction(row["accepted_mass"]) == 1
                    and Fraction(row["target_tv_upper_P192"]) <= Fraction(row["plan"]["total_tv_upper_bound"]) <= target
                    and Fraction(row["target_tv_upper_P256"]) <= target
                    and Fraction(row["proposal_tv_upper_from_component_Arb"]) <= Fraction(row["plan"]["proposal_tv_upper_bound"])
                    and row["exact_success_meets_plan_lower"] and row["ideal_target_interval_total_contains_one"]
                    and row["ideal_target_independent_max_error"] < 2e-14
                    and all(d["independent_max_error"] < 2e-14 and d["ball_total_contains_one"]
                            for d in row["component_independent_diagnostics"]),
                    f"{name} {method} target={target}: outward TV={float(Fraction(row['target_tv_upper_P256'])):.4g}")
    worker = VerifiedFiniteWork(augmented,2)
    cursor = worker.cursor(1,accuracy_bits=30)
    cursor.weights()
    cursor.advance(0)
    with ctx.workprec(cursor.precision):
        # Valid deliberately loose enclosure of the SAME exact effect.
        cursor.effect += acb_mat([[acb(arb(0,1)),0],[0,acb(arb(0,1))]])
    rebuilt = cursor.weights()
    reference = worker.cursor(1,accuracy_bits=30,initial_precision=cursor.precision)
    reference.weights(); reference.advance(0)
    expected = reference.weights()
    exp.check("P2", cursor.builds > 1 and cursor.replayed >= 1 and cursor.measured == 1
              and cursor.output == 0 and rebuilt == expected,
              f"builds={cursor.builds}, replay steps={cursor.replayed}, P={cursor.precision}")
    broken = worker.cursor(1,accuracy_bits=30,initial_precision=cursor.precision)
    broken.weights(); broken.advance(0)
    broken.effect = acb_mat([[1,0],[0,1]])  # WRONG: discard selected earlier effect
    wrong = broken.weights()
    discrepancy = sum(abs(a-b) for a,b in zip(wrong,expected))
    exp.fail_check("C1", discrepancy > 4,
                   f"wrong replay mass differs by {discrepancy}/2^30, beyond two per-weight errors")
    rows.append(dict(series="replay",rebuilds=cursor.builds,replayed_effect_steps=cursor.replayed,
        correct_weights=expected,wrong_weights=wrong,precision=cursor.precision))

    methods = ("prefix","finite_work","scalar")
    times = {m:[] for m in methods}
    samples = {m:[] for m in methods}
    for repeat in range(2):
        for seed in (624,625,626):
            offset = (repeat+seed)%3
            for method in methods[offset:]+methods[:offset]:
                rng = CountedBits(seed)
                start = time.perf_counter()
                result = VerifiedRejectionSampler(wide("rectangular"),component_method=method).sample(rng)
                elapsed = time.perf_counter()-start
                times[method].append(elapsed); samples[method].append(result)
                rows.append(dict(series="wide",repeat=repeat,seed=seed,method=method,
                    elapsed_seconds=elapsed,random_bits_requested=rng.bits,random_word_calls=rng.calls,result=result))
    exp.check("P3", all(r["total_tv_upper_bound"] <= Fraction(1,1_000_000)
              and not r["orbit_or_output_tables"] and r["attempt_cap"] is None
              and r["acceptance_enclosure_evaluations"] >= r["attempts"] for v in samples.values() for r in v)
              and all(r["prefix_vector_evaluations"] == 0 and r["component_forward_steps"] >= 63*r["attempts"] for r in samples["finite_work"])
              and all(r["prefix_vector_evaluations"] == 0 and r["component_scalar_evaluations"] >= 63*r["attempts"] for r in samples["scalar"]),
              f"attempts by method: { {m:[r['attempts'] for r in v] for m,v in samples.items()} }")
    medians = {m:median(v) for m,v in times.items()}
    exp.check("H-time",all(len(v)==6 and all(t>0 for t in v) for v in times.values()),f"medians seconds={medians}")
    path = stamp("finite_work_comparison")
    ok = exp.finish(report_path=path,rows=json_safe(rows),metadata=dict(
        medians_seconds=medians,timings_seconds=times,optional_backend="python-flint==0.9.0",
        runtime_theorem=False,native_memory_measured=False,
        scope="same requested joint TV; b=2 late backgrounds are output-invisible"))
    print(json.dumps(dict(report=str(path),ok=ok)))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        path = stamp("finite_work_comparison_failure")
        path.write_text(json.dumps(dict(error=repr(error),traceback=traceback.format_exc()),indent=2)+"\n")
        raise
