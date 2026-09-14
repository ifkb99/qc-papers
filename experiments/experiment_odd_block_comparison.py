"""Same-accuracy rejection comparison beyond the binary scalar-tail fixture.

TEMPLATE predictions registered before measurement. P1: the fixed r9,b3,t4
complete accepted laws with prefix and finite-work proposals meet the same
target TV at 1e-3 and 1e-6, with all histories, success and precision charged.
P2: supplied-wide executions remain table-free and report actual rejected
attempts, finite words and arithmetic refinement. H-time: no speed ordering
predicted. C1: deleting the tiny fixture's late backgrounds preserves its
target law must FAIL, now using outward rational intervals rather than only
the already-known floating-point discrepancy.

Only accuracy varies inside the tiny fixture. The wide comparison is a
SEPARATE supplied indexed circuit: r=3*(2^60-1), b3,t63, backgrounds at16,32,48,
reflections at21,42. It is not a width-scaling fit or a claim that the tiny
effect magnitude persists. Tiny full-law gates must pass before timing starts.
Existing r<=10,t<=4 diagnostic guards apply; no dense reference is attempted
on the wide input. Counts of matrices are not native/RSS memory measurements.
Three seeded traces per method are each repeated twice with alternating order;
these are not six independent rejection traces. Setup and every rejected
attempt are timed. The input order/index is supplied, never discovered.
"""
from __future__ import annotations

import json
import time
import traceback
from fractions import Fraction
from statistics import median

from lab import Experiment
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_odd_block_tail import fixture
from experiments.experiment_verified_rejection import run_target, target_interval, json_safe
from experiments.experiment_coherent_route_sampling import stamp
from experiments.experiment_rejection_comparison import CountedBits


def tiny():
    return VerifiedReflectionCircuit(9,4,
        {1:("x",Fraction(1,7)),3:("z",Fraction(1,5))},
        {2:(0,Fraction(1,5)),3:(1,Fraction(1,5))},block_size=3)


def wide():
    return VerifiedReflectionCircuit(3*((1 << 60)-1),63,
        {16:("x",Fraction(1,7)),32:("z",Fraction(1,5)),48:("x",Fraction(1,7))},
        {21:(0,Fraction(1,5)),42:(1,Fraction(1,5))},block_size=3)


def main():
    exp = Experiment("odd_block_comparison",doc=__doc__,exit_on_fail=False)
    exp.predict("P1","both complete finite-bit accepted laws obey all same-accuracy budgets")
    exp.predict("P2","supplied-wide execution charges setup, rejected attempts, finite bits and refinements")
    exp.predict("H-time","matched timings with no predicted method ordering")
    exp.must_fail("C1","removing the tiny odd-block late backgrounds preserves the target law")
    rows=[]
    c, f = tiny(), fixture(3)
    for target in (Fraction(1,1000),Fraction(1,1_000_000)):
        for method in ("prefix","finite_work"):
            row = run_target(c,f,target,component_method=method)
            row["series"]="tiny"
            rows.append(row)
            ok = (Fraction(row["proposal_mass"]) == Fraction(row["accepted_mass"]) == 1
                and max(Fraction(row["target_tv_upper_P192"]),Fraction(row["target_tv_upper_P256"]))
                    <= Fraction(row["plan"]["total_tv_upper_bound"]) <= target
                and Fraction(row["proposal_tv_upper_from_component_Arb"])
                    <= Fraction(row["plan"]["proposal_tv_upper_bound"])
                and row["exact_success_meets_plan_lower"] and row["ideal_target_interval_total_contains_one"]
                and row["ideal_target_independent_max_error"] < 2e-14
                and abs(row["independent_target_mass"]-1) < 2e-14
                and all(d["independent_max_error"] < 2e-14 and d["ball_total_contains_one"]
                        and abs(d["full_r_mass"]-1) < 2e-14
                        for d in row["component_independent_diagnostics"]))
            exp.check("P1",ok,f"{method} target={target}; outward TV={float(Fraction(row['target_tv_upper_P256'])):.4g}")
            if not ok:
                raise AssertionError("tiny full-law gate failed; wide timing not started")
    removed = VerifiedReflectionCircuit(9,4,{},dict(c.reflections),block_size=3)
    lower=Fraction(0)
    for g in range(c.sectors):
        for y in range(1 << c.width):
            alo,ahi=target_interval(c,g,y,192)
            blo,bhi=target_interval(removed,g,y,192)
            lower += max(Fraction(0),alo-bhi,blo-ahi)/2
    exp.fail_check("C1",lower>Fraction(1,1000),f"outward late-removal TV lower={float(lower):.8g}")
    if lower <= Fraction(1,1000):
        raise AssertionError("visibility control failed; do not benchmark an unverified fixture")

    methods=("prefix","finite_work")
    times={method:[] for method in methods}
    samples={method:[] for method in methods}
    for repeat in range(2):
        for seed in (651,652,653):
            order=methods if (repeat+seed)%2 else methods[::-1]
            for method in order:
                rng=CountedBits(seed)
                started=time.perf_counter()
                result=VerifiedRejectionSampler(wide(),component_method=method).sample(rng)
                elapsed=time.perf_counter()-started
                times[method].append(elapsed)
                samples[method].append(result)
                rows.append(dict(series="wide",repeat=repeat,seed=seed,method=method,
                    seconds=elapsed,random_bits_requested=rng.bits,random_word_calls=rng.calls,result=result))
    valid=all(s["total_tv_upper_bound"]<=Fraction(1,1_000_000)
        and not s["orbit_or_output_tables"] and s["attempt_cap"] is None and s["precision_cap"] is None
        and s["acceptance_enclosure_evaluations"]>=s["attempts"]>=1
        for values in samples.values() for s in values)
    valid &= all(s["prefix_vector_evaluations"]==0
        and s["component_forward_steps"]>=63*s["attempts"] for s in samples["finite_work"])
    exp.check("P2",valid,f"attempts={ {m:[s['attempts'] for s in ss] for m,ss in samples.items()} }")
    medians={m:median(v) for m,v in times.items()}
    exp.check("H-time",all(len(v)==6 and all(x>0 for x in v) for v in times.values()),f"median seconds={medians}")
    path=stamp("odd_block_comparison")
    ok=exp.finish(report_path=path,rows=json_safe(rows),metadata=dict(
        medians_seconds=medians,timing_seconds=times,late_removal_tv_lower=str(lower),
        runtime_theorem=False,native_memory_measured=False,optional_backend="python-flint==0.9.0",
        benchmark_scope="same certified b3 joint target; scalar shortcut rejected, no optimality claim"))
    print(json.dumps(dict(ok=ok,report=str(path))))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        path=stamp("odd_block_comparison_failure")
        path.write_text(json.dumps(dict(error=repr(error),traceback=traceback.format_exc()),indent=2)+"\n")
        raise
