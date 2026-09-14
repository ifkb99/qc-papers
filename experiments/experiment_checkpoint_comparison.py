"""Exact-law comparison for on-demand forward checkpoints.

The fixed primary input is r=9,b=3,t=4 with embedded Rx(pi/7) at W1..W4
and two reflection insertions (four histories). For every mask, initial sector, and forced output,
``VerifiedFiniteWork`` is run with checkpoint spacing None, 1, 2, and 4 at
target TV 1e-6.  Laws must agree as exact Fractions; checkpointing changes
only reconstruction work and retained-storage counters.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  All checkpoint spacings produce exactly the same complete finite-bit
      component laws, and the uncheckpointed law stays within its independent
      outward ideal target interval budget.
  P2  Extra checkpointing lowers the reported persistent forward/branch
      storage bound but increases recomputation work; no time-speedup claim is
      made.
  P3  The b=3 no-background zero-prefix fallback preserves normalization, and
      a forced low-precision replay after a selected bit preserves the exact
      law while charging replay.
  C1  Resetting the backward effect at the same measured prefix changes the
      next conditional law and is therefore not an admissible checkpoint.

The wide row is a supplied r=3*(2^60-1), b=3, t=63 schedule with no routes.
It compares forced output 0 and a few matched seeded samples.  Counters are
algorithmic storage/work records, not RSS measurements.  No production or
generic propagator is added.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from random import Random

import numpy as np

from lab import Experiment
from lab.verified_finite_work import VerifiedFiniteWork
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_rejection import VerifiedRejectionSampler
from lab.sampling_error import plan_prefix_mass_accuracy
from experiments.experiment_verified_rejection import target_interval, tv_from_intervals


TARGET = Fraction(1, 1_000_000)
PERIOD = 9
BLOCK = 3
WIDTH = 4
SPACINGS = (None, 1, 2, 4)
MAX_BYTES = 32 << 20

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all checkpoint spacings preserve exact finite-bit laws and the outward target budget")
exp.predict("P2", "checkpoint spacing trades persistent storage for recomputation work")
exp.predict("P3", "zero-prefix fallback and forced precision replay preserve normalization and charge work")
exp.must_fail("C1", "resetting the backward effect at a measured prefix preserves the next conditional law")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def report_path():
    return Path("out") / f"checkpoint_comparison_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def fixture(*, backgrounds=True, width=WIDTH, period=PERIOD):
    bg = {s: ("x", Fraction(1, 7)) for s in range(1, width + 1)} if backgrounds else {}
    refs = {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))} if width <= 4 else {}
    return VerifiedReflectionCircuit(period, width, bg, refs, block_size=BLOCK)


def component_law(circuit, mask, spacing):
    worker = VerifiedFiniteWork(circuit, mask, checkpoint_spacing=spacing)
    M, Q = circuit.sectors, 1 << circuit.width
    if circuit.width > 4 or circuit.period > 10 or M*Q*4096 > MAX_BYTES:
        raise MemoryError("only capped tiny complete-law dictionaries are allowed")
    guard((M, Q), label="tiny component law")
    law = {}
    conditional_masses = []
    stats = []
    started = time.perf_counter()
    for initial in range(M):
        conditional = Fraction(0)
        final = set()
        for output in range(Q):
            result = worker.path(initial, output=output, target_tv=TARGET)
            probability = Fraction(result["conditional_path_probability"])
            conditional += probability
            final.add(result["final_coarse_sector"])
            law[(result["final_coarse_sector"], output)] = law.get(
                (result["final_coarse_sector"], output), Fraction(0)) + probability / M
            if output == 0:
                stats.append({k: result.get(k) for k in (
                    "initial_coarse_sector", "final_coarse_sector", "forward_steps",
                    "recomputed_forward_steps", "branch_matrix_pair_constructions",
                    "peak_retained_forward_matrix_count", "stored_branch_matrix_count",
                    "stored_matrix_count_upper_bound", "max_working_precision",
                    "refinement_retries", "replayed_effect_steps",
                    "zero_approximate_block_fallbacks")})
        if conditional != 1 or len(final) != 1:
            raise AssertionError("checkpoint component conditional law lost mass or route determinism")
        conditional_masses.append(str(conditional))
    if sum(law.values()) != 1:
        raise AssertionError("checkpoint component law did not normalize")
    elapsed = time.perf_counter() - started
    return law, dict(conditional_masses=conditional_masses, stats=stats,
                     elapsed_seconds=elapsed)


def target_bound(circuit, mask, law):
    route = VerifiedRejectionSampler(circuit).component(mask)
    intervals = {(g, y): target_interval(route, g, y, 256)
                 for g in range(circuit.sectors) for y in range(1 << circuit.width)}
    plan = plan_prefix_mass_accuracy(circuit.width, TARGET)
    return dict(tv_upper=str(tv_from_intervals(law, intervals)),
                plan=json_safe(plan),
                ideal_mass_contains_one=sum(lo for lo,_ in intervals.values()) <= 1 <= sum(hi for _,hi in intervals.values()),
                tv_within_plan=tv_from_intervals(law, intervals) <= plan["total_tv_upper_bound"] <= TARGET)


def zero_prefix_control():
    zero = fixture(backgrounds=False)
    rows = []
    laws = {}
    for spacing in SPACINGS:
        law, stats = component_law(zero, 0, spacing)
        laws[str(spacing)] = law
        rows.append(dict(spacing=spacing, **stats))
    # The t=4 null checks complete-law normalization.  Its branches are not
    # yet small enough to round to an all-zero dyadic block, so use a bounded
    # no-background t=32 forced rare prefix to exercise the actual fallback.
    deep = fixture(backgrounds=False, width=32)
    deep_rows = []
    for spacing in (None, 4, 8):
        result = VerifiedFiniteWork(deep, 0, checkpoint_spacing=spacing).path(
            0, output=1, target_tv=TARGET)
        deep_rows.append(dict(spacing=spacing,
                              probability=str(result["conditional_path_probability"]),
                              zero_approximate_block_fallbacks=result["zero_approximate_block_fallbacks"],
                              forward_steps=result["forward_steps"],
                              recomputed_forward_steps=result["recomputed_forward_steps"]))
    return dict(rows=rows, laws_equal=all(laws[str(s)] == laws["None"] for s in SPACINGS),
                fallback_observed=any(
                    item.get("zero_approximate_block_fallbacks", 0) > 0
                    for row in rows for item in row["stats"]),
                masses=[row["conditional_masses"] for row in rows],
                deep_rows=deep_rows,
                deep_fallback_observed=any(row["zero_approximate_block_fallbacks"] > 0
                                           for row in deep_rows),
                deep_probabilities_equal=len({row["probability"] for row in deep_rows}) == 1)


def replay_and_reset_control(circuit):
    worker = VerifiedFiniteWork(circuit, 3, checkpoint_spacing=2)
    # Keep one fixed requested accuracy. Widen a valid E enclosure only
    # AFTER a bit is chosen; the actual production rebuild must replay it.
    from flint import ctx,arb
    cursor = worker.cursor(1, accuracy_bits=40, initial_precision=64)
    first = cursor.weights()
    cursor.advance(0)
    with ctx.workprec(64):
        for i in range(3):
            cursor.effect[i,i] += arb(0,1)
    correct = cursor.weights()
    replayed = cursor.replayed
    reference = worker.cursor(1, accuracy_bits=40, initial_precision=cursor.precision)
    reference.weights()
    reference.advance(0)
    reference_next = reference.weights()
    reset = worker.cursor(1, accuracy_bits=40, initial_precision=cursor.precision)
    reset.weights()
    reset.advance(0)
    with ctx.workprec(cursor.precision):
        reset.effect = reset.owner._identity()
    wrong = reset.weights()
    conditional_tv = sum(abs(Fraction(a,sum(correct))-Fraction(b,sum(wrong)))
                         for a,b in zip(correct,wrong))/2
    return dict(first_weights=first, correct_next_weights=correct,
                reference_next_weights=reference_next,
                forced_replay_matches_reference=correct == reference_next,
                reset_next_weights=wrong, reset_conditional_tv=str(conditional_tv),
                reset_changes=conditional_tv > Fraction(1,10000),
                forced_replayed_effect_steps=replayed,
                forced_forward_builds=cursor.builds,forced_precision=cursor.precision)


def wide_compare():
    width = 63
    period = 3 * ((1 << 60) - 1)
    bg = {s: ("x", Fraction(1, 7)) for s in range(1, width + 1)}
    circuit = VerifiedReflectionCircuit(period, width, bg, {}, block_size=BLOCK)
    output_rows = []
    for spacing in (None, 8):
        started = time.perf_counter()
        worker = VerifiedFiniteWork(circuit, 0, checkpoint_spacing=spacing)
        forced = worker.path(0, output=0, target_tv=TARGET)
        elapsed = time.perf_counter() - started
        output_rows.append(dict(spacing=spacing, forced_probability=str(forced["conditional_path_probability"]),
                                elapsed_seconds=elapsed,
                                counters={k: forced.get(k) for k in (
                                    "forward_steps", "recomputed_forward_steps",
                                    "branch_matrix_pair_constructions", "peak_retained_forward_matrix_count",
                                    "stored_branch_matrix_count", "stored_matrix_count_upper_bound",
                                    "max_working_precision", "refinement_retries","forward_builds",
                                    "total_tv_upper_bound")}))
    samples = []
    for seed in (624, 625, 626):
        a = VerifiedFiniteWork(circuit, 0, checkpoint_spacing=None).sample(Random(seed), target_tv=TARGET)
        b = VerifiedFiniteWork(circuit, 0, checkpoint_spacing=8).sample(Random(seed), target_tv=TARGET)
        samples.append(dict(seed=seed, same_output=a["output"] == b["output"],
                           same_final_sector=a["final_coarse_sector"] == b["final_coarse_sector"],
                           full_output=a["output"], checkpoint_output=b["output"],
                           full_final_sector=a["final_coarse_sector"], checkpoint_final_sector=b["final_coarse_sector"]))
    return dict(period=period, width=width, forced=output_rows, samples=samples,
                samples_match=all(row["same_output"] and row["same_final_sector"] for row in samples))


def main():
    started = time.time()
    report = {"status": "PASS", "rows": [], "controls": {}}
    try:
        if (len(SPACINGS)*4+6)*3*(1 << WIDTH)*4096*3 > MAX_BYTES:
            raise MemoryError("all retained law/report/serialization allowances exceed cap")
        circuit = fixture()
        for mask in range(4):
            baseline, baseline_stats = component_law(circuit, mask, None)
            target = target_bound(circuit, mask, baseline)
            for spacing in SPACINGS:
                law, stats = (baseline, baseline_stats) if spacing is None else component_law(circuit, mask, spacing)
                report["rows"].append(dict(mask=mask, spacing=spacing, law=law,
                                           stats=stats, target=target,
                                           exact_equal_to_none=law == baseline))
        report["controls"]["zero_prefix"] = zero_prefix_control()
        report["controls"]["replay_reset"] = replay_and_reset_control(circuit)
        report["controls"]["wide"] = wide_compare()
        p1 = all(row["exact_equal_to_none"] and row["target"]["tv_within_plan"]
                 and row["target"]["ideal_mass_contains_one"]
                 and all(Fraction(mass) == 1 for mass in row["stats"]["conditional_masses"])
                 for row in report["rows"])
        full,compact = report["controls"]["wide"]["forced"]
        f,q = full["counters"],compact["counters"]
        p2 = (full["forced_probability"] == compact["forced_probability"]
              and q["stored_matrix_count_upper_bound"] < f["stored_matrix_count_upper_bound"]
              and q["peak_retained_forward_matrix_count"] < f["peak_retained_forward_matrix_count"]
              and q["stored_branch_matrix_count"] == 0 < f["stored_branch_matrix_count"]
              and q["recomputed_forward_steps"] > 0 == f["recomputed_forward_steps"]
              and all(s["forward_steps"] <= 2*63*s["forward_builds"]
                      and s["total_tv_upper_bound"] <= TARGET for s in (f,q)))
        p2 &= (report["controls"]["wide"]["samples_match"]
               and report["controls"]["zero_prefix"]["laws_equal"])
        p3 = (report["controls"]["zero_prefix"]["deep_fallback_observed"]
              and report["controls"]["zero_prefix"]["deep_probabilities_equal"]
              and report["controls"]["replay_reset"]["reset_changes"]
              and report["controls"]["replay_reset"]["forced_replay_matches_reference"]
              and report["controls"]["replay_reset"]["forced_replayed_effect_steps"] > 0)
        exp.check("P1", p1, "all masks/spacings exact-equal and within outward target plan")
        exp.check("P2", p2, "checkpoint counters and wide matched seeded outputs recorded")
        exp.check("P3", p3, "zero fallback and forced precision replay retained")
        exp.fail_check("C1", report["controls"]["replay_reset"]["reset_changes"],
                       f"correct={report['controls']['replay_reset']['correct_next_weights']}, reset={report['controls']['replay_reset']['reset_next_weights']}")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time() - started
    path = report_path()
    if report["status"] == "FAIL":
        exp.check("P1",False,str(report.get("error")))
    ok = exp.finish(report_path=path,rows=json_safe(report["rows"]),
                    metadata=json_safe({k:v for k,v in report.items() if k != "rows"}))
    report["status"] = "PASS" if ok else "FAIL"
    print(json.dumps({"status": report["status"], "report": str(path),
                      "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
