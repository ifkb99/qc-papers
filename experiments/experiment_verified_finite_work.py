"""Complete-law audit of the verified finite-work component contraction.

The fixed exact-input fixture is r=10, b=2, t=4.  For each of the four
coherent history masks, every five initial coarse sectors and sixteen forced
QFT outputs are evaluated through ``VerifiedFiniteWork.path``.  The returned
finite-bit probabilities are exact Fractions and are assembled by the
deterministic final-sector route, rather than by a sampled histogram.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every component conditional law and complete joint law normalizes
      exactly; its outward ideal-TV bound is no larger than the declared
      finite-work plan and that plan is no larger than the requested target.
  P2  The independent full-r deterministic K(pi)-route products agree with
      the midpoint of the verified component enclosures; this is a diagnostic,
      not an interval proof of the float calculation.
  C1  Reversing a nonzero deterministic route q changes its complete joint
      law, so a route-sign shortcut is a genuine must-fail control.

There is no new generic propagator and no production change.  Exact zero
probabilities are retained by the Fraction kernel; no amplitude cutoff or
discarded branch is used.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.verified_finite_work import VerifiedFiniteWork, VerifiedScalarWork
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_verified_sampling import fixtures
from experiments.experiment_verified_rejection import target_interval, tv_from_intervals
from experiments.experiment_coherent_route_sampling import direct_joint, rx


TARGETS = (Fraction(1, 1000), Fraction(1, 1_000_000))
MAX_BYTES = 32 << 20

exp = Experiment("verified_finite_work", doc=__doc__, exit_on_fail=False)
exp.predict("P1", "density/effect and scalar component laws normalize and meet their declared outward TV budgets")
exp.predict("P2", "both verified component implementations agree with independent full-r deterministic products")
exp.must_fail("C1", "reversing a nonzero deterministic route changes the component joint law")


def guard(shape, dtype=np.float64, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def guard_entries(entries, bytes_per_entry=4096, label="retained entries"):
    if type(entries) is not int or entries < 0:
        raise ValueError("entry count must be a nonnegative integer")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError(f"{label} allowance exceeds 32 MiB")


def report_path():
    return Path("out") / f"verified_finite_work_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def add(out, key, value):
    out[key] = out.get(key, Fraction(0)) + Fraction(value)


def deterministic_float_component(comparator, verified, mask):
    """Independent full-r reference for one selected K(pi) history route."""
    M = verified.sectors
    selected = {}
    for i, (s, (q, _theta)) in enumerate(verified.reflections.items()):
        if mask & (1 << i):
            selected[s] = (q, np.pi)
    # Reuse exactly the comparator's supplied background block matrices.  Only
    # the deterministic route insertion is changed; no new state propagator
    # is introduced here.
    return direct_joint(CoherentReflectionCircuit(
        verified.period, verified.b, verified.width,
        comparator.background.defects, selected))


def component_finite_law(verified, mask, target, worker_type):
    M, Q = verified.sectors, 1 << verified.width
    guard((M, Q), label="finite-work joint law")
    guard_entries(M * Q + M + 8, label="finite-work retained law entries")
    worker = worker_type(verified, mask)
    law = {}
    conditional_rows = []
    path_stats = []
    for initial_sector in range(M):
        conditional_mass = Fraction(0)
        final_sectors = set()
        for output in range(Q):
            result = worker.path(initial_sector, output=output, target_tv=target)
            probability = Fraction(result["conditional_path_probability"])
            if probability < 0:
                raise AssertionError("finite-work path probability became negative")
            conditional_mass += probability
            final_sectors.add(result["final_coarse_sector"])
            add(law, (result["final_coarse_sector"], output),
                probability / M)
            if output == 0:
                path_stats.append(dict(initial_sector=initial_sector,
                                       final_coarse_sector=result["final_coarse_sector"],
                                       zero_approximate_block_fallbacks=result["zero_approximate_block_fallbacks"],
                                       max_working_precision=result["max_working_precision"],
                                       refinement_retries=result["refinement_retries"],
                                       forward_builds=result.get("forward_builds"),
                                       backward_child_evaluations=result.get("backward_child_evaluations"),
                                       scalar_conditional_evaluations=result.get("scalar_conditional_evaluations")))
        if conditional_mass != 1:
            raise AssertionError(f"conditional law for initial sector {initial_sector} lost mass")
        if len(final_sectors) != 1:
            raise AssertionError("deterministic route changed final sector by output")
        conditional_rows.append(dict(initial_sector=initial_sector,
                                     mass=str(conditional_mass),
                                     final_coarse_sector=next(iter(final_sectors))))
    if sum(law.values()) != 1:
        raise AssertionError("finite-work joint law did not normalize")
    # The path planner is independent of the forced output and therefore its
    # declared bound is the same for every row in this component.
    plan = worker.path(0, output=0, target_tv=target)
    return law, dict(plan={key: value for key, value in plan.items()
                           if key in ("target_tv", "total_tv_upper_bound",
                                      "oracle_tv_upper_bound",
                                      "finite_random_bits_tv_upper_bound",
                                      "weight_abs_error", "categorical_random_bits",
                                      "weight_accuracy_bits")},
                     conditional_rows=conditional_rows,
                     path_stats=path_stats)


def component_intervals(component, precision):
    M, Q = component.sectors, 1 << component.width
    return {(gamma, output): target_interval(component, gamma, output, precision)
            for gamma in range(M) for output in range(Q)}


def run_target(verified, comparator, target, worker_type):
    if verified.period != 10 or verified.width != 4:
        raise ValueError("only the frozen r10,t4 fixture is allowed")
    sampler = VerifiedRejectionSampler(verified)
    rows = []
    for mask in range(1 << sampler.k):
        component = sampler.component(mask)
        law, stats = component_finite_law(verified, mask, target, worker_type)
        intervals_192 = component_intervals(component, 192)
        intervals_256 = component_intervals(component, 256)
        tv_192 = tv_from_intervals(law, intervals_192)
        tv_256 = tv_from_intervals(law, intervals_256)
        ideal_mid = np.array([[float(sum(intervals_256[(gamma, output)]) / 2)
                               for output in range(1 << verified.width)]
                              for gamma in range(verified.sectors)])
        independent = deterministic_float_component(comparator, verified, mask)
        interval_total_lo = sum(lo for lo, _hi in intervals_256.values())
        interval_total_hi = sum(hi for _lo, hi in intervals_256.values())
        rows.append(dict(mask=mask, law_mass=str(sum(law.values())),
                         plan=json_safe(stats["plan"]),
                         conditional_rows=stats["conditional_rows"],
                         path_stats=stats["path_stats"],
                         component_tv_upper_P192=str(tv_192),
                         component_tv_upper_P256=str(tv_256),
                         plan_meets_target=Fraction(stats["plan"]["total_tv_upper_bound"]) <= target,
                         tv_meets_plan_P256=tv_256 <= stats["plan"]["total_tv_upper_bound"],
                         ideal_interval_total_contains_one=interval_total_lo <= 1 <= interval_total_hi,
                         independent_float_mass=float(independent.sum()),
                         ideal_midpoint_float_max_error=float(np.max(np.abs(ideal_mid-independent))),
                         independent_float_law=json_safe(independent.tolist())))
    return rows


def wrong_route_control(verified, comparator):
    """Use q -> -q for the one nonzero-route mask as a frozen false shortcut."""
    mask = 2  # second insertion has q=1 in the frozen fixture
    selected = {}
    for i, (s, (q, _theta)) in enumerate(verified.reflections.items()):
        if mask & (1 << i):
            selected[s] = (q, np.pi)
    wrong = {s: ((-q) % verified.sectors, theta) for s, (q, theta) in selected.items()}
    correct_law = deterministic_float_component(comparator, verified, mask)
    wrong_law = direct_joint(CoherentReflectionCircuit(
        verified.period, verified.b, verified.width,
        comparator.background.defects, wrong))
    tv = float(np.abs(correct_law - wrong_law).sum() / 2)
    return dict(mask=mask, correct_q={str(s): int(q) for s, (q, _theta) in selected.items()},
                wrong_q={str(s): int(q) for s, (q, _theta) in wrong.items()},
                wrong_route_tv=tv, correct_mass=float(correct_law.sum()),
                wrong_mass=float(wrong_law.sum()))


def zero_block_edge(target):
    """Exact-zero QFT branches must use the declared fallback, not deletion."""
    edge = VerifiedReflectionCircuit(10, 4, {}, {})
    worker = VerifiedFiniteWork(edge, 0)
    total_fallbacks = 0
    conditional_masses = []
    for initial_sector in range(edge.sectors):
        mass = Fraction(0)
        for output in range(1 << edge.width):
            result = worker.path(initial_sector, output=output, target_tv=target)
            mass += Fraction(result["conditional_path_probability"])
            total_fallbacks += result["zero_approximate_block_fallbacks"]
        conditional_masses.append(mass)
    if any(mass != 1 for mass in conditional_masses):
        raise AssertionError("exact-zero fallback edge lost conditional mass")
    return dict(conditional_masses=[str(mass) for mass in conditional_masses],
                joint_mass=str(sum(conditional_masses, Fraction(0)) / edge.sectors),
                zero_approximate_block_fallbacks=total_fallbacks,
                fallback_observed=total_fallbacks > 0)


def fixtures_with_w0():
    """The same gates plus one fixed insertion-0 Rx(pi/4), no sweep."""
    from fractions import Fraction
    backgrounds = {0: ("x", Fraction(1, 4)),
                   1: ("x", Fraction(1, 7)),
                   3: ("z", Fraction(1, 5))}
    reflections = {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}
    verified = VerifiedReflectionCircuit(10, 4, backgrounds, reflections)
    comparator = CoherentReflectionCircuit(
        10, 2, 4,
        {0: rx(np.pi / 4), 1: rx(np.pi / 7),
         3: np.diag([np.exp(-1j * np.pi / 10), np.exp(1j * np.pi / 10)])},
        {2: (0, np.pi / 5), 3: (1, np.pi / 5)})
    return verified, comparator


def main():
    started = time.time()
    report = {"status": "PASS", "rows": [], "controls": {}}
    try:
        fixture_set = [("original", *fixtures()), ("insertion0_rx_pi4", *fixtures_with_w0())]
        worker_set = [("density_effect", VerifiedFiniteWork),
                      ("scalar", VerifiedScalarWork)]
        for fixture_name, verified, comparator in fixture_set:
            for worker_name, worker_type in worker_set:
                for target in TARGETS:
                    report["rows"].append(dict(fixture=fixture_name,
                                               worker=worker_name,
                                               target_tv=str(target),
                                               components=run_target(verified, comparator, target, worker_type)))
        report["controls"]["wrong_route"] = wrong_route_control(verified, comparator)
        report["controls"]["wrong_route_fails"] = report["controls"]["wrong_route"]["wrong_route_tv"] > 1e-8
        report["controls"]["exact_zero_edge"] = zero_block_edge(TARGETS[0])
        p1 = True
        for row in report["rows"]:
            target = Fraction(row["target_tv"])
            for component in row["components"]:
                plan_bound = Fraction(component["plan"]["total_tv_upper_bound"])
                p1 &= (Fraction(component["law_mass"]) == 1
                       and plan_bound <= target
                       and Fraction(component["component_tv_upper_P192"]) <= plan_bound
                       and Fraction(component["component_tv_upper_P256"]) <= plan_bound)
        p1 &= (report["controls"]["exact_zero_edge"]["fallback_observed"]
               and Fraction(report["controls"]["exact_zero_edge"]["joint_mass"]) == 1)
        p2 = all(component["ideal_interval_total_contains_one"]
                  and abs(component["independent_float_mass"] - 1) < 2e-14
                  and component["ideal_midpoint_float_max_error"] < 2e-14
                  for row in report["rows"] for component in row["components"])
        exp.check("P1", p1, "both verified component implementations and both fixtures normalized within their plans")
        exp.check("P2", p2, "independent deterministic full-r products agree with verified ideal midpoints")
        exp.fail_check("C1", report["controls"]["wrong_route_fails"],
                       f"wrong-route TV={report['controls']['wrong_route']['wrong_route_tv']:.8g}")
        path = report_path()
        if not exp.finish(report_path=path,rows=json_safe(report["rows"]),metadata=dict(
                controls=json_safe(report["controls"]),optional_backend="python-flint==0.9.0")):
            raise AssertionError("Experiment harness failed")
        print(json.dumps(dict(ok=True,report=str(path))))
        return
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time() - started
    path = report_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(report), indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": str(path),
                      "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
