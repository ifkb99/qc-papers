"""Verified tiny-law audit for the fixed observable b=3 coherent fixture.

This uses r=9, b=3, t=4, M=3 with W1=embedded Rx(pi/7),
W3=embedded Rz(pi/5), and K2(q=0,pi/5), K3(q=1,pi/5).  The generalized
shared references enumerate the complete finite-bit laws; this file adds the
per-history outward TV checks that are needed before mixing the four histories
or applying rejection acceptance.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Prefix and finite-work component laws normalize exactly and each stays
      within its own planner budget at target TV 1e-3 and 1e-6.
  P2  Both component methods, in rectangular and norm enclosure modes, give
      normalized accepted laws inside the full sampler budget and agree with
      independent full-r products at diagnostic precision.
  P3  The b=3 late-background removal control changes the joint law, while a
      final-only W4 remains invisible.
  C1  The universal claim that every late background is removable must fail on
      this fixed b=3 input.

No scalar b=2 shortcut is used; no generic propagator or amplitude cutoff is
introduced.  Dense and retained tiny-reference allocations are guarded by the
shared r<=10,t<=4 limits.
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
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_coherent_route_sampling import direct_joint, rx
from experiments.experiment_verified_rejection import (
    component_laws, run_target, target_interval, tv_from_intervals)


TARGETS = (Fraction(1, 1000), Fraction(1, 1_000_000))
MAX_BYTES = 32 << 20
M = 3
BLOCK = 3
WIDTH = 4

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all b=3 component laws normalize and meet individual finite-work/prefix budgets")
exp.predict("P2", "rectangular and norm accepted laws meet full plans and independent diagnostics")
exp.predict("P3", "late-removal and final-only b=3 boundary controls behave distinctly")
exp.must_fail("C1", "a universal all-late-background removal claim fails on the fixed b=3 input")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def report_path():
    return Path("out") / f"verified_odd_block_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


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
    return value


def verified_fixture(mode="rectangular", *, late=True, final_only=False):
    backgrounds = {}
    if late:
        backgrounds = {1: ("x", Fraction(1, 7)),
                       3: ("z", Fraction(1, 5))}
    if final_only:
        backgrounds = {4: ("x", Fraction(1, 7))}
    refs = {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}
    return VerifiedReflectionCircuit(9, WIDTH, backgrounds, refs,
                                     enclosure_mode=mode, block_size=BLOCK)


def float_fixture(*, late=True, final_only=False):
    backgrounds = {}
    if late:
        rx3 = np.eye(BLOCK, dtype=complex)
        rx3[:2, :2] = rx(np.pi / 7)
        rz3 = np.eye(BLOCK, dtype=complex)
        rz3[0, 0] = np.exp(-1j * np.pi / 10)
        rz3[1, 1] = np.exp(1j * np.pi / 10)
        backgrounds = {1: rx3, 3: rz3}
    if final_only:
        rx3 = np.eye(BLOCK, dtype=complex)
        rx3[:2, :2] = rx(np.pi / 7)
        backgrounds = {4: rx3}
    return CoherentReflectionCircuit(
        9, BLOCK, WIDTH, backgrounds,
        {2: (0, np.pi / 5), 3: (1, np.pi / 5)})


def joint_matrix(circuit):
    guard((circuit.sectors, 1 << circuit.width), label="joint diagnostic law")
    return np.array([[circuit.joint_probability(g, y)
                      for y in range(1 << circuit.width)]
                     for g in range(circuit.sectors)])


def component_budget_audit(verified, target, method):
    sampler = VerifiedRejectionSampler(verified, component_method=method)
    plan = sampler.plan(target)
    components, _thresholds = component_laws(sampler, plan)
    rows = []
    for mask, entry in components.items():
        component = sampler.component(mask)
        intervals_192 = {(g, y): target_interval(component, g, y, 192)
                          for g in range(M) for y in range(1 << WIDTH)}
        intervals_256 = {(g, y): target_interval(component, g, y, 256)
                          for g in range(M) for y in range(1 << WIDTH)}
        tv_192 = tv_from_intervals(entry["law"], intervals_192)
        tv_256 = tv_from_intervals(entry["law"], intervals_256)
        if method == "prefix":
            component_plan = component.accuracy_plan(plan["component_target_tv"])
        else:
            from lab.sampling_error import plan_prefix_mass_accuracy
            component_plan = plan_prefix_mass_accuracy(WIDTH, plan["component_target_tv"])
        rows.append(dict(mask=mask, mass=str(sum(entry["law"].values())),
                         tv_upper_P192=str(tv_192), tv_upper_P256=str(tv_256),
                         plan=json_safe(component_plan), stats=json_safe(entry["stats"]),
                         tv_meets_plan=(tv_192 <= component_plan["total_tv_upper_bound"]
                                        and tv_256 <= component_plan["total_tv_upper_bound"]),
                         plan_meets_component_target=(component_plan["total_tv_upper_bound"]
                                                      <= plan["component_target_tv"])))
    return rows


def boundary_controls():
    late = direct_joint(float_fixture(late=True))
    removed = direct_joint(float_fixture(late=False))
    no_background = direct_joint(float_fixture(late=False))
    final_only = direct_joint(float_fixture(late=False, final_only=True))
    late_delta = np.abs(late - removed)
    final_delta = np.abs(no_background - final_only)
    return dict(late_mass=float(late.sum()), removed_mass=float(removed.sum()),
                late_removal_max_abs_error=float(late_delta.max()),
                late_removal_tv=float(late_delta.sum() / 2),
                no_background_mass=float(no_background.sum()),
                final_only_mass=float(final_only.sum()),
                final_only_max_abs_error=float(final_delta.max()),
                final_only_tv=float(final_delta.sum() / 2),
                late_removal_changes=float(late_delta.sum() / 2) > 1e-10,
                final_only_invisible=float(final_delta.max()) < 2e-14)


def main():
    started = time.time()
    report = {"status": "PASS", "rows": [], "controls": {}}
    try:
        # Conservative allowance for all retained row dictionaries, in
        # addition to the shared evaluator's transient frontier/cache guards.
        if len(TARGETS)*4*(4+16)*M*(1 << WIDTH)*4096 > MAX_BYTES:
            raise MemoryError("retained multi-row report allowance exceeded")
        # Every mode/method/target is a separate fixed input row; only target
        # TV changes within each repeated fixture.
        for mode in ("rectangular", "norm"):
            for method in ("prefix", "finite_work"):
                verified = verified_fixture(mode)
                comparator = float_fixture(late=True)
                for target in TARGETS:
                    row = run_target(verified, comparator, target,
                                     component_method=method)
                    row["enclosure_mode"] = mode
                    row["component_budget_audit"] = component_budget_audit(
                        verified, target, method)
                    report["rows"].append(row)
        report["controls"] = boundary_controls()
        p1 = True
        p2 = True
        for row in report["rows"]:
            target = Fraction(row["target_tv"])
            plan = row["plan"]
            full_plan = Fraction(plan["total_tv_upper_bound"])
            proposal_plan = Fraction(plan["proposal_tv_upper_bound"])
            p1 &= (Fraction(row["proposal_mass"]) == 1
                   and Fraction(row["accepted_mass"]) == 1
                   and Fraction(row["target_tv_upper_P192"]) <= full_plan
                   and Fraction(row["target_tv_upper_P256"]) <= full_plan
                   and Fraction(row["proposal_tv_upper_from_component_Arb"]) <= proposal_plan
                   and row["exact_success_meets_plan_lower"]
                   and row["ideal_target_interval_total_contains_one"]
                   and abs(row["independent_target_mass"] - 1) < 2e-14)
            p2 &= (full_plan <= target
                   and row["ideal_target_independent_max_error"] < 2e-14
                   and abs(row["independent_target_mass"] - 1) < 2e-14
                   and all(d["ball_total_contains_one"]
                           and abs(d["full_r_mass"] - 1) < 2e-14
                           and d["independent_max_error"] < 2e-14
                           for d in row["component_independent_diagnostics"])
                   and all(Fraction(c["mass"]) == 1
                           and c["tv_meets_plan"]
                           and c["plan_meets_component_target"]
                           for c in row["component_budget_audit"]))
        p3 = (report["controls"]["late_removal_changes"]
              and report["controls"]["final_only_invisible"]
              and abs(report["controls"]["late_mass"] - 1) < 1e-13
              and abs(report["controls"]["removed_mass"] - 1) < 1e-13
              and abs(report["controls"]["no_background_mass"] - 1) < 1e-13
              and abs(report["controls"]["final_only_mass"] - 1) < 1e-13)
        exp.check("P1", p1, "all b=3 accepted/proposal/component laws normalized within plans")
        exp.check("P2", p2, "rectangular/norm and prefix/finite-work diagnostics agree independently")
        exp.check("P3", p3, "late-removal visibility and final-only invariance checked")
        exp.fail_check("C1", report["controls"]["late_removal_changes"],
                       f"b3 late-removal TV={report['controls']['late_removal_tv']:.8g}")
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]),
                          metadata=dict(controls=json_safe(report["controls"]),
                                        elapsed_seconds=time.time()-started,
                                        optional_backend="python-flint==0.9.0",
                                        native_memory_measured=False)):
            raise AssertionError("Experiment harness failed")
        print(json.dumps(dict(ok=True, report=str(path))))
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
