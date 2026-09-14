"""Finite-TV audit of the exact-input verified coherent-component rejection law.

The frozen fixture is the r=10,t=4 exact rational-pi family used by
``experiment_verified_sampling``.  For target TV values 1e-3 and 1e-6, all
four history component laws are enumerated with the existing finite-bit
transition reference, mixed with exact Fraction history probabilities, and
then passed through every actual dyadic acceptance threshold for all 80
(final-sector, output) cells.  Arb target intervals are converted to exact
Fractions directly from binary endpoints; no float-to-Fraction path is used.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The finite proposal and accepted submeasure laws normalize exactly, and
      accepted success and target TV remain within the planned diagnostic
      budget at both target settings.
  P2  Outward Arb proposal TV is within the separate proposal budget, success
      exceeds its certified lower bound, and IDEAL Arb target/component
      midpoints agree with independent full-r float matrix products.

  C1  Skipping acceptance leaves a proposal law measurably different from the
      target despite accurate component proposals.
  C2  Acceptance-error accuracy alone cannot certify a deliberately biased
      proposal (a two-cell exact counterexample).

No production helper is modified and no wide timing claim is made.
The first agent report had an invalid upper-bound-as-lower-bound control and
a vacuous nonnegative-proposal-TV predicate; main strengthened both before
interpreting it. Original reports/logs are retained; see the research note.
Tiny enumeration uses the existing r10,t4/cache guards. Retained dictionaries
are preflighted with a conservative 4096-byte allowance per entry at these
two fixed precisions; this is not a measurement of native/RSS memory.
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
from lab.verified_prefix import binary_fraction
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_verified_sampling import exact_transition_law, fixtures
from experiments.experiment_coherent_route_sampling import direct_joint


TARGETS = (Fraction(1, 1000), Fraction(1, 1_000_000))
MAX_BYTES = 32 << 20

exp = Experiment("verified_rejection", doc=__doc__, exit_on_fail=False)
exp.predict("P1", "exact finite-bit proposal/acceptance laws normalize and meet planned diagnostics")
exp.predict("P2", "proposal TV/success meet budgets and ideal intervals match independent products")
exp.must_fail("C1", "skipping acceptance preserves proposal accuracy but changes the target law")
exp.must_fail("C2", "acceptance error alone cannot certify a biased proposal")


def guard(shape, dtype=np.float64, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 32 MiB")


def report_path():
    return Path("out") / f"verified_rejection_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key,item in value.items()}
    if isinstance(value, (list,tuple)):
        return [json_safe(item) for item in value]
    return value


def add(out, key, value):
    out[key] = out.get(key, Fraction(0)) + Fraction(value)


def target_interval(circuit, sector, output, precision):
    """Outward exact Fraction interval for p(sector,output)."""
    from flint import ctx
    with ctx.workprec(precision):
        vector = circuit.prefix_enclosure(sector, 0, circuit.width,
                                          measured=circuit.width, output=output,
                                          working_precision=precision)
        total = sum((z.real*z.real + z.imag*z.imag for z in vector)) / circuit.sectors
        lo = binary_fraction(total.lower())
        hi = binary_fraction(total.upper())
    return max(Fraction(0), lo), max(Fraction(0), hi)


def true_history_weight_bounds(sampler, precision):
    from flint import arb, ctx, fmpq
    c = sampler.circuit
    with ctx.workprec(precision):
        factors = []
        B_lo, B_hi = Fraction(1), Fraction(1)
        for _, theta in c.reflections.values():
            sine, cosine = arb(fmpq(theta.numerator, 2*theta.denominator)).sin_cos_pi()
            ss = abs(sine); cc = abs(cosine)
            slo, shi = binary_fraction(ss.lower()), binary_fraction(ss.upper())
            clo, chi = binary_fraction(cc.lower()), binary_fraction(cc.upper())
            denom_lo, denom_hi = slo + clo, shi + chi
            factors.append((slo, shi, clo, chi, denom_lo, denom_hi))
            B_lo *= denom_lo
            B_hi *= denom_hi
        bounds = []
        for mask in range(1 << sampler.k):
            lo = hi = Fraction(1)
            for i, (slo, shi, clo, chi, dlo, dhi) in enumerate(factors):
                if mask & (1 << i):
                    lo *= slo / dhi if dhi else 0
                    hi *= shi / dlo if dlo else 1
                else:
                    lo *= clo / dhi if dhi else 0
                    hi *= chi / dlo if dlo else 1
            bounds.append((lo, hi))
        return bounds, (B_lo, B_hi)


def component_laws(sampler, plan):
    thresholds = sampler.history_thresholds(plan["history_random_bits"])
    Lh = thresholds["random_bits"]
    masks = {}
    for mask in range(1 << sampler.k):
        weight = Fraction(1)
        for i, threshold in enumerate(thresholds["thresholds"]):
            hit = bool(mask & (1 << i))
            count = threshold if hit else (1 << Lh)-threshold
            weight *= Fraction(count, 1 << Lh)
        marginal = {}
        if sampler.component_method == "prefix":
            component = sampler.component(mask)
            component_plan = component.accuracy_plan(plan["component_target_tv"])
            law, stats = exact_transition_law(
                component, component_plan["absolute_accuracy_bits"],
                component_plan["random_bits"])
            for (sector, _p, output), mass in law.items():
                add(marginal, (sector, output), mass)
        else:
            from lab.verified_finite_work import VerifiedFiniteWork,VerifiedScalarWork
            from lab.sampling_error import plan_prefix_mass_accuracy,plan_conditional_weight_accuracy
            factory = VerifiedFiniteWork if sampler.component_method == "finite_work" else VerifiedScalarWork
            planner = plan_prefix_mass_accuracy if sampler.component_method == "finite_work" else plan_conditional_weight_accuracy
            component = factory(sampler.circuit,mask)
            component_plan = planner(sampler.circuit.width,plan["component_target_tv"])
            stats = dict(max_working_precision=0,refinement_retries=0,forced_paths=0)
            for initial in range(sampler.circuit.sectors):
                for output in range(1 << sampler.circuit.width):
                    forced = component.path(initial,output=output,target_tv=plan["component_target_tv"])
                    add(marginal,(forced["final_coarse_sector"],output),
                        forced["conditional_path_probability"]/sampler.circuit.sectors)
                    stats["max_working_precision"] = max(stats["max_working_precision"],forced["max_working_precision"])
                    stats["refinement_retries"] += forced["refinement_retries"]
                    stats["forced_paths"] += 1
        if sum(marginal.values()) != 1:
            raise AssertionError("component law did not normalize")
        masks[mask] = dict(weight=weight, law=marginal,
                           stats=stats, accuracy_plan=component_plan)
    if sum(row["weight"] for row in masks.values()) != 1:
        raise AssertionError("history proposal weights did not normalize")
    return masks, thresholds


def tv_from_intervals(law, intervals):
    total = Fraction(0)
    # Include cells with zero exact mass as well.  The target has all 80
    # finite cells, and omitting a zero proposal/accepted cell would hide a
    # support error in the TV audit.
    for key in set(law) | set(intervals):
        mass = law.get(key, Fraction(0))
        lo, hi = intervals[key]
        total += max(abs(mass-lo), abs(mass-hi))
    return total / 2


def tv_lower_from_intervals(law, intervals):
    return sum(max(Fraction(0), lo-law.get(key,Fraction(0)),
                   law.get(key,Fraction(0))-hi)
               for key,(lo,hi) in intervals.items())/2


def run_target(verified, comparator, target, *, component_method="prefix"):
    if (verified.period > 10 or verified.width > 4 or verified.b not in (2, 3)
            or target not in TARGETS
            or (verified.period,verified.b,verified.width) !=
               (comparator.period,comparator.b,comparator.width)):
        raise ValueError("matched tiny r<=10,t<=4,b=2/3 reference budget only")
    # Existing transition evaluator guards its sequential query cache.
    # Here count retained marginals, interval/decision tables and row data.
    retained_entries = ((1 << len(verified.reflections)) + 16) * verified.sectors * (1 << verified.width)
    if retained_entries * 4096 > MAX_BYTES:
        raise MemoryError("retained Python dictionary allowance exceeded")
    sampler = VerifiedRejectionSampler(verified,component_method=component_method)
    plan = sampler.plan(target)
    components, thresholds = component_laws(sampler, plan)
    M, Q = verified.sectors, 1 << verified.width
    guard((M, Q), label="tiny joint law")
    proposal = {}
    for row in components.values():
        for key, mass in row["law"].items():
            add(proposal, key, row["weight"]*mass)
    if sum(proposal.values()) != 1:
        raise AssertionError("proposal law did not normalize")
    accepted_submeasure = {}
    acceptance_cells = {}
    for sector in range(M):
        for output in range(Q):
            decision = sampler.acceptance(
                sector, output,
                interval_tolerance=plan["acceptance_interval_width_tolerance"],
                random_bits=plan["acceptance_random_bits"])
            acceptance = decision.probability
            acceptance_cells[(sector, output)] = dict(probability=acceptance,
                numerator_bounds=tuple(str(x) for x in decision.numerator_bounds),
                denominator_bounds=tuple(str(x) for x in decision.denominator_bounds),
                width_sum=str(decision.width_sum), refinements=decision.refinements)
            add(accepted_submeasure, (sector, output),
                proposal.get((sector, output), Fraction(0))*acceptance)
    success = sum(accepted_submeasure.values())
    if success <= 0:
        raise AssertionError("verified rejection success mass vanished")
    accepted = {key: value/success for key, value in accepted_submeasure.items()}
    if sum(accepted.values()) != 1:
        raise AssertionError("accepted law did not normalize")
    target_intervals_192 = {(g,y): target_interval(verified,g,y,192)
                            for g in range(M) for y in range(Q)}
    target_intervals_256 = {(g,y): target_interval(verified,g,y,256)
                            for g in range(M) for y in range(Q)}
    target_tv_192 = tv_from_intervals(accepted, target_intervals_192)
    target_tv_256 = tv_from_intervals(accepted, target_intervals_256)
    proposal_target_tv = tv_from_intervals(proposal, target_intervals_256)
    proposal_target_tv_lower = tv_lower_from_intervals(proposal, target_intervals_256)
    ideal_weight_bounds, B_bounds = true_history_weight_bounds(sampler, 192)
    proposal_intervals = {}
    for sector in range(M):
        for output in range(Q):
            lo = hi = Fraction(0)
            for mask, row in components.items():
                comp_lo, comp_hi = target_interval(sampler.component(mask), sector, output, 192)
                wlo, whi = ideal_weight_bounds[mask]
                lo += wlo*comp_lo
                hi += whi*comp_hi
            proposal_intervals[(sector, output)] = (lo, hi)
    proposal_tv_upper = tv_from_intervals(proposal, proposal_intervals)
    independent = direct_joint(comparator)
    ideal_target_mid = np.array([[float(sum(target_intervals_256[(g,y)])/2)
                                 for y in range(Q)] for g in range(M)])
    component_diagnostics = []
    for mask in components:
        route = sampler.component(mask)
        independent_component = direct_joint(CoherentReflectionCircuit(verified.period,verified.b,verified.width,
            comparator.background.defects,
            {s:(q,np.pi) for s,(q,_) in route.reflections.items()}))
        ideal_intervals = {(g,y):target_interval(route,g,y,192) for g in range(M) for y in range(Q)}
        ideal_mid = np.array([[float(sum(ideal_intervals[(g,y)])/2)
                              for y in range(Q)] for g in range(M)])
        component_diagnostics.append(dict(mask=mask,
            full_r_mass=float(independent_component.sum()),
            ball_total_contains_one=(sum(lo for lo,hi in ideal_intervals.values()) <= 1
                                    <= sum(hi for lo,hi in ideal_intervals.values())),
            independent_max_error=float(np.max(np.abs(ideal_mid-independent_component)))))
    accepted_float = np.zeros_like(independent)
    proposal_float = np.zeros_like(independent)
    for (g,y), mass in accepted.items():
        accepted_float[g,y] = float(mass)
    for (g,y), mass in proposal.items():
        proposal_float[g,y] = float(mass)
    return dict(target_tv=str(target), component_method=component_method,
                plan={k:str(v) if isinstance(v,Fraction) else v for k,v in plan.items()},
                history_thresholds=thresholds,
                history_weights={str(k):str(v["weight"]) for k,v in components.items()},
                component_stats={str(k):v["stats"] for k,v in components.items()},
                proposal_mass=str(sum(proposal.values())), accepted_mass=str(sum(accepted.values())),
                success_probability=str(success),
                exact_success_meets_plan_lower=success >= plan["success_probability_lower_bound"],
                target_tv_upper_P192=str(target_tv_192), target_tv_upper_P256=str(target_tv_256),
                proposal_vs_target_tv_upper_P256=str(proposal_target_tv),
                proposal_vs_target_tv_lower_P256=str(proposal_target_tv_lower),
                ideal_target_independent_max_error=float(np.max(np.abs(ideal_target_mid-independent))),
                independent_target_mass=float(independent.sum()),
                ideal_target_interval_total_contains_one=(sum(lo for lo,hi in target_intervals_256.values())
                    <= 1 <= sum(hi for lo,hi in target_intervals_256.values())),
                component_independent_diagnostics=component_diagnostics,
                independent_float_target_max_error=float(np.max(np.abs(accepted_float-independent))),
                independent_float_proposal_max_error=float(np.max(np.abs(proposal_float-independent))),
                true_B_bounds=tuple(str(x) for x in B_bounds),
                true_B2_upper=str(B_bounds[1]*B_bounds[1]),
                proposal_tv_upper_from_component_Arb=str(proposal_tv_upper),
                # JSON object keys cannot be tuple cells; retain the exact
                # cell coordinates textually without changing any arithmetic.
                acceptance_cells={f"{g},{y}": value
                                  for (g, y), value in acceptance_cells.items()})


def main():
    started = time.time()
    report = {"status":"PASS", "rows":[], "controls":{}}
    try:
        verified, comparator = fixtures()
        for target in TARGETS:
            report["rows"].append(run_target(verified, comparator, target))
        first = report["rows"][0]
        proposal_tv = Fraction(first["proposal_vs_target_tv_lower_P256"])
        # Exact two-cell biased-proposal control: zero acceptance numerical
        # error does not imply a proposal law has the required support/weights.
        biased = (Fraction(9,10), Fraction(1,10))
        ideal = (Fraction(1,2), Fraction(1,2))
        biased_tv = sum(abs(a-b) for a,b in zip(biased, ideal))/2
        report["controls"] = dict(skip_acceptance_tv_lower=str(proposal_tv),
                                   skip_acceptance_fails=proposal_tv > Fraction(1,100),
                                   biased_proposal_tv=str(biased_tv),
                                   biased_acceptance_error=Fraction(0),
                                   biased_acceptance_only_fails=biased_tv > Fraction(1,100))
        p1 = all(Fraction(row["proposal_mass"]) == 1
                 and Fraction(row["accepted_mass"]) == 1
                 and Fraction(row["target_tv_upper_P256"]) <= Fraction(row["plan"]["total_tv_upper_bound"])
                 <= Fraction(row["target_tv"])
                 and Fraction(row["target_tv_upper_P192"]) <= Fraction(row["plan"]["total_tv_upper_bound"])
                 and row["exact_success_meets_plan_lower"]
                 for row in report["rows"])
        # The independent full-r product is deliberately a float diagnostic;
        # its cellwise discrepancy includes the finite-bit verified target
        # error.  Require a finite, useful agreement at the scale of the
        # requested target rather than pretending it is an outward proof.
        p2 = all(0 <= Fraction(row["proposal_tv_upper_from_component_Arb"])
                 <= Fraction(row["plan"]["proposal_tv_upper_bound"])
                 and row["ideal_target_independent_max_error"] < 2e-14
                 and abs(row["independent_target_mass"]-1) < 2e-14
                 and row["ideal_target_interval_total_contains_one"]
                 and all(d["independent_max_error"] < 2e-14
                         and abs(d["full_r_mass"]-1) < 2e-14 and d["ball_total_contains_one"]
                         for d in row["component_independent_diagnostics"])
                 for row in report["rows"])
        exp.check("P1", p1, f"rows={len(report['rows'])}, exact laws normalized and planned bounds checked")
        exp.check("P2", p2, "Arb outward proposal/target diagnostics and independent float comparison recorded")
        exp.fail_check("C1", report["controls"]["skip_acceptance_fails"],
                       f"proposal-vs-target TV lower={float(proposal_tv):.8g}")
        exp.fail_check("C2", report["controls"]["biased_acceptance_only_fails"],
                       f"biased proposal TV={biased_tv} with acceptance error=0")
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(report["rows"]), metadata=dict(
                controls=json_safe(report["controls"]), elapsed_seconds=time.time()-started,
                optional_backend="python-flint==0.9.0", native_memory_measured=False)):
            raise AssertionError("Experiment harness failed")
        print(json.dumps(dict(ok=True, report=str(path))))
        return
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type":type(exc).__name__, "message":str(exc),
                            "traceback":traceback.format_exc()}
    report["elapsed_seconds"] = time.time()-started
    path = report_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True,
                               default=lambda value: str(value)
                               if isinstance(value, Fraction) else repr(value)),
                    encoding="utf-8")
    print(json.dumps({"status":report["status"], "report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
