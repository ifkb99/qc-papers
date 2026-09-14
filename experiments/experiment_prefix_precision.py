"""Bounded numerical audit of the coherent-route prefix sampler.

This experiment keeps the production oracle untouched.  It rounds each
queried complex prefix amplitude to a fixed binary grid, normalizes every
conditional law, and explicitly uses a uniform-within-block law when a whole
rounded block is zero.  The exact tiny transition law and final joint law are
the references; no generic propagator is introduced.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  For every precision, the exact-parent weighted conditional TV sum gives
      final TV <= sum_s weighted_TV_s <= 2 sum_s FULL_Born_TV_s.
  P2  Every rounded law normalizes; the finest tested precision improves over
      the coarsest row.  Monotonicity of round-to-grid error is not assumed.
  P3  A normalized rounded final-amplitude comparator is reported separately;
      this prefix audit does not certify production rejection acceptance.

  C1  A normalized 64-coordinate example refutes TV<=maximum amplitude error
      without the dimension factor (an independent analytic control).
  C2  Dropping rounded bins without a declared normalized fallback is not a
      probability law and cannot be treated as a certificate.

All arrays are capped for r=10,b=2,t=4. The unmodified production default is
complex128, not exact arithmetic. This is output-grid quantization of its
amplitudes, not an internal mantissa/trigonometric precision sweep. Earlier
failed hypotheses and verifier bugs are retained in timestamped reports.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from experiments.experiment_coherent_route_sampling import fixture, direct_joint


MAX_BYTES = 16 << 20
CAP_R, CAP_T = 14, 4
PRECISIONS = (0, 1, 2, 4, 6, 8, 10, 12)

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "the weighted conditional TV chain budget bounds every rounded final law")
exp.predict("P2", "precision improves the law while exact-zero and explicit fallback cases normalize")
exp.predict("P3", "normalized rounded-final amplitude comparator obeys its global L2 bound, not a rejection-implementation claim")
exp.must_fail("C1", "pointwise amplitude error alone omits global dimension")
exp.must_fail("C2", "discarding zeroed bins without normalization is not a probability law")


def guard(shape, dtype=np.float64, label="array"):
    size = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if size > MAX_BYTES:
        raise MemoryError(f"{label} allocation {size} exceeds 16 MiB")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"prefix_precision_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"prefix_precision_{stamp}_{serial}.json"
    return path


def add(dic, key, value):
    if value != 0:
        dic[key] = dic.get(key, 0.) + float(value)


def quantize(vector, precision):
    """Fixed binary real/imaginary rounding, with no tolerance deletion."""
    step = 2.0 ** (-precision)
    return (np.rint(np.real(vector) / step) * step
            + 1j * np.rint(np.imag(vector) / step) * step)


def tv_law(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ArithmeticError("nonfinite law")
    if left.shape != right.shape or np.any(left < 0) or np.any(right < 0):
        raise ValueError("TV requires matching nonnegative arrays")
    if abs(float(left.sum())-1) > 5e-11 or abs(float(right.sum())-1) > 5e-11:
        raise AssertionError(f"law not normalized: {left.sum()} {right.sum()}")
    return float(.5*np.abs(left-right).sum())


def transition(circuit, key, stage, approximate, precision, cache, allow_fallback=True):
    """Return child keys and a normalized conditional law for one block."""
    a, exponent, p, output = key
    kind, index = stage
    if kind == "background":
        children = [(a, exponent, pp, output) for pp in range(circuit.b)]
        query = [(a, exponent, index, "background", 0, 0, pp) for pp in range(circuit.b)]
    elif kind == "reflection":
        q, _ = circuit.reflections[index]
        target = (-a-q) % circuit.sectors
        if target == a:
            return [key], np.array([1.]), {"exact_zero": 0, "raw_l2": 0., "fallback": False}
        children = [(x, exponent, p, output) for x in (a, target)]
        query = [(x, exponent, index, "reflection", 0, 0, p) for x in (a, target)]
    elif kind == "qft":
        measured = index
        ys = (output, output | (1 << (measured-1)))
        children = [(a, exponent, p, y) for y in ys]
        query = [(a, exponent, circuit.width, "reflection", measured, y, p) for y in ys]
    else:
        raise ValueError(f"unknown stage {kind}")
    weights = []
    exact_zero = 0
    raw_l2 = 0.
    for a0, e0, stop, boundary, measured, y, pp in query:
        cache_key = (a0, e0, stop, boundary, measured, y)
        if cache_key not in cache:
            exact = circuit.prefix_vector(a0, e0, stop, boundary=boundary,
                                          measured=measured, output=y)
            approx = quantize(exact, precision)
            cache[cache_key] = (exact, approx)
        exact, approx = cache[cache_key]
        weights.append(float(abs((approx if approximate else exact)[pp])**2))
        if float(abs(exact[pp])**2) == 0.:
            exact_zero += 1
        raw_l2 += float(abs(approx[pp]-exact[pp])**2)
    weights = np.asarray(weights, dtype=float)
    total = float(weights.sum())
    fallback = False
    if not np.isfinite(total) or total < 0:
        raise ArithmeticError("nonfinite rounded transition mass")
    if total == 0.:
        fallback = True
        if allow_fallback:
            probs = np.full(len(children), 1./len(children))
        else:
            probs = np.zeros(len(children))
    else:
        probs = weights / total
    if allow_fallback and abs(float(probs.sum())-1) > 5e-13:
        raise AssertionError("conditional law did not normalize")
    return children, probs, {"exact_zero": exact_zero, "raw_l2": math.sqrt(raw_l2),
                             "fallback": fallback}


def deterministic_update(frontier, bit, block):
    out = {}
    for (a, e, p, y), mass in frontier.items():
        key = (a, e, (p+bit) % block, y) if e & bit else (a, e, p, y)
        add(out, key, mass)
    return out


def apply_stage(circuit, exact_frontier, approx_frontier, stage, precision, cache):
    exact_next, approx_next = {}, {}
    weighted_tv = 0.
    max_local_tv = 0.
    worst_local_parent_mass = 0.
    fallback_count = 0
    exact_zero_count = 0
    raw_l2_weighted = 0.
    for parent, mass in exact_frontier.items():
        children, pe, meta = transition(circuit, parent, stage, False, precision, cache)
        _, pr, _ = transition(circuit, parent, stage, True, precision, cache)
        weighted_tv += mass * tv_law(pe, pr)
        local_tv = tv_law(pe, pr)
        if local_tv > max_local_tv:
            max_local_tv = local_tv
            worst_local_parent_mass = mass
        exact_zero_count += meta["exact_zero"]
        raw_l2_weighted += mass * meta["raw_l2"]
    # Propagate each law from its own current frontier.  Approximate states
    # outside exact support are retained; their conditional law is still
    # explicitly normalized and can only increase the measured discrepancy.
    for parent, mass in exact_frontier.items():
        children, probs, _ = transition(circuit, parent, stage, False, precision, cache)
        for child, prob in zip(children, probs):
            add(exact_next, child, mass*prob)
    for parent, mass in approx_frontier.items():
        children, probs, meta = transition(circuit, parent, stage, True, precision, cache)
        fallback_count += int(meta["fallback"])
        for child, prob in zip(children, probs):
            add(approx_next, child, mass*prob)
    if abs(sum(exact_next.values())-1) > 5e-11 or abs(sum(approx_next.values())-1) > 5e-11:
        raise AssertionError("stage law lost normalization")
    global_tv = .5*sum(abs(exact_next.get(k,0.)-approx_next.get(k,0.))
                        for k in set(exact_next) | set(approx_next))
    return exact_next, approx_next, dict(stage=f"{stage[0]}{stage[1]}",
        weighted_conditional_tv=weighted_tv, max_local_tv=max_local_tv,
        worst_local_parent_mass=worst_local_parent_mass,
        global_prefix_tv=float(global_tv), fallback_blocks=fallback_count,
        exact_zero_children=exact_zero_count, weighted_raw_l2=raw_l2_weighted,
        exact_mass=float(sum(exact_next.values())), approx_mass=float(sum(approx_next.values())))


def apply_discard_stage(circuit, frontier, stage, precision, cache):
    """Propagate the deliberately unnormalized zero-bin control."""
    nxt, dropped = {}, 0
    for parent, mass in frontier.items():
        children, probs, meta = transition(circuit, parent, stage, True,
                                            precision, cache, allow_fallback=False)
        if meta["fallback"]:
            dropped += 1
        for child, prob in zip(children, probs):
            add(nxt, child, mass*prob)
    return nxt, dropped


def full_prefix_stats(circuit, stage, precision):
    """Reference-only post-prefix Born laws and raw amplitude L2 error.

    The table is bounded by r<=14,t<=4 and is never used by either sampler.
    Sector/work coordinates are the indexed basis supplied to the prefix
    oracle; QFT rows retain the currently fixed low exponent bits and output
    prefix.  The sqrt(M) prefix convention is divided by sqrt(M) here.
    """
    kind, index = stage
    M, Q, b, t = circuit.sectors, 1 << circuit.width, circuit.b, circuit.width
    if circuit.period > CAP_R or t > CAP_T:
        raise ValueError("full prefix reference exceeds finite cap")
    guard((M,Q,b), np.complex128, "full prefix table")
    exact_law, rounded_law = {}, {}
    raw_sq = 0.
    max_coord_error = 0.
    max_scaled_coordinate_error = 0.
    if kind in ("background", "reflection"):
        e_values = range(Q)
        args = [(e, 0) for e in e_values]
        stop, boundary, measured = index, kind, 0
    else:
        e_values = range(1 << (t-index))
        args = [(e, y) for e in e_values for y in range(1 << index)]
        stop, boundary, measured = t, "reflection", index
    for gamma in range(M):
        for exponent, output in args:
            scaled = circuit.prefix_vector(gamma, exponent, stop,
                                           boundary=boundary, measured=measured,
                                           output=output)
            vector = scaled / np.sqrt(M)
            rounded = quantize(scaled, precision) / np.sqrt(M)
            scaled_delta = quantize(scaled, precision)-scaled
            max_scaled_coordinate_error = max(max_scaled_coordinate_error,
                float(np.max(abs(scaled_delta.real))), float(np.max(abs(scaled_delta.imag))))
            raw_sq += float(np.vdot(rounded-vector, rounded-vector).real)
            max_coord_error = max(max_coord_error,
                                  float(np.max(np.abs(rounded-vector))))
            for p in range(b):
                key = (gamma, exponent, p, output)
                exact_law[key] = float(abs(vector[p])**2)
                rounded_law[key] = float(abs(rounded[p])**2)
    exact_mass = sum(exact_law.values())
    rounded_mass = sum(rounded_law.values())
    if abs(exact_mass-1) > 5e-10:
        raise AssertionError(f"reference prefix law not normalized: {stage} {exact_mass}")
    full_fallback = False
    if rounded_mass <= 0 or not np.isfinite(rounded_mass):
        # Diagnostic-only full-law counterpart of the block fallback: retain
        # a declared normalized law so TV remains meaningful. Production never
        # constructs this table and never rounds amplitudes.
        full_fallback = True
        rounded_norm = {key: 1./len(rounded_law) for key in rounded_law}
    else:
        rounded_norm = {key: value/rounded_mass for key, value in rounded_law.items()}
    born_tv = tv_law([v/exact_mass for v in exact_law.values()],
                    [rounded_norm[key] for key in exact_law])
    return dict(stage=f"{kind}{index}", exact_mass=float(exact_mass),
                rounded_raw_mass=float(rounded_mass), full_born_tv=float(born_tv),
                raw_global_l2=math.sqrt(raw_sq), full_fallback=full_fallback,
                max_coord_error=max_coord_error,
                max_scaled_real_imag_error=max_scaled_coordinate_error,
                uniform_grid_state_l2_bound=math.sqrt(2*b*Q)*2.**(-precision-1),
                table_entries=len(exact_law))


def run_precision(circuit, target, precision):
    if circuit.period > CAP_R or circuit.width > CAP_T:
        raise ValueError("precision reference exceeds r<=14,t<=4 cap")
    guard((circuit.sectors, 1 << circuit.width, circuit.b, 1 << circuit.width),
          label="precision state budget")
    exact = {(a, e, 0, 0): 1./(circuit.sectors*(1 << circuit.width))
             for a in range(circuit.sectors) for e in range(1 << circuit.width)}
    approx = dict(exact)
    discard = dict(exact)
    cache = {}
    stages = []
    for stop in range(circuit.width+1):
        if stop:
            bit = 1 << (stop-1)
            exact = deterministic_update(exact, bit, circuit.b)
            approx = deterministic_update(approx, bit, circuit.b)
            discard = deterministic_update(discard, bit, circuit.b)
        if stop in circuit.background.defects:
            stage = ("background", stop)
            exact, approx, row = apply_stage(circuit, exact, approx,
                                              stage, precision, cache)
            discard, dropped = apply_discard_stage(circuit, discard, stage, precision, cache)
            row["discard_dropped_blocks"] = dropped
            row.update(full_prefix_stats(circuit, stage, precision))
            stages.append(row)
        if stop in circuit.reflections:
            stage = ("reflection", stop)
            exact, approx, row = apply_stage(circuit, exact, approx,
                                              stage, precision, cache)
            discard, dropped = apply_discard_stage(circuit, discard, stage, precision, cache)
            row["discard_dropped_blocks"] = dropped
            row.update(full_prefix_stats(circuit, stage, precision))
            stages.append(row)
    for measured in range(1, circuit.width+1):
        stage = ("qft", measured)
        exact, approx, row = apply_stage(circuit, exact, approx,
                                          stage, precision, cache)
        discard, dropped = apply_discard_stage(circuit, discard, stage, precision, cache)
        row["discard_dropped_blocks"] = dropped
        row.update(full_prefix_stats(circuit, stage, precision))
        stages.append(row)
    final = np.zeros_like(target)
    rounded = np.zeros_like(target)
    for (a, e, p, y), mass in exact.items():
        final[a, y] += mass
    for (a, e, p, y), mass in approx.items():
        rounded[a, y] += mass
    rounded_mass = float(rounded.sum())
    if abs(rounded_mass-1) > 5e-11:
        raise AssertionError("rounded final law not normalized")
    exact_reference_tv = tv_law(target, final)
    rounded_final_tv = tv_law(target, rounded)
    weighted_sum = sum(row["weighted_conditional_tv"] for row in stages)
    full_born_sum = sum(row["full_born_tv"] for row in stages)
    local_full_holds = all(row["weighted_conditional_tv"] <=
                           2*row["full_born_tv"] + 5e-12 for row in stages)
    raw_global_l2_sum = sum(row["raw_global_l2"] for row in stages)
    final_comparator = np.zeros_like(target)
    for gamma in range(circuit.sectors):
        for y in range(1 << circuit.width):
            exact_vec = circuit.prefix_vector(gamma, 0, circuit.width,
                                              measured=circuit.width, output=y)
            rounded_vec = quantize(exact_vec, precision)
            final_comparator[gamma, y] = float(np.vdot(rounded_vec, rounded_vec).real)
    final_comparator /= float(circuit.sectors)
    comparator_mass = float(final_comparator.sum())
    if not np.isfinite(comparator_mass) or comparator_mass < 0:
        raise ArithmeticError("nonfinite rounded final amplitude mass")
    if comparator_mass == 0:
        final_comparator[:] = 1./final_comparator.size
    else:
        final_comparator /= comparator_mass
    comparator_tv = tv_law(target, final_comparator)
    return dict(precision=precision, final_tv=rounded_final_tv,
                exact_reference_tv=exact_reference_tv,
                rounded_mass=rounded_mass, weighted_tv_sum=weighted_sum,
                full_born_tv_sum=full_born_sum,
                conservative_bound=min(1.,2*full_born_sum),
                local_full_born_bound_holds=local_full_holds,
                chain_weighted_bound_holds=rounded_final_tv <= weighted_sum + 5e-12,
                bound_holds=rounded_final_tv <= weighted_sum + 5e-12 and local_full_holds
                    and all(row["full_born_tv"] <= row["raw_global_l2"]+5e-12
                            and row["raw_global_l2"] <= row["uniform_grid_state_l2_bound"]+5e-12
                            and row["max_scaled_real_imag_error"] <= 2.**(-precision-1)+5e-12
                            for row in stages),
                rounded_final_comparator_tv=comparator_tv,
                rounded_final_comparator_raw_mass=comparator_mass,
                comparator_global_l2_bound=stages[-1]["raw_global_l2"],
                fallback_blocks=sum(row["fallback_blocks"] for row in stages),
                exact_zero_children=sum(row["exact_zero_children"] for row in stages),
                max_local_tv=max(row["max_local_tv"] for row in stages),
                max_global_prefix_tv=max(row["global_prefix_tv"] for row in stages),
                raw_global_l2_sum=raw_global_l2_sum,
                max_raw_global_l2=max(row["raw_global_l2"] for row in stages),
                discard_only_final_mass=float(sum(discard.values())),
                discard_dropped_blocks=sum(row["discard_dropped_blocks"] for row in stages),
                stage_rows=stages, prefix_query_count=len(cache))


def main():
    started = time.time()
    report = {"status":"PASS", "predictions":"P1/P2/P3 with C1/C2 must-fail",
              "rows":[], "controls":{}, "reference":"existing direct_joint and prefix_vector; r10,b2,t4"}
    try:
        circuit = fixture()
        target = direct_joint(circuit)
        rows = [run_precision(circuit, target, precision) for precision in PRECISIONS]
        report["rows"] = rows
        low = next(row for row in rows if row["precision"] == 0)
        finest = rows[-1]
        guard((64,), np.complex128, "pointwise control")
        control_psi = np.ones(64)/8
        control_phi = np.concatenate((np.zeros(32),np.ones(32)/np.sqrt(32)))
        control_tv = tv_law(abs(control_psi)**2, abs(control_phi)**2)
        coordinate_error = float(np.max(abs(control_phi-control_psi)))
        report["controls"] = {
            "pointwise_only_bound": coordinate_error,
            "pointwise_only_actual_tv": control_tv,
            "pointwise_only_fails": control_tv > coordinate_error,
            "discard_only_law_mass": low["discard_only_final_mass"],
            "uniform_fallback_blocks_at_p0": low["fallback_blocks"],
        }
        p1 = all(row["bound_holds"] for row in rows)
        p2 = (all(abs(row["rounded_mass"]-1) < 5e-11 for row in rows)
              and finest["final_tv"] < low["final_tv"])
        p3 = all(row["rounded_final_comparator_tv"] <= row["comparator_global_l2_bound"] + 5e-12
                 for row in rows)
        exp.check("P1", p1 and all(row["exact_reference_tv"]<5e-12 for row in rows),
                  f"all full-law/weighted/global-L2 budgets and independent final reference hold={p1}")
        exp.check("P2", p2 and low["fallback_blocks"] >= 1,
                  f"normalized/finest-improves={p2}, explicit p0 fallback blocks={low['fallback_blocks']}")
        exp.check("P3", p3, f"rounded-final amplitude comparator/global-L2 bound={p3}")
        exp.fail_check("C1", report["controls"]["pointwise_only_fails"],
                       f"independent normalized control TV={control_tv:.6g} > max coordinate error={coordinate_error:.6g}")
        # A deliberately unnormalized discard-only law is measured as a failed
        # certificate, while the production/fallback law above remains normalized.
        discard_mass = report["controls"]["discard_only_law_mass"]
        exp.fail_check("C2", abs(discard_mass-1.) > 1e-15 and low["fallback_blocks"] >= 1,
                       f"discard-only mass {discard_mass}; explicit fallback is required")
        if not exp.finish():
            raise AssertionError("Experiment harness failed")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time()-started
    path = report_path()
    path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"status":report["status"], "report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
