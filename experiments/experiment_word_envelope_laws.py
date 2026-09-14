"""Bounded law audit of the chronological word-specific support envelope.

This is the first law-level check of C73's E/O word plan.  It compares the
existing sparse reverse instrument with ``support_mode='alphabet'`` and
``support_mode='word'`` on the same supplied gates, and checks every forced
finite-bit reverse submass against the independent full-r reference.  The
word plan is setup-only metadata; it does not retain offset sets or enumerate
reflection histories.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Alphabet and word modes have identical complex prefixes and complete
      accepted joint laws; each accepted submass normalizes to 1/(b*S).
  P2  The word setup reports finite released-offset work, zero retained
      offsets/history enumeration, and preserves the independent sequential
      output marginal on all tiny and the bounded k=9 law.
  P3  The chronological (0,1,0) word has a strictly smaller global envelope
      than the two-label alphabet bound on the unsaturated r=14,b=2 fixture.
  C1  A gamma=0-only reached-set count is a valid global envelope on that
      fixture.

This is a complex128 diagnostic, not a finite-TV certificate or timing claim.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionInput
from lab.coherent_reverse import SparseCoherentReverse
from lab.merged_prefix import MergedCoherentPrefixes


MAX_BYTES = 16 << 20
TOL = 5e-11
MAX_FORCED_ROWS = 200_000
MAX_FRONTIER_WORK = 1_000_000_000

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "alphabet and chronological-word laws/prefixes agree")
exp.predict("P2", "word setup is released, history-free, and matches sequential output")
exp.predict("P3", "the unsaturated chronological word strictly tightens the envelope")
exp.must_fail("C1", "gamma=0-only reached-set count is a valid global envelope")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"word_envelope_laws_{stamp}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def guard(shape, dtype=np.float64, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2),
                    np.exp(1j * theta / 2)]).astype(complex)


def embed(matrix, block):
    result = np.eye(block, dtype=complex)
    result[:2, :2] = matrix
    return result


def fixture(period, block, width, route_labels, *, angle=np.pi / 4):
    defects = {s: embed(rx(np.pi / 7) if s % 2 else rz(np.pi / 5), block)
               for s in range(width + 1)}
    reflections = {int(s): (int(q), float(angle))
                   for s, q in route_labels.items()}
    return CoherentReflectionInput(period, block, width, defects, reflections)


def tiny_fixtures():
    rows = []
    for block, period in ((2, 6), (3, 9)):
        for width in (0, 1, 4):
            if width == 0:
                labels = {}
            elif width == 1:
                labels = {0: 0}
            else:
                labels = {0: 0, 1: 1, 3: 0, 4: 1}
            rows.append((f"r{period}_b{block}_t{width}",
                         fixture(period, block, width, labels)))
    # M=4 is even, q-difference 2 is noninvertible modulo M, and all route
    # rotations are exactly zero.  The structural support plan must still be
    # valid even though the physical reflection is the identity.
    rows.append(("even_r12_b3_t4_noninvertible_zero", fixture(
        12, 3, 4, {0: 0, 1: 2, 3: 0, 4: 2}, angle=0.0)))
    rows.append(("even_r12_b3_t4_noninvertible_nonzero", fixture(
        12, 3, 4, {0: 0, 1: 2, 3: 0, 4: 2})))
    return rows


def unsaturated_fixture():
    # Chronological q=(0,1,0), r=14,b=2,M=7: word S=6 while the
    # two-label alphabet envelope is min(7,1+2*3)=7.
    return fixture(14, 2, 3, {0: 0, 1: 1, 3: 0})


def wide_fixture():
    return fixture(9, 3, 8, {s: s % 3 for s in range(9)})


def preflight(circuit):
    # Use the actual shared guard before any joint-law/frontier allocation.
    from experiments.experiment_coherent_route_sampling import reference_guard
    reference = reference_guard(circuit, max_width=8)
    guard((circuit.sectors, 1 << circuit.width), label="joint law")
    r, q, t = circuit.period, 1 << circuit.width, circuit.width
    estimated_frontier = (4 * circuit.sectors * circuit.sectors
                          * circuit.b * circuit.b * q
                          * max(1, t + 1 + 2 * q))
    forced_rows = circuit.sectors * circuit.b * q
    if forced_rows > MAX_FORCED_ROWS:
        raise MemoryError("forced reverse-row budget exceeded")
    if estimated_frontier > MAX_FRONTIER_WORK:
        raise MemoryError("frontier scalar-visit budget exceeded")
    return {"period": r, "block": circuit.b, "width": t,
            **reference,
            "estimated_transition_frontier_work": estimated_frontier,
            "forced_rows": forced_rows,
            "max_forced_rows": MAX_FORCED_ROWS,
            "max_frontier_work": MAX_FRONTIER_WORK}


def counted_reverse_sample(worker):
    fields = ("reverse_steps", "work_matvecs", "reflection_vector_contributions",
              "qft_matrix_pair_constructions")
    counts = dict.fromkeys(fields, 0)
    count = 0
    original = worker.attempt
    def counted(*args, **kwargs):
        nonlocal count
        result = original(*args, **kwargs)
        count += 1
        for key in fields:
            counts[key] += result[key]
        return result
    worker.attempt = counted
    try:
        result = worker.sample(np.random.default_rng(719), max_proposals=100_000)
    finally:
        worker.attempt = original
    assert count == result["rejection_proposals"]
    assert all(result[key] == value for key, value in counts.items())
    return result


def counted_prefix_sample(worker):
    from experiments.experiment_merged_prefix_sampling import expected_prefix_cost
    counts = dict.fromkeys(worker._counts, 0)
    original = worker.prefix_vector
    def counted(*args, **kwargs):
        expected = expected_prefix_cost(worker.circuit, *args, **kwargs)
        result = original(*args, **kwargs)
        for key, value in expected.items():
            assert worker.last_prefix_stats[key] == value
            counts[key] = (max(counts[key], value) if key == "peak_sector_count"
                           else counts[key] + value)
        return result
    worker.prefix_vector = counted
    try:
        result = worker.sample(np.random.default_rng(731))
    finally:
        worker.prefix_vector = original
    assert all(result[key] == value for key, value in counts.items())
    assert result["prefix_vector_evaluations"] == counts["merged_prefix_queries"]
    return result


def prefix_audit(circuit):
    from experiments.experiment_coherent_route_sampling import direct_prefix
    calls = {"count": 0}
    had = hasattr(CoherentReflectionInput, "_histories")
    original_history = getattr(CoherentReflectionInput, "_histories", None)
    def forbidden_history(*args, **kwargs):
        calls["count"] += 1
        raise AssertionError("word prefix adapter called forbidden _histories")
    CoherentReflectionInput._histories = forbidden_history
    alpha = MergedCoherentPrefixes(circuit, support_mode="alphabet")
    word = MergedCoherentPrefixes(circuit, support_mode="word")
    max_mode = max_direct = 0.0
    requests = 0
    try:
        for stop in range(circuit.width + 1):
            for boundary in ("arithmetic", "background", "reflection"):
                measured_values = (range(circuit.width + 1)
                                   if stop == circuit.width and boundary == "reflection"
                                   else (0,))
                for measured in measured_values:
                    low = min(stop, circuit.width - measured)
                    for gamma in range(circuit.sectors):
                        for exponent in range(1 << low):
                            for output in range(1 << measured):
                                a = alpha.prefix_vector(gamma, exponent, stop,
                                    boundary=boundary, measured=measured, output=output)
                                w = word.prefix_vector(gamma, exponent, stop,
                                    boundary=boundary, measured=measured, output=output)
                                ref = direct_prefix(circuit, gamma, exponent, stop,
                                    boundary, measured, output, max_width=8)
                                max_mode = max(max_mode, float(np.max(np.abs(a - w))))
                                max_direct = max(max_direct,
                                                 float(np.max(np.abs(a - ref))),
                                                 float(np.max(np.abs(w - ref))))
                                requests += 1
    finally:
        if had:
            CoherentReflectionInput._histories = original_history
        else:
            delattr(CoherentReflectionInput, "_histories")
    return {"requests": requests, "mode_max_abs_error": max_mode,
            "direct_max_abs_error": max_direct,
            "history_sentinel_calls": calls["count"],
            "alphabet_stats": alpha.stats(), "word_stats": word.stats()}


def reverse_law(circuit, mode):
    from experiments.experiment_coherent_route_sampling import direct_joint, full_gates, sequential_path
    preflight(circuit)
    guard((circuit.sectors, 1 << circuit.width), label=f"{mode} law")
    calls = {"count": 0}
    had = hasattr(CoherentReflectionInput, "_histories")
    original_history = getattr(CoherentReflectionInput, "_histories", None)
    def forbidden_history(*args, **kwargs):
        calls["count"] += 1
        raise AssertionError("word reverse adapter called forbidden _histories")
    CoherentReflectionInput._histories = forbidden_history
    worker = SparseCoherentReverse(circuit, support_mode=mode)
    M, b, Q = circuit.sectors, circuit.b, 1 << circuit.width
    proposal = np.zeros((M, Q), float)
    accepted = np.zeros((M, Q), float)
    max_submass_error = max_acceptance_error = 0.0
    max_peak = 0
    zero_terminal_rows = 0
    peak_violation_rows = 0
    forced_rows = 0
    started = time.perf_counter()
    for gamma in range(M):
        for boundary in range(b):
            for output in range(Q):
                result = worker.attempt(gamma, boundary, output=output)
                q = float(result["proposal_path_probability"])
                a = float(result["acceptance_probability"])
                sub = float(result["accepted_joint_submass"])
                proposal[gamma, output] += q / (M * b)
                accepted[gamma, output] += sub
                max_submass_error = max(max_submass_error,
                    abs(sub - q * a / (M * b)))
                max_peak = max(max_peak, int(result["peak_sector_count"]))
                if int(result["peak_sector_count"]) > int(result["support_bound"]):
                    peak_violation_rows += 1
                terminal_norm = float(result["terminal_norm"])
                terminal_num = float(result["terminal_coherent_numerator"])
                if terminal_norm == 0.0:
                    # Exact zero proposal branches are allowed and must not
                    # be converted into a numerical fallback.  The forced
                    # result must carry zero acceptance and zero numerator.
                    zero_terminal_rows += 1
                    if a != 0.0 or terminal_num != 0.0 or q != 0.0 or sub != 0.0:
                        raise AssertionError("zero terminal branch received nonzero proposal/acceptance")
                else:
                    max_acceptance_error = max(max_acceptance_error,
                        abs(a - terminal_num
                        / (float(result["support_bound"]) * terminal_norm)))
                if int(result["history_enumerations"]) != 0:
                    raise AssertionError("reverse law enumerated histories")
                forced_rows += 1
    # Check exact row cardinality and setup metadata before doing independent
    # dense/reference allocations.
    expected_rows = M * b * Q
    if forced_rows != expected_rows:
        raise AssertionError(f"forced-row count {forced_rows} != {expected_rows}")
    setup = worker.support_plan
    setup_ok = (int(setup["support_bound"]) <= int(setup["alphabet_support_bound"])
                and int(setup["word_retained_offset_entries"]) == 0)
    if mode == "word":
        setup_ok &= (int(setup["word_setup_steps"]) == len(circuit.reflections)
                     and int(setup["word_setup_copied_offsets"]) >= 0
                     and int(setup["word_setup_reflected_offsets"]) >= 0
                     and int(setup["word_setup_peak_live_offset_entries"]) >= 1)
    else:
        setup_ok &= all(int(setup[key]) == 0 for key in (
            "word_setup_steps", "word_setup_copied_offsets",
            "word_setup_reflected_offsets", "word_setup_peak_live_offset_entries"))
    if not setup_ok or peak_violation_rows:
        raise AssertionError("support setup/peak metadata invalid")
    S = int(worker.support_bound)
    worker_sample = counted_reverse_sample(worker)
    worker_sample_ok = (
        1 <= int(worker_sample["rejection_proposals"]) <= 100_000
        and int(worker_sample["history_enumerations"]) == 0
        and int(worker_sample["peak_sector_count"]) <= S
        and int(worker_sample["reverse_steps"])
            == int(worker_sample["rejection_proposals"]) * circuit.width)
    if not worker_sample_ok:
        raise AssertionError("reverse returned-sample counter audit failed")
    # Restore the sentinel before constructing the merged adapter, then audit
    # its actual one-sample counters under a fresh instance.
    if had:
        CoherentReflectionInput._histories = original_history
    else:
        delattr(CoherentReflectionInput, "_histories")
    merged_calls = {"count": 0}
    merged_had = hasattr(CoherentReflectionInput, "_histories")
    merged_original = getattr(CoherentReflectionInput, "_histories", None)
    def merged_forbidden(*args, **kwargs):
        merged_calls["count"] += 1
        raise AssertionError("word merged adapter called forbidden _histories")
    CoherentReflectionInput._histories = merged_forbidden
    merged = MergedCoherentPrefixes(circuit, support_mode=mode)
    merged_sample = counted_prefix_sample(merged)
    if merged_had:
        CoherentReflectionInput._histories = merged_original
    else:
        delattr(CoherentReflectionInput, "_histories")
    merged_sample_ok = (merged_calls["count"] == 0
                        and int(merged_sample["rejection_proposals"]) == 0
                        and int(merged_sample["history_enumerations"]) == 0
                        and int(merged_sample["peak_sector_count"]) <= S
                        and int(merged_sample["prefix_vector_evaluations"]) > 0
                        and int(merged_sample["prefix_vector_evaluations"])
                            == int(merged_sample["merged_prefix_queries"]))
    if not merged_sample_ok:
        raise AssertionError("merged returned-sample counter audit failed")
    target = direct_joint(circuit, max_width=8)
    pairs, initial = full_gates(circuit, max_width=8)
    sequential = np.array([sequential_path(pairs, initial, output=y)
                            ["conditional_path_probability"] for y in range(Q)])
    accepted_mass = float(accepted.sum())
    proposal_mass = float(proposal.sum())
    S = int(worker.support_bound)
    accepted_normalized = accepted / accepted_mass
    return {"mode": mode, "support_bound": S,
            "support_plan": worker.support_plan,
            "forced_rows": forced_rows,
            "expected_forced_rows": expected_rows,
            "accepted_mass": accepted_mass,
            "expected_accepted_mass": 1.0 / (b * S),
            "accepted_mass_error": abs(accepted_mass - 1.0 / (b * S)),
            "proposal_mass": proposal_mass,
            "proposal_mass_error": abs(proposal_mass - 1.0),
            "accepted_vs_target_tv": float(np.abs(accepted_normalized - target).sum() / 2),
            "target_mass_error": abs(float(target.sum()) - 1.0),
            "target_vs_sequential_tv": float(np.abs(target.sum(axis=0) - sequential).sum() / 2),
            "max_submass_identity_error": max_submass_error,
            "max_acceptance_identity_error": max_acceptance_error,
            "zero_terminal_rows": zero_terminal_rows,
            "peak_violation_rows": peak_violation_rows,
            "history_sentinel_calls": calls["count"],
            "merged_sample": {k: merged_sample[k] for k in (
                "output", "final_coarse_sector", "prefix_vector_evaluations",
                "rejection_proposals", "history_enumerations", "peak_sector_count")},
            "worker_sample": {k: worker_sample[k] for k in (
                "output", "final_coarse_sector", "rejection_proposals",
                "history_enumerations", "reverse_steps", "peak_sector_count")},
            "worker_sample_ok": worker_sample_ok,
            "merged_history_sentinel_calls": merged_calls["count"],
            "merged_sample_ok": merged_sample_ok,
            "setup_metadata_ok": setup_ok,
            "peak_sector_count": max_peak,
            "history_enumerations": 0,
            "worker_stats": worker._stats(),
            "seconds": time.perf_counter() - started,
            "target": target.tolist(), "accepted_normalized": accepted_normalized.tolist()}


def reached_set(sectors, qs, gamma):
    current = {int(gamma)}
    for q in reversed(tuple(int(q) % sectors for q in qs)):
        current |= {(-alpha - q) % sectors for alpha in current}
    return current


def main():
    report = {"status": "PASS", "preflight": [], "prefixes": [],
              "laws": [], "controls": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        small = tiny_fixtures()
        unsat = unsaturated_fixture()
        wide = wide_fixture()
        all_cases = small + [("unsaturated_r14_b2_t3", unsat),
                             ("wide_r9_b3_k9_t8", wide)]
        report["preflight"] = [dict(name=name, **preflight(circuit))
                                for name, circuit in all_cases]

        # Prefixes are fully enumerated on the tiny cases only; the k=9 law
        # remains bounded by the explicit dense/reference preflight.
        for name, circuit in small + [("unsaturated_r14_b2_t3", unsat)]:
            row = dict(name=name, **prefix_audit(circuit))
            report["prefixes"].append(row)
            p1 &= (row["mode_max_abs_error"] < TOL
                   and row["direct_max_abs_error"] < TOL)
            p2 &= row["history_sentinel_calls"] == 0

        for name, circuit in all_cases:
            for mode in ("alphabet", "word"):
                row = reverse_law(circuit, mode)
                row["name"] = name
                report["laws"].append(row)
                p1 &= (row["accepted_mass_error"] < TOL
                       and row["accepted_vs_target_tv"] < TOL
                       and row["target_mass_error"] < TOL)
                p2 &= (row["proposal_mass_error"] < TOL
                       and row["target_vs_sequential_tv"] < TOL
                       and row["max_submass_identity_error"] < TOL
                       and row["max_acceptance_identity_error"] < TOL
                       and row["history_enumerations"] == 0
                       and row["peak_violation_rows"] == 0
                       and row["forced_rows"] == row["expected_forced_rows"]
                       and row["setup_metadata_ok"]
                       and row["worker_sample_ok"]
                       and row["merged_sample_ok"]
                       and row["history_sentinel_calls"] == 0
                       and row["merged_history_sentinel_calls"] == 0)

        word_unsat = next(x for x in report["laws"]
                          if x["name"] == "unsaturated_r14_b2_t3" and x["mode"] == "word")
        alpha_unsat = next(x for x in report["laws"]
                           if x["name"] == "unsaturated_r14_b2_t3" and x["mode"] == "alphabet")
        p3 = (word_unsat["support_bound"] < alpha_unsat["support_bound"]
              and word_unsat["support_plan"]["word_retained_offset_entries"] == 0
              and word_unsat["support_plan"]["word_setup_steps"] == 3)

        # Meaningful must-fail: on the unsaturated fixture, the gamma=0
        # reached set has size 3 but other target sectors require 6.  This is
        # a set-level invalidity of a local shortcut, not a forced acceptance
        # overflow claim.
        qs = (0, 1, 0)
        gamma_sizes = [len(reached_set(7, qs, gamma)) for gamma in range(7)]
        gamma0_shortcut_valid = max(gamma_sizes) <= gamma_sizes[0]
        report["controls"]["gamma0_only"] = {
            "sizes_by_gamma": gamma_sizes,
            "gamma0_size": gamma_sizes[0],
            "global_max_size": max(gamma_sizes),
            "shortcut_claim": gamma0_shortcut_valid}

        exp.check("P1", p1, "prefixes, accepted submasses, and dense laws")
        exp.check("P2", p2, "released word setup, sequential marginal, and no histories")
        exp.check("P3", p3, "unsaturated word envelope is strictly smaller")
        exp.fail_check("C1", not gamma0_shortcut_valid,
                       f"gamma sizes={gamma_sizes}")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "PASS" if all(ok for _, _, ok, _ in exp._results) else "FAIL"
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "dense_payload_cap_bytes": MAX_BYTES, "tolerance": TOL,
        "wide_fixture": "r9,b3,k9,t8", "float_only": True,
        "word_offsets_retained": False})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
