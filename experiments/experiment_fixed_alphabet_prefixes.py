"""Bounded fixed-alphabet prefix and complete-law audit beyond eight routes.

This experiment uses the structural ``CoherentReflectionInput`` (which has no
history API) with nine insertions and a cyclic three-label alphabet.  A capped
test-only cache canonicalizes the irrelevant high exponent bits while retaining
all boundary, measured-prefix, and output coordinates.  The resulting merged
prefixes and complete transition laws are compared with the existing dense
full-r reference at opt-in width cap eight.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  All distinct cached complex prefixes for structural b=2 and b=3 cases
      agree with direct full-r contraction, and complete transition laws
      normalize and agree on at least the b=3 nine-insertion case.
  P2  The structural path has no history API/calls, reports bounded merged
      support and charged cache entries, and agrees with an equivalent small
      legacy k<=8 fixture.
  P3  Intermediate/end noncommuting mixers remain visible in the complete law;
      the default width guard and legacy k>8 history cap reject wider calls.
  C1  A wrong internal QFT phase changes the complete law.
  C2  A wrong reflection boundary changes the complete law.
  C3  Removing the explicit opt-in/default/legacy guards succeeds.

The dense reference and transition/frontier budgets are explicit opt-in caps;
no production cache, history enumeration, or broad timing claim is made.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit, CoherentReflectionInput
from lab.merged_prefix import MergedCoherentPrefixes


MAX_BYTES = 16 << 20
TARGET_TOL = 5e-11
MAX_CACHE_ENTRIES = 16_384

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "structural prefixes and complete laws agree with opt-in dense references")
exp.predict("P2", "structural history-free counters/cache and legacy equivalence hold")
exp.predict("P3", "intermediate effects remain visible and explicit guards reject unsupported defaults")
exp.must_fail("C1", "wrong internal QFT phase preserves the complete law")
exp.must_fail("C2", "wrong reflection boundary preserves the complete law")
exp.must_fail("C3", "default/legacy/reference guards reject no unsupported call")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"fixed_alphabet_prefixes_{stamp}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def guard_entries(entries, bytes_per_entry=1024):
    if type(entries) is not int or entries < 0 or entries > MAX_CACHE_ENTRIES:
        raise MemoryError("prefix cache entry cap exceeded")
    if entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("prefix cache payload exceeds 16 MiB")


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)]).astype(complex)


def embed(matrix, block):
    out = np.eye(block, dtype=complex)
    out[:2, :2] = matrix
    return out


def structural_fixture(block, width=8, *, route_count=9, intermediate=True):
    period = 6 if block == 2 else 9
    defects = {0: embed(rx(np.pi / 4), block)}
    if intermediate:
        for s in range(1, width + 1):
            defects[s] = embed(rx(np.pi / 7) if s % 2 else rz(np.pi / 5), block)
    reflections = {s: (i % 3, np.pi / 4) for i, s in enumerate(range(route_count))}
    return CoherentReflectionInput(period, block, width, defects, reflections)


def reference_guard_info(circuit):
    from experiments.experiment_coherent_route_sampling import reference_guard
    r, q, t = circuit.period, 1 << circuit.width, circuit.width
    # Same conservative scalar-visit estimate enforced by
    # enumerate_transitions; the simultaneously retained frontier is much
    # smaller and is reported separately by the production helper.
    estimated_frontier = (4 * (r // circuit.b) ** 2 * circuit.b ** 2 * q
                          * max(1, t + 1 + 2 * q))
    return {**reference_guard(circuit,max_width=8),
            "period": r, "block": circuit.b, "width": t,
            "estimated_transition_frontier_work": estimated_frontier,
            "two_frontier_entry_allowance_bytes": 2*r*q*1024,
            "max_width": 8}


def forbidden_history_classpatch():
    calls = {"count": 0}
    had = hasattr(CoherentReflectionInput, "_histories")
    original = getattr(CoherentReflectionInput, "_histories", None)
    def sentinel(*args, **kwargs):
        calls["count"] += 1
        raise AssertionError("structural merged path called forbidden _histories")
    CoherentReflectionInput._histories = sentinel
    return calls, had, original


def restore_history_classpatch(had, original):
    if had:
        CoherentReflectionInput._histories = original
    else:
        delattr(CoherentReflectionInput, "_histories")


def prefix_audit(supplied):
    from experiments.experiment_coherent_route_sampling import direct_prefix
    from experiments.experiment_merged_prefix_sampling import expected_prefix_cost
    merged = MergedCoherentPrefixes(supplied)
    cache = {}
    max_error = 0.0
    requests = 0
    started = time.perf_counter()
    t, M = supplied.width, supplied.sectors
    expected_unique = M * (3 * (2 ** (t + 1) - 1) + t * (2 ** t))
    guard_entries(expected_unique)
    totals = dict.fromkeys(expected_prefix_cost(merged,0,0,0),0)
    aliases,alias_error = 0,0.
    for stop in range(t + 1):
        for boundary in ("arithmetic", "background", "reflection"):
            measured_values = range(t + 1) if stop == t and boundary == "reflection" else (0,)
            for measured in measured_values:
                low = min(stop, t - measured)
                outputs = range(1 << measured)
                for gamma in range(M):
                    for exponent in range(1 << low):
                        for output in outputs:
                            key = (gamma, exponent, stop, boundary, measured, output)
                            if key not in cache:
                                # Guard before constructing/retaining the next
                                # pair of complex vectors.  The 1024-byte
                                # allowance is conservative Python/cache
                                # bookkeeping; actual ndarray payload is
                                # accounted separately below.
                                guard_entries(len(cache) + 1)
                                got = merged.prefix_vector(gamma, exponent, stop,
                                    boundary=boundary, measured=measured, output=output)
                                expected_cost = expected_prefix_cost(merged,gamma,exponent,stop,
                                    boundary=boundary,measured=measured,output=output)
                                for key_cost,value_cost in expected_cost.items():
                                    totals[key_cost] = (max(totals[key_cost],value_cost)
                                        if key_cost=="peak_sector_count" else totals[key_cost]+value_cost)
                                if not all(merged.last_prefix_stats[key_cost]==value_cost
                                           for key_cost,value_cost in expected_cost.items()):
                                    raise AssertionError("independent prefix operation accounting mismatch")
                                ref = direct_prefix(supplied, gamma, exponent, stop,
                                    boundary, measured, output, max_width=8)
                                cache[key] = (got, ref)
                            got, ref = cache[key]
                            max_error = max(max_error, float(np.max(np.abs(got - ref))))
                            requests += 1
    # Independently test selected noncanonical API requests with every ignored
    # exponent bit set. Keep their counters separate from the audited sweep.
    stats = merged.stats()
    for index,(key,(value,_)) in enumerate(cache.items()):
        gamma,exponent,stop,boundary,measured,output = key
        low = min(stop,t-measured)
        if index%97 == 0 and low<t:
            alternate = exponent | (((1 << t)-1)^((1 << low)-1))
            alias_error = max(alias_error,float(np.max(np.abs(value-merged.prefix_vector(
                gamma,alternate,stop,boundary=boundary,measured=measured,output=output)))))
            aliases += 1
    actual_bytes = len(cache) * (2 * supplied.b * np.dtype(np.complex128).itemsize)
    conservative_bytes = len(cache) * 1024
    return {"max_abs_error": max_error, "cache_entries": len(cache),
            "expected_unique_entries": expected_unique,
            "cache_conservative_bytes": conservative_bytes,
            "cache_actual_numpy_bytes": actual_bytes,
            "requests": requests, "seconds": time.perf_counter() - started,
            "cost_exact": all(stats[key]==value for key,value in totals.items()),
            "alias_requests": aliases, "alias_max_abs_error": alias_error,
            "support_bound": merged.support_bound,
            "stats": stats}


def law_audit(supplied):
    from experiments.experiment_coherent_route_sampling import direct_joint, enumerate_transitions
    from experiments.experiment_merged_prefix_sampling import expected_prefix_cost
    merged = MergedCoherentPrefixes(supplied)
    # enumerate_transitions revisits the same canonical prefix many times.
    # Install a test-only memoizer on this instance; high exponent bits are
    # irrelevant once the prefix/QFT boundary has fixed the requested low
    # bits, but all other request coordinates remain part of the key.
    cache = {}
    original_prefix = merged.prefix_vector
    def cached_prefix(sector, exponent, stop, *, boundary="reflection",
                      measured=0, output=0):
        low = min(int(stop), int(merged.width) - int(measured))
        canonical = int(exponent) & ((1 << low) - 1) if low else 0
        key = (int(sector), canonical, int(stop), str(boundary),
               int(measured), int(output))
        if key not in cache:
            guard_entries(len(cache) + 1)
            cache[key] = original_prefix(int(sector), canonical, int(stop),
                                         boundary=boundary, measured=int(measured),
                                         output=int(output))
        return cache[key]
    merged.prefix_vector = cached_prefix
    reference = CoherentReflectionInput(supplied.period, supplied.b, supplied.width,
        dict(supplied.background.defects), dict(supplied.reflections))
    calls, had, original = forbidden_history_classpatch()
    try:
        target = direct_joint(reference, max_width=8)
        enum = enumerate_transitions(merged, max_width=8)
        enum_stats = merged.stats()
        # Sample separately with NO memoizer; charge every real oracle call.
        totals = dict.fromkeys(expected_prefix_cost(merged,0,0,0),0)
        def counted_prefix(*args,**kwargs):
            expected = expected_prefix_cost(merged,*args,**kwargs)
            for key,value in expected.items():
                totals[key] = max(totals[key],value) if key=="peak_sector_count" else totals[key]+value
            return original_prefix(*args,**kwargs)
        merged.prefix_vector = counted_prefix
        sampled = merged.sample(np.random.default_rng(7410+supplied.b))
    finally:
        restore_history_classpatch(had, original)
    actual_bytes = len(cache) * supplied.b * np.dtype(np.complex128).itemsize
    return {"target_mass": float(target.sum()), "transition_mass": float(enum.sum()),
            "max_abs_error": float(np.max(np.abs(enum - target))),
            "tv_diagnostic": float(np.sum(np.abs(enum - target)) / 2),
            "cache_entries": len(cache),
            "cache_conservative_bytes": len(cache) * 1024,
            "cache_actual_numpy_bytes": actual_bytes,
            "history_api_present_before": had, "history_sentinel_calls": calls["count"],
            "stats": enum_stats, "sample": sampled,
            "uncached_sample_cost_exact": all(sampled[key]==value for key,value in totals.items())}


def wrong_control_error(cls):
    from experiments.experiment_coherent_route_sampling import direct_joint, enumerate_transitions
    reference = structural_fixture(3, width=4, route_count=4)
    target = direct_joint(reference, max_width=8)
    wrong = cls(structural_fixture(3, width=4, route_count=4))
    enum = enumerate_transitions(wrong, max_width=8)
    return float(np.max(np.abs(enum - target))), float(np.sum(np.abs(enum - target)) / 2)


class WrongPhaseMerged(MergedCoherentPrefixes):
    def prefix_vector(self, *args, **kwargs):
        import lab.merged_prefix as module
        original = module.unit_phase
        module.unit_phase = lambda n, d: original(n, d) * np.exp(0.123j)
        try:
            return super().prefix_vector(*args, **kwargs)
        finally:
            module.unit_phase = original


class WrongBoundaryMerged(MergedCoherentPrefixes):
    def prefix_vector(self, sector, exponent, stop, **kwargs):
        if kwargs.get("boundary", "reflection") == "reflection" and kwargs.get("measured", 0) == 0:
            kwargs["boundary"] = "background"
        return super().prefix_vector(sector, exponent, stop, **kwargs)


def legacy_equivalence():
    from experiments.experiment_coherent_route_sampling import direct_joint, enumerate_transitions
    structural = structural_fixture(3, width=4, route_count=3)
    legacy = CoherentReflectionCircuit(9, 3, 4,
        dict(structural.background.defects), dict(structural.reflections))
    target_s = direct_joint(structural, max_width=8)
    target_l = direct_joint(legacy, max_width=4)
    merged_s = MergedCoherentPrefixes(structural)
    law_s = enumerate_transitions(merged_s, max_width=4)
    law_l = enumerate_transitions(legacy, max_width=4)
    rng_s = np.random.default_rng(991)
    rng_l = np.random.default_rng(991)
    sample_s = merged_s.sample(rng_s)
    sample_l = legacy.sample(rng_l)
    sample_equal = all(sample_s[k] == sample_l[k]
                       for k in ("output", "final_coarse_sector", "final_within_sector"))
    return {"joint_max_abs_error": float(np.max(np.abs(target_s - target_l))),
            "merged_transition_max_abs_error": float(np.max(np.abs(law_s - law_l))),
            "merged_transition_mass_s": float(law_s.sum()),
            "merged_transition_mass_l": float(law_l.sum()),
            "seeded_sample_equal": bool(sample_equal),
            "seeded_sample_s": {k: sample_s[k] for k in ("output", "final_coarse_sector", "final_within_sector")},
            "seeded_sample_l": {k: sample_l[k] for k in ("output", "final_coarse_sector", "final_within_sector")},
            "structural_mass": float(target_s.sum()), "legacy_mass": float(target_l.sum()),
            "legacy_reflection_count": len(legacy.reflections)}


def main():
    report = {"status": "PASS", "preflight": [], "prefixes": [], "laws": [],
              "legacy": {}, "controls": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        b2 = structural_fixture(2)
        b3 = structural_fixture(3)
        report["preflight"] = [reference_guard_info(b2), reference_guard_info(b3)]
        guard((9, 9), label="small dense preflight")
        for name, fixture in (("b2_k9", b2), ("b3_k9", b3)):
            calls, had, original = forbidden_history_classpatch()
            try:
                prefix = prefix_audit(fixture)
            finally:
                restore_history_classpatch(had, original)
            prefix["name"] = name
            prefix["history_api_present_before"] = had
            prefix["history_sentinel_calls"] = calls["count"]
            report["prefixes"].append(prefix)
            p1 &= (prefix["max_abs_error"] < TARGET_TOL
                   and prefix["cache_entries"] <= MAX_CACHE_ENTRIES
                   and prefix["cache_entries"] == prefix["expected_unique_entries"]
                   and prefix["alias_requests"] > 0 and prefix["alias_max_abs_error"] < TARGET_TOL)
            p2 &= not prefix["history_api_present_before"] and prefix["history_sentinel_calls"] == 0
            p2 &= prefix["cost_exact"] and prefix["stats"]["merged_prefix_queries"]==prefix["requests"]
            p2 &= 1<=prefix["stats"]["peak_sector_count"]<=prefix["support_bound"]
            for key in ("history_enumerations","orbit_table_entries","sector_table_entries","output_table_entries"):
                p2 &= prefix["stats"][key]==0
        # Complete law is required for b3 and attempted for b2 within the same
        # explicit opt-in reference/frontier budgets.
        for name, fixture in (("b3_k9", b3), ("b2_k9", b2)):
            law = law_audit(fixture)
            law["name"] = name
            report["laws"].append(law)
            p1 &= law["max_abs_error"] < TARGET_TOL and law["tv_diagnostic"] < TARGET_TOL
            p1 &= abs(law["target_mass"] - 1) < TARGET_TOL and abs(law["transition_mass"] - 1) < TARGET_TOL
            p2 &= not law["history_api_present_before"] and law["history_sentinel_calls"] == 0
            p2 &= law["stats"]["merged_prefix_queries"]==law["cache_entries"]
            p2 &= law["uncached_sample_cost_exact"] and law["sample"]["rejection_proposals"]==0
            p2 &= law["sample"]["merged_prefix_queries"]==law["sample"]["prefix_vector_evaluations"]
            p2 &= 0<=law["sample"]["output"]<(1 << fixture.width)
            p2 &= 0<=law["sample"]["final_coarse_sector"]<fixture.sectors
            p2 &= not law["sample"]["precision_certificate"]
            for key in ("history_enumerations","orbit_table_entries","sector_table_entries","output_table_entries"):
                p2 &= law["stats"][key]==law["sample"][key]==0
        report["legacy"] = legacy_equivalence()
        p2 &= report["legacy"]["joint_max_abs_error"] < TARGET_TOL
        p2 &= report["legacy"]["merged_transition_max_abs_error"] < TARGET_TOL
        p2 &= report["legacy"]["seeded_sample_equal"]

        # Intermediate effects must remain observable with the same routes.
        from experiments.experiment_coherent_route_sampling import direct_joint
        with_background = direct_joint(b3, max_width=8)
        no_background = direct_joint(structural_fixture(3, intermediate=False), max_width=8)
        background_error = float(np.max(np.abs(with_background - no_background)))
        report["background_effect"] = {"max_abs_error": background_error}
        p3 &= background_error > 1e-5

        for label, cls in (("wrong_phase", WrongPhaseMerged), ("wrong_boundary", WrongBoundaryMerged)):
            error, tv = wrong_control_error(cls)
            report["controls"][label] = {"max_abs_error": error, "tv_diagnostic": tv}
        default_failed = legacy_failed = optin_failed = False
        try:
            from experiments.experiment_coherent_route_sampling import direct_joint
            direct_joint(b3)  # default width cap must reject t=8
        except ValueError:
            default_failed = True
        try:
            direct_joint(b3, max_width=9)  # opt-in cap is bounded at 8
        except ValueError:
            optin_failed = True
        try:
            CoherentReflectionCircuit(9, 3, 8,
                dict(b3.background.defects), dict(b3.reflections))
        except ValueError:
            legacy_failed = True
        report["controls"]["guards"] = {"default_width_rejected": default_failed,
                                          "max_width9_rejected": optin_failed,
                                          "legacy_k9_rejected": legacy_failed}
        exp.check("P1", p1, "prefix and complete-law comparisons")
        exp.check("P2", p2, "history absence/cache/legacy equivalence")
        exp.check("P3", p3, f"intermediate background error {background_error:.6g}")
        exp.fail_check("C1", report["controls"]["wrong_phase"]["tv_diagnostic"] > 1e-5,
                       "wrong internal phase changes the law")
        exp.fail_check("C2", report["controls"]["wrong_boundary"]["tv_diagnostic"] > 1e-5,
                       "wrong boundary changes the law")
        exp.fail_check("C3", default_failed and optin_failed and legacy_failed,
                       "unsupported default/opt-in/legacy calls are rejected")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
        exp.fail_check("C2", False, "exception before control")
        exp.fail_check("C3", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "PASS" if all(ok for _, _, ok, _ in exp._results) else "FAIL"
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)],
                    metadata={"max_width_optin": 8, "cache_cap_entries": MAX_CACHE_ENTRIES,
                              "target_tolerance": TARGET_TOL,
                              "native_memory_measured": False,
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
