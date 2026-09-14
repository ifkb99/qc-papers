"""Independent audit of the bounded word-specific E/O support plan.

This is an integer-only verifier.  It compares the production setup against
an independently written recurrence and checks the returned bound against
actual reverse supports.  It does not propagate amplitudes or enumerate
histories at the large-word rows.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  A separately evolved direct support recurrence equals the E/O normal
      form at every reverse prefix, and the selected bound covers it.
  P2  Its E/O final counts, alphabet bound, word envelope, and setup counters
      agree with an independent exact recurrence.
  P3  The plan snapshots/order-validates its input and does not retain a
      mutable caller list.

  C1  The cumulative offset guard rejects a genuinely growing alphabet before
      an oversized integer-set allocation.

The setup counters count candidate offset updates/copies, not native CPU
instructions or resident-memory measurements.  No claim about a full sampler
law is made here.
"""
from __future__ import annotations

import math
import time
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment
from lab.coherent_reverse import SparseCoherentReverse, word_support_plan
from lab.coherent_routes import CoherentReflectionInput
from lab.merged_prefix import MergedCoherentPrefixes


MAX_SET_ENTRIES = 65_536
# Frozen before the strengthened rerun: direct recurrence updates now charge
# each gamma/depth element and retain every prefix for comparison.  The first
# strengthened 300,000 cap was exhausted at 300,001; 500,000 is the new fixed
# ceiling, still an enumerated-operation cap rather than a native CPU/RSS
# claim.
MAX_OPS = 500_000
MAX_TINY_MODULUS = 31

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "word selected bounds cover every reverse-prefix support")
exp.predict("P2", "word E/O counts, bounds and setup counters match independent recurrence")
exp.predict("P3", "word-plan inputs are ordered, validated and snapshotted")
exp.must_fail("C1", "a genuinely growing alphabet trips the cumulative offset guard")


class Budget:
    def __init__(self, limit=MAX_OPS):
        self.limit = int(limit)
        self.operations = 0

    def add(self, amount=1):
        self.operations += int(amount)
        if self.operations > self.limit:
            raise MemoryError(f"word-support verifier operation budget exceeded: "
                              f"{self.operations}>{self.limit}")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"word_envelope_support_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard(entries, label):
    if type(entries) is not int or entries < 0 or entries > MAX_SET_ENTRIES:
        raise MemoryError(f"{label} exceeds {MAX_SET_ENTRIES} integer entries")


def exact_eo_prefixes(modulus, qs, budget):
    """Independent simultaneous E/O recurrence, including every prefix."""
    modulus = int(modulus)
    even, odd = {0}, set()
    rows = [(even, odd)]
    retained_entries = 1
    for q in reversed(tuple(int(q) % modulus for q in qs)):
        old_even, old_odd = even, odd
        old_count = len(old_even) + len(old_odd)
        guard(retained_entries + 2 * old_count,
              "independent E/O retained-prefix allocation")
        budget.add(old_count)
        even = set(old_even)
        even.update((-x - q) % modulus for x in old_odd)
        odd = set(old_odd)
        odd.update((-x - q) % modulus for x in old_even)
        guard(len(even) + len(odd), "independent E/O row")
        rows.append((even, odd))
        retained_entries += len(even) + len(odd)
        guard(retained_entries, "independent E/O retained-prefix total")
    return rows


def exact_reached(modulus, gamma, pair):
    even, odd = pair
    return ({(int(gamma) + x) % modulus for x in even}
            | {(-int(gamma) + x) % modulus for x in odd})


def direct_reverse_prefixes(modulus, gamma, qs, budget):
    """Direct A <- A union f_q(A), independent of E/O offsets."""
    modulus = int(modulus)
    reached = {int(gamma) % modulus}
    prefixes = [reached]
    retained_entries = len(reached)
    guard(retained_entries, "direct initial support")
    for q in reversed(tuple(int(q) % modulus for q in qs)):
        # Guard old plus the largest possible new set before allocating it;
        # retained prefixes are also charged because this reference keeps
        # every depth for the equality audit.
        guard(retained_entries + 2 * len(reached),
              "direct retained-prefix allocation")
        budget.add(len(reached))
        next_reached = set(reached)
        for alpha in reached:
            budget.add(1)
            next_reached.add((-alpha - q) % modulus)
        guard(len(next_reached), "direct support update")
        reached = next_reached
        prefixes.append(reached)
        retained_entries += len(reached)
        guard(retained_entries, "direct retained-prefix total")
    return prefixes


def field(value, name):
    if isinstance(value, Mapping):
        if name in value:
            return value[name]
    if hasattr(value, name):
        return getattr(value, name)
    raise AssertionError(f"word plan result lacks field {name}")


def plan_fields(result):
    names = {
        "exact_steps": "word_setup_steps",
        "copied_offsets": "word_setup_copied_offsets",
        "mapped_offsets": "word_setup_reflected_offsets",
        "peak_live_old_new": "word_setup_peak_live_offset_entries",
        "final_even_count": "word_even_offset_count",
        "final_odd_count": "word_odd_offset_count",
        "alphabet_support": "alphabet_support_bound",
        "word_envelope": "word_support_bound",
        "selected_bound": "support_bound",
        "retained_entries": "word_retained_offset_entries",
    }
    return {name: field(result, api_name) for name, api_name in names.items()}


def alphabet_bound(modulus, qs):
    labels = tuple(dict.fromkeys(int(q) % int(modulus) for q in qs))
    d, k = len(labels), len(qs)
    if d == 0:
        value = 1
    elif d == 1:
        value = 2 if k else 1
    elif d == 2:
        value = 1 + 2 * k
    else:
        value = 2 * sum(2 ** j * math.comb(d - 1, j) * math.comb(k, j)
                         for j in range(min(d - 1, k) + 1))
    return min(int(modulus), 1 << k, value)


def independent_stats(modulus, qs, budget):
    rows = exact_eo_prefixes(modulus, qs, budget)
    even, odd = rows[-1]
    word = min(int(modulus), len(even) + len(odd))
    alpha = alphabet_bound(modulus, qs)
    selected = min(alpha, word)
    old_counts = [len(e) + len(o) for e, o in rows[:-1]]
    # One copy and one mapped update for each live old offset.  The production
    # peak records old plus the two materialized new sets after collisions.
    peak = 1
    for old, new in zip(rows[:-1], rows[1:]):
        peak = max(peak, sum(map(len, old)) + sum(map(len, new)))
    return {
        "rows": rows, "final_even_count": len(even),
        "final_odd_count": len(odd), "word_envelope": word,
        "alphabet_support": alpha, "selected_bound": selected,
        "exact_steps": len(qs), "copied_offsets": sum(old_counts),
        "mapped_offsets": sum(old_counts),
        "peak_live_old_new": peak,
    }


def validate_plan(modulus, qs, gammas, budget, *, all_tiny=False):
    qs = tuple(int(q) for q in qs)
    independent = independent_stats(modulus, qs, budget)
    result = word_support_plan(modulus, qs, max_offset_entries=MAX_SET_ENTRIES)
    actual = plan_fields(result)
    expected_keys = ("exact_steps", "final_even_count", "final_odd_count",
                     "alphabet_support", "word_envelope", "selected_bound")
    counter_mismatches = sum(actual[key] != independent[key] for key in expected_keys)
    counter_mismatches += actual["copied_offsets"] != independent["copied_offsets"]
    counter_mismatches += actual["mapped_offsets"] != independent["mapped_offsets"]
    counter_mismatches += actual["peak_live_old_new"] != independent["peak_live_old_new"]
    counter_mismatches += actual["retained_entries"] != 0
    prefix_failures = 0
    direct_mismatches = 0
    scalar_only = True
    scalar_only = all(isinstance(value, (int, float, str, bool, type(None)))
                      for value in result.values())
    scalar_only &= actual["retained_entries"] == 0
    # All gamma frontiers coexist in this verifier, unlike production.
    # Bound their combined retention and an additional live next-set BEFORE
    # constructing any of them, plus the already-retained E/O references.
    per_gamma = sum(min(modulus, 1 << depth) for depth in range(len(qs)+1))
    eo_entries = sum(len(e)+len(o) for e, o in independent["rows"])
    guard(eo_entries + len(gammas)*per_gamma + 2*min(modulus, 1 << len(qs)),
          "all-gamma reference retention and temporary set")
    direct_by_gamma = {
        gamma: direct_reverse_prefixes(modulus, gamma, qs, budget)
        for gamma in gammas}
    for depth, pair in enumerate(independent["rows"]):
        # The plan is global; every tiny gamma and selected larger gamma must
        # fit the same selected bound at every reverse-prefix depth.
        for gamma in gammas:
            reached = direct_by_gamma[gamma][depth]
            budget.add(sum(map(len, pair)))
            direct_mismatches += reached != exact_reached(modulus, gamma, pair)
            prefix_failures += len(reached) > actual["selected_bound"]
    return {"modulus": int(modulus), "word": qs, "gammas": tuple(gammas),
            "actual": actual, "independent": {
                key: value for key, value in independent.items() if key != "rows"},
            "counter_mismatches": int(counter_mismatches),
            "prefix_bound_failures": int(prefix_failures),
            "direct_eo_mismatches": int(direct_mismatches),
            "scalar_only_plan": bool(scalar_only)}


def input_snapshot_audit():
    labels = [0, 1, 101, 3]
    before = word_support_plan(1009, labels, max_offset_entries=MAX_SET_ENTRIES)
    before_fields = plan_fields(before)
    labels[:] = [999, 998]
    after_fields = plan_fields(before)
    tuple_fields = plan_fields(word_support_plan(1009, (0, 1, 101, 3),
                                                 max_offset_entries=MAX_SET_ENTRIES))
    # Adapters sort insertion keys into original-time order and share one
    # complete chronological setup. Insert the same gates in reverse order.
    reflection_a = {4: (1, 0.3), 1: (0, 0.1), 3: (101, 0.2)}
    reflection_b = {3: (101, 0.2), 4: (1, 0.3), 1: (0, 0.1)}
    input_a = CoherentReflectionInput(6, 3, 4, {}, reflection_a)
    input_b = CoherentReflectionInput(6, 3, 4, {}, reflection_b)
    expected_adapter = plan_fields(word_support_plan(
        input_a.sectors, (0, 101, 1), max_offset_entries=MAX_SET_ENTRIES))
    reverse_a = SparseCoherentReverse(input_a, support_mode="word")
    reverse_b = SparseCoherentReverse(input_b, support_mode="word")
    merged_a = MergedCoherentPrefixes(input_a, support_mode="word")
    merged_b = MergedCoherentPrefixes(input_b, support_mode="word")
    adapters = (reverse_a, reverse_b, merged_a.reverse, merged_b.reverse)
    reverse_snapshot_before = plan_fields(reverse_a.support_plan)
    moved = input_a.reflections.pop(1)
    input_a.reflections[1] = moved
    reverse_snapshot_after = plan_fields(reverse_a.support_plan)
    wide_input = CoherentReflectionInput(
        3 * 1_000_003, 3, 63, {},
        {s: (s % 9, math.pi / 4) for s in range(64)})
    wide_reverse = SparseCoherentReverse(wide_input, support_mode="word")
    wide_merged = MergedCoherentPrefixes(wide_input, support_mode="word")
    try:
        SparseCoherentReverse(wide_input, support_mode="alphabet")
    except MemoryError:
        wide_alphabet_rejected = True
    else:
        wide_alphabet_rejected = False
    return {"list_tuple_equal": before_fields == tuple_fields,
            "mutation_did_not_change_result": before_fields == after_fields,
            "adapter_order_snapshot_equal": all(
                plan_fields(adapter.support_plan) == expected_adapter
                for adapter in adapters)
                and reverse_snapshot_before == reverse_snapshot_after,
            "adapter_support_mode_word": all(
                field(adapter.support_plan, "support_mode") == "word"
                for adapter in adapters),
            "adapter_bound_matches_plan": all(
                adapter.support_bound == expected_adapter["selected_bound"]
                for adapter in adapters),
            "wide_word_input_accepted": (
                wide_reverse.support_bound == 226
                and wide_merged.support_bound == 226),
            "wide_alphabet_cap_rejected": wide_alphabet_rejected,
            "invalid_inputs_rejected": invalid_input_checks()}


def invalid_input_checks():
    checks = []
    for args in ((0, (0,)), (17, "01"), (17, (0, "x"))):
        try:
            word_support_plan(*args, max_offset_entries=MAX_SET_ENTRIES)
        except (TypeError, ValueError):
            checks.append(True)
        else:
            checks.append(False)
    return all(checks)


def growing_cap_control():
    # Distinct labels force the E/O candidate count to grow; the production
    # guard must reject before crossing the frozen 65,536-entry cap.
    # Powers of three deliberately avoid the near-linear cancellation of the
    # consecutive-label word while remaining within the 64-label API limit.
    qs = tuple(3 ** j for j in range(16))
    try:
        word_support_plan(1_000_000_007, qs,
                          max_offset_entries=MAX_SET_ENTRIES)
    except MemoryError as exc:
        return {"rejected": True, "word_length": len(qs),
                "exception_type": type(exc).__name__,
                "message": str(exc)}
    except ValueError as exc:
        return {"rejected": False, "word_length": len(qs),
                "exception_type": type(exc).__name__,
                "message": str(exc)}
    return {"rejected": False, "word_length": len(qs)}


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    return value


def main():
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "operations": 0}
    p1 = p2 = p3 = True
    try:
        budget = Budget()
        # All gamma values on tiny sectors; selected labels on wider sectors.
        specs = []
        for modulus in range(1, MAX_TINY_MODULUS + 1):
            specs.append((modulus, tuple(range(modulus)), True))
        specs.extend(((101, (0, 1, 50), False),
                      (1009, (0, 1, 17, 503), False),
                      (1_000_003, (0, 1, 500_001), False)))
        for modulus, gammas, all_tiny in specs:
            qs = tuple((j * j + 3 * j + 1) % modulus for j in range(8))
            report["rows"].append(validate_plan(modulus, qs, gammas, budget,
                                                 all_tiny=all_tiny))
        p1 = all(row["prefix_bound_failures"] == 0
                 and row["direct_eo_mismatches"] == 0
                 for row in report["rows"])
        p2 = all(row["counter_mismatches"] == 0
                 and row["scalar_only_plan"]
                 for row in report["rows"])
        report["snapshot"] = input_snapshot_audit()
        p3 = all(report["snapshot"].values())
        report["cap_control"] = growing_cap_control()
        c1 = report["cap_control"]["rejected"]
        report["operations"] = budget.operations
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "all tested reverse-prefix supports fit the selected word bound")
        exp.check("P2", p2, "independent E/O counts and setup counters match")
        exp.check("P3", p3, "input ordering, snapshot and validation checks pass")
        exp.fail_check("C1", c1, "growing alphabet is rejected by the offset guard")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "max_set_entries": MAX_SET_ENTRIES,
        "max_operations": MAX_OPS,
        "operation_scope": "independent candidate/update operations, not native CPU count",
        "no_history_enumeration": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
