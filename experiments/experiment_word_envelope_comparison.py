"""Bounded comparison of alphabet and word-specific global envelopes.

TODO32/C73 supplies a word-specific E/O offset envelope.  This diagnostic
asks only whether that cheaper envelope is integrated and charged correctly
on the separated three-label input; it makes no timing-ranking or
finite-TV claim.  The C71 adjacent-label approximation is not applicable to
this widely separated fixture, and failure of its sufficient bound would not
be a lower bound on any method.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The independent chronological-word plan agrees with production's word
      plan, and its global bound is never larger than the alphabet bound.
  P2  Successful reverse/prefix samples stay within the selected global
      envelope, without history expansion, and method-specific counters equal
      the instrumented attempt/query calls.
  P3  Setup, integer offset work, complex slots, and returned-sample work are
      reported separately for every (k, method, seed) row.
  P4  The word setup is charged once per word-mode adapter and no history or
      rejection work is assigned to the merged-prefix method.

  C1  Charging only one reverse attempt must undercharge a returned sample
      that actually required retries.
  C2  Omitting nonzero word-plan setup must undercharge every nonempty word.

The wide sweep is deliberately gated by the caller: the core gate and the
independent tiny-law audit must pass before running this module.  Dense orbit
or output arrays are never allocated.  Integer-entry counts are declared
arithmetic payload, not native RSS.
"""
from __future__ import annotations

import math
import json
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit, CoherentReflectionInput
from lab.coherent_reverse import SparseCoherentReverse
from lab.merged_prefix import MergedCoherentPrefixes


PERIOD = 3 * (2**40 - 1)
BLOCK = 3
WIDTH = 16
MODULUS = PERIOD // BLOCK
LABELS = (0, 679_535_556_937, 314_159_265_359)
THETA = np.pi / 4
SEEDS = (7430, 7431, 7432)
METHODS = ("reverse_alphabet", "reverse_word",
           "prefix_alphabet", "prefix_word")
LEGACY_ANCHORS = (0, 4, 8)
MAX_PROPOSALS = 8192
ROW_SECONDS = 60.0
MAX_OFFSET_ENTRIES = 65_536
MAX_BYTES = 16 << 20
SWEEP_SECONDS = 600.0  # frozen main-run budget, checked between rows

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=True)
exp.predict("P1", "independent word E/O plans agree with production and improve or match alphabet bounds")
exp.predict("P2", "selected envelopes contain returned support and instrumented counters match charged work")
exp.predict("P3", "setup, integer, complex-slot, and method-specific sample work remain separate per row")
exp.predict("P4", "word setup is charged and merged prefixes receive no rejection/history work")
exp.must_fail("C1", "a one-attempt reverse charge must undercount a sample requiring retries")
exp.must_fail("C2", "dropping nonzero word-plan setup must undercount a nonempty word")


def stamp():
    now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"word_envelope_comparison_{now}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"word_envelope_comparison_{now}_{serial}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


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
    if isinstance(value, (float,)) and not math.isfinite(value):
        return str(value)
    return value


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2),
                    np.exp(1j * theta / 2)]).astype(complex)


def embed(block):
    result = np.eye(BLOCK, dtype=complex)
    result[:2, :2] = block
    return result


def defects():
    result = {0: embed(rx(np.pi / 4))}
    for insertion in range(1, WIDTH + 1):
        result[insertion] = embed(
            rx(np.pi / 7) if insertion % 2 else rz(np.pi / 5))
    return result


def word_for(k):
    return tuple(LABELS[i % len(LABELS)] for i in range(k))


def independent_word_plan(labels):
    """Independent E/O recurrence; returns counts, never a sector table."""
    even, odd = {0}, set()
    copied = reflected = 0
    peak_live = 1
    retained = 1
    for q in reversed(tuple(int(q) % MODULUS for q in labels)):
        old_size = len(even) + len(odd)
        # Existing sets plus two candidate images are the only live integer
        # sets.  Check before creating the next sets.
        if 3 * old_size > MAX_OFFSET_ENTRIES:
            raise MemoryError("independent word live-offset cap exceeded")
        copied += old_size
        reflected += old_size
        next_even = set(even)
        next_odd = set(odd)
        next_even.update((-value - q) % MODULUS for value in odd)
        next_odd.update((-value - q) % MODULUS for value in even)
        even, odd = next_even, next_odd
        retained += len(even) + len(odd)
        peak_live = max(peak_live, old_size + len(even) + len(odd))
        if peak_live > MAX_OFFSET_ENTRIES:
            raise MemoryError("independent word peak-offset cap exceeded")
    distinct = len(set(labels))
    if distinct == 0:
        alphabet_bound = 1
    elif distinct == 1:
        alphabet_bound = min(MODULUS, 1 << len(labels),
                             1 if not labels else 2)
    elif distinct == 2:
        alphabet_bound = min(MODULUS, 1 << len(labels), 1 + 2 * len(labels))
    else:
        alphabet_bound = min(MODULUS, 1 << len(labels),
                             2 * sum((2**j) * math.comb(distinct - 1, j)
                                      * math.comb(len(labels), j)
                                      for j in range(min(distinct - 1,
                                                         len(labels)) + 1)))
    word_bound = min(MODULUS, len(even) + len(odd))
    return dict(support_mode="word", support_bound=min(alphabet_bound, word_bound),
                alphabet_support_bound=alphabet_bound,
                word_support_bound=word_bound,
                word_even_offset_count=len(even), word_odd_offset_count=len(odd),
                word_setup_steps=len(labels),
                word_setup_copied_offsets=copied,
                word_setup_reflected_offsets=reflected,
                word_setup_peak_live_offset_entries=peak_live,
                word_setup_max_offset_entries=MAX_OFFSET_ENTRIES,
                # The helper releases its temporary sets before returning;
                # this is not a retained allocation or a cumulative charge.
                word_retained_offset_entries=0,
                independent_word_validation_retained_offset_entries=0,
                independent_word_integer_work=copied + reflected)


def fixture(k, *, legacy=False):
    reflections = {i: (LABELS[i % len(LABELS)], THETA) for i in range(k)}
    cls = CoherentReflectionCircuit if legacy else CoherentReflectionInput
    return cls(PERIOD, BLOCK, WIDTH, defects(), reflections)


def method_mode(method):
    return "word" if method.endswith("_word") else "alphabet"


def make_worker(k, method):
    if method.startswith("reverse"):
        return SparseCoherentReverse(fixture(k), support_mode=method_mode(method))
    return MergedCoherentPrefixes(fixture(k), support_mode=method_mode(method))


def worker_stats(worker):
    if hasattr(worker, "stats"):
        return dict(worker.stats())
    return dict(worker._stats())


def install_clock(worker, method, started):
    calls = {"attempts": 0, "queries": 0, "history_items": 0,
             "completed": {}, "first_completed": None, "completed_calls": 0}

    def record_completed(values):
        if calls["first_completed"] is None:
            calls["first_completed"] = dict(values)
        calls["completed_calls"] += 1
        for key, value in values.items():
            calls["completed"][key] = calls["completed"].get(key, 0) + value

    def check_clock():
        if time.perf_counter() - started > ROW_SECONDS:
            raise RuntimeError(f"row exceeded soft {ROW_SECONDS:.0f}s limit")

    if method.startswith("reverse"):
        original = worker.attempt

        def timed_attempt(*args, **kwargs):
            check_clock()
            calls["attempts"] += 1
            result = original(*args, **kwargs)
            record_completed({key: int(result.get(key, 0)) for key in (
                "reverse_steps", "work_matvecs", "reflection_vector_contributions",
                "qft_matrix_pair_constructions")})
            check_clock()
            return result

        worker.attempt = timed_attempt
    elif method == "legacy_history":
        original_prefix = worker.prefix_vector
        original_histories = worker._histories

        def timed_query(*args, **kwargs):
            check_clock()
            calls["queries"] += 1
            result = original_prefix(*args, **kwargs)
            record_completed({"prefix_vector_evaluations": 1})
            check_clock()
            return result

        def timed_histories(*args, **kwargs):
            check_clock()
            count = 0
            try:
                for item in original_histories(*args, **kwargs):
                    count += 1
                    yield item
            finally:
                calls["history_items"] += count
                check_clock()

        worker.prefix_vector = timed_query
        worker._histories = timed_histories
    else:
        original = worker.prefix_vector

        def timed_query(*args, **kwargs):
            check_clock()
            calls["queries"] += 1
            result = original(*args, **kwargs)
            last = getattr(worker, "last_prefix_stats", None) or {}
            record_completed({key: int(last.get(key, 0)) for key in (
                "merged_prefix_queries", "prefix_control_steps",
                "background_block_products", "control_block_products",
                "fixed_control_block_scales", "qft_matrix_pair_constructions",
                "reflection_block_contributions")})
            check_clock()
            return result

        worker.prefix_vector = timed_query
    return calls


def completed_sum(calls, key):
    return int(calls["completed"].get(key, 0))


def first_completed_sum(calls, keys):
    if calls["first_completed"] is None:
        return 0
    return sum(int(calls["first_completed"].get(key, 0)) for key in keys)


def work_fields(method, result, calls):
    if method.startswith("reverse"):
        keys = ("reverse_steps", "work_matvecs",
                "reflection_vector_contributions",
                "qft_matrix_pair_constructions")
        actual_attempts = int(result.get("rejection_proposals", calls["attempts"]))
        values = {key: int(result.get(key, completed_sum(calls, key)))
                  for key in keys}
        return dict(work_kind="rejection_attempts_and_reverse_steps",
                    actual_attempts=actual_attempts,
                    instrumented_attempts=int(calls["attempts"]),
                    instrumented_queries=0,
                    **values,
                    actual_reverse_work=sum(values.values()),
                    wrong_one_attempt_work=first_completed_sum(calls, keys),
                    prefix_vector_evaluations=0,
                    merged_prefix_queries=0,
                    background_block_products=0,
                    control_block_products=0,
                    completed_work={key: completed_sum(calls, key) for key in keys})
    if method == "legacy_history":
        return dict(work_kind="legacy_history_prefix_and_history_items",
                    actual_attempts=0, instrumented_attempts=0,
                    instrumented_queries=int(calls["queries"]),
                    history_items=int(calls["history_items"]),
                    reverse_steps=0, work_matvecs=0,
                    reflection_vector_contributions=0,
                    qft_matrix_pair_constructions=0,
                    actual_reverse_work=0, wrong_one_attempt_work=0,
                    prefix_vector_evaluations=int(
                        result.get("prefix_vector_evaluations", calls["queries"])),
                    merged_prefix_queries=0, background_block_products=0,
                    control_block_products=0,
                    completed_work={"prefix_vector_evaluations": calls["queries"],
                                    "history_items": calls["history_items"]})
    keys = ("merged_prefix_queries", "prefix_control_steps",
            "background_block_products", "control_block_products",
            "fixed_control_block_scales", "qft_matrix_pair_constructions",
            "reflection_block_contributions")
    values = {key: int(result.get(key, completed_sum(calls, key))) for key in keys}
    return dict(work_kind="merged_prefix_queries_and_block_products",
                actual_attempts=0,
                instrumented_attempts=0,
                instrumented_queries=int(calls["queries"]),
                reverse_steps=0,
                work_matvecs=0,
                reflection_vector_contributions=0,
                **values,
                actual_reverse_work=0,
                wrong_one_attempt_work=0,
                prefix_vector_evaluations=int(result.get(
                    "prefix_vector_evaluations", calls["queries"])),
                completed_work={key: completed_sum(calls, key) for key in keys})


def run_row(k, method, seed, *, phase):
    started = time.perf_counter()
    mode = method_mode(method) if method in METHODS else "history"
    row = dict(phase=phase, k=int(k), method=method, support_mode=mode,
               seed=int(seed), period=PERIOD, block_size=BLOCK, width=WIDTH,
               labels=list(word_for(k)), max_proposals=MAX_PROPOSALS,
               output_count=1 << WIDTH, coarse_sector_count=MODULUS,
               fine_work_dimension=BLOCK,
               target_law_cells=MODULUS * (1 << WIDTH),
               payload_cap_bytes=MAX_BYTES)
    plan = None
    worker = None
    before = {}
    after = {}
    calls = {"attempts": 0, "queries": 0, "history_items": 0,
             "completed": {}, "first_completed": None, "completed_calls": 0}
    setup_started = setup_done = sample_started = sample_done = None
    result = {}
    plan_calls = {"count": 0}
    try:
        validation_started = time.perf_counter()
        plan = independent_word_plan(word_for(k)) if method in METHODS else None
        row["independent_validation_seconds"] = time.perf_counter() - validation_started

        setup_started = time.perf_counter()
        if method in METHODS:
            import lab.coherent_reverse as reverse_module
            original_plan = reverse_module.word_support_plan

            def counted_plan(*args, **kwargs):
                plan_calls["count"] += 1
                return original_plan(*args, **kwargs)

            if mode == "word":
                reverse_module.word_support_plan = counted_plan
            try:
                worker = make_worker(k, method)
            finally:
                if mode == "word":
                    reverse_module.word_support_plan = original_plan
        else:
            worker = CoherentReflectionCircuit(
                PERIOD, BLOCK, WIDTH, defects(),
                {i: (LABELS[i % len(LABELS)], THETA) for i in range(k)})
        before = worker_stats(worker)
        setup_done = time.perf_counter()
        calls = install_clock(worker, method, started)
        rng = np.random.default_rng(seed)
        sample_started = time.perf_counter()
        def no_rebuild(*args, **kwargs):
            raise AssertionError("word plan rebuilt during sampling")
        if mode == "word":
            reverse_module.word_support_plan = no_rebuild
        try:
            if method.startswith("reverse"):
                result = worker.sample(rng, max_proposals=MAX_PROPOSALS)
            else:
                result = worker.sample(rng)
        finally:
            if mode == "word":
                reverse_module.word_support_plan = original_plan
        sample_done = time.perf_counter()
        after = worker_stats(worker)
        fields = work_fields(method, result, calls)
        production_plan = {key: after.get(key) for key in (
            "support_bound", "alphabet_support_bound", "word_support_bound",
            "word_setup_steps", "word_setup_copied_offsets",
            "word_setup_reflected_offsets", "word_setup_peak_live_offset_entries",
            "word_setup_max_offset_entries", "word_retained_offset_entries",
            "word_even_offset_count", "word_odd_offset_count")}
        if plan is not None and mode == "word":
            plan_match = all(production_plan[key] == plan[key] for key in (
                "support_bound", "alphabet_support_bound", "word_support_bound",
                "word_setup_steps", "word_setup_copied_offsets",
                "word_setup_reflected_offsets", "word_setup_peak_live_offset_entries",
                "word_setup_max_offset_entries", "word_retained_offset_entries",
                "word_even_offset_count", "word_odd_offset_count"))
            integer_setup_work = int(plan["independent_word_integer_work"])
        else:
            plan_match = True
            integer_setup_work = 0
        plan_once = (plan_calls["count"] == 1) if mode == "word" else None
        payload_slots = int(after.get(
            "working_complex_coordinate_slots_upper_bound", 0)) + int(
                after.get("supplied_gate_complex_entries", 0))
        row.update(fields, setup_seconds=setup_done-setup_started,
                   sample_seconds=sample_done-setup_done,
                   adapter_seconds=setup_done-setup_started,
                   independent_word_plan=plan,
                   production_support_plan=production_plan,
                   word_plan_matches=bool(plan_match),
                   word_setup_matches=bool(plan_match),
                   word_support_plan_calls=int(plan_calls["count"]),
                   word_support_plan_exact_once=plan_once,
                   word_setup_released=(
                       plan is None or plan.get("word_retained_offset_entries", 0) == 0),
                   charged_integer_setup_work=integer_setup_work,
                   wrong_omitted_word_setup_work=0,
                   working_complex_coordinate_slots_upper_bound=int(
                       after.get("working_complex_coordinate_slots_upper_bound", 0)),
                   supplied_gate_complex_entries=int(
                       after.get("supplied_gate_complex_entries", 0)),
                   payload_upper_bound_bytes=16*payload_slots,
                   payload_upper_bound_positive=16*payload_slots > 0,
                   payload_within_cap=0 < 16*payload_slots <= MAX_BYTES,
                   verifier_storage="constant-size running counters; no attempt log",
                   output=int(result["output"]),
                   final_coarse_sector=int(result["final_coarse_sector"]),
                   peak_sector_count=int(result.get("peak_sector_count", 0)),
                   support_bound=int(result.get("support_bound", after.get("support_bound", 0))),
                   history_enumerations=(int(calls["history_items"]) if method=="legacy_history"
                                         else int(result.get("history_enumerations", 0))),
                   rejection_proposals=int(result.get("rejection_proposals", 0)),
                   setup_stats=before, final_stats=after, successful=True,
                   counter_match=(
                       fields["instrumented_attempts"] == fields["actual_attempts"]
                       and all(fields["completed_work"].get(key, 0) == fields.get(key, 0)
                               for key in ("reverse_steps", "work_matvecs",
                                           "reflection_vector_contributions",
                                           "qft_matrix_pair_constructions"))
                       if method.startswith("reverse") else
                       fields["instrumented_queries"] == fields["prefix_vector_evaluations"]
                       and all(fields["completed_work"].get(key, 0) == fields.get(key, 0)
                               for key in ("merged_prefix_queries", "prefix_control_steps",
                                           "background_block_products", "control_block_products",
                                           "fixed_control_block_scales",
                                           "qft_matrix_pair_constructions",
                                           "reflection_block_contributions"))),
                   peak_within_bound=(
                       int(result.get("peak_sector_count", 0)) <=
                       int(result.get("support_bound", after.get("support_bound", 0)))))
        return row
    except Exception as exc:
        if worker is not None and not after:
            try:
                after = worker_stats(worker)
            except Exception:
                after = {}
        fields = work_fields(method, result, calls)
        row.update(fields,
                   setup_seconds=(None if setup_started is None else
                                  (setup_done or time.perf_counter())-setup_started),
                   sample_seconds=(None if sample_started is None else
                                   (sample_done or time.perf_counter())-sample_started),
                   adapter_seconds=(None if setup_started is None else
                                    (setup_done or time.perf_counter())-setup_started),
                   independent_word_plan=plan, setup_stats=before,
                   final_stats=after,
                   word_support_plan_calls=int(plan_calls["count"]),
                   word_support_plan_exact_once=None,
                   successful=False,
                   budget_censored=(type(exc).__name__ == "RuntimeError"
                                    and ("soft" in str(exc) or "cap exhausted" in str(exc))),
                   completed_calls=calls["completed_calls"],
                   error_type=type(exc).__name__, error=str(exc),
                   traceback=traceback.format_exc())
        return row


def main():
    rows = []
    checkpoint = stamp().with_suffix(".partial.json")
    def record_row(k, method, seed, phase):
        print(f"starting {phase} k={k} method={method} seed={seed}", flush=True)
        rows.append(run_row(k, method, seed, phase=phase))
        # Generated diagnostics, outside measured setup/sample intervals.
        # Preserve completed rows even if the next native call aborts Python.
        checkpoint.write_text(json.dumps(dict(
            status="INCOMPLETE_CHECKPOINT_NOT_A_VERDICT", python=sys.version,
            numpy=np.__version__, rows=json_safe(rows)), indent=2, allow_nan=False)+"\n")
    preflight_only = os.environ.get("WORD_ENVELOPE_PREFLIGHT_ONLY") == "1"
    expected_sweep_rows = len(METHODS) * (WIDTH + 1) * len(SEEDS)
    expected_word_sweep_rows = 2 * (WIDTH + 1) * len(SEEDS)
    metadata = dict(period=PERIOD, modulus=MODULUS, block_size=BLOCK,
                    width=WIDTH, labels=list(LABELS), theta="pi/4",
                    seeds=list(SEEDS), methods=list(METHODS),
                    max_proposals=MAX_PROPOSALS, row_seconds=ROW_SECONDS,
                    max_offset_entries=MAX_OFFSET_ENTRIES, max_bytes=MAX_BYTES,
                    timing_claim="diagnostic only; no method ranking prediction",
                    adjacent_approximation=(
                        "C71 neighboring-label sufficient bound is not applicable "
                        "to this separated-label fixture; failure is not a lower bound"),
                    preflight_status="not_run_until_tiny_law_audit_passes",
                    expected_sweep_rows=expected_sweep_rows,
                    expected_word_sweep_rows=expected_word_sweep_rows,
                    preflight_only=preflight_only, sweep_seconds=SWEEP_SECONDS,
                    censoring_policy="retain capped failures without fallback; do not pool as successes")
    try:
        # This function is intentionally not called until the parent agent has
        # approved the core/tiny-law gate.  When called, preserve every failed
        # preflight row before deciding whether to sweep.
        for method in METHODS:
            record_row(WIDTH, method, SEEDS[0], "preflight")
        preflight = list(rows)
        preflight_ok = all(row.get("successful") and row.get("counter_match")
                           and row.get("peak_within_bound") and row.get("payload_within_cap")
                           and row.get("word_plan_matches") for row in preflight)
        metadata["preflight_status"] = "passed" if preflight_ok else "failed"
        metadata["preflight_rows"] = len(preflight)
        if preflight_only:
            metadata["scientific_verdict"] = "preflight_only; full sweep not executed"
            payload = dict(name=exp.name, ok=bool(preflight_ok),
                           predictions=exp._predictions, controls=exp._controls,
                           checks=[], warnings=["full sweep intentionally omitted"],
                           rows=json_safe(rows), metadata=json_safe(metadata))
            path = stamp()
            path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
            if not preflight_ok:
                sys.exit(1)
            return
        if preflight_ok and time.perf_counter() - exp._started < 20.0:
            # Optional legacy anchors are decided before the full sweep, so
            # their inclusion is genuinely conditional on the preflight
            # budget rather than on the already-completed benchmark.
            for k in LEGACY_ANCHORS:
                for seed in SEEDS:
                    record_row(k, "legacy_history", seed, "legacy_anchor")
        if preflight_ok:
            sweep_started = time.perf_counter()
            for k in range(WIDTH + 1):
                for method in METHODS:
                    for seed in SEEDS:
                        if time.perf_counter()-sweep_started > SWEEP_SECONDS:
                            raise RuntimeError("frozen sweep budget exhausted between rows")
                        record_row(k, method, seed, "sweep")
                    print(f"completed k={k} method={method}", flush=True)
        metadata["sweep_rows"] = sum(row.get("phase") == "sweep" for row in rows)
        metadata["legacy_rows"] = sum(row.get("phase") == "legacy_anchor" for row in rows)
    except Exception as exc:
        metadata["harness_exception"] = dict(type=type(exc).__name__,
                                               error=str(exc),
                                               traceback=traceback.format_exc())

    successful = [row for row in rows if row.get("successful")]
    primary_successful = [row for row in successful if row["method"] in METHODS]
    sweep_rows = [row for row in rows if row.get("phase") == "sweep"
                  and row.get("method") in METHODS]
    successful_sweep = [row for row in sweep_rows if row.get("successful")]
    # A capped rejection run returns no sample. It is a preserved diagnostic
    # outcome, not a wrong output law and not a successful timing datum.
    accounted_sweep = (len(sweep_rows) == expected_sweep_rows and all(
        row.get("successful") or (row.get("budget_censored")
            and row.get("setup_stats") and row.get("completed_calls", 0)>0
            and row.get("traceback")) for row in sweep_rows))
    word_rows = [row for row in primary_successful if row.get("support_mode") == "word"]
    reverse_rows = [row for row in primary_successful if row["method"].startswith("reverse")]
    setup_omission_rows = [row for row in word_rows
                           if row.get("k", 0) > 0
                           and row.get("charged_integer_setup_work", 0) > 0]
    retry_rows = [row for row in reverse_rows if row.get("actual_attempts", 0) > 1]
    exp.check("P1", len(sweep_rows) == expected_sweep_rows
              and len([row for row in successful_sweep
                       if row.get("support_mode") == "word"]) == expected_word_sweep_rows
              and all(row.get("word_plan_matches")
                                             and row["support_bound"] <= row["production_support_plan"].get("alphabet_support_bound", row["support_bound"])
                                             and row.get("word_support_plan_exact_once")
                                             for row in successful_sweep
                                             if row.get("support_mode") == "word"),
              f"word plan rows={len(word_rows)}")
    exp.check("P2", accounted_sweep and bool(successful_sweep)
              and all(row.get("counter_match")
                                                       and row.get("peak_within_bound")
                                                       for row in successful_sweep),
              f"successful sweep rows={len(successful_sweep)}/{expected_sweep_rows}")
    exp.check("P3", accounted_sweep and bool(successful_sweep)
              and all(row.get("setup_seconds", 0) >= 0
                                               and row.get("working_complex_coordinate_slots_upper_bound", 0) >= 0
                                               and row.get("charged_integer_setup_work", 0) >= 0
                                               and row.get("payload_upper_bound_positive")
                                               and row.get("payload_within_cap")
                                               and row.get("word_setup_released", True)
                                               for row in successful_sweep),
              f"resource rows={len(successful_sweep)}/{expected_sweep_rows}")
    exp.check("P4", accounted_sweep and bool(successful_sweep)
              and bool(setup_omission_rows) and all(
        row.get("word_setup_matches", True)
        and row.get("history_enumerations", 0) == 0
        and (not row["method"].startswith("prefix")
             or row.get("rejection_proposals", 0) == 0)
        for row in successful_sweep),
              f"word setup rows={len(setup_omission_rows)}")
    exp.fail_check("C1", accounted_sweep
        and bool(retry_rows) and any(
        row.get("actual_reverse_work", 0) > row.get("wrong_one_attempt_work", 0)
        for row in retry_rows),
        f"retry rows={len(retry_rows)}; one-attempt reverse charge is deliberately compared")
    exp.fail_check("C2", accounted_sweep
        and bool(setup_omission_rows) and any(
        row.get("charged_integer_setup_work", 0) > row.get("wrong_omitted_word_setup_work", 0)
        for row in setup_omission_rows),
        f"nonempty word rows with charged setup={len(setup_omission_rows)}")
    metadata["row_counts"] = {"all": len(rows), "successful": len(successful),
                                "failed": len(rows)-len(successful)}
    metadata["sweep_complete"] = (len(sweep_rows) == expected_sweep_rows
                                   and len(successful_sweep) == expected_sweep_rows)
    metadata["all_planned_sweep_rows_accounted"] = bool(accounted_sweep)
    metadata["budget_censored_rows"] = sum(
        bool(row.get("budget_censored")) for row in rows)
    metadata["unexpected_failure_rows"] = sum(
        bool(row.get("successful") is False and not row.get("budget_censored"))
        for row in rows)
    metadata["controls"] = dict(retry_rows=len(retry_rows),
                                 word_setup_rows=len(setup_omission_rows))
    metadata["incremental_checkpoint"] = str(checkpoint)
    metadata["numpy"] = np.__version__
    path = stamp()
    ok = exp.finish(report_path=path, rows=json_safe(rows), metadata=json_safe(metadata))
    print(f"report: {path}", flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
