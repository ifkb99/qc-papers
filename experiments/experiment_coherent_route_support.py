"""Exact finite-time support audit for two-label coherent reflections.

For chronological reflection labels q_1,...,q_k, a reverse final-sector
support is updated in the actual order

    R <- R union { -x-q_i mod M : x in R },  i=k,...,1.

Two fixed labels therefore give at most 1+2k reached sectors, irrespective of
whether their difference is invertible.  This experiment audits that statement
and its peak-prefix version without constructing a history table.  It also
keeps the generic-label ordering control: the supplied chronological word
(0,1,0,3,0,9,0,27) has backward support 108 at M=1009,gamma=17, while
applying it in listed order has support 88.  Neither belongs under the
two-label bound.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every binary chronological q word through k=8 obeys both leaf and
      every-prefix bounds min(M,1+2d) for M=7,16,101,1009 at the requested
      starting sectors.
  P2  The noninvertible pair {2,6} modulo 16 and coinciding/fixed-point pairs
      obey the same support bound; no affine relabeling assumption is needed.
  P3  The generic backward control has support 108, contains all 16 positive
      subset-sum translations of gamma=17, and differs from the direct-listed
      order support 88.
  C1  Applying the two-label bound to the generic eight-label word must fail.
  C2  Treating chronological and reversed insertion order as the same support
      must fail on the generic control.

Only bounded integer sets are used transiently. The report stores counts and
small control sets, not all route histories or a sector table.
"""
from __future__ import annotations

import itertools
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment


MAX_BYTES = 16 << 20
MAX_K = 8

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "binary two-label words satisfy leaf and peak support bounds")
exp.predict("P2", "noninvertible and coinciding/fixed-point pairs satisfy the same bound")
exp.predict("P3", "generic backward support is 108 with positive subset translations and differs from order-forward support 88")
exp.must_fail("C1", "the two-label bound incorrectly applied to the generic eight-label word")
exp.must_fail("C2", "chronological and reversed insertion order produce the same generic support")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"coherent_route_support_{stamp}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def guard_entries(entries, bytes_per_entry=256):
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("support diagnostic payload exceeds 16 MiB")


def reverse_prefixes(M, gamma, chronological):
    reached = {int(gamma) % M}
    prefixes = [reached]
    for q in reversed(tuple(chronological)):
        reached = reached | {(-x - int(q)) % M for x in reached}
        prefixes.append(reached)
    return prefixes


def binary_pair_sweep(M, qpair, gammas):
    rows = []
    for k in range(MAX_K + 1):
        words = 1 << k
        max_leaf = max_peak = 0
        violating_words = 0
        for word in itertools.product(qpair, repeat=k):
            for gamma in gammas:
                prefixes = reverse_prefixes(M, gamma, word)
                sizes = [len(support) for support in prefixes]
                max_leaf = max(max_leaf, sizes[-1])
                max_peak = max(max_peak, max(sizes))
                if any(size > min(M, 1 + 2 * depth)
                       for depth, size in enumerate(sizes)):
                    violating_words += 1
        rows.append({"k": k, "word_count": words, "gamma_count": len(tuple(gammas)),
                     "max_leaf": max_leaf, "max_peak": max_peak,
                     "bound": min(M, 1 + 2 * k), "violating_words": violating_words})
    return rows


def fixed_points(M, q):
    return [x for x in range(M) if (-x - q) % M == x]


def generic_control():
    M, gamma = 1009, 17
    word = (0, 1, 0, 3, 0, 9, 0, 27)
    backward = reverse_prefixes(M, gamma, word)[-1]
    forward = {gamma}
    for q in word:
        forward = forward | {(-x - q) % M for x in forward}
    subset_sums = {sum(q for bit, q in zip(mask, word) if bit) % M
                   for mask in itertools.product((0, 1), repeat=len(word))}
    positive_translations = {(gamma + shift) % M for shift in subset_sums}
    return {
        "M": M, "gamma": gamma, "word": word,
        "backward_support_size": len(backward),
        "forward_listed_support_size": len(forward),
        "backward_support_min": min(backward),
        "backward_support_max": max(backward),
        "positive_subset_translation_count": len(positive_translations),
        "positive_subset_translations_contained": positive_translations <= backward,
        "two_label_bound_at_k8": min(M, 1 + 2 * len(word)),
        # Small exact sets are retained because they are the control itself.
        "positive_subset_translations": sorted(positive_translations),
        "ordering_sets_differ": backward != forward,
    }


def main():
    report = {"status": "PASS", "sweeps": [], "special_cases": [], "generic": {}}
    p1 = p2 = p3 = True
    started = time.perf_counter()
    try:
        # One word's full prefix-set list is retained: at most 2^(k+1)-1
        # labels across that list even for unrestricted reflections. Include
        # simultaneous union temporaries and compact report summaries, not
        # every discarded word/starting-sector calculation in the sweep.
        guard_entries(3 * (1 << (MAX_K+1)) + 7 * (MAX_K + 1) + 2048)
        sweep_specs = [
            ("M7_pair01_allgamma", 7, (0, 1), tuple(range(7))),
            ("M16_pair01_allgamma", 16, (0, 1), tuple(range(16))),
            ("M101_pair01_selected", 101, (0, 1), (0, 1, 17)),
            ("M1009_pair01_selected", 1009, (0, 1), (0, 1, 17)),
            ("M16_noninvertible_pair26", 16, (2, 6), tuple(range(16))),
        ]
        for name, M, qpair, gammas in sweep_specs:
            rows = binary_pair_sweep(M, qpair, gammas)
            report["sweeps"].append({"name": name, "M": M, "qpair": qpair,
                                     "gammas": gammas, "rows": rows})
            p1 &= all(row["violating_words"] == 0 for row in rows)
            if name == "M16_noninvertible_pair26":
                p2 &= all(row["violating_words"] == 0 for row in rows)

        special_specs = [("M16_coinciding_q2_fixed", 16, (2, 2), (0, 7, 15)),
                         ("M7_coinciding_q0_fixed", 7, (0, 0), (0, 1, 6))]
        for name, M, qpair, gammas in special_specs:
            rows = binary_pair_sweep(M, qpair, gammas)
            fixed = fixed_points(M, qpair[0])
            report["special_cases"].append({"name": name, "M": M,
                                              "qpair": qpair, "gammas": gammas,
                                              "fixed_points": fixed, "rows": rows})
            p2 &= all(row["violating_words"] == 0 for row in rows)
            p2 &= all(row["max_leaf"] <= 2 for row in rows)
            p2 &= fixed == ([7,15] if M==16 else [0])
            p2 &= all(all(len(support)==1 for support in
                          reverse_prefixes(M,gamma,(qpair[0],)*MAX_K))
                      for gamma in fixed)

        generic = generic_control()
        report["generic"] = generic
        p3 &= generic["backward_support_size"] == 108
        p3 &= generic["positive_subset_translation_count"] == 16
        p3 &= generic["positive_subset_translations_contained"]
        p3 &= generic["forward_listed_support_size"] == 88
        p3 &= generic["ordering_sets_differ"]

        exp.check("P1", p1, "all binary pair sweep prefixes stay within min(M,1+2d)")
        exp.check("P2", p2, "noninvertible and fixed-point special cases stay bounded")
        exp.check("P3", p3, "generic backward/forward ordering control")
        exp.fail_check("C1", generic["backward_support_size"] > generic["two_label_bound_at_k8"],
                       f"generic backward support {generic['backward_support_size']} > two-label bound {generic['two_label_bound_at_k8']}")
        exp.fail_check("C2", generic["ordering_sets_differ"],
                       f"backward size {generic['backward_support_size']} vs listed-order size {generic['forward_listed_support_size']}")
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
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "PASS" if all(ok for _, _, ok, _ in exp._results) else "FAIL"
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)],
                    metadata={"max_k": MAX_K, "payload_guard_bytes": MAX_BYTES,
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
