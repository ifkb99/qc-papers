"""Exact integer audit of a nested binary cut for two earlier phases.

This is a C79/TODO40 follow-up after the C80 checkpoint. It uses no amplitudes, NumPy arrays,
orbit table, or sampler.  For r=60, s=7, source labels x=0..59 and u=0..2,
it compares the original phase-argument word for positions (3,v), v=3..7,
with the cut formula for every w=0..7.  It then partitions compatible labels
by (k,rho,z), checks the resulting arithmetic progressions and verifies that
the phase word is constant within each class.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Every compatible label has the same original and cut phase word, and
      every refined class is a disjoint arithmetic progression with constant
      word; all compatible labels are covered exactly once.
  P2  For positions (3,6), the structural factor K*P_left is 16 at cuts 3
      and 7 but 4 at cut 6 (K=2, P_left=2); observed class counts do not
      exceed their corresponding factors.
  P3  With alpha=v2(60)=2, B=5 and the largest effective insertion gap Delta,
      the frozen largest-gap formula 2^(B-Delta) equals the minimum of all
      measured cut factors for every (3,v), v=3..7.
  C1  A control that drops the left refinement and forces P_left=1 detects the
      frozen witness a=0,60,u=x=k=rho=0,w=6: phase-3 arguments are 0 and 4,
      while phase-6 arguments are both 0.

These are exact integer partition checks only.  They do not test amplitudes,
output laws, sampling, runtime, or a multi-phase production implementation.
"""
from __future__ import annotations

import json
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment


PERIOD = 60
WIDTH = 7
L = 1 << WIDTH
US = range(3)
XS = range(PERIOD)
CUTS = range(WIDTH + 1)
SECOND_POSITIONS = range(3, WIDTH + 1)
MAX_VISITS = 5_000_000
MAX_NUMERIC_BYTES = 1 << 20
TOL = 0


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_cuts_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def original_word(u, l, positions, visits):
    visits["argument"] += len(positions)
    return tuple((int(u) + (int(l) % (1 << int(v)))) % PERIOD
                  for v in positions)


def cut_word(u, x, l, cut, positions, visits):
    visits["argument"] += len(positions)
    a = int(l) % (1 << int(cut))
    k = int(l) // (1 << int(cut))
    word = []
    for v in positions:
        if v >= cut:
            shift = 1 << (int(v) - int(cut))
            argument = int(x) - (1 << int(v)) * (k // shift)
        else:
            argument = int(u) + (a % (1 << int(v)))
        word.append(argument % PERIOD)
    return tuple(word)


def factor(cut, positions, visits):
    visits["factor"] += 1
    a = 1 << int(cut)
    k_count = L // a
    left = [int(v) for v in positions if int(v) < int(cut)]
    if not left:
        p_left = 1
    else:
        vmax = max(left)
        p_left = (1 << vmax) // __import__("math").gcd(PERIOD, 1 << vmax)
    return k_count, p_left, k_count * p_left


def witness_control(visits):
    """Return the deliberately false P_left=1 prediction's witness data."""
    positions = (3, 6)
    cut = 6
    labels = (0, 60)
    actual = [original_word(0, label, positions, visits) for label in labels]
    wrong = []
    # The invalid control retains the cut-6 right formula but collapses the
    # left phase to a single value, rather than refining by n mod P_left.
    for label in labels:
        visits["witness_label"] += 1
        visits["argument"] += len(positions)
        a = label % (1 << cut)
        k = label // (1 << cut)
        wrong.append(tuple((0 if v < cut else
                            (0 - (1 << v) * (k // (1 << (v-cut)))) % PERIOD)
                           for v in positions))
    return {"labels": labels, "actual_words": actual, "wrong_words": wrong,
            "actual_differs": actual[0] != actual[1],
            "wrong_collapses": wrong[0] == wrong[1]}


def independently_constructed_partition(u, x, cut, positions, visits):
    """Build the proposed classes from k/rho/N, independent of l scanning."""
    a_size = 1 << int(cut)
    k_count = L // a_size
    _k_count, p_left, _factor = factor(cut, positions, visits)
    constructed = defaultdict(list)
    for k in range(k_count):
        visits["construction_k"] += 1
        rho = (int(x) - a_size * k - int(u)) % PERIOD
        count_n = max(0, 1 + (a_size - 1 - rho) // PERIOD)
        for z in range(min(p_left, count_n)):
            visits["construction_class"] += 1
            start = a_size * k + rho + PERIOD * z
            count = 1 + (count_n - 1 - z) // p_left
            stride = PERIOD * p_left
            labels = []
            for t in range(count):
                visits["construction_label"] += 1
                label = start + stride * t
                if not 0 <= label < L:
                    raise AssertionError("constructed progression leaves l range")
                labels.append(label)
            key = (k, rho, z)
            constructed[key] = [(label, original_word(u, label, positions, visits))
                                 for label in labels]
    return constructed


def largest_gap_prediction(positions):
    alpha = 2  # v2(60), frozen explicitly rather than inferred from samples.
    B = max(0, WIDTH - alpha)
    effective = sorted({0, B} | {
        max(0, int(v) - alpha) for v in positions})
    gaps = [b - a for a, b in zip(effective, effective[1:])]
    delta = max(gaps, default=0)
    return {"alpha": alpha, "B": B, "effective_positions": effective,
            "Delta": delta, "predicted_min_factor": 1 << (B - delta)}


def preflight():
    cases = len(SECOND_POSITIONS) * len(CUTS) * len(US) * len(XS)
    max_compatible = (L + PERIOD - 1) // PERIOD
    relation = cases * L
    arguments = cases * max_compatible * 2 * 2 * 2  # two words, two phases
    partition = cases * max_compatible
    class_checks = cases * max_compatible * 2
    reserved = relation + arguments + partition + class_checks
    # Only small aggregate summaries and one <=128-label streaming case are
    # retained; reserve conservatively for integer summaries and JSON copies,
    # with no numeric orbit/Q-by-r payload.
    payload = 64 * 1024
    groups = len(SECOND_POSITIONS) * len(US) * len(XS)
    construction_k = groups * sum(L // (1 << cut) for cut in CUTS)
    # P_left is at most 16 for this r and v<=7; this is conservative across
    # every second insertion position and is charged before construction.
    construction_classes = groups * sum(
        (L // (1 << cut)) * min(1 << cut, 16) for cut in CUTS)
    construction_labels = groups * len(CUTS) * L
    factor_work = len(SECOND_POSITIONS) * len(CUTS) * 8
    construction_word_work = cases * max_compatible * 2 + cases
    witness_work = 64
    reserved += construction_k + construction_classes + construction_labels
    reserved += factor_work + witness_work
    reserved += construction_word_work
    if reserved > MAX_VISITS:
        raise MemoryError(f"integer visit preflight {reserved}>{MAX_VISITS}")
    if payload > MAX_NUMERIC_BYTES:
        raise MemoryError(f"numeric payload preflight {payload}>{MAX_NUMERIC_BYTES}")
    return {"cases": cases, "max_compatible_labels_per_case": max_compatible,
            "reserved_relation_visits": relation,
            "reserved_argument_visits": arguments,
            "reserved_partition_visits": partition,
            "reserved_class_check_visits": class_checks,
            "reserved_construction_k_visits": construction_k,
            "reserved_construction_class_visits": construction_classes,
            "reserved_construction_label_visits": construction_labels,
            "reserved_factor_work": factor_work,
            "reserved_construction_word_and_factor_work": construction_word_work,
            "reserved_witness_work": witness_work,
            "reserved_total_visits": reserved,
            "planned_numeric_payload_bytes": payload,
            "max_numeric_payload_bytes": MAX_NUMERIC_BYTES,
            "max_integer_visits": MAX_VISITS}


def main():
    exp = Experiment("nested_phase_cuts", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "cut words and refined partitions are exact and exhaustive")
    exp.predict("P2", "nested factor is smaller at the interior cut")
    exp.predict("P3", "largest-gap factor matches the minimum over all cuts")
    exp.must_fail("C1", "dropping left refinement misses the explicit witness")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        visits = defaultdict(int)
        summaries = {}
        factor_summary = {}
        max_class_count = 0
        max_word_classes = 0
        for second in SECOND_POSITIONS:
            positions = (3, int(second))
            for cut in CUTS:
                k_count, p_left, structural_factor = factor(cut, positions, visits)
                key = f"{positions[0]},{positions[1]}:w{cut}"
                factor_summary[key] = {
                    "K": k_count, "P_left": p_left,
                    "structural_factor": structural_factor,
                    "max_observed_classes": 0,
                    "cases": 0,
                }
                for u in US:
                    for x in XS:
                        # Charge the conservative entire case before any of
                        # its l/word/partition loops.
                        a_size = 1 << cut
                        k_count = L // a_size
                        construction_case_reserved = (
                            k_count + k_count * min(a_size, 16) + L
                            + 1 + 2 * ((L + PERIOD - 1) // PERIOD))
                        case_reserved = (L
                                         + ((L + PERIOD - 1) // PERIOD)
                                         * (2 * len(positions) * 2 + 1 + 2)
                                         + construction_case_reserved)
                        if visits["reserved"] + case_reserved > MAX_VISITS:
                            raise MemoryError("case reserve exceeds integer visit cap")
                        visits["reserved"] += case_reserved
                        classes = defaultdict(list)
                        compatible = 0
                        for l in range(L):
                            visits["relation"] += 1
                            if (int(u) + l - int(x)) % PERIOD != 0:
                                continue
                            compatible += 1
                            old = original_word(u, l, positions, visits)
                            new = cut_word(u, x, l, cut, positions, visits)
                            if old != new:
                                raise AssertionError(
                                    f"word mismatch positions={positions} cut={cut} "
                                    f"u={u} x={x} l={l}: {old}!={new}")
                            a = l % (1 << cut)
                            k = l // (1 << cut)
                            rho = (int(x) - (1 << cut) * k - int(u)) % PERIOD
                            if (a - rho) % PERIOD != 0:
                                raise AssertionError("rho does not describe compatible a")
                            n = (a - rho) // PERIOD
                            z = n % p_left
                            visits["partition"] += 1
                            classes[(k, rho, z)].append((l, old))
                        if compatible > (L + PERIOD - 1) // PERIOD:
                            raise AssertionError("too many compatible labels")
                        covered = sum(len(values) for values in classes.values())
                        if covered != compatible:
                            raise AssertionError("partition lost or duplicated labels")
                        constructed = independently_constructed_partition(
                            u, x, cut, positions, visits)
                        reference = {key: sorted(values)
                                     for key, values in classes.items()}
                        constructed = {key: sorted(values)
                                       for key, values in constructed.items()}
                        if reference != constructed:
                            raise AssertionError(
                                "independent k/rho progression construction mismatch")
                        observed_classes = len(classes)
                        if observed_classes > structural_factor:
                            raise AssertionError("partition exceeds structural factor")
                        for values in classes.values():
                            visits["class_check"] += 1
                            words = {word for _label, word in values}
                            if len(words) != 1:
                                raise AssertionError("phase word varies within class")
                            labels = sorted(label for label, _word in values)
                            if len(labels) > 1:
                                expected_stride = PERIOD * p_left
                                if any(b - a != expected_stride
                                       for a, b in zip(labels, labels[1:])):
                                    raise AssertionError("class is not one progression")
                        factor_summary[key]["max_observed_classes"] = max(
                            factor_summary[key]["max_observed_classes"], observed_classes)
                        factor_summary[key]["cases"] += 1
                        max_class_count = max(max_class_count, observed_classes)
                        max_word_classes = max(max_word_classes, len(classes))
        witness = witness_control(visits)
        c1 = bool(witness["actual_differs"] and witness["wrong_collapses"])
        targeted = {key: value for key, value in factor_summary.items()
                    if key in ("3,6:w3", "3,6:w6", "3,6:w7")}
        # The relation loop is over every frozen source-label/l value, while
        # argument/partition loops are restricted to compatible labels.
        p1 = bool(visits["relation"] == report["preflight"]["cases"] * L
                  and visits["partition"] > 0
                  and visits["class_check"] > 0
                  and visits["construction_k"] > 0
                  and visits["construction_class"] > 0
                  and visits["construction_label"] > 0
                  and visits["reserved"] <= MAX_VISITS
                  and visits["reserved"] <= report["preflight"]["reserved_total_visits"]
                  and sum(value for key, value in visits.items() if key != "reserved")
                  <= report["preflight"]["reserved_total_visits"])
        p2 = bool(targeted["3,6:w3"]["structural_factor"] == 16
                  and targeted["3,6:w6"]["K"] == 2
                  and targeted["3,6:w6"]["P_left"] == 2
                  and targeted["3,6:w6"]["structural_factor"] == 4
                  and targeted["3,6:w7"]["structural_factor"] == 16
                  and all(value["max_observed_classes"]
                          <= value["structural_factor"]
                          for value in factor_summary.values()))
        largest_gap = {}
        for second in SECOND_POSITIONS:
            visits["gap_formula_calls"] += 1
            largest_gap[f"3,{second}"] = largest_gap_prediction((3, int(second)))
        p3 = bool(all(
            min(factor_summary[f"3,{second}:w{cut}"]["structural_factor"]
                for cut in CUTS)
            == largest_gap[f"3,{second}"]["predicted_min_factor"]
            for second in SECOND_POSITIONS)
                  and visits["factor"] == (len(SECOND_POSITIONS) * len(CUTS)
                                           + report["preflight"]["cases"]))
        report.update({
            "fixture": {"r": PERIOD, "s": WIDTH, "L": L,
                        "u_values": list(US), "x_values": list(XS),
                        "phase_positions": [list((3, v)) for v in SECOND_POSITIONS],
                        "cuts": list(CUTS)},
            "factor_summary": factor_summary,
            "targeted_3_6_factors": targeted,
            "largest_gap_predictions": largest_gap,
            "witness_control": witness,
            "visits": dict(visits),
            "max_observed_class_count": max_class_count,
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
            "status": "PASS" if p1 and p2 and p3 and c1 else "FAIL",
            "python_version": sys.version,
            "platform": platform.platform(),
            "scope": "exact integer words/partitions only; no amplitudes or sampler",
        })
        exp.check("P1", p1, "exact words, coverage, and progression classes")
        exp.check("P2", p2, "interior nested factor and endpoint factors")
        exp.check("P3", p3, "largest-gap formula over all frozen cuts")
        exp.fail_check("C1", c1, "wrong P_left=1 witness detected")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before integer audit")
        exp.fail_check("C1", False, "exception before witness control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": "r60,s7,positions(3,v),v=3..7",
        "scope": "exact integer structural pilot",
        "max_integer_visits": MAX_VISITS,
        "max_numeric_payload_bytes": MAX_NUMERIC_BYTES,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
