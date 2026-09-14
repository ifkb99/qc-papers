"""Exact selector audit for nested multi-phase progression construction.

This is an integer-only follow-up to C81/TODO40.  It compares direct
``(h,k,rho,z)`` enumeration with the closed T/S formulas for coefficient-pair
visits and pre-cancellation geometric candidates.  It also checks that cuts
at phase positions and the final endpoint contain the all-cut minima.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Direct clipped/unclipped counts agree with
      C=min(H*K*b^2*P,T) and G=min(H*K*R*P,S), including r=b aliases.
  P2  Between phase positions the exact naive structural cost is
      nonincreasing, so phase-position/end cuts contain its minimum.
  C1  The largest-gap K*P selector chooses the wrong cut for the frozen
      r=60,b=3,s=7,(v1,v2)=(3,6),j=0 cost control: direct work is 78 at
      cut 6 but 60 at cut 7.

The test has no amplitudes, phases, propagators, or orbit tables.  Its cost
model counts grouping pairs, coefficient-pair visits, candidate products and
pointwise phase calls for the stated direct constructor only.

Run: uv run python -m experiments.experiment_nested_phase_selector
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment


MAX_VISITS = 100_000
MAX_PAYLOAD = 1 << 20
H = 2


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_selector_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def v2(value):
    value = int(value)
    result = 0
    while value and value % 2 == 0:
        value //= 2
        result += 1
    return result


def phase_period(r, positions, cut):
    left = [int(v) for v in positions if int(v) < int(cut)]
    if not left:
        return 1
    vmax = max(left)
    return (1 << vmax) // math.gcd(int(r), 1 << vmax)


def difference_multiplicities(r, b, ledger=None):
    values = {}
    for u in range(int(b)):
        for q in range(int(b)):
            if ledger is not None:
                ledger.charge("difference_pair", 1)
            d = (q - u) % int(r)
            values[d] = values.get(d, 0) + 1
    return tuple(sorted(values.items()))


def history_cells(r, b, s, j=0, ledger=None):
    """C values for the two frozen high histories, using j-e=u+p-q."""
    L = 1 << int(s)
    cells = []
    for h in range(H):
        if ledger is not None:
            ledger.charge("history_setup", 1)
        jh = (int(j) - L * h) % int(r)
        cells.append(int(b) * (jh // int(b)))
    return tuple(cells)


def n_count(limit, r, rho):
    return max(0, 1 + (int(limit) - 1 - int(rho)) // int(r))


def full_TS(r, b, s, cells, differences, ledger):
    """Full-L pair incidence T and distinct-support count S."""
    L = 1 << int(s)
    T = 0
    S = 0
    for c in cells:
        seen = set()
        for d, multiplicity in differences:
            ledger.charge("full_pair", 1)
            rho = (int(c) + int(d)) % int(r)
            count = n_count(L, r, rho)
            T += int(multiplicity) * count
            if count:
                seen.add(rho)
        for rho in seen:
            ledger.charge("full_support", 1)
            S += n_count(L, r, rho)
    return T, S


class Ledger:
    def __init__(self, reserved):
        self.reserved = int(reserved)
        self.total = 0
        self.counts = {}

    def charge(self, name, amount=1):
        amount = int(amount)
        if amount < 0:
            raise ValueError("negative visit charge")
        if self.total + amount > MAX_VISITS:
            raise MemoryError(f"actual {name} visits exceed frozen cap")
        if self.total + amount > self.reserved:
            raise MemoryError(f"actual {name} visits exceed preflight reserve")
        self.total += amount
        self.counts[name] = self.counts.get(name, 0) + amount


def direct_counts(r, b, s, cells, differences, positions, cut, period, ledger):
    """Direct integer (h,k,d,z) enumeration, independent of T/S formulas."""
    L = 1 << int(s)
    A = 1 << int(cut)
    K = L // A
    P = int(period)
    C = G = 0
    for c in cells:
        for k in range(K):
            # Actually enumerate the pairs rather than charging a symbolic
            # grouping bound as if that many loop visits had occurred.
            pairs = {}
            for u in range(b):
                for q in range(b):
                    ledger.charge("grouping_pair")
                    rho = (int(c) + q - u - A * k) % int(r)
                    pairs[rho] = pairs.get(rho, 0) + 1
            for rho, multiplicity in pairs.items():
                ledger.charge("rho_class", 1)
                count = n_count(A, r, rho)
                for z in range(min(P, count)):
                    ledger.charge("candidate_class", 1)
                    for _pair in range(multiplicity):
                        ledger.charge("coefficient_pair")
                        C += 1
                    G += 1
    return C, G


def cut_record(r, b, s, cells, differences, positions, cut, period, ts, ledger):
    L = 1 << int(s)
    A = 1 << int(cut)
    K = L // A
    P = int(period)
    T, S = ts
    R = min(int(r), 2 * int(b) - 1)
    ledger.charge("formula", 1)
    C_formula = min(H * K * b * b * P, T)
    G_formula = min(H * K * R * P, S)
    grouping = H * K * b * b
    C_direct, G_direct = direct_counts(
        r, b, s, cells, differences, positions, cut, P, ledger)
    phase_calls = H + len(positions) * C_direct
    products = (len(positions) + 1) * C_direct + G_direct
    # This is the explicitly stated direct-constructor scalar cost.  It is
    # not a timing model and does not charge arbitrary oracle internals.
    cost = grouping + C_direct + G_direct + phase_calls + products
    return {
        "cut": int(cut), "A": A, "K": K, "P": P,
        "grouping_pairs": grouping,
        "coefficient_pairs": C_direct,
        "candidate_count": G_direct,
        "coefficient_formula": C_formula,
        "candidate_formula": G_formula,
        "phase_calls": phase_calls,
        "products": products,
        "cost": cost,
        "structural_factor": K * P,
        "all_classes_clipped": P <= A // r,
    }


def preflight_case(r, b, s, positions):
    L = 1 << int(s)
    cells = history_cells(r, b, s)
    differences = difference_multiplicities(r, b)
    R = min(int(r), 2 * int(b) - 1)
    total = H + b * b + H * b * b * 2  # setup, differences, T/S summaries
    for cut in range(int(s) + 1):
        A = 1 << cut
        K = L // A
        P = phase_period(r, positions, cut)
        F = min(P, (A + r - 1) // r)
        # Conservative direct bounds before any execution.
        grouping = H * K * b * b
        rho_classes = H * K * R
        candidates = H * K * R * F
        pairs = H * K * b * b * F
        # phase-period calls, monotonicity comparisons, formula checks, and
        # selector summaries are reserved separately from structural loops.
        total += (grouping + rho_classes + candidates + pairs
                  + (s + 1) + (s + 4) + 3 * (s + 1) + 4)
    return total, cells, differences


def candidate_cuts(s, positions):
    return tuple(sorted({int(s), *(int(v) for v in positions)}))


def check_case(fixture, positions, ledger):
    r, b, s = fixture["r"], fixture["b"], fixture["s"]
    cells = history_cells(r, b, s, ledger=ledger)
    differences = difference_multiplicities(r, b, ledger=ledger)
    ts = full_TS(r, b, s, cells, differences, ledger)
    records = []
    for cut in range(s + 1):
        ledger.charge("phase_period", 1)
        records.append(cut_record(
            r, b, s, cells, differences, positions, cut,
            phase_period(r, positions, cut), ts, ledger))
    by_cut = {row["cut"]: row for row in records}
    candidates = candidate_cuts(s, positions)
    # Costs decrease between phase positions; a phase at the endpoint remains
    # right-side, so the endpoint itself is included as a candidate cut.
    monotone = True
    phase_set = set(int(v) for v in positions)
    for cut in range(s):
        if cut not in phase_set:
            ledger.charge("selector_comparison", 1)
            monotone &= by_cut[cut]["cost"] >= by_cut[cut + 1]["cost"]
    ledger.charge("selector_comparison", 2 * (len(records) + len(candidates)))
    all_min = min(row["cost"] for row in records)
    candidate_min = min(by_cut[cut]["cost"] for cut in candidates)
    structural_min = min(row["structural_factor"] for row in records)
    structural_candidate_min = min(by_cut[cut]["structural_factor"]
                                   for cut in candidates)
    if all_min != candidate_min or structural_min != structural_candidate_min:
        raise AssertionError("candidate cuts missed an all-cut minimum")
    return {
        "fixture": fixture,
        "positions": list(positions),
        "cells": list(cells),
        "difference_multiplicities": [list(x) for x in differences],
        "T": ts[0], "S": ts[1],
        "candidate_cuts": list(candidates),
        "records": records,
        "all_cost_min": all_min,
        "candidate_cost_min": candidate_min,
        "all_structural_min": structural_min,
        "candidate_structural_min": structural_candidate_min,
        "between_phase_cost_monotone": bool(monotone),
    }


def main():
    exp = Experiment("nested_phase_selector", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "direct clipped counts equal the closed T/S formulas")
    exp.predict("P2", "phase-position and endpoint cuts contain all structural-cost minima")
    exp.must_fail("C1", "largest-gap K*P chooses the wrong frozen direct-work cut")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = c1 = False
    ledger = None
    rows = []
    try:
        # Freeze fixtures and phase families before any measured loop.
        fixtures = (
            {"name": "r60_long", "r": 60, "b": 3, "s": 7},
            {"name": "r6_clipping", "r": 6, "b": 2, "s": 4},
            {"name": "alias_r3_b3", "r": 3, "b": 3, "s": 4},
            {"name": "alias_r2_b2", "r": 2, "b": 2, "s": 3},
        )
        phase_families = {
            "r60_long": ((), (0,), (7,), (3, 3), (3, 6)),
            "r6_clipping": ((), (0,), (2,), (4,), (2, 2)),
            "alias_r3_b3": ((), (0,), (1,), (4,), (1, 1)),
            "alias_r2_b2": ((), (0,), (1,), (3,), (1, 1)),
        }
        reserved = 0
        preflight_rows = []
        for fixture in fixtures:
            name = fixture["name"]
            for positions in phase_families[name]:
                bound, cells, differences = preflight_case(
                    fixture["r"], fixture["b"], fixture["s"], positions)
                reserved += bound
                preflight_rows.append({"fixture": name,
                                       "positions": list(positions),
                                       "reserved_visits": bound,
                                       "cells": list(cells),
                                       "difference_count": len(differences)})
        payload = 192 * 1024  # integer summaries plus retained JSON records
        report["preflight"] = {"reserved_visits": reserved, "payload_bytes": payload,
                               "max_visits": MAX_VISITS, "max_payload_bytes": MAX_PAYLOAD,
                               "rows": preflight_rows}
        if reserved > MAX_VISITS:
            raise MemoryError(f"preflight visits {reserved}>{MAX_VISITS}")
        if payload > MAX_PAYLOAD:
            raise MemoryError(f"preflight payload {payload}>{MAX_PAYLOAD}")
        ledger = Ledger(reserved)
        rows = []
        for fixture in fixtures:
            for positions in phase_families[fixture["name"]]:
                rows.append(check_case(fixture, positions, ledger))

        # Frozen wrong-cost control: structural factor prefers cut 6, while
        # the exact direct work prefers cut 7 (78 versus 60 grouping+coeff).
        # This exact case is already in rows; reuse it rather than silently
        # charging a second unreserved traversal of the frozen control.
        control = next(row for row in rows
                       if row["fixture"]["name"] == "r60_long"
                       and row["positions"] == [3, 6])
        rec6 = control["records"][6]
        rec7 = control["records"][7]
        c1 = bool(rec6["structural_factor"] < rec7["structural_factor"]
                  and rec6["grouping_pairs"] + rec6["coefficient_pairs"] == 78
                  and rec7["grouping_pairs"] + rec7["coefficient_pairs"] == 60
                  and rec6["cost"] > rec7["cost"])
        audit_records = sum(len(row["records"]) for row in rows)
        ledger.charge("audit_comparison", 3 * audit_records)
        clipping_branches = {record["all_classes_clipped"]
                             for row in rows for record in row["records"]}
        p1 = bool(all(record["coefficient_pairs"] == record["coefficient_formula"]
                      and record["candidate_count"] == record["candidate_formula"]
                      for row in rows for record in row["records"])
                  and clipping_branches == {False, True})
        p2 = bool(all(row["between_phase_cost_monotone"]
                      and row["all_cost_min"] == row["candidate_cost_min"]
                      and row["all_structural_min"] == row["candidate_structural_min"]
                      for row in rows))
        report.update({
            "status": "PASS" if p1 and p2 and c1 else "FAIL",
            "fixtures": fixtures,
            "phase_families": {key: [list(x) for x in value]
                               for key, value in phase_families.items()},
            "preflight": {"reserved_visits": reserved,
                          "payload_bytes": payload,
                          "max_visits": MAX_VISITS,
                          "max_payload_bytes": MAX_PAYLOAD,
                          "rows": preflight_rows},
            "rows": rows,
            "cost_control": control,
            "actual_visits": ledger.total,
            "actual_visit_counts": ledger.counts,
            "fixture_phase_cases": len(rows), "cut_cases": audit_records,
            "both_clipping_branches": clipping_branches == {False, True},
            "checks": {"P1": p1, "P2": p2, "C1": c1},
            "python_version": sys.version,
            "platform": platform.platform(),
            "scope": "exact integer structural counts; no amplitudes or oracles",
        })
        exp.check("P1", p1, "direct clipped counts match T/S formulas")
        exp.check("P2", p2, "candidate cuts contain all structural-cost minima")
        exp.fail_check("C1", c1, "largest-gap factor is not an actual direct-cost optimizer")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        report["completed_rows"] = rows
        if ledger is not None:
            report["actual_visits"] = ledger.total
            report["actual_visit_counts"] = ledger.counts
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before integer selector audit")
        exp.check("P2", False, "exception before integer selector audit")
        exp.fail_check("C1", False, "exception before frozen cost control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "scope": "bounded exact T/S selector audit",
        "max_named_visits": MAX_VISITS,
        "max_numeric_payload_bytes": MAX_PAYLOAD,
        "no_amplitudes": True,
        "no_production_changes": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
