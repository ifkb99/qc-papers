"""Independent exact audit of the uniform-prefix near-collision certificate.

The production helper counts

    {1 <= q < 2**d : dist_r((2**(t-d))*q, 0) <= 2*R}

with two Euclidean floor sums.  This experiment compares that result with a
small direct integer loop, then checks boundary and width-63 cases by closed
integer formulas without enumerating their exponentially many histories.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The helper count, certificate flag, and exact metadata agree with direct
      set-counting for every r=1..40, t=0..8, d<=t, R=0..4.
  P2  The stated edge formulas hold: d=0 is vacuous, r=1 collides for every
      nonempty prefix, D>=floor(r/2) covers all histories, and zero step has
      the expected count.  Width-63 dyadic and coprime full-cycle cases agree
      with independently derived integer quotients without enumeration.
  C1  A deliberately weak strict-separation test (rejecting only distances
      <2R) incorrectly certifies touching intervals at r=60,t=6,d=2,R=6;
      the helper must refuse that case. A normalized basis-state witness at
      the shared endpoint yields a NONUNIFORM prefix despite this weak test.

This is arithmetic-only.  It allocates no orbit/prefix arrays and stores only
aggregate rows plus a handful of edge cases.  The direct family sweep is
preflighted at 193600 loop terms, below the 2M bounded-test budget.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment
from lab.periodic import uniform_prefix_certificate


MAX_DIRECT_TERMS = 2_000_000
MAX_MISMATCH_EXAMPLES = 16
LIMIT = (1 << 63) - 1

def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"uniform_prefix_counts_{stamp}.json"


def distance(period: int, residue: int) -> int:
    residue %= period
    return min(residue, period - residue)


def direct_count(period: int, width: int, prefix_bits: int, radius: int) -> tuple[int, int]:
    """Return (count, loop terms), using only native integer arithmetic."""
    histories = 1 << prefix_bits
    step = 1 << (width - prefix_bits)
    threshold = 2 * radius
    count = 0
    terms = 0
    for q in range(1, histories):
        terms += 1
        if distance(period, step * q) <= threshold:
            count += 1
    return count, terms


def expected_metadata(result, period, width, prefix_bits, radius, count):
    histories = 1 << prefix_bits
    step = (1 << (width - prefix_bits)) % period
    certified = count == 0
    return (
        result["period"] == period
        and result["exponent_width"] == width
        and result["prefix_bits"] == prefix_bits
        and result["radius"] == radius
        and result["high_history_count"] == histories
        and result["step_mod_period"] == step
        and result["near_collision_count"] == count
        and result["certified_uniform"] is certified
        and result["vacuous_prefix"] is (prefix_bits == 0)
        and (result["prefix_probability_numerator"] == (1 if certified else None))
        and (result["prefix_probability_denominator"] == (histories if certified else None))
        and result["orbit_table_entries"] == 0
        and result["prefix_table_entries"] == 0
        and result["floor_sum_calls"] in (0, 2)
        and 0 <= result["euclidean_iterations"] <= 4*period.bit_length()+4
    )


def run_family():
    cases = 0
    predicted_terms = 0
    for period in range(1, 41):
        for width in range(9):
            for prefix_bits in range(width + 1):
                for radius in range(5):
                    cases += 1
                    predicted_terms += (1 << prefix_bits) - 1
    if predicted_terms > MAX_DIRECT_TERMS:
        raise MemoryError("direct integer-term preflight exceeds 2M budget")

    checked = 0
    direct_terms = 0
    max_floor_calls = 0
    max_euclidean_iterations = 0
    failures = []
    mismatch_count = 0
    for period in range(1, 41):
        for width in range(9):
            for prefix_bits in range(width + 1):
                for radius in range(5):
                    expected, terms = direct_count(period, width, prefix_bits, radius)
                    direct_terms += terms
                    result = uniform_prefix_certificate(period, width, prefix_bits, radius)
                    max_floor_calls = max(max_floor_calls, result["floor_sum_calls"])
                    max_euclidean_iterations = max(
                        max_euclidean_iterations, result["euclidean_iterations"])
                    if not expected_metadata(result, period, width, prefix_bits,
                                             radius, expected):
                        mismatch_count += 1
                        if len(failures) < MAX_MISMATCH_EXAMPLES:
                            failures.append({
                            "period": period, "width": width,
                            "prefix_bits": prefix_bits, "radius": radius,
                            "expected_count": expected, "result": result,
                        })
                    checked += 1
    return {
        "cases": cases,
        "checked": checked,
        "predicted_direct_terms": predicted_terms,
        "actual_direct_terms": direct_terms,
        "max_floor_sum_calls": max_floor_calls,
        "max_euclidean_iterations": max_euclidean_iterations,
        "failures": failures,
        "mismatch_count": mismatch_count,
    }


def edge_case(period, width, prefix_bits, radius, expected, label):
    result = uniform_prefix_certificate(period, width, prefix_bits, radius)
    ok = expected_metadata(result, period, width, prefix_bits, radius, expected)
    return {
        "label": label,
        "period": period,
        "width": width,
        "prefix_bits": prefix_bits,
        "radius": radius,
        "expected_count": expected,
        "actual_count": result["near_collision_count"],
        "certified_uniform": result["certified_uniform"],
        "floor_sum_calls": result["floor_sum_calls"],
        "euclidean_iterations": result["euclidean_iterations"],
        "ok": ok,
    }


def run_edges():
    # For d=0, H=1 and there are no nonzero q histories.
    rows = [edge_case(17, 8, 0, 4, 0, "zero_bit_vacuous")]
    # At r=1 every nonempty history is exactly at distance zero.
    rows.append(edge_case(1, 4, 1, 0, 1, "period_one_nonempty"))
    # D >= floor(r/2) covers every circular residue.
    rows.append(edge_case(9, 4, 3, 2, 7, "covering_circle"))
    # The step is zero modulo r, so every nonzero q collides.
    rows.append(edge_case(32, 6, 1, 0, 1, "zero_step"))
    # Explicit equality at distance D=2R: q=3 gives residue 48 and distance 12.
    rows.append(edge_case(60, 6, 2, 6, 1, "exact_touching_distance"))

    # Width 63: all expected counts are closed integer formulas, not loops.
    h = 1 << 63
    dyadic_period = 1 << 62
    rows.append(edge_case(dyadic_period, 63, 1, 0, 1, "width63_dyadic_zero_step"))
    rows.append(edge_case(dyadic_period, 63, 63, 0,
                          (h - 1) // dyadic_period,
                          "width63_dyadic_full_cycles"))
    coprime_period = 37
    rows.append(edge_case(coprime_period, 63, 63, 0,
                          (h - 1) // coprime_period,
                          "width63_coprime_full_cycles"))
    odd_period = (1 << 62) - 1
    rows.append(edge_case(odd_period, 63, 1, 0, 0,
                          "width63_coprime_single_step"))
    return rows


def run_control():
    period, width, prefix_bits, radius = 60, 6, 2, 6
    result = uniform_prefix_certificate(period, width, prefix_bits, radius)
    # Wrong predicate: it only rejects strict overlap (distance < 2R), so it
    # incorrectly calls equality/touching a certificate.
    histories = 1 << prefix_bits
    step = 1 << (width - prefix_bits)
    strict_overlap_count = sum(
        distance(period, step * q) < 2 * radius for q in range(1, histories)
    )
    weak_certifies = strict_overlap_count == 0
    # A separate state allowed by the support promise (not the N61 circuit):
    # for each low l choose high-history labels l+labels[h]. Histories 0,3
    # share endpoint -6==54; histories 1,2 are at their centers. Every state
    # is normalized, with displacement <=6, but Fourier interference remains.
    labels = [(-radius)%period, step%period, (2*step)%period,
              (3*step+radius)%period]
    support_distances = [distance(period, labels[h]-step*h) for h in range(histories)]
    prefix_numerators = []
    gaussian_terms = 0
    for z in range(histories):
        numerator = 0j
        for h in range(histories):
            for hp in range(histories):
                numerator += (1j)**(z*(h-hp)) * (labels[h] == labels[hp])
                gaussian_terms += 1
        assert numerator.imag == 0 and numerator.real == int(numerator.real)
        prefix_numerators.append(int(numerator.real))
    detected = (result["near_collision_count"] == 1
                and not result["certified_uniform"] and weak_certifies
                and max(support_distances) <= radius
                and prefix_numerators == [6,4,2,4]
                and sum(prefix_numerators) == histories**2)
    return {
        "case": {"period": period, "width": width,
                 "prefix_bits": prefix_bits, "radius": radius},
        "helper_count": result["near_collision_count"],
        "helper_certified": result["certified_uniform"],
        "strict_overlap_count": strict_overlap_count,
        "weak_predicate_certified": weak_certifies,
        "touching_q": 3,
        "touching_distance": 12,
        "detected_false_certificate": detected,
        "support_witness_labels_relative_to_low_input": labels,
        "support_witness_distances": support_distances,
        "prefix_numerators": prefix_numerators,
        "prefix_denominator": histories**2,
        "witness_gaussian_integer_terms": gaussian_terms,
        "witness_scope": "a distinct normalized state satisfying the supplied support promise, not the fixed phase circuit",
    }


def main():
    exp = Experiment("uniform_prefix_counts", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "floor-sum helper agrees with direct modular distance counts")
    exp.predict("P2", "edge, covering-circle, zero-step, and width-63 formulas agree")
    exp.must_fail("C1", "weak strict-separation accepts touching support intervals")
    started = time.perf_counter()
    report = {
        "status": "FAIL",
        "resource_scope": {
            "max_direct_terms": MAX_DIRECT_TERMS,
            "mismatch_examples_cap": MAX_MISMATCH_EXAMPLES,
            "probe_numeric_arrays_allocated": False,
            "memory_scope": "constant-size integer loop state and capped report examples; not total RSS",
        },
    }
    p1 = p2 = c1 = False
    try:
        family = run_family()
        edges = run_edges()
        control = run_control()
        report.update({"family": family, "edge_cases": edges, "control": control})
        p1 = (family["checked"] == family["cases"] == 9000
              and family["actual_direct_terms"] == family["predicted_direct_terms"]
              and family["mismatch_count"] == 0
              and family["actual_direct_terms"]+3 <= MAX_DIRECT_TERMS)
        p2 = (all(row["ok"] for row in edges)
              and all(isinstance(row["actual_count"], int) for row in edges))
        c1 = bool(control["detected_false_certificate"])
        report["status"] = "PASS" if p1 and p2 and c1 else "FAIL"
        exp.check("P1", p1, "all bounded direct counts and helper metadata agree")
        exp.check("P2", p2, "edge and width-63 closed integer checks hold")
        exp.fail_check("C1", c1, "strict separation rejects overlap but misses touching")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before bounded family completed")
        exp.check("P2", False, "exception before edge checks completed")
        exp.fail_check("C1", False, "exception before touching control")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "helper": "lab.periodic.uniform_prefix_certificate",
        "direct_count": "native Python integer modular distances only",
        "formula_scope": "count <=2R; nonzero count is inconclusive",
        "no_orbit_or_prefix_arrays": True,
        "width63_cases": "closed quotient formulas; no history enumeration",
        "limits": {"period": LIMIT, "width": 63, "direct_terms": MAX_DIRECT_TERMS},
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
