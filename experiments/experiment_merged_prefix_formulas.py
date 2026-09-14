"""Bounded all-prefix audit of coherent sector merging.

MergedCoherentPrefixes is checked as an amplitude oracle against the existing
C59 route expansion and its independent tiny full-r direct_prefix helper.  The
audit enumerates every valid sector, exponent, stop, boundary, measured output
prefix, and fine-work coordinate on small fixtures, while retaining only error
summaries and bounded laws.  No merged-path history list or generic propagator
is introduced here.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  every valid merged prefix vector equals C59 prefix_vector and
      independent direct_prefix, including intermediate boundaries and exact
      zero amplitudes.
  P2  merged joint probabilities agree with C59 and independent direct_joint;
      samples and exposed stats remain table-free and no-history.
  P3  endpoint, t=0/1, b=2, fixed-point, wrap, collision, and π-zero fixtures
      obey the same prefix identity.

  C1  treating an intermediate background boundary as arithmetic changes a
      nontrivial prefix.
  C2  omitting the prepared-tail factor changes an intermediate prefix.
  C3  omitting the inverse-QFT output phase changes a terminal measured prefix.

This is a bounded complex128 diagnostic, not a finite-precision certificate.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit


MAX_DENSE_BYTES = 16 << 20
TOL = 2e-10

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all valid merged prefixes match C59 and independent direct_prefix")
exp.predict("P2", "joint laws and table-free no-history stats match independent references")
exp.predict("P3", "endpoint, b=2, fixed-point, wrap, collision, and zero fixtures pass")
exp.must_fail("C1", "wrong intermediate boundary inclusion changes a prefix")
exp.must_fail("C2", "omitting the prepared tail changes an intermediate prefix")
exp.must_fail("C3", "omitting the inverse-QFT phase changes a measured prefix")


def guard(shape, dtype=np.complex128, label="array"):
    entries = math.prod(int(x) for x in shape)
    payload = entries * np.dtype(dtype).itemsize
    if entries < 0 or payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def stamp(prefix="merged_prefix_formulas"):
    now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{now}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"{prefix}_{now}_{serial}.json"
    return path


def rx(theta, block):
    out = np.eye(block, dtype=complex)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    out[:2, :2] = ((c, -1j * s), (-1j * s, c))
    return out


def rz(theta, block):
    out = np.eye(block, dtype=complex)
    out[0, 0] = np.exp(-1j * theta / 2)
    out[1, 1] = np.exp(1j * theta / 2)
    return out


def fixture_primary():
    return CoherentReflectionCircuit(
        9, 3, 4,
        {0: rx(np.pi / 4, 3), 1: rx(np.pi / 7, 3),
         3: rz(np.pi / 5, 3)},
        {0: (0, np.pi / 5), 2: (1, np.pi / 3),
         3: (0, np.pi / 7), 4: (1, np.pi / 4)})


def fixture_t0():
    return CoherentReflectionCircuit(
        9, 3, 0, {0: rx(np.pi / 6, 3)}, {0: (0, np.pi / 3)})


def fixture_t1():
    return CoherentReflectionCircuit(
        9, 3, 1, {0: rx(np.pi / 6, 3), 1: rz(np.pi / 5, 3)},
        {0: (1, np.pi / 3), 1: (0, np.pi / 4)})


def fixture_b2():
    return CoherentReflectionCircuit(
        10, 2, 2,
        {0: rx(np.pi / 5, 2), 1: rz(np.pi / 7, 2)},
        {0: (4, np.pi / 5), 1: (0, np.pi / 3),
         2: (1, np.pi / 4)})


def fixture_zero():
    # Float pi is a deterministic-route limit, not an exact zero cosine.
    return CoherentReflectionCircuit(
        9, 3, 1, {0: rx(np.pi / 4, 3), 1: rx(np.pi / 7, 3)},
        {0: (0, np.pi), 1: (1, np.pi)})


def direct_prefix(circuit, final, exponent, stop, boundary, measured, output):
    from experiments.experiment_coherent_route_sampling import direct_prefix as ref
    return np.asarray(ref(circuit, final, exponent, stop, boundary,
                          measured, output), dtype=complex)


def stats_no_history(stats):
    # Full stats expose these counters; last_prefix_stats intentionally only
    # exposes the local prefix costs.  If a forbidden counter is present it
    # must be zero, while absence itself is not a hidden history table.
    forbidden = ("history_enumerations", "orbit_table_entries",
                 "output_table_entries", "full_sector_table_entries",
                 "sector_table_entries")
    return all(stats.get(key, 0) == 0 for key in forbidden)


def run_fixture(name, circuit, MergedCoherentPrefixes, *, sample_seed):
    guard((circuit.period, circuit.b), label=f"{name} work workspace")
    merged = MergedCoherentPrefixes(circuit)
    stats = merged.stats()
    if not stats_no_history(stats):
        raise AssertionError(f"{name}: merged stats retain forbidden tables/history")
    direct_max = 0.0
    merged_max = 0.0
    zero_count = near_zero_count = 0
    query_count = 0
    wrong_boundary_delta = 0.0
    wrong_tail_delta = 0.0
    wrong_phase_delta = 0.0
    for final in range(circuit.sectors):
        for exponent in range(1 << circuit.width):
            for stop in range(circuit.width + 1):
                for boundary in ("arithmetic", "background", "reflection"):
                    # C59 permits measured prefixes only at the terminal
                    # reflection boundary; all other valid requests are m=0.
                    measured_outputs = ((0, 0),)
                    if stop == circuit.width and boundary == "reflection":
                        measured_outputs += tuple(
                            (measured, output)
                            for measured in range(1, circuit.width + 1)
                            for output in range(1 << measured))
                    for measured, output in measured_outputs:
                        reference = circuit.prefix_vector(
                            final, exponent, stop, boundary=boundary,
                            measured=measured, output=output)
                        independent = direct_prefix(
                            circuit, final, exponent, stop, boundary,
                            measured, output)
                        with patch.object(CoherentReflectionCircuit, "_histories",
                                          side_effect=AssertionError("hidden history expansion")):
                            got = merged.prefix_vector(
                                final, exponent, stop, boundary=boundary,
                                measured=measured, output=output)
                        if len(got) != circuit.b:
                            raise AssertionError(f"{name}: wrong prefix dimension")
                        merged_max = max(merged_max,
                                         float(np.max(np.abs(got-reference))))
                        direct_max = max(direct_max,
                                         float(np.max(np.abs(reference-independent))))
                        zero_count += int(np.count_nonzero(reference == 0))
                        near_zero_count += int(np.count_nonzero(np.abs(reference) < 1e-14))
                        query_count += 1
                        prefix_stats = getattr(merged, "last_prefix_stats", None)
                        if prefix_stats is None or not stats_no_history(prefix_stats):
                            raise AssertionError(f"{name}: prefix retained history/tables")
                        if not (1 <= prefix_stats["peak_sector_count"] <=
                                prefix_stats["prefix_support_bound"] <= stats["support_bound"]):
                            raise AssertionError(f"{name}: inconsistent prefix support bound")

                        # Wrong-boundary control: W at a nontrivial intermediate
                        # stop is deliberately compared with the arithmetic one.
                        if (stop == 1 and boundary == "background"
                                and circuit.width >= 1):
                            arithmetic = circuit.prefix_vector(
                                final, exponent, stop, boundary="arithmetic")
                            wrong_boundary_delta = max(
                                wrong_boundary_delta,
                                float(np.max(np.abs(got-arithmetic))))
                        # Wrong-tail control: remove the factor for unvisited
                        # controls from a nonterminal arithmetic prefix.
                        if stop == 1 and boundary == "arithmetic" and circuit.width > 1:
                            wrong = got * (2.0 ** ((circuit.width-stop) / 2))
                            wrong_tail_delta = max(
                                wrong_tail_delta,
                                float(np.max(np.abs(wrong-reference))))
                        # Wrong-phase control: compare a nonzero terminal output
                        # with the same contraction after suppressing its phase.
                        if (stop == circuit.width and boundary == "reflection"
                                and measured == circuit.width and output == 1):
                            no_phase = direct_prefix(
                                circuit, final, exponent, stop, boundary,
                                measured, 0)
                            wrong_phase_delta = max(
                                wrong_phase_delta,
                                float(np.max(np.abs(got-no_phase))))

    prefix_stats = merged.stats()
    prefix_count_ok = (prefix_stats.get("merged_prefix_queries") == query_count
                       and prefix_stats.get("history_enumerations", 0) == 0)
    shape = (circuit.sectors, 1 << circuit.width)
    guard(shape, np.float64, f"{name} law")
    target = np.zeros(shape, dtype=float)
    existing_joint = np.zeros_like(target)
    from experiments.experiment_coherent_route_sampling import direct_joint as full_joint
    independent_joint = full_joint(circuit)
    for final in range(circuit.sectors):
        for output in range(1 << circuit.width):
            with patch.object(CoherentReflectionCircuit, "_histories",
                              side_effect=AssertionError("hidden history expansion")):
                target[final, output] = float(merged.joint_probability(final, output))
            existing_joint[final, output] = float(circuit.joint_probability(final, output))
    target_mass = float(target.sum())
    existing_mass = float(existing_joint.sum())
    with patch.object(CoherentReflectionCircuit, "_histories",
                      side_effect=AssertionError("hidden history expansion")):
        sample = merged.sample(np.random.default_rng(sample_seed))
    sample_stats = merged.stats()
    if not stats_no_history(sample_stats):
        raise AssertionError(f"{name}: sample stats retain history/tables")
    sample_ok = (0 <= int(sample["output"]) < (1 << circuit.width)
                 and 0 <= int(sample["final_coarse_sector"]) < circuit.sectors
                 and sample["merged_prefix_queries"] == sample["prefix_vector_evaluations"]
                 and sample["rejection_proposals"] == 0)
    return {
        "name": name,
        "period": circuit.period,
        "block_size": circuit.b,
        "width": circuit.width,
        "prefix_queries": query_count,
        "max_merged_vs_c59": merged_max,
        "max_c59_vs_direct": direct_max,
        "exact_zero_coordinates": zero_count,
        "near_zero_coordinates": near_zero_count,
        "target_mass": target_mass,
        "existing_mass": existing_mass,
        "target_mass_error": abs(target_mass - 1.0),
        "existing_mass_error": abs(existing_mass - 1.0),
        "joint_tv": float(np.sum(np.abs(target-existing_joint)) / 2),
        "independent_joint_tv": float(np.sum(np.abs(target-independent_joint)) / 2),
        "independent_mass_error": abs(float(independent_joint.sum())-1.),
        "wrong_boundary_delta": wrong_boundary_delta,
        "wrong_tail_delta": wrong_tail_delta,
        "wrong_phase_delta": wrong_phase_delta,
        "sample_ok": sample_ok,
        "prefix_count_ok": prefix_count_ok,
        "sample": {key: sample[key] for key in
                    ("output", "final_coarse_sector") if key in sample},
        "stats": prefix_stats,
    }


def main():
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "controls": {}}
    try:
        from lab.merged_prefix import MergedCoherentPrefixes
        fixtures = (("primary_r9_b3_t4", fixture_primary()),
                    ("endpoint_t0", fixture_t0()),
                    ("endpoint_t1", fixture_t1()),
                    ("wrap_fixed_b2_r10_t2", fixture_b2()),
                    ("pi_limit_r9_b3_t1", fixture_zero()),
                    ("exact_zero_r6_b2_t2", CoherentReflectionCircuit(6,2,2,{},{})))
        for index, (name, circuit) in enumerate(fixtures):
            report["rows"].append(run_fixture(
                name, circuit, MergedCoherentPrefixes, sample_seed=701 + index))
        p1 = all(row["max_merged_vs_c59"] <= TOL
                 and row["max_c59_vs_direct"] <= TOL
                 for row in report["rows"])
        p1 &= report["rows"][-1]["exact_zero_coordinates"] > 0
        p2 = all(row["joint_tv"] <= TOL
                 and row["independent_joint_tv"] <= TOL
                 and row["independent_mass_error"] <= TOL
                 and row["target_mass_error"] <= TOL
                 and row["existing_mass_error"] <= TOL
                 and row["prefix_count_ok"]
                 and row["sample_ok"] for row in report["rows"])
        p3 = all(row["max_merged_vs_c59"] <= TOL for row in report["rows"])
        c1 = max(row["wrong_boundary_delta"] for row in report["rows"]) > 1e-5
        c2 = max(row["wrong_tail_delta"] for row in report["rows"]) > 1e-5
        c3 = max(row["wrong_phase_delta"] for row in report["rows"]) > 1e-5
        report["controls"] = dict(
            wrong_boundary_delta=max(row["wrong_boundary_delta"] for row in report["rows"]),
            wrong_tail_delta=max(row["wrong_tail_delta"] for row in report["rows"]),
            wrong_phase_delta=max(row["wrong_phase_delta"] for row in report["rows"]),
            all_controls_fail=bool(c1 and c2 and c3))
        exp.check("P1", p1, "merged/C59/direct prefix amplitudes")
        exp.check("P2", p2, "joint laws, normalization, and table-free samples")
        exp.check("P3", p3, "edge, fixed-point, wrap, collision, and zero fixtures")
        exp.fail_check("C1", c1, f"wrong-boundary delta={report['controls']['wrong_boundary_delta']:.6g}")
        exp.fail_check("C2", c2, f"wrong-tail delta={report['controls']['wrong_tail_delta']:.6g}")
        exp.fail_check("C3", c3, f"wrong-phase delta={report['controls']['wrong_phase_delta']:.6g}")
        report["status"] = "PASS" if all((p1, p2, p3, c1, c2, c3)) else "FAIL"
        path = stamp()
        if not exp.finish(report_path=path, rows=report["rows"], metadata={
                "allocation_bound_bytes": MAX_DENSE_BYTES,
                "tolerance": TOL,
                "method": "MergedCoherentPrefixes vs C59 prefix_vector/direct_prefix",
                "controls": report["controls"],
                "elapsed_seconds": time.perf_counter() - started}):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = stamp("merged_prefix_formulas_failure")
        path.write_text(json.dumps(report, indent=2, default=str) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
