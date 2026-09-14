"""Complete-law audit of the coherently merged gate-by-gate prefix oracle.

The merged-prefix implementation replaces C59's history-expanded prefix oracle
with the existing gate-by-gate sampler loop.  This experiment checks every
tiny joint (coarse-sector, output) cell against the independent full-r
reference and the existing transition-law enumerator.  During merged runs the
copied circuit's ``_histories`` method is replaced by a raising sentinel, so a
passing law cannot be explained by hidden history enumeration.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The main b=3 r=9,t=4 fixture and bounded t=0/1, b=2, and generic-q
      controls have complete normalized joint laws matching both references.
  P2  Merged sampling has zero rejection proposals and its reported prefix,
      block-product, QFT-pair, reflection, and peak-sector counters are
      equal an independently counted reached-set recurrence for every query.
  P3  Removing the noncommuting intermediate background blocks changes the
      joint law; the merged oracle retains those effects while remaining
      history-free.
  P4  Taking componentwise absolute values AFTER the full prefix contraction
      preserves the sampling law, unlike discarding phases inside it.
  C1  Wrong prefix magnitudes change the complete law.
  C2  Omitting reflection inclusion at the requested boundary changes the law.
  C3  Using a wrong internal QFT phase changes the law.

These are complex128 diagnostics, not finite-TV certificates. Dense reference
arrays are guarded and only tiny complete laws are retained.
"""
from __future__ import annotations

import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.merged_prefix import MergedCoherentPrefixes


MAX_BYTES = 16 << 20
TARGET_TOL = 5e-11

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all tiny merged/reference transition and joint laws normalize and agree")
exp.predict("P2", "merged sample counters show history-free zero-rejection execution")
exp.predict("P3", "noncommuting intermediate backgrounds change the complete joint law")
exp.predict("P4", "componentwise absolute values of completed oracle outputs preserve the law")
exp.must_fail("C1", "wrong prefix magnitudes preserve the complete law")
exp.must_fail("C2", "omitting reflection boundary inclusion preserves the complete law")
exp.must_fail("C3", "using a wrong internal QFT phase preserves the complete law")


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return Path("out") / f"merged_prefix_sampling_{stamp}.json"


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
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError("merged-prefix diagnostic payload exceeds 16 MiB")


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], complex)


def rz(theta):
    return np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)]).astype(complex)


def embed(matrix):
    out = np.eye(3, dtype=complex)
    out[:2, :2] = matrix
    return out


def main_fixture(*, intermediate=True):
    defects = {0: embed(rx(np.pi / 4))}
    if intermediate:
        defects.update({1: embed(rx(np.pi / 7)), 2: embed(rz(np.pi / 5)),
                        3: embed(rx(np.pi / 6)), 4: embed(rz(np.pi / 7))})
    reflections = {0: (0, np.pi / 4), 1: (1, np.pi / 4),
                   2: (0, np.pi / 4), 4: (1, np.pi / 4)}
    return CoherentReflectionCircuit(9, 3, 4, defects, reflections)


def case_specs():
    return [
        ("main_b3_collisions", main_fixture(), 9, 3, 4),
        ("t0_endpoint", CoherentReflectionCircuit(9, 3, 0,
             {0: embed(rx(np.pi / 4))}, {0: (0, np.pi / 4)}), 9, 3, 0),
        ("t1_endpoint", CoherentReflectionCircuit(9, 3, 1,
             {0: embed(rx(np.pi / 4)), 1: embed(rz(np.pi / 5))},
             {0: (0, np.pi / 4), 1: (1, np.pi / 4)}), 9, 3, 1),
        ("b2_control", CoherentReflectionCircuit(6, 2, 2,
             {0: np.eye(2, dtype=complex), 1: rx(np.pi / 7), 2: rz(np.pi / 5)},
             {0: (0, np.pi / 4), 2: (1, np.pi / 4)}), 6, 2, 2),
        ("generic_q_control", CoherentReflectionCircuit(9, 3, 3,
             {0: embed(rx(np.pi / 4)), 1: embed(rx(np.pi / 7)),
              2: embed(rz(np.pi / 5)), 3: embed(rx(np.pi / 6))},
             {0: (0, np.pi / 4), 1: (2, np.pi / 4), 2: (4, np.pi / 4)}), 9, 3, 3),
    ]


def forbid_histories(circuit):
    def fail(*args, **kwargs):
        raise AssertionError("merged-prefix execution called forbidden _histories")
    circuit._histories = fail


def independent_helpers():
    from experiments.experiment_coherent_route_sampling import direct_joint, enumerate_transitions
    return direct_joint, enumerate_transitions


def expected_prefix_cost(circuit, sector, exponent, stop, *, boundary="reflection",
                         measured=0, output=0):
    reached = {sector}
    cost = dict(merged_prefix_queries=1, prefix_control_steps=stop,
                background_block_products=0, control_block_products=0,
                fixed_control_block_scales=0, qft_matrix_pair_constructions=0,
                reflection_block_contributions=0, peak_sector_count=1)
    for s in range(stop,-1,-1):
        if s in circuit.reflections and (s < stop or boundary == "reflection"):
            cost["reflection_block_contributions"] += 2*len(reached)
            q = circuit.reflections[s][0]
            reached |= {(-a-q)%circuit.sectors for a in reached}
            cost["peak_sector_count"] = max(cost["peak_sector_count"],len(reached))
        if s < stop or boundary != "arithmetic":
            cost["background_block_products"] += len(reached)
        if s:
            if s-1 >= circuit.width-measured:
                cost["control_block_products"] += len(reached)
                cost["qft_matrix_pair_constructions"] += len(reached)
            else:
                cost["control_block_products"] += bool(exponent & (1 << (s-1)))*len(reached)
                cost["fixed_control_block_scales"] += len(reached)
    return cost


def complete_case(name, supplied, period, block, width):
    direct_joint, enumerate_transitions = independent_helpers()
    guard((period // block, 1 << width), label=f"{name} full law")
    reference = CoherentReflectionCircuit(period, block, width,
        dict(supplied.background.defects), dict(supplied.reflections))
    merged = MergedCoherentPrefixes(supplied)
    # Forbid hidden history enumeration on the merged-owned copy.
    forbid_histories(merged.circuit)
    target = direct_joint(reference)
    enum = enumerate_transitions(merged)
    merged_joint = np.array([[merged.joint_probability(gamma, y)
                              for y in range(1 << width)]
                             for gamma in range(period // block)])
    sample_cost = dict.fromkeys(expected_prefix_cost(merged,0,0,0),0)
    original_prefix = merged.prefix_vector
    def counted_prefix(*args, **kwargs):
        expected = expected_prefix_cost(merged,*args,**kwargs)
        for key,value in expected.items():
            sample_cost[key] = (max(sample_cost[key],value) if key == "peak_sector_count"
                                else sample_cost[key]+value)
        return original_prefix(*args,**kwargs)
    merged.prefix_vector = counted_prefix
    sample = merged.sample(np.random.default_rng(7000 + width + block))
    stats = merged.stats()
    return {
        "name": name, "period": period, "block": block, "width": width,
        "target_mass": float(target.sum()), "enum_mass": float(enum.sum()),
        "merged_joint_mass": float(merged_joint.sum()),
        "enum_max_abs_error": float(np.max(np.abs(enum - target))),
        "merged_joint_max_abs_error": float(np.max(np.abs(merged_joint - target))),
        "enum_l1_tv_diagnostic": float(np.sum(np.abs(enum - target)) / 2),
        "merged_joint_l1_tv_diagnostic": float(np.sum(np.abs(merged_joint - target)) / 2),
        "sample": {k: sample.get(k) for k in (
            "output", "final_coarse_sector", "prefix_vector_evaluations",
            "rejection_proposals", "resampled_blocks")},
        "stats": stats,
        "sample_cost_exact": all(stats[key] == value for key,value in sample_cost.items()),
        "histories_forbidden": True,
    }


class WrongMagnitudeMerged(MergedCoherentPrefixes):
    def prefix_vector(self, *args, **kwargs):
        result = super().prefix_vector(*args, **kwargs).copy()
        if kwargs.get("stop", args[2] if len(args) > 2 else 0) > 0:
            result[0] *= 1.3
        return result


class WrongBoundaryMerged(MergedCoherentPrefixes):
    def prefix_vector(self, sector, exponent, stop, **kwargs):
        # Partial-QFT queries require reflection boundary semantics by API;
        # perturb only unmeasured gate-boundary inclusion.
        if (kwargs.get("boundary", "reflection") == "reflection"
                and kwargs.get("measured", 0) == 0):
            kwargs["boundary"] = "background"
        return super().prefix_vector(sector, exponent, stop, **kwargs)


class WrongPhaseMerged(MergedCoherentPrefixes):
    def prefix_vector(self, *args, **kwargs):
        # Perturb the internal QFT phase used by the shared merged-prefix
        # contraction, then restore the module binding immediately.  This
        # exercises a wrong oracle phase without introducing a propagator.
        import lab.merged_prefix as merged_module
        original = merged_module.unit_phase
        merged_module.unit_phase = lambda numerator, denominator: (
            original(numerator, denominator) * np.exp(0.123j))
        try:
            return super().prefix_vector(*args, **kwargs)
        finally:
            merged_module.unit_phase = original


class AbsoluteOutputMerged(MergedCoherentPrefixes):
    def prefix_vector(self, *args, **kwargs):
        value = super().prefix_vector(*args, **kwargs)
        # Sampling only uses |amplitude|^2, so this is law-preserving. It does
        # NOT permit deleting the phases while the amplitude is being formed.
        return np.abs(value)


def control_error(cls):
    direct_joint, enumerate_transitions = independent_helpers()
    reference = main_fixture()
    target = direct_joint(reference)
    wrong = cls(main_fixture())
    forbid_histories(wrong.circuit)
    law = enumerate_transitions(wrong)
    return float(np.max(np.abs(law - target))), float(np.sum(np.abs(law - target)) / 2)


def main():
    report = {"status": "PASS", "cases": [], "controls": {}, "background_effect": {}}
    p1 = p2 = p3 = p4 = True
    started = time.perf_counter()
    try:
        guard_entries(256)
        for name, circuit, period, block, width in case_specs():
            row = complete_case(name, circuit, period, block, width)
            report["cases"].append(row)
            p1 &= all(row[key] < TARGET_TOL for key in (
                "enum_max_abs_error", "merged_joint_max_abs_error"))
            p1 &= all(abs(row[key] - 1) < TARGET_TOL for key in (
                "target_mass", "enum_mass", "merged_joint_mass"))
            stats = row["stats"]
            sample = row["sample"]
            p2 &= sample["rejection_proposals"] == 0
            p2 &= sample["prefix_vector_evaluations"] == stats["merged_prefix_queries"]
            p2 &= stats["history_enumerations"] == 0 and stats["orbit_table_entries"] == 0
            p2 &= stats["merged_prefix_queries"] > 0
            p2 &= row["sample_cost_exact"]
            p2 &= stats["peak_sector_count"] <= stats["support_bound"]
            p2 &= all(stats[key] == 0 for key in ("sector_table_entries","output_table_entries"))
            p2 &= 0 <= sample["output"] < (1 << width)
            p2 &= 0 <= sample["final_coarse_sector"] < period // block
            if width > 0:
                p2 &= stats["prefix_control_steps"] > 0 and stats["control_block_products"] > 0

        reference = main_fixture()
        no_intermediate = main_fixture(intermediate=False)
        direct_joint, _ = independent_helpers()
        target = direct_joint(reference)
        no_target = direct_joint(no_intermediate)
        background_error = float(np.max(np.abs(target - no_target)))
        report["background_effect"] = {"max_abs_joint_difference": background_error,
                                        "target_mass": float(target.sum()),
                                        "no_intermediate_mass": float(no_target.sum())}
        p3 &= background_error > 1e-5

        for label, cls in (("wrong_magnitude", WrongMagnitudeMerged),
                           ("wrong_boundary", WrongBoundaryMerged),
                           ("wrong_phase", WrongPhaseMerged)):
            error, tv = control_error(cls)
            report["controls"][label] = {"max_abs_error": error, "tv_diagnostic": tv}
        abs_error,abs_tv = control_error(AbsoluteOutputMerged)
        report["absolute_output_invariance"] = {"max_abs_error": abs_error, "tv_diagnostic": abs_tv}
        p4 &= abs_tv < TARGET_TOL
        exp.check("P1", p1, "complete merged/reference laws and masses")
        exp.check("P2", p2, "history-free merged sample/query counters")
        exp.check("P3", p3, f"intermediate background law difference {background_error:.6g}")
        exp.check("P4", p4, "absolute output amplitudes preserve squared magnitudes")
        exp.fail_check("C1", report["controls"]["wrong_magnitude"]["tv_diagnostic"] > 1e-5,
                       "wrong prefix magnitudes change the complete law")
        exp.fail_check("C2", report["controls"]["wrong_boundary"]["tv_diagnostic"] > 1e-5,
                       "wrong reflection boundary changes the complete law")
        exp.fail_check("C3", report["controls"]["wrong_phase"]["tv_diagnostic"] > 1e-5,
                       "wrong internal QFT phase changes the complete law")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        exp.check("P1", False, "exception before completion")
        exp.check("P2", False, "exception before completion")
        exp.check("P3", False, "exception before completion")
        exp.check("P4", False, "exception before completion")
        exp.fail_check("C1", False, "exception before control")
        exp.fail_check("C2", False, "exception before control")
        exp.fail_check("C3", False, "exception before control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report["status"] = "PASS" if all(ok for _, _, ok, _ in exp._results) else "FAIL"
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)],
                    metadata={"target_tolerance": TARGET_TOL,
                              "history_calls_forbidden": True,
                              "native_memory_measured": False,
                              "detail_report_is_first_row": True})
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
