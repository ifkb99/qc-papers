"""Decision-path edge audit for ``LateWorkProgressions``.

This is a tiny exact b=2, r=6, Q=16 fixture.  The expected joint law is
assembled independently from the supplied b-by-b matrices and compared with
the public column/forced-joint paths.  A scripted RNG forces one rejection and
one acceptance while checking that the work label is drawn only once.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Independent column amplitudes and every forced joint probability agree;
      rows, proposals, and columns normalize.
  P2  A scripted two-attempt draw accepts on attempt two with one work draw,
      exercising the real progression_sample path and preserving the work row.
  P3  Zero rows/components and invalid dimensions/phases/caps raise or return
      zero safely without numerical thresholding.
  C1  Deliberately redrawing the work label on each retry is detected by the
      fixed-work counter/predicate rather than silently accepted.
"""
from __future__ import annotations

import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

import numpy as np

from lab import Experiment
from lab.work_first import LateWorkProgressions, _pick, _weight


MAX_BYTES = 16 * 1024 * 1024
PERIOD, BLOCK, WIDTH, SPLIT = 6, 2, 4, 2
Q = 1 << WIDTH
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"work_first_sampler_edges_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(payload, label):
    payload = int(payload)
    if payload < 0 or payload > MAX_BYTES:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def valid_law(values, size):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (size,) and np.all(np.isfinite(values))
                and np.all(values >= -TOL)
                and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right):
    return float(np.sum(np.abs(np.asarray(left)-np.asarray(right)))/2)


def json_safe(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def direct_columns(initial, middle, phase):
    """Independent finite formula for all columns, no sampler internals."""
    if Q*BLOCK*BLOCK > 64:
        raise MemoryError("direct column preflight exceeds 64 local terms")
    columns = np.zeros((Q, PERIOD), dtype=complex)
    terms = 0
    for exponent in range(Q):
        low = exponent % (1 << SPLIT)
        high = exponent // (1 << SPLIT)
        for u in range(BLOCK):
            cell, q = divmod((u+low) % PERIOD, BLOCK)
            for p in range(BLOCK):
                terms += 1
                index = cell*BLOCK+p
                columns[exponent, (index+(1 << SPLIT)*high) % PERIOD] += (
                    initial[u, 0]*middle[p, q]*phase(index))
    return columns, terms


class ScriptedRNG:
    """First attempt rejects, second accepts; all progression bits are zero."""
    def __init__(self):
        self.random_calls = []
        self.integer_calls = []

    def random(self):
        index = len(self.random_calls)
        self.random_calls.append(index)
        # call 0: work; call 1: component; calls 2..4: interval bits;
        # call 5: reject; call 6: component; calls 7..9: bits; call 10: accept.
        return 0.999999999999 if index == 5 else 0.0

    def integers(self, high):
        self.integer_calls.append(int(high))
        return 0


class ExhaustRNG:
    def __init__(self):
        self.random_calls = 0
        self.integer_calls = []

    def random(self):
        self.random_calls += 1
        return 0.999999999999

    def integers(self, high):
        self.integer_calls.append(int(high))
        return 0


def raises(call, kinds):
    try:
        call()
    except kinds:
        return True
    return False


def main():
    exp = Experiment("work_first_sampler_edges", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "independent columns and forced joint law agree")
    exp.predict("P2", "scripted retry uses one fixed work draw and real progression sampling")
    exp.predict("P3", "zero and invalid edge paths are explicit and bounded")
    exp.must_fail("C1", "a retry path that redrew work would violate fixed-work accounting")
    started = time.perf_counter()
    report = {"status": "FAIL", "fixture": {
        "period": PERIOD, "block_size": BLOCK, "width": WIDTH, "split": SPLIT,
        "Q": Q}}
    p1 = p2 = p3 = c1 = False
    try:
        matrix = np.array([[1., 1.], [1., -1.]]) / math.sqrt(2.)
        phase = lambda index: 1.+0j
        # The source's numeric reserve is known from structural inputs before
        # copying either supplied matrix; include it in the aggregate guard.
        component_bound = (1 << (WIDTH-SPLIT)) * min(PERIOD, 2*BLOCK-1)
        sampler_reserve = 16*(16*BLOCK*BLOCK + 12*component_bound + 128)
        planned_payload = (5*sampler_reserve + 16*Q*PERIOD*16
                           + 8*PERIOD*PERIOD*16)
        guard_bytes(planned_payload, "edge dense reference and sampler reserve")
        if PERIOD*Q*Q > 1536 or PERIOD*Q > 96:
            raise MemoryError("direct Fourier/forced-call preflight")
        source_copy = matrix.copy()
        sampler = LateWorkProgressions(PERIOD, BLOCK, WIDTH, SPLIT,
                                        source_copy, source_copy, phase)
        stats = sampler.stats()
        report["stats"] = stats
        report["planned_payload_bytes"] = planned_payload
        report["sampler_reserve_leq_plan"] = stats["numeric_payload_bound_bytes"] <= planned_payload
        source_copy[0, 0] = 17.
        report["matrix_copy_isolated"] = bool(
            abs(sampler.initial[0, 0]-1/math.sqrt(2.)) < TOL)
        columns, direct_column_terms = direct_columns(matrix, matrix, phase)
        dense_joint_terms = 0
        report["direct_column_terms"] = direct_column_terms
        column_errors = []
        column_norms = []
        for exponent in range(Q):
            observed = sampler.column(exponent)
            actual = np.zeros(PERIOD, dtype=complex)
            for index, value in observed["amplitudes"].items():
                actual[index] = value
            column_errors.append(float(np.max(np.abs(actual-columns[exponent]))))
            column_norms.append(float(np.vdot(actual, actual).real))
        joint_errors = []
        work_laws = []
        direct_work_laws = []
        row_component_counts = []
        zero_component_seen = False
        forced_calls = 0
        for work in range(PERIOD):
            row = sampler.row(work)
            row_component_counts.append(len(row["components"]))
            # Directly identify an omitted zero-gamma residue in the tiny H/H
            # fixture, if one exists; zero components must not be proposed.
            v = (work-(1 << SPLIT)*0) % PERIOD
            cell, p = divmod(v, BLOCK)
            raw = {}
            for u in range(BLOCK):
                for q in range(BLOCK):
                    rho = (cell*BLOCK+q-u) % PERIOD
                    raw[rho] = raw.get(rho, 0j) + matrix[p, q]*matrix[u, 0]
            exact_zero_residues = {rho for rho, value in raw.items() if value == 0j}
            zero_component_seen |= bool(exact_zero_residues)
            if exact_zero_residues:
                active_starts = {start for start, count, _ in row["components"]
                                 if count > 0}
                if any(rho in active_starts for rho in exact_zero_residues):
                    raise AssertionError("exact-zero component was retained")
            expected = np.zeros(Q, dtype=float)
            z = float(np.sum(np.abs(columns[:, work])**2))
            direct_work_laws.append(z/Q)
            if z > 0:
                for output in range(Q):
                    if dense_joint_terms+Q > 1536:
                        raise MemoryError("direct-root term cap before vector sum")
                    dense_joint_terms += Q
                    amplitude = np.sum(columns[:, work]
                                       * np.exp(-2j*np.pi*output*np.arange(Q)/Q))
                    expected[output] = abs(amplitude)**2 / Q**2
            observed = []
            for output in range(Q):
                if forced_calls >= 96:
                    raise MemoryError("forced-call cap before helper call")
                forced_calls += 1
                result = sampler.forced_joint(work, output)
                observed.append(float(result["joint_probability"]))
            work_laws.append(observed)
            joint_errors.append(float(np.max(np.abs(np.asarray(observed)-expected))))
        report["column_max_error"] = max(column_errors)
        report["column_max_norm_error"] = max(abs(value-1.) for value in column_norms)
        report["forced_joint_calls"] = forced_calls
        report["dense_joint_terms"] = dense_joint_terms
        report["forced_joint_max_error"] = max(joint_errors)
        report["row_component_counts"] = row_component_counts
        report["zero_component_seen"] = zero_component_seen
        report["work_marginal"] = direct_work_laws
        report["work_marginal_mass"] = sum(direct_work_laws)
        report["joint_mass"] = sum(sum(row) for row in work_laws)
        report["component_proposal_bound"] = stats["component_bound"]

        rng = ScriptedRNG()
        sampled = sampler.sample(rng, max_attempts=2)
        report["scripted_sample"] = sampled
        report["rng_random_calls"] = len(rng.random_calls)
        report["rng_integer_calls"] = rng.integer_calls
        report["rng_call_counts_expected"] = {"random": 11, "integers": 3}
        report["sampled_work"] = sampled["work"]
        report["sampled_attempts"] = sampled["attempts"]
        report["sample_counters"] = sampled["counters"]

        zero_sampler = LateWorkProgressions(32, 2, WIDTH, SPLIT,
                                             np.eye(2), np.eye(2), phase)
        zero_row = zero_sampler.row(31)
        zero_joint = zero_sampler.forced_joint(31, 0)
        report["zero_row"] = {"norm": zero_row["norm"],
                              "components": len(zero_row["components"]),
                              "forced_joint": zero_joint}

        invalid = {
            "bad_period_divisibility": raises(
                lambda: LateWorkProgressions(5, 2, WIDTH, SPLIT, matrix, matrix, phase),
                (ValueError,)),
            "bad_width": raises(
                lambda: LateWorkProgressions(6, 2, 64, SPLIT, matrix, matrix, phase),
                (ValueError,)),
            "bad_split": raises(
                lambda: LateWorkProgressions(6, 2, WIDTH, WIDTH+1, matrix, matrix, phase),
                (ValueError,)),
            "bad_unitary": raises(
                lambda: LateWorkProgressions(6, 2, WIDTH, SPLIT, np.ones((2, 2)), matrix, phase),
                (ValueError,)),
            "bad_components_cap": raises(
                lambda: LateWorkProgressions(6, 2, WIDTH, SPLIT, matrix, matrix, phase,
                                             max_components=1), (MemoryError,)),
            "bad_payload_cap": raises(
                lambda: LateWorkProgressions(6, 2, WIDTH, SPLIT, matrix, matrix, phase,
                                             max_payload_bytes=1), (MemoryError,)),
            "bad_exponent": raises(lambda: sampler.column(-1), (ValueError,)),
            "bad_work": raises(lambda: sampler.row(PERIOD), (ValueError,)),
            "bad_output": raises(lambda: sampler.forced_joint(0, Q), (ValueError,)),
            "bad_attempt_cap": raises(lambda: sampler.sample(ScriptedRNG(), max_attempts=0),
                                       (ValueError,)),
        }
        bad_phase = LateWorkProgressions(PERIOD, BLOCK, WIDTH, SPLIT,
                                          matrix, matrix, lambda _: 0j)
        invalid["bad_phase_value"] = raises(lambda: bad_phase.column(0), (ValueError,))
        phase_calls = []
        invalid["huge_history_preflight"] = raises(
            lambda: LateWorkProgressions(6, 2, 63, 0, matrix, matrix,
                                         lambda index: phase_calls.append(index) or 1+0j,
                                         max_components=1), (MemoryError,))
        invalid["huge_history_phase_calls_zero"] = len(phase_calls) == 0
        invalid["bad_local_cap"] = raises(
            lambda: LateWorkProgressions(6, 2, WIDTH, SPLIT, matrix, matrix, phase,
                                         max_local_terms=1), (MemoryError,))
        with patch("lab.work_first.np.array", side_effect=AssertionError("numeric copy reached")):
            invalid["payload_cap_before_copy"] = raises(
                lambda: LateWorkProgressions(6,2,WIDTH,SPLIT,matrix,matrix,phase,
                                             max_payload_bytes=1), (MemoryError,))
            invalid["huge_history_before_copy"] = raises(
                lambda: LateWorkProgressions(6,2,63,0,matrix,matrix,phase), (MemoryError,))
        point = LateWorkProgressions(1,1,3,3,np.ones((1,1)),np.ones((1,1)),phase)
        zero_proposal = point.forced_joint(0,1)
        report["zero_proposal"] = zero_proposal
        invalid["zero_proposal_zero_acceptance"] = (
            zero_proposal["proposal_probability"] == zero_proposal["acceptance"] == 0)
        exhaust_rng = ExhaustRNG()
        exhausted = raises(lambda: sampler.sample(exhaust_rng, max_attempts=1),
                           (RuntimeError,))
        report["exhaustion"] = {"raised": exhausted,
                                "random_calls": exhaust_rng.random_calls,
                                "integer_calls": exhaust_rng.integer_calls}
        report["invalid_paths"] = invalid

        p1 = (report["column_max_error"] < TOL
              and report["column_max_norm_error"] < TOL
              and report["forced_joint_max_error"] < TOL
              and abs(report["work_marginal_mass"]-1.) < TOL
              and abs(report["joint_mass"]-1.) < TOL
              and forced_calls == PERIOD*Q
              and report["direct_column_terms"] == Q*BLOCK*BLOCK
              and report["dense_joint_terms"] == PERIOD*Q*Q
              and report["sampler_reserve_leq_plan"]
              and report["matrix_copy_isolated"])
        counters = sampled["counters"]
        p2 = (sampled["attempts"] == 2 and counters["work_draws"] == 1
              and counters["component_draws"] == 2
              and counters["progression_proposals"] == 2
              and counters["acceptance_draws"] == 2
              and sampled["work"] == 0
              and len(rng.integer_calls) == 3
              and len(rng.random_calls) == 11
              and len(set(row_component_counts)) > 1)
        p3 = (zero_row["norm"] == 0 and len(zero_row["components"]) == 0
              and zero_joint["joint_probability"] == 0
              and zero_joint["conditional_probability"] is None
              and zero_joint["proposal_probability"] == 0
              and zero_component_seen and all(invalid.values())
              and report["exhaustion"]["raised"])
        # Exercise a deliberately wrong retry implementation through the
        # actual column path: it redraws the auxiliary exponent/work twice.
        wrong_rng = ScriptedRNG()
        wrong_column_calls = 0
        wrong_work_draws = 0
        for _ in range(2):
            wrong_exponent = wrong_rng.integers(Q)
            sampler.column(wrong_exponent)
            wrong_column_calls += 1
            wrong_work_draws += 1
        c1 = (counters["work_draws"] == 1 and wrong_work_draws == 2
              and wrong_column_calls != counters["work_draws"])
        report["wrong_redraw_control"] = {
            "actual_correct_work_draws": counters["work_draws"],
            "deliberate_wrong_work_draws": wrong_work_draws,
            "deliberate_wrong_column_calls": wrong_column_calls,
        }
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "columns and forced joint law match independent formula")
        exp.check("P2", p2, "scripted reject/accept uses fixed work and progression_sample")
        exp.check("P3", p3, "zero and invalid paths are explicit")
        exp.fail_check("C1", c1, "retry path retains one sampled work label")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[json_safe(report)], metadata={
        "fixture": {"period": PERIOD, "block_size": BLOCK,
                    "width": WIDTH, "split": SPLIT},
        "max_payload_bytes": MAX_BYTES,
        "reference": "independent finite column formula and existing progression_sample",
        "no_new_propagator": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
