"""Tiny weighted-component Cauchy--Schwarz sampler audit.

Question: does the frozen mass-weighted proposal formula reproduce the target
law while lowering the exact rejection envelope on a row whose components have
unequal masses?

The primary fixture has Q=8, common stride 2, singleton component e=0 with
gamma=sqrt(.9), and the progression e=2,4 with gamma=sqrt(.05).  The target is
computed two ways: an explicitly expanded 8-entry row and the independent
finite geometric component formula.  Existing ``progression_sample`` is then
called on every reachable component/output path.  Path probabilities are
computed independently by finite Fourier sums and prefix conditioning,
including the kernel's gcd lift; this is not a precomputed proposal table.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The expanded-row FFT and component formulas agree; target, old proposal,
      weighted proposal, and their accepted submeasures normalize as derived.
  P2  The unequal-mass row has B^2=8/5 for the weighted envelope, versus m=2
      for the old mean-m envelope; exact-zero components are ignored,
      and equal-mass components give the same proposal and normalized law.
  P3  Enumerated calls to the existing progression kernel agree with the
      independently derived path probabilities, including gcd lifts.
  P4  The opt-in production row loop agrees with those paths for both proposal
      modes, including a forced retry and an explicit cap failure.
  C1  Applying the OLD acceptance factor to the NEW sqrt-mass proposal must
      give a normalized law different from the target; the component Fourier
      laws are deliberately different.

The split-progression case is a representation control: splitting the labels
2,4 into two singleton components preserves the coherent target but changes
the proposal/envelope.  The opt-in production hook is exercised only on this
isolated conditional row, not a physical full-joint schedule.  This is a
bounded float algebra/RNG-path audit, not a precision certificate, timing
result, or novelty claim.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.work_first import LateWorkProgressions
from lab.fourier_sampling import (
    geometric_sum,
    interval_prefix_probability,
    progression_sample,
    unit_phase,
)


Q = 8
WIDTH = 3
MAX_BYTES = 1 << 20
MAX_COMPONENT_TERMS = 4096
MAX_ROOT_TERMS = 512
MAX_FFT_ENTRIES = 256
MAX_RNG_DECISIONS = 2048
MAX_MARGINAL_QUERIES = 4096
MAX_PATH_RECORDS = 128
MAX_PRODUCTION_TERMS = 4096
TOL = 3e-12


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"component_weighting_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def json_safe(value):
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    return value


def tv(left, right):
    return float(np.sum(np.abs(np.asarray(left) - np.asarray(right))) / 2.)


def valid_law(values):
    values = np.asarray(values, dtype=float)
    return bool(values.shape == (Q,) and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) <= TOL)


def preflight():
    """Guard all arrays and finite work before the first science operation."""
    # Four frozen fixtures coexist only as 8-entry laws; retain FFT, formula,
    # proposal and accepted vectors plus scalar diagnostics and path records.
    # component_laws retains target/proposal/acceptance vectors, the per-y
    # component values, and their JSON-safe copies for all four fixtures.
    formula_float_bytes = 12 * 4 * Q * np.dtype(np.float64).itemsize
    formula_complex_bytes = 4 * Q * 3 * np.dtype(np.complex128).itemsize
    expanded_fft_bytes = 4 * Q * 3 * np.dtype(np.complex128).itemsize
    path_bytes = MAX_PATH_RECORDS * 12 * 8
    production_trace_bytes = MAX_PATH_RECORDS * 20 * 8
    production_law_bytes = 4 * 2 * Q * np.dtype(np.float64).itemsize
    json_copy_reserve = 2 * (formula_float_bytes + formula_complex_bytes
                             + expanded_fft_bytes + path_bytes)
    scalar_bytes = 16_384
    payload = (formula_float_bytes + formula_complex_bytes
               + expanded_fft_bytes + path_bytes + production_trace_bytes
               + production_law_bytes + json_copy_reserve + scalar_bytes)
    if payload > MAX_BYTES:
        raise MemoryError(f"payload preflight {payload}>{MAX_BYTES}")
    fixture_component_counts = (2, 2, 2, 3)
    component_terms = sum(Q * count for count in fixture_component_counts)
    # Distinct normalized progression kernels are (count,stride)=(1,2) and
    # (2,2); cache each independently derived finite Fourier table once.
    root_terms = Q * (1 + 2)
    fft_entries = len(fixture_component_counts) * Q
    # Each reachable output path has WIDTH' interval decisions and, for stride
    # two, one gcd lift. This is a conservative call-bound, not an RNG table.
    rng_decisions = sum(Q * count * (2 + 1) for count in fixture_component_counts)
    marginal_queries = sum(Q * count * 4 for count in fixture_component_counts)
    if component_terms > MAX_COMPONENT_TERMS:
        raise MemoryError("component-term preflight exceeds cap")
    if root_terms > MAX_ROOT_TERMS:
        raise MemoryError("independent root-term preflight exceeds cap")
    if fft_entries > MAX_FFT_ENTRIES:
        raise MemoryError("FFT-entry preflight exceeds cap")
    if rng_decisions > MAX_RNG_DECISIONS:
        raise MemoryError("RNG-decision preflight exceeds cap")
    if marginal_queries > MAX_MARGINAL_QUERIES:
        raise MemoryError("marginal-query preflight exceeds cap")
    return {
        "max_numeric_payload_bytes": MAX_BYTES,
        "planned_numeric_payload_bytes": payload,
        "payload_components": {
            "formula_float_arrays": formula_float_bytes,
            "formula_complex_values": formula_complex_bytes,
            "expanded_fft_temporaries": expanded_fft_bytes,
            "path_record_reserve": path_bytes,
            "production_trace_reserve": production_trace_bytes,
            "production_law_vectors": production_law_bytes,
            "json_copy_reserve": json_copy_reserve,
            "scalar_reserve": scalar_bytes,
        },
        "component_term_cap": MAX_COMPONENT_TERMS,
        "root_term_cap": MAX_ROOT_TERMS,
        "fft_entry_cap": MAX_FFT_ENTRIES,
        "rng_decision_cap": MAX_RNG_DECISIONS,
        "marginal_query_cap": MAX_MARGINAL_QUERIES,
        "path_record_cap": MAX_PATH_RECORDS,
        "production_term_cap": MAX_PRODUCTION_TERMS,
        "planned_component_terms": component_terms,
        "planned_root_terms": root_terms,
        "planned_fft_entries": fft_entries,
        "planned_rng_decisions": rng_decisions,
        "planned_marginal_queries": marginal_queries,
    }


def component_f(c, y):
    return (complex(c["gamma"])
            * unit_phase(-int(y) * int(c["start"]), Q)
            * geometric_sum(int(c["count"]), -int(y) * int(c["stride"]), Q))


def active_components(components):
    return [c for c in components if complex(c["gamma"]) != 0j]


def expanded_row(components, counters, caps):
    if counters["expanded_terms"] + sum(int(c["count"]) for c in components) > caps["expanded_terms"]:
        raise MemoryError("expanded row term cap before allocation")
    row = np.zeros(Q, dtype=np.complex128)
    for c in components:
        for n in range(int(c["count"])):
            counters["expanded_terms"] += 1
            index = int(c["start"]) + int(c["stride"]) * n
            if not 0 <= index < Q:
                raise ValueError("synthetic progression leaves Q-register")
            row[index] += complex(c["gamma"])
    return row


def component_laws(components, counters, caps):
    active = active_components(components)
    z = math.fsum(int(c["count"]) * abs(complex(c["gamma"]))**2
                  for c in active)
    if z <= 0 or not math.isfinite(z):
        raise ValueError("empty synthetic row")
    target_num = np.zeros(Q, dtype=float)
    incoherent = np.zeros(Q, dtype=float)
    values = [[0j for _ in active] for _ in range(Q)]
    for y in range(Q):
        for index, c in enumerate(active):
            if counters["component_terms"] >= caps["component_terms"]:
                raise MemoryError("component formula term cap before evaluation")
            counters["component_terms"] += 1
            value = component_f(c, y)
            values[y][index] = value
        coherent = sum(values[y], 0j)
        target_num[y] = abs(coherent)**2
        incoherent[y] = math.fsum(abs(value)**2 for value in values[y])
    target = target_num / (Q * z)
    old_proposal = incoherent / (Q * z)
    lambdas = [int(c["count"]) * abs(complex(c["gamma"]))**2 / z
               for c in active]
    root_sum = math.fsum(math.sqrt(value) for value in lambdas)
    weighted_proposal = np.zeros(Q, dtype=float)
    weighted_s = np.zeros(Q, dtype=float)
    old_accept = np.zeros(Q, dtype=float)
    weighted_accept = np.zeros(Q, dtype=float)
    m = len(active)
    for y in range(Q):
        weighted_s[y] = math.fsum(
            abs(values[y][i])**2 / math.sqrt(lambdas[i])
            for i in range(m) if lambdas[i] > 0)
        weighted_proposal[y] = weighted_s[y] / (Q * z * root_sum)
        if incoherent[y] > 0:
            old_accept[y] = target_num[y] / (m * incoherent[y])
        if weighted_s[y] > 0:
            weighted_accept[y] = target_num[y] / (root_sum * weighted_s[y])
    return {
        "active": active,
        "z": z,
        "values": values,
        "target_num": target_num,
        "target": target,
        "old_proposal": old_proposal,
        "weighted_proposal": weighted_proposal,
        "lambdas": lambdas,
        "B": root_sum,
        "B2": root_sum * root_sum,
        "old_accept": old_accept,
        "weighted_accept": weighted_accept,
        "old_accepted": old_proposal * old_accept,
        "weighted_accepted": weighted_proposal * weighted_accept,
    }


class ScriptedPathRNG:
    """Minimal RNG script for one progression_sample decision path."""

    def __init__(self, random_values, lift):
        self.random_values = list(random_values)
        self.lift = int(lift)
        self.random_calls = 0
        self.integer_calls = 0

    def random(self):
        if not self.random_values:
            raise AssertionError("progression_sample requested an extra random draw")
        self.random_calls += 1
        value = float(self.random_values.pop(0))
        if not 0. <= value < 1.:
            raise AssertionError("scripted random value outside [0,1)")
        return value

    def integers(self, high):
        self.integer_calls += 1
        if high <= self.lift:
            raise AssertionError("invalid scripted gcd lift")
        return self.lift


def direct_progression_distribution(count, stride, counters, caps):
    """Independent finite Fourier law and low-bit conditional paths."""
    g = math.gcd(int(stride), Q)
    reduced_dim = Q // g
    beta = int(stride) // g
    inv = pow(beta, -1, reduced_dim) if reduced_dim > 1 else 0
    reduced = np.zeros(reduced_dim, dtype=float)
    for output in range(reduced_dim):
        amplitude_terms = []
        for n in range(int(count)):
            if counters["root_terms"] >= caps["root_terms"]:
                raise MemoryError("independent root-term cap before summand")
            counters["root_terms"] += 1
            amplitude_terms.append(np.exp(-2j * np.pi * output * n / reduced_dim))
        amplitude = sum(amplitude_terms) / math.sqrt(int(count))
        reduced[output] = abs(amplitude)**2 / reduced_dim
    paths = []
    bits = reduced_dim.bit_length() - 1
    for reduced_output, probability in enumerate(reduced):
        # A contiguous finite sum is exactly zero here iff the nonzero
        # Fourier character completes an integer number of cycles.  Do not
        # discard merely-small positive probabilities.
        if reduced_output != 0 and (int(count) * reduced_output) % reduced_dim == 0:
            probability = 0.
        if probability == 0.:
            continue
        prefix_probability = 1.
        prefix = 0
        random_values = []
        path_probability = 1.
        for bit in range(bits):
            child = []
            mask = (1 << (bit + 1)) - 1
            for candidate in (0, 1):
                candidate_prefix = prefix | (candidate << bit)
                marginal = sum(
                    reduced[z] for z in range(reduced_dim)
                    if (z & mask) == candidate_prefix)
                child.append(float(marginal / prefix_probability))
            chosen = (reduced_output >> bit) & 1
            branch = child[chosen]
            if branch <= 0.:
                raise AssertionError("independent path selected zero branch")
            random_values.append(0.5 * child[0] if chosen == 0
                                 else child[0] + 0.5 * child[1])
            path_probability *= branch
            prefix |= chosen << bit
            prefix_probability *= branch
        if abs(path_probability - probability) > 1e-12:
            raise AssertionError("prefix path probability does not telescope")
        for lift in range(g):
            output = ((inv * reduced_output) % reduced_dim
                      if reduced_dim > 1 else 0) + reduced_dim * lift
            paths.append({"reduced": reduced_output, "lift": lift,
                          "output": output,
                          "probability": probability / g,
                          "random_values": random_values,
                          "random_calls": bits,
                          "integer_calls": int(g > 1)})
    return paths


def enumerate_progression_paths(component, paths, counters, caps):
    if counters["path_records"] + len(paths) > caps["path_records"]:
        raise MemoryError("path-record cap before progression calls")
    observed = []
    for path in paths:
        planned = path["random_calls"] + path["integer_calls"]
        if counters["rng_decisions"] + planned > caps["rng_decisions"]:
            raise MemoryError("RNG decision cap before progression call")
        planned_queries = 2 * path["random_calls"]
        if counters["marginal_queries"] + planned_queries > caps["marginal_queries"]:
            raise MemoryError("marginal-query cap before progression call")
        rng = ScriptedPathRNG(path["random_values"], path["lift"])
        actual = progression_sample(
            WIDTH, int(component["start"]), int(component["stride"]),
            int(component["count"]), 0, 1, rng)
        counters["rng_decisions"] += rng.random_calls + rng.integer_calls
        counters["path_records"] += 1
        if actual["output"] != path["output"]:
            raise AssertionError("progression_sample path output mismatch")
        if rng.random_calls != path["random_calls"] or rng.integer_calls != path["integer_calls"]:
            raise AssertionError("progression_sample RNG decision count mismatch")
        if actual["marginal_queries"] != 2 * path["random_calls"]:
            raise AssertionError("progression marginal query count mismatch")
        counters["marginal_queries"] += int(actual["marginal_queries"])
        observed.append({"component_start": int(component["start"]),
                         "component_count": int(component["count"]),
                         "component_stride": int(component["stride"]),
                         "output": int(actual["output"]),
                         "independent_probability": path["probability"],
                         "random_calls": rng.random_calls,
                         "integer_calls": rng.integer_calls,
                         "marginal_queries": actual["marginal_queries"]})
    return observed


def ensure_production(total, increments, caps, label):
    for key, amount in increments.items():
        if total.get(key, 0) + int(amount) > caps[key]:
            raise MemoryError(
                f"{label} exceeds production {key} cap: "
                f"{total.get(key, 0)}+{amount}>{caps[key]}")


def add_counts(total, counters):
    for key, value in (counters or {}).items():
        if isinstance(value, (int, np.integer)):
            total[key] = total.get(key, 0) + int(value)


def production_fourier_audit(model, row, formula, production_totals, caps):
    """Compare every production proposal/acceptance coordinate to formulas."""
    results = {}
    m = len(row["components"])
    for mode in ("mass", "root_mass"):
        planned_prepare = {"component_weight_terms": m}
        if mode == "root_mass":
            planned_prepare["component_sqrt_terms"] = m
        ensure_production(production_totals, planned_prepare, caps,
                          f"{mode} proposal preparation")
        counters = model._new_counters()
        prepared = model._prepare_proposal(row, mode, counters)
        add_counts(production_totals, counters)
        proposal = np.zeros(Q, dtype=float)
        acceptance = np.zeros(Q, dtype=float)
        accepted = np.zeros(Q, dtype=float)
        for output in range(Q):
            increments = {"fourier_component_terms": m}
            if mode == "root_mass":
                increments["weighted_ratio_divisions"] = m
            ensure_production(production_totals, increments, caps,
                              f"{mode} Fourier output")
            before = dict(counters)
            observed = model._fourier(row, output, counters, prepared=prepared)
            delta = {key: value - before.get(key, 0)
                     for key, value in counters.items()
                     if isinstance(value, (int, np.integer))}
            add_counts(production_totals, delta)
            proposal[output] = float(observed["proposal_probability"])
            acceptance[output] = float(observed["acceptance"])
            accepted[output] = proposal[output] * acceptance[output]
        expected_proposal = (formula["old_proposal"] if mode == "mass"
                             else formula["weighted_proposal"])
        expected_acceptance = (formula["old_accept"] if mode == "mass"
                               else formula["weighted_accept"])
        if np.max(np.abs(proposal - expected_proposal)) > TOL:
            raise AssertionError(f"{mode} production proposal mismatch")
        if np.max(np.abs(acceptance - expected_acceptance)) > TOL:
            raise AssertionError(f"{mode} production acceptance mismatch")
        normalized = accepted / accepted.sum()
        if not valid_law(normalized) or np.max(np.abs(normalized - formula["target"])) > TOL:
            raise AssertionError(f"{mode} production normalized accepted law mismatch")
        results[mode] = {
            "proposal": proposal,
            "acceptance": acceptance,
            "accepted": accepted,
            "normalized_accepted": normalized,
            "expected_attempts": float(prepared["expected_attempts"]),
            "counters": dict(counters),
        }
    return results


def component_selector_random(weights, index):
    total = math.fsum(weights)
    before = math.fsum(weights[:index])
    return (before + 0.5 * weights[index]) / total


def production_row_call(model, row, formula, mode, component_index, path,
                        production_totals, caps, *, accept, max_attempts=1,
                        partial_retry=False):
    """Force one actual _sample_row path, retaining counters on failure."""
    raw_weights = [int(c["count"]) * abs(complex(c["gamma"]))**2
                   for c in formula["active"]]
    if mode == "root_mass":
        raw_weights = [math.sqrt(weight) for weight in raw_weights]
    component_random = component_selector_random(raw_weights, component_index)
    output = int(path["output"])
    alpha = float((formula["old_accept"] if mode == "mass"
                   else formula["weighted_accept"])[output])
    if partial_retry:
        if max_attempts != 2 or not 0. < alpha < 1.:
            raise AssertionError("partial retry requires two nontrivial attempts")
        acceptance_randoms = [(alpha + 1.) / 2., alpha / 2.]
    elif accept:
        if alpha <= 0:
            raise AssertionError("cannot force acceptance at zero probability")
        acceptance_randoms = [alpha / 2.] * max_attempts
    else:
        acceptance_randoms = [(alpha + 1.) / 2.] * max_attempts
        if acceptance_randoms[0] >= 1.:
            raise AssertionError("cannot force rejection of certain acceptance")
    random_values = []
    for _attempt in range(max_attempts):
        random_values.append(component_random)
        random_values.extend(path["random_values"])
        random_values.append(acceptance_randoms[_attempt])
    planned = {
        "component_weight_terms": len(row["components"]),
        "fourier_component_terms": len(row["components"]) * max_attempts,
        "progression_proposals": max_attempts,
        "marginal_queries": 2 * path["random_calls"] * max_attempts,
        "component_draws": max_attempts,
        "acceptance_draws": max_attempts,
        "production_calls": 1,
        "rng_random_draws": (1 + path["random_calls"] + 1) * max_attempts,
        "rng_integer_draws": path["integer_calls"] * max_attempts,
    }
    if mode == "root_mass":
        planned["component_sqrt_terms"] = len(row["components"])
        planned["weighted_ratio_divisions"] = len(row["components"]) * max_attempts
    ensure_production(production_totals, planned, caps,
                      f"{mode} production row call")
    production_totals["production_calls"] += 1
    rng = ScriptedPathRNG(random_values, path["lift"])
    counters = model._new_counters()
    try:
        result = model._sample_row(row, rng, counters,
                                   max_attempts=max_attempts, proposal=mode)
        failed = False
    except RuntimeError:
        result = None
        failed = True
    add_counts(production_totals, counters)
    production_totals["rng_random_draws"] += rng.random_calls
    production_totals["rng_integer_draws"] += rng.integer_calls
    if rng.random_calls != len(random_values):
        raise AssertionError("production row consumed an unexpected random path")
    if rng.integer_calls != path["integer_calls"] * max_attempts:
        raise AssertionError("production row consumed an unexpected lift path")
    if not failed:
        if result["output"] != output or result["proposal_mode"] != mode:
            raise AssertionError("production row returned the wrong forced path")
        expected = (len(row["components"]) if mode == "mass" else formula["B2"])
        if abs(float(result["expected_attempts"]) - expected) > TOL:
            raise AssertionError("production expected-attempt field mismatch")
    return {
        "mode": mode, "component_index": component_index,
        "output": output, "independent_probability": float(path["probability"]),
        "acceptance": alpha, "accepted": not failed,
        "attempts": int(result["attempts"]) if result else max_attempts,
        "counters": dict(counters),
    }


def main():
    exp = Experiment("component_weighting", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "expanded FFT and weighted component formulas agree")
    exp.predict("P2", "weighted B2 and normalized laws match the derivation")
    exp.predict("P3", "existing progression RNG paths match independent probabilities")
    exp.predict("P4", "production row loop matches both proposal modes and cap paths")
    exp.must_fail("C1", "new proposal with old acceptance is not target law")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = p4 = c1 = False
    try:
        report["preflight"] = preflight()
        caps = {"expanded_terms": 128, "component_terms": MAX_COMPONENT_TERMS,
                "root_terms": MAX_ROOT_TERMS,
                "fft_entries": MAX_FFT_ENTRIES,
                "rng_decisions": MAX_RNG_DECISIONS,
                "marginal_queries": MAX_MARGINAL_QUERIES,
                "path_records": MAX_PATH_RECORDS}
        production_caps = {
            "production_calls": 128,
            "component_weight_terms": MAX_PRODUCTION_TERMS,
            "component_sqrt_terms": MAX_PRODUCTION_TERMS,
            "weighted_ratio_divisions": MAX_PRODUCTION_TERMS,
            "fourier_component_terms": MAX_PRODUCTION_TERMS,
            "progression_proposals": 128,
            "marginal_queries": MAX_MARGINAL_QUERIES,
            "component_draws": 128,
            "acceptance_draws": 128,
            "rng_random_draws": 512,
            "rng_integer_draws": 128,
        }
        counters = defaultdict(int)
        fixtures = {
            "unequal": [
                {"start": 0, "count": 1, "stride": 2, "gamma": math.sqrt(.9)},
                {"start": 2, "count": 2, "stride": 2, "gamma": math.sqrt(.05)},
            ],
            "unequal_with_exact_zero": [
                {"start": 0, "count": 1, "stride": 2, "gamma": math.sqrt(.9)},
                {"start": 2, "count": 2, "stride": 2, "gamma": math.sqrt(.05)},
                {"start": 6, "count": 1, "stride": 2, "gamma": 0j},
            ],
            "equal": [
                {"start": 0, "count": 1, "stride": 2, "gamma": 1/math.sqrt(2)},
                {"start": 2, "count": 1, "stride": 2, "gamma": 1/math.sqrt(2)},
            ],
            "split": [
                {"start": 0, "count": 1, "stride": 2, "gamma": math.sqrt(.9)},
                {"start": 2, "count": 1, "stride": 2, "gamma": math.sqrt(.05)},
                {"start": 4, "count": 1, "stride": 2, "gamma": math.sqrt(.05)},
            ],
        }
        laws = {}
        path_rows = {}
        fft_errors = {}
        path_errors = {}
        progression_cache = {}
        for name, components in fixtures.items():
            # Expanded row allocation is guarded by the total payload preflight
            # and again immediately before the finite support loop.
            expanded = expanded_row(components, counters, caps)
            if counters["fft_entries"] + Q > caps["fft_entries"]:
                raise MemoryError("FFT entry cap before transform")
            counters["fft_entries"] += Q
            fft_law = np.abs(np.fft.fft(expanded))**2 / (Q * np.sum(np.abs(expanded)**2))
            formula = component_laws(components, counters, caps)
            fft_errors[name] = float(np.max(np.abs(fft_law - formula["target"])))
            if fft_errors[name] > TOL:
                raise AssertionError(f"expanded FFT mismatch for {name}")
            active_path_rows = []
            actual_by_output = np.zeros(Q, dtype=float)
            for component in formula["active"]:
                cache_key = (int(component["count"]), int(component["stride"]))
                if cache_key not in progression_cache:
                    progression_cache[cache_key] = direct_progression_distribution(
                        *cache_key, counters, caps)
                component_paths = enumerate_progression_paths(
                    component, progression_cache[cache_key], counters, caps)
                active_path_rows.extend(component_paths)
                weight = (int(component["count"])
                          * abs(complex(component["gamma"]))**2 / formula["z"])
                for path in component_paths:
                    actual_by_output[path["output"]] += weight * path["independent_probability"]
            path_rows[name] = active_path_rows
            expected_component_mixture = np.zeros(Q, dtype=float)
            for component in formula["active"]:
                weight = (int(component["count"]) * abs(complex(component["gamma"]))**2
                          / formula["z"])
                cache_key = (int(component["count"]), int(component["stride"]))
                for path in progression_cache[cache_key]:
                    expected_component_mixture[path["output"]] += weight * path["probability"]
            path_errors[name] = float(np.max(np.abs(actual_by_output
                                                   - expected_component_mixture)))
            path_errors[name] = max(path_errors[name], float(np.max(
                np.abs(actual_by_output - formula["old_proposal"]))))
            if path_errors[name] > TOL:
                raise AssertionError(f"path mixture mismatch for {name}")
            laws[name] = formula

        unequal = laws["unequal"]
        zero = laws["unequal_with_exact_zero"]
        equal = laws["equal"]
        split = laws["split"]
        # Exact-zero support is semantically discarded, not treated as a
        # fallback component or a division-by-zero proposal weight.
        zero_target_error = float(np.max(np.abs(zero["target"] - unequal["target"])))
        equal_proposal_error = tv(equal["old_proposal"], equal["weighted_proposal"])
        equal_accept_error = tv(
            equal["weighted_accepted"] / equal["weighted_accepted"].sum(), equal["target"])
        split_target_error = float(np.max(np.abs(split["target"] - unequal["target"])))
        split_proposal_tv = tv(split["weighted_proposal"], unequal["weighted_proposal"])
        wrong_accepted = unequal["weighted_proposal"] * unequal["old_accept"]
        wrong_normalized = wrong_accepted / wrong_accepted.sum()
        wrong_tv = tv(wrong_normalized, unequal["target"])
        p1 = bool(all(fft_errors[name] <= TOL for name in fixtures)
                  and all(valid_law(laws[name]["target"]) for name in fixtures)
                  and all(valid_law(laws[name]["old_proposal"]) for name in fixtures)
                  and all(valid_law(laws[name]["weighted_proposal"]) for name in fixtures)
                  and all(valid_law(laws[name]["old_accepted"]
                                    / laws[name]["old_accepted"].sum())
                          for name in fixtures)
                  and all(valid_law(laws[name]["weighted_accepted"]
                                    / laws[name]["weighted_accepted"].sum())
                          for name in fixtures)
                  and all(np.max(np.abs(laws[name][key]
                                        / laws[name][key].sum()
                                        - laws[name]["target"])) <= TOL
                          for name in fixtures
                          for key in ("old_accepted", "weighted_accepted"))
                  and all(abs(float(laws[name]["old_accepted"].sum())
                              - 1/len(laws[name]["active"])) < 1e-12
                          for name in fixtures)
                  and abs(float(unequal["weighted_accepted"].sum())
                          - 1/unequal["B2"]) < 1e-12)
        p2 = bool(abs(unequal["B2"] - 8/5) < 1e-12
                  and abs(unequal["old_accepted"].sum() - .5) < 1e-12
                  and zero_target_error <= TOL
                  and equal_proposal_error <= TOL
                  and equal_accept_error <= TOL
                  and split_target_error <= TOL
                  and split_proposal_tv > 1e-8
                  and abs(float(split["weighted_accepted"].sum())
                          - 1/split["B2"]) < 1e-12)
        p3 = bool(all(path_errors[name] <= TOL for name in fixtures)
                  and counters["path_records"] <= MAX_PATH_RECORDS
                  and counters["rng_decisions"] <= MAX_RNG_DECISIONS)
        # Exercise the actual opt-in conditional production loop on the
        # synthetic row, without constructing a full physical circuit.
        production_model = LateWorkProgressions(
            8, 1, WIDTH, 0, np.ones((1, 1), dtype=np.complex128),
            np.ones((1, 1), dtype=np.complex128), lambda _index: 1.+0j,
            max_local_terms=1_000_000, max_components=4096,
            max_payload_bytes=MAX_BYTES)
        production_row = {
            "work": 0,
            "components": tuple((int(c["start"]), int(c["count"]),
                                  complex(c["gamma"]))
                                 for c in unequal["active"]),
            "stride": 2, "norm": float(unequal["z"]),
        }
        production_totals = defaultdict(int)
        production_fourier_unequal = production_fourier_audit(
            production_model, production_row, unequal,
            production_totals, production_caps)
        equal_row = {
            "work": 0,
            "components": tuple((int(c["start"]), int(c["count"]),
                                  complex(c["gamma"]))
                                 for c in equal["active"]),
            "stride": 2, "norm": float(equal["z"]),
        }
        production_fourier_equal = production_fourier_audit(
            production_model, equal_row, equal,
            production_totals, production_caps)
        production_fourier = {"unequal": production_fourier_unequal,
                              "equal": production_fourier_equal}
        production_path_laws = {}
        production_traces = {}
        production_rejection_traces = {}
        production_accepted_submasses = {}
        production_rejected_masses = {}
        production_path_errors = {}
        cap_failure_seen = {}
        for mode in ("mass", "root_mass"):
            actual_proposal = np.zeros(Q, dtype=float)
            traces = []
            rejection_traces = []
            accepted_submass = np.zeros(Q, dtype=float)
            rejected_mass = 0.
            mode_weights = [int(c["count"]) * abs(complex(c["gamma"]))**2
                            for c in unequal["active"]]
            if mode == "root_mass":
                mode_weights = [math.sqrt(value) for value in mode_weights]
            mode_total = math.fsum(mode_weights)
            for component_index, component in enumerate(unequal["active"]):
                cache_key = (int(component["count"]), int(component["stride"]))
                component_weight = mode_weights[component_index] / mode_total
                for path in progression_cache[cache_key]:
                    alpha = float((unequal["old_accept"] if mode == "mass"
                                   else unequal["weighted_accept"])[path["output"]])
                    trace = production_row_call(
                        production_model, production_row, unequal, mode,
                        component_index, path, production_totals, production_caps,
                        accept=alpha > 0.)
                    traces.append(trace)
                    actual_proposal[path["output"]] += (
                        component_weight * path["probability"])
                    if trace["accepted"] != (alpha > 0.):
                        raise AssertionError("production acceptance branch mismatch")
                    accepted_submass[path["output"]] += (
                        component_weight * path["probability"] * alpha)
                    rejected_mass += (component_weight * path["probability"]
                                      * (1. - alpha))
                    if alpha < 1.:
                        rejection = production_row_call(
                            production_model, production_row, unequal, mode,
                            component_index, path, production_totals, production_caps,
                            accept=False)
                        rejection_traces.append(rejection)
                        if rejection["accepted"]:
                            raise AssertionError("forced rejection branch unexpectedly accepted")
            production_path_laws[mode] = actual_proposal
            production_traces[mode] = traces
            production_rejection_traces[mode] = rejection_traces
            production_accepted_submasses[mode] = accepted_submass
            production_rejected_masses[mode] = rejected_mass
            expected_proposal = (unequal["old_proposal"] if mode == "mass"
                                 else unequal["weighted_proposal"])
            production_path_errors[mode] = float(np.max(
                np.abs(actual_proposal - expected_proposal)))
            cap_failure_seen[mode] = bool(rejection_traces) and all(
                not trace["accepted"] for trace in rejection_traces)

        retry_case = None
        for component_index, component in enumerate(unequal["active"]):
            cache_key = (int(component["count"]), int(component["stride"]))
            for path in progression_cache[cache_key]:
                if 0. < float(unequal["old_accept"][path["output"]]) < 1.:
                    retry_case = (component_index, path)
                    break
            if retry_case is not None:
                break
        if retry_case is None:
            raise AssertionError("no nontrivial acceptance path for retry control")
        retry_traces = {}
        for mode in ("mass", "root_mass"):
            component_index, path = retry_case
            alpha = float((unequal["old_accept"] if mode == "mass"
                           else unequal["weighted_accept"])[path["output"]])
            if not 0. < alpha < 1.:
                continue
            retry_traces[mode] = production_row_call(
                production_model, production_row, unequal, mode,
                component_index, path, production_totals, production_caps,
                accept=True, max_attempts=2, partial_retry=True)
            if retry_traces[mode]["attempts"] != 2:
                raise AssertionError("partial retry did not consume two attempts")
        p4 = bool(all(production_path_errors[mode] <= TOL
                      for mode in ("mass", "root_mass"))
                  and all(valid_law(production_fourier[fixture][mode]["normalized_accepted"])
                          and np.max(np.abs(
                              production_fourier[fixture][mode]["normalized_accepted"]
                              - laws[fixture]["target"])) <= TOL
                          for fixture in ("unequal", "equal")
                          for mode in ("mass", "root_mass"))
                  and all(cap_failure_seen.values())
                  and all(np.max(np.abs(
                      production_accepted_submasses[mode]
                      - (unequal["old_accepted"] if mode == "mass"
                         else unequal["weighted_accepted"]))) <= TOL
                          and abs(production_rejected_masses[mode]
                                  - (1. - 1./(len(unequal["active"])
                                             if mode == "mass" else unequal["B2"]))) <= TOL
                          for mode in ("mass", "root_mass"))
                  and set(retry_traces) == {"mass", "root_mass"})
        c1 = bool(wrong_tv > 1e-6)
        p1 = bool(p1 and all(valid_law(production_fourier[fixture][mode][
            "normalized_accepted"])
                             for fixture in ("unequal", "equal")
                             for mode in ("mass", "root_mass")))
        report.update({
            "fixtures": {name: components for name, components in fixtures.items()},
            "fft_max_errors": fft_errors,
            "path_max_errors": path_errors,
            "laws": {
                name: {key: value for key, value in law.items()
                       if key not in ("active", "values")}
                for name, law in laws.items()},
            "path_rows": path_rows,
            "controls": {
                "zero_target_error": zero_target_error,
                "equal_proposal_tv": equal_proposal_error,
                "equal_normalized_accepted_tv": equal_accept_error,
                "split_target_error": split_target_error,
                "split_weighted_proposal_tv": split_proposal_tv,
                "wrong_new_proposal_old_acceptance_tv": wrong_tv,
                "production_path_errors": production_path_errors,
                "production_cap_failure_seen": cap_failure_seen,
                "production_rejected_masses": production_rejected_masses,
                "production_retry_traces": retry_traces,
            },
            "expected_attempts": {
                name: {"old_m": len(law["active"]),
                       "old": float(len(law["active"])),
                       "weighted_B2": float(law["B2"])}
                for name, law in laws.items()},
            "counters": dict(counters),
            "production_fourier": production_fourier,
            "production_path_laws": production_path_laws,
            "production_traces": production_traces,
            "production_rejection_traces": production_rejection_traces,
            "production_accepted_submasses": production_accepted_submasses,
            "production_counters": dict(production_totals),
            "checks": {"P1": p1, "P2": p2, "P3": p3, "P4": p4, "C1": c1},
            "status": "PASS" if p1 and p2 and p3 and p4 and c1 else "FAIL",
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
            "production_conditional_loop_executed": True,
            "no_physical_joint_fixture": True,
        })
        exp.check("P1", p1, "expanded FFT/formula laws and mass normalization")
        exp.check("P2", p2, "weighted envelope, zero/equal/split controls")
        exp.check("P3", p3, "independent progression RNG paths")
        exp.check("P4", p4, "actual production row modes, retry, and cap failure")
        exp.fail_check("C1", c1, "wrong acceptance on new sqrt-mass proposal")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3", "P4"):
            exp.check(name, False, "exception before component audit")
        exp.fail_check("C1", False, "exception before wrong-acceptance control")
    report["elapsed_seconds"] = time.perf_counter() - started
    report = json_safe(report)
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": "Q=8 common stride2 weighted component row",
        "reference": "expanded row FFT + independent finite path probabilities",
        "production_conditional_loop_executed": True,
        "no_physical_joint_fixture": True,
        "max_numeric_payload_bytes": MAX_BYTES,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
