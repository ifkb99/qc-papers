"""Bounded end-to-end audit of the TODO41 nested-phase sampler.

This is the first experiment that calls the released ``NestedPhaseProgressions``
helper.  The complete tiny laws come from the existing ER literal-column
reference and FFTs; the helper rows are expanded only one row at a time.  The
frozen r=60 family uses the ER ``direct_column`` multi-insertion path for
empty, initial, duplicate, endpoint, and three-phase schedules.  Tiny edge
fixtures use the existing arbitrary-period literal engine, with an analytical
folding for the r=6 two-oracle revival.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Selected auto and explicit cuts reproduce every streamed row amplitude,
      row norm, work marginal, and complete normalized FFT joint law.
  P2  Empty/initial/duplicate/end schedules, the r=b alias, non-divisor
      period, exact-zero rows, and the revived cancellation row retain their
      literal laws; cache_right does not change the law.
  P3  For both mass and root_mass, the forced proposal probability multiplied
      by its acceptance and the documented envelope factor reproduces the full
      conditional law (and therefore the complete joint law).
  C1  Dropping the final early phase from the named multiple-phase schedule
      changes a literal amplitude and its FFT law.
  C2  The preserved NC near-cancellation control has tiny amplitude error but
      an unstable conditional law after normalization; no positive cutoff is
      allowed to hide it.

Run:  uv run python -m experiments.experiment_nested_phase_sampler

This is a bounded complex128/float64 integration audit.  It does not certify
finite-bit sampling or RNG frequencies; the larger returned-sample RNG audit
remains a subsequent TODO41 check.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.fourier_sampling import unit_phase


MAX_BYTES = 16 << 20
MAX_WORK = 2_000_000
TOL = 3e-10


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_sampler_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def safe(value):
    if isinstance(value, complex):
        return [float(value.real), float(value.imag)]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [safe(v) for v in value]
    return value


class Budget:
    """Frozen cumulative named-work ledger, charged before each operation."""

    def __init__(self, cap):
        self.cap = int(cap)
        self.used = defaultdict(int)
        self.observed = defaultdict(int)

    def reserve(self, name, amount):
        amount = int(amount)
        if amount < 0 or sum(self.used.values()) + amount > self.cap:
            raise MemoryError(f"{name} exceeds frozen cumulative work cap")
        self.used[name] += amount

    def observe(self, name, amount):
        amount = int(amount)
        self.observed[name] += amount


def counter_delta(before, after):
    keys = set(before) | set(after)
    return {key: int(after.get(key, 0)) - int(before.get(key, 0))
            for key in keys if int(after.get(key, 0)) != int(before.get(key, 0))}


def reserve_operation(budget, kind, limits):
    for key, amount in limits.items():
        budget.reserve(f"{kind}:{key}", amount)


def check_operation_counters(budget, kind, before, counters, limits):
    delta = counter_delta(before, counters)
    for key, value in delta.items():
        if value < 0 or value > int(limits.get(key, 0)):
            raise AssertionError(f"{kind} counter {key}={value} exceeds reserve "
                                 f"{limits.get(key, 0)}")
        budget.observe(f"{kind}:{key}", value)
    return delta


def reconcile_budget(budget):
    """Require every observed counter category to fit its prior reserve."""
    result = {}
    for reserved_key, actual in budget.observed.items():
        reserved = int(budget.used.get(reserved_key, 0))
        result[reserved_key] = {"reserved": reserved, "observed": int(actual),
                                "ok": int(actual) <= reserved}
        if actual > reserved:
            raise AssertionError(f"observed {reserved_key} exceeds reserve")
    return result


def guard_bytes(amount, label):
    if int(amount) > MAX_BYTES:
        raise MemoryError(f"{label} {amount}>{MAX_BYTES}")


def preflight():
    """Freeze setup, references, rows, laws, controls, and report reserve."""
    # The largest retained object is one Q-by-r reference plus one Q-row and
    # two FFT vectors.  Rows and references are streamed between calls.
    fixtures = {
        "r60_empty": (60, 3, 4, 3),
        "r60_initial": (60, 3, 4, 3),
        "r60_duplicate": (60, 3, 4, 3),
        "r60_endpoint": (60, 3, 4, 3),
        "r60_multiple": (60, 3, 4, 3),
        "r60_zero": (60, 3, 4, 3),
        "r9_nondivisor": (9, 3, 4, 3),
        "r6_revival": (6, 2, 2, 2),
        "r3_alias": (3, 3, 3, 2),
    }
    retained = max((1 << t) * r * 16 for r, _b, t, _s in fixtures.values())
    payload = retained + 2 * max(1 << t for _r, _b, t, _s in fixtures.values()) * 16
    payload += 2 * 64 * 64 * 16 + 256 * 1024
    guard_bytes(payload, "nested sampler aggregate numeric payload")

    # Reserve every named category before execution.  Helper calls are row
    # streamed, while the full law is built from each expanded row and never
    # from a Q-by-r helper archive.
    reference = sum((1 << t) * (r * (b * b + 6))
                    for r, b, t, _s in fixtures.values())
    constructors = 2 * len(fixtures) + 8
    rows = 520_000
    accepted = 520_000
    forced = 20_000
    fft = sum(r * (1 << t) * 3 for r, _b, t, _s in fixtures.values())
    fft_input_entries = sum(r * (1 << t) * 2 * (6 if name == "r60_multiple" else 1)
                            for name, (r, _b, t, _s) in fixtures.items()) + 2 * 16 * 60
    selector_setup = sum((1 << (t-s)) * b*b +
                         (1 << (t-s)) * min(r, 2*b-1) + 64
                         for r, b, t, s in fixtures.values())
    # At most four phase schedule entries are inspected per high history in
    # the frozen families; oracle calls are charged separately by operations.
    geometry_setup_phase_terms = sum((1 << (t-s)) * 4
                                     for _r, _b, t, s in fixtures.values())
    controls = 4096
    phase_metadata = {
        "r60_empty": (), "r60_initial": (0,), "r60_duplicate": (2, 2),
        "r60_endpoint": (3,), "r60_multiple": (0, 2, 3),
        "r60_zero": (), "r9_nondivisor": (1, 1),
        "r6_revival": (1, 2), "r3_alias": (),
    }
    total = (reference + constructors * 2 * 64 + rows + accepted + forced
             + fft + fft_input_entries + selector_setup
             + geometry_setup_phase_terms + controls)
    if total > MAX_WORK:
        raise MemoryError(f"nested sampler preflight {total}>{MAX_WORK}")
    return {
        "fixtures": fixtures,
        "phase_schedule_metadata": {
            name: {"positions": list(phase_metadata[name]),
                   "candidate_cut_bound": len(set(phase_metadata[name]) | {s})}
            for name, (_r, _b, _t, s) in fixtures.items()
        },
        "planned_numeric_payload_bytes": payload,
        "payload_components": {
            "largest_reference": retained,
            "largest_row_and_fft": 2 * max(1 << t for _r, _b, t, _s in fixtures.values()) * 16,
            "matrix_and_scratch": 2 * 64 * 64 * 16,
            "scalar_reserve": 256 * 1024,
        },
        "reserved_reference_work": reference,
        "reserved_reference_setup_work": sum(2 * b**3 for _r, b, _t, _s in fixtures.values()),
        "reserved_control_reference_work": 2 * 3**3 + 8 * 3 * (3 * 3 + 6),
        "reserved_constructor_work": constructors * 2 * 64,
        "reserved_row_work": rows,
        "reserved_accepted_law_work": accepted,
        "reserved_forced_work": forced,
        "reserved_fft_work": fft,
        "reserved_fft_input_entries": fft_input_entries,
        "reserved_selector_setup_bound": selector_setup,
        "reserved_geometry_setup_phase_terms": geometry_setup_phase_terms,
        "reserved_control_work": controls,
        "reserved_total_work": total,
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_total_work": MAX_WORK,
        "failed_report_reserve": 1,
    }


def h2():
    return np.asarray([[1., 1.], [1., -1.]], dtype=np.complex128) / math.sqrt(2.)


def identity_phase(_index):
    return 1. + 0j


def parity(index):
    return -1. + 0j if int(index) & 1 else 1. + 0j


def f2(index):
    return (1j) ** (int(index) & 1)


def g60(index):
    return unit_phase(pow(2, int(index) % 60, 61), 61)


def g9(index):
    if not 0 <= int(index) < 9:
        raise AssertionError("r9 oracle received a noncanonical index")
    return unit_phase(int(index), 9)


def f3(index):
    return unit_phase(int(index), 3)


def operation_limits(stats, sampler, kind, q):
    """Structural per-call reserves, before invoking the production method."""
    b = int(sampler.b)
    if kind == "column":
        return {"column_local_terms": b*b,
                "phase_queries": int(stats["phase_queries_per_draw_bound"]),
                "early_phase_queries": int(stats.get("column_phase_products_bound", 0)),
                "late_phase_queries": b*b,
                "column_phase_products": int(stats["column_phase_products_bound"])}
    if kind == "row":
        return {
            "selector_pair_terms": int(stats.get("selector_pair_bound", 0)),
            "selector_residue_terms": int(stats.get("selector_residue_bound", 0)),
            "selector_candidate_terms": int(stats.get("selector_candidate_bound", 0)),
            "phase_partition_terms": len(stats.get("phase_positions", ())),
            "row_local_terms": int(stats.get("row_pair_bound", stats["local_term_bound"])),
            "group_pair_visits": int(stats.get("row_pair_bound", 0)),
            "coefficient_pair_visits": int(stats.get("coefficient_pair_bound", 0)),
            "geometric_components": int(stats["component_bound"]),
            "phase_queries": int(stats["phase_queries_per_draw_bound"]),
            "early_phase_queries": int(stats.get("phase_queries_per_draw_bound", 0)),
            "late_phase_queries": int(stats.get("phase_queries_per_draw_bound", 0)),
            "row_phase_products": int(stats.get("row_phase_products_bound", 0)),
            "norm_terms": int(stats["component_bound"]),
            "right_cache_entries": int(stats.get("right_cache_live_entries_bound", 0)
                                        * max(1, int(sampler.L))
                                        * int(sampler.H)),
        }
    if kind == "accepted":
        return {"fourier_component_terms": int(stats["component_bound"]),
                "component_weight_terms": int(stats["component_bound"]),
                "component_sqrt_terms": int(stats["component_bound"]),
                "weighted_ratio_divisions": int(stats["component_bound"]),
                "phase_queries": 0,
                "row_phase_products": 0}
    if kind == "forced":
        limits = operation_limits(stats, sampler, "row", q)
        for key, value in operation_limits(stats, sampler, "accepted", q).items():
            limits[key] = max(limits.get(key, 0), value)
        return limits
    raise ValueError(f"unknown operation kind {kind}")


def expand(row, q, budget):
    components = tuple(row.get("components", ()))
    stride = int(row["stride"])
    values = np.zeros(int(q), dtype=np.complex128)
    seen = set()
    for item in components:
        if len(item) != 3:
            raise AssertionError("row component is not a three-tuple")
        start, count = map(int, item[:2])
        amplitude = complex(item[2])
        if count < 1:
            raise AssertionError("row component has nonpositive count")
        budget.reserve("row_expansion_terms", count)
        for offset in range(count):
            exponent = start + stride * offset
            if not 0 <= exponent < q or exponent in seen:
                raise AssertionError("row components overlap or leave register")
            seen.add(exponent)
            values[exponent] = amplitude
    return values


def law(values):
    values = np.asarray(values, dtype=np.complex128)
    probabilities = np.abs(values) ** 2
    total = float(probabilities.sum())
    if not np.isfinite(total) or total < 0:
        raise AssertionError("nonfinite law mass")
    return (None, total) if total == 0. else (probabilities / total, total)


def reference_columns(name, *, period, block, width, split, initial, middle,
                      phase, positions, budget, early_phases=None):
    """Use existing literal engines; never construct a second propagator."""
    q = 1 << width
    counters = defaultdict(int)
    if period == 60 and block == 3:
        from experiments.experiment_earlier_phase_cycles import direct_column
        columns = np.empty((q, period), dtype=np.complex128)
        for exponent in range(q):
            limits = {"direct_block_products": (period // block) * block * block,
                      "direct_phase_queries": period * (len(positions) + 1),
                      "modular_pow_queries": period * (len(positions) + 1),
                      "direct_shift_ops": len(positions) + 3}
            reserve_operation(budget, "reference", limits)
            before = dict(counters)
            columns[exponent] = direct_column(
                exponent, positions[0] if positions else 0, initial, middle,
                counters, split=split, early_insertions=positions)
            check_operation_counters(budget, "reference", before, counters, limits)
        return columns
    from experiments.experiment_earlier_phase_sampler_edges import literal_columns
    # Arbitrary tiny periods are covered through the existing literal engine.
    if name == "r6_revival":
        # In this r=6 fixture parity(u+l mod 2) times i^(u+l mod 4)
        # equals (-i)^(u+l mod 2), so two phases fold analytically to one
        # phase at the split endpoint.
        folded = lambda index: parity(index) * f2(index)
        early_split, early = split, folded
    elif name == "r9_nondivisor":
        early_split, early = 1, lambda index: g9(index) * g9(index)
    elif positions:
        early_split, early = positions[0], g9
    else:
        early_split, early = 0, identity_phase
    limits = {"column_calls": q, "column_local_terms": q * (period // block) * block * block,
              "early_phase_queries": q * block, "late_phase_queries": q * period,
              "shift_operations": q}
    reserve_operation(budget, "reference_literal", limits)
    edge_ref = __import__("experiments.experiment_earlier_phase_sampler_edges",
                          fromlist=["REFERENCE_LEDGER"])
    old_ledger, old_bounds = edge_ref.REFERENCE_LEDGER, edge_ref.REFERENCE_BOUNDS
    actual = defaultdict(int)
    edge_ref.REFERENCE_LEDGER, edge_ref.REFERENCE_BOUNDS = actual, limits
    try:
        columns = literal_columns(period, block, width, split, initial, middle,
                                  phase, early_split, early)
    finally:
        edge_ref.REFERENCE_LEDGER, edge_ref.REFERENCE_BOUNDS = old_ledger, old_bounds
    for key, value in actual.items():
        if value > limits.get(key, 0):
            raise AssertionError(f"reference literal {key} exceeds frozen reserve")
        budget.observe(f"reference_literal:{key}", value)
    return columns


def fixture_specs():
    w0 = np.asarray(np.eye(3), dtype=np.complex128)
    w1 = np.asarray(np.eye(3), dtype=np.complex128)
    # The work-block family is the same frozen C79 physical fixture used by ER.
    from experiments.experiment_clean_orbit_output import work_block
    c0, c1 = work_block(math.pi / 4), work_block(-math.pi / 10)
    h = h2()
    return {
        "r60_empty": dict(period=60, block=3, width=4, split=3,
            initial=c0, middle=c1, phase=g60, positions=(), early_phases=()),
        "r60_initial": dict(period=60, block=3, width=4, split=3,
            initial=c0, middle=c1, phase=g60, positions=(0,),
            early_phases=((0, g60),)),
        "r60_duplicate": dict(period=60, block=3, width=4, split=3,
            initial=c0, middle=c1, phase=g60, positions=(2, 2),
            early_phases=((2, g60), (2, g60))),
        "r60_endpoint": dict(period=60, block=3, width=4, split=3,
            initial=c0, middle=c1, phase=g60, positions=(3,),
            early_phases=((3, g60),)),
        "r60_multiple": dict(period=60, block=3, width=4, split=3,
            initial=c0, middle=c1, phase=g60, positions=(0, 2, 3),
            early_phases=((0, g60), (2, g60), (3, g60))),
        "r60_zero": dict(period=60, block=3, width=4, split=3,
            initial=w0, middle=w1, phase=g60, positions=(), early_phases=()),
        "r9_nondivisor": dict(period=9, block=3, width=4, split=3,
            initial=w0, middle=w1, phase=identity_phase, positions=(1, 1),
            early_phases=((1, g9), (1, g9))),
        "r6_revival": dict(period=6, block=2, width=2, split=2,
            initial=h, middle=h, phase=identity_phase, positions=(1, 2),
            early_phases=((1, parity), (2, f2))),
        "r3_alias": dict(period=3, block=3, width=3, split=2,
            initial=w0, middle=w1, phase=identity_phase, positions=(),
            early_phases=()),
    }


def make_sampler(cls, spec, *, cut, cache_right, budget):
    budget.reserve("constructor", 2 * 64)
    return cls(spec["period"], spec["block"], spec["width"], spec["split"],
               spec["initial"], spec["middle"], spec["phase"],
               early_phases=spec["early_phases"], cut=cut,
               cache_right=cache_right, max_local_terms=250_000,
               max_components=4096, max_payload_bytes=MAX_BYTES)


def verify_columns(sampler, spec, reference, budget):
    q, period = 1 << spec["width"], spec["period"]
    stats = sampler.stats()
    max_error = max_norm_error = 0.
    for exponent in range(q):
        limits = operation_limits(stats, sampler, "column", q)
        reserve_operation(budget, "column", limits)
        before = {}
        result = sampler.column(exponent)
        counters = result["counters"]
        check_operation_counters(budget, "column", before, counters, limits)
        values = np.zeros(period, dtype=np.complex128)
        for index, amplitude in result["amplitudes"].items():
            if not 0 <= int(index) < period:
                raise AssertionError("column amplitude leaves physical work register")
            values[int(index)] = complex(amplitude)
        max_error = max(max_error, float(np.max(np.abs(values-reference[exponent]))))
        max_norm_error = max(max_norm_error,
                             abs(float(np.sum(np.abs(values)**2))-1.))
    return {"max_column_amplitude_error": max_error,
            "max_column_norm_error": max_norm_error}


def compare_rows(name, spec, sampler, reference, budget, *, accepted_modes=True):
    q = 1 << spec["width"]
    period = spec["period"]
    row_amplitude_error = row_norm_error = 0.
    joint_error = 0.
    joint = np.zeros((q, period), dtype=float)
    selected_cuts = set()
    zero_rows = []
    component_counts = []
    revival_amplitude = None
    accepted_laws = {"mass": np.zeros_like(joint),
                     "root_mass": np.zeros_like(joint)}
    for work in range(period):
        stats = sampler.stats()
        row_limits = operation_limits(stats, sampler, "row", q)
        reserve_operation(budget, "row", row_limits)
        before = {}
        row = sampler.row(work)
        check_operation_counters(budget, "row", before, row["counters"], row_limits)
        component_counts.append(len(row.get("components", ())))
        if not isinstance(row.get("cut"), (int, np.integer)):
            raise AssertionError(f"{name} row does not expose selected cut")
        selected_cuts.add(int(row["cut"]))
        expanded = expand(row, q, budget)
        if name == "r6_revival" and work == 1:
            revival_amplitude = complex(expanded[0])
        expected = reference[:, work]
        row_amplitude_error = max(row_amplitude_error,
                                  float(np.max(np.abs(expanded - expected))))
        expected_norm = float(np.sum(np.abs(expected) ** 2))
        row_norm = float(np.sum(np.abs(expanded) ** 2))
        reported_norm = float(row["norm"])
        row_norm_error = max(row_norm_error, abs(row_norm - expected_norm),
                             abs(row_norm - reported_norm))
        if row_norm == 0.:
            zero_rows.append(work)
        budget.reserve("fft_input_entries", 2 * q)
        spectrum = np.fft.fft(expanded) / q
        joint[:, work] = np.abs(spectrum) ** 2
        reference_spectrum = np.fft.fft(expected) / q
        joint_error = max(joint_error,
                          float(np.max(np.abs(joint[:, work]
                                              - np.abs(reference_spectrum) ** 2))))
        if accepted_modes:
            for mode in ("mass", "root_mass"):
                accepted_limits = operation_limits(stats, sampler, "accepted", q)
                accepted_limits["fourier_component_terms"] *= q
                accepted_limits["weighted_ratio_divisions"] *= q
                if mode == "mass":
                    accepted_limits["component_sqrt_terms"] = 0
                    accepted_limits["weighted_ratio_divisions"] = 0
                reserve_operation(budget, "accepted", accepted_limits)
                before = dict(row["counters"])
                prepared = sampler._prepare_proposal(row, mode, row["counters"])
                proposal_total = accepted_total = 0.
                for output in range(q):
                    result = sampler._fourier(row, output, row["counters"],
                                              prepared=prepared)
                    proposal = float(result["proposal_probability"])
                    acceptance = float(result["acceptance"])
                    proposal_total += proposal
                    accepted_total += proposal * acceptance
                    factor = (float(result["expected_attempts"])
                              if mode == "root_mass"
                              else float(len(row.get("components", ()))))
                    accepted_laws[mode][output, work] = factor * proposal * acceptance
                    target = (float(q * np.abs(reference_spectrum[output])**2 / expected_norm)
                              if expected_norm else 0.)
                    if abs(accepted_laws[mode][output, work] - target) > TOL:
                        raise AssertionError(f"{name} {mode} q*acceptance mismatch")
                check_operation_counters(budget, "accepted", before,
                                         row["counters"], accepted_limits)
                if abs(proposal_total - (1. if row_norm else 0.)) > TOL:
                    raise AssertionError(f"{name} {mode} proposal law does not normalize")
                expected_inverse = (1./float(prepared["expected_attempts"])
                                    if prepared["expected_attempts"] else 0.)
                if abs(accepted_total - expected_inverse) > TOL:
                    raise AssertionError(f"{name} {mode} accepted mass mismatch")
                if work == 0:
                    forced_limits = operation_limits(stats, sampler, "forced", q)
                    reserve_operation(budget, "forced", forced_limits)
                    forced = sampler.forced_joint(work, 0, proposal=mode)
                    forced_target = float(np.abs(reference_spectrum[0])**2)
                    if abs(float(forced["joint_probability"]) - forced_target) > TOL:
                        raise AssertionError(f"{name} public forced_joint mismatch")
                    check_operation_counters(budget, "forced", {},
                                             forced["counters"], forced_limits)
    accepted_errors = None
    if accepted_modes:
        denominator = joint.sum(axis=0, keepdims=True)
        conditional = np.divide(joint, denominator,
                                out=np.zeros_like(joint), where=denominator != 0)
        accepted_errors = {
            mode: float(np.max(np.abs(values - conditional)))
            for mode, values in accepted_laws.items()
        }
    # The accepted law is checked against forced conditional values above; its
    # joint sum is the same law after multiplying by the row work marginal.
    return {
        "max_row_amplitude_error": row_amplitude_error,
        "max_row_norm_error": row_norm_error,
        "max_joint_probability_error": joint_error,
        "work_marginal_total": float(joint.sum()),
        "selected_cuts": sorted(selected_cuts),
        "zero_rows": zero_rows,
        "accepted_q_acceptance_max_error": accepted_errors,
        "component_counts": component_counts,
        "revival_work1_exponent0": (None if revival_amplitude is None
                                     else [float(revival_amplitude.real),
                                           float(revival_amplitude.imag)]),
    }, joint


def near_cancellation_control(cls, budget):
    """Reproduce NC's retained tiny-positive conditional-law failure."""
    # This is the actual r=3 NC edge fixture: its row is positive but nearly
    # canceled in complex128.  The existing literal engine is the reference;
    # this check deliberately keeps every tiny positive helper coefficient.
    from experiments.experiment_nested_phase_edges import f3, f3_matrix, identity
    from experiments.experiment_earlier_phase_sampler_edges import literal_columns
    w03 = f3_matrix()
    old_middle = np.diag([1., 1j, -1j]) @ w03
    initial_e = np.diag([f3(j) for j in range(3)]) @ w03
    budget.reserve("nc_reference_setup", 2 * 3**3)
    budget.reserve("nc_reference_work", 8 * 3 * (3 * 3 + 6))
    reference = literal_columns(3, 3, 3, 2, initial_e, old_middle,
                                identity, 1, f3)[:, 0]
    exact = np.asarray(reference, dtype=np.complex128)
    spec = dict(period=3, block=3, width=3, split=2, initial=w03,
                middle=old_middle, phase=identity,
                early_phases=((0, f3), (1, f3)))
    sampler = make_sampler(cls, spec, cut=0, cache_right=False, budget=budget)
    limits = operation_limits(sampler.stats(), sampler, "row", 8)
    reserve_operation(budget, "nc_row", limits)
    before = {}
    row = sampler.row(0)
    check_operation_counters(budget, "nc_row", before, row["counters"], limits)
    expanded = expand(row, 8, budget)
    perturbed = np.asarray(expanded, dtype=np.complex128)
    exact_mass = float(np.sum(np.abs(exact) ** 2))
    perturbed_mass = float(np.sum(np.abs(perturbed) ** 2))
    budget.reserve("nc_fft_input_entries", 16)
    exact_law = np.abs(np.fft.fft(exact)) ** 2
    exact_law /= exact_law.sum()
    perturbed_law = np.abs(np.fft.fft(perturbed)) ** 2
    perturbed_law /= perturbed_law.sum()
    return {
        "amplitude_error": float(np.max(np.abs(exact - perturbed))),
        "exact_mass": exact_mass,
        "perturbed_mass": perturbed_mass,
        "conditional_tv": float(np.sum(np.abs(exact_law - perturbed_law)) / 2.),
        "work_weighted_joint_l1_half": float(np.sum(np.abs(
            exact_law*exact_mass/8-perturbed_law*perturbed_mass/8))/2),
    }


def main():
    exp = Experiment("nested_phase_sampler", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "streamed selected-cut rows reproduce amplitudes, norms, and full FFT laws")
    exp.predict("P2", "all required edge schedules and cache modes preserve their literal laws")
    exp.predict("P3", "mass and root_mass accepted q*acceptance laws reproduce conditionals")
    exp.must_fail("C1", "dropping a named early phase changes the independent literal law")
    exp.must_fail("C2", "near-canceled positive mass remains conditionally precision-unstable")

    report = {"status": "FAIL"}
    budget = Budget(MAX_WORK)
    budget.reserve("failed_report_reserve", 1)
    p1 = p2 = c1 = c2 = False
    p3 = True
    rows = {}
    try:
        report["preflight"] = preflight()
        specs = fixture_specs()
        from lab.work_first import NestedPhaseProgressions

        for name, spec in specs.items():
            q = 1 << spec["width"]
            budget.reserve("reference_setup", 2 * spec["block"] ** 3)
            reference = reference_columns(name, budget=budget, **spec)
            if reference.shape != (q, spec["period"]):
                raise AssertionError(f"{name} literal shape mismatch")
            column_norms = np.sum(np.abs(reference) ** 2, axis=1)
            if np.max(np.abs(column_norms - 1.)) > TOL:
                raise AssertionError(f"{name} literal columns are not normalized")

            auto_sampler = make_sampler(NestedPhaseProgressions, spec,
                                        cut="auto", cache_right=False,
                                        budget=budget)
            column_detail = verify_columns(auto_sampler, spec, reference, budget)
            cuts = (("auto", 0, spec["split"])
                    if name == "r60_multiple" else ("auto",))
            cache_modes = (False, True) if name == "r60_multiple" else (False,)
            for cache in cache_modes:
                for cut in cuts:
                    sampler = (auto_sampler if cut == "auto" and not cache else
                               make_sampler(NestedPhaseProgressions, spec,
                                            cut=cut, cache_right=cache,
                                            budget=budget))
                    detail, joint = compare_rows(name, spec, sampler, reference,
                                                 budget, accepted_modes=cut == "auto")
                    key = f"{name}/cut={cut}/cache_right={cache}"
                    detail["stats"] = safe(sampler.stats())
                    detail["cache_right"] = cache
                    detail["requested_cut"] = cut
                    detail["column_norm_error"] = float(np.max(np.abs(column_norms - 1.)))
                    detail.update(column_detail)
                    rows[key] = detail
                    if (detail["max_row_amplitude_error"] > TOL
                            or detail["max_row_norm_error"] > TOL
                            or detail["max_joint_probability_error"] > TOL
                            or abs(detail["work_marginal_total"] - 1.) > TOL):
                        raise AssertionError(f"{key} law mismatch")
                    if cut == "auto":
                        p3 = p3 and all(value <= TOL for value in
                                        (detail["accepted_q_acceptance_max_error"] or {}).values())

        # Deliberate science control: omit the final early phase in the named
        # multiple-phase schedule, while retaining the production row support.
        control_spec = dict(specs["r60_multiple"])
        control_spec["early_phases"] = control_spec["early_phases"][:-1]
        control_spec["positions"] = control_spec["positions"][:-1]
        control_sampler = make_sampler(NestedPhaseProgressions, control_spec,
                                       cut="auto", cache_right=False, budget=budget)
        control_reference = reference_columns("r60_multiple", budget=budget,
                                              **specs["r60_multiple"])
        wrong_reference = reference_columns("r60_multiple", budget=budget,
                                            **control_spec)
        control_error = float(np.max(np.abs(control_reference - wrong_reference)))
        control_stats = control_sampler.stats()
        control_limits = operation_limits(control_stats, control_sampler, "row",
                                          1 << control_spec["width"])
        for work in range(control_spec["period"]):
            reserve_operation(budget, "control_row", control_limits)
            before = {}
            wrong_row = control_sampler.row(work)
            check_operation_counters(budget, "control_row", before,
                                     wrong_row["counters"], control_limits)
            wrong_values = expand(wrong_row, 1 << control_spec["width"], budget)
            control_error = max(control_error,
                                float(np.max(np.abs(wrong_values
                                                    - control_reference[:, work]))))
        budget.reserve("control_fft_input_entries", 2 * (1 << control_spec["width"])
                       * control_spec["period"])
        control_law_error = float(np.max(np.abs(
            np.abs(np.fft.fft(control_reference, axis=0)) ** 2
            - np.abs(np.fft.fft(wrong_reference, axis=0)) ** 2)))
        c1 = control_error > TOL and control_law_error > TOL
        report["wrong_phase_control"] = {"amplitude_error": control_error,
                                          "law_unnormalized_error": control_law_error}

        nc = near_cancellation_control(NestedPhaseProgressions, budget)
        c2 = (nc["amplitude_error"] < TOL and nc["exact_mass"] > 0.
              and nc["perturbed_mass"] > 0. and nc["conditional_tv"] > 1e-3)
        report["near_cancellation_control"] = nc
        p1 = all(value["max_row_amplitude_error"] <= TOL
                 and value["max_row_norm_error"] <= TOL
                 and value["max_joint_probability_error"] <= TOL
                 for value in rows.values())
        revival = rows["r6_revival/cut=auto/cache_right=False"]
        p2 = (p1 and bool(rows) and any(value["zero_rows"] for value in rows.values())
              and "r3_alias/cut=auto/cache_right=False" in rows
              and "r9_nondivisor/cut=auto/cache_right=False" in rows
              and revival["revival_work1_exponent0"] is not None
              and abs(complex(*revival["revival_work1_exponent0"])
                      - (0.5 + 0.5j)) <= TOL)
        cached_key = "r60_multiple/cut=auto/cache_right=True"
        plain_key = "r60_multiple/cut=auto/cache_right=False"
        cache_equal = all(abs(float(rows[cached_key][field])
                              - float(rows[plain_key][field])) <= TOL
                          for field in ("max_row_amplitude_error",
                                        "max_row_norm_error",
                                        "max_joint_probability_error",
                                        "work_marginal_total"))
        p2 = p2 and cache_equal
        report["edge_witnesses"] = {
            "r6_revival": revival["revival_work1_exponent0"],
            "r3_alias_present": "r3_alias/cut=auto/cache_right=False" in rows,
            "r9_nondivisor_present": "r9_nondivisor/cut=auto/cache_right=False" in rows,
            "cached_auto_law_equal": cache_equal,
        }
        report["budget_used"] = dict(budget.used)
        report["budget_total_used"] = sum(budget.used.values())
        report["counter_reconciliation"] = reconcile_budget(budget)
        report["rows"] = rows
        report["status"] = "PASS" if all((p1, p2, p3, c1, c2)) else "FAIL"
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()
        report["budget_used"] = dict(budget.used)
        report["budget_total_used"] = sum(budget.used.values())
        report["rows"] = rows
    exp.check("P1", p1, "all selected cuts matched streamed literal rows and FFT laws")
    exp.check("P2", p2, "edge schedules, alias, zero rows, and cache mode were exercised")
    exp.check("P3", p3, "both proposal modes matched accepted q*acceptance law")
    exp.fail_check("C1", c1, "omitting the final multiple-phase insertion changes amplitude/law")
    exp.fail_check("C2", c2, "positive near-canceled mass exposes conditional precision instability")
    report["observed_counters"] = dict(budget.observed)
    report["observed_counter_checksum"] = sum(budget.observed.values())
    report["budget_semantics"] = "budget_used/total are reservations; observed counters are separate"
    return exp.finish(report_path=report_path(), rows=report.get("rows", {}),
                      metadata={"result": report, "no_timing_claim": True,
                                "rng_distribution_audit_deferred": True})


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
