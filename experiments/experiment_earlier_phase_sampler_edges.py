"""Bounded edge and RNG-path checks for the opt-in earlier-phase sampler.

This test is deliberately small: an independent indexed-column FFT oracle is
used for r=6,b=2,L=4 cancellation revival, while a separate r=3,b=1,width=3
fixture exercises stride six, gcd lifting, retries, and exhaustion through the
actual sampler RNG path.  It is not a generic propagator or a sampler claim.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1 cycle/dual/auto rows and forced joints agree with literal columns,
     including the r=6,b=2 cancellation-revival row.
  P2 aliases (r=b), arbitrary stride, zero rows, noncallable/cap inputs, and
     callable-but-invalid phase queries take the documented paths.
  P3 scripted actual sampling preserves work across a rejection/retry and
     raises on a finite cap; the r=3 stride-six path records two fair bits and
     a gcd lift rather than silently using a unit-stride shortcut.
  C1 an intentionally dropped early phase or wrong stride is not accepted by
     the independent literal joint law checks.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment


MAX_BYTES = 16 << 20
MAX_WORK = 250_000
TOL = 2e-10
REFERENCE_LEDGER = None
REFERENCE_BOUNDS = None
SETUP_LEDGER = None


def reference_charge(key, amount=1):
    if REFERENCE_LEDGER is not None:
        if REFERENCE_LEDGER[key]+amount > REFERENCE_BOUNDS[key]:
            raise MemoryError(f"reference {key} exceeds budget before operation")
        REFERENCE_LEDGER[key] += amount


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"earlier_phase_sampler_edges_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value} exceeds 16 MiB")
    return value


def h2():
    return np.asarray([[1., 1.], [1., -1.]], dtype=np.complex128) / math.sqrt(2.)


def literal_columns(period, block, width, split, initial, middle,
                    phase, early_split, early_phase):
    """Independent indexed schedule, retaining no physical orbit labels."""
    Q, L = 1 << width, 1 << split
    columns = np.zeros((Q, period), dtype=np.complex128)
    ledger = REFERENCE_LEDGER
    for exponent in range(Q):
        reference_charge("column_calls")
        low, high = exponent % L, exponent // L
        early = low % (1 << early_split)
        vector = np.zeros(period, dtype=np.complex128)
        for u in range(block):
            reference_charge("early_phase_queries")
            before = (u + early) % period
            after = (before + low - early) % period
            vector[after] += initial[u, 0] * early_phase(before)
        mixed = np.zeros(period, dtype=np.complex128)
        for cell in range(period // block):
            start = cell * block
            reference_charge("column_local_terms", block*block)
            mixed[start:start + block] = middle @ vector[start:start + block]
        for index in np.flatnonzero(mixed):
            reference_charge("late_phase_queries")
            mixed[index] *= phase(int(index))
        reference_charge("shift_operations")
        columns[exponent] = np.roll(mixed, L * high)
    return columns


def literal_joint(columns):
    Q = columns.shape[0]
    reference_charge("fft_calls")
    reference_charge("fft_input_entries", int(columns.size))
    amplitudes = np.fft.fft(columns, axis=0) / Q
    return np.abs(amplitudes) ** 2


def valid_joint(values):
    values = np.asarray(values, dtype=float)
    return bool(values.ndim == 2 and np.all(np.isfinite(values))
                and np.min(values) >= -TOL
                and abs(float(values.sum()) - 1.) < TOL)


def tv(left, right):
    return float(np.abs(np.asarray(left) - np.asarray(right)).sum() / 2.)


class ScriptedRNG:
    """Small deterministic RNG exposing every actual primitive draw."""

    def __init__(self, random_values, integer_values):
        self.random_values = list(random_values)
        self.integer_values = list(integer_values)
        self.random_calls = 0
        self.integer_calls = 0

    def random(self):
        if not self.random_values:
            raise AssertionError("scripted random stream exhausted")
        self.random_calls += 1
        value = float(self.random_values.pop(0))
        if not 0 <= value < 1:
            raise AssertionError("scripted real outside RNG contract [0,1)")
        return value

    def integers(self, high):
        if not self.integer_values:
            raise AssertionError("scripted integer stream exhausted")
        value = int(self.integer_values.pop(0))
        if not 0 <= value < int(high):
            raise AssertionError(f"scripted integer {value} outside [0,{high})")
        self.integer_calls += 1
        return value


def preflight():
    # Literal matrices for r=6 and r=3, joint laws, and a bounded collection
    # of actual path records.  No full Q-by-r family is allocated.
    numeric = (8 * 6 * 16 + 8 * 3 * 16 + 2 * 8 * 6 * 8
               + 2 * 8 * 3 * 8 + 256 * 1024)
    path_count_bound = 8 * (3 + 3 + 2) * 4
    attempt_call_bound = 2 * path_count_bound
    retained_path_bytes = path_count_bound * 64
    reference_bounds = {
        "column_calls": 38,
        "column_local_terms": 152,
        "early_phase_queries": 44,
        "late_phase_queries": 124,
        "shift_operations": 38,
        "fft_calls": 5,
        "fft_input_entries": 124,
        "control_root_terms": 2*16*16,
    }
    guard_bytes(numeric + retained_path_bytes,
                "earlier-phase sampler edge aggregate payload")
    # <=128 row/forced/invalid calls at <=50 named units, plus <=512
    # attempts at <=40 units and <=7 RNG primitives each. These conservative
    # executable budgets precede measurement; the category cap is unchanged.
    operation_bound = (128*50 + attempt_call_bound*(40+7)
                       + sum(reference_bounds.values()) + 16*2*2**3)
    if operation_bound > MAX_WORK:
        raise MemoryError("earlier-phase sampler edge work cap")
    return {"payload_bytes": numeric + retained_path_bytes,
            "payload_components": {"base_numeric": numeric,
                                    "retained_path_records": retained_path_bytes},
            "path_count_bound": path_count_bound,
            "attempt_call_bound": attempt_call_bound,
            "reference_operation_bounds": reference_bounds,
            "operation_bound": operation_bound,
            "max_payload_bytes": MAX_BYTES, "max_work": MAX_WORK}


def import_sampler():
    # Keep the import inside the released run: the production API is opt-in
    # and was intentionally not guessed while it was under construction.
    from lab.work_first import EarlierPhaseProgressions
    class AuditedEarlierPhaseProgressions(EarlierPhaseProgressions):
        def __init__(self, *args, **kwargs):
            self.audit_counter_dicts = []
            super().__init__(*args, **kwargs)

        def _new_counters(self):
            counters = super()._new_counters()
            self.audit_counter_dicts.append(counters)
            return counters
    return AuditedEarlierPhaseProgressions


NONOVERLAP_COUNTER_KEYS = (
    "column_local_terms", "row_local_terms", "phase_queries",
    "fourier_component_terms", "progression_proposals", "marginal_queries",
    "work_draws", "component_draws", "acceptance_draws",
)


def counter_units(counters):
    """Count total work without double-counting early/late phase subcounts."""
    return sum(int(counters.get(key, 0)) for key in NONOVERLAP_COUNTER_KEYS)


def audit_units(sampler):
    return sum(counter_units(counters)
               for counters in getattr(sampler, "audit_counter_dicts", ()))


def reserve_call(sampler, ledger, kind, *, max_attempts=1):
    stats = sampler.stats()
    pair_bound = max(stats.get("cycle_row_pair_visit_bound", 0),
                     stats.get("dual_row_pair_visit_bound", 0))
    phase_bound = stats.get("phase_queries_per_draw_bound", 0)
    component_bound = stats.get("component_bound", 0)
    local_bound = stats.get("local_term_bound", 0)
    if kind == "row":
        bound = pair_bound + phase_bound + component_bound
    elif kind == "forced":
        bound = pair_bound + phase_bound + component_bound * 2
    elif kind == "sample":
        bound = (local_bound + pair_bound + phase_bound
                 + max_attempts * (component_bound + 2*sampler.width + 5))
    elif kind == "column":
        bound = sampler.b**2 + phase_bound
    else:
        raise ValueError(f"unknown audited call kind {kind}")
    if ledger["reserved_units"] + bound > MAX_WORK:
        raise MemoryError("audited per-call reservation exceeds frozen work cap")
    ledger["reserved_units"] += int(bound)
    ledger["call_counts"][kind] = ledger["call_counts"].get(kind, 0) + 1
    return bound


def sampler_kwargs(cls, *, period, block, width, split, initial, middle,
                   phase, early_split, early_phase, cover, **extra):
    options = dict(max_local_terms=250_000,
                   max_components=4096, max_payload_bytes=MAX_BYTES)
    options.update(extra)
    if SETUP_LEDGER["constructor_calls"] >= 16:
        raise MemoryError("constructor count cap before call")
    SETUP_LEDGER["constructor_calls"] += 1
    result = cls(period, block, width, split, initial, middle, phase,
                 early_split=early_split, early_phase=early_phase,
                 cover=cover, **options)
    # Successful constructors run the existing two dense b-by-b checks.
    SETUP_LEDGER["unitary_check_product_terms"] += 2*block**3
    return result


def row_components(row):
    components = row.get("components", ())
    parsed = []
    for item in components:
        if len(item) != 3:
            raise AssertionError("row component must remain a three-tuple")
        parsed.append(item)
    return parsed


def main():
    global REFERENCE_LEDGER, REFERENCE_BOUNDS, SETUP_LEDGER
    exp = Experiment("earlier_phase_sampler_edges", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "literal joint law agrees for cycle, dual, and auto covers")
    exp.predict("P2", "aliases, zero rows, and invalid/cap paths are guarded")
    exp.predict("P3", "actual RNG paths preserve work and exercise gcd lift/retry")
    exp.must_fail("C1", "dropping early phase or stride does not pass joint checks")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        SETUP_LEDGER = {"constructor_calls": 0,"unitary_check_product_terms": 0}
        EarlierPhaseProgressions = import_sampler()
        H2 = h2()
        audit_ledger = {"reserved_units": 0, "call_counts": {}}
        audited_samplers = []
        reference_ledger = {"column_calls": 0, "column_local_terms": 0,
                            "early_phase_queries": 0, "late_phase_queries": 0,
                            "shift_operations": 0, "fft_calls": 0,
                            "fft_input_entries": 0, "control_root_terms": 0}
        REFERENCE_LEDGER = reference_ledger
        REFERENCE_BOUNDS = report["preflight"]["reference_operation_bounds"]
        parity = lambda j: complex(1. if int(j) % 2 == 0 else -1.)
        # Cancellation-revival fixture: old rho=0 cancels, early phase revives it.
        columns = literal_columns(6, 2, 2, 2, H2, H2, parity, 1, parity)
        joint = literal_joint(columns)
        if not valid_joint(joint):
            raise AssertionError("literal r6 joint law invalid")
        cover_rows = {}
        forced_errors = {}
        row_component_counts = {}
        forced_work_sums = {}
        forced_counter_totals = {}
        for cover in ("cycle", "dual", "auto"):
            sampler = sampler_kwargs(EarlierPhaseProgressions, period=6, block=2,
                                     width=2, split=2, initial=H2, middle=H2,
                                     phase=parity, early_split=1,
                                     early_phase=parity, cover=cover)
            audited_samplers.append(sampler)
            row_component_counts[cover] = {}
            forced_counter_totals[cover] = {}
            max_error = 0.
            work_sums = np.zeros(6, dtype=float)
            for work in range(6):
                reserve_call(sampler, audit_ledger, "row")
                row = sampler.row(work)
                comps = row_components(row)
                row_component_counts[cover][str(work)] = len(comps)
                for output in range(4):
                    reserve_call(sampler, audit_ledger, "forced")
                    forced = sampler.forced_joint(work, output)
                    got = forced["joint_probability"]
                    for key, value in forced["counters"].items():
                        if key in NONOVERLAP_COUNTER_KEYS:
                            forced_counter_totals[cover][key] = (
                                forced_counter_totals[cover].get(key, 0) + value)
                    work_sums[work] += float(got)
                    max_error = max(max_error, abs(float(got) - joint[output, work]))
            cover_rows[cover] = sampler.stats()
            forced_errors[cover] = max_error
            forced_work_sums[cover] = work_sums.tolist()
        # At j=1 the canceled old rho=0 class must be represented after the
        # early phase; this catches an implementation that only refines old
        # nonzero coefficients.
        revival_counts = {cover: row_component_counts[cover]["1"]
                          for cover in row_component_counts}

        # Aliasing endpoint r=b, and an independent zero-row search.
        identity2 = np.eye(2, dtype=np.complex128)
        alias_sampler = sampler_kwargs(EarlierPhaseProgressions, period=2, block=2,
                                       width=1, split=1, initial=identity2,
                                       middle=identity2, phase=lambda j: 1.+0j,
                                       early_split=0, early_phase=lambda j: 1.+0j,
                                       cover="auto")
        audited_samplers.append(alias_sampler)
        alias_rows = []
        for j in range(2):
            reserve_call(alias_sampler, audit_ledger, "row")
            alias_rows.append(alias_sampler.row(j))
        alias_columns = literal_columns(2, 2, 1, 1, identity2, identity2,
                                        lambda j: 1.+0j, 0,
                                        lambda j: 1.+0j)
        alias_joint = literal_joint(alias_columns)
        alias_error = 0.
        for work in range(2):
            for output in range(2):
                reserve_call(alias_sampler, audit_ledger, "forced")
                forced = alias_sampler.forced_joint(work, output)
                alias_error = max(alias_error,
                                  abs(float(forced["joint_probability"])
                                      - alias_joint[output, work]))
        zero_sampler = sampler_kwargs(EarlierPhaseProgressions, period=6, block=2,
                                      width=2, split=2, initial=identity2,
                                      middle=identity2, phase=lambda j: 1.+0j,
                                      early_split=1,
                                      early_phase=lambda j: 1.+0j, cover="auto")
        audited_samplers.append(zero_sampler)
        zero_rows = []
        for work in range(6):
            reserve_call(zero_sampler, audit_ledger, "row")
            row = zero_sampler.row(work)
            if row.get("norm", 0.) == 0:
                zero_rows.append(work)
                reserve_call(zero_sampler, audit_ledger, "forced")
                forced = zero_sampler.forced_joint(work, 0)
                if forced["work_probability"] != 0 or forced["joint_probability"] != 0:
                    raise AssertionError("zero row returned nonzero forced mass")

        # Invalid booleans, invalid phase, cover, and structural caps must
        # reject before invoking a phase oracle or copying matrices.
        phase_calls = {"n": 0}
        def counted_phase(j):
            phase_calls["n"] += 1
            return 1.+0j
        invalid = {}
        for name, kwargs in (
            ("bool_early_split", {"early_split": True}),
            ("bad_cover", {"cover": "wrong"}),
            ("bad_phase", {"phase": lambda j: 2.+0j}),
            ("bad_early_phase", {"early_phase": lambda j: 2.+0j}),
            ("noncallable_phase", {"phase": None}),
            ("noncallable_early_phase", {"early_phase": None}),
        ):
            constructor_rejected = False
            query_rejected = False
            try:
                base = dict(period=6, block=2, width=2, split=2,
                            initial=H2, middle=H2, phase=counted_phase,
                            early_split=1, early_phase=counted_phase, cover="auto")
                base.update(kwargs)
                candidate = sampler_kwargs(EarlierPhaseProgressions, **base)
                audited_samplers.append(candidate)
                if name in ("bad_phase", "bad_early_phase"):
                    try:
                        reserve_call(candidate,audit_ledger,"column")
                        candidate.column(0)
                    except ValueError:
                        query_rejected = True
            except (ValueError, MemoryError):
                constructor_rejected = True
            invalid[name] = constructor_rejected
            if name in ("bad_phase", "bad_early_phase"):
                invalid[f"{name}_query_rejected"] = query_rejected
                invalid[f"{name}_before_copy"] = constructor_rejected
        old_array = np.array
        copy_trap = {"called": False}
        def trap_array(*args, **kwargs):
            copy_trap["called"] = True
            raise AssertionError("matrix copy occurred before cap rejection")
        np.array = trap_array
        try:
            for name, key in (("noncallable_phase_before_copy","phase"),
                              ("noncallable_early_before_copy","early_phase")):
                kwargs = dict(period=6,block=2,width=2,split=2,initial=H2,
                              middle=H2,phase=counted_phase,early_split=1,
                              early_phase=counted_phase,cover="auto")
                kwargs[key] = None
                try:
                    sampler_kwargs(EarlierPhaseProgressions,**kwargs)
                except ValueError:
                    invalid[name] = not copy_trap["called"]
                else:
                    invalid[name] = False
            try:
                sampler_kwargs(EarlierPhaseProgressions, period=6, block=2,
                               width=2, split=2, initial=H2, middle=H2,
                               phase=counted_phase, early_split=1,
                               early_phase=counted_phase, cover="auto",
                               max_components=1)
            except MemoryError:
                invalid["component_cap_before_copy"] = not copy_trap["called"]
            else:
                invalid["component_cap_before_copy"] = False
        finally:
            np.array = old_array

        # Actual RNG transition fixture: r=3,b=1, width=3, split=2,
        # early_split=1 gives stride r*P=6, reduced dimension four, lift two.
        one = np.ones((1, 1), dtype=np.complex128)
        rng_sampler = sampler_kwargs(EarlierPhaseProgressions, period=3, block=1,
                                     width=3, split=2, initial=one, middle=one,
                                     phase=lambda j: 1.+0j, early_split=1,
                                     early_phase=parity, cover="cycle")
        audited_samplers.append(rng_sampler)
        rng_columns = literal_columns(3, 1, 3, 2, one, one,
                                      lambda j: 1.+0j, 1, parity)
        rng_joint = literal_joint(rng_columns)
        rng_rows = []
        for work in range(3):
            reserve_call(rng_sampler, audit_ledger, "row")
            rng_rows.append(rng_sampler.row(work))
        m_by_work = [len(row_components(row)) for row in rng_rows]
        if m_by_work != [3, 3, 2]:
            raise AssertionError(f"unexpected r3 component multiplicities {m_by_work}")
        path_count = sum(m_by_work[seed % 3] for seed in range(8)) * 4 * 2
        if path_count > MAX_WORK:
            raise MemoryError("full RNG path enumeration exceeds frozen work cap")
        accepted_mass = np.zeros((8, 3), dtype=float)
        rejected_mass = np.zeros((8, 3), dtype=float)
        rng_paths = []
        rejection_calls = 0
        failed_sample_calls = 0
        path_rng_shapes = set()
        actual_rng_random_calls = 0
        actual_rng_integer_calls = 0
        for seed in range(8):
            work = seed % 3
            m = m_by_work[work]
            for component in range(m):
                component_draw = (component + .5) / m
                for k in range(4):
                    # interval_path appends bit<<step, so the first scripted
                    # draw is the least-significant bit.
                    bits = (k & 1, (k >> 1) & 1)
                    for lift in range(2):
                        output = (3 * k) % 4 + 4 * lift
                        acceptance = 64. * float(rng_joint[output, work]) / (m*m)
                        if not -TOL <= acceptance <= 1.+TOL:
                            raise AssertionError("invalid analytic acceptance")
                        acceptance = min(1., max(0., acceptance))
                        path_weight = 1. / (64. * m)
                        accepted_mass[output, work] += path_weight * acceptance
                        rejected_mass[output, work] += path_weight * (1.-acceptance)
                        bit_values = [.25 if bit == 0 else .75 for bit in bits]
                        rng = ScriptedRNG(
                            [.5, component_draw, *bit_values,
                             acceptance/2. if acceptance > 0 else .5],
                            [seed, lift])
                        got = None
                        try:
                            reserve_call(rng_sampler, audit_ledger, "sample")
                            got = rng_sampler.sample(rng, max_attempts=1)
                        except RuntimeError:
                            failed_sample_calls += 1
                            if acceptance > 0:
                                raise
                        path_rng_shapes.add((rng.random_calls, rng.integer_calls))
                        actual_rng_random_calls += rng.random_calls
                        actual_rng_integer_calls += rng.integer_calls
                        if acceptance > 0:
                            if got is None or got["work"] != work or got["output"] != output:
                                raise AssertionError("accepted scripted path changed work/output")
                        elif got is not None:
                            raise AssertionError("zero-acceptance path returned a sample")
                        rng_paths.append({"seed": seed, "work": work,
                                          "component": component, "k": k,
                                          "lift": lift, "output": output,
                                          "acceptance": acceptance,
                                          "accepted": acceptance > 0,
                                          "random_calls": rng.random_calls,
                                          "integer_calls": rng.integer_calls})
                        if acceptance < 1.-TOL:
                            rejection_calls += 1
                            reject_rng = ScriptedRNG(
                                [.5, component_draw, *bit_values,
                                 (1.+acceptance)/2.], [seed, lift])
                            try:
                                reserve_call(rng_sampler, audit_ledger, "sample")
                                rng_sampler.sample(reject_rng, max_attempts=1)
                            except RuntimeError:
                                failed_sample_calls += 1
                                path_rng_shapes.add((reject_rng.random_calls,
                                                     reject_rng.integer_calls))
                                actual_rng_random_calls += reject_rng.random_calls
                                actual_rng_integer_calls += reject_rng.integer_calls
                                pass
                            else:
                                raise AssertionError("positive reject probability was not rejected")

        work_probs = np.asarray(m_by_work, dtype=float) / 8.
        accepted_conditionals = np.zeros_like(accepted_mass)
        restored_joint = np.zeros_like(accepted_mass)
        for work, m in enumerate(m_by_work):
            accepted_conditionals[:, work] = accepted_mass[:, work] / (work_probs[work]/m)
            restored_joint[:, work] = accepted_conditionals[:, work] * work_probs[work]
        restored_error = float(np.max(np.abs(restored_joint-rng_joint)))
        restored_norm = float(restored_joint.sum())
        accepted_total = float(accepted_mass.sum())
        rejected_total = float(rejected_mass.sum())

        # Exact partial-acceptance retry: e=0,y=0 has A=1/9 in this fixture.
        retry_acceptance = 64.*float(rng_joint[0, 0])/(m_by_work[0]**2)
        if abs(retry_acceptance - 1./9.) > TOL:
            raise AssertionError(f"frozen retry acceptance changed: {retry_acceptance}")
        retry_rng = ScriptedRNG(
            [.5, 1./6., .25, .25, 5./9., 1./6., .25, .25, 1./18.],
            [0, 0, 0])
        reserve_call(rng_sampler, audit_ledger, "sample", max_attempts=2)
        retry = rng_sampler.sample(retry_rng, max_attempts=2)
        actual_rng_random_calls += retry_rng.random_calls
        actual_rng_integer_calls += retry_rng.integer_calls
        if retry["attempts"] != 2 or retry["work"] != 0 or retry["output"] != 0:
            raise AssertionError("retry path did not preserve the frozen work/output")
        exhaust_rng = ScriptedRNG([.5, 1./6., .25, .25, 5./9.], [0, 0])
        try:
            reserve_call(rng_sampler, audit_ledger, "sample")
            rng_sampler.sample(exhaust_rng, max_attempts=1)
        except RuntimeError:
            exhausted = True
            failed_sample_calls += 1
        else:
            exhausted = False
        actual_rng_random_calls += exhaust_rng.random_calls
        actual_rng_integer_calls += exhaust_rng.integer_calls
        actual_attempt_calls = path_count + rejection_calls + 2 + 1
        if actual_attempt_calls > 512:
            raise MemoryError("actual scripted attempt calls exceed preflight")

        # Independent restored joint law, plus deliberate wrong controls.
        rng_sampler_joint = np.zeros_like(rng_joint)
        rng_forced_counters = {}
        for work in range(3):
            for output in range(8):
                reserve_call(rng_sampler, audit_ledger, "forced")
                forced = rng_sampler.forced_joint(work, output)
                rng_sampler_joint[output, work] = forced["joint_probability"]
                for key, value in forced["counters"].items():
                    if key in NONOVERLAP_COUNTER_KEYS:
                        rng_forced_counters[key] = rng_forced_counters.get(key, 0) + value
        rng_forced_error = float(np.max(np.abs(rng_sampler_joint-rng_joint)))
        rng_work_sums = rng_sampler_joint.sum(axis=0)
        dropped = literal_joint(literal_columns(3, 1, 3, 2, one, one,
                                                lambda j: 1.+0j, 1,
                                                lambda j: 1.+0j))
        # Genuine wrong-stride control: use a width-4 row with count>1 and
        # expand its components with stride r instead of the declared 6.
        wide_columns = literal_columns(3, 1, 4, 3, one, one,
                                       lambda j: 1.+0j, 1, parity)
        wide_joint = literal_joint(wide_columns)
        wide_sampler = sampler_kwargs(EarlierPhaseProgressions, period=3, block=1,
                                      width=4, split=3, initial=one, middle=one,
                                      phase=lambda j: 1.+0j, early_split=1,
                                      early_phase=parity, cover="cycle")
        audited_samplers.append(wide_sampler)
        wrong_stride_joint = np.zeros_like(wide_joint)
        correct_stride_joint = np.zeros_like(wide_joint)
        wrong_stride_rows = 0
        for work in range(3):
            reserve_call(wide_sampler, audit_ledger, "row")
            row = wide_sampler.row(work)
            comps = row_components(row)
            if any(count > 1 for _start, count, _amp in comps):
                wrong_stride_rows += 1
            norm = row["norm"]
            for output in range(16):
                coherent = 0j
                correct_coherent = 0j
                for start, count, amplitude in comps:
                    reference_charge("control_root_terms",2*count)
                    geom = sum(np.exp(-2j*np.pi*output*(start+3*n)/16)
                               for n in range(count))
                    coherent += amplitude*geom
                    correct_geom = sum(np.exp(-2j*np.pi*output
                                               *(start+row["stride"]*n)/16)
                                       for n in range(count))
                    correct_coherent += amplitude*correct_geom
                wrong_stride_joint[output, work] = abs(coherent)**2/(16*16)
                correct_stride_joint[output, work] = abs(correct_coherent)**2/(16*16)
        wrong_stride_half_l1 = tv(wide_joint.ravel(), wrong_stride_joint.ravel())
        correct_stride = tv(wide_joint.ravel(), correct_stride_joint.ravel())
        wrong_stride_norm = float(wrong_stride_joint.sum())
        correct_stride_norm = float(correct_stride_joint.sum())
        # Wrong stride changes the input norm through spurious collisions.
        # Report that defect; TV requires normalizing the altered law first.
        wrong_stride = tv(wide_joint, wrong_stride_joint/wrong_stride_norm)
        audited_counter_total = sum(audit_units(sampler) for sampler in audited_samplers)
        actual_counter_total = (audited_counter_total + actual_rng_random_calls
                                + actual_rng_integer_calls + sum(reference_ledger.values())
                                + SETUP_LEDGER["unitary_check_product_terms"])
        actual_work_within_bound = actual_counter_total <= MAX_WORK
        actual_operation_bound_ok = actual_counter_total <= report["preflight"]["operation_bound"]
        reference_bounds_ok = all(reference_ledger[key] <= value
            for key,value in REFERENCE_BOUNDS.items())
        reservation_within_bound = audited_counter_total <= audit_ledger["reserved_units"] <= MAX_WORK
        c1 = (tv(rng_joint.ravel(), dropped.ravel()) > 1e-3
              and wrong_stride > 1e-3 and correct_stride < TOL
              and abs(correct_stride_norm-1.) < TOL)
        p1 = (reference_bounds_ok
              and all(error < TOL for error in forced_errors.values())
              and all(count >= 2 for count in revival_counts.values())
              and valid_joint(joint))
        p2 = (alias_error < TOL
              and invalid["bad_phase_query_rejected"]
              and invalid["bad_early_phase_query_rejected"]
              and all(invalid[name] for name in
                  ("bool_early_split", "bad_cover", "noncallable_phase",
                   "noncallable_early_phase", "component_cap_before_copy"))
              and invalid["noncallable_phase_before_copy"]
              and invalid["noncallable_early_before_copy"]
              and bool(zero_rows)
              and valid_joint(alias_joint)
              and all(abs(row["norm"]-1.) < TOL for row in alias_rows))
        p3 = (len(rng_paths) == path_count and retry["attempts"] == 2
              and retry["counters"]["work_draws"] == 1 and exhausted
              and rejection_calls > 0
              and restored_error < TOL and abs(restored_norm - 1.) < TOL
              and abs(accepted_total + rejected_total - 1.) < TOL
              and path_rng_shapes == {(5, 2)}
              and rng_forced_error < TOL
              and abs(float(rng_sampler_joint.sum()) - 1.) < TOL
              and actual_work_within_bound and actual_operation_bound_ok
              and reservation_within_bound)
        report.update({
            "fixture": {"r6_b2": {"period": 6, "block": 2, "width": 2,
                                    "split": 2, "early_split": 1},
                        "rng": {"period": 3, "block": 1, "width": 3,
                                "split": 2, "early_split": 1}},
            "forced_joint_max_errors": forced_errors,
            "forced_work_sums": forced_work_sums,
            "forced_counter_totals": forced_counter_totals,
            "row_component_counts": row_component_counts,
            "cancellation_revival_counts": revival_counts,
            "cover_stats": cover_rows,
            "zero_rows": zero_rows,
            "invalid_and_cap_checks": invalid,
            "alias_max_error": alias_error,
            "rng_path_count": len(rng_paths),
            "rng_path_expected_count": path_count,
            "rng_rejection_calls": rejection_calls,
            "rng_path_draw_shapes": sorted(path_rng_shapes),
            "accepted_mass": accepted_mass.tolist(),
            "rejected_mass": rejected_mass.tolist(),
            "accepted_conditional_laws": accepted_conditionals.tolist(),
            "restored_joint": restored_joint.tolist(),
            "restored_joint_max_error": restored_error,
            "restored_joint_norm": restored_norm,
            "rng_paths_sample": rng_paths[:16],
            "retry": {"work": retry["work"], "attempts": retry["attempts"],
                      "counters": retry["counters"],
                      "random_calls": retry_rng.random_calls,
                      "integer_calls": retry_rng.integer_calls},
            "exhaustion_raised": exhausted,
            "rng_forced_max_error": rng_forced_error,
            "rng_forced_work_sums": rng_work_sums.tolist(),
            "rng_forced_joint_norm": float(rng_sampler_joint.sum()),
            "rng_forced_counters": rng_forced_counters,
            "actual_counter_total": actual_counter_total,
            "audited_counter_total": audited_counter_total,
            "actual_attempt_calls": actual_attempt_calls,
            "failed_sample_calls": failed_sample_calls,
            "actual_rng_random_calls": actual_rng_random_calls,
            "actual_rng_integer_calls": actual_rng_integer_calls,
            "accepted_total": accepted_total,
            "rejected_total": rejected_total,
            "accepted_plus_rejected": accepted_total + rejected_total,
            "actual_work_within_bound": actual_work_within_bound,
            "actual_operation_bound_ok": actual_operation_bound_ok,
            "reservation_within_bound": reservation_within_bound,
            "audit_ledger": audit_ledger,
            "reference_ledger": reference_ledger,
            "setup_ledger": dict(SETUP_LEDGER),
            "reference_bounds_ok": reference_bounds_ok,
            "wrong_stride_tv": wrong_stride,
            "wrong_stride_unnormalized_half_l1": wrong_stride_half_l1,
            "correct_stride_tv": correct_stride,
            "wrong_stride_norm": wrong_stride_norm,
            "correct_stride_norm": correct_stride_norm,
            "wrong_stride_rows_with_count_gt1": wrong_stride_rows,
            "dropped_early_phase_tv": tv(rng_joint.ravel(), dropped.ravel()),
            "checks": {"P1": p1, "P2": p2, "P3": p3, "C1": c1},
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "literal and all cover joint laws agree")
        exp.check("P2", p2, "edge/cap/zero-row paths are guarded")
        exp.check("P3", p3, "actual RNG path preserves work and exhausts honestly")
        exp.fail_check("C1", c1, "wrong early phase/stride controls differ")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": "bounded r6/b2 cancellation and r3/b1 stride-six RNG path",
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_work": MAX_WORK,
        "reference": "independent indexed literal columns and actual EarlierPhaseProgressions RNG path",
        "no_generic_propagator": True,
        "no_production_claim": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
