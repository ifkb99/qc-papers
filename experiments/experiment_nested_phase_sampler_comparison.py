"""Bounded returned-sample comparison for the nested-phase work sampler.

This freezes the first matched returned-sample comparison requested by TODO41.
The physical state is the NC/C81 fixture with ``r=60,b=3,t=13,s=12`` and
early phases at controls 5 and 11.  Six construction configurations are
compared in both proposal modes: the three fixed cuts, auto without a cache,
auto with right-phase caching, and auto with both right caching and one
explicit 60-entry phase table.  Each configuration/mode receives four sample
requests with ``max_attempts=1024``.  Failed calls are retained with their
``last_counters`` and are never redrawn.

The helper's four forced output marginals are compared with the existing
``lab.semiclassical.sequential_path`` instrument using supplied dense 60-by-60
branches.  The sequential reference is only a small matched output check; its
dense setup, matrix products, QR shapes and cubic work are charged separately
from scalar progression/Fourier work.  No Q-by-r archive, timing result or
finite-bit sampling certificate is produced.

Run only after review and confirmation of this frozen design:
``uv run python -m experiments.experiment_nested_phase_sampler_comparison``.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import time
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.fourier_sampling import unit_phase
from lab.semiclassical import sequential_path
from lab.work_first import NestedPhaseProgressions
from experiments.experiment_clean_orbit_output import work_block


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT = 13, 12
Q, L, H = 1 << WIDTH, 1 << SPLIT, 1 << (WIDTH - SPLIT)
POSITIONS = (5, 11)
OUTPUTS = (0, 1, Q // 2, Q - 1)
PROPOSALS = ("mass", "root_mass")
MAX_ATTEMPTS = 1024
SAMPLE_CALLS = 4
TOL = 3e-10
MAX_BYTES = 16 << 20
MAX_CATEGORY_TERMS = 80_000_000
MAX_DENSE_TERMS = 300_000_000

# Six construction configurations x two proposal modes.  The last row is the
# one explicit phase-table baseline; all configurations use the same physical
# phase, and only this row pays table construction/storage/lookups separately.
CONSTRUCTION_CONFIGS = (
    ("cut5", 5, True, False),
    ("cut11", 11, False, False),
    ("cut12", 12, True, False),
    ("auto", "auto", False, False),
    ("auto_cache", "auto", True, False),
    ("auto_cache_table", "auto", True, True),
)


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"nested_phase_sampler_comparison_{stamp}.json"
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
    """Separate cumulative scalar and dense ledgers, all reserved first."""

    def __init__(self, cap, *, name):
        self.cap = int(cap)
        self.name = name
        self.reserved = defaultdict(int)
        self.used = defaultdict(int)

    def reserve(self, key, amount=1):
        amount = int(amount)
        if amount < 0 or self.reserved[key] + amount > self.cap:
            raise MemoryError(f"{self.name} category {key} exceeds cap")
        self.reserved[key] += amount

    def charge(self, key, amount=1):
        amount = int(amount)
        if amount < 0 or self.used[key] + amount > self.reserved[key]:
            raise MemoryError(f"{self.name} actual {key} exceeds preflight")
        self.used[key] += amount


def guard_bytes(value, label):
    value = int(value)
    if value < 0 or value > MAX_BYTES:
        raise MemoryError(f"{label} payload {value}>{MAX_BYTES}")
    return value


def g_phase(index):
    return unit_phase(pow(BASE, int(index) % PERIOD, N), N)


def table_phase(table, lookup_counter, index):
    lookup_counter["phase_table_lookups"] += 1
    return table[int(index) % PERIOD]


def aggregate(dst, src):
    for key, value in (src or {}).items():
        if isinstance(value, (int, np.integer)):
            dst[str(key)] += int(value)


def phase_period(cut):
    left = [v for v in POSITIONS if v < int(cut)]
    if not left:
        return 1
    return (1 << max(left)) // math.gcd(PERIOD, 1 << max(left))


def row_bounds(cut):
    """C81 conservative bounds for this fixed r=60,b=3,L=4096 fixture."""
    A = 1 << int(cut)
    K = L // A
    P = phase_period(cut)
    R = min(PERIOD, 2 * BLOCK - 1)
    F = min(P, (A + PERIOD - 1) // PERIOD)
    base = H * BLOCK * BLOCK
    ceil_length = (L + PERIOD - 1) // PERIOD
    C = min(base * K * P, base * ceil_length)
    G = min(H * K * min(A, R * F), H * R * ceil_length, Q)
    return {
        "group_pair_visits": base * K,
        "coefficient_pair_visits": C,
        "geometric_components": G,
        "phase_calls": H + len(POSITIONS) * C,
        "phase_products": (len(POSITIONS) + 1) * C + G,
        "norm_terms": G,
        "K": K, "P": P, "F": F,
    }


def auto_bounds():
    """Uniform auto bounds obtained from min safe (group+C) over candidates."""
    base = H * BLOCK * BLOCK
    ceil_length = (L + PERIOD - 1) // PERIOD
    safe = []
    for cut in sorted({*POSITIONS, SPLIT}):
        K, P = L // (1 << cut), phase_period(cut)
        safe.append(base * K + min(base * K * P, base * ceil_length))
    pair = min(safe)
    coeff = pair - base
    # G<=C and the full-L support bound are both valid.  The first is the
    # selected auto bound; the latter keeps the component cap independent.
    G = min(coeff, H * min(PERIOD, 2 * BLOCK - 1) * ceil_length, Q)
    return {
        "group_pair_visits": pair,
        "coefficient_pair_visits": coeff,
        "geometric_components": G,
        "phase_calls": H + len(POSITIONS) * coeff,
        "phase_products": (len(POSITIONS) + 1) * coeff + G,
        "norm_terms": G,
    }


def config_bound(cut, cache):
    bound = auto_bounds() if cut == "auto" else row_bounds(cut)
    result = dict(bound)
    if cache:
        # The implementation retains a length-b right-product vector per
        # (h,k) slice.  This is a scalar phase cost and a separate live cache.
        K = L // (1 << (SPLIT if cut == "auto" else int(cut)))
        # Auto's selected K is not known before a row; the uniform pair bound
        # safely bounds H*K*b by group_pair_bound//b.
        cache_pairs = (bound["group_pair_visits"] // BLOCK)
        result["phase_calls"] += len(POSITIONS) * cache_pairs
        result["phase_products"] += len(POSITIONS) * cache_pairs
        result["right_cache_entries"] = cache_pairs
    else:
        result["right_cache_entries"] = 0
    result["cut"] = cut
    result["cache_right"] = bool(cache)
    return result


def preflight():
    """Freeze aggregate scalar work, dense reference work and numeric payload."""
    scalar = defaultdict(int)
    dense = defaultdict(int)
    rows_per_config = PERIOD
    for _name, cut, cache, _table in CONSTRUCTION_CONFIGS:
        bound = config_bound(cut, cache)
        for mode in PROPOSALS:
            operation_rows = rows_per_config * (1 + len(OUTPUTS)) + SAMPLE_CALLS
            scalar["helper_constructor_calls"] += 1
            scalar["helper_constructor_matrix_entries"] += 2 * BLOCK * BLOCK
            scalar["helper_constructor_unitary_terms"] += 2 * BLOCK ** 3
            scalar["row_local_terms"] += operation_rows * (
                bound["group_pair_visits"] + bound["coefficient_pair_visits"])
            scalar["geometric_components"] += operation_rows * bound["geometric_components"]
            scalar["phase_calls"] += operation_rows * bound["phase_calls"]
            scalar["phase_products"] += operation_rows * bound["phase_products"]
            scalar["norm_terms"] += operation_rows * bound["norm_terms"]
            scalar["integer_metadata_terms"] += operation_rows * (4+len(POSITIONS))
            scalar["geometry_setup_phase_terms"] += len(POSITIONS)*(3 if cut == "auto" else 1)
            if cut == "auto":
                scalar["selector_pair_terms"] += operation_rows * H * BLOCK ** 2
                scalar["selector_residue_terms"] += operation_rows * H * min(
                    PERIOD, 2 * BLOCK - 1)
                scalar["selector_candidate_terms"] += operation_rows * (len(POSITIONS) + 1)
            # Four public forced queries per work row.  Their Fourier terms are
            # charged separately from row construction.
            scalar["forced_calls"] += rows_per_config * len(OUTPUTS)
            scalar["fourier_component_terms"] += (
                rows_per_config * len(OUTPUTS) * bound["geometric_components"])
            # Four returned-sample requests, each with the frozen cap.  This
            # is a reservation, not a claim that every call will exhaust.
            scalar["sample_calls"] += SAMPLE_CALLS
            scalar["sample_attempts"] += SAMPLE_CALLS * MAX_ATTEMPTS
            scalar["fourier_component_terms"] += (
                SAMPLE_CALLS * MAX_ATTEMPTS * bound["geometric_components"])
            # Proposal preparation and root-mode weighted divisions use the
            # same component-sized scalar workspace and are charged into the
            # Fourier component category.
            scalar["fourier_component_terms"] += operation_rows * bound["geometric_components"]
            if mode == "root_mass":
                scalar["fourier_component_terms"] += (
                    operation_rows + rows_per_config * len(OUTPUTS) + SAMPLE_CALLS * MAX_ATTEMPTS
                ) * bound["geometric_components"]
            scalar["sample_progression_marginal_queries"] += (
                SAMPLE_CALLS * MAX_ATTEMPTS * 2 * WIDTH)
            scalar["sample_progression_proposals"] += SAMPLE_CALLS * MAX_ATTEMPTS
            scalar["sample_work_draws"] += SAMPLE_CALLS
            scalar["sample_component_draws"] += SAMPLE_CALLS * MAX_ATTEMPTS
            scalar["sample_acceptance_draws"] += SAMPLE_CALLS * MAX_ATTEMPTS
            scalar["sample_column_local_terms"] += SAMPLE_CALLS * BLOCK ** 2
            scalar["phase_calls"] += SAMPLE_CALLS * (
                len(POSITIONS) * BLOCK + BLOCK ** 2)
            scalar["phase_products"] += SAMPLE_CALLS * (
                len(POSITIONS) * BLOCK + BLOCK ** 2)
            scalar["right_cache_entries"] += (
                operation_rows * bound["right_cache_entries"])
        if _table:
            scalar["phase_table_values"] += PERIOD
            scalar["phase_table_storage_bytes"] += PERIOD * 16
            scalar["phase_table_modular_pow_queries"] += PERIOD
            cb = config_bound(cut, cache)
            per_mode = (rows_per_config * 5 * cb["phase_calls"]
                        + SAMPLE_CALLS *
                        (len(POSITIONS) * BLOCK + BLOCK ** 2
                         + cb["phase_calls"]))
            scalar["phase_table_lookups"] += 2 * per_mode

    # The dense reference uses four forced target/control marginals and two
    # target returned draws.  All calls use the same supplied branches.
    seq_calls = 4 + 4 + 2
    dense["sequential_calls"] = seq_calls
    # 2 validation + 2 forward + 4 backward matrix products per full step.
    # Unitarity terms below are an overlapping diagnostic subcount, not extra
    # products to add to this count. QR shapes are reported separately.
    dense["sequential_matmul_terms"] = seq_calls * 8 * WIDTH * PERIOD ** 3
    dense["sequential_qr_combined_entries"] = seq_calls * WIDTH * 2 * PERIOD ** 2
    dense["sequential_qr_factor_entries"] = seq_calls * WIDTH * PERIOD ** 2
    dense["sequential_branch_entries"] = 4 * WIDTH * PERIOD ** 2
    dense["sequential_shift_entries"] = 4 * WIDTH * PERIOD ** 2
    dense["sequential_repeated_block_entries"] = 2 * PERIOD ** 2
    dense["sequential_setup_matmul_units"] = (1+2*WIDTH) * PERIOD ** 3
    dense["sequential_shift_entry_writes"] = 2 * WIDTH * PERIOD
    dense["sequential_unitarity_terms"] = seq_calls * 2 * WIDTH * PERIOD ** 3

    # Retain one reference output vector, six-by-two 60-entry work summaries,
    # one phase table, one sampler row and sequential_path's documented QR
    # scratch estimate.  This is numeric payload, not RSS or Python headers.
    payload_parts = {
        "sequential_path_scratch": 16 * (5 * WIDTH + 16) * PERIOD ** 2,
        "sequential_branch_matrices": (
            dense["sequential_branch_entries"]
            + dense["sequential_shift_entries"]
            + dense["sequential_repeated_block_entries"]) * 16,
        "reference_output_vectors": 4 * Q * 8,
        "work_summary_vectors": len(CONSTRUCTION_CONFIGS) * len(PROPOSALS) * 4 * PERIOD * 8,
        "one_sampler_row_and_lists": 16 * (64 * BLOCK ** 2 + 32 * 690 + 256),
        "phase_table": PERIOD * 16,
        "report_scalars": 512 * 1024,
    }
    payload = guard_bytes(sum(payload_parts.values()), "comparison aggregate numeric payload")
    scalar_total = sum(scalar.values())
    dense_total = sum(dense.values())
    if scalar_total > MAX_CATEGORY_TERMS:
        raise MemoryError(f"scalar preflight {scalar_total}>{MAX_CATEGORY_TERMS}")
    if dense_total > MAX_DENSE_TERMS:
        raise MemoryError(f"dense preflight {dense_total}>{MAX_DENSE_TERMS}")
    if any(value > MAX_CATEGORY_TERMS for value in scalar.values()):
        raise MemoryError("scalar category exceeds frozen cap")
    if any(value > MAX_DENSE_TERMS for value in dense.values()):
        raise MemoryError("dense category exceeds frozen cap")
    return {
        "scalar_bounds": dict(scalar), "dense_bounds": dict(dense),
        "scalar_total": scalar_total, "dense_total": dense_total,
        "payload_parts": payload_parts, "planned_numeric_payload_bytes": payload,
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_scalar_category": MAX_CATEGORY_TERMS,
        "max_dense_category": MAX_DENSE_TERMS,
        "configs": [dict(name=n, cut=c, cache_right=ca, table=tb)
                    for n, c, ca, tb in CONSTRUCTION_CONFIGS],
        "proposal_modes": list(PROPOSALS), "sample_calls_per_mode": SAMPLE_CALLS,
        "max_attempts": MAX_ATTEMPTS,
    }


def repeated_block(block):
    out = np.zeros((PERIOD, PERIOD), dtype=np.complex128)
    for cell in range(PERIOD // BLOCK):
        lo = cell * BLOCK
        out[lo:lo + BLOCK, lo:lo + BLOCK] = block
    return out


def shift_matrix(power):
    out = np.zeros((PERIOD, PERIOD), dtype=np.complex128)
    for source in range(PERIOD):
        out[(source + int(power)) % PERIOD, source] = 1.
    return out


def sequential_pairs(phase_diag, late_phase_diag, W1, *, include_early=True):
    """Build dense branches for sequential_path; this is setup, not a propagator."""
    identity = np.eye(PERIOD, dtype=np.complex128)
    after = [identity] * WIDTH
    if include_early:
        after[4] = phase_diag
        after[10] = phase_diag
    # The caller already supplied D_late @ W1. Multiplying W1 again would
    # silently change the physical circuit used as the comparison reference.
    after[11] = late_phase_diag
    pairs = []
    for index in range(WIDTH):
        shift = shift_matrix(1 << index)
        pairs.append((after[index], after[index] @ shift))
    return pairs


def sequential_forced(pairs, initial, outputs, dense_actual, dense_budget):
    result = {}
    for output in outputs:
        dense_budget.charge("sequential_calls")
        dense_budget.charge("sequential_matmul_terms", 8 * WIDTH * PERIOD ** 3)
        dense_budget.charge("sequential_unitarity_terms", 2 * WIDTH * PERIOD ** 3)
        dense_budget.charge("sequential_qr_combined_entries", WIDTH * 2 * PERIOD ** 2)
        dense_budget.charge("sequential_qr_factor_entries", WIDTH * PERIOD ** 2)
        dense_actual["sequential_calls"] += 1
        dense_actual["sequential_matmul_terms"] += 8 * WIDTH * PERIOD ** 3
        dense_actual["sequential_unitarity_terms"] += 2 * WIDTH * PERIOD ** 3
        dense_actual["sequential_qr_combined_entries"] += WIDTH * 2 * PERIOD ** 2
        dense_actual["sequential_qr_factor_entries"] += WIDTH * PERIOD ** 2
        value = sequential_path(pairs, initial, output=int(output),
                                max_payload_bytes=MAX_BYTES)
        skipped_terms = 4*(WIDTH-value["steps"])*PERIOD**3
        dense_actual["sequential_matmul_terms"] -= skipped_terms
        dense_budget.used["sequential_matmul_terms"] -= skipped_terms
        result[int(output)] = float(value["conditional_path_probability"])
    return result


def sequential_draws(pairs, initial, dense_actual, dense_budget):
    rows = []
    for seed in (4105, 4111):
        dense_budget.charge("sequential_calls")
        dense_budget.charge("sequential_matmul_terms", 8 * WIDTH * PERIOD ** 3)
        dense_budget.charge("sequential_unitarity_terms", 2 * WIDTH * PERIOD ** 3)
        dense_budget.charge("sequential_qr_combined_entries", WIDTH * 2 * PERIOD ** 2)
        dense_budget.charge("sequential_qr_factor_entries", WIDTH * PERIOD ** 2)
        dense_actual["sequential_calls"] += 1
        dense_actual["sequential_matmul_terms"] += 8 * WIDTH * PERIOD ** 3
        dense_actual["sequential_unitarity_terms"] += 2 * WIDTH * PERIOD ** 3
        dense_actual["sequential_qr_combined_entries"] += WIDTH * 2 * PERIOD ** 2
        dense_actual["sequential_qr_factor_entries"] += WIDTH * PERIOD ** 2
        value = sequential_path(pairs, initial, rng=np.random.default_rng(seed),
                                max_payload_bytes=MAX_BYTES)
        assert value["steps"] == WIDTH  # full-step named count used above
        rows.append({"seed": seed, "output": int(value["output"]),
                     "path_probability": float(value["conditional_path_probability"]),
                     "steps": int(value["steps"])})
    return rows


def sampler_instance(spec, *, phase_oracle, budget, actual):
    name, cut, cache, use_table = spec
    budget.charge("helper_constructor_calls")
    budget.charge("helper_constructor_matrix_entries", 2 * BLOCK * BLOCK)
    budget.charge("helper_constructor_unitary_terms", 2 * BLOCK ** 3)
    actual["helper_constructor_calls"] += 1
    actual["helper_constructor_matrix_entries"] += 2 * BLOCK * BLOCK
    actual["helper_constructor_unitary_terms"] += 2 * BLOCK ** 3
    W0, W1 = work_block(math.pi / 4), work_block(-math.pi / 10)
    sampler = NestedPhaseProgressions(
        PERIOD, BLOCK, WIDTH, SPLIT, W0, W1, phase_oracle,
        early_phases=((5, phase_oracle), (11, phase_oracle)),
        cut=cut, cache_right=cache, max_local_terms=1_000_000,
        max_components=4096, max_payload_bytes=MAX_BYTES)
    setup = sampler.stats()["geometry_setup_phase_terms"]
    budget.charge("geometry_setup_phase_terms", setup)
    actual["geometry_setup_phase_terms"] += setup
    return sampler


def capture_operation(actual, budget, counters):
    # Production counter names are retained in the report, while the frozen
    # preflight uses shared categories for row/phase/Fourier work.
    category = {
        "column_local_terms": "sample_column_local_terms",
        "row_local_terms": "row_local_terms",
        "phase_queries": "phase_calls",
        "early_phase_queries": "phase_calls",
        "late_phase_queries": "phase_calls",
        "fourier_component_terms": "fourier_component_terms",
        "component_weight_terms": "fourier_component_terms",
        "component_sqrt_terms": "fourier_component_terms",
        "weighted_ratio_divisions": "fourier_component_terms",
        "column_phase_products": "phase_products",
        "row_phase_products": "phase_products",
        "geometric_components": "geometric_components",
        "norm_terms": "norm_terms",
        "right_cache_entries": "right_cache_entries",
        "selector_pair_terms": "selector_pair_terms",
        "selector_residue_terms": "selector_residue_terms",
        "selector_candidate_terms": "selector_candidate_terms",
        "phase_partition_terms": "integer_metadata_terms",
        "progression_proposals": "sample_progression_proposals",
        "marginal_queries": "sample_progression_marginal_queries",
        "work_draws": "sample_work_draws",
        "component_draws": "sample_component_draws",
        "acceptance_draws": "sample_acceptance_draws",
    }
    for key, value in (counters or {}).items():
        if isinstance(value, (int, np.integer)):
            actual[str(key)] += int(value)
            # phase_queries already equals early_phase_queries plus
            # late_phase_queries in the production counter dictionary.
            if str(key) in ("early_phase_queries", "late_phase_queries",
                            "group_pair_visits", "coefficient_pair_visits"):
                continue
            budget.charge(category.get(str(key), str(key)), int(value))


def row_summaries(sampler, mode, scalar_actual, scalar_budget):
    work_rows = []
    marginals = {int(output): 0. for output in OUTPUTS}
    for work in range(PERIOD):
        row = sampler.row(work)
        scalar_actual["integer_metadata_terms"] += 4
        scalar_budget.charge("integer_metadata_terms", 4)
        prepared = sampler._prepare_proposal(row, mode, row["counters"])
        capture_operation(scalar_actual, scalar_budget, row["counters"])
        components = tuple(row["components"])
        m = len(components)
        E = float(prepared["expected_attempts"])
        stride = int(row["stride"])
        gcd_bits = math.gcd(stride, Q).bit_length() - 1
        work_mass = float(row["norm"]) / Q
        expected_fourier = m * E
        expected_progression = 2 * (WIDTH - gcd_bits) * E
        work_rows.append({
            "work": work, "work_mass": work_mass, "m": m,
            "root_envelope_E": E, "expected_fourier_terms": expected_fourier,
            "expected_progression_marginal_queries": expected_progression,
            "stride": stride, "gcd_stride_Q_bits": gcd_bits,
            "cut": int(row["cut"]), "norm": float(row["norm"]),
            "construction_counters": {key: int(row["counters"][key]) for key in (
                "row_local_terms", "group_pair_visits", "coefficient_pair_visits",
                "phase_queries", "row_phase_products", "selector_pair_terms",
                "selector_residue_terms", "selector_candidate_terms", "phase_partition_terms")},
        })
        for output in OUTPUTS:
            # Public forced_joint is intentionally used for the matched query;
            # it rebuilds the row and its counters are retained separately.
            forced = sampler.forced_joint(work, output, proposal=mode)
            capture_operation(scalar_actual, scalar_budget, forced["counters"])
            marginals[int(output)] += float(forced["joint_probability"])
    total = sum(item["work_mass"] for item in work_rows)
    if not math.isfinite(total) or abs(total-1.) > TOL:
        raise ArithmeticError("invalid work-mass total")
    normalized = [item["work_mass"] / total for item in work_rows]
    for item, value in zip(work_rows, normalized, strict=True):
        item["normalized_work_mass"] = value
    return work_rows, marginals, total


def sample_requests(sampler, mode, seed_base, scalar_actual, scalar_budget):
    requests = []
    for index in range(SAMPLE_CALLS):
        seed = int(seed_base + index)
        rng = np.random.default_rng(seed)
        row = {"request": index, "seed": seed, "proposal": mode,
               "max_attempts": MAX_ATTEMPTS}
        try:
            result = sampler.sample(rng, max_attempts=MAX_ATTEMPTS, proposal=mode)
        except Exception as exc:
            counters = sampler.last_counters or {}
            capture_operation(scalar_actual, scalar_budget, counters)
            row.update({"status": "failed", "exception": repr(exc),
                        "counters": safe(dict(counters))})
        else:
            counters = result.get("counters", sampler.last_counters or {})
            capture_operation(scalar_actual, scalar_budget, counters)
            row.update({"status": "success", "work": int(result["work"]),
                        "output": int(result["output"]),
                        "attempts": int(result["attempts"]),
                        "component_count": int(result["component_count"]),
                        "counters": safe(dict(counters))})
        requests.append(row)
    return requests


def main():
    exp = Experiment("nested_phase_sampler_comparison", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "forced helper marginals match the dense sequential reference")
    exp.predict("P2", "all rows normalize and expose the frozen C80 cost quantities")
    exp.predict("P3", "bounded returned calls retain successes and exhausted failures")
    exp.predict("P4", "the fixed phase-table baseline has the same matched output")
    exp.must_fail("C1", "omitting both early phases changes a checked output marginal")
    started = time.perf_counter()
    report = {"status": "FAIL"}
    scalar_actual = defaultdict(int)
    dense_actual = defaultdict(int)
    scalar_budget = Budget(MAX_CATEGORY_TERMS, name="scalar")
    dense_budget = Budget(MAX_DENSE_TERMS, name="dense")
    p1 = p2 = p3 = p4 = c1 = False
    rows, sample_rows = {}, {}
    report["rows"], report["sample_requests"] = rows, sample_rows
    try:
        pre = preflight()
        report["preflight"] = pre
        for key, value in pre["scalar_bounds"].items():
            scalar_budget.reserve(key, value)
        for key, value in pre["dense_bounds"].items():
            dense_budget.reserve(key, value)

        table = tuple(complex(g_phase(index)) for index in range(PERIOD))
        scalar_actual["phase_table_values"] = PERIOD
        scalar_actual["phase_table_storage_bytes"] = PERIOD * 16
        scalar_actual["phase_table_modular_pow_queries"] = PERIOD
        table_lookup_actual = defaultdict(int)
        scalar_budget.charge("phase_table_values", PERIOD)
        scalar_budget.charge("phase_table_storage_bytes", PERIOD * 16)
        scalar_budget.charge("phase_table_modular_pow_queries", PERIOD)

        W0 = repeated_block(work_block(math.pi / 4))
        W1 = repeated_block(work_block(-math.pi / 10))
        initial = W0[:, 0]
        D = np.diag(np.asarray(table, dtype=np.complex128))
        late = D @ W1
        dense_budget.charge("sequential_setup_matmul_units", (1+2*WIDTH) * PERIOD ** 3)
        dense_budget.charge("sequential_branch_entries", 4 * WIDTH * PERIOD ** 2)
        dense_budget.charge("sequential_shift_entries", 4 * WIDTH * PERIOD ** 2)
        dense_budget.charge("sequential_repeated_block_entries", 2 * PERIOD ** 2)
        dense_budget.charge("sequential_shift_entry_writes", 2 * WIDTH * PERIOD)
        dense_actual["sequential_setup_matmul_units"] += (1+2*WIDTH) * PERIOD ** 3
        dense_actual["sequential_branch_entries"] += 4 * WIDTH * PERIOD ** 2
        dense_actual["sequential_shift_entries"] += 4 * WIDTH * PERIOD ** 2
        dense_actual["sequential_repeated_block_entries"] += 2 * PERIOD ** 2
        dense_actual["sequential_shift_entry_writes"] += 2 * WIDTH * PERIOD
        target_pairs = sequential_pairs(D, late, W1, include_early=True)
        control_pairs = sequential_pairs(D, late, W1, include_early=False)
        target_reference = sequential_forced(
            target_pairs, initial, OUTPUTS, dense_actual, dense_budget)
        control_reference = sequential_forced(
            control_pairs, initial, OUTPUTS, dense_actual, dense_budget)
        sequential_samples = sequential_draws(
            target_pairs, initial, dense_actual, dense_budget)

        helper_differences = []
        for config_index, spec in enumerate(CONSTRUCTION_CONFIGS):
            name, cut, cache, use_table = spec
            oracle = (lambda index, values=table, counters=table_lookup_actual:
                      table_phase(values, counters, index)) if use_table else g_phase
            for mode_index, mode in enumerate(PROPOSALS):
                sampler = sampler_instance(
                    spec, phase_oracle=oracle, budget=scalar_budget,
                    actual=scalar_actual)
                # Constructor bounds are frozen before matrices are copied;
                # this report records the resulting structural contract.
                row_data, marginals, work_total = row_summaries(
                    sampler, mode, scalar_actual, scalar_budget)
                key = f"{name}/proposal={mode}"
                output_errors = {
                    str(output): abs(marginals[int(output)] - target_reference[int(output)])
                    for output in OUTPUTS}
                helper_differences.extend(output_errors.values())
                rows[key] = {
                    "construction": {"name": name, "cut": cut,
                                     "cache_right": cache, "phase_table": use_table},
                    "proposal": mode,
                    "stats": safe(sampler.stats()),
                    "work_rows": row_data,
                    "work_mass_sum": work_total,
                    "output_marginals": {str(k): float(v) for k, v in marginals.items()},
                    "target_reference_marginals": {str(k): float(v)
                                                    for k, v in target_reference.items()},
                    "max_forced_marginal_error": max(output_errors.values()),
                }
                # Common seed exponents/work draws make these four requests
                # comparable; no seedwise or statistical speedup is predicted.
                seed_base = 510000 + mode_index * 10
                sample_rows[key] = sample_requests(
                    sampler, mode, seed_base, scalar_actual, scalar_budget)

        scalar_actual["phase_table_lookups"] = table_lookup_actual[
            "phase_table_lookups"]
        scalar_budget.charge("phase_table_lookups",
                             scalar_actual["phase_table_lookups"])

        control_differences = {
            str(output): abs(target_reference[int(output)] - control_reference[int(output)])
            for output in OUTPUTS}
        c1 = max(control_differences.values()) > 1e-8
        p1 = (bool(rows) and max(helper_differences, default=float("inf")) < TOL
              and all(math.isfinite(v) for v in helper_differences))
        p2 = all(
            math.isfinite(float(item["work_mass"])) and item["work_mass"] >= 0
            and math.isfinite(float(item["root_envelope_E"]))
            and 1-1e-12 <= float(item["root_envelope_E"]) <= item["m"]*(1+1e-12)
            and math.isfinite(float(item["expected_fourier_terms"]))
            and math.isfinite(float(item["expected_progression_marginal_queries"]))
            for detail in rows.values() for item in detail["work_rows"])
        all_requests = [item for values in sample_rows.values() for item in values]
        successes = [item for item in all_requests if item["status"] == "success"]
        failures = [item for item in all_requests if item["status"] == "failed"]
        p3 = (len(all_requests) == len(CONSTRUCTION_CONFIGS) * len(PROPOSALS) * SAMPLE_CALLS
              and all(1 <= item["attempts"] <= MAX_ATTEMPTS for item in successes)
              and all("counters" in item for item in all_requests)
              and all(item["counters"]["work_draws"] == 1 for item in all_requests)
              and all("cap exhausted" in item["exception"]
                      and item["counters"]["acceptance_draws"] == MAX_ATTEMPTS
                      for item in failures)
              and all(0 <= item["work"] < PERIOD and 0 <= item["output"] < Q
                      for item in successes)
              and len(successes) + len(failures) == len(all_requests))
        table_key = "auto_cache_table/proposal=mass"
        p4 = (all(rows[f"auto_cache_table/proposal={mode}"]["max_forced_marginal_error"] < TOL
                  for mode in PROPOSALS)
              and scalar_actual["phase_table_values"] == PERIOD
              and scalar_actual["phase_table_storage_bytes"] == PERIOD * 16
              and scalar_actual["phase_table_lookups"] > 0)

        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "width": WIDTH, "split": SPLIT, "Q": Q, "L": L, "H": H,
                        "positions": list(POSITIONS)},
            "rows": rows,
            "sample_requests": sample_rows,
            "sequential": {
                "target_reference_marginals": target_reference,
                "omitted_early_phase_reference_marginals": control_reference,
                "omission_control_abs_differences": control_differences,
                "returned_draws": sequential_samples,
                "branch_dimension": PERIOD,
                "branch_count": WIDTH * 2,
                "qr_shapes": {"forward_combined": [PERIOD, 2 * PERIOD],
                              "forward_factor": [PERIOD, PERIOD]},
            },
            "scalar_actual": dict(scalar_actual),
            "dense_actual": dict(dense_actual),
            "successful_sample_count": len(successes),
            "failed_sample_count": len(failures),
            "elapsed_seconds": time.perf_counter() - started,
            "scope": "bounded returned-sample work comparison; no timing claim",
        })
        summary = {}
        for key, detail in rows.items():
            totals = defaultdict(int)
            for request in sample_rows[key]:
                aggregate(totals, request["counters"])
            succeeded = sum(request["status"] == "success" for request in sample_rows[key])
            weighted = lambda field: sum(item["normalized_work_mass"]*item[field]
                                         for item in detail["work_rows"])
            construction = {field: sum(item["normalized_work_mass"]*
                                      item["construction_counters"][field]
                                      for item in detail["work_rows"])
                            for field in detail["work_rows"][0]["construction_counters"]}
            summary[key] = dict(
                requested=SAMPLE_CALLS, returned=succeeded, exhausted=SAMPLE_CALLS-succeeded,
                expected_uncapped_attempts=weighted("root_envelope_E"),
                expected_uncapped_fourier_terms=weighted("expected_fourier_terms"),
                expected_uncapped_marginal_queries=weighted("expected_progression_marginal_queries"),
                expected_row_construction=construction,
                all_request_counters=dict(totals),
                counters_per_request={k: v/SAMPLE_CALLS for k, v in totals.items()},
                counters_per_return_including_failures={k: v/succeeded if succeeded else None
                                                       for k, v in totals.items()})
        report["comparison_summary"] = summary
        report["cost_semantics"] = (
            "C80 envelopes give uncapped mathematical expectations from float row masses. "
            "Observed counters include failed requests; four seeds are not a runtime estimate. "
            "Dense product/QR counts and scalar geometric-query counts are different operations. "
            "Unitary-validation terms are a subset of dense matmul terms.")
        report["status"] = "PASS" if p1 and p2 and p3 and p4 and c1 else "FAIL"
        exp.check("P1", p1, "helper all-work forced marginals match sequential reference")
        exp.check("P2", p2, "row masses, C80 envelopes and expected query costs are finite")
        exp.check("P3", p3, "all bounded sample requests retain success/failure records")
        exp.check("P4", p4, "explicit phase-table baseline preserves matched output")
        exp.fail_check("C1", c1, "omitting early phases changes the checked reference output")
    except Exception as exc:
        report["exception"] = repr(exc)
        report["traceback"] = traceback.format_exc()
        exp.check("P1", False, "exception before matched output comparison")
        exp.check("P2", False, "exception before row metrics")
        exp.check("P3", False, "exception before returned samples")
        exp.check("P4", False, "exception before phase-table baseline")
        exp.fail_check("C1", False, "exception before omission control")

    report["scalar_actual"] = dict(scalar_actual)
    report["dense_actual"] = dict(dense_actual)
    report["scalar_reserved"] = dict(scalar_budget.reserved)
    report["scalar_budget_used"] = dict(scalar_budget.used)
    report["dense_reserved"] = dict(dense_budget.reserved)
    report["dense_budget_used"] = dict(dense_budget.used)
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                    "width": WIDTH, "split": SPLIT},
        "reference": "existing lab.semiclassical.sequential_path with supplied 60x60 branches",
        "same_physical_state": True,
        "returned_sample_calls_per_configuration": SAMPLE_CALLS,
        "max_attempts": MAX_ATTEMPTS,
        "failed_calls_retained": True,
        "no_q_by_r_archive": True,
        "no_timing_claim": True,
        "float_diagnostic": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
