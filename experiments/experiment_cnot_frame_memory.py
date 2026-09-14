"""Bounded allocation audit for TODO43's lazy CNOT coordinate frame.

This integrates a worker proposal with coordinator-reviewed accounting. It compares the
existing propagator with and without ``affine_frame=True`` on the same output,
then audits three specialized CNOT-boundary baselines.  The boundary kernels
are not alternative propagators: they begin with an already-built coefficient
dictionary and contain no nonlinear branching, truncation, or tracing.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  On a fixed 4^6-term nonlinear suffix, adding a nonempty CNOT prefix
      preserves support; physical results agree BETWEEN modes at each depth. The lazy
      frame has a lower traced whole-call peak than the unchanged loop because
      the unchanged CNOT step retains old, relabelled, and coefficient-filtered
      dictionaries simultaneously.
  P2  CNOT depth changes allocation churn but not support cardinality.  For D
      CNOTs over S retained terms the unchanged loop performs 2DS CNOT-boundary
      dictionary insertions, while the frame performs 2D frame-vector XORs.
  P3  Streaming a framed result allocates no support-sized collection.  Asking
      for ``dict(result.final_terms)`` explicitly adds one physical-key table;
      that materialization must be charged together with the retained backing
      coordinate dictionary and frame.
  P4  The strongest inexpensive pure-boundary comparators are (a) in-place
      pair swaps with an S-reference key snapshot and (b) one batched affine
      materialization.  Both return the same physical dictionary, but neither
      carries a frame through an interleaved nonlinear gate.  A packed uint64
      key array is recorded as a representation-boundary lower comparator and
      is not a same-API dictionary propagator.
  P5  A bounded nonlinear/physical-key sentinel with max_weight and trace_plus
      agrees exactly between modes, including physical final keys and trace
      events.
  P6  At D=128, an existing-engine suffix plus one batched CNOT materialization
      returns the same full physical dictionary. Its whole-call peak, and the
      framed-plus-materialization whole-call peak, are open comparisons: either
      can erase the advantage over the unchanged engine. Setup and output
      allocation are inside each trace; prebuilt circuit objects are outside.
  C1  The one-dictionary final-memory property must fail after explicit
      physical-key materialization: the added table must have positive
      support-sized retained bytes.

No native/RSS or timing claim is made.  ``tracemalloc`` measures traced Python
allocations for these bounded calls.  Shallow retained layouts are itemized
separately; allocator arenas, Circuit storage, interpreter state, and RSS are
outside that accounting.
"""
from __future__ import annotations

import gc
import json
import sys
import tracemalloc
from pathlib import Path

import numpy as np

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESEARCH_ROOT))

from circuits import Circuit
from lab.harness import Experiment
from perm_pps import FramedTerms, propagate_perm


OUT = RESEARCH_ROOT / "out"
M = 6
N = 3 * M
DEPTHS = (0, 1, 8, 32, 128)


def cnot_schedule(depth: int, n: int = N) -> list[tuple[int, int]]:
    """A prefix of one fixed deterministic schedule; depth is the only sweep."""
    ops = []
    for step in range(depth):
        control = (7 * step + 1) % n
        target = (11 * step + 5) % n
        if control == target:
            target = (target + 1) % n
        ops.append((control, target))
    return ops


def fixture(depth: int) -> tuple[Circuit, int]:
    """Forward CNOT prefix, then disjoint Toffolis; reverse support is 4^M."""
    qc = Circuit(N)
    for control, target in cnot_schedule(depth):
        qc.cnot(control, target)
    zmask = 0
    for block in range(M):
        a, b, target = 3 * block, 3 * block + 1, 3 * block + 2
        qc.toffoli(a, b, target)
        zmask |= 1 << target
    return qc, zmask


def traced_propagation(qc: Circuit, zmask: int, framed: bool):
    gc.collect()
    tracemalloc.start()
    result = propagate_perm(qc, zmask, affine_frame=framed)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, current, peak


def dict_layout(terms: dict[int, float]) -> dict[str, int]:
    return {
        "entries": len(terms),
        "table_bytes": sys.getsizeof(terms),
        "unique_key_bytes": sum(sys.getsizeof(key) for key in terms),
        "unique_value_bytes": sum(
            sys.getsizeof(value) for value in {id(v): v for v in terms.values()}.values()
        ),
    }


def framed_layout(view: FramedTerms) -> dict[str, object]:
    backing = view._terms
    frame = view._frame
    vectors = {id(v): v for v in (*frame.rows, *frame.inverse_columns)}
    return {
        "view_bytes": sys.getsizeof(view),
        "view_dict_bytes": sys.getsizeof(view.__dict__),
        "backing": dict_layout(backing),
        "frame_object_bytes": sys.getsizeof(frame),
        "frame_dict_bytes": sys.getsizeof(frame.__dict__),
        "frame_list_bytes": sys.getsizeof(frame.rows) + sys.getsizeof(frame.inverse_columns),
        "frame_unique_vector_bytes": sum(sys.getsizeof(v) for v in vectors.values()),
        "frame_vector_count": len(vectors),
    }


def incremental_access(view: FramedTerms, materialize: bool) -> tuple[int, int, object]:
    """Traced bytes allocated after the already-resident framed result exists."""
    gc.collect()
    tracemalloc.start()
    tracemalloc.reset_peak()
    if materialize:
        value = dict(view)
    else:
        value = sum((key + 1) * coefficient for key, coefficient in view.items())
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return current, peak, value


def relabel(key: int, control: int, target: int) -> int:
    return key ^ (((key >> target) & 1) << control)


def boundary_rebuild(seed: dict[int, float], ops: list[tuple[int, int]]):
    terms = dict(seed)
    for control, target in reversed(ops):
        terms = {relabel(key, control, target): value for key, value in terms.items()}
    return terms


def boundary_in_place(seed: dict[int, float], ops: list[tuple[int, int]]):
    """Relabel one involution using one key-reference list and pair swaps.

    Pair-complete supports only swap values.  Singleton pairs pop/reinsert a
    key and may accumulate dummy slots or trigger a dict resize; this is why
    measured peak and final table size are reported rather than inferred.
    """
    terms = dict(seed)
    for control, target in reversed(ops):
        keys = list(terms)
        control_bit = 1 << control
        for key in keys:
            if ((key >> target) & 1) and not (key & control_bit):
                partner = key ^ control_bit
                if partner in terms:
                    terms[key], terms[partner] = terms[partner], terms[key]
                else:
                    terms[partner] = terms.pop(key)
        for key in keys:
            if ((key >> target) & 1) and (key & control_bit):
                partner = key ^ control_bit
                if partner not in terms:
                    terms[partner] = terms.pop(key)
    return terms


def physical_from_rows(key: int, rows: list[int]) -> int:
    return sum(((row & key).bit_count() & 1) << q for q, row in enumerate(rows))


def boundary_batched(seed: dict[int, float], ops: list[tuple[int, int]], n: int = N):
    rows = [1 << q for q in range(n)]
    for control, target in reversed(ops):
        rows[control] ^= rows[target]
    return {physical_from_rows(key, rows): value for key, value in seed.items()}


def boundary_packed(seed: dict[int, float], ops: list[tuple[int, int]]):
    keys = np.fromiter(seed, dtype=np.uint64, count=len(seed))
    values = np.fromiter(seed.values(), dtype=np.float64, count=len(seed))
    for control, target in reversed(ops):
        keys ^= ((keys >> np.uint64(target)) & np.uint64(1)) << np.uint64(control)
    return keys, values


def boundary_packed_chunked(seed: dict[int, float], ops: list[tuple[int, int]]):
    """Same packed boundary, with RHS temporaries bounded to 256 entries."""
    keys = np.fromiter(seed, dtype=np.uint64, count=len(seed))
    values = np.fromiter(seed.values(), dtype=np.float64, count=len(seed))
    for control, target in reversed(ops):
        for start in range(0, len(keys), 256):
            chunk = keys[start:start + 256]
            chunk ^= ((chunk >> np.uint64(target)) & np.uint64(1)) << np.uint64(control)
    return keys, values


def traced_complete_boundary(qc, zmask, *, ops=None, materialize=False):
    """Include propagation and final physical dictionary in ONE trace.

    ops supplied: qc is the fixed nonlinear suffix, processed by the EXISTING
    default propagator, then one batched linear remap replaces final_terms.
    ops absent: qc is the complete circuit, processed by the framed mode.
    This comparator applies only to a CNOT prefix with no intermediate filters.
    """
    gc.collect()
    tracemalloc.start()
    result = propagate_perm(qc, zmask, affine_frame=ops is None)
    if ops is not None:
        result.final_terms = boundary_batched(result.final_terms, ops, qc.n)
        result.n_terms.extend([len(result.final_terms)] * len(ops))
    elif materialize:
        # Retain the backing result while the requested dict is consumed, as a
        # caller normally would. This measures the more demanding boundary.
        physical = dict(result.final_terms)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, current, peak, physical if ops is None and materialize else None


def traced_boundary(function, seed, ops):
    gc.collect()
    tracemalloc.start()
    result = function(seed, ops)
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return result, current, peak


exp = Experiment("experiment_cnot_frame_memory", doc=__doc__)
exp.predict("P1", "fixed nonlinear support and physical law; lower framed whole-call peak for D>0")
exp.predict("P2", "support fixed while CNOT insertion churn is 2*D*S versus 2*D frame XORs")
exp.predict("P3", "streaming has sub-support extra allocation; dict(view) adds a physical-key table")
exp.predict("P4", "specialized in-place, batched and packed boundaries preserve the physical law")
exp.predict("P5", "physical max_weight and trace_plus sentinel agrees exactly")
exp.predict("P6", "complete batched and frame-plus-materialization calls return the same law; memory ordering is open")
exp.must_fail("C1", "explicit dict(view) must violate the one-dictionary final-memory property")

rows = []
whole_results = {}
for depth in DEPTHS:
    qc, zmask = fixture(depth)
    legacy, legacy_current, legacy_peak = traced_propagation(qc, zmask, False)
    framed, framed_current, framed_peak = traced_propagation(qc, zmask, True)
    physical = dict(framed.final_terms)
    same = (
        physical == legacy.final_terms
        and framed.n_terms == legacy.n_terms
        and framed.expectation == legacy.expectation
        and framed.frame_cnot_updates == depth
    )
    support = len(physical)
    rows.append({
        "kind": "whole_propagation",
        "depth": depth,
        "support": support,
        "same_physical_result": same,
        "legacy_current_traced_bytes": legacy_current,
        "legacy_peak_traced_bytes": legacy_peak,
        "framed_current_traced_bytes": framed_current,
        "framed_peak_traced_bytes": framed_peak,
        "derived_unchanged_loop_cnot_entry_insertions": 2 * depth * support,
        "derived_frame_vector_xors": 2 * depth,
    })
    whole_results[depth] = (legacy, framed, legacy_peak, framed_peak)

support_counts = {row["support"] for row in rows}
same_all = all(row["same_physical_result"] for row in rows)
peak_advantage = all(
    row["framed_peak_traced_bytes"] < row["legacy_peak_traced_bytes"]
    for row in rows if row["depth"] > 0
)
exp.check("P1", same_all and support_counts == {4 ** M} and peak_advantage,
          f"supports={sorted(support_counts)}, same={same_all}, lower_peak={peak_advantage}")
exp.check("P2", support_counts == {4 ** M},
          f"D=128: insertions={2*128*(4**M)}, frame_xors={2*128}")

_, final_view, _, _ = whole_results[DEPTHS[-1]]
assert isinstance(final_view.final_terms, FramedTerms)
stream_current, stream_peak, stream_checksum = incremental_access(final_view.final_terms, False)
material_current, material_peak, materialized = incremental_access(final_view.final_terms, True)
backing_layout = framed_layout(final_view.final_terms)
material_layout = dict_layout(materialized)
same_value_ids = {id(v) for v in materialized.values()} <= {
    id(v) for v in final_view.final_terms._terms.values()
}
rows.append({
    "kind": "final_access",
    "depth": DEPTHS[-1],
    "support": len(final_view.final_terms),
    "stream_incremental_current_bytes": stream_current,
    "stream_incremental_peak_bytes": stream_peak,
    "stream_checksum": stream_checksum,
    "materialize_incremental_current_bytes": material_current,
    "materialize_incremental_peak_bytes": material_peak,
    "backing_and_frame_shallow_layout": backing_layout,
    "materialized_shallow_layout": material_layout,
    "materialized_values_shared_with_backing": same_value_ids,
})
exp.check("P3", stream_peak < material_peak and same_value_ids,
          f"stream_peak={stream_peak}, material_peak={material_peak}, shared_values={same_value_ids}")
exp.fail_check("C1", material_current >= material_layout["table_bytes"] > 0,
               f"additional_current={material_current}, table={material_layout['table_bytes']}")

seed_qc, seed_mask = fixture(0)
seed = dict(propagate_perm(seed_qc, seed_mask).final_terms)
ops = cnot_schedule(DEPTHS[-1])
reference, rebuild_current, rebuild_peak = traced_boundary(boundary_rebuild, seed, ops)
in_place, in_place_current, in_place_peak = traced_boundary(boundary_in_place, seed, ops)
batched, batched_current, batched_peak = traced_boundary(boundary_batched, seed, ops)
(packed_keys, packed_values), packed_current, packed_peak = traced_boundary(boundary_packed, seed, ops)
packed = dict(zip(map(int, packed_keys), map(float, packed_values)))
(chunk_keys, chunk_values), chunk_current, chunk_peak = traced_boundary(boundary_packed_chunked, seed, ops)
chunked = dict(zip(map(int, chunk_keys), map(float, chunk_values)))
baseline_same = reference == in_place == batched == packed == chunked
rows.append({
    "kind": "specialized_cnot_boundary",
    "depth": DEPTHS[-1],
    "support": len(seed),
    "same_physical_result": baseline_same,
    "rebuild_current_traced_bytes": rebuild_current,
    "rebuild_peak_traced_bytes": rebuild_peak,
    "in_place_current_traced_bytes": in_place_current,
    "in_place_peak_traced_bytes": in_place_peak,
    "in_place_final_table_bytes": sys.getsizeof(in_place),
    "batched_current_traced_bytes": batched_current,
    "batched_peak_traced_bytes": batched_peak,
    "packed_current_traced_bytes": packed_current,
    "packed_peak_traced_bytes": packed_peak,
    "packed_payload_bytes": packed_keys.nbytes + packed_values.nbytes,
    "packed_chunked_current_traced_bytes": chunk_current,
    "packed_chunked_peak_traced_bytes": chunk_peak,
    "packed_chunk_entries": 256,
    "scope": "prebuilt coefficient dictionary; CNOT-only boundary; packed form lacks dict API/nonlinear merge",
})
exp.check("P4", baseline_same and len(reference) == 4 ** M,
          f"same={baseline_same}, peaks rebuild/in-place/batched/packed="
          f"{rebuild_peak}/{in_place_peak}/{batched_peak}/{packed_peak}")

full_qc, full_mask = fixture(DEPTHS[-1])
complete_batched, batch_current, batch_peak, _ = traced_complete_boundary(
    seed_qc, seed_mask, ops=ops)
complete_framed, complete_current, complete_peak, complete_physical = traced_complete_boundary(
    full_qc, full_mask, materialize=True)
complete_same = (complete_batched.final_terms == complete_physical == reference
                 and complete_batched.n_terms == complete_framed.n_terms
                 and complete_batched.expectation == complete_framed.expectation)
rows.append({
    "kind": "complete_materialized_comparison",
    "depth": DEPTHS[-1],
    "support": len(reference),
    "same_physical_result": complete_same,
    "batched_current_traced_bytes": batch_current,
    "batched_peak_traced_bytes": batch_peak,
    "framed_plus_materialization_current_traced_bytes": complete_current,
    "framed_plus_materialization_peak_traced_bytes": complete_peak,
    "scope": "same full physical dictionary; batching specialized to this CNOT-prefix fixture; original circuit objects excluded in all modes",
})
exp.check("P6", complete_same,
          f"complete batched peak={batch_peak}; frame plus requested dict peak={complete_peak}")

sentinel = (Circuit(7)
            .cnot(0, 4)
            .toffoli(1, 2, 4)
            .cnot(3, 1)
            .toffoli(0, 5, 2)
            .cnot(4, 6)
            .toffoli(2, 3, 6)
            .x(1))
sentinel_args = dict(zmask=(1 << 2) | (1 << 4) | (1 << 6), delta=0.10,
                     max_weight=4, trace_plus=(3,))
plain_sentinel = propagate_perm(sentinel, **sentinel_args)
frame_sentinel = propagate_perm(sentinel, **sentinel_args, affine_frame=True)
sentinel_same = (
    dict(frame_sentinel.final_terms) == plain_sentinel.final_terms
    and frame_sentinel.n_terms == plain_sentinel.n_terms
    and frame_sentinel.trace_events == plain_sentinel.trace_events
    and frame_sentinel.expectation == plain_sentinel.expectation
)
sentinel_nontrivial = (
    0 < len(plain_sentinel.final_terms) < (1 << sentinel.n)
    and any(before > after for _, _, before, after in plain_sentinel.trace_events)
)
rows.append({
    "kind": "nonlinear_weight_trace_sentinel",
    "support": len(plain_sentinel.final_terms),
    "same_physical_result": sentinel_same,
    "nontrivial": sentinel_nontrivial,
    "n_terms": plain_sentinel.n_terms,
    "trace_events": plain_sentinel.trace_events,
})
exp.check("P5", sentinel_same and sentinel_nontrivial,
          f"same={sentinel_same}, support={len(plain_sentinel.final_terms)}, "
          f"trace={plain_sentinel.trace_events}")

exp.finish(
    report_path=OUT / "cnot_frame_memory_main.json",
    rows=rows,
    metadata={
        "native_memory_measured": False,
        "timing_claim": False,
        "fixture_terms": 4 ** M,
        "max_cnot_depth": DEPTHS[-1],
        "counts_derived_not_instrumented": ["CNOT entry insertions", "frame XORs"],
        "exclusions": ["prebuilt circuit objects", "preexisting reference results", "interpreter and allocator arenas", "RSS"],
    },
)
