"""Conditional Walsh queries by constant-workspace itinerary replay.

Frozen contract: out/agent-board/workers/A050923a6a5524fc3/reconciled_contract.md.
This module supplies the candidate, shared-lazy C89 baseline and cold-child
measurement interface. A separate versioned harness validates the candidate
against the frozen blind native gate reference before running performance.
No full outside average, intermediate Pauli support or general memory win is
implied. Arithmetic and the four-state range DP are reused from C89.
"""
from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass
import json
import resource
import statistics
import time
import tracemalloc

from experiments.experiment_controlled_intervals import (
    compile_intervals, macro_image, prefix_walsh,
)


@dataclass(frozen=True, slots=True)
class NativeWords(Sequence):
    """An immutable O(m+log q)-bit view of native a=1 power-of-two words."""

    q: int
    x: int = 0
    u: int = 1
    controls: bool = False
    reverse: bool = False

    def __len__(self):
        return self.q

    def __getitem__(self, index):
        if not isinstance(index, int):
            raise TypeError("integer indexing only; no materialized slices")
        if index < 0:
            index += self.q
        if not 0 <= index < self.q:
            raise IndexError(index)
        i = self.q - 1 - index if self.reverse else index
        return self.u * ((self.x >> i) & 1) if self.controls else 1 << i


def validate(width, modulus, constants, t, h, controls, queries=()):
    if not 1 <= width <= 4096 or len(constants) != len(controls) or len(constants) > 256:
        raise ValueError("bounded equal-length constant/control sequences required")
    size = 1 << width
    if not 0 < modulus < size or not 0 <= t < size or h not in (0, 1):
        raise ValueError("invalid fixed scratch/modulus")
    if any(not 0 <= c < size for c in constants) or any(e not in (0, 1) for e in controls):
        raise ValueError("invalid constant/control")
    if any(not (0 <= a < 2*size and 0 <= c < 2*size) for a, c in queries):
        raise ValueError("mask out of range")


def local_span(width, modulus, constant, t, h, control, z, *, omit_wrap=False):
    """Next local ordinary-translation boundary; constant many residues only.

    omit_wrap is a deliberately wrong, bounded correctness control. It is never
    selected by either benchmark method.
    """
    size, half = 1 << width, 1 << (width-1)
    b = z % size
    r = modulus - 2*(t & modulus)
    s, d = control*(constant - 2*(t & constant)), t+h
    span = size-b  # lane endpoint also handles the residue zero
    for raw in (half, r, r+half, r-s, r-s+half):
        cut = raw % size
        if cut > b:
            span = min(span, cut-b)
    if not omit_wrap:
        for raw in (-d-s+r, -d-s):
            cut = raw % size
            if cut > b:
                span = min(span, cut-b)
    return span


def segment_at(width, modulus, constants, t, h, controls, lo, *, wrong_min=False):
    """Replay the itinerary of lo while retaining every prior distance cap."""
    span, z = (1 << (width+1))-lo, lo
    for c, e in zip(constants, controls):
        cap = local_span(width, modulus, c, t, h, e, z)
        span = cap if wrong_min else min(span, cap)
        z = macro_image(width, modulus, c, t, h, e, z)
    return lo, lo+span, z-lo


def stream_pieces(width, modulus, constants, t, h, controls, stats=None):
    """Yield maximal final translations using one pending segment triple.

    The unmerged itinerary has <=1+15q pieces. Since all retained points follow
    the same local translation at each step, each emitted interval is exact.
    Adjacent equal shifts are merged before the Walsh DP, matching C89 _merge.
    stats, if supplied, contains only two counters, never a history or table.
    """
    validate(width, modulus, constants, t, h, controls)
    size, lo, raw_count, emitted = 1 << (width+1), 0, 0, 0
    pending = None
    while lo < size:
        piece = segment_at(width, modulus, constants, t, h, controls, lo)
        _, hi, shift = piece
        if not lo < hi <= size or not 0 <= lo+shift < hi+shift <= size:
            raise AssertionError("invalid streamed interval")
        raw_count += 1
        if raw_count > 1+15*len(constants):
            raise AssertionError("itinerary bound violated")
        if pending is not None and pending[2] == shift:
            pending = (pending[0], hi, shift)
        else:
            if pending is not None:
                emitted += 1
                yield pending
            pending = piece
        lo = hi
    if pending is not None:
        emitted += 1
        yield pending
    if stats is not None:
        stats.update(raw_pieces=raw_count, emitted_pieces=emitted)


def contract_pieces(width, pieces, queries):
    """Retain k integer numerators, using the same C89 DP for both methods."""
    sums = [0] * len(queries)
    for lo, hi, shift in pieces:
        for j, (alpha, beta) in enumerate(queries):
            sums[j] += (prefix_walsh(width+1, hi, shift, alpha, beta)
                        - prefix_walsh(width+1, lo, shift, alpha, beta))
    return tuple(sums)


def evaluate(method, width, modulus, constants, t, h, controls, queries):
    validate(width, modulus, constants, t, h, controls, queries)
    stats = {}
    if method == "stream":
        pieces = stream_pieces(width, modulus, constants, t, h, controls, stats)
    elif method == "table":
        pieces = compile_intervals(width, modulus, constants, t, h, controls)
        stats["emitted_pieces"] = len(pieces)
    elif method == "noop":
        return (0, 0, 0, 0), stats
    else:
        raise ValueError("unknown method")
    return contract_pieces(width, pieces, queries), stats


def frozen_inputs(q):
    if q not in (1, 3, 7, 15, 31, 63, 127):
        raise ValueError("q not in the frozen campaign")
    width = 129
    modulus = (1 << 128)-1
    t = (1 << 128) | int("96"*16, 16)
    x = int("69"*16, 16)
    constants = NativeWords(q)
    controls = NativeWords(q, x=x, u=1, controls=True)
    pairs = (
        ((128,), (129, 0)),
        ((128, 17, 73), (129, 5, 96)),
        ((128, 3, 64, 111), (129, 29, 89, 127)),
        ((128, 11, 47, 101, 127), (129, 1, 32, 67, 120)),
    )
    queries = tuple((sum(1 << i for i in a), sum(1 << i for i in c)) for a, c in pairs)
    return width, modulus, constants, t, 1, controls, queries


def whole_operation(method, q):
    # Inputs and output retention are part of each complete timed invocation.
    return evaluate(method, *frozen_inputs(q))


def measure_child(method, mode, q):
    """Cold process; parent supplies wall timeout, this supplies CPU limit.

    Task allocation begins after common imports, before packed inputs/masks and
    lazy views. JSON serialization and tracemalloc bookkeeping are not mixed
    into its task-allocation peak. RSS and timing are untraced separate modes.
    """
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    if mode == "allocation":
        tracemalloc.start()
        result, stats = whole_operation(method, q)
        current, peak = tracemalloc.get_traced_memory()
        tracer_bytes = tracemalloc.get_tracemalloc_memory()
        tracemalloc.stop()
        metric = dict(current_bytes=current, peak_bytes=peak, tracer_bytes=tracer_bytes)
    elif mode == "rss":
        result, stats = whole_operation(method, q)
        metric = dict(maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    elif mode == "timing":
        before = time.perf_counter()
        result, stats = whole_operation(method, q)
        warmup = time.perf_counter()-before
        # >=10 ms per batch is well above the performance counter resolution.
        repeats = max(1, min(1000, int(0.01/max(warmup, 1e-9))+1))
        batches = []
        for _ in range(3):
            before = time.perf_counter()
            for _ in range(repeats):
                result, stats = whole_operation(method, q)
            batches.append(time.perf_counter()-before)
        per_call = [b/repeats for b in batches]
        metric = dict(repetitions=repeats, batch_seconds=batches,
                      median_seconds=statistics.median(per_call),
                      range_seconds=[min(per_call), max(per_call)],
                      clock_resolution=time.get_clock_info("perf_counter").resolution)
    else:
        raise ValueError("unknown metric")
    return dict(method=method, mode=mode, q=q, numerators=result, stats=stats,
                metric=metric, cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime
                + resource.getrusage(resource.RUSAGE_SELF).ru_stime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", choices=("stream", "table", "noop"), required=True)
    parser.add_argument("--mode", choices=("allocation", "rss", "timing"), required=True)
    parser.add_argument("--q", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(measure_child(args.child, args.mode, args.q)))
