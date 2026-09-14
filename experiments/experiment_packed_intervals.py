"""Compact C89 interval-table comparator for the C94 conditional evaluator.

Each record stores high endpoint and signed shift; low is the previous high.
Two composition buffers and a separate packed final-image validation buffer
are charged. Shared native macro construction and range DP remain unchanged.
The versioned task harness supplies pre-execution predictions, native controls
and isolated measurements. No strongest-implementation or general CNOT claim.
"""
from __future__ import annotations

import argparse
import json
import resource
import statistics
import time
import tracemalloc

from experiments.experiment_controlled_intervals import macro_intervals
from experiments.experiment_interval_stream import (
    contract_pieces, evaluate as existing_evaluate, frozen_inputs, validate,
)


class PackedIntervals:
    """Flat exact records, without an O(q) Python-int or tuple container."""

    __slots__ = ("word_bytes", "size", "data", "count", "last_hi", "last_shift")

    def __init__(self, width):
        self.word_bytes = (width+2+7)//8
        self.size = 1 << (width+1)
        self.data = bytearray()
        self.count = 0
        self.last_hi = 0
        self.last_shift = None

    def append(self, lo, hi, shift):
        if lo != self.last_hi or not lo < hi <= self.size:
            raise AssertionError("invalid packed input partition")
        if not 0 <= lo+shift < hi+shift <= self.size:
            raise AssertionError("invalid packed image interval")
        w = self.word_bytes
        if self.count and shift == self.last_shift:
            pos = (self.count-1)*2*w
            self.data[pos:pos+w] = hi.to_bytes(w, "little")
        else:
            self.data.extend(hi.to_bytes(w, "little"))
            self.data.extend(shift.to_bytes(w, "little", signed=True))
            self.count += 1
        self.last_hi, self.last_shift = hi, shift

    def __iter__(self):
        w, lo = self.word_bytes, 0
        view = memoryview(self.data)
        for offset in range(0, len(view), 2*w):
            hi = int.from_bytes(view[offset:offset+w], "little")
            shift = int.from_bytes(view[offset+w:offset+2*w], "little", signed=True)
            yield lo, hi, shift
            lo = hi


def validate_packed_images(table):
    """Validate the image tiling using a packed in-place heapsort.

    A transient (image_lo,image_hi) byte buffer replaces C89's sorted tuple
    list. Heapsort uses constant many records, O(P log P) comparisons, and no
    recursive stack or Python-int array. The original domain order is kept.
    """
    if table.last_hi != table.size:
        raise AssertionError("incomplete packed domain")
    w, data = table.word_bytes, bytearray()
    for lo, hi, shift in table:
        data.extend((lo+shift).to_bytes(w, "little"))
        data.extend((hi+shift).to_bytes(w, "little"))
    stride = 2*w
    view = memoryview(data)

    def left(index):
        offset = stride*index
        return int.from_bytes(view[offset:offset+w], "little")

    def swap(a, b):
        a, b = a*stride, b*stride
        temporary = bytes(view[a:a+stride])
        view[a:a+stride] = view[b:b+stride]
        view[b:b+stride] = temporary

    def sift(start, stop):
        root = start
        while 2*root+1 < stop:
            child = 2*root+1
            if child+1 < stop and left(child) < left(child+1):
                child += 1
            if left(root) >= left(child):
                break
            swap(root, child)
            root = child

    for start in range(table.count//2-1, -1, -1):
        sift(start, table.count)
    for stop in range(table.count-1, 0, -1):
        swap(0, stop)
        sift(0, stop)
    cursor = 0
    for offset in range(0, len(view), stride):
        lo = int.from_bytes(view[offset:offset+w], "little")
        hi = int.from_bytes(view[offset+w:offset+stride], "little")
        if lo != cursor or not lo < hi <= table.size:
            raise AssertionError("packed images overlap or have gaps")
        cursor = hi
    if cursor != table.size:
        raise AssertionError("incomplete packed image")


def compile_packed(width, modulus, constants, t, h, controls):
    validate(width, modulus, constants, t, h, controls)
    current = PackedIntervals(width)
    current.append(0, 1 << (width+1), 0)
    peak_composition_payload = len(current.data)
    for step, (constant, control) in enumerate(zip(constants, controls), 1):
        local = macro_intervals(width, modulus, constant, t, h, control)
        after = PackedIntervals(width)
        for lo, hi, shift in current:
            image_lo, image_hi = lo+shift, hi+shift
            for after_lo, after_hi, after_shift in local:
                left, right = max(image_lo, after_lo), min(image_hi, after_hi)
                if left < right:
                    after.append(left-shift, right-shift, shift+after_shift)
        if after.count > 1+15*step:
            raise AssertionError("packed composition count bound violated")
        peak_composition_payload = max(peak_composition_payload,
                                       len(current.data)+len(after.data))
        current = after
        # Do not retain the previous destination through the next allocation.
        del after
    validate_packed_images(current)
    return current, peak_composition_payload


def evaluate(method, width, modulus, constants, t, h, controls, queries):
    if method != "packed":
        return existing_evaluate(method, width, modulus, constants, t, h, controls, queries)
    validate(width, modulus, constants, t, h, controls, queries)
    pieces, peak_payload = compile_packed(width, modulus, constants, t, h, controls)
    values = contract_pieces(width, pieces, queries)
    return values, dict(emitted_pieces=pieces.count, record_bytes=2*pieces.word_bytes,
                        final_payload_bytes=len(pieces.data),
                        peak_composition_payload_bytes=peak_payload,
                        validation_payload_bytes=len(pieces.data))


def whole_operation(method, q):
    return evaluate(method, *frozen_inputs(q))


def measure_child(method, mode, q):
    # Same trace boundaries and untraced timing protocol as the frozen C94
    # campaign. The unchanged C94 evaluator handles stream/noop invocations.
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
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return dict(method=method, mode=mode, q=q, numerators=result, stats=stats,
                metric=metric, cpu_seconds=usage.ru_utime+usage.ru_stime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child", choices=("packed", "stream", "noop"), required=True)
    parser.add_argument("--mode", choices=("allocation", "rss", "timing"), required=True)
    parser.add_argument("--q", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(measure_child(args.child, args.mode, args.q)))
