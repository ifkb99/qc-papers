"""Exhaustive audit of the exact dyadic categorical-kernel primitives.

This bounded diagnostic targets the C60 sampling-kernel boundary, not the
floating prefix oracle.  It exhausts every L-bit random word for tiny rational
weight vectors and independently computes the ceil-CDF bin counts.  The
unbiasedness of a supplied random-bit source is an assumption; seeded RNG
frequencies are not used as proof.

P1: round_dyadic implements exact nearest, ties-to-even signed rounding for
    positive, negative, exact, and half-grid Fraction inputs.
P2: integer_cdf_index agrees word-by-word with the strict rational CDF rule,
    and integer_cdf_counts agrees with exhaustive enumeration over a bounded
    sweep; its TV error is at most (d-1)/2^L for tested laws.
P3: zero-weight/all-zero and invalid-negative cases have explicit semantics,
    with all-zero fallback tested on a divisible uniform word space; scripted
    uniform_integer consumes exactly the needed rejection/no-draw words.
C1: reducing random words modulo the number of bins is not a valid weighted
    categorical kernel.
C2: deleting a genuinely all-zero approximate block and renormalizing changes
    the global law; filtering exact zero-weight bins while preserving labels is
    tested as valid rather than treated as a failure.

The helper under test is imported only when the experiment runs.  No generic
propagator or production sampler is reimplemented here.  All arrays are tiny
and the exhaustive word count is capped before construction.
"""
from __future__ import annotations

import json
import itertools
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


MAX_WORDS = 1 << 8
MAX_BINS = 8
MAX_BYTES = 16 << 20


def report_path(prefix: str = "dyadic_kernel") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_words(random_bits: int, bins: int) -> int:
    if not isinstance(random_bits, int) or not 0 <= random_bits <= 8:
        raise ValueError("random_bits outside bounded exhaustive fixture")
    if not isinstance(bins, int) or not 1 <= bins <= MAX_BINS:
        raise ValueError("bin count outside bounded exhaustive fixture")
    words = 1 << random_bits
    if words * bins * 8 > MAX_BYTES:
        raise MemoryError("exhaustive categorical workspace exceeds cap")
    return words


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def expected_positive_counts(weights: list[Fraction], random_bits: int) -> tuple[int, ...]:
    words = guard_words(random_bits, len(weights))
    if any(weight < 0 for weight in weights) or not any(weights):
        raise ValueError("expected-count reference requires nonnegative nonzero weights")
    total = sum(weights, Fraction(0))
    scale = 1 << random_bits
    previous = 0
    counts = []
    cumulative = Fraction(0)
    for weight in weights:
        cumulative += weight
        boundary = min(words, ceil_fraction(Fraction(scale) * cumulative / total))
        counts.append(boundary - previous)
        previous = boundary
    return tuple(counts)


def expected_index(weights: list[Fraction], random_word: int,
                   random_bits: int) -> int:
    words = guard_words(random_bits, len(weights))
    if any(weight < 0 for weight in weights):
        raise ValueError("expected-index reference requires nonnegative weights")
    if not any(weights):
        weights = [Fraction(1)] * len(weights)
    total = sum(weights, Fraction(0))
    cutoff = Fraction(random_word * total, words)
    cumulative = Fraction(0)
    for index, weight in enumerate(weights):
        cumulative += weight
        if cutoff < cumulative:
            return index
    raise AssertionError("strict exact CDF reference did not cover word")


def integer_weights(weights: list[Fraction]) -> list[int]:
    """Clear exact denominators without changing any categorical ratios."""
    denominator = math.lcm(*(value.denominator for value in weights))
    return [int(value * denominator) for value in weights]


def as_strings(values):
    return [str(value) for value in values]


class ScriptedRng:
    """Deterministic API harness; not evidence about RNG unbiasedness."""
    def __init__(self, words):
        self.words = list(words)
        self.calls = []

    def getrandbits(self, bits):
        self.calls.append(bits)
        if not self.words:
            raise AssertionError("scripted RNG exhausted")
        return self.words.pop(0)


def main() -> None:
    # Import after predictions are declared and only after the parent helper is
    # available.  This keeps the standalone artifact honest during handoff.
    from lab.verified_prefix import (integer_cdf_counts, integer_cdf_index,
                                     round_dyadic, uniform_integer)

    exp = Experiment("dyadic_kernel", doc=__doc__)
    exp.predict("P1", "dyadic nearest-even signed rounding matches exact Fraction cases")
    exp.predict("P2", "strict-CDF index/count enumeration and TV bound agree word-by-word")
    exp.predict("P3", "zero/negative/fallback and exact uniform-integer edge cases behave explicitly")
    exp.must_fail("C1", "weighted modulo random-word mapping must disagree with exact CDF")
    exp.must_fail("C2", "deleting an all-zero approximate block must change the global law")
    started = time.perf_counter()
    rows: list[dict] = []

    exp.section("P1 signed nearest-even dyadic rounding")
    rounding_cases = [
        (Fraction(0), 2, 0),
        (Fraction(1, 8), 2, 0),       # +0.5 grid tie -> even 0
        (Fraction(3, 8), 2, 2),       # +1.5 grid tie -> even 2
        (Fraction(-1, 8), 2, 0),      # negative tie
        (Fraction(-3, 8), 2, -2),
        (Fraction(5, 16), 2, 1),      # non-tie
        (Fraction(-5, 16), 2, -1),
        (Fraction(7, 3), 0, 2),
        (Fraction(-7, 3), 0, -2),
    ]
    rounding_failures = []
    for value, bits, expected in rounding_cases:
        got = round_dyadic(value, bits)
        if got != expected:
            rounding_failures.append((str(value), bits, expected, got))
    exp.check("P1", not rounding_failures,
              f"cases={len(rounding_cases)}, failures={rounding_failures}")
    rows.append(dict(series="rounding", cases=len(rounding_cases),
                     failures=rounding_failures))

    exp.section("P2 exhaustive strict-CDF counts and TV bound")
    categorical_cases = [
        ([Fraction(1), Fraction(1), Fraction(2)], 4),
        ([Fraction(1), Fraction(1), Fraction(1)], 2),
        ([Fraction(2), Fraction(0), Fraction(3), Fraction(1)], 3),
        ([Fraction(1, 3), Fraction(2, 5), Fraction(7, 11)], 5),
    ]
    categorical_failures = []
    max_tv = Fraction(0)
    for weights, bits in categorical_cases:
        words = guard_words(bits, len(weights))
        expected = expected_positive_counts(weights, bits)
        integer = integer_weights(weights)
        got_counts = tuple(integer_cdf_counts(integer, bits))
        observed = [0] * len(weights)
        bad_words = []
        for word in range(words):
            index = integer_cdf_index(integer, word, bits)
            reference_index = expected_index(weights, word, bits)
            if not 0 <= index < len(weights):
                bad_words.append((word, index))
            elif index != reference_index:
                bad_words.append((word, index, reference_index))
            else:
                observed[index] += 1
        if tuple(observed) != expected or got_counts != expected or bad_words:
            categorical_failures.append(dict(weights=as_strings(weights), bits=bits,
                                              expected=expected, observed=observed,
                                              helper_counts=got_counts,
                                              bad_words=bad_words))
        total = sum(weights, Fraction(0))
        ideal = [weight / total for weight in weights]
        empirical = [Fraction(count, words) for count in observed]
        tv = sum(abs(a - b) for a, b in zip(empirical, ideal)) / 2
        max_tv = max(max_tv, tv)
        bound = Fraction(len(weights) - 1, words)
        if tv > bound:
            categorical_failures.append(dict(weights=as_strings(weights), bits=bits,
                                              tv=str(tv), bound=str(bound),
                                              kind="TV-bound"))
        rows.append(dict(series="categorical", weights=as_strings(weights), bits=bits,
                         expected_counts=list(expected), observed_counts=observed,
                         helper_counts=list(got_counts), tv=str(tv),
                         tv_bound=str(bound)))
    # Bounded exhaustive sweep: every integer weight vector in {0,1,2,3}^d,
    # d=2..4, every L=0..4, and every word in the finite random-word space.
    sweep_cases = 0
    sweep_words = 0
    sweep_failures = []
    for bins in range(2, 5):
        for raw in itertools.product(range(4), repeat=bins):
            weights = [Fraction(value) for value in raw]
            integer = list(raw)
            for bits in range(5):
                words = guard_words(bits, bins)
                reference_weights = weights if any(weights) else [Fraction(1)] * bins
                expected_counts = expected_positive_counts(reference_weights, bits)
                helper_counts = tuple(integer_cdf_counts(integer, bits))
                observed = [0] * bins
                for word in range(words):
                    got = integer_cdf_index(integer, word, bits)
                    want = expected_index(weights, word, bits)
                    sweep_words += 1
                    if got != want:
                        sweep_failures.append((raw, bits, word, got, want))
                    observed[got] += 1
                sweep_cases += 1
                if tuple(observed) != expected_counts or helper_counts != expected_counts:
                    sweep_failures.append((raw, bits, "counts",
                                           tuple(observed), helper_counts, expected_counts))
    rows.append(dict(series="bounded_exhaustive_sweep", bins="2..4", weights="0..3",
                     random_bits="0..4", cases=sweep_cases, words=sweep_words,
                     failures=sweep_failures[:20], failure_count=len(sweep_failures)))
    exp.check("P2", not categorical_failures and not sweep_failures,
              f"cases={len(categorical_cases)}+{sweep_cases}, max_tv={max_tv}, "
              f"word_count={sweep_words}, failures={categorical_failures}")

    exp.section("P3 zero, negative, and fallback semantics")
    zero_weights = [0, 0, 0, 0]
    zero_bits = 4
    zero_words = guard_words(zero_bits, len(zero_weights))
    zero_counts = [0] * len(zero_weights)
    zero_invalid = []
    for word in range(zero_words):
        index = integer_cdf_index(zero_weights, word, zero_bits)
        if not 0 <= index < len(zero_weights):
            zero_invalid.append((word, index))
        else:
            zero_counts[index] += 1
    rejected = []
    for call in (lambda: integer_cdf_index([1, -1], 0, 2),
                 lambda: integer_cdf_counts([1, -1], 2)):
        try:
            call()
            rejected.append(False)
        except (ValueError, ArithmeticError):
            rejected.append(True)
    negative_rejected = all(rejected)
    exp.check("P3", not zero_invalid and zero_counts == [4, 4, 4, 4]
              and negative_rejected,
              f"all_zero_counts={zero_counts}, invalid={zero_invalid}, "
              f"negative_rejected={negative_rejected}")
    filtered_zero_weights = [0, 3, 1]
    filtered_counts = tuple(integer_cdf_counts(filtered_zero_weights, 4))
    filtered_expected = expected_positive_counts(
        [Fraction(value) for value in filtered_zero_weights], 4)
    scripted_reject = ScriptedRng([3, 2])
    draw_reject = uniform_integer(3, scripted_reject)
    scripted_accept = ScriptedRng([2])
    draw_accept = uniform_integer(3, scripted_accept)
    scripted_one = ScriptedRng([])
    draw_one = uniform_integer(1, scripted_one)
    uniform_ok = (draw_reject == 2 and scripted_reject.calls == [2, 2]
                  and draw_accept == 2 and scripted_accept.calls == [2]
                  and draw_one == 0 and scripted_one.calls == [])
    exp.check("P3", filtered_counts == filtered_expected and uniform_ok,
              f"filtered_zero_counts={filtered_counts}, "
              f"uniform_scripted={uniform_ok}")
    rows.append(dict(series="edge_cases", all_zero_counts=zero_counts,
                     negative_rejected=negative_rejected,
                     filtered_zero_counts=list(filtered_counts),
                     filtered_zero_expected=list(filtered_expected),
                     uniform_scripted=uniform_ok,
                     uniform_calls=[scripted_reject.calls,
                                   scripted_accept.calls, scripted_one.calls]))

    exp.section("must-fail biased controls")
    weights = [Fraction(1), Fraction(3)]
    bits = 4
    exact = expected_positive_counts(weights, bits)
    integer = integer_weights(weights)
    modulo = [0, 0]
    for word in range(1 << bits):
        modulo[word % len(weights)] += 1
    exp.fail_check("C1", tuple(modulo) != exact,
                   f"weighted modulo counts={modulo}, exact={exact}")

    # A genuinely all-zero approximate block cannot be silently removed: its
    # ideal global mass remains part of normalization.  Filtering an exact
    # zero coordinate inside a block is tested above as valid behavior.
    ideal_block_masses = [Fraction(3), Fraction(1)]
    deleted_block_masses = [Fraction(0), Fraction(1)]
    ideal_block_law = [x / sum(ideal_block_masses) for x in ideal_block_masses]
    deleted_block_law = [x / sum(deleted_block_masses) for x in deleted_block_masses]
    deletion_tv = sum(abs(a - b) for a, b in zip(ideal_block_law,
                                                  deleted_block_law)) / 2
    exp.fail_check("C2", deletion_tv > 0,
                   f"ideal block law={as_strings(ideal_block_law)}, "
                   f"deleted-block law={as_strings(deleted_block_law)}, "
                   f"TV={deletion_tv}")
    rows.append(dict(series="must_fail", weighted_modulo=modulo,
                     weighted_exact=list(exact),
                     deleted_all_zero_block_tv=str(deletion_tv)))

    report = report_path()
    exp.finish(report_path=report, rows=rows,
               metadata=dict(max_words=MAX_WORDS, max_bins=MAX_BINS,
                             max_dense_bytes=MAX_BYTES,
                             rng_unbiasedness="assumed, not empirically certified",
                             cdf_rule="strict word*total < cumulative*2^L",
                             elapsed_seconds=time.perf_counter() - started))
    print(f"report: {report}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("dyadic_kernel_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
