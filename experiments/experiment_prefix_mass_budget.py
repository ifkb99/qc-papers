"""Exact binary-tree test of the unnormalized prefix-mass error budget.

For a depth-t binary law, let the exact child weights at every prefix be
the *unnormalized* ideal prefix masses.  Replace each nonnegative weight by a
nonnegative dyadic value with absolute error at most eta=2^-p, normalize the
two children at that prefix, and use an exact L-bit inverse-CDF kernel.  An
all-zero pair uses the declared uniform fallback.  The proposed law-level
budget is

    TV <= 2*(2^t-1)*eta + t*2^-L.

P1 tests this bound on one fixed t=5 law while only p varies.  The law has
exactly-zero leaves and a positive 2^-128 leaf, so all-zero local prefixes
and cancellation/rare-mass behavior are exercised.  C1 records that a
conditional TV error can be 1/2 on a tiny positive prefix even
though its absolute weighted deficit is tiny.  C2 records why an uncharged
global rescaling before absolute rounding is invalid: scaling a two-bin law
by 2^-128 at p=20 makes both approximate weights zero and invokes uniform
fallback, giving TV 1/4.

All scientific arithmetic is exact Fraction or integer arithmetic.  The
random-word space is enumerated only at fixed L=12; no generic propagator or
floating approximation is used.  Reports state bounded operation counts,
not process RSS.
"""
from __future__ import annotations

import json
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

from lab import Experiment


TREE_DEPTH = 5
FIXED_RANDOM_BITS = 12
SWEEP_MAX_P = 12
CONTROL_P = 20
TINY_EXPONENT = 128
MAX_BYTES = 16 << 20


def report_path(prefix: str = "prefix_mass_budget") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_tree(depth: int, random_bits: int) -> tuple[int, int]:
    if type(depth) is not int or not 0 <= depth <= TREE_DEPTH:
        raise ValueError("depth outside bounded binary-tree fixture")
    if type(random_bits) is not int or not 0 <= random_bits <= FIXED_RANDOM_BITS:
        raise ValueError("random-bit width outside bounded fixture")
    leaves = 1 << depth
    words = 1 << random_bits
    # Words are ITERATED, not stored. Guard the retained tiny level/frontier
    # entries with a conservative allowance for this fixed 128-bit fixture.
    # This is neither a NumPy payload formula nor a measured RSS claim.
    if (4*leaves+16) * 1024 > MAX_BYTES:
        raise MemoryError("prefix-mass diagnostic exceeds allocation cap")
    return leaves, words


def floor_dyadic(value: Fraction, bits: int) -> Fraction:
    if not isinstance(value, Fraction) or value < 0 or type(bits) is not int:
        raise ValueError("floor_dyadic requires a nonnegative Fraction")
    if not 0 <= bits <= CONTROL_P:
        raise ValueError("dyadic precision outside bounded fixture")
    scaled = value * (1 << bits)
    return Fraction(scaled.numerator // scaled.denominator, 1 << bits)


def strict_cdf_counts(integer_weights: tuple[int, int], random_bits: int) -> tuple[int, int]:
    """Enumerate the exact strict integer CDF at every L-bit word."""
    if len(integer_weights) != 2 or any(type(w) is not int or w < 0
                                        for w in integer_weights):
        raise ValueError("require two nonnegative integer weights")
    _, words = guard_tree(1, random_bits)
    weights = integer_weights if any(integer_weights) else (1, 1)
    total = sum(weights)
    counts = [0, 0]
    for word in range(words):
        cutoff = word * total
        if cutoff < weights[0] * words:
            counts[0] += 1
        else:
            counts[1] += 1
    return tuple(counts)


def prefix_levels(leaf_weights: list[Fraction], depth: int) -> list[list[Fraction]]:
    guard_tree(depth, FIXED_RANDOM_BITS)
    if len(leaf_weights) != 1 << depth or any(w < 0 for w in leaf_weights):
        raise ValueError("invalid exact leaf law")
    levels: list[list[Fraction]] = []
    for prefix_depth in range(depth + 1):
        width = 1 << prefix_depth
        block = 1 << (depth - prefix_depth)
        levels.append([sum(leaf_weights[start:start + block], Fraction(0))
                       for start in range(0, len(leaf_weights), block)])
        if len(levels[-1]) != width:
            raise AssertionError("prefix level shape mismatch")
    return levels


def approximate_law(levels: list[list[Fraction]], precision: int,
                    random_bits: int) -> tuple[list[Fraction], dict]:
    """Enumerate the complete normalized approximate tree law."""
    depth = len(levels) - 1
    guard_tree(depth, random_bits)
    eta = Fraction(1, 1 << precision)
    current = {0: Fraction(1)}
    all_zero_nodes = 0
    local_failures = []
    for level in range(depth):
        child_mass: dict[int, Fraction] = {}
        for prefix, parent_probability in current.items():
            true_children = (levels[level + 1][2 * prefix],
                             levels[level + 1][2 * prefix + 1])
            approximate = tuple(floor_dyadic(value, precision)
                                for value in true_children)
            if any(value - rounded < 0 or value - rounded >= eta
                   for value, rounded in zip(true_children, approximate)):
                local_failures.append((level, prefix, true_children, approximate))
            integer = tuple(int(value * (1 << precision))
                            for value in approximate)
            if not any(integer):
                all_zero_nodes += 1
            counts = strict_cdf_counts(integer, random_bits)
            if sum(counts) != 1 << random_bits:
                local_failures.append((level, prefix, "cdf", counts))
            for branch, count in enumerate(counts):
                child = 2 * prefix + branch
                child_mass[child] = child_mass.get(child, Fraction(0)) \
                    + parent_probability * Fraction(count, 1 << random_bits)
        if sum(child_mass.values(), Fraction(0)) != 1:
            local_failures.append((level, "normalization", sum(child_mass.values())))
        current = child_mass
    leaves = [current.get(index, Fraction(0)) for index in range(1 << depth)]
    return leaves, dict(all_zero_nodes=all_zero_nodes,
                        local_failures=local_failures,
                        final_mass=sum(leaves, Fraction(0)))


def total_variation(left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum(abs(a - b) for a, b in zip(left, right)) / 2


def fixture_leaves() -> list[Fraction]:
    tiny = Fraction(1, 1 << TINY_EXPONENT)
    zero_indices = {3, 7, 16, 17, 24}
    tiny_index = 31
    ordinary = [index for index in range(1 << TREE_DEPTH)
                if index not in zero_indices and index != tiny_index]
    ordinary_mass = (Fraction(1) - tiny) / len(ordinary)
    leaves = [Fraction(0) for _ in range(1 << TREE_DEPTH)]
    for index in ordinary:
        leaves[index] = ordinary_mass
    leaves[tiny_index] = tiny
    if sum(leaves, Fraction(0)) != 1:
        raise AssertionError("fixture leaf law is not normalized")
    return leaves


def main() -> None:
    exp = Experiment("prefix_mass_budget", doc=__doc__)
    exp.predict("P1", "exact normalized binary-tree law obeys the absolute mass budget")
    exp.must_fail("C1", "tiny weighted mass does not imply a small conditional relative error")
    exp.must_fail("C2", "uncharged pre-rounding rescaling can violate the original budget")
    started = time.perf_counter()
    rows: list[dict] = []

    exp.section("P1 fixed depth-5 law; sweep only dyadic precision p")
    guard_tree(TREE_DEPTH, FIXED_RANDOM_BITS)
    ideal = fixture_leaves()
    levels = prefix_levels(ideal, TREE_DEPTH)
    assert any(mass == 0 for level in levels[:-1] for mass in level)
    from lab.sampling_error import prefix_mass_error_budget
    failures = []
    worst_ratio = Fraction(0)
    for precision in range(SWEEP_MAX_P + 1):
        approximate, details = approximate_law(levels, precision, FIXED_RANDOM_BITS)
        eta = Fraction(1, 1 << precision)
        tv = total_variation(approximate, ideal)
        bound = 2 * ((1 << TREE_DEPTH) - 1) * eta \
            + TREE_DEPTH * Fraction(1, 1 << FIXED_RANDOM_BITS)
        ratio = tv / bound if bound else Fraction(0)
        worst_ratio = max(worst_ratio, ratio)
        implemented_bound = prefix_mass_error_budget(TREE_DEPTH,eta,FIXED_RANDOM_BITS)["total_tv_upper_bound"]
        if (details["local_failures"] or details["final_mass"] != 1 or tv > bound
                or implemented_bound != min(Fraction(1),bound)):
            failures.append(dict(precision=precision, tv=str(tv), bound=str(bound),
                                 details=details))
        row = dict(series="precision_sweep", precision=precision, eta=str(eta),
                   tv=str(tv), bound=str(bound), ratio=str(ratio),
                   all_zero_nodes=details["all_zero_nodes"],
                   approximate_mass=str(details["final_mass"]),
                   local_failure_count=len(details["local_failures"]))
        rows.append(row)
    exp.check("P1", not failures,
              f"depth={TREE_DEPTH}, leaves={len(ideal)}, L={FIXED_RANDOM_BITS}, "
              f"p=0..{SWEEP_MAX_P}, worst_tv_to_bound={worst_ratio}, "
              f"failures={failures}")

    exp.section("must-fail relative-error and rescaling controls")
    tiny = Fraction(1, 1 << TINY_EXPONENT)
    tiny_rounded = floor_dyadic(tiny, CONTROL_P)
    relative_error = (tiny - tiny_rounded) / tiny
    weighted_deficit = tiny - tiny_rounded
    tiny_counts = strict_cdf_counts((0,0),FIXED_RANDOM_BITS)
    tiny_kernel = [Fraction(count,1 << FIXED_RANDOM_BITS) for count in tiny_counts]
    conditional_tv = total_variation(tiny_kernel,[Fraction(1),Fraction(0)])
    exp.fail_check("C1", conditional_tv == Fraction(1,2)
                   and tiny*conditional_tv <= weighted_deficit,
                   f"tiny_mass={tiny}, rounded={tiny_rounded}, "
                   f"conditional_TV={conditional_tv}, weighted_deficit={weighted_deficit}")

    scale = Fraction(1, 1 << TINY_EXPONENT)
    root_ideal = [Fraction(3, 4), Fraction(1, 4)]
    scaled = [scale * value for value in root_ideal]
    rounded = [floor_dyadic(value, CONTROL_P) for value in scaled]
    counts = strict_cdf_counts(tuple(int(value * (1 << CONTROL_P)) for value in rounded),
                               FIXED_RANDOM_BITS)
    rescaled_law = [Fraction(count, 1 << FIXED_RANDOM_BITS) for count in counts]
    rescaled_tv = total_variation(rescaled_law, root_ideal)
    rescaled_bound = 2 * (2 - 1) * Fraction(1, 1 << CONTROL_P) \
        + Fraction(1, 1 << FIXED_RANDOM_BITS)
    exp.fail_check("C2", rescaled_tv > rescaled_bound,
                   f"ideal={root_ideal}, scale={scale}, rounded={rounded}, "
                   f"fallback_counts={counts}, law={rescaled_law}, "
                   f"TV={rescaled_tv}, original_bound={rescaled_bound}")
    rows.append(dict(series="controls", tiny_mass=str(tiny),
                     tiny_relative_error=str(relative_error),
                     tiny_conditional_tv=str(conditional_tv),
                     tiny_weighted_deficit=str(weighted_deficit),
                     rescaling_ideal=[str(x) for x in root_ideal],
                     rescaling_factor=str(scale),
                     rescaling_rounded=[str(x) for x in rounded],
                     rescaling_law=[str(x) for x in rescaled_law],
                     rescaling_tv=str(rescaled_tv),
                     original_bound=str(rescaled_bound)))

    report = report_path()
    exp.finish(report_path=report, rows=rows,
               metadata=dict(depth=TREE_DEPTH, fixed_random_bits=FIXED_RANDOM_BITS,
                             precision_sweep=f"0..{SWEEP_MAX_P}",
                             tiny_exponent=TINY_EXPONENT,
                             arithmetic="exact Fraction and integer arithmetic only",
                             allocation_bound="bounded retained-entry allowance, not word-array allocation/RSS",
                             elapsed_seconds=time.perf_counter() - started))
    print(f"report: {report}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("prefix_mass_budget_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
