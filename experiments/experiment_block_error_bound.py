"""Exact-rational audit of the TODO24 block-kernel error bound.

P_s and R_s below are FULL basis-state Born laws.  Their block marginals are
only used to define the within-block kernels.  The experiment checks the
derived inequality

    sum_B P(B) TV(p(.|B), r(.|B)) <= 2 TV(P, R),

including an explicitly normalized fallback when R(B)=0.  It also checks the
one-step Markov-kernel triangle bound and the exact rational budgets in
lab.sampling_error.  No circuit propagator or floating-point amplitude is
used.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The full-law weighted conditional inequality holds exactly for all
      deterministic rational partition cases, including zero approximate
      blocks.
  P2  The one-step approximate kernel bound and helper budgets hold exactly
      under their stated oracle/random-bit assumptions.
  C1  Replacing full-law TV by block-mass TV must fail: equal block masses can
      hide completely different within-block laws.
  C2  Omitting the R(B)=0 fallback contribution must fail on a rare block.

This is a bounded theorem/constant audit, not a numerical sampler certificate.
"""
from __future__ import annotations

import random
from fractions import Fraction
from datetime import datetime, timezone
from pathlib import Path

from lab import Experiment
from lab.sampling_error import plan_prefix_accuracy, prefix_error_budget, rejection_error_budget


def tv(a, b):
    if len(a) != len(b) or sum(a) != 1 or sum(b) != 1 or min(a) < 0 or min(b) < 0:
        raise ValueError("TV requires matching normalized nonnegative full laws")
    return sum(abs(x-y) for x, y in zip(a, b)) / 2


def block_masses(law, blocks):
    return tuple(sum(law[i] for i in block) for block in blocks)


def conditional_bound(P, R, blocks):
    """Return weighted ideal conditional TV and the exact 2TV bound."""
    lhs = Fraction(0)
    for block in blocks:
        p_mass = sum(P[i] for i in block)
        r_mass = sum(R[i] for i in block)
        if p_mass == 0:
            continue
        if r_mass == 0:
            # Measure the actual uniform fallback, not its worst-case bound.
            lhs += sum(abs(P[i] - p_mass/Fraction(len(block))) for i in block) / 2
            continue
        lhs += sum(abs(P[i] - p_mass*R[i]/r_mass) for i in block) / 2
    full_tv = tv(P, R)
    return lhs, 2*full_tv


def apply_kernel(Q, R, blocks):
    """Apply R's normalized within-block law while retaining Q block mass."""
    result = [Fraction(0)] * len(Q)
    for block in blocks:
        q_mass = sum(Q[i] for i in block)
        r_mass = sum(R[i] for i in block)
        if r_mass:
            for i in block:
                result[i] = q_mass * R[i] / r_mass
        else:
            share = q_mass / len(block)
            for i in block:
                result[i] = share
    return result


def random_law(rng, size, total=120):
    weights = [rng.randrange(total + 1) for _ in range(size)]
    if not any(weights):
        weights[rng.randrange(size)] = 1
    denominator = sum(weights)
    return [Fraction(x, denominator) for x in weights]


def random_partition(rng, size, block_count):
    cuts = sorted(rng.sample(range(1, size), block_count - 1))
    edges = [0] + cuts + [size]
    return [tuple(range(edges[j], edges[j+1])) for j in range(block_count)]


def main():
    exp = Experiment("block_error_bound", doc=__doc__)
    exp.predict("P1", "full-law weighted conditional TV is at most twice full-law TV")
    exp.predict("P2", "kernel triangle and rational budget constants hold")
    exp.must_fail("C1", "block-mass-only TV falsely certifies equal-mass opposite conditionals")
    exp.must_fail("C2", "zero approximate block is discarded without fallback mass")

    # Exact rational randomized partitions: deterministic seed, no floating
    # arithmetic, and no amplitude/state arrays.
    rng = random.Random(20260924)
    max_slack = Fraction(0)
    chain_slack = Fraction(0)
    cases = 0
    for size in range(2, 13):
        for block_count in range(1, min(size, 5) + 1):
            for _ in range(8):
                blocks = random_partition(rng, size, block_count)
                P = random_law(rng, size)
                R = random_law(rng, size)
                lhs, bound = conditional_bound(P, R, blocks)
                if lhs > bound:
                    raise AssertionError(("conditional inequality", size, blocks, P, R, lhs, bound))
                max_slack = max(max_slack, bound-lhs)

                # Construct Pprev and Pnext with identical block masses, then
                # perturb the prior law before applying R's kernel.
                Pprev = random_law(rng, size)
                Pnext = list(Pprev)
                for block in blocks:
                    mass = sum(Pprev[i] for i in block)
                    local = random_law(rng, len(block))
                    for i, weight in zip(block, local):
                        Pnext[i] = mass * weight
                Qprev = [(Pprev[i] + P[i]) / 2 for i in range(size)]
                Qnext = apply_kernel(Qprev, R, blocks)
                rhs = tv(Qprev, Pprev) + 2*tv(Pnext, R)
                actual = tv(Qnext, Pnext)
                if actual > rhs:
                    raise AssertionError(("chain inequality", size, blocks, actual, rhs))
                chain_slack = max(chain_slack, rhs-actual)
                cases += 1
    exp.check("P1", cases == sum(8*min(size, 5) for size in range(2, 13))
              and max_slack >= 0,
              f"exact rational partitions={cases}, max_bound_slack={max_slack}")
    exp.check("P2", chain_slack >= 0,
              f"exact one-step triangle cases={cases}, max_chain_slack={chain_slack}")

    # The clarified negative control: same block masses but full-law TV one.
    blocks = [(0, 1), (2, 3)]
    P = [Fraction(1, 2), 0, Fraction(1, 2), 0]
    R = [0, Fraction(1, 2), 0, Fraction(1, 2)]
    block_tv = tv(block_masses(P, blocks), block_masses(R, blocks))
    full_tv = tv(P, R)
    weighted, full_bound = conditional_bound(P, R, blocks)
    exp.fail_check("C1", weighted > 2*block_tv and full_tv == 1
                   and weighted <= full_bound,
                   f"block_mass_TV={block_tv}, full_TV={full_tv}, "
                   f"weighted_conditional_TV={weighted}, full_bound={full_bound}")

    # A zero approximate block: omitting its fallback term would claim zero,
    # while the exact rare-block contribution is positive and charged by full TV.
    eps = Fraction(1, 100)
    P0 = [eps, 0, 1-eps, 0]
    R0 = [0, 0, 1, 0]
    weighted0, bound0 = conditional_bound(P0, R0, blocks)
    omitted_zero_term = Fraction(0)
    exp.fail_check("C2", weighted0 > omitted_zero_term
                   and weighted0 <= bound0
                   and tv(P0, R0) == eps,
                   f"zero_block_weight={weighted0}, omitted={omitted_zero_term}, "
                   f"full_TV={tv(P0, R0)}, bound={bound0}")

    budget = prefix_error_budget(4, 2, [Fraction(1, 1024)], random_bits=4)
    planned = plan_prefix_accuracy(4, 2, 3, Fraction(1, 10))
    exp.check("P2", budget["scaled_coordinate_to_state_l2_factor"] == 8
              and budget["oracle_tv_upper_bound"] == Fraction(1, 64)
              and budget["finite_random_bits_tv_upper_bound"] == Fraction(1, 16)
              and budget["total_tv_upper_bound"] == Fraction(5, 64)
              and planned["total_tv_upper_bound"] <= Fraction(1, 10),
              f"budget={budget}, planner_total={planned['total_tv_upper_bound']}")
    wide = plan_prefix_accuracy(63, 2, 135, Fraction(1,1_000_000))
    for draw in range(50):
        P = random_law(rng, 8)
        approximate_law = random_law(rng, 8)
        mass = Fraction(draw, 50)
        C = Fraction(8)
        accepted_l1 = sum(abs(mass*q-p/C) for p,q in zip(P, approximate_law))
        if tv(P,approximate_law) > rejection_error_budget(C,accepted_l1):
            raise AssertionError("normalized accepted-law error exceeds proved bound")
    exp.check("P2", wide["total_tv_upper_bound"] <= Fraction(1,1_000_000),
              f"50 exact rejection-measure cases; wide oracle accuracy bits="
              f"{wide['absolute_accuracy_bits']}, random bits={wide['random_bits']}")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"block_error_bound_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    exp.finish(report_path=path, metadata={
        "exact_rational_cases": cases, "max_dimension": 12,
        "zero_block_actual_weighted_tv": str(weighted0),
        "zero_block_full_tv": str(tv(P0, R0)),
        "arithmetic": "exact Fraction; seeded finite partition audit",
        "full_law_not_block_marginal": True,
        "requires_uniform_oracle_contract": True,
        "wide_accuracy_plan": {k:str(v) if isinstance(v,Fraction) else v
                               for k,v in wide.items()},
        "rejection_measure_cases": 50,
    })
    print(f"report: {path}")


if __name__ == "__main__":
    main()
