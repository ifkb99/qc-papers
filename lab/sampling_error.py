"""Exact-rational error budgets, NOT a certificate for complex128 execution.

The prefix budget is conditional on a UNIFORM error bound for every real and
imaginary coordinate of the sqrt(M)-scaled post-prefix amplitude function.
Input representation, trig evaluation and contraction rounding must all be
included in that bound. Observed errors, output quantization bits, and a
machine's mantissa size do not establish the required oracle contract.

The kernel uses exact normalized squared weights, or a specified normalized
within-block fallback at an all-zero approximate block. An optional bound
covers exact inverse-CDF comparison to uniform finite random bits. It does
not automatically cover numpy.random plus floating-point division/CDF sums.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt


def _dimensions(width, block_size):
    if (type(width) is not int or not 0 <= width <= 63
            or type(block_size) is not int or not 1 <= block_size <= 64):
        raise ValueError("require integer 0<=width<=63 and 1<=block_size<=64")


def _rational(value, *, positive=False):
    if type(value) not in (int, Fraction):
        raise ValueError("error/accuracy bounds must be exact int or Fraction, not floats")
    result = Fraction(value)
    if result < 0 or positive and result == 0:
        raise ValueError("invalid nonnegative error or positive accuracy")
    return result


def prefix_error_budget(width, block_size, coordinate_errors, *, random_bits=None):
    """Conservative TV budget under an explicit scaled-amplitude oracle promise.

    One coordinate_errors entry per STOCHASTIC BLOCK UPDATE, not per amplitude
    query. If every real and imaginary scaled coordinate has error <=eta_s,
    global state L2 error is <=sqrt(2*b*2^t)*eta_s: M cancels. The block-kernel
    induction gives TV <=2*sum_s eps_s, without a lower block-mass bound.
    Ceiling the square root gives an exactly representable upper bound.

    Initial labels are sampled exactly. With L random bits and exact CDF
    boundaries, at most (d-1)/2^L TV is added per update on d<=max(b,2) bins.
    No arrays scale with period, exponent strings, or output dimension.
    """
    _dimensions(width, block_size)
    if not isinstance(coordinate_errors, (tuple, list)) or len(coordinate_errors) > 256:
        raise ValueError("supply at most 256 explicit per-update error bounds")
    errors = tuple(_rational(v) for v in coordinate_errors)
    if random_bits is not None and (type(random_bits) is not int or not 1 <= random_bits <= 4096):
        raise ValueError("random_bits must be None (ideal) or an integer in [1,4096]")
    dimension_factor = 2 * block_size * (1 << width)
    root = isqrt(dimension_factor)
    root += root * root < dimension_factor
    oracle_tv = 2 * root * sum(errors, Fraction(0))
    random_tv = (Fraction(len(errors) * (max(block_size, 2) - 1), 1 << random_bits)
                 if random_bits is not None else Fraction(0))
    return dict(oracle_tv_upper_bound=min(Fraction(1), oracle_tv),
                finite_random_bits_tv_upper_bound=min(Fraction(1), random_tv),
                total_tv_upper_bound=min(Fraction(1), oracle_tv + random_tv),
                scaled_coordinate_to_state_l2_factor=root,
                updates=len(errors), width=width, block_size=block_size,
                random_bits=random_bits, arithmetic="exact rational budget",
                requires_uniform_oracle_error_bound=True,
                certifies_existing_float_sampler=False)


def plan_prefix_accuracy(width, block_size, updates, target_tv):
    """Allocate half the target to oracle error and half to ideal finite bits.

    Returns an ABSOLUTE coordinate accuracy request 2^-p, not a sufficient
    machine mantissa. Producing this oracle with verified input/trig/contraction
    error remains a separate implementation task.
    """
    _dimensions(width, block_size)
    if type(updates) is not int or not 1 <= updates <= 256:
        raise ValueError("require 1..256 block updates")
    target = _rational(target_tv, positive=True)
    if target >= 1:
        raise ValueError("target TV must be strictly below one")
    p = random_bits = 1
    factor = isqrt(2 * block_size * (1 << width))
    factor += factor * factor < 2 * block_size * (1 << width)
    while Fraction(2 * updates * factor, 1 << p) > target / 2:
        p += 1
        if p > 4096:
            raise ValueError("accuracy request exceeds the bounded planner")
    while Fraction(updates * (max(block_size, 2) - 1), 1 << random_bits) > target / 2:
        random_bits += 1
        if random_bits > 4096:
            raise ValueError("random-bit request exceeds the bounded planner")
    error = Fraction(1, 1 << p)
    result = prefix_error_budget(width, block_size, [error] * updates,
                                 random_bits=random_bits)
    return dict(**result, requested_real_imag_abs_error=error,
                absolute_accuracy_bits=p, target_tv=target,
                machine_mantissa_bits=None)


def rejection_error_budget(envelope, accepted_mass_l1_error):
    """TV <= C*zeta for a perturbed accepted subprobability measure.

    The exact accepted measure is p/C with p normalized. If the approximate
    accepted measure s is nonnegative and ||s-p/C||_1<=zeta, its normalized
    law has TV<=C*zeta. If s has zero mass, C*zeta>=1 and any declared fallback
    is covered only by the trivial bound. Proposal and acceptance errors must
    BOTH be included in zeta; this routine does not estimate them.
    """
    C = _rational(envelope, positive=True)
    error = _rational(accepted_mass_l1_error)
    if C < 1:
        raise ValueError("a normalized rejection target needs envelope >=1")
    return min(Fraction(1), C * error)


def prefix_mass_error_budget(width, weight_error, random_bits):
    """Binary output tree, ABSOLUTE unnormalized child-mass oracle errors.

    Unlike prefix_error_budget these are probability errors, not amplitude
    errors. Ideal-parent-weighted conditional TV is <=2*eta at each node;
    there are 2^t-1 internal nodes, counted without allocating a tree.
    Covered normalized all-zero fallback; exact finite-bit inverse CDF.
    """
    _dimensions(width,2)
    eta = _rational(weight_error)
    if type(random_bits) is not int or not 1 <= random_bits <= 4096:
        raise ValueError("require 1..4096 categorical random bits")
    oracle = 2*((1 << width)-1)*eta
    finite = Fraction(width,1 << random_bits)
    return dict(total_tv_upper_bound=min(Fraction(1),oracle+finite),
                oracle_tv_upper_bound=min(Fraction(1),oracle),
                finite_random_bits_tv_upper_bound=min(Fraction(1),finite),
                weight_abs_error=eta, categorical_random_bits=random_bits,
                requires_uniform_unnormalized_weight_error=True)


def plan_prefix_mass_accuracy(width, target_tv):
    _dimensions(width,2)
    target = _rational(target_tv,positive=True)
    if target >= 1:
        raise ValueError("target TV must be below one")
    bits = random_bits = 1
    while Fraction(2*((1 << width)-1),1 << bits) > target/2:
        bits += 1
        if bits > 4096:
            raise ValueError("weight accuracy exceeds bounded planner")
    while Fraction(width,1 << random_bits) > target/2:
        random_bits += 1
        if random_bits > 4096:
            raise ValueError("random bits exceed bounded planner")
    return dict(**prefix_mass_error_budget(width,Fraction(1,1 << bits),random_bits),
                weight_accuracy_bits=bits,target_tv=target, machine_mantissa_bits=None)


def plan_conditional_weight_accuracy(width,target_tv):
    """Normalized binary weights at every prefix: TV<=2*t*eta+t*2^-L.

    This STRONGER oracle promise cannot be inferred by dividing inaccurate
    prefix masses by their sum. The b=2 scalar component supplies it directly.
    """
    _dimensions(width,2)
    target = _rational(target_tv,positive=True)
    if target >= 1:
        raise ValueError("target TV must be below one")
    bits = random_bits = 1
    while Fraction(2*width,1 << bits) > target/2:
        bits += 1
        if bits > 4096:
            raise ValueError("weight accuracy exceeds bounded planner")
    while Fraction(width,1 << random_bits) > target/2:
        random_bits += 1
        if random_bits > 4096:
            raise ValueError("random bits exceed bounded planner")
    oracle,finite = Fraction(2*width,1 << bits),Fraction(width,1 << random_bits)
    return dict(total_tv_upper_bound=oracle+finite,oracle_tv_upper_bound=oracle,
        finite_random_bits_tv_upper_bound=finite,weight_abs_error=Fraction(1,1 << bits),
        categorical_random_bits=random_bits,weight_accuracy_bits=bits,target_tv=target,
        requires_uniform_normalized_weight_error=True,machine_mantissa_bits=None)
