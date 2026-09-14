"""Finite-sum Fourier samplers for phase-modulated intervals/progressions.

No vector of length 2^width is allocated. All phases enter as rational turns
and are reduced with Python integers before float64 trigonometry. Mathematical
identities are exact; this implementation is not arbitrary-precision sampling.
These are specialized closed-form distributions, not a circuit propagator.
"""
from __future__ import annotations
import math
import numpy as np


def _ints(*values):
    return all(isinstance(v, (int, np.integer)) for v in values)


def _sinpi_ratio(numerator, denominator):
    """sin(pi*numerator/denominator), with exact zeros and reduced angles."""
    a = int(numerator) % (2*denominator)
    if a == 0 or a == denominator:
        return 0.
    sign = -1. if a > denominator else 1.
    a %= denominator
    return sign * float(np.sin(np.pi * (min(a, denominator-a)/denominator)))


def unit_phase(numerator, denominator):
    """exp(2*pi*i*numerator/denominator) with integer modular reduction."""
    if not _ints(numerator, denominator) or denominator < 1:
        raise ValueError("phase must have integer numerator and positive denominator")
    denominator = int(denominator)
    a = int(numerator) % denominator
    if a == 0:
        return 1.+0j
    if 2*a == denominator:
        return -1.+0j
    if a > denominator//2:
        a -= denominator
    angle = 2*np.pi*(a/denominator)
    return complex(np.cos(angle), np.sin(angle))


def geometric_sum(length, numerator, denominator):
    """Sum exp(2*pi*i*numerator*m/denominator), 0<=m<length."""
    if not _ints(length, numerator, denominator) or length < 0 or denominator < 1:
        raise ValueError("invalid finite geometric sum")
    length, numerator, denominator = int(length), int(numerator), int(denominator)
    if length == 0:
        return 0j
    numerator %= denominator
    if numerator == 0:
        return complex(length)
    ratio = _sinpi_ratio(length*numerator, denominator) / _sinpi_ratio(numerator, denominator)
    return ratio * unit_phase((length-1)*numerator, 2*denominator)


def geometric_norm2(length, numerator, denominator):
    """Squared modulus of the same finite sum; phase-free evaluation."""
    if not _ints(length, numerator, denominator) or length < 0 or denominator < 1:
        raise ValueError("invalid finite geometric norm")
    length, numerator, denominator = int(length), int(numerator), int(denominator)
    if length == 0:
        return 0.
    numerator %= denominator
    if numerator == 0:
        return float(length)**2
    ratio = _sinpi_ratio(length*numerator, denominator) / _sinpi_ratio(numerator, denominator)
    return ratio*ratio


def interval_prefix_probability(width, count, phase_num, phase_den, bits, prefix):
    """Probability of LOW Fourier output bits for a phase-modulated interval.

    Input amplitudes are exp(-2*pi*i*phase*m)/sqrt(count), 0<=m<count, in
    dimension T=2^width. Fourier sign is negative (the circuit's inverse QFT).
    Summing unmeasured outputs forces input differences to be multiples of
    d=T/2^bits. The d residue classes have only two possible lengths.
    """
    if (not _ints(width, count, phase_num, phase_den, bits, prefix)
            or not 0 <= width <= 63 or not 0 <= bits <= width or phase_den < 1
            or not 1 <= count <= 1 << int(width) or not 0 <= prefix < 1 << int(bits)):
        raise ValueError("invalid interval/prefix parameters")
    width, count, bits, prefix = int(width), int(count), int(bits), int(prefix)
    phase_num, phase_den = int(phase_num), int(phase_den)
    M, d = 1 << bits, 1 << (width-bits)
    quotient, extra = divmod(count, d)
    numerator = prefix*phase_den + phase_num*d*M
    denominator = M*phase_den
    total = ((d-extra)*geometric_norm2(quotient, numerator, denominator)
             + extra*geometric_norm2(quotient+1, numerator, denominator))
    return total / count / M


def interval_path(width, count, phase_num, phase_den, *, rng=None, output=None):
    """Sample or force a Fourier outcome using two scalar marginals per bit."""
    interval_prefix_probability(width, count, phase_num, phase_den, 0, 0)
    width = int(width)
    if (rng is None) == (output is None):
        raise ValueError("supply exactly one of rng or output")
    if output is not None and (not _ints(output) or not 0 <= output < 1 << width):
        raise ValueError("forced Fourier output out of range")
    prefix, probability = 0, 1.
    for step in range(width):
        weights = np.array([interval_prefix_probability(width, count, phase_num, phase_den,
                                                        step+1, prefix | (bit << step))
                            for bit in (0, 1)])
        total = float(weights.sum())
        if not np.all(np.isfinite(weights)) or np.min(weights) < 0 or total <= 0:
            raise ArithmeticError("invalid interval conditional probabilities")
        probs = weights/total
        bit = int(rng.random() >= probs[0]) if output is None else (int(output) >> step) & 1
        prefix |= bit << step
        probability *= float(probs[bit])
        if probs[bit] == 0:
            return dict(output=int(output), path_probability=0., marginal_queries=2*(step+1))
    return dict(output=prefix, path_probability=probability, marginal_queries=2*width)


def progression_sample(width, start, stride, count, phase_num, phase_den, rng):
    """Fourier sample of a nonwrapping, phase-modulated arithmetic progression.

    Input labels start+stride*m must lie inside [0,2^width). Their amplitudes
    are exp(-2*pi*i*phase*label)/sqrt(count). The start contributes only a
    global Fourier amplitude phase; it remains relevant in coherent sums.
    """
    if (not _ints(width, start, stride, count, phase_num, phase_den)
            or not 0 <= width <= 63 or start < 0 or stride < 1 or count < 1
            or phase_den < 1):
        raise ValueError("invalid progression parameters")
    width, start, stride, count = int(width), int(start), int(stride), int(count)
    L = 1 << width
    if start+stride*(count-1) >= L:
        raise ValueError("progression must not wrap or leave the input register")
    g = math.gcd(stride, L)
    T, beta = L//g, stride//g
    reduced = interval_path(T.bit_length()-1, count, int(phase_num)*stride,
                            int(phase_den), rng=rng)
    residue = (pow(beta, -1, T)*reduced["output"]) % T if T > 1 else 0
    lift = int(rng.integers(g)) if g > 1 else 0
    return dict(output=residue+T*lift, marginal_queries=reduced["marginal_queries"],
                gcd_lifts=g, reduced_dimension=T)
