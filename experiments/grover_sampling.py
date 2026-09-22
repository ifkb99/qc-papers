"""Exact output sampling for the explicit path-predicate Grover family.

Supports 2..128 qubits and 0..8 iterations using integer arithmetic.
See claims/C117.md and notes/GS-structured-grover-sampling.md for scope.
The sampling law assumes independent unbiased input random bits; the seeded
convenience wrapper supplies reproducible PRNG output.
"""
from __future__ import annotations

from fractions import Fraction
import random


def fib_pair(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a, b


def validate(n, steps):
    if type(n) is not int or not 2 <= n <= 128:
        raise ValueError('n must be an integer in2..128')
    if type(steps) is not int or not 0 <= steps <= 8:
        raise ValueError('steps must be an integer in0..8')


def integer_coefficients(n, marked, steps):
    N = 1 << n
    A, B, Q = 1, 0, 1
    for _ in range(steps):
        A, B, Q = (N-4*marked)*A-2*marked*B, N*(2*A+B), N*Q
    return A, B, Q


def rolling_environment(u, v, power):
    return ((1, u, v), (u, power*u, power*v), (v, power*v, power*v))


def contract_environment(E):
    """Generic bond-three MPS suffix contraction: sum_s T_s E T_s^T."""
    Ts = (((1, 0, 0), (0, 1, 1), (0, 1, 0)),
          ((0, 0, 0), (0, 1, -1), (0, 1, 0)))
    return tuple(tuple(sum(T[i][a]*E[a][b]*T[j][b]
                           for T in Ts for a in range(3) for b in range(3))
                       for j in range(3)) for i in range(3))


def branch_row(row, site, bit):
    c, z, o = row
    sign = 1 if bit == 0 else -1
    return (c if bit == 0 else 0, 0 if site == 0 else z+o, sign*z)


def weight(row, E):
    return sum(row[i]*E[i][j]*row[j] for i in range(3) for j in range(3))


def exact_randbelow(total, getrandbits):
    if type(total) is not int or total < 1:
        raise ValueError('total must be a positive integer')
    if total == 1:
        return 0
    k = (total-1).bit_length()
    while True:
        value = getrandbits(k)
        if type(value) is not int or not 0 <= value < (1 << k):
            raise ValueError('getrandbits violated its integer range contract')
        if value < total:
            return value


def choose_bit(w0, w1, getrandbits):
    if type(w0) is not int or type(w1) is not int or min(w0, w1) < 0 or w0+w1 == 0:
        raise ValueError('branch weights must be nonnegative with positive sum')
    if w0 == 0:
        return 1
    if w1 == 0:
        return 0
    return int(exact_randbelow(w0+w1, getrandbits) >= w0)


class StructuredGroverSampler:
    def __init__(self, n, steps=3, method='rolling'):
        validate(n, steps)
        if method not in ('rolling', 'cached'):
            raise ValueError('unknown contraction method')
        self.n, self.steps, self.method = n, steps, method
        self.marked, following = fib_pair(n)
        self.start_pair = (following, self.marked)
        A, B, Q = integer_coefficients(n, self.marked, steps)
        self.initial = (A*(1 << n), B, 0)
        self.denominator = (Q*(1 << n))**2
        self.cached = None
        if method == 'cached':
            self.cached = [((1, 1, 1),)*3]
            for _ in range(n-1):
                self.cached.append(contract_environment(self.cached[-1]))

    def environments(self):
        u, v = self.start_pair
        power = 1 << (self.n-1)
        for remaining in range(self.n-1, -1, -1):
            yield (self.cached[remaining] if self.cached is not None
                   else rolling_environment(u, v, power))
            u, v, power = v, u-v, power >> 1

    def prefix_numerator(self, prefix, length):
        if type(length) is not int or not 0 <= length <= self.n:
            raise ValueError('invalid prefix length')
        if type(prefix) is not int or not 0 <= prefix < (1 << length):
            raise ValueError('invalid prefix value')
        if length == 0:
            # Verify the actual root contraction rather than assuming normalization.
            E = next(self.environments())
            return sum(weight(branch_row(self.initial, 0, b), E) for b in (0, 1))
        row = self.initial
        for i, E in enumerate(self.environments()):
            row = branch_row(row, i, (prefix >> i) & 1)
            if i+1 == length:
                return weight(row, E)
        raise AssertionError('unreachable valid prefix')

    def prefix_probability(self, prefix, length):
        return Fraction(self.prefix_numerator(prefix, length), self.denominator)

    def sample_one(self, getrandbits, *, prefix=0, length=0):
        if self.prefix_numerator(prefix, length) == 0:
            raise ValueError('cannot condition on a zero-probability prefix')
        row, result = self.initial, prefix
        for i, E in enumerate(self.environments()):
            if i < length:
                row = branch_row(row, i, (prefix >> i) & 1)
                continue
            rows = (branch_row(row, i, 0), branch_row(row, i, 1))
            weights = tuple(weight(r, E) for r in rows)
            bit = choose_bit(*weights, getrandbits)
            row = rows[bit]
            result |= bit << i
        return result

    def sample(self, count, getrandbits, *, prefix=0, length=0):
        if type(count) is not int or count < 0:
            raise ValueError('count must be a nonnegative integer')
        # Validate conditioning even when count is zero.
        if self.prefix_numerator(prefix, length) == 0:
            raise ValueError('cannot condition on a zero-probability prefix')
        return [self.sample_one(getrandbits, prefix=prefix, length=length) for _ in range(count)]


def seeded_samples(n, steps=3, count=16, method='rolling', seed=20260921):
    sampler = StructuredGroverSampler(n, steps, method)
    rng = random.Random(seed)
    return sampler.sample(count, rng.getrandbits)
