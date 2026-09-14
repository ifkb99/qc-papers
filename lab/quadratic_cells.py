"""Exact signed Walsh queries by streamed affine-support quadratic cells.

This is a restricted stabilizer/quadratic-form method, not a new simulation
principle. It follows the actual X/CNOT/Toffoli circuit without extracting an
exponential Boolean polynomial or storing all cells. Supplied cut qubits
partition the observable's output coordinates before reverse propagation.
Unsupported intermediate cells raise, even if another cut or a later
cancellation could make the requested answer simple.

The output is ONE full-space physical Walsh coefficient, normalized by 2**n.
It is not a sampled quantum output, a complete coefficient dictionary, or an
arbitrary biased-product expectation. Arithmetic is exact Python integers.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from operator import index


def _bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


class QuadraticCellEscape(ValueError):
    """The supplied partition fails the intermediate quadratic-cell promise."""

    def __init__(self, reason, *, branch=None, reverse_step=None):
        self.reason = reason
        self.branch = branch
        self.reverse_step = reverse_step
        super().__init__(f"{reason}; branch={branch}, reverse_step={reverse_step}")


class _Quadratic:
    """Boolean q(y): constant, linear mask, symmetric zero-diagonal polar rows."""

    def __init__(self, d, linear=0, constant=0):
        self.d = d
        self.linear = linear
        self.constant = constant
        self.rows = [0] * d

    def copy(self):
        other = _Quadratic(self.d, self.linear, self.constant)
        other.rows = self.rows.copy()
        return other

    def value(self, y):
        value = self.constant ^ ((self.linear & y).bit_count() & 1)
        for i in _bits(y):
            value ^= (self.rows[i] & y & ~((1 << (i + 1)) - 1)).bit_count() & 1
        return value

    def derivative(self, direction):
        linear = 0
        for i in _bits(direction):
            linear ^= self.rows[i]
        return self.value(direction) ^ self.constant, linear

    def add_product(self, a0, alpha, b0, beta):
        """XOR (a0+alpha.y)(b0+beta.y), using y_i^2=y_i."""
        self.constant ^= a0 & b0
        self.linear ^= (beta if a0 else 0) ^ (alpha if b0 else 0) ^ (alpha & beta)
        for i in _bits(alpha):
            self.rows[i] ^= beta
        for i in _bits(beta):
            self.rows[i] ^= alpha
        # The two diagonal toggles cancel, even where alpha and beta overlap.

    def gauss_sum(self, extra_linear=0):
        """Return sum_y (-1)**(q(y)+extra_linear.y), plus eliminated pairs.

        Eliminate each hyperbolic pair using sum_ab (-1)**(ab+aA+bB)
        = 2*(-1)**(AB). No truth vector or floating arithmetic is used.
        """
        work = self.copy()
        work.linear ^= extra_linear
        active = (1 << self.d) - 1
        pairs = 0
        while active:
            i = next((i for i in _bits(active) if work.rows[i] & active), None)
            if i is None:
                if work.linear & active:
                    return 0, pairs
                return (-1 if work.constant else 1) << (pairs + active.bit_count()), pairs
            j = (work.rows[i] & active & -(work.rows[i] & active)).bit_length() - 1
            active ^= (1 << i) | (1 << j)
            alpha, beta = work.rows[i] & active, work.rows[j] & active
            a0, b0 = (work.linear >> i) & 1, (work.linear >> j) & 1
            work.linear &= active
            for k in _bits(active):
                work.rows[k] &= active
            work.rows[i] = work.rows[j] = 0
            work.add_product(a0, alpha, b0, beta)
            pairs += 1
        return (-1 if work.constant else 1) << pairs, pairs


def _solve_target(rows, target):
    """Solve Vv=e_target from embedding ROWS V, or return None.

    The embedding is injective by construction; ordinary GF(2) elimination
    charges temporary basis masks. It is performed afresh when needed.
    """
    basis = {}
    for i, original in enumerate(rows):
        mask, rhs = original, int(i == target)
        while mask:
            pivot = mask.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (mask, rhs)
                break
            previous, value = basis[pivot]
            mask ^= previous
            rhs ^= value
        if not mask and rhs:
            return None
    value = 0
    for pivot in sorted(basis):
        mask, rhs = basis[pivot]
        if rhs ^ ((mask & value).bit_count() & 1):
            value |= 1 << pivot
    return value


class _Cell:
    """Indicator of x=offset+V y, with sign (-1)**q(y); V stays injective."""

    def __init__(self, n, cut, branch, observable):
        cut_set = set(cut)
        self.rows = [0] * n
        self.offset = sum(((branch >> j) & 1) << q for j, q in enumerate(cut))
        free = 0
        for q in range(n):
            if q not in cut_set:
                self.rows[q] = 1 << free
                free += 1
        linear = 0
        for q in _bits(observable):
            linear ^= self.rows[q]
        self.quadratic = _Quadratic(free, linear, (observable & self.offset).bit_count() & 1)

    def reverse_gate(self, op):
        if op[0] == "x":
            self.offset ^= 1 << op[1]
            return "x"
        if op[0] == "cnot":
            control, target = op[1:]
            self.rows[target] ^= self.rows[control]
            self.offset ^= ((self.offset >> control) & 1) << target
            return "cnot"
        a, b, target = op[1:]
        a0, b0 = (self.offset >> a) & 1, (self.offset >> b) & 1
        alpha, beta = self.rows[a], self.rows[b]
        if not alpha or not beta or alpha == beta:
            # T restricted to this cell is affine. Transport its embedding;
            # other cells can acquire different tangent spaces.
            self.rows[target] ^= (beta if a0 else 0) ^ (alpha if b0 else 0) ^ (alpha & beta)
            self.offset ^= (a0 & b0) << target
            return "toffoli_affine"
        direction = _solve_target(self.rows, target)
        if direction is None:
            raise QuadraticCellEscape("non-affine cell image")
        # The cell is invariant. Its sign must also remain quadratic.
        delta, linear = self.quadratic.derivative(direction)
        if linear == 0:
            mu, nu = 0, 0
        elif linear == alpha:
            mu, nu = 1, 0
        elif linear == beta:
            mu, nu = 0, 1
        elif linear == alpha ^ beta:
            mu, nu = 1, 1
        else:
            raise QuadraticCellEscape("cubic sign on invariant cell")
        gamma = delta ^ (mu & a0) ^ (nu & b0) ^ mu ^ nu
        if gamma:
            self.quadratic.add_product(a0, alpha, b0, beta)
        return "toffoli_quadratic"

    def signed_sum(self, query):
        linear = 0
        for q in _bits(query):
            linear ^= self.rows[q]
        value, pairs = self.quadratic.gauss_sum(linear)
        if (query & self.offset).bit_count() & 1:
            value = -value
        return value, pairs


@dataclass(frozen=True)
class QuadraticCellResult:
    coefficient: Fraction
    branches: int
    free_dimension: int
    counters: dict[str, int]


def quadratic_cell_walsh(circuit, observable: int, query: int, *,
                         cut_qubits=(), max_branches=65536,
                         max_qubits=1024) -> QuadraticCellResult:
    """One exact physical coefficient of U^dagger Z^observable U.

    Supplied cut qubits refer to output-side coordinates where the reverse
    walk starts. Every one of 2**k branches replays the complete circuit and
    computes a signed Gauss sum; no order/oracle/polynomial is supplied free.
    Only one O(n**2)-bit cell, solver and Gauss workspace are live at a time,
    plus integer scalars and Python objects. Circuit storage is additional.
    Gate scan/replay and arithmetic costs, including O(n**3) dense binary
    elimination per needed solve/query, are charged for EVERY branch.

    On a certificate failure, raises with the first encountered branch and
    reverse step; a partial sum is never returned. No truncation, trace_plus,
    implicit cut search, sampling, or finite-precision fallback is performed.
    k=n becomes exhaustive basis enumeration; this is not a worst-case
    polynomial-time simulator. A selected coefficient is also the real
    amplitude <query|H^n U^dagger Z^observable U H^n|0>.
    """
    n = index(circuit.n)
    if not 0 <= n <= index(max_qubits):
        raise ValueError("qubit bound exceeded")
    observable, query = index(observable), index(query)
    if any(mask < 0 or mask.bit_length() > n for mask in (observable, query)):
        raise ValueError("observable/query outside circuit")
    cut = tuple(sorted({index(q) for q in cut_qubits}))
    if any(q < 0 or q >= n for q in cut):
        raise ValueError("cut qubit outside circuit")
    branches = 1 << len(cut)
    if branches > index(max_branches):
        raise ValueError("branch budget exceeded before cell allocation")
    if not circuit.is_classical():
        raise ValueError("expected an X/CNOT/Toffoli circuit")
    for op in circuit.logical:
        if len(set(op[1:])) != len(op) - 1 or any(q < 0 or q >= n for q in op[1:]):
            raise ValueError("gate operands must be distinct and in range")
    counts = dict(x=0, cnot=0, toffoli_affine=0, toffoli_quadratic=0,
                  gauss_pairs=0, branch_initializations=0)
    total = 0
    for branch in range(branches):
        cell = _Cell(n, cut, branch, observable)
        counts["branch_initializations"] += 1
        for reverse_step, op in enumerate(reversed(circuit.logical), 1):
            try:
                kind = cell.reverse_gate(op)
            except QuadraticCellEscape as exc:
                raise QuadraticCellEscape(exc.reason, branch=branch,
                                          reverse_step=reverse_step) from None
            counts[kind] += 1
        value, pairs = cell.signed_sum(query)
        total += value
        counts["gauss_pairs"] += pairs
        # Release the old cell before constructing the next; no branch list.
        del cell
    return QuadraticCellResult(Fraction(total, 1 << n), branches,
                               n - len(cut), counts)
