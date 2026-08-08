"""Toffoli-compiled modular exponentiation -- the permutation-form counterpart
to the Fourier-compiled `modexp.ModExp`.

Same logical map |x>|1> -> |x>|a^x mod N>, but every gate is Toffoli / CNOT / X,
so the circuit is a permutation matrix *gate by gate*. That is the property
under test (C3/C6 in ABSTRACT.md): Z-type Pauli strings should stay Z-type under
back-propagation, so PPS support should collapse, where the Fourier compilation
of the identical map blows up.

Arithmetic: Cuccaro et al. (quant-ph/0410184) ripple-carry adder, taken mod 2^m
by dropping the carry-out. Constants are loaded into a scratch register by
(controlled) X gates, so a single uncontrolled adder serves for the controlled
case -- adding zero is a no-op. Modular reduction follows the standard
add / subtract-N / conditional-restore pattern.

Layout (qubit 0 = least significant), m = n+1:
    b   : m  accumulator; b[m-1] is the sign bit after subtraction
    t   : m  scratch operand (holds the constant currently being added)
    x   : n  multiplicand; ends holding a^exp mod N
    c0  : 1  adder carry-in ancilla
    anc : 1  sign flag
    exp : n_exp  exponent register
"""
from __future__ import annotations
import math
from circuits import Circuit


class ToffoliModExp:
    def __init__(self, N: int, a: int, n_exp: int | None = None):
        if math.gcd(a, N) != 1:
            raise ValueError(f"a={a} must be coprime to N={N}")
        self.N = N
        self.a = a % N
        self.n = N.bit_length()
        self.m = self.n + 1
        self.n_exp = n_exp if n_exp is not None else 2 * self.n

        m, n = self.m, self.n
        self.b = list(range(0, m))
        self.t = list(range(m, 2 * m))
        self.x = list(range(2 * m, 2 * m + n))
        self.c0 = 2 * m + n
        self.anc = 2 * m + n + 1
        self.exp = list(range(self.anc + 1, self.anc + 1 + self.n_exp))
        self.n_qubits = self.anc + 1 + self.n_exp

    # -- Cuccaro ripple-carry, taken mod 2^m (carry-out dropped) ------------
    @staticmethod
    def _maj(qc: Circuit, c: int, b: int, a: int):
        qc.cnot(a, b); qc.cnot(a, c); qc.toffoli(c, b, a)

    @staticmethod
    def _uma(qc: Circuit, c: int, b: int, a: int):
        qc.toffoli(c, b, a); qc.cnot(a, c); qc.cnot(c, b)

    def _add_t_into_b(self, qc: Circuit):
        """b += t (mod 2^m). Restores t and c0."""
        m, t, b, c0 = self.m, self.t, self.b, self.c0
        self._maj(qc, c0, b[0], t[0])
        for i in range(1, m):
            self._maj(qc, t[i - 1], b[i], t[i])
        for i in range(m - 1, 0, -1):
            self._uma(qc, t[i - 1], b[i], t[i])
        self._uma(qc, c0, b[0], t[0])

    # -- constant loading into the scratch register ------------------------
    def _load(self, qc: Circuit, value: int, ctrls: tuple[int, ...] = ()):
        """XOR `value` into t, optionally controlled. Self-inverse."""
        v = value % (1 << self.m)
        for j in range(self.m):
            if not ((v >> j) & 1):
                continue
            if len(ctrls) == 0:
                qc.x(self.t[j])
            elif len(ctrls) == 1:
                qc.cnot(ctrls[0], self.t[j])
            elif len(ctrls) == 2:
                qc.toffoli(ctrls[0], ctrls[1], self.t[j])
            else:
                raise ValueError("at most two controls supported")

    def _add_const(self, value: int, ctrls: tuple[int, ...] = ()) -> Circuit:
        """Sub-circuit: b += value (mod 2^m), controlled on `ctrls`."""
        qc = Circuit(self.n_qubits)
        self._load(qc, value, ctrls)
        self._add_t_into_b(qc)
        self._load(qc, value, ctrls)      # unload (self-inverse)
        return qc

    # -- doubly-controlled modular addition  b <- (b + c) mod N ------------
    def cc_add_mod(self, qc: Circuit, c1: int, c2: int, c: int):
        msb = self.b[self.m - 1]
        add_c = self._add_const(c, (c1, c2))
        qc.extend(add_c)                                  # b += c
        qc.extend(self._add_const(self.N).inverse())      # b -= N
        qc.cnot(msb, self.anc)                            # anc = 1 iff negative
        qc.extend(self._add_const(self.N, (self.anc,)))   # restore if negative
        qc.extend(add_c.inverse())                        # b -= c
        qc.x(msb)
        qc.cnot(msb, self.anc)                            # uncompute anc
        qc.x(msb)
        qc.extend(add_c)                                  # b += c

    # -- controlled multiply-accumulate  b <- b + a*x mod N ----------------
    def cmult_mod(self, ctrl: int, a: int) -> Circuit:
        qc = Circuit(self.n_qubits)
        for i, xq in enumerate(self.x):
            self.cc_add_mod(qc, ctrl, xq, (a * (1 << i)) % self.N)
        return qc

    # -- controlled  |x> -> |a*x mod N> ------------------------------------
    def u_a(self, ctrl: int, a: int) -> Circuit:
        qc = Circuit(self.n_qubits)
        qc.extend(self.cmult_mod(ctrl, a))
        for xq, bq in zip(self.x, self.b[: self.n]):
            qc.cswap(ctrl, xq, bq)
        qc.extend(self.cmult_mod(ctrl, pow(a, -1, self.N)).inverse())
        return qc

    # -- full modular exponentiation ---------------------------------------
    def build(self, init_x: bool = True) -> Circuit:
        qc = Circuit(self.n_qubits)
        if init_x:
            qc.x(self.x[0])
        for i, eq in enumerate(self.exp):
            qc.extend(self.u_a(eq, pow(self.a, 1 << i, self.N)))
        return qc

    def build_shor(self) -> Circuit:
        qc = Circuit(self.n_qubits)
        qc.x(self.x[0])
        for eq in self.exp:
            qc.h(eq)
        for i, eq in enumerate(self.exp):
            qc.extend(self.u_a(eq, pow(self.a, 1 << i, self.N)))
        qc.qft(self.exp, inverse=True)
        return qc

    def order(self) -> int:
        r, v = 1, self.a % self.N
        while v != 1:
            v = (v * self.a) % self.N
            r += 1
        return r

    def summary(self) -> str:
        return (f"N={self.N} a={self.a} r={self.order()} | n={self.n} "
                f"n_exp={self.n_exp} qubits={self.n_qubits} "
                f"[b={self.b} t={self.t} x={self.x} c0={self.c0} "
                f"anc={self.anc} exp={self.exp}]")
