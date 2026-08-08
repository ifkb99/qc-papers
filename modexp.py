"""Modular exponentiation |x>|1> -> |x>|a^x mod N>, Beauregard construction.

Reference: Beauregard, "Circuit for Shor's algorithm using 2n+3 qubits"
(quant-ph/0205095). Arithmetic is done in Fourier space (Draper), so constant
addition is a layer of phase rotations rather than a Toffoli ripple.

Register layout (qubit 0 = least significant):
    b    : n+1 qubits   accumulator, top bit is the overflow/sign bit
    x    : n   qubits   multiplicand; ends holding a^exp mod N
    anc  : 1   qubit    comparison ancilla, must return to |0>
    exp  : n_exp qubits exponent register

Note on PPS character: the Fourier adder emits *small-angle* Rz/CP rotations
(theta ~ 2*pi*c*2^j/2^m), whereas the Toffoli path in circuits.py emits
theta = pi/4 exactly. Same logical arithmetic, opposite branching profile.
That contrast is thread 3 in NOTES.md.
"""
from __future__ import annotations
import math
import numpy as np
from circuits import Circuit

TWO_PI = 2 * np.pi


class ModExp:
    def __init__(self, N: int, a: int, n_exp: int | None = None):
        if math.gcd(a, N) != 1:
            raise ValueError(f"a={a} must be coprime to N={N}")
        self.N = N
        self.a = a % N
        self.n = N.bit_length()
        self.m = self.n + 1                       # width of the b register
        self.n_exp = n_exp if n_exp is not None else 2 * self.n

        self.b = list(range(self.m))
        self.x = list(range(self.m, self.m + self.n))
        self.anc = self.m + self.n
        self.exp = list(range(self.anc + 1, self.anc + 1 + self.n_exp))
        self.n_qubits = self.anc + 1 + self.n_exp

    # -- Fourier-space constant addition on the b register ------------------
    def _theta(self, c: int, j: int) -> float:
        return TWO_PI * (c % (1 << self.m)) * (1 << j) / (1 << self.m)

    def phi_add(self, qc: Circuit, c: int, sign: int = 1):
        for j in range(self.m):
            qc.rz(self.b[j], sign * self._theta(c, j))

    def c_phi_add(self, qc: Circuit, ctrl: int, c: int, sign: int = 1):
        for j in range(self.m):
            qc.cphase(ctrl, self.b[j], sign * self._theta(c, j))

    def cc_phi_add(self, qc: Circuit, c1: int, c2: int, c: int, sign: int = 1):
        for j in range(self.m):
            qc.ccphase(c1, c2, self.b[j], sign * self._theta(c, j))

    def qft_b(self, qc: Circuit, inverse=False):
        qc.qft(self.b, inverse=inverse)

    # -- doubly-controlled modular addition  b <- (b + c) mod N -------------
    def cc_phi_add_mod(self, qc: Circuit, c1: int, c2: int, c: int):
        msb = self.b[self.m - 1]
        self.cc_phi_add(qc, c1, c2, c, +1)      # b += c
        self.phi_add(qc, self.N, -1)            # b -= N
        self.qft_b(qc, inverse=True)
        qc.cnot(msb, self.anc)                  # anc = 1 iff it went negative
        self.qft_b(qc)
        self.c_phi_add(qc, self.anc, self.N, +1)   # add N back if so
        self.cc_phi_add(qc, c1, c2, c, -1)      # b -= c  (to test the ancilla)
        self.qft_b(qc, inverse=True)
        qc.x(msb)
        qc.cnot(msb, self.anc)                  # uncompute anc
        qc.x(msb)
        self.qft_b(qc)
        self.cc_phi_add(qc, c1, c2, c, +1)      # b += c

    # -- controlled multiply-accumulate  b <- b + a*x mod N ------------------
    def cmult_mod(self, ctrl: int, a: int) -> Circuit:
        qc = Circuit(self.n_qubits)
        self.qft_b(qc)
        for i, xq in enumerate(self.x):
            self.cc_phi_add_mod(qc, ctrl, xq, (a * (1 << i)) % self.N)
        self.qft_b(qc, inverse=True)
        return qc

    # -- controlled  |x> -> |a*x mod N> -------------------------------------
    def u_a(self, ctrl: int, a: int) -> Circuit:
        qc = Circuit(self.n_qubits)
        qc.extend(self.cmult_mod(ctrl, a))                 # b = a*x mod N
        for xq, bq in zip(self.x, self.b[: self.n]):
            qc.cswap(ctrl, xq, bq)                         # swap x and b
        a_inv = pow(a, -1, self.N)
        qc.extend(self.cmult_mod(ctrl, a_inv).inverse())   # clear b
        return qc

    # -- full modular exponentiation ----------------------------------------
    def build(self, init_x: bool = True) -> Circuit:
        qc = Circuit(self.n_qubits)
        if init_x:
            qc.x(self.x[0])                                # x <- |1>
        for i, eq in enumerate(self.exp):
            qc.extend(self.u_a(eq, pow(self.a, 1 << i, self.N)))
        return qc

    def build_shor(self) -> Circuit:
        """Full order-finding skeleton: H on exponent, modexp, inverse QFT."""
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
                f"[b={self.b} x={self.x} anc={self.anc} exp={self.exp}]")
