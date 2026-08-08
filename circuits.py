"""Circuits expressed as sequences of Pauli rotations exp(-i*theta*sigma/2).

Clifford gates land on theta = +/- pi/2 and do not branch under PPS.
T / T-dagger land on theta = +/- pi/4 and branch with weight cos = sin = 1/sqrt(2),
which is the maximum possible. That contrast is the point of the experiment.
"""
from __future__ import annotations
import numpy as np
from pauli import rotation_matrix

PI = np.pi


class Circuit:
    """Ordered list of (sigma, theta); gates[0] is applied to the state first."""

    def __init__(self, n: int):
        self.n = n
        self.gates: list[tuple[tuple[int, int], float]] = []
        # Parallel trace at the *logical* level. Named reversible gates record
        # themselves rather than their Clifford+T decomposition, so a circuit
        # built only from X/CNOT/Toffoli can be replayed as a permutation of
        # basis states (see walsh.py). Anything else records ("nonclassical",).
        self.logical: list[tuple] = []
        self._suppress = 0

    def _log(self, op: tuple):
        if self._suppress == 0:
            self.logical.append(op)

    def _atomic(self, op: tuple):
        """Record `op` and suppress logging from its decomposition."""
        self._log(op)
        outer = self._suppress
        self._suppress = outer + 1
        return outer

    # -- primitive Pauli rotations -----------------------------------------
    def rot(self, sigma: tuple[int, int], theta: float):
        self._log(("nonclassical",))
        self.gates.append((sigma, theta)); return self

    def rz(self, q: int, theta: float): return self.rot((0, 1 << q), theta)
    def rx(self, q: int, theta: float): return self.rot((1 << q, 0), theta)

    # -- named gates, each a product of Pauli rotations ---------------------
    def t(self, q: int):    return self.rz(q, PI / 4)
    def tdg(self, q: int):  return self.rz(q, -PI / 4)
    def s(self, q: int):    return self.rz(q, PI / 2)
    def sdg(self, q: int):  return self.rz(q, -PI / 2)
    def z(self, q: int):    return self.rz(q, PI)

    def x(self, q: int):
        outer = self._atomic(("x", q))
        self.rx(q, PI)
        self._suppress = outer
        return self

    def h(self, q: int):
        # H = Rz(pi/2) Rx(pi/2) Rz(pi/2) up to global phase
        self.rz(q, PI / 2); self.rx(q, PI / 2); self.rz(q, PI / 2); return self

    def cnot(self, c: int, t: int):
        # CNOT = Rzx(-pi/2) then Rz_c(pi/2), Rx_t(pi/2), up to global phase
        outer = self._atomic(("cnot", c, t))
        self.rot(((1 << t), (1 << c)), -PI / 2)   # sigma = Z_c X_t
        self.rz(c, PI / 2)
        self.rx(t, PI / 2)
        self._suppress = outer
        return self

    def toffoli(self, a: int, b: int, c: int):
        """Standard Clifford+T decomposition: 7 T gates, 6 CNOTs, 2 H."""
        outer = self._atomic(("toffoli", a, b, c))
        self.h(c)
        self.cnot(b, c); self.tdg(c)
        self.cnot(a, c); self.t(c)
        self.cnot(b, c); self.tdg(c)
        self.cnot(a, c); self.t(b); self.t(c)
        self.h(c)
        self.cnot(a, b); self.t(a); self.tdg(b)
        self.cnot(a, b)
        self._suppress = outer
        return self

    def cphase(self, a: int, b: int, theta: float):
        """diag(1,1,1,e^{i theta}) as three Z-type rotations (up to phase)."""
        self.rz(a, theta / 2)
        self.rz(b, theta / 2)
        self.rot((0, (1 << a) | (1 << b)), -theta / 2)
        return self

    def swap(self, a: int, b: int):
        self.cnot(a, b); self.cnot(b, a); self.cnot(a, b); return self

    def cswap(self, c: int, a: int, b: int):
        self.cnot(b, a); self.toffoli(c, a, b); self.cnot(b, a); return self

    def ccphase(self, c1: int, c2: int, t: int, theta: float):
        """Doubly-controlled phase. Phase applied only when c1=c2=t=1."""
        self.cphase(c2, t, theta / 2)
        self.cnot(c1, c2)
        self.cphase(c2, t, -theta / 2)
        self.cnot(c1, c2)
        self.cphase(c1, t, theta / 2)
        return self

    def inverse(self) -> "Circuit":
        inv = Circuit(self.n)
        inv.gates = [(sigma, -th) for sigma, th in reversed(self.gates)]
        # X / CNOT / Toffoli are involutions, so reversing the trace inverts it.
        inv.logical = list(reversed(self.logical))
        return inv

    def extend(self, other: "Circuit"):
        self.gates += other.gates
        self.logical += other.logical
        return self

    def is_classical(self) -> bool:
        """True if every logical op is a basis-state permutation."""
        return all(op[0] in ("x", "cnot", "toffoli") for op in self.logical)

    def qft(self, qubits: list[int], inverse=False, swaps=True):
        """QFT on `qubits`. With swaps=True the Fourier index is in standard
        bit order, which is what Fourier-space arithmetic requires (bit
        reversal does not commute with addition)."""
        m = len(qubits)
        seq = []
        for j in range(m):
            seq.append(("h", qubits[j], 0.0))
            for k in range(j + 1, m):
                seq.append(("cp", (qubits[k], qubits[j]), PI / (2 ** (k - j))))

        def do_seq(s):
            for kind, q, th in s:
                if kind == "h":
                    self.h(q)
                else:
                    self.cphase(q[0], q[1], th)

        def do_swaps():
            if not swaps:
                return
            for j in range(m // 2):
                self.swap(qubits[j], qubits[m - 1 - j])

        # forward is U = Body . Swaps (swaps applied to the state first): the
        # swap network cancels the bit reversal the H/CP body puts on the
        # *input* index, leaving the plain DFT.
        if inverse:
            do_seq([(kind, q, -th) for kind, q, th in reversed(seq)])
            do_swaps()
        else:
            do_swaps()
            do_seq(seq)
        return self

    # -- bookkeeping --------------------------------------------------------
    def n_nonclifford(self) -> int:
        return sum(1 for _, th in self.gates
                   if abs(abs(np.sin(th)) - 1) > 1e-9 and abs(np.sin(th)) > 1e-9)

    def to_unitary(self) -> np.ndarray:
        U = np.eye(2 ** self.n, dtype=complex)
        for sigma, th in self.gates:
            U = rotation_matrix(sigma, th, self.n) @ U
        return U

    def apply(self, state: np.ndarray) -> np.ndarray:
        for sigma, th in self.gates:
            state = rotation_matrix(sigma, th, self.n) @ state
        return state


# ---------------------------------------------------------------------------
# Reversible arithmetic: Cuccaro et al. (quant-ph/0410184) ripple-carry adder
# ---------------------------------------------------------------------------

def _maj(qc: Circuit, c: int, b: int, a: int):
    qc.cnot(a, b); qc.cnot(a, c); qc.toffoli(c, b, a)


def _uma(qc: Circuit, c: int, b: int, a: int):
    qc.toffoli(c, b, a); qc.cnot(a, c); qc.cnot(c, b)


def ripple_adder(nbits: int) -> tuple[Circuit, dict]:
    """b += a on nbits-wide registers. Qubit layout:
         0            : carry-in ancilla (must start |0>)
         1..nbits     : a
         nbits+1..2n  : b
         2*nbits+1    : carry-out
    Uses 2*nbits Toffolis => 14*nbits T gates.
    """
    n_q = 2 * nbits + 2
    qc = Circuit(n_q)
    c0 = 0
    a = [1 + i for i in range(nbits)]
    b = [1 + nbits + i for i in range(nbits)]
    z = 2 * nbits + 1

    _maj(qc, c0, b[0], a[0])
    for i in range(1, nbits):
        _maj(qc, a[i - 1], b[i], a[i])
    qc.cnot(a[nbits - 1], z)
    for i in range(nbits - 1, 0, -1):
        _uma(qc, a[i - 1], b[i], a[i])
    _uma(qc, c0, b[0], a[0])

    return qc, {"c0": c0, "a": a, "b": b, "z": z, "nbits": nbits}


# ---------------------------------------------------------------------------
# Brickwork circuits: the family the PPS resource framework was derived on
# (Gharibyan et al. Eq. 10 -- Clifford ZZ layer, random-angle X layer)
# ---------------------------------------------------------------------------

def brickwork(n: int, depth: int, rng, theta_x=None, correlated=False) -> Circuit:
    qc = Circuit(n)
    for _ in range(depth):
        for i in range(n - 1):
            qc.rot(((0), (1 << i) | (1 << (i + 1))), -PI / 2)   # ZZ, Clifford
        for i in range(n):
            if theta_x is not None:
                th = theta_x
            elif correlated:
                th = PI / 4
            else:
                th = rng.uniform(-PI / 4, PI / 4)
            qc.rx(i, th)
    return qc
