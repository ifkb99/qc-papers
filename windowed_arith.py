"""Windowed / table-lookup modular exponentiation -- the constructions TODO 12
tests the C29 involution criterion against.

Both classes consume the exponent register `w` bits at a time instead of one bit
at a time, and both compute the same map |x>|1> -> |x>|a^x mod N>. They differ
in what happens *off the valid subspace*, which is exactly where the Walsh
pullback lives (the pullback sums over ALL basis states, scratch included).

  WindowedModExp  -- Gidney-style TABLE LOOKUP. For each exponent window j:
                     look up T[j] = a^(j*2^(kw)) mod N into a scratch register,
                     multiply x by the looked-up *register*, unlook it up. The
                     lookup and the multiply happen for every j INCLUDING j = 0.

  SelectModExp    -- SELECT-multiply. For each nonzero window value j, apply the
                     constant multiplier u_a(a^(j*2^(kw))) controlled on
                     [window == j]; j = 0 is skipped, because multiplying by 1
                     "does nothing". At w = 1 this is exactly ToffoliModExp.

The distinction that matters for C29: in the identity tail (a^(2^(kw)) = 1 mod N,
which happens for every window past v2(r) when r is a power of two) the lookup
table is all-ones, so WindowedModExp's block is a FIXED permutation that never
reads the window register, whereas SelectModExp's block is V controlled on
OR(window) -- the same involution V, but gated by a nonlinear function of the
window bits. See NOTES.md SS WD.

Layout note: both subclass ToffoliModExp and reuse its Cuccaro adder, constant
loader and (for SelectModExp) u_a, then splice their own scratch registers in
front of the exponent register.
"""
from __future__ import annotations

import math

from circuits import Circuit
from toffoli_arith import ToffoliModExp


# ---------------------------------------------------------------------------
# multi-controlled AND, used to activate on a window value
# ---------------------------------------------------------------------------

def and_into(qc: Circuit, ctrls: list[int], out: int, work: list[int]):
    """out ^= AND(ctrls). Needs len(ctrls) - 2 clean qubits in `work`.
    Self-inverse: calling it twice with the same arguments is the identity."""
    k = len(ctrls)
    if k == 0:
        qc.x(out)
    elif k == 1:
        qc.cnot(ctrls[0], out)
    elif k == 2:
        qc.toffoli(ctrls[0], ctrls[1], out)
    else:
        if len(work) < k - 2:
            raise ValueError(f"need {k - 2} work qubits, got {len(work)}")
        qc.toffoli(ctrls[0], ctrls[1], work[0])
        for i in range(2, k - 1):
            qc.toffoli(work[i - 2], ctrls[i], work[i - 1])
        qc.toffoli(work[k - 3], ctrls[k - 1], out)
        for i in range(k - 2, 1, -1):
            qc.toffoli(work[i - 2], ctrls[i], work[i - 1])
        qc.toffoli(ctrls[0], ctrls[1], work[0])


class _WindowedBase(ToffoliModExp):
    """Shared layout: ToffoliModExp's registers, plus activation ancillas."""

    def __init__(self, N: int, a: int, n_exp: int | None = None, w: int = 2,
                 extra: int = 0):
        if w < 1:
            raise ValueError("window width must be >= 1")
        super().__init__(N, a, n_exp)
        self.w = w
        nxt = self.anc + 1
        self.s = list(range(nxt, nxt + extra))          # lookup output (may be [])
        nxt += extra
        # lk[-1] is the activation bit; lk[:-1] are AND-ladder work qubits.
        self.lk = list(range(nxt, nxt + max(1, w - 1)))
        nxt += len(self.lk)
        self.exp = list(range(nxt, nxt + self.n_exp))
        self.n_qubits = nxt + self.n_exp

    # -- window bookkeeping -------------------------------------------------
    def windows(self) -> list[list[int]]:
        return [self.exp[i:i + self.w] for i in range(0, self.n_exp, self.w)]

    def window_base(self, k: int) -> int:
        """The multiplier a^(2^(k*w)) mod N that window k raises to the j-th
        power. Equal to 1 exactly when the window lies in the identity tail."""
        return pow(self.a, 1 << (k * self.w), self.N)

    def is_tail_window(self, k: int) -> bool:
        return self.window_base(k) == 1

    def _activate(self, qc: Circuit, window: list[int], j: int):
        """Toggle the activation bit on [window == j]. Self-inverse."""
        act = self.lk[-1]
        work = self.lk[:-1]
        flips = [q for p, q in enumerate(window) if not ((j >> p) & 1)]
        for q in flips:
            qc.x(q)
        and_into(qc, window, act, work)
        for q in flips:
            qc.x(q)


# ---------------------------------------------------------------------------
# Design B -- table lookup (Gidney-style)
# ---------------------------------------------------------------------------

class WindowedModExp(_WindowedBase):
    """Exponent windowed by QROM lookup + register-register modular multiply.

    Extra register beyond ToffoliModExp: `s`, n qubits, holding the looked-up
    multiplier. `t` stays the Cuccaro operand register; s is copied into it.
    """

    def __init__(self, N: int, a: int, n_exp: int | None = None, w: int = 2):
        super().__init__(N, a, n_exp, w, extra=N.bit_length())

    # -- QROM: s ^= table[j] where j is the window value --------------------
    def lookup(self, qc: Circuit, window: list[int], table: list[int]):
        """s ^= table[window]. Self-inverse, so the same call unlooks it up.

        Exactly one j matches, so the net permutation is s ^= table[j]. When
        every entry of `table` is equal, that permutation does not depend on
        the window register at all -- the point of this whole file.
        """
        if len(table) != (1 << len(window)):
            raise ValueError("table size must be 2^len(window)")
        for j, val in enumerate(table):
            if val == 0:
                continue
            self._activate(qc, window, j)
            for p in range(self.n):
                if (val >> p) & 1:
                    qc.cnot(self.lk[-1], self.s[p])
            self._activate(qc, window, j)

    # -- b += s (mod 2^m), controlled on ctrl -------------------------------
    def _add_s(self, ctrl: int) -> Circuit:
        qc = Circuit(self.n_qubits)
        for p in range(self.n):
            qc.toffoli(ctrl, self.s[p], self.t[p])
        self._add_t_into_b(qc)
        for p in range(self.n):
            qc.toffoli(ctrl, self.s[p], self.t[p])
        return qc

    # -- b <- (b + s) mod N, controlled on ctrl -----------------------------
    def c_add_s_mod(self, qc: Circuit, ctrl: int):
        """Same add / subtract-N / conditional-restore shape as
        ToffoliModExp.cc_add_mod, with the constant replaced by register s."""
        msb = self.b[self.m - 1]
        add_s = self._add_s(ctrl)
        qc.extend(add_s)                                  # b += s
        qc.extend(self._add_const(self.N).inverse())      # b -= N
        qc.cnot(msb, self.anc)                            # anc = 1 iff negative
        qc.extend(self._add_const(self.N, (self.anc,)))   # restore if negative
        qc.extend(add_s.inverse())                        # b -= s
        qc.x(msb)
        qc.cnot(msb, self.anc)                            # uncompute anc
        qc.x(msb)
        qc.extend(add_s)                                  # b += s

    # -- b += x * T[window] mod N ------------------------------------------
    def mult_acc(self, window: list[int], table: list[int]) -> Circuit:
        """One shifted lookup per bit of x: the entry (T[j] * 2^i) mod N is
        pre-reduced classically, so the accumulator never leaves [0, N)."""
        qc = Circuit(self.n_qubits)
        for i, xq in enumerate(self.x):
            shifted = [(v * (1 << i)) % self.N for v in table]
            self.lookup(qc, window, shifted)
            self.c_add_s_mod(qc, xq)
            self.lookup(qc, window, shifted)
        return qc

    # -- one windowed block: x <- x * T[window] mod N -----------------------
    def window_block(self, window: list[int], k: int) -> Circuit:
        base = self.window_base(k)
        table = [pow(base, j, self.N) for j in range(1 << len(window))]
        tinv = [pow(v, -1, self.N) for v in table]
        qc = Circuit(self.n_qubits)
        qc.extend(self.mult_acc(window, table))            # b = x*T mod N
        for xq, bq in zip(self.x, self.b[: self.n]):
            qc.swap(xq, bq)                                # x <-> b
        qc.extend(self.mult_acc(window, tinv).inverse())   # clear b
        return qc

    def build(self, init_x: bool = True) -> Circuit:
        qc = Circuit(self.n_qubits)
        if init_x:
            qc.x(self.x[0])
        for k, window in enumerate(self.windows()):
            qc.extend(self.window_block(window, k))
        return qc

    def summary(self) -> str:
        return (f"windowed N={self.N} a={self.a} r={self.order()} w={self.w} | "
                f"n_exp={self.n_exp} qubits={self.n_qubits} "
                f"[b={self.b} t={self.t} x={self.x} c0={self.c0} "
                f"anc={self.anc} s={self.s} lk={self.lk} exp={self.exp}]")


# ---------------------------------------------------------------------------
# Design A -- select-multiply ("skip the j = 0 case")
# ---------------------------------------------------------------------------

class SelectModExp(_WindowedBase):
    """Exponent windowed by selecting a constant multiplier per window value.

    No lookup register: the multiplier is a compile-time constant per branch.

    `skip_zero` is the whole experiment. With skip_zero=True (the natural
    optimisation -- multiplying by 1 "does nothing", so don't emit the branch)
    the identity-tail block fires for j != 0 only, i.e. V gated on OR(window).
    With skip_zero=False the j = 0 branch is emitted too; in the tail all
    2^w branches are the SAME block u_a(.,1), exactly one of them fires for any
    window value, and the block degenerates to an UNCONTROLLED V -- which is
    what the lookup design does. Same arithmetic either way.
    """

    def __init__(self, N: int, a: int, n_exp: int | None = None, w: int = 2,
                 skip_zero: bool = True):
        self.skip_zero = skip_zero
        super().__init__(N, a, n_exp, w, extra=0)

    def window_block(self, window: list[int], k: int) -> Circuit:
        base = self.window_base(k)
        qc = Circuit(self.n_qubits)
        for j in range(1 if self.skip_zero else 0, 1 << len(window)):
            self._activate(qc, window, j)
            qc.extend(self.u_a(self.lk[-1], pow(base, j, self.N)))
            self._activate(qc, window, j)
        return qc

    def build(self, init_x: bool = True) -> Circuit:
        qc = Circuit(self.n_qubits)
        if init_x:
            qc.x(self.x[0])
        for k, window in enumerate(self.windows()):
            qc.extend(self.window_block(window, k))
        return qc

    def summary(self) -> str:
        return (f"select N={self.N} a={self.a} r={self.order()} w={self.w} "
                f"skip_zero={self.skip_zero} | "
                f"n_exp={self.n_exp} qubits={self.n_qubits} "
                f"[b={self.b} t={self.t} x={self.x} c0={self.c0} "
                f"anc={self.anc} lk={self.lk} exp={self.exp}]")


# ---------------------------------------------------------------------------
# verification helpers
# ---------------------------------------------------------------------------

def replay(circuit: Circuit, y: int) -> int:
    """Image of ONE basis state under a classical circuit, O(gates).

    walsh.classical_permutation replays all 2^n states at once; this is the
    single-state version, cheap enough to gate correctness on circuits far too
    wide to hold a permutation array for.
    """
    for op in circuit.logical:
        if op[0] == "x":
            y ^= 1 << op[1]
        elif op[0] == "cnot":
            y ^= ((y >> op[1]) & 1) << op[2]
        elif op[0] == "toffoli":
            y ^= (((y >> op[1]) & 1) & ((y >> op[2]) & 1)) << op[3]
        else:
            raise ValueError(f"non-classical op {op[0]!r}")
    return y


def read(y: int, qubits: list[int]) -> int:
    return sum(((y >> q) & 1) << i for i, q in enumerate(qubits))


def verify_modexp(me, scratch_clean: bool = True) -> bool:
    """x holds a^e mod N for every e, and every scratch register returns to 0.

    Uses `replay` rather than the state vector: these circuits run to ~25
    qubits and 10^4 gates, where sv.run is hopeless but replay is instant.
    """
    qc = me.build()
    scratch = list(me.b) + list(me.t) + [me.c0, me.anc] + list(me.lk) \
        + list(getattr(me, "s", []))
    for e in range(1 << me.n_exp):
        y0 = 0
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1:
                y0 |= 1 << eq
        y = replay(qc, y0)
        if read(y, me.x) != pow(me.a, e, me.N):
            return False
        if read(y, me.exp) != e:
            return False
        if scratch_clean and any((y >> q) & 1 for q in scratch):
            return False
    return True
