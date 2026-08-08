"""Overridable modexp builders -- the machinery behind the reduction-variant
experiments, generalising the `Variant` subclass that was copy-pasted through
experiment_reduction2 / affstruct / resid1 / resid2.

A *wrap* is a named conjugation applied around every cc_add_mod call:

    @register_wrap("my_wrap")
    def my_wrap(me, qc):
        qc.toffoli(me.b[me.m - 1], me.t[0], me.anc)

    me = build_modexp(N=7, a=3, n_exp=4, wraps=["my_wrap"])

Because the wrap is applied symmetrically before and after the reduction and
acts on scratch qubits that are 0 on the valid subspace, correctness of
a^e mod N is preserved -- but ALWAYS confirm with verify_correctness(); the
statevector check is cheap at experiment sizes and the assumption is exactly
the kind that fails silently.

Known wraps (NOTES.md SS L2 / SS RS): the residue story is a function of which
*monomials* couple the msb, not how many wraps are stacked -- msb_t0_anc and
msb_t0_t1 share the monomial msb&t0 and cap density at 3/4 (C34); adding
msb_t1_anc introduces a second monomial and removes the cap.
"""
from __future__ import annotations

import numpy as np

from circuits import Circuit
from toffoli_arith import ToffoliModExp
import statevec as sv

WRAPS: dict[str, callable] = {}


def register_wrap(name: str):
    def deco(fn):
        if name in WRAPS:
            raise ValueError(f"wrap {name!r} already registered")
        WRAPS[name] = fn
        return fn
    return deco


@register_wrap("msb_t0_anc")
def _wrap_msb_t0_anc(me, qc):
    """anc ^= msb & t0 -- the original structure-breaking wrap (v4)."""
    qc.toffoli(me.b[me.m - 1], me.t[0], me.anc)


@register_wrap("msb_t0_t1")
def _wrap_msb_t0_t1(me, qc):
    """t1 ^= msb & t0 -- same monomial, different target (v5)."""
    qc.toffoli(me.b[me.m - 1], me.t[0], me.t[1])


@register_wrap("msb_t1_anc")
def _wrap_msb_t1_anc(me, qc):
    """anc ^= msb & t1 -- an INDEPENDENT monomial; removes the 3/4 cap."""
    qc.toffoli(me.b[me.m - 1], me.t[1], me.anc)


class VariantModExp(ToffoliModExp):
    """ToffoliModExp with named wraps conjugated around every reduction."""

    def __init__(self, N: int, a: int, n_exp: int | None = None,
                 wraps: tuple[str, ...] = ()):
        unknown = [w for w in wraps if w not in WRAPS]
        if unknown:
            raise ValueError(f"unknown wraps {unknown}; "
                             f"registered: {sorted(WRAPS)}")
        self.wraps = tuple(wraps)
        super().__init__(N, a, n_exp)

    def _apply_wraps(self, qc: Circuit):
        for name in self.wraps:
            WRAPS[name](self, qc)

    def cc_add_mod(self, qc: Circuit, c1: int, c2: int, c: int):
        self._apply_wraps(qc)
        super().cc_add_mod(qc, c1, c2, c)
        self._apply_wraps(qc)


def build_modexp(N: int, a: int, n_exp: int | None = None,
                 wraps: tuple[str, ...] = (),
                 compilation: str = "toffoli"):
    """Builder returning a modexp object with .build() / .build_shor().

    compilation="toffoli" supports wraps; "fourier" (Beauregard) does not --
    its reduction is phase-space and has no cc_add_mod hook.
    """
    if compilation == "toffoli":
        return VariantModExp(N=N, a=a, n_exp=n_exp, wraps=wraps)
    if compilation == "fourier":
        if wraps:
            raise ValueError("wraps are only supported for the Toffoli "
                             "compilation")
        from modexp import ModExp
        return ModExp(N, a, n_exp=n_exp)
    raise ValueError(f"unknown compilation {compilation!r}")


def verify_correctness(me, tol: float = 1e-8) -> bool:
    """Statevector check: the built circuit computes a^e mod N for every e.
    Run this for EVERY new variant before trusting any measurement of it."""
    for e in range(1 << me.n_exp):
        qc = Circuit(me.n_qubits)
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1:
                qc.x(eq)
        qc.extend(me.build())
        j, amp = sv.peak(sv.run(qc))
        if sv.read_register(j, me.x) != pow(me.a, e, me.N) or abs(amp - 1) > tol:
            return False
    return True
