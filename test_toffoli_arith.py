"""Layer-by-layer verification of the Toffoli-compiled modular arithmetic.

Mirrors test_modexp.py so the two compilations are held to the same standard.
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from toffoli_arith import ToffoliModExp
import statevec as sv

FAILED = []
def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra else ''}")
    if not cond: FAILED.append(name)


def prep(me, b0=0, x0=None, ctrls=(), exps=()):
    """Circuit that loads the given classical values into registers."""
    qc = Circuit(me.n_qubits)
    for j in range(me.m):
        if (b0 >> j) & 1: qc.x(me.b[j])
    if x0 is not None:
        for j in range(me.n):
            if (x0 >> j) & 1: qc.x(me.x[j])
    for q in ctrls: qc.x(q)
    for q in exps: qc.x(q)
    return qc


def clean(me, idx):
    """t, c0, anc all back to |0>?"""
    return (sv.read_register(idx, me.t) == 0
            and not ((idx >> me.c0) & 1)
            and not ((idx >> me.anc) & 1))


print("\n[A] ripple adder mod 2^m  (b += t, carry-out dropped)")
for N in (5, 7):
    me = ToffoliModExp(N=N, a=3 if N == 7 else 2, n_exp=1)
    M = 1 << me.m
    ok, why = True, ""
    for t0 in range(M):
        for b0 in range(M):
            qc = prep(me, b0=b0)
            for j in range(me.m):
                if (t0 >> j) & 1: qc.x(me.t[j])
            me._add_t_into_b(qc)
            j, amp = sv.peak(sv.run(qc))
            got = sv.read_register(j, me.b)
            if got != (b0 + t0) % M or abs(amp - 1) > 1e-8 or ((j >> me.c0) & 1):
                ok = False; why = f"t={t0} b={b0} got={got} want={(b0+t0)%M}"; break
        if not ok: break
    check(f"m={me.m}: b += t mod 2^m on all {M}x{M} pairs, c0 restored", ok, why)

print("\n[B] controlled constant addition  b += c (mod 2^m)")
me = ToffoliModExp(N=7, a=3, n_exp=1)
M = 1 << me.m
ok, why = True, ""
for c in range(M):
    for b0 in range(0, M, 3):
        for ctrl in (0, 1):
            qc = prep(me, b0=b0, ctrls=(me.exp[0],) if ctrl else ())
            qc.extend(me._add_const(c, (me.exp[0],)))
            j, amp = sv.peak(sv.run(qc))
            got = sv.read_register(j, me.b)
            want = (b0 + c) % M if ctrl else b0
            if got != want or abs(amp - 1) > 1e-8 or not clean(me, j):
                ok = False; why = f"c={c} b={b0} ctrl={ctrl} got={got} want={want}"; break
        if not ok: break
    if not ok: break
check("controlled _add_const, scratch restored", ok, why)

print("\n[C] inverse gives subtraction")
ok, why = True, ""
for c in range(M):
    for b0 in range(M):
        qc = prep(me, b0=b0)
        qc.extend(me._add_const(c).inverse())
        j, amp = sv.peak(sv.run(qc))
        got = sv.read_register(j, me.b)
        if got != (b0 - c) % M or abs(amp - 1) > 1e-8 or not clean(me, j):
            ok = False; why = f"c={c} b={b0} got={got} want={(b0-c)%M}"; break
    if not ok: break
check("_add_const(...).inverse() == subtraction", ok, why)

print("\n[D] doubly-controlled modular addition  b <- (b+c) mod N")
for N, a in ((5, 2), (7, 3)):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    ok, why = True, ""
    for c in range(N):
        for b0 in range(N):
            for ctrl in (0, 1):
                cs = (me.exp[0], me.x[0]) if ctrl else ()
                qc = prep(me, b0=b0, ctrls=cs)
                me.cc_add_mod(qc, me.exp[0], me.x[0], c)
                j, amp = sv.peak(sv.run(qc))
                got = sv.read_register(j, me.b)
                want = (b0 + c) % N if ctrl else b0
                if got != want or abs(amp - 1) > 1e-8 or not clean(me, j):
                    ok = False
                    why = f"N={N} c={c} b={b0} ctrl={ctrl} got={got} want={want}"
                    break
            if not ok: break
        if not ok: break
    check(f"N={N}: cc_add_mod, all ancillas restored", ok, why)

print("\n[E] controlled multiply  |x> -> |a*x mod N>")
for N, a in ((5, 2), (7, 3)):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    ok, why = True, ""
    for x0 in range(N):
        for ctrl in (0, 1):
            qc = prep(me, x0=x0, ctrls=(me.exp[0],) if ctrl else ())
            qc.extend(me.u_a(me.exp[0], a))
            j, amp = sv.peak(sv.run(qc))
            got = sv.read_register(j, me.x)
            want = (a * x0) % N if ctrl else x0
            bclean = sv.read_register(j, me.b) == 0
            if got != want or abs(amp - 1) > 1e-8 or not clean(me, j) or not bclean:
                ok = False
                why = f"N={N} x={x0} ctrl={ctrl} got={got} want={want} bclean={bclean}"
                break
        if not ok: break
    check(f"N={N} a={a}: u_a correct, b/t/c0/anc cleared", ok, why)

print("\n[F] full modular exponentiation  |e> -> |a^e mod N>")
for N, a, ne in ((5, 2, 2), (7, 3, 2)):
    me = ToffoliModExp(N=N, a=a, n_exp=ne)
    ok, why = True, ""
    for e in range(1 << ne):
        qc = prep(me, exps=[me.exp[i] for i in range(ne) if (e >> i) & 1])
        qc.extend(me.build())
        j, amp = sv.peak(sv.run(qc))
        got = sv.read_register(j, me.x)
        want = pow(a, e, N)
        if got != want or abs(amp - 1) > 1e-8:
            ok = False; why = f"e={e} got={got} want={want} amp={amp:.4f}"; break
    check(f"N={N} a={a} r={me.order()} ({me.n_qubits} qubits): a^e mod N", ok, why)

print("\n[G] cross-check: Toffoli and Fourier compilations agree")
from modexp import ModExp
for N, a, ne in ((5, 2, 2), (7, 3, 2)):
    tof = ToffoliModExp(N=N, a=a, n_exp=ne)
    fou = ModExp(N=N, a=a, n_exp=ne)
    agree = True
    for e in range(1 << ne):
        qt = prep(tof, exps=[tof.exp[i] for i in range(ne) if (e >> i) & 1])
        qt.extend(tof.build())
        jt, _ = sv.peak(sv.run(qt))
        qf = Circuit(fou.n_qubits)
        for i in range(ne):
            if (e >> i) & 1: qf.x(fou.exp[i])
        qf.extend(fou.build())
        jf, _ = sv.peak(sv.run(qf))
        if sv.read_register(jt, tof.x) != sv.read_register(jf, fou.x):
            agree = False; break
    check(f"N={N} a={a}: both compilations give identical a^e mod N", agree)

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
