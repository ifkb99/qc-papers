"""Layer-by-layer verification of the modular exponentiation generator."""
from __future__ import annotations
import numpy as np, math
from circuits import Circuit
from modexp import ModExp
import statevec as sv

FAILED = []
def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra else ''}")
    if not cond: FAILED.append(name)


print("\n[A] ccphase / cswap primitives")
qc = Circuit(3); qc.ccphase(0, 1, 2, 0.7)
U = qc.to_unitary(); tgt = np.eye(8, dtype=complex); tgt[7, 7] = np.exp(0.7j)
ph = U[0, 0] / tgt[0, 0]
check("ccphase = phase only on |111>", np.allclose(U, ph * tgt, atol=1e-9))

qc = Circuit(3); qc.cswap(0, 1, 2)
def cswap_f(i):
    if not (i & 1): return i
    b1, b2 = (i >> 1) & 1, (i >> 2) & 1
    return (i & 1) | (b2 << 1) | (b1 << 2)
tgt = np.zeros((8, 8), dtype=complex)
for i in range(8): tgt[cswap_f(i), i] = 1
U = qc.to_unitary(); ph = U[cswap_f(0), 0] / 1
check("cswap (Fredkin)", np.allclose(U, ph * tgt, atol=1e-9))

print("\n[B] statevec matches dense to_unitary")
rng = np.random.default_rng(0)
qc = Circuit(4)
for _ in range(25):
    k = rng.integers(0, 4)
    q = int(rng.integers(0, 4))
    if k == 0: qc.h(q)
    elif k == 1: qc.t(q)
    elif k == 2: qc.rx(q, float(rng.uniform(0, np.pi)))
    else:
        t = int(rng.integers(0, 4))
        if t != q: qc.cnot(q, t)
psi_dense = qc.to_unitary() @ sv.basis(4, 0)
psi_fast = sv.run(qc)
check("statevec == dense", np.allclose(psi_dense, psi_fast, atol=1e-10),
      f"max diff {np.max(np.abs(psi_dense-psi_fast)):.2e}")

print("\n[C] Fourier constant adder on the b register (mod 2^m)")
me = ModExp(N=7, a=3, n_exp=1)
M = 1 << me.m
ok = True
for c in range(M):
    for b0 in range(M):
        qc = Circuit(me.n_qubits)
        for j in range(me.m):
            if (b0 >> j) & 1: qc.x(me.b[j])
        me.qft_b(qc); me.phi_add(qc, c); me.qft_b(qc, inverse=True)
        psi = sv.run(qc); j, amp = sv.peak(psi)
        if sv.read_register(j, me.b) != (b0 + c) % M or abs(amp - 1) > 1e-8:
            ok = False; break
    if not ok: break
check(f"phi_add correct on all {M}x{M} (c, b) pairs", ok)

print("\n[D] doubly-controlled modular addition  b <- (b+c) mod N")
for N, a in ((5, 2), (7, 3)):
    me = ModExp(N=N, a=a, n_exp=1)
    ok, why = True, ""
    for c in range(N):
        for b0 in range(N):
            for ctrl in (0, 1):
                qc = Circuit(me.n_qubits)
                if ctrl:
                    qc.x(me.x[0]); qc.x(me.exp[0])
                for j in range(me.m):
                    if (b0 >> j) & 1: qc.x(me.b[j])
                me.qft_b(qc)
                me.cc_phi_add_mod(qc, me.exp[0], me.x[0], c)
                me.qft_b(qc, inverse=True)
                psi = sv.run(qc); j, amp = sv.peak(psi)
                got = sv.read_register(j, me.b)
                want = (b0 + c) % N if ctrl else b0
                if got != want or abs(amp - 1) > 1e-8 or (j >> me.anc) & 1:
                    ok = False; why = f"N={N} c={c} b={b0} ctrl={ctrl} got={got} want={want}"
                    break
            if not ok: break
        if not ok: break
    check(f"N={N}: cc_phi_add_mod + ancilla restored", ok, why)

print("\n[E] controlled multiply  |x> -> |a*x mod N>")
for N, a in ((5, 2), (7, 3), (15, 7)):
    me = ModExp(N=N, a=a, n_exp=1)
    ok, why = True, ""
    for x0 in range(N):
        for ctrl in (0, 1):
            qc = Circuit(me.n_qubits)
            if ctrl: qc.x(me.exp[0])
            for j in range(me.n):
                if (x0 >> j) & 1: qc.x(me.x[j])
            qc.extend(me.u_a(me.exp[0], a))
            psi = sv.run(qc); j, amp = sv.peak(psi)
            got = sv.read_register(j, me.x)
            want = (a * x0) % N if ctrl else x0
            bclean = sv.read_register(j, me.b) == 0 and not ((j >> me.anc) & 1)
            if got != want or abs(amp - 1) > 1e-8 or not bclean:
                ok = False
                why = f"N={N} a={a} x={x0} ctrl={ctrl} got={got} want={want} clean={bclean}"
                break
        if not ok: break
    check(f"N={N} a={a}: u_a correct, b and anc cleared", ok, why)

print("\n[F] full modular exponentiation  |e> -> |a^e mod N>")
for N, a in ((5, 2), (7, 3), (15, 7), (15, 2)):
    me = ModExp(N=N, a=a, n_exp=3)
    ok, why = True, ""
    for e in range(1 << me.n_exp):
        qc = Circuit(me.n_qubits)
        for i, eq in enumerate(me.exp):
            if (e >> i) & 1: qc.x(eq)
        qc.extend(me.build())
        psi = sv.run(qc); j, amp = sv.peak(psi)
        got = sv.read_register(j, me.x)
        want = pow(a, e, N)
        if got != want or abs(amp - 1) > 1e-8:
            ok = False; why = f"N={N} a={a} e={e} got={got} want={want} amp={amp:.4f}"
            break
    check(f"N={N} a={a} (r={me.order()}, {me.n_qubits} qubits): a^e mod N", ok, why)

print("\n[G] period structure survives into superposition")
me = ModExp(N=15, a=7, n_exp=4)
qc = me.build_shor()
psi = sv.run(qc)
probs = np.abs(psi) ** 2
exp_probs = np.zeros(1 << me.n_exp)
for idx in np.nonzero(probs > 1e-12)[0]:
    exp_probs[sv.read_register(idx, me.exp)] += probs[idx]
top = np.argsort(exp_probs)[::-1][:6]
r = me.order()
print(f"  N=15 a=7, true order r={r}, exponent register = {me.n_exp} qubits")
print(f"  {'y':>4} {'prob':>10} {'y*r/2^t':>10}")
for y in top:
    print(f"  {y:4d} {exp_probs[y]:10.5f} {y * r / (1 << me.n_exp):10.3f}")
peaks = [y for y in top if exp_probs[y] > 0.05]
near_int = all(abs(y * r / (1 << me.n_exp) - round(y * r / (1 << me.n_exp))) < 1e-6
               for y in peaks)
check("QFT peaks sit at y*r/2^t ~ integer", near_int and len(peaks) >= 2,
      f"peaks={peaks}")
check("total probability on exponent register = 1",
      abs(exp_probs.sum() - 1) < 1e-8, f"{exp_probs.sum():.10f}")

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
