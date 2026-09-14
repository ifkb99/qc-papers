"""Correctness tests. Everything downstream depends on these passing."""
from __future__ import annotations
import numpy as np, itertools, random
from pauli import pauli_mult, i_sigma_p, commutes, to_matrix
from circuits import Circuit, ripple_adder, brickwork
from pps import propagate, exact_expectation

rng = np.random.default_rng(0)
FAILED = []

def check(name, cond, extra=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' -- ' + extra) if extra and not cond else ''}")
    if not cond: FAILED.append(name)

def same_up_to_phase(A, B, tol=1e-9):
    idx = np.unravel_index(np.argmax(np.abs(B)), B.shape)
    if abs(B[idx]) < tol: return False
    ph = A[idx] / B[idx]
    return abs(abs(ph) - 1) < tol and np.allclose(A, ph * B, atol=tol)

def perm_matrix(n, f):
    M = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i in range(2 ** n): M[f(i), i] = 1
    return M

# --------------------------------------------------------------------------
print("\n[1] Pauli algebra vs dense matrices (n=3, all pairs)")
n = 3
bad = 0
for x1, z1, x2, z2 in itertools.product(range(8), repeat=4):
    p, q = (x1, z1), (x2, z2)
    r, k = pauli_mult(p, q)
    lhs = to_matrix(p, n) @ to_matrix(q, n)
    rhs = (1j ** k) * to_matrix(r, n)
    if not np.allclose(lhs, rhs): bad += 1
check("pauli_mult phase convention", bad == 0, f"{bad} mismatches")

bad = 0
for x1, z1, x2, z2 in itertools.product(range(8), repeat=4):
    p, q = (x1, z1), (x2, z2)
    if commutes(p, q): continue
    r, s = i_sigma_p(p, q)
    lhs = 1j * to_matrix(p, n) @ to_matrix(q, n)
    if not np.allclose(lhs, s * to_matrix(r, n)): bad += 1
check("i*sigma*P is Hermitian with real sign", bad == 0, f"{bad} mismatches")

bad = 0
for x1, z1, x2, z2 in itertools.product(range(8), repeat=4):
    p, q = (x1, z1), (x2, z2)
    A, B = to_matrix(p, n), to_matrix(q, n)
    if np.allclose(A @ B, B @ A) != commutes(p, q): bad += 1
check("symplectic commutation test", bad == 0, f"{bad} mismatches")

# --------------------------------------------------------------------------
print("\n[2] Gate decompositions into Pauli rotations")
qc = Circuit(1); qc.h(0)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
check("H = Rz(pi/2) Rx(pi/2) Rz(pi/2)", same_up_to_phase(qc.to_unitary(), H))

qc = Circuit(1); qc.x(0)
check("X = Rx(pi)", same_up_to_phase(qc.to_unitary(), np.array([[0, 1], [1, 0]], dtype=complex)))

qc = Circuit(1); qc.t(0)
check("T = Rz(pi/4)", same_up_to_phase(qc.to_unitary(),
      np.diag([1, np.exp(1j * np.pi / 4)]).astype(complex)))

qc = Circuit(1); qc.s(0)
check("S = Rz(pi/2)", same_up_to_phase(qc.to_unitary(), np.diag([1, 1j]).astype(complex)))

qc = Circuit(2); qc.cnot(0, 1)
check("CNOT(0,1) decomposition",
      same_up_to_phase(qc.to_unitary(), perm_matrix(2, lambda i: i ^ 2 if i & 1 else i)))

qc = Circuit(2); qc.cnot(1, 0)
check("CNOT(1,0) decomposition",
      same_up_to_phase(qc.to_unitary(), perm_matrix(2, lambda i: i ^ 1 if i & 2 else i)))

qc = Circuit(3); qc.toffoli(0, 1, 2)
check("Toffoli = Clifford+T (7 T gates)",
      same_up_to_phase(qc.to_unitary(),
                       perm_matrix(3, lambda i: i ^ 4 if (i & 1) and (i & 2) else i)))
check("Toffoli T-count is 7", qc.n_nonclifford() == 7, f"got {qc.n_nonclifford()}")

# --------------------------------------------------------------------------
print("\n[3] Ripple-carry adder computes a+b")
for nbits in (2, 3):
    qc, lay = ripple_adder(nbits)
    U = qc.to_unitary()
    ok = True
    for a in range(2 ** nbits):
        for b in range(2 ** nbits):
            idx = sum(((a >> i) & 1) << lay["a"][i] for i in range(nbits)) \
                + sum(((b >> i) & 1) << lay["b"][i] for i in range(nbits))
            out = U[:, idx]
            j = int(np.argmax(np.abs(out)))
            if abs(abs(out[j]) - 1) > 1e-9: ok = False; break
            a_out = sum(((j >> lay["a"][i]) & 1) << i for i in range(nbits))
            b_out = sum(((j >> lay["b"][i]) & 1) << i for i in range(nbits))
            carry = (j >> lay["z"]) & 1
            if (a_out, b_out + (carry << nbits)) != (a, a + b): ok = False; break
        if not ok: break
    check(f"{nbits}-bit adder correct on all {4**nbits} inputs", ok)
    check(f"{nbits}-bit adder T-count = 14*nbits", qc.n_nonclifford() == 14 * nbits,
          f"got {qc.n_nonclifford()}")

# --------------------------------------------------------------------------
print("\n[4] PPS with no truncation reproduces exact expectation values")
py = random.Random(3)
worst = 0.0
for trial in range(6):
    n = py.randint(2, 4)
    qc = Circuit(n)
    for _ in range(py.randint(4, 10)):
        kind = py.choice(["h", "t", "cnot", "rx", "rz"])
        q = py.randrange(n)
        if kind == "cnot" and n > 1:
            t = py.randrange(n)
            while t == q: t = py.randrange(n)
            qc.cnot(q, t)
        elif kind == "h": qc.h(q)
        elif kind == "t": qc.t(q)
        elif kind == "rx": qc.rx(q, py.uniform(0, np.pi))
        else: qc.rz(q, py.uniform(0, np.pi))
    obs = {(0, 1 << py.randrange(n)): 1.0}
    got = propagate(qc, obs, delta=0.0).expectation
    want = exact_expectation(qc, obs)
    worst = max(worst, abs(got - want))
check("PPS(delta=0) == dense expectation", worst < 1e-9, f"max err {worst:.2e}")

# --------------------------------------------------------------------------
print("\n[4b] REGRESSION: theta = pi gates (X/Y/Z) must negate, not no-op")
# The original random-circuit test never emitted an X/Y/Z gate, so a bug that
# treated every sin(theta)=0 gate as the identity went undetected.
for name, build, obs_q in (("X then Z_0", lambda c: c.x(0), 0),
                           ("Z then X_0", lambda c: c.z(0), 0)):
    qc = Circuit(1); build(qc)
    P = (0, 1) if name.endswith("Z_0") else (1, 0)
    got = propagate(qc, {P: 1.0}, delta=0.0).expectation
    want = exact_expectation(qc, {P: 1.0})
    check(f"{name}: PPS == dense", abs(got - want) < 1e-9, f"{got:+.6f} vs {want:+.6f}")

py2 = random.Random(11)
worst = 0.0
for _ in range(8):
    n = py2.randint(1, 4)
    qc = Circuit(n)
    for _ in range(py2.randint(3, 12)):
        kind = py2.choice(["x", "z", "h", "t", "cnot"])
        q = py2.randrange(n)
        if kind == "cnot" and n > 1:
            t = py2.randrange(n)
            while t == q: t = py2.randrange(n)
            qc.cnot(q, t)
        elif kind == "x": qc.x(q)
        elif kind == "z": qc.z(q)
        elif kind == "h": qc.h(q)
        else: qc.t(q)
    obs = {(0, 1 << py2.randrange(n)): 1.0}
    worst = max(worst, abs(propagate(qc, obs, delta=0.0).expectation
                           - exact_expectation(qc, obs)))
check("random circuits INCLUDING X/Z gates", worst < 1e-9, f"max err {worst:.2e}")

print("\n[5] PPS on a real adder circuit vs dense")
qc, lay = ripple_adder(2)
obs = {(0, 1 << lay["b"][0]): 1.0}
got = propagate(qc, obs, delta=0.0).expectation
want = exact_expectation(qc, obs)
check("adder: PPS(delta=0) == dense", abs(got - want) < 1e-9,
      f"{got:.10f} vs {want:.10f}")

print("\n[6] lazy CNOT frame preserves physical operators and filter boundaries")
from perm_pps import propagate_perm, FramedTerms

# Reverse propagation sees CNOT before Toffoli. Both the transformed target
# predicate and one inverse-column offset differ from their internal key bits.
qc = Circuit(3).toffoli(1, 2, 0).cnot(0, 1)
framed = propagate_perm(qc, 6, affine_frame=True)
expected = {7: .5, 5: .5, 3: .5, 1: -.5}
check("nonlinear interleaving uses physical predicate and inverse offsets",
      dict(framed.final_terms) == expected)
unitary = qc.to_unitary()
dense = unitary.conj().T @ to_matrix((0, 6), 3) @ unitary
represented = sum(value * to_matrix((0, key), 3)
                  for key, value in framed.final_terms.items())
check("full framed observable equals independent dense conjugation",
      np.max(np.abs(dense - represented)) < 1e-12)
check("lazy output exposes physical-key lookup, streaming items and values",
      isinstance(framed.final_terms, FramedTerms)
      and framed.final_terms[5] == .5
      and set(framed.final_terms) == set(expected)
      and sum(framed.final_terms.values()) == framed.expectation)

qc = Circuit(2).cnot(0, 1)
check("first CNOT still applies coefficient threshold",
      not propagate_perm(qc, 2, delta=1.1, affine_frame=True).final_terms)
check("trace uses physical bit after CNOT",
      not propagate_perm(qc, 2, trace_plus=[0], affine_frame=True).final_terms)
qc.cnot(0, 1)
check("physical weight truncation is not delayed across a CNOT cancellation",
      not propagate_perm(qc, 2, max_weight=1, affine_frame=True).final_terms)
check("reused |+> control contracts only at its final reverse use",
      dict(propagate_perm(qc, 2, trace_plus=[0],
                          affine_frame=True).final_terms) == {2: 1.})
check("physical X sign after a deferred CNOT",
      dict(propagate_perm(Circuit(2).x(0).cnot(0, 1), 2,
                          affine_frame=True).final_terms) == {3: -1.})
capped = propagate_perm(Circuit(3).toffoli(0, 1, 2), 4,
                        trace_plus=[0], max_terms=2, affine_frame=True)
check("term cap is checked before a shrinking trace",
      capped.hit_cap and capped.n_max == 4 and not capped.trace_events)
wide = propagate_perm(Circuit(65).cnot(64, 0), 1, affine_frame=True)
check("frame and physical lookup retain bits beyond a machine word",
      wide.final_terms[1 | (1 << 64)] == 1.)
check("control: treating internal frame keys as physical keys gives wrong operator",
      dict(propagate_perm(Circuit(2).cnot(0, 1), 2,
                          affine_frame=True).final_terms) != {2: 1.})

print("\n[7] exact quadratic-cell queries preserve signs and reject escaped cells")
from fractions import Fraction
from lab.quadratic_cells import quadratic_cell_walsh, QuadraticCellEscape

qc = Circuit(5).toffoli(0, 1, 3).toffoli(3, 2, 4)
check("signed branch cancellation is exact",
      quadratic_cell_walsh(qc, 16, 17, cut_qubits=(0,)).coefficient == 0)
for cut, reason in (((), "cubic sign on invariant cell"),
                    ((3,), "non-affine cell image")):
    escaped = False
    try:
        quadratic_cell_walsh(qc, 16, 17, cut_qubits=cut)
    except QuadraticCellEscape as exc:
        escaped = exc.reason == reason and exc.reverse_step == 2
    check(f"unsupported intermediate cell rejects for cut={cut}", escaped)

qc = Circuit(3).x(0).cnot(2, 0).toffoli(0, 1, 2)
unitary = qc.to_unitary()
dense = unitary.conj().T @ to_matrix((0, 4), 3) @ unitary
cell_error = 0.0
negative = False
for query in range(8):
    expected = np.trace(to_matrix((0, query), 3) @ dense).real / 8
    for cut in ((), (1, 2), (0, 1, 2)):
        value = quadratic_cell_walsh(qc, 4, query, cut_qubits=cut).coefficient
        cell_error = max(cell_error, abs(float(value) - expected))
        negative |= value < 0
check("no-cut, nonparallel and singleton partitions equal dense coefficients",
      cell_error < 1e-12 and negative, f"max err {cell_error:.2e}")
wide = quadratic_cell_walsh(Circuit(65).x(64).cnot(64, 0),
                            1, 1 | (1 << 64))
check("exact query retains physical sign and masks beyond a machine word",
      wide.coefficient == Fraction(-1))
capped = False
try:
    quadratic_cell_walsh(Circuit(65), 1, 1, cut_qubits=range(65))
except ValueError as exc:
    capped = "branch budget exceeded" in str(exc)
check("branch budget rejects before exponential cell construction", capped)

print("\n" + ("ALL TESTS PASSED" if not FAILED else f"FAILURES: {FAILED}"))
raise SystemExit(bool(FAILED))
