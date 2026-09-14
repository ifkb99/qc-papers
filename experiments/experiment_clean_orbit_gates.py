"""Clean physical-orbit mixers and reflection rotations, not just basis maps.

Derived BEFORE measuring (TODO 33 / C74). With supplied exact r=b*M and
x=a^(bm+p), U:|x,0> -> |x*a^-p,p> exposes only the small coordinate p.
U^-1 (I tensor W) U implements repeated W and erases p coherently.
Conjugating unit inversion by U implements cell reflection. One Hadamard /
controlled-involution / Hadamard gadget, Rz and its INVERSE exponentiate it
with a second clean ancilla. No full discrete logarithm is evaluated.

P1: all-basis COMPLEX matrices of the capped coordinate/permutation fixtures
    equal their declared full-space extensions, with no per-column phases.
P2: every clean orbit column of repeated W and exp(-i theta R/2) matches an
    independently assembled indexed matrix, including coherent amplitudes.
C1: omitting coordinate uncomputation must leave nonzero scratch probability.
C2: replacing true controlled inversion with minus that branch gives the
    opposite reflection rotation: probabilities alone can miss this defect.
C3: plain physical inversion must disagree with cell reflection (p != 0).

Frozen fixtures (N,a,b,r)=(7,3,3,6),(13,2,3,12),(15,2,2,4).
This is EXPONENTIAL tiny truth-table gate synthesis for validation only,
not an implementation of a scalable arithmetic compiler. Reuses existing
Circuit/statevec and the existing two-level Pauli fixture; no propagator.
Every gate list, batch payload and gate-times-state-entry count is capped.
The arithmetic existence argument is separate from measured fixture costs.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python /usr/bin/python3.12
     --with 'numpy==2.4.6' python -u -m experiments.experiment_clean_orbit_gates
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import traceback

import numpy as np

from circuits import Circuit
from experiments.experiment_physical_phase import (
    controlled_rx_fixture, two_level_rotation_fixture)
from lab import Experiment
import statevec

MAX_BYTES = 16 * 1024 * 1024
MAX_GATES = 30000
MAX_ENTRIES_UPDATED = 100_000_000
TOL = 3e-10


@dataclass(frozen=True)
class Fixture:
    N: int
    a: int
    b: int
    r: int

    def __post_init__(self):
        if (not all(type(v) is int for v in (self.N, self.a, self.b, self.r))
                or not 3 <= self.N <= 15 or not 1 < self.b <= 3
                or not 1 <= self.r <= 12 or self.r % self.b
                or not 1 <= self.a < self.N or math.gcd(self.N, self.a) != 1):
            raise ValueError("outside frozen small physical fixture caps")
        # This tiny order audit is setup/reference work, not order discovery
        # hidden in a purported large-instance algorithm.
        if pow(self.a, self.r, self.N) != 1 or any(
                pow(self.a, d, self.N) == 1 for d in range(1, self.r)):
            raise ValueError("supplied order is not exact")

    @property
    def M(self): return self.r // self.b
    @property
    def nw(self): return (self.N - 1).bit_length()
    @property
    def np(self): return (self.b - 1).bit_length()
    @property
    def nb(self): return self.nw + self.np
    @property
    def work(self): return list(range(self.nw))
    @property
    def fine(self): return list(range(self.nw, self.nb))

    def table(self):
        return {pow(self.a, self.M*p, self.N): p for p in range(self.b)}


def report_path(prefix="clean_orbit_gates"):
    return Path("out") / (prefix + "_" + datetime.now(timezone.utc).strftime(
        "%Y%m%dT%H%M%S%fZ") + ".json")


def gate_guard(qc, additional=0):
    if qc.n > 10 or len(qc.gates) + additional > MAX_GATES:
        raise MemoryError("finite fixture gate/qubit budget exceeded")


def append(qc, part):
    if part.n > qc.n:
        raise ValueError("fixture cannot discard occupied qubits")
    gate_guard(qc, len(part.gates))
    qc.extend(part)


def pattern_phase(qc, qubits, value, angle):
    """Exact exp(i angle |value><value|), including identity Pauli term."""
    k = len(qubits)
    if len(set(qubits)) != k or not 0 <= value < 1 << k:
        raise ValueError("invalid fixed pattern")
    gate_guard(qc, 1 << k)
    for subset in range(1 << k):
        z = sum(1 << q for bit, q in enumerate(qubits) if subset >> bit & 1)
        sign = -1 if (subset & value).bit_count() & 1 else 1
        qc.rot((0, z), -2*angle*sign/(1 << k))


def pair_rx(qc, qubits, left, right, angle):
    """Existing two-level fixture, with its named-CNOT global phase removed."""
    k = len(qubits)
    # Guard before the existing helper allocates its Pauli list.
    gate_guard(qc, (1 << (k-1)) + 6*(k-1) + 1)
    part = two_level_rotation_fixture(left, right, angle, qubits)
    append(qc, part)
    # Circuit.cnot = exp(-i*pi/4) CNOT, and the fixture uses 2*(d-1).
    qc.rot((0, 0), -np.pi*((left ^ right).bit_count()-1))


def pair_swap(qc, qubits, left, right):
    pair_rx(qc, qubits, left, right, np.pi)
    # exp(i*pi/2*(P_left+P_right)): half the diagonal Pauli coefficients
    # cancel EXACTLY between the two projectors. Emit their sum directly,
    # not two larger expansions. This is algebraic cancellation, no cutoff.
    k, diff = len(qubits), left ^ right
    gate_guard(qc, 1 << (k-1))
    for subset in range(1 << k):
        if (subset & diff).bit_count() & 1:
            continue
        z = sum(1 << q for bit, q in enumerate(qubits) if subset >> bit & 1)
        sign = -1 if (subset & left).bit_count() & 1 else 1
        qc.rot((0, z), -np.pi*sign/(1 << (k-1)))


def emit_permutation(qc, qubits, mapping):
    """Capped fixture-only cycle synthesis; NOT a polynomial bit-size compiler."""
    if len(qubits) > 6 or len(mapping) != 1 << len(qubits):
        raise ValueError("tiny permutation cap")
    if sorted(mapping) != list(range(len(mapping))):
        raise ValueError("map is not a permutation")
    seen = set()
    for start in range(len(mapping)):
        if start in seen:
            continue
        cycle, x = [], start
        while x not in seen:
            seen.add(x)
            cycle.append(x)
            x = mapping[x]
        # Chronological (c0,c1),(c0,c2),... realizes c_i -> c_(i+1).
        for other in cycle[1:]:
            pair_swap(qc, qubits, cycle[0], other)


def coordinate_parts(f):
    table = f.table()
    # Full-domain completion: x outside N or not in table has p0(x)=0.
    pvalues = [table.get(pow(x, f.M, f.N), 0) if x < f.N else 0
               for x in range(1 << f.nw)]
    c = Circuit(f.nb)
    for x, p in enumerate(pvalues):
        for bit, target in enumerate(f.fine):
            if p >> bit & 1:
                gate_guard(c, 1 << f.nw)
                controlled_rx_fixture(c, target, f.work,
                                      [(x >> j) & 1 for j in range(f.nw)], np.pi)
                pattern_phase(c, f.work, x, np.pi/2)
    cmapping = [x ^ (pvalues[x & ((1 << f.nw)-1)] << f.nw)
                for x in range(1 << f.nb)]
    multipliers = [pow(f.a, -p, f.N) for p in range(f.b)]
    fmapping = []
    for z in range(1 << f.nb):
        x, p = z & ((1 << f.nw)-1), z >> f.nw
        y = x*multipliers[p] % f.N if x < f.N and p < f.b else x
        fmapping.append(y | (p << f.nw))
    strip = Circuit(f.nb)
    emit_permutation(strip, list(range(f.nb)), fmapping)
    u = Circuit(f.nb)
    append(u, c)
    append(u, strip)
    setup = dict(subgroup_table_entries=f.b, physical_domain_scan=1 << f.nw,
                 full_coordinate_map_entries=2*(1 << f.nb),
                 table_modular_powers=f.b, fine_exponentiations=f.N,
                 inverse_multiplier_powers=f.b, supplied_order_checks=f.r,
                 coordinate_compute_gates=len(c.gates), strip_gates=len(strip.gates),
                 coordinate_gates=len(u.gates),
                 cost_scope="setup scalar/table counters; NOT bit runtime or RSS")
    return c, strip, u, cmapping, fmapping, setup


def fine_matrix(b, theta):
    w = np.eye(b, dtype=complex)
    for left, right, angle in [(0, 1, theta)] + (
            [(1, 2, theta/2)] if b == 3 else []):
        v = np.eye(b, dtype=complex)
        v[left, left] = v[right, right] = np.cos(angle/2)
        v[left, right] = v[right, left] = -1j*np.sin(angle/2)
        w = v @ w
    w[1, :] *= np.exp(1j*theta/3)
    return w


def fine_gate(f, theta):
    qc = Circuit(f.nb)
    pair_rx(qc, f.fine, 0, 1, theta)
    if f.b == 3:
        pair_rx(qc, f.fine, 1, 2, theta/2)
    pattern_phase(qc, f.fine, 1, theta/3)
    return qc


def mixer(f, theta, *, uncompute=True):
    c, strip, u, _, _, _ = coordinate_parts(f)
    qc = Circuit(f.nb)
    append(qc, u)
    append(qc, fine_gate(f, theta))
    append(qc, strip.inverse())
    if uncompute:
        append(qc, c.inverse())
    return qc


def inversion_map(f):
    return [pow(x, -1, f.N) if x < f.N and math.gcd(x, f.N) == 1 else x
            for x in range(1 << f.nw)]


def stripped_reflection(f, theta, *, wrong_control_phase=False):
    """Rotate unit inversion while already in (y,p) coordinates."""
    flag = f.nb
    qc = Circuit(flag+1)
    gadget = Circuit(qc.n).h(flag)
    inv = inversion_map(f)
    mapping = list(range(1 << f.nw)) + [x | (1 << f.nw) for x in inv]
    emit_permutation(gadget, f.work + [flag], mapping)
    if wrong_control_phase:
        # Replace controlled V by controlled (-V); still unitary, same basis
        # probabilities, but its eigenvalue labels (hence rotation) reverse.
        pattern_phase(gadget, [flag], 1, np.pi)
    gadget.h(flag)
    append(qc, gadget)
    qc.rz(flag, theta)
    append(qc, gadget.inverse())
    return qc


def reflection(f, theta, *, wrong_control_phase=False):
    _, _, u, _, _, _ = coordinate_parts(f)
    qc = Circuit(f.nb+1)
    append(qc, u)
    append(qc, stripped_reflection(f, theta, wrong_control_phase=wrong_control_phase))
    append(qc, u.inverse())
    return qc


def controlled_multiplier(qc, f, power, control):
    constant = pow(f.a, power, f.N)
    plain = [x*constant % f.N if x < f.N else x for x in range(1 << f.nw)]
    mapping = list(range(1 << f.nw)) + [x | (1 << f.nw) for x in plain]
    emit_permutation(qc, f.work + [control], mapping)


def checked_run(qc, psi, budget):
    gate_guard(qc)
    # Conservative simultaneous ndarray payload: caller initial/reference,
    # result, Pauli gather/rotation temporaries and index/sign arrays. Gate
    # tuples/Python-object overhead separately bounded by MAX_GATES, not RSS.
    payload = 16*(12*psi.size + 8*psi.shape[0])
    updates = len(qc.gates)*psi.size
    if payload > MAX_BYTES or budget[0]+updates > MAX_ENTRIES_UPDATED:
        raise MemoryError("fixture payload or cumulative gate-entry budget")
    budget[0] += updates
    return statevec.run(qc, psi=psi), dict(payload_bound_bytes=payload,
        gate_count=len(qc.gates), state_entries=psi.size, gate_entry_updates=updates)


def checked_columns(qc, labels, budget):
    # Guard before allocation, not merely before propagation.
    size, columns = 1 << qc.n, len(labels)
    if 16*(12*size*columns + 8*size) > MAX_BYTES:
        raise MemoryError("batch pre-allocation cap")
    psi = np.zeros((size, columns), dtype=complex)
    psi[labels, np.arange(columns)] = 1
    return checked_run(qc, psi, budget)


def physical_output(f, width, theta, *, omit_second=False):
    if (f.N, f.a, f.b, f.r) != (13, 2, 3, 12) or width != 3:
        raise ValueError("full-output test is frozen to N13 / width3")
    offset = f.nb+1
    qc = Circuit(offset+width).x(0)
    for control in range(offset, qc.n):
        qc.h(control)
    append(qc, mixer(f, np.pi/4))
    for i, control in enumerate(range(offset, qc.n)):
        controlled_multiplier(qc, f, 1 << i, control)
        # Adjacent clean mixer/reflection conjugations share a coordinate
        # transform: (U^-1 V U)(U^-1 W U)=U^-1 V W U. No arithmetic moves.
        _, _, u, _, _, _ = coordinate_parts(f)
        append(qc, u)
        append(qc, fine_gate(f, (-1)**i*np.pi/(5+i)))
        if i == 0:
            append(qc, stripped_reflection(f, np.pi/3))
        elif i == 1 and not omit_second:
            append(qc, stripped_reflection(f, theta))
        append(qc, u.inverse())
    qc.qft(list(range(offset, qc.n)), inverse=True)
    # Fixed width was checked before building or allocating.
    if 16*(20*(1 << qc.n)) > MAX_BYTES:
        raise MemoryError("full-output allocation cap")
    psi, costs = checked_run(qc, statevec.basis(qc.n, 0), [0])
    probabilities = np.sum(abs(psi.reshape(1 << width, 1 << offset))**2, axis=1)
    idx = np.arange(len(psi))
    leakage = float(np.sum(abs(psi[(idx >> f.nw) & ((1 << (f.np+1))-1) != 0])**2))
    return dict(probabilities=probabilities.tolist(), leakage=leakage,
                norm=float(np.vdot(psi, psi).real), qubits=qc.n, **costs)


def main():
    exp = Experiment("clean_orbit_gates", doc=__doc__)
    exp.predict("P1", "full-space coordinate/permutation complex matrices match")
    exp.predict("P2", "clean orbit gate columns, including phases, match indexed action")
    exp.must_fail("C1", "omitting scratch erasure leaks probability")
    exp.must_fail("C2", "a wrong controlled global phase reverses reflection rotation")
    exp.must_fail("C3", "plain inversion differs from cell reflection")
    rows, budget = [], [0]
    try:
        for params in ((7, 3, 3, 6), (13, 2, 3, 12), (15, 2, 2, 4)):
            f = Fixture(*params)
            c, strip, u, cmapping, fmapping, setup = coordinate_parts(f)
            primitive_errors = []
            for name, qc, mapping in (("compute", c, cmapping), ("strip", strip, fmapping)):
                actual, cost = checked_columns(qc, list(range(1 << f.nb)), budget)
                expected = np.zeros_like(actual)
                expected[mapping, np.arange(len(mapping))] = 1
                error = float(np.max(abs(actual-expected)))
                primitive_errors.append(error)
                rows.append(dict(series=name, fixture=params, error=error, **cost))
                del actual, expected
            exp.check("P1", max(primitive_errors) < TOL, f"{params}: {primitive_errors}")
            orbit = [pow(f.a, j, f.N) for j in range(f.r)]
            clean_columns = np.array(orbit)
            raw_r = np.zeros((f.r, f.r), complex)
            for j in range(f.r):
                m, p = divmod(j, f.b)
                raw_r[f.b*((-m) % f.M)+p, j] = 1
            for angle in (0., np.pi/7, -np.pi/3):
                for name, qc, target in (
                        ("mixer", mixer(f, angle), np.kron(np.eye(f.M), fine_matrix(f.b, angle))),
                        ("reflection", reflection(f, angle), np.cos(angle/2)*np.eye(f.r)
                         - 1j*np.sin(angle/2)*raw_r)):
                    actual, cost = checked_columns(qc, clean_columns, budget)
                    expected = np.zeros_like(actual)
                    expected[clean_columns, :] = target
                    error = float(np.max(abs(actual-expected)))
                    leakage = float(np.max(np.sum(abs(actual[1 << f.nw:, :])**2, axis=0)))
                    exp.check("P2", error < TOL and leakage < TOL,
                              f"{params} {name} theta/pi={angle/np.pi:g}: err={error:.3g} leak={leakage:.3g}")
                    rows.append(dict(series=name, fixture=params, theta=float(angle),
                                     complex_error=error, leakage=leakage, setup=setup, **cost))
                    del actual, expected
            # Explicit witnesses, not a predicate about a changed circuit's label.
            dirty, cost = checked_columns(mixer(f, np.pi/3, uncompute=False), clean_columns, budget)
            dirty_leak = float(np.max(np.sum(abs(dirty[1 << f.nw:, :])**2, axis=0)))
            del dirty
            bad, badcost = checked_columns(reflection(f, np.pi/3, wrong_control_phase=True), clean_columns, budget)
            true = np.cos(np.pi/6)*np.eye(f.r)-1j*np.sin(np.pi/6)*raw_r
            wrong = np.cos(np.pi/6)*np.eye(f.r)+1j*np.sin(np.pi/6)*raw_r
            bad_error = float(np.max(abs(bad[clean_columns, :]-true)))
            wrong_error = float(np.max(abs(bad[clean_columns, :]-wrong)))
            plain_mismatches = sum(pow(x, -1, f.N) != orbit[f.b*((-(j//f.b)) % f.M)+j%f.b]
                                   for j, x in enumerate(orbit))
            # At M=2 R is identity, so the bad rotation differs by GLOBAL phase
            # only; only N13 M4 is a meaningful coherent wrong-sign witness.
            overlap = np.vdot(true, bad[clean_columns, :]) / f.r
            phase = overlap/abs(overlap) if abs(overlap) else 1
            projective_error = float(np.max(abs(bad[clean_columns, :]-phase*true)))
            exp.fail_check("C1", dirty_leak > .01, f"{params}: scratch={dirty_leak:g}")
            exp.fail_check("C3", plain_mismatches > 0, f"{params}: mismatches={plain_mismatches}")
            if f.M > 2:
                exp.fail_check("C2", projective_error > .1 and wrong_error < TOL,
                               f"{params}: phase-invariant error={projective_error:g}, opposite-sign={wrong_error:.3g}")
            rows.append(dict(series="controls", fixture=params, dirty_leakage=dirty_leak,
                             wrong_rotation_error=bad_error, opposite_sign_error=wrong_error,
                             projective_error=projective_error, plain_mismatches=plain_mismatches,
                             dirty_cost=cost, wrong_cost=badcost))
            del bad
        exp.check("P1", budget[0] <= MAX_ENTRIES_UPDATED,
                  f"charged cumulative Pauli-gate x complex-entry updates={budget[0]}")
    except Exception as exc:
        rows.append(dict(series="exception", error=repr(exc), traceback=traceback.format_exc()))
        exp.check("execution", False, repr(exc))
    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(numpy=np.__version__,
        max_bytes=MAX_BYTES, max_gates=MAX_GATES, max_gate_entry_updates=MAX_ENTRIES_UPDATED,
        charged_gate_entry_updates=budget[0], amplitude_tolerance=TOL,
        synthesis="exponential tiny truth-table fixture, not scalable arithmetic",
        references="full-space basis permutations and independent clean indexed matrices",
        timings="no timing comparison; no hardware stability conclusion"))
    print(f"report: {path}")


if __name__ == "__main__":
    main()
