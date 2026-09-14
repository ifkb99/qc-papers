"""Independent review probes, 2026-09-09; does not alter paper claims.

Predictions before measurement: the Walsh identity survives all probes, but
the universal global-peak formula may fail outside the measured arithmetic
family; odd multiplicative order need not be retained by each output bit;
parity confinement alone does not force exactly equal support-sector sizes.
The single-Toffoli control must fail to refute the peak formula.

Run: uv run python -m experiments.experiment_independent_review
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from circuits import Circuit
from lab import Experiment
from modexp import ModExp
from pauli import to_matrix
from perm_pps import propagate_perm
from pps import propagate
from toffoli_arith import ToffoliModExp
from walsh import classical_permutation, pullback_coefficients, wht


def emit(n, ops):
    circuit = Circuit(n)
    ends = []
    for name, *args in ops:
        getattr(circuit, name)(*args)
        ends.append(len(circuit.gates))
    return circuit, ends


def peak_probe(n, ops, target):
    circuit, ends = emit(n, ops)
    rot = propagate(circuit, {(0, 1 << target): 1.0})
    perm = propagate_perm(circuit, 1 << target)
    m = int(np.argmax(rot.n_terms)) + 1
    gate_index = len(circuit.gates) - m
    op_index = next(i for i, end in enumerate(ends) if end > gate_index)
    if ops[op_index][0] != "toffoli":
        return None
    c = ops[op_index][-1]
    boundary_sets = []
    for step, count in enumerate(perm.n_terms, 1):
        if count == perm.n_max:
            suffix, _ = emit(n, ops[-step:])
            boundary_sets.append(set(propagate_perm(suffix, 1 << target).final_terms))
    predicted = sorted({2 * perm.n_max - sum(not (z & (1 << c)) for z in s)
                        for s in boundary_sets})
    return circuit, rot, perm, m, c, predicted


def main():
    exp = Experiment("independent_review", doc=__doc__)
    exp.predict("P1", "a small circuit can refute the universal global-peak formula")
    exp.predict("P2", "the counterexample still obeys the final Walsh identity")
    exp.predict("P3", "a nonconstant output bit can lose the odd factor of the order")
    exp.predict("P4", "the useful-support fraction need not be exactly 2^-(alpha+1)")
    exp.predict("P5", "the supplied Fourier and Toffoli implementations have different final supports")
    exp.predict("P6", "averaging an actual identity-tail involution is an idempotent, nontrivial map")
    exp.must_fail("C1", "one Toffoli must fail to refute the peak formula")
    exp.must_fail("C2", "averaging an order-three cycle must fail idempotence")

    control = peak_probe(3, [("toffoli", 0, 1, 2)], 2)
    exp.fail_check("C1", control[1].n_max in control[-1],
                   f"one Toffoli: perm={control[2].n_max}, rot={control[1].n_max}, "
                   f"formula={control[-1]}")

    rng = np.random.default_rng(20260909)
    found = None
    for trial in range(1000):
        n = 4
        ops = []
        for _ in range(int(rng.integers(2, 13))):
            a, b, c = map(int, rng.choice(n, 3, replace=False))
            ops.append(("toffoli", a, b, c) if rng.random() < .8
                       else ("cnot", a, b))
        target = int(rng.integers(n))
        result = peak_probe(n, ops, target)
        if result is not None and result[1].n_max not in result[-1]:
            found = (ops, target, result)
            break
    exp.check("P1", found is not None, f"searched {trial + 1} seeded circuits")
    if found:
        ops, target, result = found
        # A deletion-minimal witness also violates the weaker factor-two bound.
        changed = True
        while changed and result[1].n_max > 2 * result[2].n_max:
            changed = False
            for index in range(len(ops)):
                candidate_ops = ops[:index] + ops[index + 1:]
                if not candidate_ops:
                    continue
                candidate = peak_probe(n, candidate_ops, target)
                if candidate is not None and candidate[1].n_max > 2 * candidate[2].n_max:
                    ops, result = candidate_ops, candidate
                    changed = True
                    break
        circuit, rot, perm, m, c, predicted = result
        exp.log(f"ops={ops}; observable=Z{target}; peak target={c}")
        exp.log(f"perm peak={perm.n_max}; rot peak={rot.n_max}; formula values={predicted}")
        coeffs = pullback_coefficients(circuit, target)
        err = max(abs(coeffs[z] - perm.final_terms.get(z, 0.))
                  for z in range(1 << circuit.n))
        suffix = Circuit(circuit.n)
        suffix.gates = circuit.gates[-m:]
        unitary = suffix.to_unitary()
        op = unitary.conj().T @ to_matrix((0, 1 << target), circuit.n) @ unitary
        dense_count = sum(abs(np.trace(to_matrix((x, z), circuit.n) @ op)
                              / (1 << circuit.n)) > 1e-10
                          for x, z in itertools.product(range(1 << circuit.n), repeat=2))
        exp.check("P2", err < 1e-12 and dense_count == rot.n_max,
                  f"final Walsh error={err}; dense peak={dense_count}")
    else:
        exp.check("P2", False, "no counterexample found")

    bit_example = None
    for modulus in range(5, 100, 2):
        for base in range(2, modulus):
            if math.gcd(modulus, base) != 1:
                continue
            orbit = [1]
            value = base % modulus
            while value != 1:
                orbit.append(value)
                value = value * base % modulus
            order = len(orbit)
            if order & (order - 1) == 0:
                continue
            for bit in range(modulus.bit_length()):
                seq = [(v >> bit) & 1 for v in orbit]
                if len(set(seq)) < 2:
                    continue
                period = next(d for d in range(1, order + 1) if order % d == 0
                              and all(seq[k] == seq[k % d] for k in range(order)))
                if period & (period - 1) == 0:
                    bit_example = modulus, base, order, bit, period, seq
                    break
            if bit_example:
                break
        if bit_example:
            break
    exp.check("P3", bit_example is not None, repr(bit_example))
    if bit_example:
        modulus, base, order, bit, period, seq = bit_example
        counts = []
        for width in (4, 8, 12):
            chi = 1 - 2 * np.array([seq[i % order] for i in range(1 << width)])
            counts.append(int(np.count_nonzero(wht(chi))))
        exp.log(f"nonconstant bit Walsh counts at widths 4,8,12: {counts}")

    any_unequal = False
    for modulus, base, width, alpha in [(7, 6, 2, 1), (7, 6, 3, 1), (5, 2, 3, 2)]:
        me = ToffoliModExp(N=modulus, a=base, n_exp=width)
        coeffs = pullback_coefficients(me.build(), me.x[0])
        supp = np.flatnonzero(coeffs != 0)
        exp_mask = sum(1 << q for q in me.exp)
        useful = int(np.count_nonzero((supp & exp_mask) == 0))
        differs = useful * (1 << (alpha + 1)) != len(supp)
        any_unequal |= differs
        exp.log(f"N={modulus},a={base},t={width}: useful={useful}/{len(supp)} "
                f"({useful / len(supp):.10f}); claimed={2. ** -(alpha + 1)}")
    exp.check("P4", any_unequal, "compare counts with exact integer equality")

    fourier = ModExp(N=5, a=2, n_exp=2)
    toffoli = ToffoliModExp(N=5, a=2, n_exp=2)
    fourier_result = propagate(fourier.build(), {(0, 1 << fourier.x[0]): 1.0})
    fourier_support = sum(abs(v) > 1e-10 for v in fourier_result.final_terms.values())
    toffoli_support = int(np.count_nonzero(pullback_coefficients(toffoli.build(), toffoli.x[0])))
    exp.check("P5", fourier_support != toffoli_support,
              f"N=5,a=2,t=2: Fourier q={fourier.n_qubits}, S={fourier_support}; "
              f"Toffoli q={toffoli.n_qubits}, S={toffoli_support}")

    me = ToffoliModExp(N=7, a=6, n_exp=1)
    block = classical_permutation(me.u_a(me.exp[0], 1))
    work = np.arange(1 << me.exp[0])
    images = block[work | (1 << me.exp[0])] & ((1 << me.exp[0]) - 1)
    chi = 1. - 2 * ((work >> me.x[0]) & 1)
    once = (chi + chi[images]) / 2
    twice = (once + once[images]) / 2
    error = float(np.max(np.abs(once - twice)))
    exp.check("P6", error == 0 and bool(np.any(once != chi)),
              f"averaging error={error}; changed entries={np.count_nonzero(once != chi)}")
    cycle, chi = np.array([1, 2, 0]), np.array([1., -1., 1.])
    once = (chi + chi[cycle]) / 2
    twice = (once + once[cycle]) / 2
    error = float(np.max(np.abs(once - twice)))
    exp.fail_check("C2", error > 0, f"order-three averaging error={error}")
    exp.finish()


if __name__ == "__main__":
    main()
