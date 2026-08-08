"""Permutation-native Pauli propagation.

Standard PPS decomposes a Toffoli into Clifford+T and propagates through the
rotations. Those rotations include Hadamards, so the operator leaves the
diagonal mid-circuit even though the Toffoli as a whole is a permutation --
the peak Pauli count is therefore much larger than the final one, and is not
a Walsh quantity.

Treating X / CNOT / Toffoli as ATOMIC permutation gates fixes both problems.
Conjugating a Z-string by a permutation gives a diagonal operator, so the
expansion stays purely Z-type at every step, and the count after k gates is
exactly the Walsh sparsity of the k-gate prefix. Peak memory is then

    N_max = max over prefixes of (Walsh sparsity of that prefix)

which is exact and computable, closing the gap Fable identified.

Conjugation rules for Z^z (all derived from (-1)^{popcount(perm(y) & z)}):

  X(q)            : Z^z -> (-1)^{z_q} Z^z                       1 term
  CNOT(c,t)       : Z^z -> Z^{z XOR (z_t << c)}                 1 term
  Toffoli(a,b,c)  : z_c = 0 -> Z^z                              1 term
                    z_c = 1 -> 1/2 [ Z^z + Z^{z^a} + Z^{z^b}
                                     - Z^{z^a^b} ]              4 terms
  (last line uses (-1)^{y_a y_b} = 1/2(1 + (-1)^{y_a} + (-1)^{y_b}
                                        - (-1)^{y_a + y_b}))
"""
from __future__ import annotations
import numpy as np


class PermPPSResult:
    def __init__(self):
        self.n_terms: list[int] = []
        self.n_max: int = 0
        self.final_terms: dict[int, float] = {}
        self.expectation: float = 0.0
        self.hit_cap: bool = False


def propagate_perm(circuit, zmask: int, delta: float = 0.0,
                   max_terms: int = 4_000_000,
                   max_weight: int | None = None) -> PermPPSResult:
    """Back-propagate the Z-type observable Z^zmask through a classical
    reversible circuit, keeping only Z-strings (keys are z bitmasks).

    `delta`      : coefficient truncation (drop |c| < delta)
    `max_weight` : weight truncation (drop Pauli weight popcount(z) > k), the
                   other standard PPS knob. For permutation circuits the Pauli
                   weight of Z^z is popcount(z), which is exactly the Fourier
                   DEGREE of that Walsh coefficient -- so weight truncation here
                   is literally low-degree Fourier truncation of the pulled-back
                   Boolean function.
    """
    if not circuit.is_classical():
        bad = next(op[0] for op in circuit.logical
                   if op[0] not in ("x", "cnot", "toffoli"))
        raise ValueError(f"circuit is not a permutation (found {bad!r})")

    terms: dict[int, float] = {zmask: 1.0}
    res = PermPPSResult()

    for op in reversed(circuit.logical):
        new: dict[int, float] = {}

        if op[0] == "x":
            q = op[1]
            for z, c in terms.items():
                new[z] = new.get(z, 0.0) + (-c if (z >> q) & 1 else c)

        elif op[0] == "cnot":
            cq, tq = op[1], op[2]
            for z, c in terms.items():
                zz = z ^ (((z >> tq) & 1) << cq)
                new[zz] = new.get(zz, 0.0) + c

        else:                                   # toffoli(a, b, c)
            a, b, cq = op[1], op[2], op[3]
            ba, bb = 1 << a, 1 << b
            for z, c in terms.items():
                if not ((z >> cq) & 1):
                    new[z] = new.get(z, 0.0) + c
                    continue
                h = 0.5 * c
                for zz, s in ((z, h), (z ^ ba, h), (z ^ bb, h), (z ^ ba ^ bb, -h)):
                    new[zz] = new.get(zz, 0.0) + s

        if delta > 0:
            new = {z: v for z, v in new.items() if abs(v) >= delta}
        else:
            new = {z: v for z, v in new.items() if abs(v) > 1e-13}
        if max_weight is not None:
            new = {z: v for z, v in new.items() if z.bit_count() <= max_weight}
        terms = new

        res.n_terms.append(len(terms))
        if len(terms) > max_terms:
            res.hit_cap = True
            break

    res.n_max = max(res.n_terms) if res.n_terms else 0
    res.final_terms = terms
    res.expectation = float(sum(terms.values()))   # <0|Z^z|0> = 1 for all z
    return res
