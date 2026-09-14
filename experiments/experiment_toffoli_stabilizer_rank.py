"""Does C83 recognize more one-cell Toffoli closures than a generic tableau?

Predictions before measurement: generic commutator rank <=2 is exactly C83's
acceptance on the same real affine-quadratic input. Signed joint probability
selects an exact physical Clifford replacement without a split/merge. A
single added target-neighbor phase edge causes rank3 and cubic rejection;
discarding the global minus of the p=1 case must fail a vector comparison.

This is a bounded correctness corpus and algebraic baseline, not a runtime
or allocated-byte benchmark. Stim prepares the same state and supplies exact
Pauli expectations; its tableau storage/construction is additional to the
constant-workspace rank scan. Dense vectors are tiny validation only.
Run: uv run python -m experiments.experiment_toffoli_stabilizer_rank
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from random import Random

import numpy as np
import stim

from lab.harness import Experiment
from lab.quadratic_cells import _Cell, QuadraticCellEscape
from pauli import to_matrix


def rank_certificate(tableau, a, b, t):
    """Only 3n direct tableau entries; <=3 three-bit basis columns retained."""
    basis = [0, 0, 0]
    reads = 0
    for i in range(len(tableau)):
        pa, pb, pt = (tableau.z_output_pauli(i, q) for q in (a, b, t))
        reads += 3
        column = int(pa in (1, 2)) | (int(pb in (1, 2)) << 1) | (int(pt in (2, 3)) << 2)
        for pivot in (2, 1, 0):
            if (column >> pivot) & 1:
                if basis[pivot]:
                    column ^= basis[pivot]
                else:
                    basis[pivot] = column
                    break
    columns = [v for v in basis if v]
    rank = len(columns)
    # Select independent ORIGINAL conditions (not arbitrary combinations).
    rows = [sum(((v >> j) & 1) << i for i, v in enumerate(columns)) for j in range(3)]
    chosen = ()
    if rank == 1:
        chosen = (next(j for j, row in enumerate(rows) if row),)
    elif rank == 2:
        chosen = next((i, j) for i, j in combinations(range(3), 2)
                      if rows[i] and rows[j] and rows[i] != rows[j])
    return rank, chosen, reads


def projector_probability(sim, n, a, b, t):
    total = 1
    for subset in range(1, 8):
        p = stim.PauliString(n)
        if subset & 1:
            p[a] = "Z"
        if subset & 2:
            p[b] = "Z"
        if subset & 4:
            p[t] = "X"
        total += (-1)**subset.bit_count() * sim.peek_observable_expectation(p)
    return Fraction(total, 8)


def replacement(rank, chosen, probability):
    if rank == 3:
        return None
    if probability == 0:
        return "I"
    assert probability == Fraction(1, 1 << rank)
    if rank == 0:
        return "-I"
    if rank == 1:
        return ("Za", "Zb", "Xt")[chosen[0]]
    return {(0, 1): "CZab", (0, 2): "CNOTat", (1, 2): "CNOTbt"}[chosen]


def direct_toffoli(psi, a, b, t):
    index = np.arange(len(psi))
    perm = index ^ (((index >> a) & (index >> b) & 1) << t)
    return psi[perm]


def replacement_vector(psi, label, a, b, t):
    n = (len(psi)-1).bit_length()
    if label == "I":
        return psi.copy()
    if label == "-I":
        return -psi
    if label in ("Za", "Zb", "Xt"):
        return to_matrix({"Za": (0, 1 << a), "Zb": (0, 1 << b), "Xt": (1 << t, 0)}[label], n) @ psi
    idx = np.arange(len(psi))
    if label == "CZab":
        return (1-2*((idx >> a) & (idx >> b) & 1)) * psi
    control = a if label == "CNOTat" else b
    return psi[idx ^ (((idx >> control) & 1) << t)]


def cell_vector(cell):
    q = cell.quadratic
    psi = np.zeros(1 << len(cell.rows), dtype=complex)
    for y in range(1 << q.d):
        x = cell.offset
        for bit, row in enumerate(cell.rows):
            x ^= ((row & y).bit_count() & 1) << bit
        psi[x] = (-1)**q.value(y) / np.sqrt(1 << q.d)
    return psi


def fixture(n, cut, branch, linear, edges, affine=()):
    cell = _Cell(n, cut, branch, 0)
    free = [q for q in range(n) if q not in cut]
    prep = stim.Circuit()
    prep.append("I", [n-1])
    for j, q in enumerate(cut):
        if (branch >> j) & 1:
            prep.append("X", [q])
    for q in free:
        prep.append("H", [q])
    cell.quadratic.linear = linear
    for i, q in enumerate(free):
        if (linear >> i) & 1:
            prep.append("Z", [q])
    for i, j in edges:
        cell.quadratic.rows[i] ^= 1 << j
        cell.quadratic.rows[j] ^= 1 << i
        prep.append("CZ", [free[i], free[j]])
    for op in affine:
        cell.reverse_gate(op)
        prep.append("X" if op[0] == "x" else "CNOT", list(op[1:]))
    return cell, prep


def main():
    exp = Experiment("toffoli_stabilizer_rank", doc=__doc__, exit_on_fail=True)
    exp.predict("P1", "same input cell and Stim tableau; generic rank acceptance equals C83 and signed physical replacements are exact")
    exp.predict("P2", "adding only target-spectator phase edge at fixed n=5 changes rank2 to rank3 and forces C83 cubic escape")
    exp.must_fail("C1", "rank3 Toffoli output cannot be another stabilizer: input overlap is3/4")
    exp.must_fail("C2", "discarding p=1 global minus fails signed vector equality despite identical tableau state")
    exp.must_fail("C3", "accepting invariant support without phase test, or curved image without geometry test, must fail")

    n, a, b, t = 5, 0, 1, 2
    cases = []
    for edge in (0, 1):
        cases.append((f"edge{edge}", *fixture(n, (), 0, 1 << t,
                                             [(t, 3)] if edge else [])))
    # Explicit witnesses collectively exercise all eight exact replacements.
    cases.extend([
        ("inactive", *fixture(n, (a,), 0, 0, [])),
        ("minus_identity", *fixture(n, (a, b), 3, 1, [])),
        ("Za", *fixture(n, (b,), 1, 1 << 1, [])),
        ("Zb", *fixture(n, (a,), 1, 1 << 1, [])),
        ("Xt", *fixture(n, (a, b, t), 3, 0, [])),
        ("CNOTat", *fixture(n, (b, t), 1, 0, [])),
        ("CNOTbt", *fixture(n, (a, t), 1, 0, [])),
        ("curved", *fixture(n, (t,), 0, 0, [])),
    ])
    rng = Random(20260912)
    for trial in range(24):
        cut = tuple(sorted(rng.sample(range(n), trial % 3)))
        d = n-len(cut)
        edges = [(i, j) for i, j in combinations(range(d), 2) if rng.randrange(2)]
        affine = []
        for step in range(5):
            c, target = rng.sample(range(n), 2)
            affine.append(("cnot", c, target))
        affine.append(("x", rng.randrange(n)))
        cases.append((f"random{trial}", *fixture(n, cut, rng.randrange(1 << len(cut)),
                                                rng.randrange(1 << d), edges, affine)))

    rows = []
    for label, cell, prep in cases:
        tableau = stim.Tableau.from_circuit(prep)
        sim = stim.TableauSimulator()
        sim.do(prep)
        psi = cell_vector(cell)
        reference = tableau.to_state_vector(endian="little").astype(complex)
        same_input = bool(abs(abs(np.vdot(psi, reference))-1) < 2e-7)
        rank, chosen, reads = rank_certificate(tableau, a, b, t)
        probability = projector_probability(sim, n, a, b, t)
        out = direct_toffoli(psi, a, b, t)
        updated = deepcopy(cell)
        reason = None
        try:
            updated.reverse_gate(("toffoli", a, b, t))
            cell_accept = True
            cell_error = float(np.max(np.abs(cell_vector(updated)-out)))
        except QuadraticCellEscape as exc:
            cell_accept, cell_error, reason = False, None, exc.reason
        gate = replacement(rank, chosen, probability)
        error = None if gate is None else float(np.max(np.abs(replacement_vector(psi, gate, a, b, t)-out)))
        ok = same_input and (rank <= 2) == cell_accept == (probability != Fraction(1, 8))
        ok &= (error is None or error < 1e-12) and (cell_error is None or cell_error < 1e-12)
        exp.check("P1", ok and reads == 3*n, f"{label}: rank={rank}, p={probability}, replacement={gate}, C83={reason or 'accept'}")
        rows.append(dict(label=label, rank=rank, probability=str(probability), replacement=gate,
                         cell_accept=cell_accept, reason=reason, same_input=same_input,
                         replacement_error=error, cell_error=cell_error, tableau_reads=reads,
                         input_output_overlap=float(np.vdot(psi, out).real)))

    exp.check("P2", rows[0]["rank"] == 2 and rows[1]["rank"] == 3
              and rows[1]["reason"] == "cubic sign on invariant cell", str(rows[:2]))
    exp.fail_check("C1", abs(rows[1]["input_output_overlap"]-.75) < 1e-12,
                   "distinct stabilizer overlap cannot lie between1/sqrt(2) and1")
    minus = next(row for row in rows if row["label"] == "minus_identity")
    exp.fail_check("C2", minus["replacement"] == "-I" and abs(minus["input_output_overlap"]+1) < 1e-12,
                   "omitted global minus changes the input/output overlap from-1 to+1")
    curved = next(row for row in rows if row["label"] == "curved")
    exp.fail_check("C3", not rows[1]["cell_accept"] and curved["reason"] == "non-affine cell image",
                   "invariant cubic and curved images both rejected")
    labels = {row["replacement"] for row in rows if row["replacement"] is not None}
    exp.check("P1", labels == {"I", "-I", "Za", "Zb", "Xt", "CZab", "CNOTat", "CNOTbt"},
              f"all eight signed replacements exercised: {sorted(labels)}")
    exp.finish(report_path="out/toffoli_stabilizer_rank_main.json", rows=rows,
               metadata={"stim_version": stim.__version__, "n": n,
                         "rank_scan": "3n tableau entries, <=3 three-bit basis columns",
                         "extra_costs": "tableau and Stim construction; 7 Pauli expectations for signed probability; tiny dense vectors for verification only",
                         "novelty": "generic stabilizer measurement rank matches C83 acceptance; no recognition advantage established"})


if __name__ == "__main__":
    main()
