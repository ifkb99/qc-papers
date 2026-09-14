"""Does a circuit-derived quadratic-cell certificate beat Walsh enumeration?

PREDICTIONS, WRITTEN BEFORE MEASURING.
P1: Freeze m=4, n=14 and a CNOT prefix. Vary only the count p of shared-control
    cubic-generating Toffolis. Four signed physical coefficients agree between
    streamed cells, existing sparse PPS, and a circuit-verified closed formula.
P2: Two cells suffice at every p; exact Walsh support can still grow. Measure
    actual Python allocation for ONE negative-coefficient query, charging
    construction. Whether cells beat the closed formula is an OPEN comparison.
P3: At m=32 a charged two-cell calculation agrees with the same formula without
    enumeration. This is a restricted known-stabilizer class, not novelty or
    a generic circuit result. No sparse run or numerical-floor comparison there.
C1: Squaring/absolutizing cell contributions before their coherent sum must fail
    on the zero coefficient obtained by cancellation between the two cells.
C2: The no-cut quadratic promise must fail once the first cubic gate is present.

All methods receive the same actual Circuit and physical masks. The closed
formula first recognizes its supported gate schedule; it does not receive a
free polynomial oracle. Traces exclude the prebuilt Circuit/interpreter/RSS.
Only short tracemalloc runs: no native timing sweep while TODO34 is unresolved.
"""
from __future__ import annotations

from fractions import Fraction
import gc
from pathlib import Path
import tracemalloc

from circuits import Circuit
from lab.harness import Experiment
from lab.quadratic_cells import quadratic_cell_walsh, QuadraticCellEscape, _Cell
from perm_pps import propagate_perm

OUT = Path(__file__).resolve().parents[1] / "out"


def fixture(m, p, depth=19):
    n, w, target = 3*m+2, 3*m, 3*m+1
    qc = Circuit(n)
    for step in range(depth):
        control, to = (7*step+1) % n, (11*step+5) % n
        if control == to:
            to = (to+1) % n
        qc.cnot(control, to)
    for i in range(p):
        qc.toffoli(w, 3*i+2, 3*i)
    for i in range(m):
        qc.toffoli(3*i, 3*i+1, target)
    return qc, 1 << target, w


def base_to_physical(qc, mask):
    end = next((j for j, op in enumerate(qc.logical) if op[0] != "cnot"), len(qc.logical))
    for j in range(end-1, -1, -1):
        _, control, target = qc.logical[j]
        mask ^= ((mask >> target) & 1) << control
    return mask


def recognized_coefficient(qc, observable, query):
    """Strong special-family baseline, recognizing/charging the actual gates.

    Before the CNOT prefix, f=t+sum_i u_i*v_i+w*sum_(i<p) a_i*v_i.
    Summing each u_i and a_i forces v_i=z_ui and z_ai=w*z_ui (or 0
    for i>=p). Each surviving block contributes 4*(-1)**(z_ui*z_vi).
    The sum over w, and its Fourier sign, MUST remain coherent.
    """
    if qc.n < 5 or (qc.n-2) % 3:
        raise ValueError("unrecognized register layout")
    m, w, target = (qc.n-2)//3, qc.n-2, qc.n-1
    if observable != 1 << target:
        raise ValueError("unrecognized observable")
    prefix = 0
    visits = 0
    for op in qc.logical:
        visits += 1
        if op[0] != "cnot":
            break
        _, control, to = op
        query ^= ((query >> to) & 1) << control
        prefix += 1
    p = len(qc.logical)-prefix-m
    if not 0 <= p <= m:
        raise ValueError("unrecognized nonlinear schedule")
    for i in range(p):
        visits += 1
        if qc.logical[prefix+i] != ("toffoli", w, 3*i+2, 3*i):
            raise ValueError("unrecognized cubic gate")
    for i in range(m):
        visits += 1
        if qc.logical[prefix+p+i] != ("toffoli", 3*i, 3*i+1, target):
            raise ValueError("unrecognized quadratic gate")
    if not ((query >> target) & 1):
        return Fraction(0), visits
    sign = sum(((query >> (3*i)) & 1)*((query >> (3*i+1)) & 1) for i in range(m)) & 1
    total = 0
    for control_value in (0, 1):
        compatible = all(((query >> (3*i+2)) & 1) ==
                         (control_value*((query >> (3*i)) & 1) if i < p else 0)
                         for i in range(m))
        if compatible:
            total += (-1) ** (sign ^ (control_value*((query >> w) & 1)))
    return Fraction(total, 1 << (m+1)), visits


def measure_query(method, qc, observable, query, cut):
    gc.collect()
    tracemalloc.start()
    if method == "sparse":
        result = propagate_perm(qc, observable, max_terms=100000)
        if result.hit_cap:
            raise AssertionError("bounded sparse reference hit cap")
        value = Fraction(result.final_terms.get(query, 0.0))
        counters = {"peak_terms": result.n_max, "final_terms": len(result.final_terms)}
    elif method == "cells":
        result = quadratic_cell_walsh(qc, observable, query, cut_qubits=(cut,))
        value = result.coefficient
        counters = dict(result.counters, branches=result.branches)
    else:
        result = recognized_coefficient(qc, observable, query)
        value = result[0]
        counters = {"recognized_gate_visits": result[1]}
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return value, dict(current_bytes=current, peak_bytes=peak, counters=counters)


exp = Experiment("experiment_quadratic_cell_memory", doc=__doc__)
exp.predict("P1", "all signed queries agree; fixed n and CNOT prefix, only cubic gate count varies")
exp.predict("P2", "two streamed cells suffice; allocation comparison with strong formula remains open")
exp.predict("P3", "large two-cell exact queries match circuit-recognized formula, without enumerating support")
exp.must_fail("C1", "incoherent branch magnitudes cannot reproduce a coherently canceled coefficient")
exp.must_fail("C2", "no-cut quadratic certificate rejects a cubic intermediate sign")

rows = []
for p in range(5):
    qc, observable, cut = fixture(4, p)
    bases = [observable, observable | 3, observable | (1 << cut),
             observable | 1 | 4 | (1 << cut)]
    values = []
    for base in bases:
        query = base_to_physical(qc, base)
        cell = quadratic_cell_walsh(qc, observable, query, cut_qubits=(cut,))
        formula, _ = recognized_coefficient(qc, observable, query)
        sparse = propagate_perm(qc, observable)
        value = Fraction(sparse.final_terms.get(query, 0.))
        values.append(str(cell.coefficient))
        exp.check("P1", cell.coefficient == formula == value,
                  f"p={p}, query={query}, coefficient={cell.coefficient}")
    query = base_to_physical(qc, observable | 3)
    measurements = {}
    for method in ("sparse", "cells", "formula"):
        value, measurement = measure_query(method, qc, observable, query, cut)
        measurements[method] = dict(measurement, coefficient=str(value))
    exp.check("P2", measurements["cells"]["counters"]["branches"] == 2
              and len({m["coefficient"] for m in measurements.values()}) == 1,
              f"p={p}; peaks=" + str({key: v["peak_bytes"] for key, v in measurements.items()}))
    rows.append(dict(kind="fixed_width", m=4, n=qc.n, cubic_gates=p,
                     values=values, measurements=measurements))

qc, observable, cut = fixture(32, 32)
query = base_to_physical(qc, observable | 3)
large = {}
for method in ("cells", "formula"):
    value, measurement = measure_query(method, qc, observable, query, cut)
    large[method] = dict(measurement, coefficient=str(value))
exp.check("P3", large["cells"]["coefficient"] == large["formula"]["coefficient"]
          == str(Fraction(-1, 1 << 33)), str(large))
rows.append(dict(kind="large_restricted", m=32, n=qc.n, cubic_gates=32,
                 measurements=large, sparse_run=False,
                 support_count_derived=(1 << 66) - 3*(1 << 32)))

qc, observable, cut = fixture(1, 1, depth=0)
query = observable | (1 << cut)
contributions = []
for branch in (0, 1):
    cell = _Cell(qc.n, (cut,), branch, observable)
    for op in reversed(qc.logical):
        cell.reverse_gate(op)
    contributions.append(Fraction(cell.signed_sum(query)[0], 1 << qc.n))
exp.fail_check("C1", sum(contributions) == 0 and sum(map(abs, contributions)) > 0,
               f"signed branches={list(map(str, contributions))}")
try:
    quadratic_cell_walsh(qc, observable, query)
except QuadraticCellEscape as exc:
    exp.fail_check("C2", exc.reason == "cubic sign on invariant cell", str(exc))
else:
    exp.fail_check("C2", False, "unsupported cubic circuit was accepted")

exp.finish(report_path=OUT / "quadratic_cell_memory_main.json", rows=rows,
           metadata={"output": "one exact signed full-space physical Walsh coefficient",
                     "memory": "tracemalloc, excludes prebuilt Circuit and RSS",
                     "novelty": "quadratic/stabilizer decomposition is known",
                     "controls": list(map(str, contributions))})
