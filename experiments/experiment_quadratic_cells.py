"""Bounded exact audit of streamed quadratic-cell Walsh queries (C83).

Integrated from the independently authored board draft, with strengthened
all-cut success and signed nonparallel-cell assertions. Predictions precede
measurement; this is a correctness corpus, not a scaling sweep.
Run: uv run python -m experiments.experiment_quadratic_cells
"""
from __future__ import annotations
from fractions import Fraction
from random import Random
from circuits import Circuit
from walsh import classical_permutation
from lab.harness import Experiment
from lab.quadratic_cells import _Cell, _Quadratic, QuadraticCellEscape, quadratic_cell_walsh

def direct_q_value(q: _Quadratic, y: int) -> int:
    # Evaluate the upper triangle independently of the implementation helper.
    v = q.constant ^ ((q.linear & y).bit_count() & 1)
    for i in range(q.d):
        for j in range(i + 1, q.d):
            v ^= ((q.rows[i] >> j) & 1) & ((y >> i) & 1) & ((y >> j) & 1)
    return v


def exhaustive_q(q: _Quadratic, extra: int = 0) -> int:
    return sum(-1 if (direct_q_value(q, y) ^ ((extra & y).bit_count() & 1)) else 1
               for y in range(1 << q.d))

def exact_wht(a: list[int]) -> list[int]:
    out = a[:]
    h = 1
    while h < len(out):
        for base in range(0, len(out), 2 * h):
            for j in range(base, base + h):
                x, y = out[j], out[j + h]
                out[j], out[j + h] = x + y, x - y
        h *= 2
    return out

def replay_walsh(circuit: Circuit, observable: int) -> list[int]:
    perm = classical_permutation(circuit)
    signs = [-1 if ((int(perm[y]) & observable).bit_count() & 1) else 1
             for y in range(1 << circuit.n)]
    return exact_wht(signs)

def reference_coeff(circuit: Circuit, observable: int, query: int) -> Fraction:
    return Fraction(replay_walsh(circuit, observable)[query], 1 << circuit.n)

exp = Experiment("quadratic_cells_independent_audit", doc=__doc__, exit_on_fail=True)
exp.predict("P1", "gauss_sum agrees with exhaustive integer Boolean sums for seeded quadratic forms through d=5, including affine products with repeated variables")
exp.predict("P2", "successful streamed coefficients on mixed X/CNOT/Toffoli circuits equal exact replay+integer-WHT coefficients, including negative and zero values")
exp.predict("P3", "cutting a Toffoli control and target produces different tangent subspaces with exact signed transport")
exp.must_fail("C1", "the cubic two-Toffoli no-cut case and curved cut[3] case must raise rather than silently approximate; querymask17 must exhibit coherent cancellation unlike an incoherent branch sum")

exp.section("P1 quadratic Gauss kernel and repeated-bit products")
rng = Random(81273)
q_rows = 0
gauss_bad = []
for d in range(0, 6):
    for _ in range(32 if d else 1):
        q = _Quadratic(d, rng.randrange(1 << d), rng.randrange(2))
        for i in range(d):
            for j in range(i + 1, d):
                if rng.randrange(2):
                    q.rows[i] |= 1 << j
                    q.rows[j] |= 1 << i
        for extra in range(1 << min(d, 3)):
            got, _ = q.gauss_sum(extra)
            want = exhaustive_q(q, extra)
            q_rows += 1
            if got != want:
                gauss_bad.append((d, q.linear, q.constant, q.rows, extra, got, want))

product_bad = []
for d in range(1, 6):
    for _ in range(40):
        a0, b0 = rng.randrange(2), rng.randrange(2)
        alpha, beta = rng.randrange(1 << d), rng.randrange(1 << d)
        q = _Quadratic(d)
        q.add_product(a0, alpha, b0, beta)
        for y in range(1 << d):
            want = (a0 ^ ((alpha & y).bit_count() & 1)) & (b0 ^ ((beta & y).bit_count() & 1))
            if direct_q_value(q, y) != want:
                product_bad.append((d, a0, alpha, b0, beta, y, direct_q_value(q, y), want))
                break
exp.check("P1", not gauss_bad and not product_bad,
          f"{q_rows} quadratic/extras; 200 overlapping affine products; gauss_bad={len(gauss_bad)} product_bad={len(product_bad)}")

exp.section("P2 exact circuit coefficients from replay and integer WHT")
cases = []
aff = Circuit(5).x(0).cnot(0, 2).cnot(2, 4).x(3).cnot(4, 1)
cases.extend(("affine", aff, 1, 1, cut) for cut in ((), (0,), (2,), (0, 2), tuple(range(5))))
cases.append(("sign", Circuit(1).x(0), 1, 1, ()))
quad = Circuit(5).x(1).toffoli(0, 1, 3).cnot(3, 4)
cases.extend(("quadratic", quad, 16, query, cut)
             for query in (0, 1, 16, 17)
             for cut in ((), (0,), (1,), (0, 1), tuple(range(5))))
# Seeded actual mixed circuits, with both a difficult no-cut certificate and
# an exhaustive all-output cut certificate.
mix_rng = Random(9917)
for idx in range(8):
    n = 4 + (idx % 4)
    c = Circuit(n)
    c.x(idx % n)
    c.cnot(0, 1)
    c.toffoli(0, 1, 2)
    for _ in range(2 + idx % 3):
        kind = mix_rng.randrange(3)
        if kind == 0:
            c.x(mix_rng.randrange(n))
        elif kind == 1:
            a, b = mix_rng.sample(range(n), 2)
            c.cnot(a, b)
        else:
            a, b, t = mix_rng.sample(range(n), 3)
            c.toffoli(a, b, t)
    obs = 1 << ((idx + 2) % n)
    query = (1 << (idx % n)) | (1 << ((idx + 1) % n))
    cases.extend(((f"mixed{idx}-none", c, obs, query, ()),
                  (f"mixed{idx}-all", c, obs, query, tuple(range(n)))))
cubic = Circuit(5).toffoli(0, 1, 3).toffoli(3, 2, 4)
circuit_pass, circuit_escapes, circuit_negative, circuit_zero = 0, 0, 0, 0
successful_rows = []
unexpected_full_cut_escapes = 0
for label, circuit, obs, query, cut in cases:
    want = reference_coeff(circuit, obs, query)
    try:
        got = quadratic_cell_walsh(circuit, obs, query, cut_qubits=cut)
    except QuadraticCellEscape as exc:
        circuit_escapes += 1
        unexpected_full_cut_escapes += len(cut) == circuit.n
        successful_rows.append(dict(label=label, cut=list(cut), status="escape", reason=exc.reason))
        continue
    circuit_pass += 1
    ok = got.coefficient == want
    if got.coefficient < 0:
        circuit_negative += 1
    if got.coefficient == 0:
        circuit_zero += 1
    if not ok:
        successful_rows.append(dict(label=label, cut=list(cut), status="wrong", got=str(got.coefficient), want=str(want)))
    else:
        successful_rows.append(dict(label=label, cut=list(cut), status="pass", got=str(got.coefficient), want=str(want), branches=got.branches))
exp.check("P2", circuit_pass > 0 and not unexpected_full_cut_escapes and not [r for r in successful_rows if r["status"] == "wrong"] and circuit_negative > 0 and circuit_zero > 0,
          f"{circuit_pass} exact passes, {circuit_escapes} honest escapes, negative={circuit_negative}, zero={circuit_zero}")

exp.section("P3 nonparallel affine cut transport")
transport = Circuit(3).toffoli(0, 1, 2)
branch_rows = []
transport_images_ok = True
transport_tangents = []
for branch in (0, 1):
    old = _Cell(3, (1, 2), branch, 4)
    new = _Cell(3, (1, 2), branch, 4)
    before_rows, before_offset = old.rows[:], old.offset
    new.reverse_gate(("toffoli", 0, 1, 2))
    branch_rows.append((tuple(new.rows), new.offset, new.signed_sum(0)[0]))
    tangent = set()
    for y in range(2):
        x = before_offset
        for q, row in enumerate(before_rows):
            if row and ((row & y).bit_count() & 1):
                x ^= 1 << q
        image = x ^ ((((x >> 0) & 1) & ((x >> 1) & 1)) << 2)
        expected = new.offset
        for q, row in enumerate(new.rows):
            if row and ((row & y).bit_count() & 1):
                expected ^= 1 << q
        tangent.add(expected ^ new.offset)
        transport_images_ok &= image == expected and direct_q_value(new.quadratic, y) == direct_q_value(old.quadratic, y)
    transport_tangents.append(tangent)
transport_ok = transport_tangents[0] != transport_tangents[1]
transport_coeff = quadratic_cell_walsh(transport, 4, 4, cut_qubits=(1, 2))
transport_want = reference_coeff(transport, 4, 4)
exp.check("P3", transport_ok and transport_images_ok and transport_coeff.coefficient == transport_want,
          f"tangent subspaces differ={transport_ok}, pointwise images={transport_images_ok}, branch rows={[r[0] for r in branch_rows]}, coefficient={transport_coeff.coefficient}, reference={transport_want}")

exp.section("C1 deliberate rejection and coherent cancellation controls")
escape_no_cut = False
escape_cut3 = False
try:
    quadratic_cell_walsh(cubic, 16, 17, cut_qubits=())
except QuadraticCellEscape as exc:
    escape_no_cut = exc.reason in ("non-affine cell image", "cubic sign on invariant cell")
try:
    quadratic_cell_walsh(cubic, 16, 17, cut_qubits=(3,))
except QuadraticCellEscape:
    escape_cut3 = True
coherent = quadratic_cell_walsh(cubic, 16, 17, cut_qubits=(0,))
perm = classical_permutation(cubic)
perm_wht = replay_walsh(cubic, 16)
incoherent = 0
for branch in (0, 1):
    subtotal = 0
    for y in range(1 << cubic.n):
        if ((y >> 0) & 1) != branch:
            continue
        subtotal += -1 if ((int(perm[y]) & 16).bit_count() ^ ((y & 17).bit_count() & 1)) else 1
    incoherent += abs(subtotal)
coherent_ok = coherent.coefficient == 0 and perm_wht[17] == 0 and incoherent != 0
exp.fail_check("C1", escape_no_cut and escape_cut3 and coherent_ok,
               f"no-cut escape={escape_no_cut}, cut3 escape={escape_cut3}, cut0 coefficient={coherent.coefficient}, WHT={perm_wht[17]}, incoherent absolute branch sum={incoherent}")

rows = dict(gauss_forms=q_rows, circuit_cases=len(cases), circuit_successes=circuit_pass,
            circuit_escapes=circuit_escapes, negative_coefficients=circuit_negative,
            zero_coefficients=circuit_zero, explicit_cubic_cut0=str(coherent.coefficient),
            explicit_cubic_no_cut_escape=escape_no_cut, explicit_cubic_cut3_escape=escape_cut3,
            explicit_cubic_incoherent_branch_sum=incoherent,
            transport_branch_rows=[list(r[0]) for r in branch_rows],
            transport_tangents=[sorted(t) for t in transport_tangents],
            unexpected_full_cut_escapes=unexpected_full_cut_escapes, successful_rows=successful_rows)
exp.finish(report_path="out/quadratic_cells_main.json",
           rows=rows, metadata={"reference": "walsh.classical_permutation + integer WHT; exhaustive integer forms"})
