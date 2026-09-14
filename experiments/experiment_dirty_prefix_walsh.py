"""Discover a nonzero full-space dirty three-macro Walsh discriminator.

Predictions before execution: N7,a1,n_exp1,q3, output C=b0+f. The derived
all-n dirty-label witness guarantees a nonzero chronological coefficient
distinguishing each of reversed order and omitted incoming carry. A common
mask is possible but not guaranteed; record joint or split discovery honestly.
All gate labels must match C89, and P/J symmetries force forbidden-mask zeros.
Frozen discovered masks require a separate independent interval verification.

Only 16384 labels, no size sweep or performance claim. Existing walsh.wht
uses float64 but every butterfly here is exact: integer inputs +/-1 and all
intermediates have magnitude <=16384<2^53. Fractions normalize integer sums.
No new propagator, statevector, or production helper is introduced.

Run: uv run python -u -X faulthandler -m experiments.experiment_dirty_prefix_walsh
"""
from fractions import Fraction
from itertools import product

import numpy as np

from circuits import Circuit
from experiments.experiment_controlled_intervals import macro_image, pack
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp
from walsh import classical_permutation, wht


def parity_sign(value):
    return -1 if int(value).bit_count() & 1 else 1


def run():
    exp = Experiment("dirty_prefix_walsh", doc=__doc__)
    exp.predict("P1", "forward and reverse Circuit maps agree with C89 on every input")
    exp.predict("P2", "both gate maps commute with exact P and J bit flips")
    exp.predict("P3", "exact Walsh sums obey Parseval, balance, and all forbidden-mask zeros")
    exp.predict("P4", "proved dirty label has true endpoint (9,1), reverse and no-h (6,1)")
    exp.predict("P5", "separate nonzero chronological order/carry masks exist; common mask is open")
    exp.predict("P6", "selected transform sums agree with direct integer character sums")
    exp.must_fail("C1", "reversed chronological order disagrees on a selected nonzero coefficient")
    exp.must_fail("C2", "omitting incoming carry disagrees on a selected nonzero coefficient")

    model = ToffoliModExp(7, 1, n_exp=1)
    assert (model.n, model.m, model.n_qubits) == (3, 4, 14)
    size, count = 1 << model.m, 1 << model.n_qubits
    constants = (1, 2, 4)
    schedules = ((0, 1, 2), (2, 1, 0))
    maps, gate_counts = [], []
    for order in schedules:
        circuit = Circuit(model.n_qubits)
        for j in order:
            model.cc_add_mod(circuit, model.exp[0], model.x[j], constants[j])
        assert len(circuit.logical) <= 600
        maps.append(classical_permutation(circuit))
        gate_counts.append(len(circuit.logical))

    no_h = np.empty(count, dtype=np.int64)
    map_errors = [0, 0]
    for t, h, u, x, z in product(range(size), (0, 1), (0, 1), range(8), range(2*size)):
        label = pack(model, t, h, u, x, z)
        for version, order in enumerate(schedules):
            after = z
            for j in order:
                after = macro_image(model.m, model.N, constants[j], t, h, u*((x >> j) & 1), after)
            map_errors[version] += int(maps[version][label]) != pack(model, t, h, u, x, after)
        after = z
        for j in schedules[0]:
            after = macro_image(model.m, model.N, constants[j], t, 0, u*((x >> j) & 1), after)
        # Wrong arithmetic still preserves the physical input carry wire.
        no_h[label] = pack(model, t, h, u, x, after)
    maps.append(no_h)
    exp.check("P1", map_errors == [0, 0], str(map_errors))

    labels = np.arange(count, dtype=np.int64)
    p_mask = (1 << (model.m-1)) | (1 << model.anc)
    j_mask = (size-1) | ((size-1) << model.m) | (1 << model.c0) | (1 << model.anc)
    symmetry_errors = [[int(np.count_nonzero(arr[labels ^ mask] != (arr ^ mask)))
                        for mask in (p_mask, j_mask)] for arr in maps[:2]]
    exp.check("P2", symmetry_errors == [[0, 0], [0, 0]], str(symmetry_errors))

    output_mask = 1 | (1 << model.anc)
    spectra = []
    transform_ok = True
    for arr in maps:
        signs = np.fromiter((parity_sign(y & output_mask) for y in arr), dtype=np.int64, count=count)
        transformed = wht(signs)
        integer = transformed.astype(np.int64)
        transform_ok &= bool(np.array_equal(transformed, integer))
        transform_ok &= sum(int(v)**2 for v in integer) == count**2
        transform_ok &= int(integer[0]) == 0
        spectra.append(integer)
    allowed = np.array([parity_sign((a ^ output_mask) & p_mask) == 1 and
                        parity_sign((a ^ output_mask) & j_mask) == 1 for a in range(count)])
    forbidden_errors = [int(np.count_nonzero(s[~allowed])) for s in spectra[:2]]
    exp.check("P3", transform_ok and int(allowed.sum()) == 4096 and forbidden_errors == [0, 0],
              f"integer/Parseval/balance={transform_ok}; allowed={allowed.sum()}; forbidden={forbidden_errors}")

    witness = pack(model, 2, 1, 1, 7, 0)
    expected = [pack(model, 2, 1, 1, 7, 9+size), pack(model, 2, 1, 1, 7, 6+size)]
    observed = [int(arr[witness]) for arr in maps]
    exp.check("P4", observed == [expected[0], expected[1], expected[1]], str(observed))

    nonzero = (spectra[0] != 0) & allowed
    order_candidates = np.flatnonzero(nonzero & (spectra[0] != spectra[1]))
    carry_candidates = np.flatnonzero(nonzero & (spectra[0] != spectra[2]))
    joint_candidates = np.intersect1d(order_candidates, carry_candidates)
    exp.check("P5", len(order_candidates) > 0 and len(carry_candidates) > 0,
              f"order={len(order_candidates)}; carry={len(carry_candidates)}; common={len(joint_candidates)}")
    selections = {}
    for name, candidates in (("order", order_candidates), ("carry", carry_candidates), ("joint", joint_candidates)):
        if len(candidates):
            selections[name] = int(candidates[0])
    rows, direct_ok = [], True
    for alpha in sorted(set(selections.values())):
        values = [int(s[alpha]) for s in spectra]
        direct = [sum(parity_sign((a & alpha) ^ (int(b) & output_mask))
                      for a, b in enumerate(arr)) for arr in maps]
        direct_ok &= values == direct
        row = dict(input_mask=alpha, output_mask=output_mask,
                   roles=[name for name, mask in selections.items() if mask == alpha],
                   numerator=values, denominator=count,
                   coefficients=[str(Fraction(v, count)) for v in values])
        rows.append(row)
        exp.log(str(row))
    exp.check("P6", direct_ok and bool(rows), str(selections))
    a_order, a_carry = selections.get("order"), selections.get("carry")
    exp.fail_check("C1", a_order is not None and spectra[0][a_order] != 0 and
                   spectra[0][a_order] != spectra[1][a_order], str(selections))
    exp.fail_check("C2", a_carry is not None and spectra[0][a_carry] != 0 and
                   spectra[0][a_carry] != spectra[2][a_carry], str(selections))
    exp.finish(report_path="out/dirty_prefix_walsh_report.json", rows=rows,
               metadata=dict(N=7, a=1, n=3, q=3, n_exp=1, constants=constants,
                             labels=count, gates=gate_counts, map_errors=map_errors,
                             symmetry_errors=symmetry_errors, selections=selections,
                             candidates=dict(order=len(order_candidates), carry=len(carry_candidates),
                                             joint=len(joint_candidates)),
                             explicit_label=dict(input=witness, outputs=observed),
                             branch="joint" if len(joint_candidates) else "split",
                             arithmetic="exact integer butterflies represented in float64; Fraction normalization",
                             scope="mask discovery only; separate frozen-mask interval verification required"))


if __name__ == "__main__":
    run()
