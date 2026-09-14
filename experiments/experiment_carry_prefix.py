"""Do actual Cuccaro logical prefixes satisfy the known carry-correlation certificate?

PREDICTIONS WRITTEN BEFORE MEASURING. Fix width3; vary only logical prefix
length, including partial MAJ/UMA. For zero, single-wire and three fixed
multi-wire observables, both exact carry methods must equal every coefficient
from existing full-space gate replay and Walsh transform. Each wire needs
at most three original/carry symbols. Negative coefficients must occur.
Controls: erasing the carry-in, forgetting signs, or recognizing an altered
gate must fail. This validates a known Wallen baseline, not a speed claim.

Budget: width3 gives 256 basis labels and at most 12*256*20 scalar queries
per backend; one reference vector at a time. No width or memory pilot.
The existing float64 WHT is exact here: its integer intermediates have
absolute value <=256. All final comparisons use exact Fraction values.

Run: uv run python -m experiments.experiment_carry_prefix
"""
from fractions import Fraction

import numpy as np

from circuits import Circuit, ripple_adder
from lab.carry_prefix import compile_ripple_prefix
from lab.harness import Experiment
from walsh import classical_permutation, wht


def run():
    exp = Experiment("carry_prefix", doc=__doc__)
    exp.predict("P1", "both exact backends equal all requested signed FWHT coefficients at every prefix")
    exp.predict("P2", "every prefix compiles with <=3 symbols per wire; partial UMA reaches3")
    exp.predict("P3", "first MAJ majority has full-carry cubic-parity coefficient -1/2")
    exp.must_fail("C1", "clean-carry shortcut gives1 where the full-space coefficient is0")
    exp.must_fail("C2", "magnitude-only extraction loses a negative coefficient")
    exp.must_fail("C3", "changed logical gate must be rejected before using carry rewrites")
    width = 3
    full, regs = ripple_adder(width)
    n, dim = full.n, 1 << full.n
    masks = [0] + [1 << q for q in range(n)] + [3, (1 << n)-1, (1 << regs["a"][1]) | (1 << regs["b"][2])]
    prefix = Circuit(n)
    rows = []
    errors = negatives = comparisons = maximum = 0
    for steps in range(len(full.logical)+1):
        perm = classical_permutation(prefix)
        row_errors = row_negatives = row_max = 0
        for obs in masks:
            cert = compile_ripple_prefix(prefix, obs)
            row_max = max(row_max, cert.max_wire_symbols)
            signs = np.fromiter((1-2*((int(y) & obs).bit_count() & 1) for y in perm), dtype=np.int64)
            integers = wht(signs)
            assert np.array_equal(integers, integers.astype(np.int64))
            for query, numerator in enumerate(integers):
                expected = Fraction(int(numerator), dim)
                actual = cert.coefficient(query)
                reference = cert.coefficient(query, backend="transfer")
                row_errors += int(actual != expected or reference != expected)
                row_negatives += int(expected < 0)
                comparisons += 1
        errors += row_errors
        negatives += row_negatives
        maximum = max(maximum, row_max)
        rows.append(dict(steps=steps, observables=len(masks), queries_per_observable=dim,
                         mismatches=row_errors, negative_coefficients=row_negatives,
                         max_wire_symbols=row_max))
        exp.log(rows[-1])
        if steps < len(full.logical):
            op = full.logical[steps]
            getattr(prefix, op[0])(*op[1:])
    exp.check("P1", errors == 0, f"{comparisons} coefficient triples; mismatches={errors}")
    exp.check("P2", maximum == 3, f"maximum={maximum}")

    first = Circuit(n)
    for op in full.logical[:3]:
        getattr(first, op[0])(*op[1:])
    majority = compile_ripple_prefix(first, 1 << regs["a"][0])
    cubic = (1 << regs["a"][0]) | (1 << regs["b"][0]) | 1
    negative = majority.coefficient(cubic)
    exp.check("P3", negative == Fraction(-1, 2) and negatives > 0,
              f"first majority={negative}; negative reference entries={negatives}")

    obs = 1 << regs["b"][0]
    missing_carry = obs | (1 << regs["a"][0])
    cert = compile_ripple_prefix(full, obs)
    clean_numerator = sum(1-2*(((int(y) & obs).bit_count() ^ (x & missing_carry).bit_count()) & 1)
                          for x, y in enumerate(perm) if not x & 1)
    clean = Fraction(clean_numerator, dim//2)
    unconditioned = cert.coefficient(missing_carry)
    exp.fail_check("C1", clean == 1 and unconditioned == 0 and cert.coefficient(missing_carry | 1) == 1,
                   f"clean={clean}, full-space={unconditioned}, carry included={cert.coefficient(missing_carry | 1)}")
    exp.fail_check("C2", abs(negative) != negative, f"abs={abs(negative)} versus signed={negative}")
    changed = Circuit(n).cnot(regs["a"][0], regs["b"][1])
    rejected = False
    try:
        compile_ripple_prefix(changed, obs)
    except ValueError as exc:
        rejected = "unrecognized Cuccaro logical gate" in str(exc)
    exp.fail_check("C3", rejected, "altered first gate rejected")
    exp.finish(report_path="out/carry_prefix_report.json", rows=rows,
               metadata=dict(width=width, qubits=n, coefficient_triples=comparisons,
                             strongest_baseline="Wallen A84 Sections3.1-3.3, already applied to every prefix",
                             numerical_scope="FWHT integer magnitude<=256, represented exactly by float64",
                             output="selected 2^-N-normalized signed full-space Walsh coefficients"))


if __name__ == "__main__":
    run()
