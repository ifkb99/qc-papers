"""Can two actual modular macros be averaged without a scratch truth table?

PREDICTIONS BEFORE EXECUTION: for the b0-to-b0 Walsh family, two actual
dirty macros reduce to a top-bit autocorrelation at D+s1+s2-r. A signed
carry DP with states in [-3,4] sums all scratch bits and incoming carry;
eight common-control/multiplicand assignments complete the full scalar.
Native N3,a2 has Q(b0+t2+h+x0 ; b0)=-1/16. All128 outside characters
must match actual gate replay. Missing h makes that probe zero. This family
is order independent, so the separate N3,a1,t2,h0,u1 conditional probe
must distinguish chronological/reversed macros: -1/8 versus +1/8.

Fixed correctness fixtures only. The general proof is fixed TWO macros,
one specified endpoint-mask family. No peak-memory, speed, novelty, growing
prefix or output-sampling claim. Negative-carry arithmetic is checked at
N9,c=(3,6),m5 by explicit integer sums, without a large circuit truth table.

Run: uv run python -u -X faulthandler -m experiments.experiment_shared_control_average
"""
from fractions import Fraction
from itertools import product

from circuits import Circuit
from experiments.experiment_controlled_intervals import macro_image, pack
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp
from walsh import classical_permutation


def sign(value):
    return -1 if value.bit_count() & 1 else 1


def enabled_average(width, modulus, constants, enables, t_mask, h_mask):
    """Exact E_(t,h,b,f) chi_t chi_h chi_b0 chi_b0_after for two macros."""
    if not 2 <= width <= 4096 or len(constants) != 2 or len(enables) != 2:
        raise ValueError("bounded width and exactly two macros required")
    size = 1 << width
    if not 0 < modulus < size or not modulus & 1:
        raise ValueError("positive odd modulus below word range required")
    if any(not 0 <= c < modulus for c in constants):
        raise ValueError("constant outside modulus")
    if any(e not in (0, 1) for e in enables) or h_mask not in (0, 1):
        raise ValueError("binary enables/carry mask required")
    if not 0 <= t_mask < size:
        raise ValueError("scratch mask outside word")
    c1, c2 = constants
    e1, e2 = enables
    kappa = (e1*c1 + e2*c2 - modulus) % size
    state = {0: 1, 1: sign(h_mask)}
    for j in range(width):
        coefficient = (1 + 2*((modulus >> j) & 1)
                       - 2*e1*((c1 >> j) & 1) - 2*e2*((c2 >> j) & 1))
        row = {}
        for carry, weight in state.items():
            for w, t in product((0, 1), repeat=2):
                total = w + coefficient*t + ((kappa >> j) & 1) + carry
                next_carry = total // 2  # Floor division is required for negatives.
                if not -3 <= next_carry <= 4:
                    raise AssertionError("derived carry bound violated")
                phase = ((t_mask >> j) & 1)*t
                if j == width-1:
                    phase ^= w ^ (total & 1)
                row[next_carry] = row.get(next_carry, 0) + sign(phase)*weight
        state = row
    prefactor = -sign((e1*(c1 & 1)) ^ (e2*(c2 & 1)))
    return Fraction(prefactor*sum(state.values()), 1 << (2*width+1))


def shared_average(width, modulus, constants, t_mask, h_mask, x_mask, u_mask):
    if not 0 <= x_mask < 4 or u_mask not in (0, 1):
        raise ValueError("two multiplicand bits and one common control required")
    numerator = Fraction(0)
    for u, x0, x1 in product((0, 1), repeat=3):
        x = x0 + 2*x1
        numerator += sign((x & x_mask).bit_count() % 2 ^ (u*u_mask)) * enabled_average(
            width, modulus, constants, (u*x0, u*x1), t_mask, h_mask)
    return numerator / 8


def two_macros(model, reverse=False):
    constants = (model.a, 2*model.a % model.N)
    circuit = Circuit(model.n_qubits)
    indices = (1, 0) if reverse else (0, 1)
    for j in indices:
        model.cc_add_mod(circuit, model.exp[0], model.x[j], constants[j])
    return constants, circuit


def direct_coefficient(images, alpha, beta):
    return Fraction(sum(sign((x & alpha).bit_count() % 2 ^
                             (int(y) & beta).bit_count() % 2)
                        for x, y in enumerate(images)), len(images))


def triangular_ac(width, shift):
    size = 1 << width
    residue = shift % size
    return 1 - Fraction(4*min(residue, size-residue), size)


def run():
    exp = Experiment("shared_control_average", doc=__doc__)
    exp.predict("P1", "all128 full-space outside-character coefficients match actual N3a2 gates")
    exp.predict("P2", "nonzero carry-sensitive full coefficient is -1/16")
    exp.predict("P3", "two-macro conditional b0 correlation equals triangular autocorrelation")
    exp.predict("P4", "signed carry DP includes exact negative-carry arithmetic")
    exp.predict("P5", "conditional chronology probe is -1/8 forward and +1/8 reversed")
    exp.must_fail("C1", "fixing incoming h to zero annihilates the full nonzero probe")
    exp.must_fail("C2", "disabled-as-identity gives zero for the chronological conditional probe")
    exp.predict("P6", "the b0 endpoint family itself is order independent")

    model = ToffoliModExp(3, 2, n_exp=1)
    constants, circuit = two_macros(model)
    assert model.n_qubits == 11 and len(circuit.logical) == 216
    actual = classical_permutation(circuit)
    errors = 0
    histogram = {}
    for tm, hm, xm, um in product(range(8), range(2), range(4), range(2)):
        alpha = (1 | (tm << model.m) | (hm << model.c0) |
                 (xm << (2*model.m)) | (um << model.exp[0]))
        predicted = shared_average(model.m, model.N, constants, tm, hm, xm, um)
        expected = direct_coefficient(actual, alpha, 1)
        errors += predicted != expected
        histogram[str(expected)] = histogram.get(str(expected), 0) + 1
    exp.check("P1", errors == 0, f"128 full-space coefficients; errors={errors}")
    probe = shared_average(3, 3, constants, 4, 1, 1, 0)
    exp.check("P2", probe == Fraction(-1, 16), str(probe))

    conditional_errors = 0
    order_errors = 0
    wrong_h_sum = 0
    for t, h, e1, e2 in product(range(8), range(2), range(2), range(2)):
        r = 3-2*(t & 3)
        s1, s2 = (e1*(constants[0]-2*(t & constants[0])),
                  e2*(constants[1]-2*(t & constants[1])))
        formula = -sign(e2)*triangular_ac(3, t+h+s1+s2-r)
        total = reversed_total = 0
        for z in range(16):
            label = pack(model, t, h, 1, e1+2*e2, z)
            total += sign((z & 1) ^ (int(actual[label]) & 1))
            after = macro_image(3, 3, constants[1], t, h, e2, z)
            after = macro_image(3, 3, constants[0], t, h, e1, after)
            reversed_total += sign((z & 1) ^ (after & 1))
        conditional_errors += Fraction(total, 16) != formula
        order_errors += total != reversed_total
    exp.check("P3", conditional_errors == 0, f"64 conditional fibers; errors={conditional_errors}")
    exp.check("P6", order_errors == 0, f"64 reversed word-map correlations; errors={order_errors}")
    for t, h, u, x, z in product(range(8), range(2), range(2), range(4), range(16)):
        after = z
        for j, c in enumerate(constants):
            after = macro_image(3, 3, c, t, 0, u*((x >> j) & 1), after)
        phase = (z & 1) ^ (after & 1) ^ ((t >> 2) & 1) ^ h ^ (x & 1)
        wrong_h_sum += sign(phase)
    wrong_h = Fraction(wrong_h_sum, 2048)
    exp.fail_check("C1", wrong_h == 0 and wrong_h != probe, f"wrong={wrong_h}; actual={probe}")

    negative_errors = 0
    for tm, hm, e1, e2 in product((0, 16, 31), range(2), range(2), range(2)):
        total = 0
        for t, h, w in product(range(32), range(2), range(32)):
            shift = t+h + e1*(3-2*(t & 3)) + e2*(6-2*(t & 6)) - (9-2*(t & 9))
            phase = ((t & tm).bit_count() % 2) ^ (h*hm) ^ (w >> 4) ^ (((w+shift) % 32) >> 4)
            total += sign(phase)
        expected = Fraction(-sign(e1)*total, 2048)
        negative_errors += enabled_average(5, 9, (3, 6), (e1, e2), tm, hm) != expected
    exp.check("P4", negative_errors == 0, f"24 signed sums with negative carry coefficients; errors={negative_errors}")

    witness_model = ToffoliModExp(3, 1, n_exp=1)
    cs, forward_circuit = two_macros(witness_model)
    _, reverse_circuit = two_macros(witness_model, reverse=True)
    forward, reverse = classical_permutation(forward_circuit), classical_permutation(reverse_circuit)
    sums = [0, 0, 0]
    for x, z in product(range(4), range(16)):
        label = pack(witness_model, 2, 0, 1, x, z)
        wrong = z
        for j, c in enumerate(cs):
            if (x >> j) & 1:
                wrong = macro_image(3, 3, c, 2, 0, 1, wrong)
        phase = (z & 13).bit_count() % 2 ^ x.bit_count() % 2
        for j, output in enumerate((int(forward[label]), int(reverse[label]), wrong)):
            sums[j] += sign(phase ^ (output & 1))
    values = [Fraction(v, 64) for v in sums]
    exp.check("P5", values[:2] == [Fraction(-1, 8), Fraction(1, 8)], str(values))
    exp.fail_check("C2", values[2] == 0 and values[2] != values[0], str(values))
    rows = [{"full_probe": str(probe), "wrong_h": str(wrong_h),
                 "chronology_forward_reverse_disabled_identity": list(map(str, values)),
                 "full_coefficient_histogram": histogram}]
    metadata = {"scope": "two-macro b0 endpoint family; conditional separate chronology witness",
                    "exact": "integer/Fraction", "outside_coefficients": 128,
                    "peak_bytes": "not measured", "negative_carry_checks": 24}
    exp.finish(report_path="out/shared_control_average_report.json", rows=rows, metadata=metadata)


if __name__ == "__main__":
    run()
