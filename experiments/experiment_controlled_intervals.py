"""Do shared-state controlled-add prefixes admit exact interval contraction?

PREDICTIONS BEFORE EXECUTION: conditioned cc_add_mod is an exchange of at most
16 translated intervals on z=f*2^width+b. Bijective composition adds cuts, so
q macros have <=1+15q pieces. A four-state range Walsh DP gives each exact
conditional coefficient. Actual N5,a2 prefixes q1,2,3 must agree on all full-
space labels and two fixed Walsh probes. At q1 their derived values are1/4
and-1/8. Other prefix values are not predicted. Wrong incoming-carry,
disabled-control and reordered shared-state references must disagree.

Known interval-exchange/digit-DP application. Conditional O(q²) construction,
O(q*width) exact contraction, with growing integer bit sizes. Full-space outer
average is exponential in unchanged scratch/control variables. No simulation
breakthrough, general matchgate extension or memory advantage is asserted.
Fixed14-qubit tiny family; 3*16384 actual gate labels. Separate 4-bit range-DP
exhaustive test. No statevector/unitary or long native permutation job.

Run: uv run python -u -X faulthandler -m experiments.experiment_controlled_intervals
"""
from fractions import Fraction
from itertools import product

from circuits import Circuit
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp
from walsh import classical_permutation


def macro_image(width, modulus, constant, t, h, control, z):
    """C86 word identity, not gate replay; all scratch/flag inputs are retained."""
    size = 1 << width
    b, f = z % size, z // size
    delta_n = modulus - 2*(t & modulus)
    a = control*(constant - 2*(t & constant))
    g = f ^ (((b+a-delta_n) % size) >> (width-1))
    after_b = (b+t+h+a-(1-g)*delta_n) % size
    after_f = g ^ 1 ^ (((b-(1-g)*delta_n) % size) >> (width-1))
    return after_b + size*after_f


def _merge(pieces):
    answer = []
    for lo, hi, shift in pieces:
        if answer and answer[-1][1] == lo and answer[-1][2] == shift:
            answer[-1] = (answer[-1][0], hi, shift)
        else:
            answer.append((lo, hi, shift))
    return answer


def _validate_exchange(pieces, size):
    cursor = 0
    for lo, hi, shift in pieces:
        if lo != cursor or not lo < hi or not 0 <= lo+shift < hi+shift <= size:
            raise AssertionError("invalid interval partition or output range")
        cursor = hi
    if cursor != size:
        raise AssertionError("incomplete interval domain")
    cursor = 0
    for lo, hi in sorted((lo+s, hi+s) for lo, hi, s in pieces):
        if lo != cursor:
            raise AssertionError("interval images overlap or have gaps")
        cursor = hi
    if cursor != size:
        raise AssertionError("incomplete interval image")


def macro_intervals(width, modulus, constant, t, h, control):
    size, half = 1 << width, 1 << (width-1)
    dn = modulus-2*(t & modulus)
    a, d = control*(constant-2*(t & constant)), t+h
    cuts = {0, size}
    cuts.update(v % size for v in (dn-a, half+dn-a, dn, half+dn,
                                  half, -d-a+dn, -d-a))
    cuts = sorted(cuts)
    pieces = []
    for flag in (0, 1):
        for lo, hi in zip(cuts, cuts[1:]):
            start = flag*size+lo
            image = macro_image(width, modulus, constant, t, h, control, start)
            pieces.append((start, flag*size+hi, image-start))
    pieces = _merge(pieces)
    _validate_exchange(pieces, 2*size)
    if len(pieces) > 16:
        raise AssertionError("derived one-macro bound violated")
    return pieces


def compose_intervals(first, second):
    """Return second after first; two full interval tables, no history tree."""
    result = []
    for lo, hi, shift in first:
        for after_lo, after_hi, after_shift in second:
            left, right = max(lo+shift, after_lo), min(hi+shift, after_hi)
            if left < right:
                result.append((left-shift, right-shift, shift+after_shift))
    return _merge(result)


def compile_intervals(width, modulus, constants, t, h, controls):
    """Experimental conditional representation, not a production circuit recognizer."""
    if not 1 <= width <= 4096 or len(constants) != len(controls) or len(constants) > 256:
        raise ValueError("bounded equal-length constant/control lists required")
    size = 1 << width
    if not 0 < modulus < size or not 0 <= t < size or h not in (0, 1):
        raise ValueError("invalid fixed scratch/modulus")
    if any(not 0 <= c < size for c in constants) or any(e not in (0, 1) for e in controls):
        raise ValueError("invalid constant or classical control")
    pieces = [(0, 2*size, 0)]
    for step, (constant, control) in enumerate(zip(constants, controls), 1):
        pieces = compose_intervals(pieces, macro_intervals(width, modulus, constant, t, h, control))
        if len(pieces) > 1+15*step:
            raise AssertionError("derived composition bound violated")
    _validate_exchange(pieces, 2*size)
    return pieces


def prefix_walsh(bits, limit, shift, alpha, beta):
    """Sum x<limit of chi_alpha(x) chi_beta(x+shift mod2^bits)."""
    size = 1 << bits
    if not 0 <= limit <= size or not 0 <= alpha < size or not 0 <= beta < size:
        raise ValueError("threshold/mask out of range")
    if not limit:
        return 0
    shift %= size
    state = {(0, 0): 1}  # addition carry, subtraction borrow for x-limit
    for i in range(bits):
        row = {}
        for (carry, borrow), weight in state.items():
            for x in (0, 1):
                total = x+((shift >> i) & 1)+carry
                y = total & 1
                after = (total >> 1, int(x-((limit >> i) & 1)-borrow < 0))
                phase = (((alpha >> i) & 1)*x) ^ (((beta >> i) & 1)*y)
                row[after] = row.get(after, 0) + (-weight if phase else weight)
        state = row
    return sum(weight for (_, borrow), weight in state.items()
               if limit == size or borrow == 1)


def interval_coefficient(width, pieces, input_mask, output_mask):
    bits = width+1
    numerator = sum(prefix_walsh(bits, hi, shift, input_mask, output_mask)
                    - prefix_walsh(bits, lo, shift, input_mask, output_mask)
                    for lo, hi, shift in pieces)
    return Fraction(numerator, 1 << bits)


def interval_image(pieces, z):
    return next(z+s for lo, hi, s in pieces if lo <= z < hi)


def actual_prefix(model, length):
    circuit = Circuit(model.n_qubits)
    constants = [(model.a*(1 << i)) % model.N for i in range(length)]
    for i, constant in enumerate(constants):
        model.cc_add_mod(circuit, model.exp[0], model.x[i], constant)
    return constants, circuit


def pack(model, t, h, u, x, z):
    size = 1 << model.m
    return ((z % size) | (t << model.m) | (x << (2*model.m)) |
            (h << model.c0) | ((z // size) << model.anc) | (u << model.exp[0]))


def run():
    exp = Experiment("controlled_intervals", doc=__doc__)
    exp.predict("P1", "all actual14-qubit prefix maps agree with interval composition")
    exp.predict("P2", "all4-bit translations/mask pairs/thresholds agree with direct signed sums")
    exp.predict("P3", "six full-space probes agree between outer interval average and gate replay")
    exp.predict("P4", "independently derived q1 probe values are1/4 and-1/8")
    exp.predict("P5", "all conditional prefixes obey the linear interval count bound")
    exp.must_fail("C1", "omitting the incoming carry annihilates both nonzero q1 probes")
    exp.must_fail("C2", "fixing the active conjunction to zero annihilates the negative q1 probe")
    exp.must_fail("C3", "reordering two actual shared-state macros changes the output flag")

    dp_errors = dp_count = 0
    for shift, alpha, beta in product(range(-8, 8), range(16), range(16)):
        expected = 0
        for limit in range(17):
            if limit:
                x = limit-1
                phase = (x & alpha).bit_count() + (((x+shift) % 16) & beta).bit_count()
                expected += -1 if phase % 2 else 1
            dp_errors += prefix_walsh(4, limit, shift, alpha, beta) != expected
            dp_count += 1
    exp.check("P2", dp_errors == 0, f"{dp_count} prefix sums; errors={dp_errors}")

    model = ToffoliModExp(5, 2, n_exp=1)
    assert (model.m, model.n_qubits) == (4, 14)
    size = 1 << model.m
    # Both input probes query b0,b3,h,f; their scratch masks differ.
    z_mask = 1 | (1 << (model.m-1)) | (1 << model.m)
    t_masks = (4, 7)
    input_masks = [1 | (1 << (model.m-1)) | (tm << model.m) |
                   (1 << model.c0) | (1 << model.anc) for tm in t_masks]
    rows = []
    all_map_errors = all_coefficient_errors = 0
    bound_ok = True
    wrong_h, wrong_e = [Fraction(0), Fraction(0)], Fraction(0)
    shared_witness = None
    for length in (1, 2, 3):
        constants, circuit = actual_prefix(model, length)
        assert len(circuit.logical) <= 200*length
        actual = classical_permutation(circuit)
        sums = [Fraction(0), Fraction(0)]
        map_errors = maximum_pieces = 0
        fibers = 0
        for t, h, u, x in product(range(size), (0, 1), (0, 1), range(1 << model.n)):
            controls = [u*((x >> i) & 1) for i in range(length)]
            pieces = compile_intervals(model.m, model.N, constants, t, h, controls)
            maximum_pieces = max(maximum_pieces, len(pieces))
            fibers += 1
            for z in range(2*size):
                label = pack(model, t, h, u, x, z)
                expected = pack(model, t, h, u, x, interval_image(pieces, z))
                map_errors += int(actual[label]) != expected
            value = interval_coefficient(model.m, pieces, z_mask, 1)
            for j, tm in enumerate(t_masks):
                sign = -1 if ((t & tm).bit_count()+h) % 2 else 1
                sums[j] += sign*value
            if length == 1:
                no_h = compile_intervals(model.m, model.N, constants, t, 0, controls)
                no_e = compile_intervals(model.m, model.N, constants, t, h, [0])
                h_value = interval_coefficient(model.m, no_h, z_mask, 1)
                e_value = interval_coefficient(model.m, no_e, z_mask, 1)
                for j, tm in enumerate(t_masks):
                    sign = -1 if ((t & tm).bit_count()+h) % 2 else 1
                    wrong_h[j] += sign*h_value
                wrong_e += (-1 if ((t & t_masks[1]).bit_count()+h) % 2 else 1)*e_value
        coefficients = [value/fibers for value in sums]
        expected = [Fraction(sum(-1 if ((label & mask).bit_count() + (int(after) & 1)) % 2 else 1
                                 for label, after in enumerate(actual)), len(actual))
                    for mask in input_masks]
        errors = sum(a != b for a, b in zip(coefficients, expected))
        all_map_errors += map_errors
        all_coefficient_errors += errors
        bound_ok &= maximum_pieces <= 1+15*length
        if length == 1:
            wrong_h = [value/fibers for value in wrong_h]
            wrong_e /= fibers
        if length == 2:
            label = pack(model, 0, 1, 1, 3, 0)
            observed = int(actual[label])
            bad = compile_intervals(model.m, model.N, list(reversed(constants)), 0, 1, [1, 1])
            incorrect = pack(model, 0, 1, 1, 3, interval_image(bad, 0))
            shared_witness = dict(input=label, actual=observed, reordered=incorrect,
                                  expected=pack(model, 0, 1, 1, 3, 3),
                                  expected_reordered=pack(model, 0, 1, 1, 3, size+3))
        row = dict(prefix=length, constants=constants, logical_gates=len(circuit.logical),
                   inputs=len(actual), fibers=fibers, max_pieces=maximum_pieces,
                   map_errors=map_errors, coefficients=[str(x) for x in coefficients],
                   expected=[str(x) for x in expected])
        rows.append(row)
        exp.log(str(row))
    exp.check("P1", all_map_errors == 0, f"49152 full-space labels; errors={all_map_errors}")
    exp.check("P3", all_coefficient_errors == 0, f"six exact full-space probes; errors={all_coefficient_errors}")
    exp.check("P4", rows[0]["coefficients"] == ["1/4", "-1/8"], str(rows[0]))
    exp.check("P5", bound_ok, str([r["max_pieces"] for r in rows]))
    exp.fail_check("C1", wrong_h == [0, 0], f"wrong={wrong_h}; correct={rows[0]['coefficients']}")
    exp.fail_check("C2", wrong_e == 0, f"wrong={wrong_e}; correct={rows[0]['coefficients'][1]}")
    exp.fail_check("C3", shared_witness["actual"] == shared_witness["expected"] and
                   shared_witness["reordered"] == shared_witness["expected_reordered"] and
                   shared_witness["actual"] != shared_witness["reordered"], str(shared_witness))
    exp.finish(report_path="out/controlled_intervals_report.json", rows=rows,
               metadata=dict(dp_checks=dp_count, exact="integer/Fraction", shared_witness=shared_witness,
                             conditioning="all t,h,u,x values explicitly averaged; no hidden free oracle",
                             scope="bounded known interval-exchange application; no memory/speed/novelty claim"))


if __name__ == "__main__":
    run()
