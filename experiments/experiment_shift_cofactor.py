"""Can final-shift cofactors contract growing dirty prefixes with polynomial memory?

PREDICTIONS BEFORE EXECUTION: the derived construction conditions on the
2^(q+1)-1 signed totals S=(T_lowq XOR X)-T_lowq. Within each affine cell,
one addition carry reveals prefix shifts while a running count/flag state
follows chronological (or inverse chronological) macros. Exact signed sums
and unsigned masses must match actual N7,a1 gate maps at q=1,2,3. Target is
input b_top+h, output b0+f. At q3 forward/reverse are -15/512,-13/512;
omitting arithmetic h makes the coefficient zero. No q2 value is predicted.

This implements the target-specific O(2^q poly(m,q)) integer algorithm with
polynomial working storage, using conservative vertical strips and exact
integer/Fraction arithmetic. It is a correctness prototype, not a measured
time/memory win, general endpoint API, or arbitrary-order construction.
Caps: n<=64,q<=8, 200000 cells and 20000 DP states; run under timeout60s.
No candidate loop enumerates scratch words, b labels or flag histories.

Run: uv run python -u -X faulthandler -m experiments.experiment_shift_cofactor
"""
from fractions import Fraction
from itertools import combinations, product

from circuits import Circuit
from experiments.experiment_controlled_intervals import macro_image, pack
from lab.harness import Experiment
from toffoli_arith import ToffoliModExp
from walsh import classical_permutation


class BudgetExceeded(RuntimeError):
    pass


def sign(bit):
    return -1 if bit & 1 else 1


def stats_new():
    return dict(cells=0, open_cells=0, point_cells=0, peak_lines=0,
                peak_events=0, line_pairs=0, peak_dp_states=0,
                digit_transitions=0, max_integer_T_span=0)


def forms(n, q, total_shift, scratch_top, carry_in, reverse):
    """All possible affine top arguments, indexed by time/count/prefix carry."""
    half, modulus = 1 << n, (1 << n)-1
    delta = scratch_top*half+carry_in
    first, second = [], []
    for j in range(q):
        if reverse:
            i = q-1-j
            after = [total_shift-(total_shift % (1 << i))+v*(1 << i) for v in (0, 1)]
            before = [total_shift-(total_shift % (1 << (i+1)))+v*(1 << (i+1)) for v in (0, 1)]
        else:
            after = [total_shift % (1 << (j+1))-v*(1 << (j+1)) for v in (0, 1)]
            before = [total_shift % (1 << j)-v*(1 << j) for v in (0, 1)]
        first.append([[(j+2*k+2, j*delta-(k+1)*modulus+s) for s in after]
                      for k in range(j+1)])
        second.append([[(j+2*k, j*delta-k*modulus+s) for s in before]
                       for k in range(j+2)])
    return first, second


def cells(n, first, second, stats):
    """Disjoint integer cells with two original-line inequalities and T bounds."""
    half, size = 1 << n, 1 << (n+1)
    lines = {(0, 0), (0, size)}  # b + slope*T = intercept
    for table in (first, second):
        for time in table:
            for choices in time:
                for slope, offset in choices:
                    high = size+slope*half
                    for ell in range(offset//half-1, (offset+high)//half+2):
                        intercept = ell*half-offset
                        if 0 <= intercept <= high:
                            lines.add((slope, intercept))
    lines = sorted(lines)
    stats['peak_lines'] = max(stats['peak_lines'], len(lines))
    if len(lines) > 2000:
        raise BudgetExceeded('more than 2000 affine lines')
    events = {Fraction(0), Fraction(half)}
    for (a, c), (b, d) in combinations(lines, 2):
        stats['line_pairs'] += 1
        if a != b:
            at = Fraction(c-d, a-b)
            if 0 < at < half and 0 <= c-a*at <= size:
                events.add(at)
    events = sorted(events)
    stats['peak_events'] = max(stats['peak_events'], len(events))

    def emit(lo, hi, lower, upper, sample_t, sample_b, kind):
        stats['cells'] += 1
        stats[kind] += 1
        stats['max_integer_T_span'] = max(stats['max_integer_T_span'], hi-lo)
        if stats['cells'] > 200000:
            raise BudgetExceeded('more than 200000 cells')
        return lo, hi, lower, upper, sample_t, sample_b

    for at in events:
        if at.denominator == 1 and at < half:
            t = int(at)
            heights = sorted({c-a*t for a, c in lines if 0 <= c-a*t <= size})
            for low, high in zip(heights, heights[1:]):
                if low < high:
                    yield emit(t, t+1, (0, low), (0, high), at,
                               Fraction(low+high, 2), 'point_cells')
    for left, right in zip(events, events[1:]):
        lo, hi = left.numerator//left.denominator+1, -(-right.numerator//right.denominator)
        if lo >= hi:
            continue
        sample_t = (left+right)/2
        ordered = sorted(lines, key=lambda line: line[1]-line[0]*sample_t)
        start, stop = ordered.index((0, 0)), ordered.index((0, size))
        for lower, upper in zip(ordered[start:stop], ordered[start+1:stop+1]):
            sample_b = (lower[1]+upper[1]-(lower[0]+upper[0])*sample_t)/2
            yield emit(lo, hi, lower, upper, sample_t, sample_b, 'open_cells')


def cell_sum(n, q, total_shift, carry_in, reverse, first, second, cell, stats):
    """Signed target numerator and unsigned physical-input mass in one cell."""
    lo, hi, lower, upper, sample_t, sample_b = cell
    half, width = 1 << n, n+1
    tests = []
    for table in (first, second):
        tests.append([[[int((sample_b+a*sample_t+c)//half) & 1 for a, c in choices]
                       for choices in time] for time in table])
    first_sign, second_sign = tests
    # State: two T borrows, two linear-comparison carries, prefix carry, k, f.
    common = sign((total_shift & 1) ^ (((q+1) & 1)*carry_in))
    if reverse:
        state = {(0, 0, 0, 0, 0, k, f): (common*sign(k+f), 1)
                 for k in range(q+1) for f in (0, 1)}
    else:
        state = {(0, 0, 0, 0, 0, 0, f): (common, 1) for f in (0, 1)}
    columns = max(width, lo.bit_length(), hi.bit_length(), lower[1].bit_length(),
                  upper[1].bit_length())+2
    for i in range(columns):
        after = {}
        for (bl, bh, cl, cu, nu, k, f), (weight, mass) in state.items():
            for b, t in product((0, 1) if i < width else (0,),
                                (0, 1) if i < n else (0,)):
                stats['digit_transitions'] += 1
                next_nu, next_k, next_f = nu, k, f
                if i < q:
                    digit_total = t+((total_shift >> i) & 1)+nu
                    next_nu = digit_total//2
                    if i == q-1 and next_nu != int(total_shift < 0):
                        continue
                    if reverse:
                        j = q-1-i
                        if not 0 <= k <= j+1:
                            continue
                        g = f ^ 1 ^ second_sign[j][k][next_nu]
                        next_k = k-1+g
                        if not 0 <= next_k <= j:
                            continue
                        next_f = g ^ first_sign[j][next_k][nu]
                    else:
                        g = f ^ first_sign[i][k][next_nu]
                        next_k = k+1-g
                        next_f = g ^ 1 ^ second_sign[i][next_k][nu]
                    if i == q-1:
                        next_nu = 0  # All prefix-carry tests are now resolved.
                new = (int(t-((lo >> i) & 1)-bl < 0),
                       int(t-((hi >> i) & 1)-bh < 0),
                       (b+lower[0]*t-((lower[1] >> i) & 1)+cl)//2,
                       (b+upper[0]*t-((upper[1] >> i) & 1)+cu)//2,
                       next_nu, next_k, next_f)
                phase = (b if i in (0, n) else 0) ^ (t*(q & 1) if i == 0 else 0)
                old_weight, old_mass = after.get(new, (0, 0))
                after[new] = (old_weight+sign(phase)*weight, old_mass+mass)
        state = after
        stats['peak_dp_states'] = max(stats['peak_dp_states'], len(state))
        if len(state) > 20000:
            raise BudgetExceeded('more than 20000 DP states')
    result, count = 0, 0
    for (bl, bh, cl, cu, _, k, f), (weight, mass) in state.items():
        if bl == 0 and bh == 1 and cl >= 0 and cu < 0 and (not reverse or k == 0):
            result += weight if reverse else sign(k+f)*weight
            count += mass
    return result, count


def cofactor(n, q, total_shift, reverse, stats):
    if not 2 <= n <= 64 or not 1 <= q <= min(n, 8) or abs(total_shift) >= 1 << q:
        raise ValueError('bounded native n>=2, q<=n and signed difference required')
    numerator = mass = 0
    for p, h in product((0, 1), repeat=2):
        first, second = forms(n, q, total_shift, p, h, reverse)
        for cell in cells(n, first, second, stats):
            value, count = cell_sum(n, q, total_shift, h, reverse, first, second, cell, stats)
            numerator += value
            mass += count
    return numerator, mass


def gate_reference(model, q, reverse):
    circuit = Circuit(model.n_qubits)
    order = reversed(range(q)) if reverse else range(q)
    for j in order:
        model.cc_add_mod(circuit, model.exp[0], model.x[j], 1 << j)
    assert len(circuit.logical) <= 200*q
    images = classical_permutation(circuit)
    sectors = {s: [0, 0] for s in range(-(1 << q)+1, 1 << q)}
    full_sum = 0
    for label, after in enumerate(images):
        phase = ((label >> model.n) & 1) ^ ((label >> model.c0) & 1)
        phase ^= (int(after) & 1) ^ ((int(after) >> model.anc) & 1)
        weight = sign(phase)
        full_sum += weight
        if (label >> model.exp[0]) & 1:
            t = (label >> model.m) & ((1 << q)-1)
            x = (label >> (2*model.m)) & ((1 << q)-1)
            total_shift = (t ^ x)-t
            sectors[total_shift][0] += weight
            sectors[total_shift][1] += 1
    spectators = 1 << (model.n-q)
    assert all(v % spectators == 0 for pair in sectors.values() for v in pair)
    sectors = {s: tuple(v//spectators for v in pair) for s, pair in sectors.items()}
    return sectors, Fraction(full_sum, len(images)), len(circuit.logical)


def run():
    exp = Experiment('shift_cofactor', doc=__doc__)
    exp.predict('P1', 'all signed shift-sector sums match actual forward/reverse gates')
    exp.predict('P2', 'all sector masses equal the independently counted signed-difference fibers')
    exp.predict('P3', 'full common-enable normalization agrees with every gate scalar')
    exp.predict('P4', 'q1 is zero; q3 is -15/512 forward, -13/512 reversed')
    exp.predict('P5', 'both nontrivial open-strip and integer-slice digit contractions are exercised')
    exp.must_fail('C1', 'reversal changes the nonzero q3 scalar')
    exp.must_fail('C2', 'omitting arithmetic h annihilates the nonzero q3 scalar')
    model = ToffoliModExp(7, 1, n_exp=1)
    stats, rows = stats_new(), []
    sum_errors = mass_errors = full_errors = 0
    values = {}
    try:
        for q, reverse in product((1, 2, 3), (False, True)):
            actual, full, gates = gate_reference(model, q, reverse)
            candidate = {}
            for total_shift in range(-(1 << q)+1, 1 << q):
                pair = cofactor(3, q, total_shift, reverse, stats)
                candidate[total_shift] = pair
                expected_mass = (1 << (2*model.m-q+2))*((1 << q)-abs(total_shift))
                sum_errors += pair[0] != actual[total_shift][0]
                mass_errors += pair[1] != actual[total_shift][1] or pair[1] != expected_mass
            value = Fraction((1 << q)*candidate[0][0]+sum(v[0] for v in candidate.values()),
                             1 << (2*model.m+q+3))
            values[q, reverse] = value
            full_errors += value != full
            row = dict(q=q, reverse=reverse, logical_gates=gates, coefficient=str(value),
                       gate_coefficient=str(full),
                       sectors=[dict(S=s, candidate=candidate[s], actual=actual[s]) for s in candidate])
            rows.append(row)
            exp.log(f'q={q}, reverse={reverse}, value={value}, gates={full}, cells={stats["cells"]}')
    except BudgetExceeded as error:
        exp.check('P1', False, f'construction stopped at declared cap: {error}')
        exp.finish(report_path='out/shift_cofactor_report.json', rows=rows,
                   metadata=dict(stats=stats, incomplete=True, reason=str(error)))
        return
    exp.check('P1', sum_errors == 0, f'50 exact signed sectors; errors={sum_errors}')
    exp.check('P2', mass_errors == 0, f'50 exact unsigned masses; errors={mass_errors}')
    exp.check('P3', full_errors == 0, f'six full-space scalars; errors={full_errors}')
    exp.check('P4', values[1, False] == values[1, True] == 0 and
              values[3, False] == Fraction(-15, 512) and values[3, True] == Fraction(-13, 512),
              str(values))
    exp.check('P5', stats['open_cells'] > 0 and stats['point_cells'] > 0 and
              stats['max_integer_T_span'] > 1, str(stats))
    exp.fail_check('C1', values[3, False] != 0 and values[3, False] != values[3, True], str(values))
    wrong_sum = 0
    size = 1 << model.m
    for t, h, u, x, z in product(range(size), (0, 1), (0, 1), range(8), range(2*size)):
        after = z
        for j in range(3):
            after = macro_image(model.m, model.N, 1 << j, t, 0, u*((x >> j) & 1), after)
        wrong_sum += sign(((z >> model.n) & 1) ^ h ^ (after & 1) ^ (after >> model.m))
    wrong = Fraction(wrong_sum, 1 << model.n_qubits)
    exp.fail_check('C2', wrong == 0 and values[3, False] != wrong, f'wrong={wrong}, true={values[3, False]}')
    exp.finish(report_path='out/shift_cofactor_report.json', rows=rows,
               metadata=dict(stats=stats, incomplete=False, arithmetic='exact integers/Fractions',
                             wrong_no_h=str(wrong), scope='bounded algorithm correctness; no resource benchmark',
                             target='input btop+h; output b0+f; full dirty-space average'))


if __name__ == '__main__':
    run()
