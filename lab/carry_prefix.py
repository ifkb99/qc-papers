"""Exact selected Walsh queries for recognized Cuccaro logical prefixes.

Known carry-correlation machinery (Wallen, HUT-TCS-A84, Sections3.1--3.3),
applied to the actual intermediate wires. This is a restricted certificate,
not a new general simulator. All input bits, including c0 and z, are uniform
in the full-space Walsh sum. Circuit.logical is the authoritative gate trace;
prefixes inside its compiled Pauli rotations are outside this API.
"""
from dataclasses import dataclass
from fractions import Fraction
from operator import index


def _schedule(width):
    for i in range(width):
        a, b, c = 1+i, 1+width+i, 0 if i == 0 else i
        yield ("cnot", a, b), None
        yield ("cnot", a, c), None
        yield ("toffoli", c, b, a), ("carry", i+1)
    yield ("cnot", width, 2*width+1), None
    for i in reversed(range(width)):
        a, b, c = 1+i, 1+width+i, 0 if i == 0 else i
        yield ("toffoli", c, b, a), ("input", a)
        yield ("cnot", a, c), None
        yield ("cnot", c, b), None


def _wallen_clean(width, alpha, beta, gamma):
    """Known MSB-first automaton; gamma includes c_width via a dummy top bit."""
    state, negative, exponent = 0, False, 0
    for i in range(width, -1, -1):
        letter = (((gamma >> i) & 1) << 2) | (((alpha >> i) & 1) << 1) | ((beta >> i) & 1)
        if state == 0:
            if letter == 4:
                state = 1
            elif letter != 0:
                return Fraction(0)
        else:
            exponent += 1
            negative ^= letter in (3, 7)
            state = int(letter in (0, 3, 5, 6))
    return Fraction(-1 if negative else 1, 1 << exponent)


def _carry_dp(width, alpha, beta, gamma, delta0):
    """Independent simple integer transfer reference, including both c0 values."""
    current = [1, -1 if delta0 else 1]
    for i in range(width):
        next_row = [0, 0]
        for carry in (0, 1):
            for a in (0, 1):
                for b in (0, 1):
                    after = (a+b+carry) >> 1
                    phase = (((alpha >> i) & 1)*a) ^ (((beta >> i) & 1)*b) ^ (((gamma >> (i+1)) & 1)*after)
                    next_row[after] += (-1 if phase else 1)*current[carry]
        current = next_row
    return Fraction(sum(current), 1 << (2*width+1))


@dataclass(frozen=True)
class CarryPrefix:
    """Certificate returned by compile_ripple_prefix; construct through that API."""
    width: int
    steps: int
    input_mask: int
    carry_mask: int
    max_wire_symbols: int

    def coefficient(self, query, *, backend="wallen"):
        """One exact 2^-N-normalized coefficient; no free cut/search/oracle."""
        n = 2*self.width+2
        query = index(query)
        if query < 0 or query.bit_length() > n:
            raise ValueError("query outside circuit")
        if backend not in ("wallen", "transfer"):
            raise ValueError("unknown carry backend")
        linear = self.input_mask ^ query
        if (linear >> (n-1)) & 1:
            return Fraction(0)
        wordmask = (1 << self.width)-1
        alpha = (linear >> 1) & wordmask
        beta = (linear >> (self.width+1)) & wordmask
        delta0 = linear & 1
        if backend == "transfer":
            return _carry_dp(self.width, alpha, beta, self.carry_mask, delta0)
        # Joint complement maps carry-in1 to carry-in0 and complements every
        # later carry. Full c0 average is either clean correlation or zero.
        parity = (alpha.bit_count()+beta.bit_count()+self.carry_mask.bit_count()+delta0) & 1
        if parity:
            return Fraction(0)
        return _wallen_clean(self.width, alpha, beta, self.carry_mask)


def compile_ripple_prefix(circuit, observable, *, max_width=4096):
    """Recognize the exact logical schedule and compile one observed parity.

    Each wire has <=3 symbolic terms, original bits or carry bits. Matching
    MAJ/UMA Toffolis use the proved carry identity; unrelated circuits raise.
    Compilation stores O(width) small symbol sets with O(log width)-bit IDs,
    then returns two O(width)-bit phase masks. Supplied Circuit storage is
    additional. Query extraction/automaton make O(width) bit visits; Python
    arbitrary-width shift/bit-operation costs and Fraction output are extra.
    No exponential term list, branch partition or truth table is constructed.
    """
    n, observable = index(circuit.n), index(observable)
    if n < 4 or n % 2 or not 1 <= (n-2)//2 <= index(max_width):
        raise ValueError("expected bounded Cuccaro register layout")
    if observable < 0 or observable.bit_length() > n:
        raise ValueError("observable outside circuit")
    width = (n-2)//2
    if len(circuit.logical) > 6*width+1:
        raise ValueError("trace exceeds Cuccaro schedule")
    wires = [frozenset((i,)) for i in range(n)]
    maximum = 1
    schedule = _schedule(width)
    for step, op in enumerate(circuit.logical, 1):
        expected, rewrite = next(schedule)
        if op != expected:
            raise ValueError(f"unrecognized Cuccaro logical gate at step {step}")
        if rewrite is None:
            _, control, target = op
            wires[target] = wires[target] ^ wires[control]
        else:
            kind, symbol = rewrite
            wires[op[3]] = frozenset((n+symbol-1 if kind == "carry" else symbol,))
        maximum = max(maximum, len(wires[op[-1]]))
        if maximum > 3:
            raise AssertionError("proved constant wire-expression support violated")
    input_mask = carry_mask = 0
    for q in range(n):
        if (observable >> q) & 1:
            for symbol in wires[q]:
                if symbol < n:
                    input_mask ^= 1 << symbol
                else:
                    carry_mask ^= 1 << (symbol-n+1)
    return CarryPrefix(width, len(circuit.logical), input_mask, carry_mask, maximum)
