"""Verified exact-input specialization of C59's finite prefix contraction.

Optional dependency: python-flint==0.9.0 (Arb/Acb). Inputs are rational
multiples of pi, NOT rounded matrices or angles. b=2 (default) or opt-in b=3;
known indexed period divisible by b; width<=63; at most eight coherent
reflections. Backgrounds rotate fine labels 0/1 and fix the rest. No orbit table.

The adaptive oracle and integer kernel implement the C60 error contract.
Certification trusts FLINT enclosures, exact Python integers, and independent
unbiased bits from rng.getrandbits. A seeded PRNG is only a test driver.
The mathematical sampler has no precision/attempt cap. A diagnostic oracle
cap raises; conditioning on successful capped runs is NOT certified.

Do not call concurrently in threads sharing flint.ctx. Working precision is
restored, but changing the process-wide context is not a thread-local lock.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from types import MappingProxyType

from lab.coherent_routes import _contract_route_component
from lab.sampling_error import plan_prefix_accuracy, prefix_error_budget


def _bits(bits):
    if type(bits) is not int or not 0 <= bits <= 4096:
        raise ValueError("require integer bits in [0,4096]")


def round_dyadic(value, bits):
    """Nearest integer to value*2^bits, with exact ties to even."""
    _bits(bits)
    if type(value) not in (int, Fraction):
        raise ValueError("round an exact int/Fraction, never a float")
    scaled = Fraction(value) * (1 << bits)
    floor, remainder = divmod(scaled.numerator, scaled.denominator)
    twice = 2 * remainder
    return floor + (twice > scaled.denominator or
                    twice == scaled.denominator and floor % 2 != 0)


def _weights(weights):
    if (not isinstance(weights, (tuple, list)) or not 1 <= len(weights) <= 64
            or any(type(w) is not int or w < 0 for w in weights)):
        raise ValueError("require 1..64 nonnegative integer weights")
    # This is a declared normalized fallback, never discarded probability.
    return tuple(weights) if any(weights) else (1,) * len(weights)


def integer_cdf_index(weights, random_word, random_bits):
    """Exact inverse CDF at word/2^L; all-zero weights mean uniform bins."""
    _bits(random_bits)
    scale = 1 << random_bits
    if type(random_word) is not int or not 0 <= random_word < scale:
        raise ValueError("random word outside its stated bit range")
    weights = _weights(weights)
    cutoff, cumulative = random_word * sum(weights), 0
    for i, w in enumerate(weights):
        cumulative += w
        if cutoff < cumulative * scale:
            return i
    raise AssertionError("exact finite CDF must cover every input word")


def integer_cdf_counts(weights, random_bits):
    """Exact bin counts of this kernel, without enumerating its random words."""
    _bits(random_bits)
    weights = _weights(weights)
    total, scale, cumulative, previous = sum(weights), 1 << random_bits, 0, 0
    counts = []
    for weight in weights:
        cumulative += weight
        endpoint = (scale * cumulative + total - 1) // total
        counts.append(endpoint - previous)
        previous = endpoint
    return tuple(counts)


def _word(rng, bits):
    value = rng.getrandbits(bits)
    if type(value) is not int or not 0 <= value < 1 << bits:
        raise ValueError("getrandbits returned an invalid word")
    return value


def uniform_integer(bound, rng):
    """Exact uniform integer under independent unbiased input bits; no cap."""
    if type(bound) is not int or bound < 1:
        raise ValueError("bound must be a positive integer")
    if bound == 1:
        return 0
    bits = (bound - 1).bit_length()
    while True:
        value = _word(rng, bits)
        if value < bound:
            return value


def binary_fraction(value):
    """Lossless extraction of an exact finite Arb value, including negatives."""
    if not value.is_finite() or not value.is_exact():
        raise ValueError("binary extraction requires an exact finite Arb value")
    mantissa, exponent = (int(v) for v in value.man_exp())
    return (Fraction(mantissa << exponent) if exponent >= 0
            else Fraction(mantissa, 1 << -exponent))


def _rational_phase(numerator, denominator):
    from flint import acb, fmpq
    return acb(fmpq(2 * (numerator % denominator), denominator)).exp_pi_i()


def _rotation_coefficients(theta):
    from flint import acb, arb, fmpq
    sine, cosine = arb(fmpq(theta.numerator, 2 * theta.denominator)).sin_cos_pi()
    return acb(cosine), acb(0, -sine)


def _apply_background(vector, axis, cosine, minus_i_sine):
    """The exact-input two-level rotation, identity on all other fine labels."""
    result = list(vector)
    a, b = result[:2]
    c, d = cosine, minus_i_sine
    result[:2] = ([c*a+d*b, d*a+c*b] if axis == "x"
                  else [(c+d)*a, (c-d)*b])
    return result


@dataclass(frozen=True)
class DyadicPrefix:
    coordinates: tuple[tuple[int, int], ...]
    accuracy_bits: int
    working_precision: int
    refinements: int
    max_ball_radius: Fraction

    @property
    def weights(self):
        return tuple(re * re + im * im for re, im in self.coordinates)

    @property
    def coordinate_error_bound(self):
        return Fraction(1, 1 << self.accuracy_bits)


class VerifiedReflectionCircuit:
    """Exact rational-pi Rx/Rz background and K_q reflection rotations.

    backgrounds: insertion -> ('x' or 'z', theta/pi as int or Fraction).
    block_size=2 by default; opt-in 3 embeds that rotation on labels 0/1
    and fixes label 2. This is not an arbitrary exact-unitary input format.
    reflections: insertion -> (integer q, theta/pi as int or Fraction).
    W then K at each insertion, including 0 and width. Every input rational
    numerator/denominator is capped at 2048 bits before exact mod-4 reduction.
    Width/history caps bound structure, not the refinement working precision.
    enclosure_mode='rectangular' retains ordinary Acb error propagation.
    Opt-in 'norm' transports prior component error with exact step norm<=1,
    evaluating each new local residual on an exact dyadic midpoint. It can
    reduce required refinement while adding bookkeeping; faster is not promised.
    """
    def __init__(self, period, width, backgrounds, reflections, *, enclosure_mode="rectangular",
                 block_size=2):
        if (type(block_size) is not int or block_size not in (2, 3)
                or type(period) is not int or not block_size <= period <= (1 << 63) - 1
                or period % block_size or type(width) is not int or not 0 <= width <= 63
                or not isinstance(backgrounds, dict) or not isinstance(reflections, dict)
                or len(reflections) > 8):
            raise ValueError("require block size 2/3 dividing int64 period, width<=63, and <=8 reflections")
        if enclosure_mode not in ("rectangular", "norm"):
            raise ValueError("enclosure_mode must be rectangular or norm")

        def angle(value):
            if type(value) not in (int, Fraction):
                raise ValueError("angles must be exact rational multiples of pi")
            value = Fraction(value)
            if max(abs(value.numerator).bit_length(), value.denominator.bit_length()) > 2048:
                raise ValueError("rational angle input exceeds 2048-bit cap")
            return value % 4

        def validate(entries, background):
            result = {}
            for s, pair in sorted(entries.items(), key=lambda kv: str(kv[0])):
                if (type(s) is not int or not 0 <= s <= width
                        or not isinstance(pair, (tuple, list)) or len(pair) != 2):
                    raise ValueError("invalid exact-gate insertion")
                kind, theta = pair
                if background:
                    if kind not in ("x", "z"):
                        raise ValueError("background must be Rx or Rz")
                else:
                    if type(kind) is not int or abs(kind).bit_length() > 2048:
                        raise ValueError("reflection q must be a bounded exact integer")
                    kind %= period // block_size
                result[s] = (kind, angle(theta))
            return MappingProxyType(dict(sorted(result.items())))

        self.period, self.sectors, self.b, self.width = period, period // block_size, block_size, width
        self.enclosure_mode = enclosure_mode
        self.backgrounds = validate(backgrounds, True)
        self.reflections = validate(reflections, False)

    def _label(self, sector, exponent, stop, boundary, measured, output):
        if (any(type(v) is not int for v in (sector, exponent, stop, measured, output))
                or not 0 <= sector < self.sectors
                or not 0 <= exponent < 1 << self.width
                or not 0 <= stop <= self.width
                or boundary not in ("arithmetic", "background", "reflection")
                or not 0 <= measured <= self.width or not 0 <= output < 1 << measured
                or measured and (stop != self.width or boundary != "reflection")):
            raise ValueError("invalid gate-prefix or partial-QFT amplitude label")

    def _enclosure(self, sector, exponent, stop, boundary, measured, output):
        from flint import acb
        result = [acb(0) for _ in range(self.b)]
        for coefficient, vector in self._history_enclosures(
                sector, exponent, stop, boundary, measured, output):
            result = [old + coefficient * v for old, v in zip(result, vector)]
        return tuple(result)

    def _history_enclosures(self, sector, exponent, stop, boundary, measured, output,
                            *, route_only=False):
        """Stream the shared contraction's components at the current precision.

        route_only is for a circuit whose listed reflections are deterministic
        J routes: select all of them, with coefficient one. This does not drop
        any coherent history in the ordinary exact-input circuit.
        """
        # Rebuild every rounded scalar at the current ctx precision. Nothing
        # computed at a previous precision is reused in the refined expression.
        from flint import acb, arb, fmpq

        phase, rotation = _rational_phase, _rotation_coefficients

        matrices = {s: (axis, *rotation(theta))
                    for s, (axis, theta) in self.backgrounds.items()}
        terms = [(s, q, *rotation(theta)) for s, (q, theta) in self.reflections.items()
                 if s < stop or s == stop and boundary == "reflection"]
        inverse_sqrt_two = 1 / arb(2).sqrt()
        tail_scale = inverse_sqrt_two ** (self.width - stop)

        def shift(power, alpha, vector):
            result = [acb(0) for _ in range(self.b)]
            for p, value in enumerate(vector):
                wrap, target = divmod(p + power, self.b)
                result[target] = phase(alpha * wrap, self.sectors) * value
            return result

        def defect(s, vector):
            if s not in matrices:
                return vector
            axis, c, d = matrices[s]  # d = -i sin(theta/2)
            return _apply_background(vector, axis, c, d)

        masks = ((1 << len(terms)) - 1,) if route_only else range(1 << len(terms))
        for mask in masks:
            coefficient, selected = acb(1), {}
            for i, (s, q, c, d) in enumerate(terms):
                chosen = mask & (1 << i)
                if not route_only:
                    coefficient *= d if chosen else c
                if chosen:
                    selected[s] = q
            # Each component factor is an exact unitary or contraction in l2:
            # W, T/sqrt(2), (I+zT)/2, and the unused-control scalar tail.
            # Apply each to an EXACT dyadic center and charge the fresh local
            # box by a Euclidean radius bound. Prior error travels
            # with norm <=1, not through repeatedly re-boxed complex phases.
            norm_radius = arb(0)
            def norm_step(operation, vector, tag):
                nonlocal norm_radius
                local = tuple(acb(z) for z in operation(vector))
                # 2*b real coordinates: l2 radius <=ceil(sqrt(2*b))*max.
                # Integer ceiling keeps the old b=2 factor exactly 2 and
                # raises it to 3 for b=3. rad() is exact binary data.
                radius_factor = isqrt(2*self.b)
                radius_factor += radius_factor**2 < 2*self.b
                local_radius = radius_factor * max((radius for z in local
                                        for radius in (z.real.rad(), z.imag.rad())), default=arb(0))
                norm_radius = (norm_radius + local_radius).upper()
                return [z.mid() for z in local]
            vector = _contract_route_component(
                self.width, self.b, self.sectors, sector, exponent, stop, boundary,
                measured, output, selected, shift=shift, defect=defect, phase=phase,
                inverse_sqrt_two=inverse_sqrt_two, tail_scale=tail_scale,
                apply_step=norm_step if self.enclosure_mode == "norm" else None)
            if self.enclosure_mode == "norm":
                # Convert the one component-vector l2 enclosure to coordinate
                # boxes ONCE, before the short coherent history summation.
                vector = [acb(arb(z.real, norm_radius), arb(z.imag, norm_radius)) for z in vector]
            # Evaluate all finite histories. No threshold or uncertain-zero test.
            yield coefficient, tuple(vector)

    def prefix_enclosure(self, sector, exponent, stop, *, boundary="reflection",
                         measured=0, output=0, working_precision=128):
        """Diagnostic fixed-precision ball, not a sampler with a failure cap."""
        self._label(sector, exponent, stop, boundary, measured, output)
        if type(working_precision) is not int or working_precision < 16:
            raise ValueError("require working precision >=16")
        from flint import ctx
        with ctx.workprec(working_precision):
            return self._enclosure(sector, exponent, stop, boundary, measured, output)

    def prefix_dyadic(self, sector, exponent, stop, *, accuracy_bits,
                      boundary="reflection", measured=0, output=0, max_precision=None):
        """Uniform absolute real/imag error <=2^-p for every valid label.

        Deterministic per-label refinement; no trajectory-dependent cache.
        Ball radius <=eta/2 plus exact midpoint rounding <=eta/2. Analytic
        finite expression with no small-amplitude divisions converges at every
        label, including zeros. No uniform runtime bound is inferred from a
        measured precision. max_precision is DIAGNOSTIC ONLY and may raise.
        """
        self._label(sector, exponent, stop, boundary, measured, output)
        _bits(accuracy_bits)
        if max_precision is not None and (type(max_precision) is not int or max_precision < 16):
            raise ValueError("invalid diagnostic precision cap")
        from flint import ctx
        precision, attempts = max(64, accuracy_bits + 16), 0
        tolerance = Fraction(1, 1 << (accuracy_bits + 1))
        while max_precision is None or precision <= max_precision:
            attempts += 1
            with ctx.workprec(precision):
                vector = self._enclosure(sector, exponent, stop, boundary, measured, output)
                scalars = tuple(z for v in vector for z in (v.real, v.imag))
                if all(z.is_finite() for z in scalars):
                    radii = tuple(binary_fraction(z.rad().upper()) for z in scalars)
                    if max(radii) <= tolerance:
                        integers = tuple(round_dyadic(binary_fraction(z.mid()), accuracy_bits)
                                         for z in scalars)
                        return DyadicPrefix(tuple(zip(integers[::2], integers[1::2])),
                                            accuracy_bits, precision, attempts - 1, max(radii))
            precision *= 2
        raise ArithmeticError("diagnostic precision cap exhausted; no conditioned-law certificate")

    def accuracy_plan(self, target_tv):
        if type(target_tv) not in (int, Fraction) or not 0 < target_tv < 1:
            raise ValueError("target TV must be an exact rational strictly between 0 and 1")
        updates = len(self.backgrounds) + len(self.reflections) + self.width
        if updates:
            return plan_prefix_accuracy(self.width, self.b, updates, target_tv)
        return dict(**prefix_error_budget(0, self.b, [], random_bits=1),
                    absolute_accuracy_bits=1, target_tv=Fraction(target_tv),
                    requested_real_imag_abs_error=Fraction(1, 2), machine_mantissa_bits=None)

    def sample(self, rng, *, target_tv=Fraction(1, 1_000_000)):
        """Covered gate-by-gate law under verified arithmetic/unbiased-bit assumptions.

        Returns a bound for the complete algorithm, not a path-specific
        certificate. No order/index discovery, certified rejection, or PRNG
        discrepancy bound is supplied by this specialization.
        """
        plan = self.accuracy_plan(target_tv)
        accuracy_bits, random_bits = plan["absolute_accuracy_bits"], plan["random_bits"]
        sector = uniform_integer(self.sectors, rng)
        exponent = _word(rng, self.width)
        p, evaluations, updates, refinements, max_working, zero_blocks = 0, 0, 0, 0, 0, 0

        def query(alpha, stop, **kw):
            nonlocal evaluations, refinements, max_working
            result = self.prefix_dyadic(alpha, exponent, stop, accuracy_bits=accuracy_bits, **kw)
            evaluations += 1
            refinements += result.refinements
            max_working = max(max_working, result.working_precision)
            return result.weights

        def draw(weights):
            nonlocal updates, zero_blocks
            updates += 1
            zero_blocks += not any(weights)
            return integer_cdf_index(weights, _word(rng, random_bits), random_bits)

        for stop in range(self.width + 1):
            if stop and exponent & (1 << (stop - 1)):
                p = (p + (1 << (stop - 1))) % self.b
            if stop in self.backgrounds:
                p = draw(query(sector, stop, boundary="background"))
            if stop in self.reflections:
                target = (-sector - self.reflections[stop][0]) % self.sectors
                if target != sector:
                    labels = (sector, target)
                    sector = labels[draw([query(alpha, stop)[p] for alpha in labels])]
        output = 0
        for measured in range(1, self.width + 1):
            children = (output, output | (1 << (measured - 1)))
            output = children[draw([query(sector, self.width, measured=measured, output=y)[p]
                                    for y in children])]
        return dict(output=output, final_coarse_sector=sector, final_within_sector=p,
                    total_tv_upper_bound=plan["total_tv_upper_bound"], target_tv=Fraction(target_tv),
                    absolute_accuracy_bits=accuracy_bits, categorical_random_bits=random_bits,
                    max_working_precision=max_working, refinement_retries=refinements,
                    prefix_vector_evaluations=evaluations, resampled_blocks=updates,
                    zero_approximate_block_fallbacks=zero_blocks,
                    block_updates_upper_bound=plan["updates"],
                    history_count_upper_bound=1 << len(self.reflections),
                    method="verified exact-input sparse-block prefix sampling",
                    arithmetic_backend="python-flint/Arb/Acb",
                    enclosure_mode=self.enclosure_mode,
                    requires_independent_unbiased_bits=True,
                    certifies_existing_float_sampler=False, precision_cap=None,
                    orbit_or_output_tables=False)
