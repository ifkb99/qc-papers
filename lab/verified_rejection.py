"""Finite-TV coherent-component rejection using verified proposals AND acceptance.

Same supplied exact-input contract as VerifiedReflectionCircuit, b=2/3. This
is not certification of C59's float finite-work instrument. Default component
proposals reuse C61's prefix sampler. Opt-in component_method='finite_work'
uses C64's verified density/effect implementation of the same finite-work
identity; component_method='scalar' uses C65's stronger b=2 specialization.
Every attempt draws a new history and initial sector.
Only (final sector, exponent output), NOT the proposal's internal work label,
has the returned target-law certificate. See C63 for accepted-measure proof.

Arithmetic refinement and rejection are uncapped. Backend/resource failures
do not certify the law conditioned on successful runs. No real-probability
equality test, minimum proposal probability, or relative amplitude error is
required. Ideal independent unbiased bits remain an explicit assumption.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction, _word, _bits


def _positive(value):
    if type(value) not in (int, Fraction) or value <= 0:
        raise ValueError("require a positive exact int/Fraction")
    return Fraction(value)


def _accuracy_bits(error):
    bits = 1
    while Fraction(1, 1 << bits) > error:
        bits += 1
        if bits > 4096:
            raise ValueError("accuracy exceeds the bounded planner")
    return bits


def conservative_threshold(numerator_lower, denominator_upper, random_bits):
    """Exact downward Bernoulli grid. Requires valid bounds on 0<=P<=D.

    Clamps implement nonnegative lower P and the known acceptance bound one.
    When D_upper=0, both true quantities are zero and acceptance is zero.
    This routine alone cannot establish that supplied bounds are valid.
    """
    _bits(random_bits)
    if any(type(x) not in (int, Fraction) for x in (numerator_lower, denominator_upper)):
        raise ValueError("acceptance endpoints must be exact")
    if denominator_upper < 0:
        raise ValueError("upper denominator must be nonnegative")
    if denominator_upper == 0:
        if numerator_lower > 0:
            raise ValueError("inconsistent zero denominator")
        return 0
    ratio = min(Fraction(1), max(Fraction(0), Fraction(numerator_lower)) /
                Fraction(denominator_upper))
    scaled = ratio * (1 << random_bits)
    return scaled.numerator // scaled.denominator


class _RoutedComponent(VerifiedReflectionCircuit):
    """One normalized J-route circuit, global K(pi) phases omitted exactly.

    The inherited block sampler uses the same partitions as K(pi); replacing
    K(pi) by J changes each prefix only by a global phase. Thus all queried
    Born weights and the C61 kernel proof are unchanged. Only one component
    contraction is evaluated, even with many selected routes.
    """
    def _enclosure(self, sector, exponent, stop, boundary, measured, output):
        return next(self._history_enclosures(sector, exponent, stop, boundary,
                                             measured, output, route_only=True))[1]


@dataclass(frozen=True)
class AcceptanceThreshold:
    threshold: int
    random_bits: int
    numerator_bounds: tuple[Fraction, Fraction]
    denominator_bounds: tuple[Fraction, Fraction]
    width_sum: Fraction
    working_precision: int
    refinements: int

    @property
    def probability(self):
        return Fraction(self.threshold, 1 << self.random_bits)


class VerifiedRejectionSampler:
    def __init__(self, circuit, *, component_method="prefix"):
        if type(circuit) is not VerifiedReflectionCircuit:
            raise ValueError("require the verified exact-input circuit, not a float/subclass oracle")
        if component_method not in ("prefix","finite_work","scalar"):
            raise ValueError("component_method must be prefix, finite_work or scalar")
        if component_method == "scalar" and circuit.b != 2:
            raise ValueError("scalar component requires block_size=2")
        self.component_method = component_method
        # Snapshot the declared gates; caller mutation cannot change a running
        # sampler's input contract. No orbit, output, or history-vector tables.
        self.circuit = VerifiedReflectionCircuit(circuit.period, circuit.width,
            dict(circuit.backgrounds), dict(circuit.reflections),
            enclosure_mode=circuit.enclosure_mode, block_size=circuit.b)
        self.k = len(circuit.reflections)
        self.envelope_upper = 1 << self.k

    def component(self, mask):
        if type(mask) is not int or not 0 <= mask < 1 << self.k:
            raise ValueError("history mask outside finite expansion")
        c = self.circuit
        routes = {s: (q, 1) for i, (s, (q, _)) in enumerate(c.reflections.items())
                  if mask & (1 << i)}
        return _RoutedComponent(c.period, c.width, dict(c.backgrounds), routes,
                                enclosure_mode=c.enclosure_mode, block_size=c.b)

    def plan(self, target_tv):
        target = _positive(target_tv)
        if target >= 1:
            raise ValueError("target TV must be below one")
        C, Q = self.envelope_upper, 1 << self.circuit.width
        component_target = target / (8 * C)
        history_bits = _accuracy_bits(target / (16 * C * self.k)) if self.k else 1
        history_error = Fraction(2 * self.k, 1 << history_bits)
        acceptance_bits = _accuracy_bits(target / (4 * C))
        interval_tolerance = target / (4 * Q)
        # Validate before doing any science/RNG work for unsupported requests.
        _accuracy_bits(interval_tolerance)
        if self.component_method == "prefix":
            self.circuit.accuracy_plan(component_target)
        elif self.component_method == "finite_work":
            from lab.sampling_error import plan_prefix_mass_accuracy
            plan_prefix_mass_accuracy(self.circuit.width,component_target)
        else:
            from lab.sampling_error import plan_conditional_weight_accuracy
            plan_conditional_weight_accuracy(self.circuit.width,component_target)
        proposal_error = component_target + history_error
        total = 2 * C * proposal_error + Q * interval_tolerance + Fraction(C, 1 << acceptance_bits)
        return dict(target_tv=target, total_tv_upper_bound=total,
                    component_target_tv=component_target,
                    history_random_bits=history_bits, history_tv_upper_bound=history_error,
                    proposal_tv_upper_bound=proposal_error,
                    acceptance_random_bits=acceptance_bits,
                    acceptance_interval_width_tolerance=interval_tolerance,
                    envelope_upper_bound=C,
                    success_probability_lower_bound=(1-total)/C,
                    expected_attempts_upper_bound=Fraction(C)/(1-total))

    def history_thresholds(self, random_bits):
        """Fixed dyadic independent choices, total TV <=2*k*2^-L.

        Each exact probability is |sin|/(|cos|+|sin|), denominator >=1.
        Refine its interval to width <=2^-L, then floor its lower endpoint.
        This is a bounded-error approximation, not an exact binary expansion.
        """
        _bits(random_bits)
        from flint import arb, ctx, fmpq
        precision, retries = max(64, random_bits + 16), 0
        tolerance = Fraction(1, 1 << random_bits)
        while True:
            with ctx.workprec(precision):
                bounds = []
                for _, theta in self.circuit.reflections.values():
                    sine, cosine = arb(fmpq(theta.numerator, 2*theta.denominator)).sin_cos_pi()
                    probability = abs(sine) / (abs(cosine) + abs(sine))
                    if not probability.is_finite():
                        break
                    bounds.append((binary_fraction(probability.lower()),
                                   binary_fraction(probability.upper())))
                if len(bounds) == self.k and all(hi-lo <= tolerance for lo, hi in bounds):
                    thresholds = tuple(conservative_threshold(lo, 1, random_bits) for lo, _ in bounds)
                    return dict(thresholds=thresholds, random_bits=random_bits, bounds=tuple(bounds),
                                working_precision=precision if self.k else 0, refinements=retries)
            precision *= 2
            retries += 1

    def acceptance_enclosure(self, sector, output, *, working_precision=128):
        """Enclose P=M*p and D=M*B^2*q, streamed, without dividing by q."""
        c = self.circuit
        c._label(sector, 0, c.width, "reflection", c.width, output)
        if type(working_precision) is not int or working_precision < 16:
            raise ValueError("require working precision >=16")
        from flint import acb, arb, ctx
        with ctx.workprec(working_precision):
            coherent, mixture, B = [acb(0) for _ in range(c.b)], arb(0), arb(0)
            for coefficient, vector in c._history_enclosures(
                    sector, 0, c.width, "reflection", c.width, output):
                magnitude = abs(coefficient)
                B += magnitude
                mixture += magnitude * sum((z.real*z.real+z.imag*z.imag for z in vector), arb(0))
                coherent = [old+coefficient*z for old, z in zip(coherent, vector)]
            P = sum((z.real*z.real+z.imag*z.imag for z in coherent), arb(0))
            return P, B * mixture

    def acceptance(self, sector, output, *, interval_tolerance, random_bits):
        tolerance = _positive(interval_tolerance)
        _bits(random_bits)
        precision, retries = max(64, _accuracy_bits(tolerance)+16), 0
        while True:
            P, D = self.acceptance_enclosure(sector, output, working_precision=precision)
            # Endpoints must be extracted at the evaluation precision too;
            # default-context rounding must not set the refinement floor.
            from flint import ctx
            with ctx.workprec(precision):
                if P.is_finite() and D.is_finite():
                    pl, pu = binary_fraction(P.lower()), binary_fraction(P.upper())
                    dl, du = binary_fraction(D.lower()), binary_fraction(D.upper())
                    width = pu-pl+du-dl
                    if du < 0 or pu < 0 or pl > du:
                        raise ArithmeticError("invalid enclosure for nonnegative P<=D")
                    if width <= tolerance:
                        threshold = conservative_threshold(pl, du, random_bits)
                        return AcceptanceThreshold(threshold, random_bits, (pl, pu), (dl, du),
                                                   width, precision, retries)
            precision *= 2
            retries += 1

    def sample(self, rng, *, target_tv=Fraction(1, 1_000_000)):
        plan = self.plan(target_tv)
        histories = self.history_thresholds(plan["history_random_bits"])
        attempts, queries, acceptance_evaluations = 0, 0, 0
        refinements, max_working = histories["refinements"], histories["working_precision"]
        updates, zero_blocks, forward_steps, replay_steps, backward_evaluations, scalar_evaluations = 0, 0, 0, 0, 0, 0
        while True:
            attempts += 1
            mask = sum((1 << i) for i, threshold in enumerate(histories["thresholds"])
                       if _word(rng, histories["random_bits"]) < threshold)
            if self.component_method == "prefix":
                proposed = self.component(mask).sample(rng,target_tv=plan["component_target_tv"])
            else:
                from lab.verified_finite_work import VerifiedFiniteWork,VerifiedScalarWork
                factory = VerifiedFiniteWork if self.component_method == "finite_work" else VerifiedScalarWork
                proposed = factory(self.circuit,mask).sample(rng,target_tv=plan["component_target_tv"])
            queries += proposed["prefix_vector_evaluations"]
            refinements += proposed["refinement_retries"]
            updates += proposed["resampled_blocks"]
            zero_blocks += proposed["zero_approximate_block_fallbacks"]
            forward_steps += proposed.get("forward_steps",0)
            replay_steps += proposed.get("replayed_effect_steps",0)
            backward_evaluations += proposed.get("backward_child_evaluations",0)
            scalar_evaluations += proposed.get("scalar_conditional_evaluations",0)
            decision = self.acceptance(proposed["final_coarse_sector"], proposed["output"],
                interval_tolerance=plan["acceptance_interval_width_tolerance"],
                random_bits=plan["acceptance_random_bits"])
            acceptance_evaluations += decision.refinements + 1
            refinements += decision.refinements
            max_working = max(max_working, proposed["max_working_precision"], decision.working_precision)
            if _word(rng, decision.random_bits) < decision.threshold:
                return dict(**plan, output=proposed["output"],
                    final_coarse_sector=proposed["final_coarse_sector"],
                    attempts=attempts, prefix_vector_evaluations=queries,
                    acceptance_enclosure_evaluations=acceptance_evaluations,
                    refinement_retries=refinements, max_working_precision=max_working,
                    resampled_blocks=updates, zero_approximate_block_fallbacks=zero_blocks,
                    history_count_upper_bound=1 << self.k, component_history_count=1,
                    method="verified finite-TV rejection",
                    component_method=self.component_method,
                    component_enclosure_mode=(self.circuit.enclosure_mode if self.component_method == "prefix"
                                              else "rectangular scalar/vector" if self.component_method == "scalar"
                                              else "rectangular density/effect"),
                    component_scalar_evaluations=scalar_evaluations,
                    component_forward_steps=forward_steps, component_replayed_effect_steps=replay_steps,
                    component_backward_child_evaluations=backward_evaluations,
                    enclosure_mode=self.circuit.enclosure_mode, orbit_or_output_tables=False,
                    requires_independent_unbiased_bits=True, certifies_existing_float_sampler=False,
                    precision_cap=None, attempt_cap=None)
