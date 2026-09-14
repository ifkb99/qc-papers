"""Finite-TV reverse-vector rejection for one exact-input work component.

This is C68's pure reverse trajectory with its norm kept UNNORMALIZED.
Children use the existing inverse-QFT operator; only their nonnegative
dyadic masses are normalized, by the exact integer CDF. There are no forward
states, branch tables, orbit/output arrays, vector renormalizations or cutoffs.
C64's mass-tree bound and C63's one-sided accepted-mass inequality certify
the composition. This is not a coherent-history outer acceptance algorithm.

Optional python-flint==0.9.0, convergent exact-input enclosures, independent
unbiased bits and a nonconcurrent precision context remain trust assumptions.
Precision and rejection are uncapped. Resource failures do not certify a
success-conditioned output law. Matrix entry counts are not native RSS.
"""
from __future__ import annotations

from fractions import Fraction

from lab.sampling_error import plan_prefix_mass_accuracy
from lab.semiclassical import _qft_branch_operators
from lab.verified_finite_work import VerifiedFiniteWork, _adjoint
from lab.verified_prefix import (binary_fraction, round_dyadic, integer_cdf_index,
    integer_cdf_counts, uniform_integer, _rational_phase, _word, _bits)
from lab.verified_rejection import (AcceptanceThreshold, conservative_threshold,
                                    _positive, _accuracy_bits)


def _norm_squared(vector):
    from flint import arb
    # General Arb exponentiation may be indeterminate across zero; these
    # finite polynomial expressions must use multiplication, as in C63.
    value = arb(0)
    for i in range(vector.nrows()):
        z = vector[i,0]
        value += z.real*z.real + z.imag*z.imag
    return value


class VerifiedReverseWork:
    """Sample one normalized deterministic history, final fine work traced.

    Initial sector is drawn exactly once per returned sample; EACH rejection
    attempt draws a fresh uniform final-boundary label j in [0,b). Keeping j
    fixed through rejection would generally return a different output law.
    """
    def __init__(self, circuit, mask):
        # Reuse the exact input snapshot, route semantics and branch builder.
        # Never call its full-history _build or create its density cursor.
        self.builder = VerifiedFiniteWork(circuit, mask)
        self.circuit, self.mask = self.builder.circuit, mask

    def plan(self, target_tv):
        target = _positive(target_tv)
        if target >= 1:
            raise ValueError("target TV must be below one")
        b, Q = self.circuit.b, 1 << self.circuit.width
        proposal_target = target/(8*b)
        proposal = plan_prefix_mass_accuracy(self.circuit.width, proposal_target)
        tolerance = target/(4*b*Q)
        _accuracy_bits(tolerance)  # validate before drawing any randomness
        acceptance_bits = _accuracy_bits(target/(4*b))
        proposal_error = proposal["total_tv_upper_bound"]
        beta = b*(2*proposal_error + Q*tolerance + Fraction(1,1 << acceptance_bits))
        return dict(target_tv=target, total_tv_upper_bound=beta,
            proposal_target_tv=proposal_target, proposal_tv_upper_bound=proposal_error,
            proposal_plan=proposal, acceptance_interval_width_tolerance=tolerance,
            acceptance_random_bits=acceptance_bits, envelope_upper_bound=b,
            success_probability_lower_bound=(1-beta)/b,
            expected_attempts_upper_bound=Fraction(b)/(1-beta),
            machine_mantissa_bits=None)

    def cursor(self, initial_sector, boundary_j, *, accuracy_bits, initial_precision=None):
        return _ReverseCursor(self, initial_sector, boundary_j, accuracy_bits, initial_precision)

    def attempt(self, initial_sector, boundary_j, *, rng=None, output=None,
                target_tv=Fraction(1,1_000_000), initial_precision=None):
        """Propose a complete path and compute its finite acceptance threshold.

        Does NOT perform the acceptance coin or return an accepted-law sample.
        Forced output returns its exact finite proposal probability conditional
        on initial_sector AND boundary_j, before their external mixture.
        """
        if (rng is None) == (output is None):
            raise ValueError("supply rng or forced output, exclusively")
        t, b = self.circuit.width, self.circuit.b
        if output is not None and (type(output) is not int or not 0 <= output < 1 << t):
            raise ValueError("forced output outside width")
        plan = self.plan(target_tv)
        proposal = plan["proposal_plan"]
        cursor = self.cursor(initial_sector,boundary_j,
            accuracy_bits=proposal["weight_accuracy_bits"],initial_precision=initial_precision)
        probability = Fraction(1)
        for step in range(t):
            weights = cursor.weights()
            L = proposal["categorical_random_bits"]
            if output is None:
                bit = integer_cdf_index(weights,_word(rng,L),L)
            else:
                bit = (output >> step)&1
                probability *= Fraction(integer_cdf_counts(weights,L)[bit],1 << L)
            cursor.advance(bit)
        decision = cursor.acceptance(
            interval_tolerance=plan["acceptance_interval_width_tolerance"],
            random_bits=plan["acceptance_random_bits"])
        result = dict(**plan,output=cursor.output,initial_coarse_sector=initial_sector,
            final_coarse_sector=self.builder._sector_at(initial_sector,t),boundary_j=boundary_j,
            acceptance_probability=decision.probability,
            acceptance_threshold=decision.threshold,
            acceptance_numerator_bounds=decision.numerator_bounds,
            acceptance_denominator_bounds=decision.denominator_bounds,
            acceptance_width_sum=decision.width_sum,
            acceptance_enclosure_evaluations=cursor.acceptance_evaluations,
            initial_vector_builds=cursor.builds,
            acceptance_inner_products=cursor.acceptance_evaluations,
            reverse_child_evaluations=cursor.evaluations,
            reverse_vector_matvecs=cursor.matvecs,
            replayed_vector_steps=cursor.replayed,
            branch_matrix_pair_constructions=cursor.branch_builds,
            refinement_retries=cursor.builds-1,max_working_precision=cursor.precision,
            zero_approximate_block_fallbacks=cursor.zero_blocks,
            resampled_blocks=t,forward_steps=0,stored_forward_matrix_count=0,
            stored_branch_matrix_count=0,retained_work_vector_count_upper_bound=4,
            stored_work_vector_scalar_count_upper_bound=4*b,
            working_matrix_scalar_count_upper_bound=24*b*b+8*b,
            work_block_size=b,orbit_or_output_tables=False,
            method="verified unnormalized reverse-vector proposal with boundary rejection",
            vector_normalization=False,requires_independent_unbiased_bits=True,
            certifies_existing_float_sampler=False,precision_cap=None,attempt_cap=None,
            initial_precision=initial_precision)
        if output is not None:
            result["proposal_path_probability"] = probability
        return result

    def sample(self, rng, *, target_tv=Fraction(1,1_000_000)):
        """Accepted (final coarse sector, output) with the declared TV bound."""
        plan = self.plan(target_tv)
        initial = uniform_integer(self.circuit.sectors,rng)
        totals = {key:0 for key in (
            "acceptance_enclosure_evaluations","initial_vector_builds",
            "acceptance_inner_products","reverse_child_evaluations",
            "reverse_vector_matvecs","replayed_vector_steps",
            "branch_matrix_pair_constructions","refinement_retries",
            "zero_approximate_block_fallbacks","resampled_blocks")}
        attempts, precision = 0, 0
        while True:
            attempts += 1
            j = uniform_integer(self.circuit.b,rng)
            proposed = self.attempt(initial,j,rng=rng,target_tv=target_tv)
            for key in totals:
                totals[key] += proposed[key]
            precision = max(precision,proposed["max_working_precision"])
            if _word(rng,plan["acceptance_random_bits"]) < proposed["acceptance_threshold"]:
                # j belongs to the proposal, not the target's final work law.
                return dict(**plan,output=proposed["output"],initial_coarse_sector=initial,
                    final_coarse_sector=proposed["final_coarse_sector"],attempts=attempts,
                    **totals,max_working_precision=precision,
                    forward_steps=0,stored_forward_matrix_count=0,stored_branch_matrix_count=0,
                    retained_work_vector_count_upper_bound=4,
                    stored_work_vector_scalar_count_upper_bound=4*self.circuit.b,
                    working_matrix_scalar_count_upper_bound=24*self.circuit.b**2+8*self.circuit.b,
                    work_block_size=self.circuit.b,orbit_or_output_tables=False,
                    method="verified reverse-vector boundary rejection",
                    requires_independent_unbiased_bits=True,vector_normalization=False,
                    certifies_existing_float_sampler=False,precision_cap=None,attempt_cap=None)


class _ReverseCursor:
    def __init__(self, owner, initial, boundary, bits, initial_precision):
        _bits(bits)
        if type(initial) is not int or not 0 <= initial < owner.circuit.sectors:
            raise ValueError("initial sector outside circuit")
        if type(boundary) is not int or not 0 <= boundary < owner.circuit.b:
            raise ValueError("boundary label outside work block")
        if initial_precision is not None and (type(initial_precision) is not int or initial_precision < 16):
            raise ValueError("initial precision must be None or integer>=16, never a cap")
        self.owner,self.initial,self.boundary,self.bits = owner,initial,boundary,bits
        self.precision = max(64,bits+16) if initial_precision is None else initial_precision
        self.output,self.measured,self.builds,self.replayed = 0,0,0,0
        self.evaluations,self.acceptance_evaluations,self.zero_blocks = 0,0,0
        self.matvecs,self.branch_builds = 0,0
        self.vector,self.initial_vector,self._pending = None,None,None
        self._rebuild()

    def _operators(self, depth):
        b = self.owner.builder
        i = self.owner.circuit.width-1-depth
        branches = b._branch(i,b._sector_at(self.initial,i),b._identity())
        self.branch_builds += 1
        phase = _rational_phase(-(self.output % (1 << depth)),1 << (depth+1))
        return _qft_branch_operators(*branches,phase)

    def _rebuild(self):
        from flint import acb_mat,ctx
        self.vector,self.initial_vector,self._pending = None,None,None
        with ctx.workprec(self.precision):
            b = self.owner.circuit.b
            self.vector = acb_mat([[int(i == self.boundary)] for i in range(b)])
            builder = self.owner.builder
            self.initial_vector = builder._background(0,builder._identity())*acb_mat(
                [[int(i == 0)] for i in range(b)])
            for depth in range(self.measured):
                K = self._operators(depth)[(self.output >> depth)&1]
                self.vector = _adjoint(K)*self.vector
                self.replayed += 1
                self.matvecs += 1
        self.builds += 1

    def weights(self):
        if self.measured >= self.owner.circuit.width:
            raise ValueError("all outputs already measured")
        if self._pending is not None:
            return self._pending[0]
        from flint import ctx
        tolerance = Fraction(1,1 << (self.bits+1))
        while True:
            with ctx.workprec(self.precision):
                children = tuple(_adjoint(K)*self.vector for K in self._operators(self.measured))
                self.matvecs += 2
                self.evaluations += 1
                masses = tuple(_norm_squared(v) for v in children)
                if all(z.is_finite() for z in masses):
                    if any(z.upper() < 0 for z in masses):
                        raise ArithmeticError("invalid nonnegative reverse mass")
                    if max(binary_fraction(z.rad().upper()) for z in masses) <= tolerance:
                        weights = tuple(max(0,round_dyadic(binary_fraction(z.mid()),self.bits)) for z in masses)
                        self.zero_blocks += not any(weights)
                        self._pending = weights,children
                        return weights
            self.precision *= 2
            self._rebuild()

    def advance(self, bit):
        if type(bit) is not int or bit not in (0,1) or self._pending is None:
            raise ValueError("query weights then advance by one exact bit")
        self.vector = self._pending[1][bit]
        self.output |= bit << self.measured
        self.measured += 1
        self._pending = None

    def acceptance(self, *, interval_tolerance, random_bits):
        if self.measured != self.owner.circuit.width:
            raise ValueError("terminal acceptance requires every output bit")
        tolerance = _positive(interval_tolerance)
        _bits(random_bits)
        from flint import ctx
        retries = 0
        while True:
            with ctx.workprec(self.precision):
                N = _norm_squared(_adjoint(self.initial_vector)*self.vector)
                D = _norm_squared(self.vector)
                self.acceptance_evaluations += 1
                if N.is_finite() and D.is_finite():
                    nl,nu = binary_fraction(N.lower()),binary_fraction(N.upper())
                    dl,du = binary_fraction(D.lower()),binary_fraction(D.upper())
                    if nu < 0 or du < 0 or nl > du:
                        raise ArithmeticError("invalid reverse acceptance enclosure for 0<=N<=D")
                    width = nu-nl+du-dl
                    if width <= tolerance:
                        return AcceptanceThreshold(conservative_threshold(nl,du,random_bits),
                            random_bits,(nl,nu),(dl,du),width,self.precision,retries)
            self.precision *= 2
            retries += 1
            self._rebuild()
