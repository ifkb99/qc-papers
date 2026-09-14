"""Locally certified, bounded-grid reverse quantum trajectories.

Exact-input branches/QFT operators are shared with C69. Their dyadic
enclosures are followed by exact Gaussian-integer vector updates, exact CDFs
and projective rounding on a fixed grid. Local joint classical/quantum trace
errors telescope with depth, not output-tree size. Terminal boundary
rejection is charged too. No forward table or physical rare-norm division.

Optional convergent Arb/Acb and independent unbiased bits remain assumptions.
Operator refinement and rejection are uncapped; resource-conditioned success
is not certified. Bounded vector integers are not a native RSS or complete
backend bit-complexity claim. This is one deterministic route component.
"""
from __future__ import annotations

from fractions import Fraction

from lab.verified_reverse_work import VerifiedReverseWork
from lab.verified_finite_work import _adjoint
from lab.semiclassical import _qft_branch_operators
from lab.verified_prefix import (binary_fraction, round_dyadic, integer_cdf_counts,
    integer_cdf_index, uniform_integer, _rational_phase, _word, _bits)
from lab.verified_rejection import _positive, _accuracy_bits, conservative_threshold


def integer_norm(vector):
    return sum(re*re+im*im for re,im in vector)


def integer_matvec(matrix, vector):
    return tuple((sum(ar*br-ai*bi for (ar,ai),(br,bi) in zip(row,vector)),
                  sum(ar*bi+ai*br for (ar,ai),(br,bi) in zip(row,vector)))
                 for row in matrix)


def canonical_quantize(vector, bits):
    """Round after max-real-coordinate scaling; return (vector, was_zero).

    For a zero child, return a declared basis vector. Such a child has exact
    zero sampling weight and can only be visited by forced diagnostics.
    """
    _bits(bits)
    if (bits < 1 or not isinstance(vector,(tuple,list)) or not 1 <= len(vector) <= 64
            or any(not isinstance(z,(tuple,list)) or len(z)!=2 or
                   any(type(x) is not int for x in z) for z in vector)):
        raise ValueError("require a nonempty bounded Gaussian-integer vector and positive grid bits")
    scale = max(abs(x) for z in vector for x in z)
    if scale == 0:
        return tuple((int(i==0)*(1 << bits),0) for i in range(len(vector))),True
    result = tuple(tuple(round_dyadic(Fraction(x,scale),bits) for x in z) for z in vector)
    return result,False


def integer_overlap_numerator(left, right):
    re = sum(ar*br+ai*bi for (ar,ai),(br,bi) in zip(left,right))
    im = sum(ar*bi-ai*br for (ar,ai),(br,bi) in zip(left,right))
    return re*re+im*im


class VerifiedQuantizedReverseWork(VerifiedReverseWork):
    def plan(self, target_tv):
        delta = _positive(target_tv)
        if delta >= 1:
            raise ValueError("target TV must be below one")
        b,t = self.circuit.b,self.circuit.width
        p = 4
        while True:
            eta = Fraction(1,1 << p)
            nu,kappa,mu = 2*b*eta,4*b*eta,4*b*eta
            operator = 2*(2+nu)*nu
            if nu < 1 and b*(t*(operator+kappa)+mu) <= delta/2:
                break
            p += 1
            if p >= 4096:
                raise ValueError("grid accuracy exceeds bounded planner")
        L = _accuracy_bits(delta/(8*b*t)) if t else 1
        La = _accuracy_bits(delta/(4*b))
        local = operator+kappa+Fraction(2,1 << L)
        accepted_error = t*local+mu+Fraction(1,1 << La)
        beta = b*accepted_error
        child_coordinate_bound = 2*b*((1 << (p+1))+1)*(1 << p)
        return dict(target_tv=delta,total_tv_upper_bound=beta,
            grid_bits=p,grid_spacing=eta,operator_grid_bits=p+1,
            operator_coordinate_error=eta/2,stacked_operator_error_upper_bound=nu,
            operator_cq_error_upper_bound=operator,state_trace_error_upper_bound=kappa,
            initial_state_trace_error_upper_bound=mu,local_cq_error_upper_bound=local,
            proposal_tv_upper_bound=min(Fraction(1),t*local/2),
            accepted_measure_l1_upper_bound=accepted_error,
            categorical_random_bits=L,acceptance_random_bits=La,envelope_upper_bound=b,
            success_probability_lower_bound=(1-beta)/b,
            expected_attempts_upper_bound=Fraction(b)/(1-beta),
            state_coordinate_magnitude_upper_bound=1 << p,
            initial_coordinate_magnitude_upper_bound=(1 << (p+1))+1,
            initial_coordinate_bits_upper_bound=p+2,
            child_coordinate_magnitude_upper_bound=child_coordinate_bound,
            state_coordinate_bits_upper_bound=p+1,
            child_coordinate_bits_upper_bound=child_coordinate_bound.bit_length(),
            machine_mantissa_bits=None)

    def cursor(self, initial_sector, boundary_j, *, target_tv=Fraction(1,1_000_000),
               initial_precision=None):
        return _QuantizedCursor(self,initial_sector,boundary_j,self.plan(target_tv),initial_precision)

    def attempt(self, initial_sector, boundary_j, *, rng=None, output=None,
                target_tv=Fraction(1,1_000_000),initial_precision=None):
        if (rng is None) == (output is None):
            raise ValueError("supply rng or forced output, exclusively")
        t,b = self.circuit.width,self.circuit.b
        if output is not None and (type(output) is not int or not 0 <= output < 1 << t):
            raise ValueError("forced output outside width")
        cursor = self.cursor(initial_sector,boundary_j,target_tv=target_tv,
                             initial_precision=initial_precision)
        plan = cursor.plan
        probability = Fraction(1)
        for depth in range(t):
            weights = cursor.weights()
            L = plan["categorical_random_bits"]
            if output is None:
                bit = integer_cdf_index(weights,_word(rng,L),L)
            else:
                bit = (output >> depth)&1
                probability *= Fraction(integer_cdf_counts(weights,L)[bit],1 << L)
            cursor.advance(bit)
        threshold,N,D = cursor.acceptance(plan["acceptance_random_bits"])
        result = dict(**plan,output=cursor.output,boundary_j=boundary_j,
            initial_coarse_sector=initial_sector,
            final_coarse_sector=self.builder._sector_at(initial_sector,t),
            acceptance_threshold=threshold,
            acceptance_probability=Fraction(threshold,1 << plan["acceptance_random_bits"]),
            acceptance_integer_numerator=N,acceptance_integer_denominator=D,
            operator_enclosure_evaluations=cursor.oracle_evaluations,
            branch_matrix_pair_constructions=cursor.branch_builds,
            reverse_vector_matvecs=cursor.matvecs,vector_compressions=cursor.compressions,
            refinement_retries=cursor.refinements,max_working_precision=cursor.precision,
            max_state_coordinate_bits=cursor.max_state_bits,
            max_initial_coordinate_bits=cursor.initial_coordinate_bits,
            max_child_coordinate_bits=cursor.max_child_bits,
            max_norm_weight_bits=cursor.max_weight_bits,
            max_acceptance_integer_bits=max(N.bit_length(),D.bit_length()),
            diagnostic_zero_child_fallbacks=cursor.zero_children,
            resampled_blocks=t,forward_steps=0,replayed_vector_steps=0,
            stored_forward_matrix_count=0,stored_branch_matrix_count=0,
            retained_work_vector_count_upper_bound=4,
            working_scalar_count_upper_bound=24*b*b+8*b,work_block_size=b,
            working_scalar_units="complex-coordinate slots; two integers per Gaussian-integer slot; excludes variable bit storage and backend scratch",
            orbit_or_output_tables=False,requires_independent_unbiased_bits=True,
            certifies_existing_float_sampler=False,precision_cap=None,attempt_cap=None,
            method="verified locally quantized reverse-vector boundary rejection",
            vector_normalization="exact projective scaling and declared dyadic rounding")
        if output is not None:
            result["proposal_path_probability"] = probability
        return result

    def sample(self,rng,*,target_tv=Fraction(1,1_000_000)):
        plan = self.plan(target_tv)
        initial = uniform_integer(self.circuit.sectors,rng)
        totals = {key:0 for key in ("operator_enclosure_evaluations","branch_matrix_pair_constructions",
            "reverse_vector_matvecs","vector_compressions","refinement_retries","resampled_blocks",
            "diagnostic_zero_child_fallbacks")}
        maxima = {key:0 for key in ("max_working_precision","max_state_coordinate_bits","max_initial_coordinate_bits",
            "max_child_coordinate_bits","max_norm_weight_bits","max_acceptance_integer_bits")}
        attempts = 0
        while True:
            attempts += 1
            j = uniform_integer(self.circuit.b,rng)
            proposed = self.attempt(initial,j,rng=rng,target_tv=target_tv)
            for key in totals:
                totals[key] += proposed[key]
            for key in maxima:
                maxima[key] = max(maxima[key],proposed[key])
            if _word(rng,plan["acceptance_random_bits"]) < proposed["acceptance_threshold"]:
                return dict(**plan,**totals,**maxima,attempts=attempts,
                    output=proposed["output"],initial_coarse_sector=initial,
                    final_coarse_sector=proposed["final_coarse_sector"],
                    forward_steps=0,replayed_vector_steps=0,
                    stored_forward_matrix_count=0,stored_branch_matrix_count=0,
                    retained_work_vector_count_upper_bound=4,
                    working_scalar_count_upper_bound=24*self.circuit.b**2+8*self.circuit.b,
                    working_scalar_units=proposed["working_scalar_units"],
                    work_block_size=self.circuit.b,orbit_or_output_tables=False,
                    requires_independent_unbiased_bits=True,certifies_existing_float_sampler=False,
                    precision_cap=None,attempt_cap=None,
                    method="verified locally quantized reverse boundary rejection")


class _QuantizedCursor:
    def __init__(self,owner,initial,boundary,plan,initial_precision):
        if type(initial) is not int or not 0 <= initial < owner.circuit.sectors:
            raise ValueError("initial sector outside circuit")
        if type(boundary) is not int or not 0 <= boundary < owner.circuit.b:
            raise ValueError("boundary outside work block")
        if initial_precision is not None and (type(initial_precision) is not int or initial_precision < 16):
            raise ValueError("initial precision must be None or exact integer>=16, never a cap")
        self.owner,self.initial,self.boundary,self.plan = owner,initial,boundary,plan
        self.precision = max(64,plan["grid_bits"]+16) if initial_precision is None else initial_precision
        self.output,self.measured,self._pending = 0,0,None
        self.oracle_evaluations,self.branch_builds,self.refinements = 0,0,0
        self.matvecs,self.compressions,self.zero_children = 0,0,0
        self.max_state_bits,self.max_child_bits,self.max_weight_bits = plan["grid_bits"]+1,0,0
        self.vector = tuple((int(i==boundary)*(1 << plan["grid_bits"]),0)
                            for i in range(owner.circuit.b))
        self.initial_vector = tuple(row[0] for row in self._round_matrices(self._initial)[0])
        self.initial_coordinate_bits = max(abs(x).bit_length() for z in self.initial_vector for x in z)
        if integer_norm(self.initial_vector) == 0:
            raise ArithmeticError("initial-state rounding violated its nonzero bound")

    def _initial(self):
        from flint import acb_mat
        b = self.owner.builder
        return (b._background(0,b._identity())*acb_mat(
            [[int(i==0)] for i in range(self.owner.circuit.b)]),)

    def _operators(self):
        b = self.owner.builder
        i = self.owner.circuit.width-1-self.measured
        pair = b._branch(i,b._sector_at(self.initial,i),b._identity())
        self.branch_builds += 1
        phase = _rational_phase(-self.output,1 << (self.measured+1))
        return tuple(_adjoint(K) for K in _qft_branch_operators(*pair,phase))

    def _round_matrices(self,factory):
        from flint import ctx
        tolerance = self.plan["grid_spacing"]/4
        while True:
            with ctx.workprec(self.precision):
                matrices = factory()
                self.oracle_evaluations += 1
                coordinates = [x for A in matrices for row in range(A.nrows())
                               for col in range(A.ncols()) for x in (A[row,col].real,A[row,col].imag)]
                if all(x.is_finite() and binary_fraction(x.rad().upper()) <= tolerance for x in coordinates):
                    p = self.plan["operator_grid_bits"]
                    return tuple(tuple(tuple((round_dyadic(binary_fraction(A[i,j].real.mid()),p),
                                              round_dyadic(binary_fraction(A[i,j].imag.mid()),p))
                                             for j in range(A.ncols())) for i in range(A.nrows())) for A in matrices)
            self.precision *= 2
            self.refinements += 1
            # Past rounded operators and current integer state NEVER change.

    def weights(self):
        if self.measured >= self.owner.circuit.width:
            raise ValueError("all output bits already measured")
        if self._pending is not None:
            return self._pending[0]
        operators = self._round_matrices(self._operators)
        children = tuple(integer_matvec(A,self.vector) for A in operators)
        weights = tuple(integer_norm(v) for v in children)
        self.matvecs += 2
        if sum(weights) == 0:
            raise ArithmeticError("rounded stacked instrument violated nu<1")
        self.max_child_bits = max(self.max_child_bits,max(abs(x).bit_length() for v in children for z in v for x in z))
        self.max_weight_bits = max(self.max_weight_bits,max(x.bit_length() for x in weights))
        self._pending = weights,children
        return weights

    def advance(self,bit):
        if type(bit) is not int or bit not in (0,1) or self._pending is None:
            raise ValueError("query weights then advance by an exact bit")
        self.vector,zero = canonical_quantize(self._pending[1][bit],self.plan["grid_bits"])
        self.compressions += 1
        self.zero_children += zero
        self.max_state_bits = max(self.max_state_bits,max(abs(x).bit_length() for z in self.vector for x in z))
        self.output |= bit << self.measured
        self.measured += 1
        self._pending = None

    def acceptance(self,random_bits):
        if self.measured != self.owner.circuit.width:
            raise ValueError("terminal acceptance requires all output bits")
        _bits(random_bits)
        N = integer_overlap_numerator(self.initial_vector,self.vector)
        D = integer_norm(self.initial_vector)*integer_norm(self.vector)
        if not 0 <= N <= D or D == 0:
            raise ArithmeticError("invalid normalized integer overlap")
        return conservative_threshold(N,D,random_bits),N,D
