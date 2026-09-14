"""Verified density/effect arithmetic for C56's sequential finite-work identity.

Exact-input b=2/3 deterministic component of C59, not an arbitrary matrix or
coherent-history sampler. Forward densities and backward effects remain
UNNORMALIZED with their physical scale. Only nonnegative dyadic CHILD
PROBABILITIES are normalized by the integer kernel. No QR, effect rescaling,
amplitude cutoff, prefix/output table, or redraw on precision refinement.

The complete-law bound uses unnormalized prefix masses, not relative branch
accuracy. Optional python-flint==0.9.0, exact integers, independent unbiased
bits and nonconcurrent precision context are the C61 trust assumptions.
Default and optional lower starting precisions are uncapped. Resource errors
do not certify a success-conditioned law. See C64 for proof/cost boundary.
"""
from __future__ import annotations

from fractions import Fraction
from weakref import proxy

from lab.semiclassical import _qft_branch_operators, _density_forward_step
from lab.sampling_error import plan_prefix_mass_accuracy
from lab.verified_prefix import (VerifiedReflectionCircuit, binary_fraction,
    round_dyadic, integer_cdf_index, integer_cdf_counts, uniform_integer,
    _rational_phase, _rotation_coefficients, _apply_background, _word, _bits)


def _adjoint(A):
    return A.conjugate().transpose()


class VerifiedFiniteWork:
    """Sample one normalized history component, with its initial sector uniform.

    mask uses the original sorted reflection insertions. Reflection angles
    weight the OUTER history selection, not this normalized J-route component.
    Both initial and final sectors can be reported; within-sector work is traced.
    """
    def __init__(self, circuit, mask, *, checkpoint_spacing=None):
        if (type(circuit) is not VerifiedReflectionCircuit or type(mask) is not int
                or not 0 <= mask < 1 << len(circuit.reflections)):
            raise ValueError("require exact-input circuit and valid history mask")
        if checkpoint_spacing is not None and (type(checkpoint_spacing) is not int
                                               or not 1 <= checkpoint_spacing <= 63):
            raise ValueError("checkpoint spacing must be None or an exact integer in [1,63]")
        self.checkpoint_spacing = checkpoint_spacing
        self.circuit = VerifiedReflectionCircuit(circuit.period,circuit.width,
            dict(circuit.backgrounds),dict(circuit.reflections),
            enclosure_mode=circuit.enclosure_mode, block_size=circuit.b)
        self.mask = mask
        self.routes = {s:q for i,(s,(q,_)) in enumerate(circuit.reflections.items()) if mask & (1 << i)}

    def _identity(self):
        from flint import acb_mat
        c = self.circuit
        return acb_mat([[int(i == j) for j in range(c.b)] for i in range(c.b)])

    def _background(self,s,I):
        from flint import acb_mat
        c = self.circuit
        if s not in c.backgrounds:
            return I
        axis,theta = c.backgrounds[s]
        co,si = _rotation_coefficients(theta)
        columns = [_apply_background([I[i,j] for i in range(c.b)],axis,co,si)
                   for j in range(c.b)]
        return acb_mat([[columns[j][i] for j in range(c.b)] for i in range(c.b)])

    def _initial_density(self,I):
        from flint import acb_mat
        initial = self._background(0,I)*acb_mat([[int(i == 0)] for i in range(self.circuit.b)])
        return initial*_adjoint(initial)

    def _branch(self,i,sector,I):
        """Shared exact branch construction; caller supplies its pre-control sector."""
        from flint import acb_mat
        c = self.circuit
        shift = acb_mat(c.b,c.b)
        for p in range(c.b):
            wraps,target = divmod(p+(1 << i),c.b)
            shift[target,p] = _rational_phase(sector*wraps,c.sectors)
        background = self._background(i+1,I)
        return background,background*shift

    def _sector_at(self,initial,stop):
        current = initial
        for insertion,q in self.routes.items():
            if insertion <= stop:
                current = (-current-q) % self.circuit.sectors
        return current

    def _build(self, initial_sector):
        """Full reference build at current precision, preserving operation order."""
        c = self.circuit
        I = self._identity()
        rho = self._initial_density(I)
        states,branches = [rho],[]
        current = initial_sector
        if 0 in self.routes:
            current = (-current-self.routes[0]) % c.sectors
        for i in range(c.width):
            B0,B1 = self._branch(i,current,I)
            branches.append((B0,B1))
            rho = _density_forward_step(B0,B1,rho,multiply=lambda A,B:A*B,adjoint=_adjoint)
            states.append(rho)
            if i+1 in self.routes:
                current = (-current-self.routes[i+1]) % c.sectors
        return branches,states,current

    def cursor(self, initial_sector, *, accuracy_bits, initial_precision=None):
        return _FiniteWorkCursor(self,initial_sector,accuracy_bits,initial_precision)

    def path(self, initial_sector, *, rng=None, output=None,
             target_tv=Fraction(1,1_000_000), initial_precision=None):
        """One sampled path, or exact finite-kernel probability of a forced output.

        Forced output is a diagnostic, not a sampler conditioned on success.
        It follows the SAME per-prefix precision schedule as sampled execution.
        """
        if (rng is None) == (output is None):
            raise ValueError("supply rng or a forced output, exclusively")
        t = self.circuit.width
        if output is not None and (type(output) is not int or not 0 <= output < 1 << t):
            raise ValueError("forced output outside circuit width")
        plan = plan_prefix_mass_accuracy(t,target_tv)
        cursor = self.cursor(initial_sector,accuracy_bits=plan["weight_accuracy_bits"],
                             initial_precision=initial_precision)
        probability = Fraction(1)
        for step in range(t):
            weights = cursor.weights()
            if output is None:
                bit = integer_cdf_index(weights,_word(rng,plan["categorical_random_bits"]),
                                         plan["categorical_random_bits"])
            else:
                bit = (output >> step)&1
                count = integer_cdf_counts(weights,plan["categorical_random_bits"])[bit]
                probability *= Fraction(count,1 << plan["categorical_random_bits"])
            cursor.advance(bit)
        result = dict(**plan,output=cursor.output,initial_coarse_sector=initial_sector,
            final_coarse_sector=cursor.final_sector,forward_builds=cursor.builds,
            forward_steps=cursor.forward_steps,backward_child_evaluations=cursor.evaluations,
            replayed_effect_steps=cursor.replayed,refinement_retries=cursor.builds-1,
            max_working_precision=cursor.precision,zero_approximate_block_fallbacks=cursor.zero_blocks,
            resampled_blocks=t,prefix_vector_evaluations=0,orbit_or_output_tables=False,
            stored_matrix_count_upper_bound=3*t+16,history_count_upper_bound=1,
            work_block_size=self.circuit.b,
            stored_matrix_scalar_count_upper_bound=(3*t+16)*self.circuit.b**2,
            method="verified finite-work unnormalized density/effect sampler",
            enclosure_mode="rectangular density/effect",
            requires_independent_unbiased_bits=True,certifies_existing_float_sampler=False,
            precision_cap=None,initial_precision=initial_precision,
            effect_normalization=False,forward_normalization=False)
        result.update(checkpoint_spacing=self.checkpoint_spacing,
            recomputed_forward_steps=cursor.recomputed_forward_steps,
            branch_matrix_pair_constructions=cursor.branch_builds,
            peak_retained_forward_matrix_count=cursor.peak_forward_matrices,
            stored_branch_matrix_count=2*t if self.checkpoint_spacing is None else 0,
            stored_branch_matrices_are_persistent=self.checkpoint_spacing is None)
        if self.checkpoint_spacing is not None:
            k = self.checkpoint_spacing
            result["stored_matrix_count_upper_bound"] = (t+k-1)//k + min(t,k) + 24
            result["stored_matrix_scalar_count_upper_bound"] = self.circuit.b**2*result["stored_matrix_count_upper_bound"]
            result["method"] = "verified finite-work checkpoint/recompute sampler"
        if output is not None:
            result["conditional_path_probability"] = probability
        return result

    def sample(self,rng,*,target_tv=Fraction(1,1_000_000)):
        # Validate target before consuming random bits.
        plan_prefix_mass_accuracy(self.circuit.width,target_tv)
        return self.path(uniform_integer(self.circuit.sectors,rng),rng=rng,target_tv=target_tv)


class _FiniteWorkCursor:
    def __init__(self, owner, initial_sector, accuracy_bits, initial_precision):
        _bits(accuracy_bits)
        if type(initial_sector) is not int or not 0 <= initial_sector < owner.circuit.sectors:
            raise ValueError("initial sector outside circuit")
        if initial_precision is not None and (type(initial_precision) is not int or initial_precision < 16):
            raise ValueError("initial precision must be None or integer >=16, never a cap")
        self.owner,self.initial_sector,self.bits = owner,initial_sector,accuracy_bits
        self.precision = max(64,accuracy_bits+16) if initial_precision is None else initial_precision
        self.measured,self.output,self.builds,self.replayed,self.evaluations,self.zero_blocks = 0,0,0,0,0,0
        self.forward_steps,self.recomputed_forward_steps,self.branch_builds,self.peak_forward_matrices = 0,0,0,0
        self._pending = None
        self._rebuild()

    def _rebuild(self):
        from flint import acb_mat,ctx
        # Discard OWN obsolete checkpoints before rebuilding; otherwise peak
        # live storage includes old and new O(t) arrays during assignment.
        self.branches,self.states,self._pending = [],[],None
        with ctx.workprec(self.precision):
            if self.owner.checkpoint_spacing is None:
                self.branches,self.states,self.final_sector = self.owner._build(self.initial_sector)
                self.forward_steps += self.owner.circuit.width
                self.branch_builds += self.owner.circuit.width
                self.peak_forward_matrices = max(self.peak_forward_matrices,len(self.states))
            else:
                self.branches = _OnDemandBranches(self)
                self.states = _CheckpointStates(self)
                self.final_sector = self.owner._sector_at(self.initial_sector,self.owner.circuit.width)
            b = self.owner.circuit.b
            self.effect = acb_mat([[int(i == j) for j in range(b)] for i in range(b)])
            for step in range(self.measured):
                phase = _rational_phase(-(self.output % (1 << step)),1 << (step+1))
                K = _qft_branch_operators(*self.branches[self.owner.circuit.width-1-step],phase)[(self.output >> step)&1]
                self.effect = _adjoint(K)*self.effect*K
                self.replayed += 1
        self.builds += 1
        self._pending = None

    def weights(self):
        """Two dyadic masses with uniform abs error<=2^-p; no small division."""
        if self.measured >= self.owner.circuit.width:
            raise ValueError("all outputs already measured")
        if self._pending is not None:
            return self._pending[0]
        from flint import ctx
        tolerance = Fraction(1,1 << (self.bits+1))
        while True:
            with ctx.workprec(self.precision):
                i = self.owner.circuit.width-1-self.measured
                phase = _rational_phase(-self.output,1 << (self.measured+1))
                children = tuple(_adjoint(K)*self.effect*K
                    for K in _qft_branch_operators(*self.branches[i],phase))
                values = tuple((child*self.states[i]).trace() for child in children)
                self.evaluations += 1
                if all(z.is_finite() for z in values):
                    if any(not z.imag.contains(0) or z.real.upper() < 0 for z in values):
                        raise ArithmeticError("invalid nonnegative real prefix-mass enclosure")
                    if max(binary_fraction(z.real.rad().upper()) for z in values) <= tolerance:
                        weights = tuple(max(0,round_dyadic(binary_fraction(z.real.mid()),self.bits)) for z in values)
                        self.zero_blocks += not any(weights)
                        self._pending = weights,children
                        return weights
            self.precision *= 2
            self._rebuild()

    def advance(self,bit):
        if type(bit) is not int or bit not in (0,1) or self._pending is None:
            raise ValueError("query weights then advance by one exact bit")
        self.effect = self._pending[1][bit]
        self.output |= bit << self.measured
        self.measured += 1
        self._pending = None


class _OnDemandBranches:
    """No persistent branch matrices or coarse-sector table."""
    def __init__(self,cursor):
        # Avoid cursor -> source -> cursor cycles retaining native matrices
        # across completed paths until a later cyclic-GC collection.
        self.cursor = proxy(cursor)

    def __getitem__(self,i):
        c = self.cursor
        if type(i) is not int or not 0 <= i < c.owner.circuit.width:
            raise IndexError("branch index outside width")
        c.branch_builds += 1
        return c.owner._branch(i,c.owner._sector_at(c.initial_sector,i),c.owner._identity())


class _CheckpointStates:
    """Store block starts and reconstruct one descending block at a time.

    Same arithmetic and precision as the full forward recurrence; no state
    approximation. Internal cursor access is descending, so each interior
    state is reconstructed at most once per precision level. Rebuilding
    starts anew and preserves the already measured output bits via E replay.
    """
    def __init__(self,cursor):
        self.cursor = proxy(cursor)
        self.spacing = cursor.owner.checkpoint_spacing
        self.boundaries,self.local = {},[]
        self.start = None
        rho = cursor.owner._initial_density(cursor.owner._identity())
        t = cursor.owner.circuit.width
        for i in range(t):
            if i % self.spacing == 0:
                self.boundaries[i] = rho
            B0,B1 = cursor.branches[i]
            rho = _density_forward_step(B0,B1,rho,multiply=lambda A,B:A*B,adjoint=_adjoint)
            cursor.forward_steps += 1
            cursor.peak_forward_matrices = max(cursor.peak_forward_matrices,len(self.boundaries)+1)

    def __getitem__(self,i):
        c = self.cursor
        t = c.owner.circuit.width
        if type(i) is not int or not 0 <= i < t:
            raise IndexError("state index outside width")
        start = i-i % self.spacing
        if start != self.start:
            # Release obsolete OWN local arrays before reconstructing.
            self.local = []
            self.start = start
            rho = self.boundaries[start]
            self.local.append(rho)
            for j in range(start,min(start+self.spacing,t)-1):
                B0,B1 = c.branches[j]
                rho = _density_forward_step(B0,B1,rho,multiply=lambda A,B:A*B,adjoint=_adjoint)
                self.local.append(rho)
                c.forward_steps += 1
                c.recomputed_forward_steps += 1
            # The first local entry aliases its retained boundary matrix.
            c.peak_forward_matrices = max(c.peak_forward_matrices,
                len(self.boundaries)+len(self.local)-1)
        return self.local[i-start]


class VerifiedScalarWork(VerifiedFiniteWork):
    """Stronger b=2 component baseline: all i>=1 shifts are scalar on p.

    All W_s, s>=1, commute to the final traced work subsystem, in their
    original mutual order. W0 remains before the only nonscalar control.
    Conditional binary weights are normalized ANALYTICALLY, not by a tiny
    numerical prefix mass. No density/effect checkpoints are needed.
    """
    def __init__(self, circuit, mask):
        if type(circuit) is not VerifiedReflectionCircuit or circuit.b != 2:
            raise ValueError("scalar component requires block_size=2")
        super().__init__(circuit, mask)

    def _sectors(self,initial):
        current = initial
        if 0 in self.routes:
            current = (-current-self.routes[0]) % self.circuit.sectors
        labels = []
        for i in range(self.circuit.width):
            labels.append(current)
            if i+1 in self.routes:
                current = (-current-self.routes[i+1]) % self.circuit.sectors
        return labels,current

    def conditional_weights(self,initial_sector,measured,output,*,accuracy_bits):
        """Fixed-label normalized branch weights, verified then dyadic-rounded."""
        c = self.circuit
        _bits(accuracy_bits)
        if (any(type(v) is not int for v in (initial_sector,measured,output))
                or not 0 <= initial_sector < c.sectors or not 0 <= measured < c.width
                or not 0 <= output < 1 << measured):
            raise ValueError("invalid scalar conditional label")
        # Compute the one needed pre-control sector, using only O(k) routes;
        # do not rebuild an O(t) label list at every bit.
        i = c.width-1-measured
        alpha = initial_sector
        for s,q in self.routes.items():
            if s <= i:
                alpha = (-alpha-q) % c.sectors
        from flint import acb,ctx
        precision,retries = max(64,accuracy_bits+16),0
        while True:
            with ctx.workprec(precision):
                phase = _rational_phase(-output,1 << (measured+1))
                if i:
                    z = phase*_rational_phase(alpha*(1 << (i-1)),c.sectors)
                    vectors = (((1+z)/2,),((1-z)/2,))
                else:
                    psi = (acb(1),acb(0))
                    if 0 in c.backgrounds:
                        axis,theta = c.backgrounds[0]
                        co,si = _rotation_coefficients(theta)
                        psi = (co,si) if axis == "x" else (co+si,acb(0))
                    shifted = (_rational_phase(alpha,c.sectors)*psi[1],psi[0])
                    vectors = tuple(tuple((a+sign*phase*b)/2 for a,b in zip(psi,shifted)) for sign in (1,-1))
                values = tuple(sum(z.real*z.real+z.imag*z.imag for z in vector) for vector in vectors)
                if all(z.is_finite() for z in values):
                    if max(binary_fraction(z.rad().upper()) for z in values) <= Fraction(1,1 << (accuracy_bits+1)):
                        weights = tuple(max(0,round_dyadic(binary_fraction(z.mid()),accuracy_bits)) for z in values)
                        return weights,precision,retries
            precision *= 2
            retries += 1

    def path(self,initial_sector,*,rng=None,output=None,target_tv=Fraction(1,1_000_000),initial_precision=None):
        from lab.sampling_error import plan_conditional_weight_accuracy
        if initial_precision is not None:
            raise ValueError("scalar mode has a fixed per-label starting precision")
        c = self.circuit
        if (type(initial_sector) is not int or not 0 <= initial_sector < c.sectors
                or (rng is None) == (output is None)
                or output is not None and (type(output) is not int or not 0 <= output < 1 << c.width)):
            raise ValueError("invalid scalar sample/forced-output request")
        plan = plan_conditional_weight_accuracy(c.width,target_tv)
        prefix,max_precision,retries,zero_blocks = 0,0,0,0
        probability = Fraction(1)
        for j in range(c.width):
            weights,P,R = self.conditional_weights(initial_sector,j,prefix,accuracy_bits=plan["weight_accuracy_bits"])
            max_precision,retries = max(max_precision,P),retries+R
            zero_blocks += not any(weights)
            L = plan["categorical_random_bits"]
            if output is None:
                bit = integer_cdf_index(weights,_word(rng,L),L)
            else:
                bit = (output >> j)&1
                probability *= Fraction(integer_cdf_counts(weights,L)[bit],1 << L)
            prefix |= bit << j
        _,final = self._sectors(initial_sector)
        result = dict(**plan,output=prefix,initial_coarse_sector=initial_sector,final_coarse_sector=final,
            max_working_precision=max_precision,refinement_retries=retries,
            zero_approximate_block_fallbacks=zero_blocks,resampled_blocks=c.width,
            prefix_vector_evaluations=0,scalar_conditional_evaluations=c.width+retries,
            method="verified b=2 scalar routed-component sampler",enclosure_mode="rectangular scalar/vector",
            requires_independent_unbiased_bits=True,certifies_existing_float_sampler=False,
            precision_cap=None,orbit_or_output_tables=False,late_backgrounds_removed=True)
        if output is not None:
            result["conditional_path_probability"] = probability
        return result

    def sample(self,rng,*,target_tv=Fraction(1,1_000_000)):
        from lab.sampling_error import plan_conditional_weight_accuracy
        plan_conditional_weight_accuracy(self.circuit.width,target_tv)
        return self.path(uniform_integer(self.circuit.sectors,rng),rng=rng,target_tv=target_tv)
