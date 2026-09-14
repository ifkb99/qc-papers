"""Sparse backward sector vectors for the supplied coherent-reflection family.

Histories merge by adding complex vectors at equal reached sector labels.
For two fixed reflection labels the support is bounded by the word ball of
two involutions; no orbit or history table is needed. The shared work shift
and QFT branch helpers preserve original gate order. Terminal rejection tests
the coherently uniform INITIAL sector boundary, not a dephased mixture.

This is an opt-in complex128 diagnostic implementation, NOT a finite-TV
certificate. Mathematical real-arithmetic identities do not certify float
rounding, RNG discrepancy, norm underflow or a resource-conditioned output.
Input gates/order/index are supplied; no generic circuit propagator is added.
"""
from __future__ import annotations

import math
import numpy as np

from lab.coherent_routes import CoherentReflectionCircuit, CoherentReflectionInput
from lab.fourier_sampling import _ints,unit_phase
from lab.semiclassical import _qft_branch_operators


def route_support_bound(sectors,qs):
    """Global bound over ALL final sectors, including intermediate supports."""
    if (not _ints(sectors) or sectors < 1 or not isinstance(qs,(tuple,list))
            or len(qs)>64 or any(not _ints(q) for q in qs)):
        raise ValueError("require a positive sector count and at most 64 integer route labels")
    qs = tuple(int(q)%int(sectors) for q in qs)
    kinds = len(set(qs))
    k = len(qs)
    # Two orientations of the integer L1 ball in d-1 coefficient dimensions.
    # Choose j nonzero coordinates, their signs, then positive magnitudes
    # whose sum is <=k. Collisions in Z_M only reduce the cover.
    cover = (1 if kinds==0 else 2 if kinds==1 else 1+2*k if kinds==2 else
             2*sum((1 << j)*math.comb(kinds-1,j)*math.comb(k,j)
                   for j in range(min(kinds-1,k)+1)))
    return min(int(sectors),1 << k,cover)


def sparse_norm(vector):
    return math.fsum(float(np.vdot(v,v).real) for v in vector.values())


def word_support_plan(sectors, qs, *, max_offset_entries=65_536):
    """Global envelope from a CHRONOLOGICAL reflection word, not one sector.

    Count affine maps alpha->alpha+E or -alpha+O, without history expansion.
    Return scalar metadata only; temporary integer sets are released at setup.
    The explicit cap covers old AND new offset entries before each update;
    entries are not native bytes. Reflected/copied counts are set operations,
    not bit-runtime or Python allocator certificates.
    """
    alphabet = route_support_bound(sectors,qs)  # shared bounded validation
    if not _ints(max_offset_entries) or not 1<=max_offset_entries<=65_536:
        raise ValueError("word setup requires an integer live-offset cap in [1,65536]")
    sectors = int(sectors)
    even,odd = {0},set()
    copied = reflected = 0
    peak = 1
    for q in reversed(tuple(int(q)%sectors for q in qs)):
        old_size = len(even)+len(odd)
        # Old sets plus at most two old_size-entry new sets. No temporary
        # Cartesian products or new-image sets are materialized.
        if 3*old_size>max_offset_entries:
            raise MemoryError("word envelope live-offset preflight cap exceeded")
        next_even,next_odd = set(even),set(odd)
        copied += old_size
        for o in odd:
            next_even.add((-o-q)%sectors)
        for e in even:
            next_odd.add((-e-q)%sectors)
        reflected += old_size
        peak = max(peak,old_size+len(next_even)+len(next_odd))
        even,odd = next_even,next_odd
    word = min(sectors,len(even)+len(odd))
    return dict(support_mode="word",support_bound=min(alphabet,word),
        alphabet_support_bound=alphabet,word_support_bound=word,
        word_even_offset_count=len(even),word_odd_offset_count=len(odd),
        word_setup_steps=len(qs),word_setup_copied_offsets=copied,
        word_setup_reflected_offsets=reflected,word_setup_peak_live_offset_entries=peak,
        word_setup_max_offset_entries=int(max_offset_entries),
        word_retained_offset_entries=0)


def _support_plan(sectors,qs,mode):
    if not isinstance(mode,str) or mode not in ("alphabet","word"):
        raise ValueError("support_mode must be 'alphabet' or 'word'")
    if mode=="word":
        return word_support_plan(sectors,qs)
    bound = route_support_bound(sectors,qs)
    return dict(support_mode="alphabet",support_bound=bound,
        alphabet_support_bound=bound,word_support_bound=None,
        word_even_offset_count=0,word_odd_offset_count=0,
        word_setup_steps=0,word_setup_copied_offsets=0,
        word_setup_reflected_offsets=0,word_setup_peak_live_offset_entries=0,
        word_setup_max_offset_entries=0,word_retained_offset_entries=0)


class SparseCoherentReverse:
    """One sparse pure reverse trajectory, with a global boundary envelope.

    Legacy input retains k<=8. Explicit CoherentReflectionInput allows k<=64
    with a 16-MiB arithmetic-payload and 65,536-sector cover preflight cap.
    Fixed d distinct labels give polynomial support in k. Zero amplitudes are
    not tolerance-pruned; dictionary keys track an algebraic support cover.
    Opt-in support_mode='word' tightens the GLOBAL envelope by one bounded
    integer-set setup; the default alphabet envelope and sampling loop remain.
    """
    def __init__(self,circuit,*,support_mode="alphabet"):
        if not isinstance(circuit,CoherentReflectionInput):
            raise ValueError("require a supplied coherent-reflection input")
        # Public input dictionaries may have been reinserted out of order.
        # Match the snapshot constructor's original-time sort exactly.
        self.support_plan = _support_plan(circuit.sectors,
            tuple(circuit.reflections[s][0] for s in sorted(circuit.reflections)),support_mode)
        self.support_bound = self.support_plan["support_bound"]
        structural = not isinstance(circuit,CoherentReflectionCircuit)
        slots = (24*circuit.b*circuit.b+8*self.support_bound*circuit.b
                 +(len(circuit.background.defects)+1)*circuit.b*circuit.b)
        if (structural or support_mode=="word") and (self.support_bound > 65_536 or 16*slots > (16 << 20)):
            raise MemoryError("structural reverse cover exceeds diagnostic payload budget")
        # Snapshot inputs; do not enumerate _histories or build r-sized gates.
        cls = CoherentReflectionInput if structural else CoherentReflectionCircuit
        self.circuit = cls(circuit.period,circuit.b,circuit.width,
            circuit.background.defects,circuit.reflections)
        self.terms = {s:(q,c.conjugate(),d.conjugate()) for s,q,c,d in self.circuit._terms}

    def _stats(self):
        c = self.circuit
        return dict(**self.support_plan,
            mathematical_mean_attempts=c.b*self.support_bound,
            route_label_count=len(set(q for q,theta in c.reflections.values())),
            coherent_rotation_count=len(c.reflections),
            history_count_upper_bound=1 << len(c.reflections),history_enumerations=0,
            orbit_table_entries=0,output_table_entries=0,full_sector_table_entries=0,
            work_block_size=c.b,
            working_complex_coordinate_slots_upper_bound=24*c.b*c.b+8*self.support_bound*c.b,
            supplied_gate_complex_entries=(len(c.background.defects)+1)*c.b*c.b,
            precision_certificate=False,arithmetic="complex128 diagnostic",
            order_and_index_discovery_included=False,
            method="sparse coherent reverse instrument with global boundary rejection")

    def _reflection_adjoint(self,vector,s,cost):
        if s not in self.terms:
            return vector
        q,c,d = self.terms[s]
        out = {}
        for alpha,v in vector.items():
            target = (-alpha-q)%self.circuit.sectors
            for label,coefficient in ((alpha,c),(target,d)):
                contribution = coefficient*v
                if label in out:
                    out[label] += contribution
                else:
                    out[label] = contribution
        cost["reflection_vector_contributions"] += 2*len(vector)
        cost["peak_sector_count"] = max(cost["peak_sector_count"],len(out))
        if len(out)>self.support_bound:
            raise ArithmeticError("algebraic support bound violated")
        return out

    def _background_adjoint(self,vector,s,cost):
        W = self.circuit.background.defects.get(s,self.circuit.background.identity)
        adjoint = W.conj().T
        cost["work_matvecs"] += len(vector)
        return {alpha:adjoint@v for alpha,v in vector.items()}

    def attempt(self,final_sector,boundary_j,*,output=None,rng=None):
        c = self.circuit
        if (rng is None)==(output is None):
            raise ValueError("supply rng or forced output exclusively")
        c._validate_label(final_sector)
        if not _ints(boundary_j) or not 0<=boundary_j<c.b:
            raise ValueError("fine boundary outside work block")
        if output is not None and (not _ints(output) or not 0<=output<1 << c.width):
            raise ValueError("forced output outside width")
        initial = np.zeros(c.b,dtype=complex)
        initial[int(boundary_j)] = 1.
        vector = {int(final_sector):initial}
        cost = dict(peak_sector_count=1,reflection_vector_contributions=0,
                    work_matvecs=0,qft_matrix_pair_constructions=0,
                    max_relative_completeness_error=0.,forced_zero_prefixes=0)
        prefix,probability = 0,1.
        for depth in range(c.width):
            i = c.width-1-depth
            parent_norm = sparse_norm(vector)
            # Forward B_z = R_(i+1) W_(i+1) T_i^z. Reverse these factors.
            common = self._reflection_adjoint(vector,i+1,cost)
            common = self._background_adjoint(common,i+1,cost)
            children = ({},{})
            phase = unit_phase(-prefix,1 << (depth+1))
            for alpha,v in common.items():
                K = _qft_branch_operators(c.background.identity,
                    c.background.shift(1 << i,alpha),phase)
                for bit in (0,1):
                    children[bit][alpha] = K[bit].conj().T@v
            cost["work_matvecs"] += 2*len(common)
            cost["qft_matrix_pair_constructions"] += len(common)
            weights = tuple(sparse_norm(child) for child in children)
            total = sum(weights)
            if not all(math.isfinite(w) and w>=0 for w in weights):
                raise ArithmeticError("invalid floating child norms")
            if parent_norm>0:
                cost["max_relative_completeness_error"] = max(
                    cost["max_relative_completeness_error"],abs(total/parent_norm-1))
            if total==0:
                if output is None or probability!=0:
                    raise ArithmeticError("zero floating proposal mass; no sampled fallback")
                cost["forced_zero_prefixes"] += 1
                bit = (int(output)>>depth)&1
            else:
                bit = CoherentReflectionCircuit._draw(weights,rng) if output is None else (int(output)>>depth)&1
                probability *= weights[bit]/total
            prefix |= bit << depth
            vector = children[bit]
        # Initial forward boundary is R0 W0 applied to coherent uniform alpha.
        vector = self._reflection_adjoint(vector,0,cost)
        vector = self._background_adjoint(vector,0,cost)
        D = sparse_norm(vector)
        amplitude = complex(math.fsum(v[0].real for v in vector.values()),
                            math.fsum(v[0].imag for v in vector.values()))
        N = float(abs(amplitude)**2)
        incoherent = math.fsum(float(abs(v[0])**2) for v in vector.values())
        acceptance = N/(self.support_bound*D) if D else 0.
        if not math.isfinite(acceptance) or not 0<=acceptance<=1+1e-10:
            raise ArithmeticError("sparse coherent boundary envelope violated numerically")
        acceptance = min(1.,acceptance)
        return dict(**self._stats(),**cost,output=prefix,
            final_coarse_sector=int(final_sector),boundary_j=int(boundary_j),
            proposal_path_probability=probability,acceptance_probability=acceptance,
            accepted_joint_submass=probability*acceptance/(c.sectors*c.b),
            terminal_norm=D,terminal_coherent_numerator=N,
            terminal_incoherent_numerator=incoherent,terminal_coherent_amplitude=amplitude,
            terminal_sector_count=len(vector),reverse_steps=c.width)

    def sample(self,rng,*,max_proposals=100_000):
        if not _ints(max_proposals) or not 1<=max_proposals<=10_000_000:
            raise ValueError("invalid explicit diagnostic rejection cap")
        totals = dict(reflection_vector_contributions=0,work_matvecs=0,
                      qft_matrix_pair_constructions=0,reverse_steps=0)
        peak,drift = 0,0.
        for attempts in range(1,int(max_proposals)+1):
            gamma = int(rng.integers(self.circuit.sectors))
            j = int(rng.integers(self.circuit.b))
            proposed = self.attempt(gamma,j,rng=rng)
            for key in totals:
                totals[key] += proposed[key]
            peak = max(peak,proposed["peak_sector_count"])
            drift = max(drift,proposed["max_relative_completeness_error"])
            if float(rng.random())<proposed["acceptance_probability"]:
                return dict(**self._stats(),**totals,output=proposed["output"],
                    final_coarse_sector=gamma,rejection_proposals=attempts,
                    peak_sector_count=peak,max_relative_completeness_error=drift,
                    max_proposals=int(max_proposals))
        raise RuntimeError("sparse coherent rejection cap exhausted; no fallback sample")
