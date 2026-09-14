"""Coherently merged backward blocks as the existing gate-by-gate oracle.

This is an opt-in complex128 diagnostic for C59's supplied indexed family.
It batches all fine-boundary columns through C71's sparse adjoint helpers,
preserving intermediate-prefix phases and endpoint inclusion. The original
C59 sampling loop is reused, with no histories, rejection or dense orbit.
There is no finite-TV certificate or approximation/truncation parameter.
"""
from __future__ import annotations

import math
import numpy as np

from lab.coherent_routes import CoherentReflectionCircuit, CoherentReflectionInput
from lab.coherent_reverse import SparseCoherentReverse, route_support_bound
from lab.fourier_sampling import unit_phase
from lab.semiclassical import _qft_branch_operators


class MergedCoherentPrefixes:
    """Batched b-by-b reverse blocks with at most S reached coarse labels.

    Legacy input caps are unchanged. Explicit CoherentReflectionInput allows
    k<=64 subject to a 16-MiB arithmetic-payload preflight. Fixed d labels
    give polynomial support; two labels retain S<=1+2k. Batching costs
    O(S*b^2) working complex entries, not the vector sampler's O(S*b).
    Counters are diagnostic mutable state, not a concurrent-worker API.
    """
    _draw = staticmethod(CoherentReflectionCircuit._draw)

    def __init__(self, circuit, *, support_mode="alphabet"):
        # One shared envelope setup and one owned gate snapshot. No prefix
        # query rebuilds E/O sets; the full-word bound covers all subwords.
        self.reverse = SparseCoherentReverse(circuit,support_mode=support_mode)
        if not isinstance(circuit, CoherentReflectionCircuit) or support_mode=="word":
            bound = self.reverse.support_bound
            slots = ((24+6*bound)*circuit.b*circuit.b
                     +(len(circuit.background.defects)+1)*circuit.b*circuit.b)
            if bound > 65_536 or 16*slots > (16 << 20):
                raise MemoryError("structural prefix cover exceeds diagnostic payload budget")
        self.circuit = self.reverse.circuit  # owns the single input snapshot
        for name in ("period", "b", "width", "sectors", "background", "reflections"):
            setattr(self, name, getattr(self.circuit, name))
        self.support_bound = self.reverse.support_bound
        self.reset_stats()

    def reset_stats(self):
        self._counts = dict(merged_prefix_queries=0, prefix_control_steps=0,
                           background_block_products=0, control_block_products=0,
                           fixed_control_block_scales=0, qft_matrix_pair_constructions=0,
                           reflection_block_contributions=0, peak_sector_count=0)
        self.last_prefix_stats = None

    def stats(self):
        return dict(**self._counts, period=self.period, block_size=self.b, width=self.width,
                    **self.reverse.support_plan,
                    coherent_rotation_count=len(self.reflections),
                    history_count_upper_bound=1 << len(self.reflections),
                    history_enumerations=0, orbit_table_entries=0,
                    sector_table_entries=0, output_table_entries=0,
                    working_complex_coordinate_slots_upper_bound=
                        (24+6*self.support_bound)*self.b*self.b,
                    supplied_gate_complex_entries=(len(self.background.defects)+1)*self.b*self.b,
                    precision_certificate=False, arithmetic="complex128 diagnostic",
                    order_and_index_discovery_included=False,
                    history_storage="none; reached-sector b-by-b adjoint blocks")

    def prefix_vector(self, sector, exponent, stop, *, boundary="reflection",
                      measured=0, output=0):
        self.circuit._validate_prefix(sector, exponent, stop, boundary, measured, output)
        included = tuple(q for s, (q, theta) in self.reflections.items()
                         if s < stop or s == stop and boundary == "reflection")
        bound = min(route_support_bound(self.sectors, included),self.support_bound)
        blocks = {int(sector): self.background.identity.copy()}
        # The shared helpers act columnwise on these matrices. Their local
        # matvec/contribution counters count BLOCK operations here, not b-vectors.
        local = dict(peak_sector_count=1, reflection_vector_contributions=0, work_matvecs=0)
        control_products = fixed_scales = qft_pairs = 0
        for s in range(int(stop), 0, -1):
            if s < stop or boundary == "reflection":
                blocks = self.reverse._reflection_adjoint(blocks, s, local)
            if s < stop or boundary != "arithmetic":
                blocks = self.reverse._background_adjoint(blocks, s, local)
            i = s-1
            children = {}
            if i >= self.width-measured:
                phase = unit_phase(-int(output)*(1 << i), 1 << self.width)
                for alpha, block in blocks.items():
                    pair = _qft_branch_operators(self.background.identity,
                        self.background.shift(1 << i, alpha), phase)
                    children[alpha] = pair[0].conj().T @ block
                qft_pairs += len(blocks)
                control_products += len(blocks)
            else:
                shifted = bool(int(exponent) & (1 << i))
                for alpha, block in blocks.items():
                    if shifted:
                        block = self.background.shift(1 << i, alpha).conj().T @ block
                    children[alpha] = block / math.sqrt(2)
                control_products += int(shifted)*len(blocks)
                fixed_scales += len(blocks)
            blocks = children
        if stop > 0 or boundary == "reflection":
            blocks = self.reverse._reflection_adjoint(blocks, 0, local)
        if stop > 0 or boundary != "arithmetic":
            blocks = self.reverse._background_adjoint(blocks, 0, local)
        if local["peak_sector_count"] > bound:
            raise ArithmeticError("prefix-specific algebraic support bound violated")
        # Column j is F_prefix^dag |sector,j>. Its coherent initial overlap
        # is conjugated to return the FORWARD amplitude in C59's convention.
        result = np.array([complex(math.fsum(v[0,j].real for v in blocks.values()),
                                   -math.fsum(v[0,j].imag for v in blocks.values()))
                           for j in range(self.b)])
        result *= math.exp2(-(self.width-stop)/2)
        if not np.all(np.isfinite(result)):
            raise ArithmeticError("nonfinite merged prefix amplitude; no fallback")
        cost = dict(merged_prefix_queries=1, prefix_control_steps=int(stop),
                    background_block_products=local["work_matvecs"],
                    control_block_products=control_products,
                    fixed_control_block_scales=fixed_scales,
                    qft_matrix_pair_constructions=qft_pairs,
                    reflection_block_contributions=local["reflection_vector_contributions"],
                    peak_sector_count=local["peak_sector_count"])
        self.last_prefix_stats = dict(**cost, prefix_support_bound=bound,
                                      included_reflections=len(included),
                                      terminal_sector_count=len(blocks))
        for key, value in cost.items():
            self._counts[key] = (max(self._counts[key], value) if key == "peak_sector_count"
                                 else self._counts[key]+value)
        return result

    def joint_probability(self, sector, output):
        return CoherentReflectionCircuit.joint_probability(self, sector, output)

    def sample(self, rng):
        self.reset_stats()
        result = CoherentReflectionCircuit.sample(self, rng)
        result["method"] = "existing gate-by-gate sampler with coherently merged prefix blocks"
        return result
