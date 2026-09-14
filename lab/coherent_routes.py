"""Few coherent cell-reflection rotations on a supplied indexed orbit.

J_q=D_q R_cell is a Hermitian unitary, with coarse route alpha -> -alpha-q.
At insertion s apply exp(-i theta_s J_q/2) AFTER the repeated block W_s.
The finite history expansion is charged explicitly. No orbit, sector or
exponent-output table is constructed. This is NOT order/index discovery,
a generic physical few-qubit gate model, or a precision certificate.

Two samplers share the same capped, original-order amplitude contraction:
coherent-component rejection and the sparse-block gate-by-gate rule of
Bravyi, Gosset and Liu (arXiv:2112.08499, Algorithm 2). The latter needs prefix
amplitudes, not merely complete output probabilities. No state propagator
or singular-value/amplitude truncation is introduced here.
"""
from __future__ import annotations

import math
import numpy as np

from lab.fourier_sampling import _ints, unit_phase
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import sequential_path


def _route_initial_sector(final, selected, sectors):
    current = int(final)
    for q in reversed(tuple(selected.values())):
        current = (-current-q) % sectors
    return current


def _contract_route_component(width, b, sectors, sector, exponent, stop,
                              boundary, measured, output, selected, *,
                              shift, defect, phase, inverse_sqrt_two, tail_scale,
                              apply_step=None):
    """Shared finite-product ordering for float and verified scalar arithmetic.

    Callbacks only supply scalar/vector arithmetic, not a second propagation
    algorithm. This evaluates one routed component's unnormalized amplitude.
    Optional apply_step(operation, vector, tag) can attach verified arithmetic
    at complete linear-contraction boundaries; it must preserve the operation.
    """
    current = _route_initial_sector(sector, selected, sectors)
    if apply_step is None:
        apply_step = lambda operation, vector, tag: operation(vector)
    vector = [1]+[0]*(b-1)
    if stop > 0 or boundary != "arithmetic":
        vector = apply_step(lambda v: defect(0, v), vector, ("background", 0))
    if 0 in selected:
        current = (-current-selected[0]) % sectors
    for i in range(stop):
        def control_step(vector):
            if i >= width-measured:
                shifted = shift(1 << i, current, vector)
                z = phase(-int(output)*(1 << i), 1 << width)
                return [(v+z*w)/2 for v,w in zip(vector,shifted)]
            if int(exponent) & (1 << i):
                vector = shift(1 << i, current, vector)
            return [v*inverse_sqrt_two for v in vector]
        vector = apply_step(control_step, vector, ("control", i))
        s = i+1
        if s < stop or boundary != "arithmetic":
            vector = apply_step(lambda v: defect(s, v), vector, ("background", s))
        if s in selected:
            current = (-current-selected[s]) % sectors
    return apply_step(lambda v: [a*tail_scale for a in v], vector, ("tail", stop))


class CoherentReflectionInput:
    """Supplied indexed gates only, r<=int64 max, b<=64, t<=63, k<=64.

    defects: insertion -> supplied b-by-b unitary, as in PeriodicOrbitCircuit.
    reflections: insertion -> (integer q, finite real theta). There is at
    most one reflection rotation at an insertion; s=0 and s=t are allowed.
    Unitarity of the reflection is algebraic, not verified by an r-sized
    numerical matrix. The physical/indexed interpretation is a supplied promise.
    This input-only class has NO history enumerator or sampler. Larger inputs
    must be consumed by explicitly resource-bounded merged adapters.
    """
    def __init__(self, period, block_size, width, defects, reflections):
        if (not _ints(width) or not 0 <= width <= 63
                or not isinstance(reflections, dict) or len(reflections) > 64):
            raise ValueError("require bounded width and at most 64 supplied rotations")
        for s, value in reflections.items():
            if (not _ints(s) or not 0 <= s <= width
                    or not isinstance(value, (tuple, list)) or len(value) != 2
                    or not _ints(value[0])
                    or not isinstance(value[1], (int, float, np.integer, np.floating))
                    or not math.isfinite(float(value[1]))):
                raise ValueError("reflection requires a valid insertion, integer q and finite angle")
        self.background = PeriodicOrbitCircuit(period, block_size, width, defects)
        self.period, self.b, self.width = int(period), int(block_size), int(width)
        self.sectors = self.background.sectors
        self.reflections = {
            int(s): (int(q) % self.sectors, float(theta))
            for s, (q, theta) in sorted(reflections.items())}
        self._terms = tuple((s, q, complex(math.cos(theta/2)),
                             complex(-1j*math.sin(theta/2)))
                            for s, (q, theta) in self.reflections.items())
        self.coefficient_l1 = math.prod(abs(c)+abs(d) for _, _, c, d in self._terms)

    def _validate_label(self, sector):
        if not _ints(sector) or not 0 <= sector < self.sectors:
            raise ValueError("sector outside supplied indexed orbit")

    def _validate_prefix(self, sector, exponent, stop, boundary, measured, output):
        """Shared request contract for history and merged prefix oracles."""
        self._validate_label(sector)
        if (not _ints(exponent, stop, measured, output)
                or not 0 <= exponent < 1 << self.width
                or not 0 <= stop <= self.width
                or boundary not in ("arithmetic", "background", "reflection")
                or not 0 <= measured <= self.width or not 0 <= output < 1 << measured
                or measured and (stop != self.width or boundary != "reflection")):
            raise ValueError("invalid gate-prefix or partial-QFT amplitude request")


class CoherentReflectionCircuit(CoherentReflectionInput):
    """Legacy history methods retain their intentional k<=8 input cap."""

    def __init__(self, period, block_size, width, defects, reflections):
        if not isinstance(reflections, dict) or len(reflections) > 8:
            raise ValueError("history circuit requires at most eight supplied rotations")
        super().__init__(period, block_size, width, defects, reflections)

    def _histories(self, stop, boundary):
        # Also guard an accidental unbound-method call on structural input.
        if not isinstance(self, CoherentReflectionCircuit) or len(self._terms) > 8:
            raise ValueError("history enumeration requires a capped history circuit")
        terms = tuple(v for v in self._terms
                      if v[0] < stop or v[0] == stop and boundary == "reflection")
        for mask in range(1 << len(terms)):
            coefficient, selected = 1+0j, {}
            for i, (s, q, c, d) in enumerate(terms):
                chosen = bool(mask & (1 << i))
                coefficient *= d if chosen else c
                if chosen:
                    selected[s] = q
            # Exact zero only: no tolerance-selected history deletion.
            if coefficient != 0:
                yield coefficient, selected

    def _initial_sector(self, final, selected):
        return _route_initial_sector(final, selected, self.sectors)

    def _component_vector(self, sector, exponent, stop, boundary, measured, output, selected):
        result = _contract_route_component(
            self.width, self.b, self.sectors, sector, exponent, stop, boundary,
            measured, output, selected,
            shift=lambda power,alpha,v: self.background.shift(power,alpha) @ v,
            defect=lambda s,v: self.background.defects.get(s,self.background.identity) @ v,
            phase=unit_phase, inverse_sqrt_two=1/math.sqrt(2),
            tail_scale=math.exp2(-(self.width-stop)/2))
        return np.asarray(result,dtype=complex)

    def prefix_vector(self, sector, exponent, stop, *, boundary="reflection",
                      measured=0, output=0):
        """sqrt(M) times the work-p vector for a gate/QFT prefix amplitude.

        stop controls have acted. At that insertion boundary is arithmetic
        (before W), background (after W, before K), or reflection (after K).
        measured>0 means that many LOW inverse-QFT output bits have been
        selected; remaining LOW exponent bits are fixed by exponent. This
        mode requires the complete arithmetic/work circuit. Amplitudes are
        unnormalized with respect to the selected outputs; phases are retained.
        """
        self._validate_prefix(sector, exponent, stop, boundary, measured, output)
        result = np.zeros(self.b, complex)
        for coefficient, selected in self._histories(stop, boundary):
            result += coefficient*self._component_vector(
                sector, exponent, stop, boundary, measured, output, selected)
        return result

    def joint_probability(self, sector, output):
        """p(final coarse sector, y), not the output marginal p(y)."""
        vector = self.prefix_vector(sector, 0, self.width,
                                    measured=self.width, output=output)
        return float(np.vdot(vector, vector).real)/self.sectors

    @staticmethod
    def _draw(weights, rng):
        weights = np.asarray(weights, dtype=float)
        total = float(weights.sum())
        if np.any(weights < 0) or not np.all(np.isfinite(weights)) or total <= 0:
            raise ArithmeticError("nonpositive/nonfinite gate-block weight; no biased fallback")
        # Select only a positive bin even if the final CDF rounds below one.
        positive = np.flatnonzero(weights > 0)
        cutoff = float(rng.random())*total
        cumulative = 0.
        for index in positive:
            cumulative += float(weights[index])
            if cutoff < cumulative:
                return int(index)
        return int(positive[-1])

    def sample(self, rng):
        """Sparse-block gate-by-gate sample using charged prefix amplitudes.

        Arithmetic is a monomial basis update. A repeated W changes only p;
        a reflection rotation changes alpha only within its known two-cycle.
        The terminating inverse QFT uses one two-outcome block at a time.
        The current classical work/exponent labels are sampled basis values,
        NOT a replacement for quantum phases: all phases live in the oracle.
        """
        sector = int(rng.integers(self.sectors))
        exponent = int(rng.integers(1 << self.width))
        p, evaluations, updates = 0, 0, 0
        for stop in range(self.width+1):
            if stop and exponent & (1 << (stop-1)):
                p = (p+(1 << (stop-1))) % self.b
            if stop in self.background.defects:
                vector = self.prefix_vector(sector, exponent, stop, boundary="background")
                p = self._draw(np.abs(vector)**2, rng)
                evaluations += 1
                updates += 1
            if stop in self.reflections:
                target = (-sector-self.reflections[stop][0]) % self.sectors
                if target != sector:
                    candidates = (sector, target)
                    weights = [abs(self.prefix_vector(a, exponent, stop)[p])**2
                               for a in candidates]
                    sector = candidates[self._draw(weights, rng)]
                    evaluations += 2
                    updates += 1
        output = 0
        for measured in range(1, self.width+1):
            children = (output, output | (1 << (measured-1)))
            weights = [abs(self.prefix_vector(sector, exponent, self.width,
                                             measured=measured, output=y)[p])**2
                       for y in children]
            output = children[self._draw(weights, rng)]
            evaluations += 2
            updates += 1
        return dict(output=output, final_coarse_sector=sector, final_within_sector=p,
                    method="sparse-block gate-by-gate prefix amplitudes",
                    prefix_vector_evaluations=evaluations, resampled_blocks=updates,
                    rejection_proposals=0, **self.stats())

    def _component_sample(self, final, selected, rng):
        current = self._initial_sector(final, selected)
        initial = self.background.defects.get(0, self.background.identity)[:, 0]
        if 0 in selected:
            current = (-current-selected[0]) % self.sectors
        pairs = []
        for i in range(self.width):
            W = self.background.defects.get(i+1, self.background.identity)
            pairs.append((W, W @ self.background.shift(1 << i, current)))
            if i+1 in selected:
                current = (-current-selected[i+1]) % self.sectors
        return sequential_path(pairs, initial, rng=rng)

    def rejection_weights(self, sector, output):
        """Joint target/proposal and acceptance, streamed over ALL histories."""
        self._validate_label(sector)
        if not _ints(output) or not 0 <= output < 1 << self.width:
            raise ValueError("output outside exponent register")
        coherent = np.zeros(self.b, complex)
        incoherent = 0.
        for coefficient, selected in self._histories(self.width, "reflection"):
            vector = self._component_vector(sector, 0, self.width, "reflection",
                                            self.width, output, selected)
            coherent += coefficient*vector
            incoherent += abs(coefficient)*float(np.vdot(vector, vector).real)
        numerator = float(np.vdot(coherent, coherent).real)
        denominator = self.coefficient_l1*incoherent
        acceptance = numerator/denominator if denominator else 0.
        if (not math.isfinite(acceptance) or acceptance < 0 or acceptance > 1+1e-10
                or denominator == 0 and numerator != 0):
            raise ArithmeticError("coherent rejection envelope violated numerically")
        return dict(joint_probability=numerator/self.sectors,
                    proposal_joint_probability=incoherent/(self.sectors*self.coefficient_l1),
                    acceptance_probability=min(1., acceptance))

    def sample_rejection(self, rng, *, max_proposals=100_000):
        """Independent proposal mechanism; cap exhaustion raises, never falls back.

        The uncapped exact-math mean is B^2 with B=sum_h |c_h|. Each attempt
        resamples BOTH final sector and history; neither is fixed outside the
        rejection loop. Every acceptance test evaluates all nonzero histories.
        """
        if not _ints(max_proposals) or not 1 <= max_proposals <= 10_000_000:
            raise ValueError("invalid rejection attempt cap")
        for attempt in range(1, int(max_proposals)+1):
            sector = int(rng.integers(self.sectors))
            selected = {}
            for s, q, c, d in self._terms:
                if self._draw((abs(c), abs(d)), rng):
                    selected[s] = q
            path = self._component_sample(sector, selected, rng)
            weights = self.rejection_weights(sector, path["output"])
            if float(rng.random()) < weights["acceptance_probability"]:
                return dict(output=path["output"], final_coarse_sector=sector,
                            method="coherent-component rejection", rejection_proposals=attempt,
                            **weights, **self.stats())
        raise RuntimeError("coherent rejection attempt cap exhausted; no fallback sample")

    def stats(self):
        return dict(period=self.period, block_size=self.b, width=self.width,
                    coherent_rotation_count=len(self._terms),
                    history_count_upper_bound=1 << len(self._terms),
                    coefficient_l1=self.coefficient_l1,
                    mathematical_mean_rejection_proposals=self.coefficient_l1**2,
                    orbit_table_entries=0, sector_table_entries=0, output_table_entries=0,
                    owned_matrix_payload_bytes=self.background.stats()["owned_block_payload_bytes"],
                    history_storage="streamed; O(k) routing integers and O(b) amplitude scalars",
                    setup_and_precision_costs_included=False)
