"""Opt-in C78/C79/C81 work-first sampling for restricted mixer/phase schedules.

Supplied known/indexed period r=b*M; initial W0|0> in the first b-cell.
The literal pre-final column is U^(L*h) D_g W1^rep U^l W0|0>, e=l+L*h.
No intervening late work gates are allowed. Common terminal work unitaries
can be stripped under the trace. This is not a generic circuit propagator.
EarlierPhaseProgressions additionally permits one pointwise phase before
the late mixer, using complementary disjoint progression covers.
NestedPhaseProgressions permits several diagonal insertions without mixers.

No orbit, exponent, output or expanded-progression arrays are constructed.
The phase callable must be deterministic and pointwise unit modulus; only
queried values are validated. Its internal costs/storage are the caller's
responsibility. Arithmetic identities are exact; float64 arithmetic and
finite RNG precision are NOT certified. Exhaustion raises: ignoring failed
calls and retaining only successes can bias the latent work distribution.
"""
from __future__ import annotations

import math
import numpy as np

from lab.fourier_sampling import geometric_sum, progression_sample, unit_phase


def _integer(value):
    return isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))


def _counters():
    return dict(column_local_terms=0, row_local_terms=0, phase_queries=0,
                fourier_component_terms=0, progression_proposals=0,
                marginal_queries=0, work_draws=0, component_draws=0,
                acceptance_draws=0, component_weight_terms=0,
                component_sqrt_terms=0, weighted_ratio_divisions=0)


def _proposal_mode(value):
    if not isinstance(value, str) or value not in ("mass", "root_mass"):
        raise ValueError("proposal must be mass or root_mass")
    return value


def _checked_phase(oracle, index, counters):
    counters["phase_queries"] += 1
    value = complex(oracle(int(index)))
    if (not math.isfinite(value.real) or not math.isfinite(value.imag)
            or abs(abs(value)-1.) > 1e-12):
        raise ValueError("queried pointwise phase must be finite and unit modulus")
    return value


def _weight(value):
    weight = float(abs(value)**2)
    if not math.isfinite(weight) or (value != 0 and weight == 0):
        raise ArithmeticError("nonfinite or underflowed nonzero amplitude weight")
    return weight


def _pick(weights, rng):
    """Finite categorical draw; ignore exact zeros, never tiny positive weights."""
    total = math.fsum(weights)
    if not math.isfinite(total) or total <= 0 or any(w < 0 for w in weights):
        raise ArithmeticError("invalid categorical weights")
    threshold = float(rng.random()) * total
    cumulative, last = 0., None
    for index, weight in enumerate(weights):
        if weight == 0:
            continue
        last = index
        cumulative += weight
        if threshold < cumulative:
            return index
    # Last-bit rounding in the cumulative sum, not a rejection fallback.
    return last


class LateWorkProgressions:
    """Bounded small-block, few-late-history FLOAT sampler; supplied order/index.

    width <= 63, block_size <= 64. Caps are checked before matrix copies and
    structural loops. phase(j) receives an orbit INDEX, not a physical residue.
    Matrices use output rows/input columns and are copied read-only.
    """

    def __init__(self, period, block_size, width, split, initial, middle, phase,
                 *, max_local_terms=1_000_000, max_components=4096,
                 max_payload_bytes=16 << 20):
        if (not all(_integer(x) for x in (period, block_size, width, split))
                or not 1 <= period < 1 << 63 or not 1 <= block_size <= 64
                or period % block_size or not 0 <= split <= width <= 63):
            raise ValueError("require indexed r=b*M, 1<=b<=64, 0<=split<=width<=63")
        if (not all(_integer(x) and x > 0 for x in
                    (max_local_terms, max_components, max_payload_bytes))
                or max_local_terms > 1_000_000 or max_components > 4096
                or max_payload_bytes > 16 << 20):
            raise ValueError("invalid local-term/component/numeric-payload caps")
        if (not callable(phase) or not all(isinstance(w, np.ndarray)
                and w.shape == (block_size, block_size) for w in (initial, middle))):
            raise ValueError("supply two b-by-b ndarrays and a pointwise phase callable")
        self.period, self.b, self.width, self.split = map(int, (period, block_size, width, split))
        self.Q, self.L, self.H = 1 << self.width, 1 << self.split, 1 << (self.width-self.split)
        (self.component_bound, self.local_term_bound,
         self.numeric_payload_bound) = self._structural_bounds()
        if (self.component_bound > max_components or self.local_term_bound > max_local_terms
                or self.numeric_payload_bound > max_payload_bytes):
            raise MemoryError("work-first structural preflight exceeds requested cap")
        matrices = []
        for value in (initial, middle):
            block = np.array(value, dtype=np.complex128, copy=True)
            if (not np.all(np.isfinite(block)) or not np.allclose(
                    block.conj().T @ block, np.eye(self.b), rtol=0, atol=1e-12)):
                raise ValueError("work blocks must be finite unitary matrices")
            block.flags.writeable = False
            matrices.append(block)
        self.initial, self.middle = matrices
        self.phase = phase

    def _structural_bounds(self):
        component_bound = self.H * min(self.period, 2*self.b-1)
        local_bound = (self.H+1)*self.b*self.b
        # Numeric payload, not Python object headers or arbitrary oracle storage.
        payload = 16*(16*self.b*self.b + 12*component_bound + 128)
        return component_bound, local_bound, payload

    def _new_counters(self):
        return _counters()

    def _phase(self, index, counters):
        return _checked_phase(self.phase, index, counters)

    def _column(self, exponent, counters):
        low, high = exponent % self.L, exponent // self.L
        before = {}
        for u in range(self.b):
            cell, q = divmod((u+low) % self.period, self.b)
            for p in range(self.b):
                counters["column_local_terms"] += 1
                amplitude = complex(self.initial[u, 0]*self.middle[p, q])
                index = cell*self.b+p
                before[index] = before.get(index, 0j)+amplitude
        result = {}
        for index, amplitude in before.items():
            if amplitude != 0:
                result[(index+self.L*high) % self.period] = amplitude*self._phase(index, counters)
        norm = math.fsum(_weight(a) for a in result.values())
        if not math.isfinite(norm) or abs(norm-1.) > 1e-10:
            raise ArithmeticError("coherent column is not normalized")
        return result

    def column(self, exponent):
        """Bounded coherent column diagnostic; no orbit table."""
        if not _integer(exponent) or not 0 <= exponent < self.Q:
            raise ValueError("exponent outside register")
        counters = self._new_counters()
        return dict(amplitudes=self._column(int(exponent), counters), counters=counters)

    def _row(self, work, counters):
        components = []
        for high in range(self.H):
            v = (work-self.L*high) % self.period
            cell, p = divmod(v, self.b)
            residues = {}
            for u in range(self.b):
                for q in range(self.b):
                    counters["row_local_terms"] += 1
                    rho = (cell*self.b+q-u) % self.period
                    amplitude = complex(self.middle[p, q]*self.initial[u, 0])
                    residues[rho] = residues.get(rho, 0j)+amplitude
            phase = self._phase(v, counters)
            for rho, amplitude in sorted(residues.items()):
                count = max(0, 1+(self.L-1-rho)//self.period)
                if count and amplitude != 0:
                    components.append((self.L*high+rho, count, phase*amplitude))
        norm = math.fsum(count*_weight(a) for _, count, a in components)
        return dict(work=work, components=tuple(components), stride=self.period,
                    norm=norm, counters=counters)

    def row(self, work):
        """Merged disjoint input progressions; empty rows have norm zero."""
        if not _integer(work) or not 0 <= work < self.period:
            raise ValueError("work index outside supplied period")
        return self._row(int(work), self._new_counters())

    def _prepare_proposal(self, row, mode, counters):
        """Raw weights avoid normalizing tiny component masses before sqrt.

        Internal rows already omit exact-zero components. A nonzero mass
        that underflows raises through _weight; no diagnostic cutoff is used.
        The original numeric preflight has room for this one component list.
        """
        mode = _proposal_mode(mode)
        weights = []
        for _, count, amplitude in row["components"]:
            counters["component_weight_terms"] += 1
            omega = count*_weight(amplitude)
            if not math.isfinite(omega) or omega <= 0:
                raise ArithmeticError("internal component must have positive finite mass")
            if mode == "root_mass":
                counters["component_sqrt_terms"] += 1
                omega = math.sqrt(omega)
            weights.append(omega)
        scale = math.fsum(weights) if mode == "root_mass" else row["norm"]
        expected = ((scale*scale/row["norm"] if row["norm"] else 0.)
                    if mode == "root_mass" else float(len(weights)))
        if weights and (not math.isfinite(scale) or not math.isfinite(expected)
                        or not 1-1e-12 <= expected <= len(weights)*(1+1e-12)):
            raise ArithmeticError("invalid component proposal normalization/envelope")
        return dict(mode=mode, weights=weights, scale=scale, expected_attempts=expected)

    def _fourier(self, row, output, counters, *, prepared=None):
        values = []
        for start, count, amplitude in row["components"]:
            counters["fourier_component_terms"] += 1
            values.append(amplitude*unit_phase(-output*start, self.Q)
                          * geometric_sum(count, -output*row["stride"], self.Q))
        coherent = _weight(sum(values, 0j))
        m, norm = len(values), row["norm"]
        if norm == 0:
            return dict(conditional_probability=None, proposal_probability=0., acceptance=0.,
                        expected_attempts=0.)
        if prepared is not None and prepared["mode"] == "root_mass":
            quotients = []
            for value, root_weight in zip(values, prepared["weights"], strict=True):
                counters["weighted_ratio_divisions"] += 1
                square = _weight(value)
                quotient = square/root_weight
                if not math.isfinite(quotient) or (square != 0 and quotient == 0):
                    raise ArithmeticError("nonfinite or underflowed weighted Fourier term")
                quotients.append(quotient)
            diagonal = math.fsum(quotients)
            denominator = prepared["scale"]*diagonal
            proposal_probability = diagonal/self.Q/prepared["scale"]
            if diagonal != 0 and proposal_probability == 0:
                raise ArithmeticError("underflowed nonzero weighted proposal probability")
            expected = prepared["expected_attempts"]
        else:
            diagonal = math.fsum(_weight(a) for a in values)
            denominator = m*diagonal
            proposal_probability = diagonal/self.Q/norm
            expected = float(m)
        acceptance = coherent/denominator if denominator > 0 else 0.
        if (not math.isfinite(acceptance) or not 0 <= acceptance <= 1+1e-12
                or not math.isfinite(diagonal) or not math.isfinite(denominator)
                or not math.isfinite(proposal_probability)
                or (denominator == 0 and coherent != 0)):
            raise ArithmeticError("invalid coherent rejection envelope")
        return dict(conditional_probability=coherent/self.Q/norm,
                    proposal_probability=proposal_probability,
                    acceptance=min(1., acceptance), expected_attempts=expected)

    def forced_joint(self, work, output, *, proposal="mass"):
        """Evaluate one joint work/output probability; summing work is NOT free."""
        _proposal_mode(proposal)
        if not _integer(output) or not 0 <= output < self.Q:
            raise ValueError("Fourier output outside register")
        row = self.row(work)
        prepared = (self._prepare_proposal(row, proposal, row["counters"])
                    if proposal == "root_mass" else None)
        result = self._fourier(row, int(output), row["counters"], prepared=prepared)
        weight = row["norm"]/self.Q
        result.update(work=int(work), output=int(output), work_probability=weight,
                      joint_probability=weight*(result["conditional_probability"] or 0.),
                      component_count=len(row["components"]), proposal_mode=proposal,
                      counters=row["counters"])
        return result

    def sample(self, rng, *, max_attempts=10_000, proposal="mass"):
        """Draw work once, discard seed exponent, retry ONLY conditional output.

        Default mass weights have mathematical mean m. Opt-in root_mass uses
        the C59 weighted envelope, with mean (sum sqrt(component mass))^2/Z.
        Neither finite attempts nor floating arithmetic is an accuracy certificate.
        """
        _proposal_mode(proposal)
        if not _integer(max_attempts) or not 1 <= max_attempts <= 100_000:
            raise ValueError("max_attempts must be an integer in [1,100000]")
        counters = self._new_counters()
        column = self._column(int(rng.integers(self.Q)), counters)
        labels = tuple(column)
        counters["work_draws"] += 1
        work = labels[_pick([_weight(column[j]) for j in labels], rng)]
        # No dependence on the auxiliary exponent survives this point.
        del column, labels
        row = self._row(work, counters)
        return self._sample_row(row, rng, counters, max_attempts=max_attempts, proposal=proposal)

    def _sample_row(self, row, rng, counters, *, max_attempts, proposal="mass"):
        """Shared conditional loop; private hook for isolated row-law tests."""
        if row["norm"] <= 0 or not row["components"]:
            raise ArithmeticError("sampled a zero-weight work row")
        prepared = self._prepare_proposal(row, proposal, counters)
        weights = prepared["weights"]
        for attempt in range(1, int(max_attempts)+1):
            counters["component_draws"] += 1
            start, count, _ = row["components"][_pick(weights, rng)]
            counters["progression_proposals"] += 1
            candidate = progression_sample(self.width, start, row["stride"], count, 0, 1, rng)
            counters["marginal_queries"] += candidate["marginal_queries"]
            law = self._fourier(row, candidate["output"], counters, prepared=prepared)
            counters["acceptance_draws"] += 1
            if rng.random() < law["acceptance"]:
                return dict(work=row["work"], output=candidate["output"], attempts=attempt,
                            component_count=len(row["components"]), counters=counters,
                            proposal_mode=prepared["mode"],
                            expected_attempts=prepared["expected_attempts"])
        raise RuntimeError("conditional Fourier rejection cap exhausted; do not discard failed calls")

    def stats(self):
        """Structural bounds, not measured RSS, native FLOPs or phase bit cost."""
        return dict(component_bound=self.component_bound, local_term_bound=self.local_term_bound,
                    matrix_payload_bytes=self.initial.nbytes+self.middle.nbytes,
                    numeric_payload_bound_bytes=self.numeric_payload_bound,
                    python_container_entry_bound=8*(self.component_bound+self.b*self.b)+128,
                    unitary_check_scalar_product_terms=2*self.b**3,
                    phase_queries_per_draw_bound=self.b*self.b+self.H,
                    proposal_modes=("mass", "root_mass"), default_proposal="mass",
                    proposal_preparation_terms_bound=self.component_bound,
                    root_mass_sqrt_terms_bound=self.component_bound,
                    root_mass_divisions_per_attempt_bound=self.component_bound,
                    orbit_table_entries=0, exponent_table_entries=0,
                    output_table_entries=0, expanded_progression_entries=0,
                    order_and_index_discovery_included=False,
                    precision="complex128/float64; no global finite-bit certificate",
                    phase_oracle_storage_and_bit_cost_included=False)


class EarlierPhaseProgressions(LateWorkProgressions):
    """C79 specialization: one extra pointwise phase after early_split controls.

    Literal columns are U^(L*h) D_g W1^rep U^(l-a) D_f U^a W0|0>,
    where a=l mod 2^early_split. Only this schedule is accepted. The two
    covers share the original C78 sampler; their row stride may differ.
    Auto chooses a conservative local-pair-visit bound, not measured time
    or the smallest realized nonzero component count. No orbit/Q-sized
    structure is constructed. Caps include both the column and chosen row.
    """

    def __init__(self, period, block_size, width, split, initial, middle, phase,
                 *, early_split, early_phase, cover="auto",
                 max_local_terms=1_000_000, max_components=4096,
                 max_payload_bytes=16 << 20):
        if (not _integer(early_split) or not _integer(split)
                or not 0 <= early_split <= split):
            raise ValueError("require 0<=early_split<=split")
        if not callable(early_phase):
            raise ValueError("supply a pointwise early phase callable")
        if not isinstance(cover, str) or cover not in ("auto", "cycle", "dual"):
            raise ValueError("cover must be auto, cycle or dual")
        self.early_split = int(early_split)
        self.early_phase = early_phase
        self.requested_cover = cover
        super().__init__(period, block_size, width, split, initial, middle, phase,
                         max_local_terms=max_local_terms, max_components=max_components,
                         max_payload_bytes=max_payload_bytes)

    def _structural_bounds(self):
        self.A = 1 << self.early_split
        self.P = self.A // math.gcd(self.period, self.A)
        self.K = self.L // self.A
        # At the final pre-mixer endpoint, f is fixed at the mixer source.
        self.effective_P = 1 if self.early_split == self.split else self.P
        self.cycle_classes_bound = min(self.effective_P,
                                      (self.L+self.period-1)//self.period)
        base_pairs = self.H*self.b*self.b
        self.cycle_pair_bound = base_pairs*(1+self.cycle_classes_bound)
        self.dual_pair_bound = base_pairs*self.K
        self.cover = ("cycle" if self.cycle_pair_bound <= self.dual_pair_bound
                      else "dual") if self.requested_cover == "auto" else self.requested_cover
        residues = min(self.period, 2*self.b-1)
        occupied_bound = self.H*residues*((self.L+self.period-1)//self.period)
        if self.cover == "cycle":
            components = self.H*min(self.L, residues*self.cycle_classes_bound)
            row_pairs = self.cycle_pair_bound
            self.row_early_query_bound = base_pairs*self.cycle_classes_bound
            self.row_stride = self.period*self.effective_P
        else:
            components = self.H*self.K*min(self.A, residues)
            row_pairs = self.dual_pair_bound
            self.row_early_query_bound = self.H*self.K*self.b
            self.row_stride = self.period
        components = min(components, occupied_bound, self.Q)
        # Includes copied matrices/unitary-check temps, the local pair lists,
        # a column and its weights, retained row tuples and Fourier temporaries.
        payload = 16*(32*self.b*self.b + 16*components + 128)
        return components, row_pairs+self.b*self.b, payload

    def _new_counters(self):
        counters = super()._new_counters()
        counters.update(early_phase_queries=0, late_phase_queries=0)
        return counters

    def _phase(self, index, counters):
        counters["late_phase_queries"] += 1
        return super()._phase(index, counters)

    def _early(self, index, counters):
        counters["early_phase_queries"] += 1
        return _checked_phase(self.early_phase, int(index) % self.period, counters)

    def _column(self, exponent, counters):
        low, high = exponent % self.L, exponent // self.L
        prefix = low % self.A
        before = {}
        for u in range(self.b):
            initial = complex(self.initial[u, 0])*self._early(u+prefix, counters)
            cell, q = divmod((u+low) % self.period, self.b)
            for p in range(self.b):
                counters["column_local_terms"] += 1
                index = cell*self.b+p
                amplitude = initial*self.middle[p, q]
                before[index] = before.get(index, 0j)+amplitude
        result = {}
        for index, amplitude in before.items():
            if amplitude != 0:
                result[(index+self.L*high) % self.period] = amplitude*self._phase(index, counters)
        norm = math.fsum(_weight(a) for a in result.values())
        if not math.isfinite(norm) or abs(norm-1.) > 1e-10:
            raise ArithmeticError("coherent early-phase column is not normalized")
        return result

    def _append(self, components, start, count, amplitude):
        if count and amplitude != 0:
            if len(components) >= self.component_bound:
                raise MemoryError("earlier-phase component bound exceeded before append")
            if not 0 <= start <= start+self.row_stride*(count-1) < self.Q:
                raise ArithmeticError("earlier-phase component leaves exponent register")
            components.append((start, count, complex(amplitude)))

    def _row(self, work, counters):
        components = []
        for high in range(self.H):
            w = (work-self.L*high) % self.period
            cell, p = divmod(w, self.b)
            c = cell*self.b
            late = self._phase(w, counters)
            if self.cover == "cycle":
                pairs = {}
                for u in range(self.b):
                    for q in range(self.b):
                        counters["row_local_terms"] += 1
                        rho = (c+q-u) % self.period
                        pairs.setdefault(rho, []).append(
                            (u, complex(self.middle[p, q]*self.initial[u, 0])))
                for rho, local_pairs in sorted(pairs.items()):
                    total_count = max(0, 1+(self.L-1-rho)//self.period)
                    for z in range(min(self.effective_P, total_count)):
                        start_low = rho+self.period*z
                        count = 1+(total_count-1-z)//self.effective_P
                        gamma = 0j
                        for u, amplitude in local_pairs:
                            counters["row_local_terms"] += 1
                            gamma += amplitude*self._early(u+(start_low % self.A), counters)
                        self._append(components, self.L*high+start_low, count, late*gamma)
            else:
                for k in range(self.K):
                    shift = self.A*k
                    residues = {}
                    for q in range(self.b):
                        early = self._early(c+q-shift, counters)
                        for u in range(self.b):
                            counters["row_local_terms"] += 1
                            rho = (c+q-shift-u) % self.period
                            amplitude = complex(self.middle[p, q]*self.initial[u, 0])*early
                            residues[rho] = residues.get(rho, 0j)+amplitude
                    for rho, amplitude in sorted(residues.items()):
                        count = max(0, 1+(self.A-1-rho)//self.period)
                        self._append(components, self.L*high+shift+rho, count, late*amplitude)
        norm = math.fsum(count*_weight(a) for _, count, a in components)
        if not math.isfinite(norm):
            raise ArithmeticError("nonfinite earlier-phase row norm")
        return dict(work=work, components=tuple(components), stride=self.row_stride,
                    norm=norm, cover=self.cover, counters=counters)

    def stats(self):
        result = super().stats()
        result.update(early_split=self.early_split, requested_cover=self.requested_cover,
                      cover=self.cover, generic_cycle_period=self.P,
                      effective_cycle_period=self.effective_P,
                      cycle_classes_bound=self.cycle_classes_bound,
                      remaining_histories=self.K, row_stride=self.row_stride,
                      cycle_row_pair_visit_bound=self.cycle_pair_bound,
                      dual_row_pair_visit_bound=self.dual_pair_bound,
                      phase_queries_per_draw_bound=(self.b+self.b*self.b+self.H
                                                    +self.row_early_query_bound),
                      early_phase_queries_per_draw_bound=self.b+self.row_early_query_bound,
                      late_phase_queries_per_draw_bound=self.b*self.b+self.H,
                      python_container_entry_bound=12*(self.component_bound+self.b*self.b)+128,
                      row_local_terms_semantics="local pair visits; cycle counts grouping and phase evaluation",
                      cover_selection="smaller safe pair-visit bound before construction; ties choose cycle")
        return result


class NestedPhaseProgressions(LateWorkProgressions):
    """C81 nested-phase rows in the unchanged work-first conditional loop.

    early_phases is a list/tuple of at most 64 (position, callable) pairs;
    each phase follows that many ascending controls, before the late mixer.
    Duplicate/initial/end positions are permitted. Pointwise phases must be
    deterministic; only queried values are checked. cut='auto' minimizes
    grouping plus coefficient pair visits, with weights (1,1), using the
    exact finite-truncated T count. Ties prefer the larger cut. This metric
    excludes phase costs, cancellation, proposal costs and rejection counts.

    cache_right caches only the right-phase product for each local q in one
    (high,k) slice. It does not change the selection metric. last_counters
    retains work from the most recent operation, including exhausted calls;
    callers must aggregate each operation before starting the next one.
    """

    def __init__(self, period, block_size, width, split, initial, middle, phase,
                 *, early_phases=(), cut="auto", cache_right=False,
                 max_local_terms=1_000_000, max_components=4096,
                 max_payload_bytes=16 << 20):
        if (not _integer(split) or not isinstance(early_phases, (list, tuple))
                or len(early_phases) > 64):
            raise ValueError("supply at most 64 insertion pairs and an integer split")
        if any(not isinstance(item, (list, tuple)) or len(item) != 2
               or not _integer(item[0]) or not 0 <= item[0] <= split
               or not callable(item[1]) for item in early_phases):
            raise ValueError("each insertion requires 0<=position<=split and a callable")
        if not (isinstance(cut, str) and cut == "auto") and not (
                _integer(cut) and 0 <= cut <= split):
            raise ValueError("cut must be auto or an integer in [0,split]")
        if not isinstance(cache_right, bool):
            raise ValueError("cache_right must be a bool")
        self.early_phases = tuple(sorted(((int(v), f) for v, f in early_phases),
                                         key=lambda item: item[0]))
        self.requested_cut, self.cache_right = cut, cache_right
        self.last_counters = None
        super().__init__(period, block_size, width, split, initial, middle, phase,
                         max_local_terms=max_local_terms, max_components=max_components,
                         max_payload_bytes=max_payload_bytes)

    def _cut_parameters(self, cut):
        left = max((v for v, _ in self.early_phases if v < cut), default=0)
        period = (1 << left)//math.gcd(self.period, 1 << left)
        A = 1 << cut
        return A, self.L//A, period

    def _structural_bounds(self):
        self.candidate_cuts = tuple(sorted({self.split} | {v for v, _ in self.early_phases}))
        cuts = self.candidate_cuts if self.requested_cut == "auto" else (int(self.requested_cut),)
        self.cut_parameters = {cut: self._cut_parameters(cut) for cut in cuts}
        self.geometry_setup_phase_terms = len(cuts)*len(self.early_phases)
        base = self.H*self.b*self.b
        ceil_length = (self.L+self.period-1)//self.period
        residues = min(self.period, 2*self.b-1)
        occupied = self.H*residues*ceil_length
        self.selector_pair_bound = base if self.requested_cut == "auto" else 0
        self.selector_setup_bound = (base+self.H*residues+len(self.candidate_cuts)
                                     if self.requested_cut == "auto" else 0)
        if self.requested_cut == "auto":
            # Actual selected g+C <= min_w (g_w+C_bound_w). Also g>=base,
            # G<=C and every component occupies a distinct compatible label.
            bounds = []
            for cut in self.candidate_cuts:
                _, K, P = self.cut_parameters[cut]
                bounds.append(base*K + min(base*K*P, base*ceil_length))
            self.row_pair_bound = min(bounds)
            self.coefficient_pair_bound = self.row_pair_bound-base
            self.group_pair_bound = self.row_pair_bound
            components = min(self.coefficient_pair_bound, occupied, self.Q)
        else:
            A, K, P = self.cut_parameters[int(self.requested_cut)]
            self.group_pair_bound = base*K
            self.coefficient_pair_bound = min(base*K*P, base*ceil_length)
            self.row_pair_bound = self.group_pair_bound+self.coefficient_pair_bound
            F = min(P, (A+self.period-1)//self.period)
            components = min(self.H*K*min(A, residues*F), occupied, self.Q)
        d = len(self.early_phases)
        # O(b^2) pair groups, O(b) right cache, O(d) schedule/candidates;
        # row/tuple conversion, proposal and Fourier lists coexist. Reserve
        # numeric slots conservatively; Python headers/oracle storage excluded.
        payload = 16*(40*self.b*self.b + 20*components + 16*d + 256)
        self.row_phase_query_bound = self.H+d*self.coefficient_pair_bound
        self.row_phase_product_bound = (d+1)*self.coefficient_pair_bound+components
        if self.cache_right:
            cache_pairs = self.group_pair_bound//self.b
            self.row_phase_query_bound += d*cache_pairs
            self.row_phase_product_bound += d*cache_pairs
        local = self.selector_setup_bound+self.row_pair_bound+self.b*self.b+d
        return components, local, payload

    def _new_counters(self):
        counters = super()._new_counters()
        counters.update(selector_pair_terms=0, selector_residue_terms=0,
                        selector_candidate_terms=0, group_pair_visits=0,
                        coefficient_pair_visits=0, geometric_components=0,
                        early_phase_queries=0, late_phase_queries=0,
                        column_phase_products=0, row_phase_products=0,
                        right_cache_entries=0, norm_terms=0, phase_partition_terms=0)
        self.last_counters = counters
        return counters

    def _phase(self, index, counters):
        counters["late_phase_queries"] += 1
        return super()._phase(index, counters)

    def _early(self, oracle, index, counters):
        counters["early_phase_queries"] += 1
        return _checked_phase(oracle, int(index) % self.period, counters)

    def _column(self, exponent, counters):
        low, high = exponent % self.L, exponent // self.L
        before = {}
        for u in range(self.b):
            initial = complex(self.initial[u, 0])
            for v, oracle in self.early_phases:
                initial *= self._early(oracle, u+low % (1 << v), counters)
                counters["column_phase_products"] += 1
            cell, q = divmod((u+low) % self.period, self.b)
            for p in range(self.b):
                counters["column_local_terms"] += 1
                index = cell*self.b+p
                before[index] = before.get(index, 0j)+initial*self.middle[p, q]
        result = {}
        for index, amplitude in before.items():
            if amplitude != 0:
                result[(index+self.L*high) % self.period] = amplitude*self._phase(index, counters)
                counters["column_phase_products"] += 1
        norm = math.fsum(_weight(a) for a in result.values())
        if not math.isfinite(norm) or abs(norm-1.) > 1e-10:
            raise ArithmeticError("coherent nested-phase column is not normalized")
        return result

    def _select_cut(self, work, counters):
        if self.requested_cut != "auto":
            return int(self.requested_cut), None
        T = S = 0
        for high in range(self.H):
            w = (work-self.L*high) % self.period
            c = self.b*(w//self.b)
            residues = set()
            for u in range(self.b):
                for q in range(self.b):
                    counters["selector_pair_terms"] += 1
                    rho = (c+q-u) % self.period
                    T += max(0, 1+(self.L-1-rho)//self.period)
                    residues.add(rho)
            for rho in residues:
                counters["selector_residue_terms"] += 1
                S += max(0, 1+(self.L-1-rho)//self.period)
        best = None
        for cut in self.candidate_cuts:
            counters["selector_candidate_terms"] += 1
            _, K, P = self.cut_parameters[cut]
            grouping = self.H*K*self.b*self.b
            coefficients = min(grouping*P, T)
            geometric = min(self.H*K*min(self.period, 2*self.b-1)*P, S)
            record = (grouping+coefficients, -cut, grouping, coefficients, geometric)
            if best is None or record < best:
                best = record
        return -best[1], dict(T=T, S=S, grouping=best[2], coefficients=best[3],
                              geometric=best[4], objective=best[0])

    def _row(self, work, counters):
        cut, selection = self._select_cut(work, counters)
        A, K, P = self.cut_parameters[cut]
        stride = self.period*P
        left, right = [], []
        for item in self.early_phases:
            counters["phase_partition_terms"] += 1
            (left if item[0] < cut else right).append(item)
        components = []
        for high in range(self.H):
            w = (work-self.L*high) % self.period
            c, p = self.b*(w//self.b), w % self.b
            late = self._phase(w, counters)
            for k in range(K):
                right_values = None
                if self.cache_right:
                    right_values = []
                    for q in range(self.b):
                        product = 1.+0j
                        for v, oracle in right:
                            argument = c+q-(1 << v)*(k//(1 << (v-cut)))
                            product *= self._early(oracle, argument, counters)
                            counters["row_phase_products"] += 1
                        right_values.append(product)
                        counters["right_cache_entries"] += 1
                pairs = {}
                for u in range(self.b):
                    for q in range(self.b):
                        counters["row_local_terms"] += 1
                        counters["group_pair_visits"] += 1
                        rho = (c+q-A*k-u) % self.period
                        pairs.setdefault(rho, []).append(
                            (u, q, complex(self.middle[p, q]*self.initial[u, 0])))
                for rho, local_pairs in sorted(pairs.items()):
                    N = max(0, 1+(A-1-rho)//self.period)
                    for z in range(min(P, N)):
                        a = rho+self.period*z
                        gamma = 0j
                        counters["geometric_components"] += 1
                        for u, q, amplitude in local_pairs:
                            counters["row_local_terms"] += 1
                            counters["coefficient_pair_visits"] += 1
                            product = right_values[q] if self.cache_right else 1.+0j
                            phases = left if self.cache_right else self.early_phases
                            for v, oracle in phases:
                                argument = (u+a % (1 << v) if v < cut else
                                            c+q-(1 << v)*(k//(1 << (v-cut))))
                                product *= self._early(oracle, argument, counters)
                                counters["row_phase_products"] += 1
                            gamma += amplitude*product
                            counters["row_phase_products"] += 1
                        gamma *= late
                        counters["row_phase_products"] += 1
                        if gamma != 0:
                            count = 1+(N-1-z)//P
                            start = self.L*high+A*k+a
                            if len(components) >= self.component_bound:
                                raise MemoryError("nested-phase component cap before append")
                            if not 0 <= start <= start+stride*(count-1) < self.Q:
                                raise ArithmeticError("nested-phase component leaves register")
                            components.append((start, count, complex(gamma)))
        counters["norm_terms"] += len(components)
        norm = math.fsum(count*_weight(a) for _, count, a in components)
        if not math.isfinite(norm):
            raise ArithmeticError("nonfinite nested-phase row norm")
        return dict(work=work, components=tuple(components), stride=stride, norm=norm,
                    cut=cut, selection=selection, cache_right=self.cache_right,
                    counters=counters)

    def stats(self):
        result = super().stats()
        d = len(self.early_phases)
        result.update(phase_positions=tuple(v for v, _ in self.early_phases),
                      requested_cut=self.requested_cut, candidate_cuts=self.candidate_cuts,
                      cache_right=self.cache_right, selector_pair_bound=self.selector_pair_bound,
                      selector_setup_bound=self.selector_setup_bound,
                      geometry_setup_phase_terms=self.geometry_setup_phase_terms,
                      phase_input_entries=d, row_phase_partition_terms_bound=d,
                      selector_residue_bound=self.H*min(self.period, 2*self.b-1)
                          if self.requested_cut == "auto" else 0,
                      selector_candidate_bound=len(self.candidate_cuts)
                          if self.requested_cut == "auto" else 0,
                      row_pair_bound=self.row_pair_bound,
                      coefficient_pair_bound=self.coefficient_pair_bound,
                      phase_queries_per_draw_bound=d*self.b+self.b*self.b+self.row_phase_query_bound,
                      column_phase_products_bound=d*self.b+self.b*self.b,
                      row_phase_products_bound=self.row_phase_product_bound,
                      phase_schedule_entries=d, right_cache_live_entries_bound=self.b
                          if self.cache_right else 0,
                      python_container_entry_bound=16*(self.component_bound+self.b*self.b+d)+256,
                      cover_selection="exact grouping+coefficient visits, weights (1,1); larger cut on ties",
                      selector_setup="H*b^2 pair terms, H*R residue terms, <=d+1 candidates; no phase calls",
                      row_local_terms_semantics="grouping plus coefficient visits; selector counted separately",
                      local_term_bound_semantics="column+row pairs, selector pairs/residues/candidates, phase partition",
                      phase_preparation="<=64 entries; stable position sort plus reported geometry phase inspections")
        return result
