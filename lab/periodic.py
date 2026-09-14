"""Known-order orbit defects via conserved or deterministically routed sectors.

For r=b*M and orbit labels j=b*m+p, measure U^b's eigenphase alpha uniformly
at the start. Every repeated b-by-b orbit block and every arithmetic power
preserves that sector. The remaining work dimension is b, even with multiple
noncommuting defects. This is standard symmetry reduction plus sequential
instrument contraction, not order discovery or a new factoring algorithm.
RoutedOrbitCircuit also permits supplied cell-character phases that permute
the sectors bijectively rather than conserving each one separately.
"""
from __future__ import annotations
import numpy as np
from lab.fourier_sampling import _ints, unit_phase
from lab.semiclassical import sequential_path


def lightcone_cover(period, width, support, radius):
    """Exact integer cover for a localized kick in a co-moving orbit frame.

    Suppose the pre-kick state is 2^(-width/2) sum_l |l>|phi_l>, with each
    normalized phi_l supported on j=l+delta (mod period), |delta|<=radius.
    A kick equal to identity outside the supplied orbit labels has hit mass
    at most the fraction of l whose residue lies within radius of support.
    For ANY such unitary kick, output TV^2 is at most min(1,4*hit_mass).

    Return exact numerator/denominator integers, not a floating-point error
    certificate. No orbit/prefix array or state propagation occurs. Circular
    interval merging costs O(d log d) integer operations and O(d) storage;
    integer bit complexity and verification of the support promise are extra.
    Caps match the periodic helper's known-order/index input contract.
    """
    if (not _ints(period, width, radius)
            or not 1 <= period <= np.iinfo(np.int64).max
            or not 0 <= width <= 63 or radius < 0
            or not isinstance(support, (list, tuple, np.ndarray))
            or isinstance(support, np.ndarray) and support.ndim != 1
            or not 0 <= len(support) <= 64):
        raise ValueError("invalid orbit cover dimensions or supplied support")
    period, width, radius = int(period), int(width), int(radius)
    if any(not _ints(p) or not 0 <= p < period for p in support):
        raise ValueError("support labels must be known orbit indices")
    labels = tuple(int(p) for p in support)
    if len(set(labels)) != len(labels):
        raise ValueError("support labels must be distinct")
    intervals = []
    span = 2*radius+1
    if labels and span >= period:
        intervals.append((0, period))
    else:
        for p in labels:
            lo = (p-radius) % period
            hi = lo+span
            if hi <= period:
                intervals.append((lo, hi))
            else:
                intervals.extend(((lo, period), (0, hi-period)))
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(hi, merged[-1][1]))
        else:
            merged.append((lo, hi))
    orbit_labels = sum(hi-lo for lo, hi in merged)
    length = 1 << width
    traversals, remainder = divmod(length, period)
    count = traversals*orbit_labels + sum(max(0, min(hi, remainder)-lo)
                                         for lo, hi in merged)
    return dict(period=period, prefix_width=width, radius=radius,
                support_size=len(labels), intervals=merged,
                covered_orbit_labels=orbit_labels, covered_prefix_labels=count,
                hit_probability_numerator=count, hit_probability_denominator=length,
                tv_squared_numerator=min(length, 4*count),
                tv_squared_denominator=length, orbit_table_entries=0,
                prefix_table_entries=0, count_arithmetic="exact integers")


def _nonnegative_floor_sum(n, modulus, multiplier, offset):
    """Sum floor((multiplier*q+offset)/modulus), q=0..n-1, and loop count.

    Standard Euclidean lattice-point recursion; compare AtCoder Library's
    CC0 floor_sum_unsigned (internal_math.hpp). Python integers here are
    exact, not that implementation's modulo-2^64 output. Internal nonnegative
    arguments only; callers validate/cap their input before reaching this.
    """
    total, iterations = 0, 0
    while True:
        iterations += 1
        quotient, multiplier = divmod(multiplier, modulus)
        total += quotient*n*(n-1)//2
        quotient, offset = divmod(offset, modulus)
        total += quotient*n
        height = multiplier*n + offset
        if height < modulus:
            return total, iterations
        n, offset = divmod(height, modulus)
        modulus, multiplier = multiplier, modulus


def uniform_prefix_certificate(period, width, prefix_bits, radius):
    """Sufficient exact support certificate for LOW inverse-QFT output bits.

    Supplied promise: before the terminating inverse QFT, the normalized
    conditional work state for exponent e has support in e+[-radius,radius]
    modulo known period. Controls initially have uniform amplitudes. Finite
    orbit-displacement mixers and any orbit-diagonal phases can establish
    this promise; this arithmetic-only helper does NOT inspect a circuit.

    Put H=2^prefix_bits and L=2^(width-prefix_bits). Count q in [1,H) with
    circular distance of L*q modulo period <=2*radius. No such q implies
    disjoint work supports for the high-input histories at each fixed low
    input, hence a uniform LOW-output prefix. A nonzero count is INCONCLUSIVE,
    not proof of nonuniformity. The zero-bit prefix is a vacuous certificate.
    Sampling the remaining output still requires correct prefix conditioning.

    Two Euclidean floor sums count the near-collisions in O(log period)
    integer iterations and constant many live integer variables, with no
    orbit/prefix enumeration. Integer bit arithmetic and verification of the
    supplied order/support promise are additional costs. Existing sampler
    defaults and propagation are unchanged; this returns no numerical law.
    """
    limit = (1 << 63)-1
    if (not _ints(period, width, prefix_bits, radius)
            or any(isinstance(v, (bool, np.bool_)) for v in (period, width, prefix_bits, radius))
            or not 1 <= period <= limit or not 0 <= width <= 63
            or not 0 <= prefix_bits <= width or not 0 <= radius <= limit):
        raise ValueError("invalid uniform-prefix dimensions or radius")
    period, width, prefix_bits, radius = map(int, (period, width, prefix_bits, radius))
    histories = 1 << prefix_bits
    step = (1 << (width-prefix_bits)) % period
    distance = 2*radius
    calls, iterations = 0, 0
    if histories == 1:
        collisions = 0
    elif distance >= period//2:
        collisions = histories-1
    else:
        # Count residues <=D and >=r-D, subtracting q=0. The two base
        # floor sums cancel; the intervals are disjoint in this branch.
        above_low, loops_low = _nonnegative_floor_sum(
            histories, period, step, period-distance-1)
        above_high, loops_high = _nonnegative_floor_sum(
            histories, period, step, distance)
        collisions = histories-1-above_low+above_high
        calls, iterations = 2, loops_low+loops_high
    if not 0 <= collisions < histories:
        raise ArithmeticError("invalid exact near-collision count")
    certified = collisions == 0
    return dict(period=period, exponent_width=width, prefix_bits=prefix_bits,
                radius=radius, high_history_count=histories,
                step_mod_period=step, near_collision_count=collisions,
                certified_uniform=certified, vacuous_prefix=(prefix_bits == 0),
                prefix_probability_numerator=1 if certified else None,
                prefix_probability_denominator=histories if certified else None,
                floor_sum_calls=calls, euclidean_iterations=iterations,
                orbit_table_entries=0, prefix_table_entries=0,
                arithmetic="exact Python integers; not a floating-state certificate")


class PeriodicOrbitCircuit:
    """Ascending controlled powers with repeated-block defects after given bits.

    defects maps insertion index s in [0,width] to one supplied b-by-b unitary,
    applied after s controls (0 means before arithmetic). At each insertion
    the same block acts on EVERY consecutive group of b orbit basis indices.
    Computational work labels are not these indices; the mapping/order is an
    external input. Orbit periodicity does not imply physical few-qubit locality.
    """
    def __init__(self, period, block_size, width, defects):
        if (not _ints(period, block_size, width)
                or not 1 <= period <= np.iinfo(np.int64).max
                or not 1 <= block_size <= 64 or period % block_size
                or not 0 <= width <= 63 or not isinstance(defects, dict)
                or len(defects) > width+1):
            raise ValueError("require known r divisible by block size, and bounded width/support")
        b, width = int(block_size), int(width)
        for s, W in defects.items():
            if (not _ints(s) or not 0 <= s <= width
                    or not isinstance(W, np.ndarray) or W.shape != (b, b)):
                raise ValueError("invalid defect insertion or supplied block shape")
        self.period, self.b, self.width = int(period), b, width
        self.sectors = self.period//b
        self.identity = np.eye(b, dtype=complex)
        self.defects = {}
        for s, W in defects.items():
            block = np.array(W, dtype=complex, copy=True)
            if (not np.all(np.isfinite(block)) or not np.allclose(
                    block.conj().T @ block, self.identity, rtol=0, atol=1e-12)):
                raise ValueError("periodic defect blocks must be finite unitaries")
            block.flags.writeable = False
            self.defects[int(s)] = block

    def shift(self, power, sector):
        """U^power within one sector, using integer phases (no float squaring)."""
        if (not _ints(power, sector) or power < 0 or not 0 <= sector < self.sectors):
            raise ValueError("invalid shift power or sector")
        power, sector = int(power), int(sector)
        result = np.zeros((self.b, self.b), dtype=complex)
        for p in range(self.b):
            wraps, target = divmod(p+power, self.b)
            result[target, p] = unit_phase(sector*wraps, self.sectors)
        return result

    def sector_path(self, sector, *, rng=None, output=None):
        if not _ints(sector) or not 0 <= sector < self.sectors:
            raise ValueError("coarse phase outside sector range")
        if (rng is None) == (output is None):
            raise ValueError("supply exactly one of rng or forced output")
        if output is not None and (not _ints(output) or not 0 <= output < 1 << self.width):
            raise ValueError("forced output outside exponent register")
        pairs = []
        for i in range(self.width):
            W = self.defects.get(i+1, self.identity)
            pairs.append((W, W @ self.shift(1 << i, sector)))
        initial = self.defects.get(0, self.identity)[:, 0]
        result = sequential_path(pairs, initial, rng=rng, output=output)
        result.update(coarse_eigenphase=int(sector), sector_probability=1/self.sectors,
                      joint_latent_output_probability=result["conditional_path_probability"]/self.sectors)
        return result

    def forced_joint(self, sector, output):
        """Joint p(alpha,y), NOT the marginal p(y); alpha is a COARSE phase."""
        if output is None:
            raise ValueError("forced output is required")
        return self.sector_path(sector, output=output)

    def sample(self, rng):
        sector = int(rng.integers(self.sectors))
        return self.sector_path(sector, rng=rng)

    def localized_kick_bound(self, split, support):
        """Bound replacing ONE extra localized kick by identity.

        The kick is after all background gates at insertion ``split``. Each
        earlier repeated block changes an orbit index by at most b-1, so the
        radius is m*(b-1). Identity blocks are conservatively counted too.
        The supplied support promise is NOT verified by this method; no kick
        matrix is an input. For an exact background sampler the returned
        rational TV^2 bound certifies omission in mathematical arithmetic.
        The float sampler's own numerical error must be added separately.
        """
        if not _ints(split) or not 0 <= split <= self.width:
            raise ValueError("kick insertion outside circuit")
        mixers = sum(s <= split for s in self.defects)
        result = lightcone_cover(self.period, int(split), support, mixers*(self.b-1))
        result.update(pre_kick_mixer_count=mixers, block_size=self.b,
                      kick_after_background_at_split=True,
                      verified_kick_support=False,
                      numerical_sampler_error_included=False)
        return result

    def stats(self):
        return dict(period=self.period, block_size=self.b, width=self.width,
                    sector_count=self.sectors, supplied_defect_count=len(self.defects),
                    owned_block_payload_bytes=(len(self.defects)+1)*self.b*self.b*16,
                    sampled_work_dimension=self.b, orbit_table_entries=0,
                    early_vector_entries=0, output_table_entries=0,
                    order_and_index_discovery_included=False,
                    path_storage="O(width*block_size^2) scalars; not process peak")


class RoutedOrbitCircuit(PeriodicOrbitCircuit):
    """Repeated blocks interleaved with known cell-character phase gates.

    At insertion s, apply D_q W_s, where D_q|b*m+p> equals
    exp(2*pi*i*q*m/M)|b*m+p>, M=r/b. ``sector_shifts`` maps s to integer q.
    D_q moves coarse alpha to alpha-q mod M WITHOUT branching it. Every
    initial sector follows a deterministic bijection, so initial uniform
    coarse dephasing is still valid after the final work trace. It is not
    necessary that each individual sector be conserved.

    Generic coherent sums of different routes are outside this contract.
    Order, indexing and the character/block description are supplied. Such
    an orbit-character phase is not automatically a physical one-qubit gate.
    This reuses the finite-work sequential instrument, not a new propagator.
    """
    def __init__(self, period, block_size, width, defects, sector_shifts):
        if (not _ints(width) or not 0 <= width <= 63
                or not isinstance(sector_shifts, dict) or len(sector_shifts) > width+1
                or any(not _ints(s, q) or not 0 <= s <= width
                       for s, q in sector_shifts.items())):
            raise ValueError("routing phases need bounded insertion indices and integer charges")
        super().__init__(period, block_size, width, defects)
        self.sector_shifts = {int(s): int(q) % self.sectors
                              for s, q in sector_shifts.items()}

    def sector_path(self, sector, *, rng=None, output=None):
        if not _ints(sector) or not 0 <= sector < self.sectors:
            raise ValueError("initial coarse phase outside sector range")
        if (rng is None) == (output is None):
            raise ValueError("supply exactly one of rng or forced output")
        if output is not None and (not _ints(output) or not 0 <= output < 1 << self.width):
            raise ValueError("forced output outside exponent register")
        initial_sector = int(sector)
        current = (initial_sector-self.sector_shifts.get(0, 0)) % self.sectors
        initial = self.defects.get(0, self.identity)[:, 0]
        pairs = []
        for i in range(self.width):
            W = self.defects.get(i+1, self.identity)
            pairs.append((W, W @ self.shift(1 << i, current)))
            current = (current-self.sector_shifts.get(i+1, 0)) % self.sectors
        result = sequential_path(pairs, initial, rng=rng, output=output)
        result.update(coarse_eigenphase=initial_sector,
                      final_coarse_eigenphase=current,
                      sector_probability=1/self.sectors,
                      joint_latent_output_probability=result["conditional_path_probability"]/self.sectors,
                      sector_routing="deterministic cell-character shifts")
        return result

    def localized_kick_bound(self, split, support):
        # D_q is diagonal in the orbit basis, hence adds no orbit displacement.
        result = super().localized_kick_bound(split, support)
        result["diagonal_routing_phases_add_zero_displacement"] = True
        return result

    def stats(self):
        result = super().stats()
        result.update(routing_gate_count=len(self.sector_shifts),
                      routing_integer_count=len(self.sector_shifts),
                      individual_sector_conservation_required=False,
                      initial_sector_mixture_valid=True)
        return result
