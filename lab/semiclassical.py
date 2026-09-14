"""Conditional instruments for a terminating inverse QFT.

The later sequential_path API contracts forward unmeasured states and backward
measurement effects; unlike the original permutation path APIs, it preserves
arbitrary noncommuting time order. See its separate input/allocation contract.

These are work-state updates, not a replacement Pauli propagator. A pair of
permutations (P0,P1) is extracted from an existing controlled arithmetic block;
measuring its |+> control after phase feedback gives
M_b=(P0+(-1)^b exp(i phase) P1)/2. Interference must be retained.

Pairs are supplied in increasing exponent-bit order. Measurements run in the
reverse order; the first measured bit is the least-significant output bit.
This defines an iterative instrument sequence. Equivalence to a coherent
controlled-block circuit followed by inverse QFT additionally requires the
reordered arithmetic actions to commute on the reachable work subspace. This
holds for valid modular multiplication, NOT arbitrary full-scratch branch maps.
The dense work vector is exponential in work width. Full-distribution
enumeration is exponential in exponent width; sample() follows one path only.
"""
from __future__ import annotations
import numpy as np


def validate_inputs(pairs, initial):
    initial = np.asarray(initial, dtype=complex)
    if initial.ndim != 1 or not initial.size or not np.all(np.isfinite(initial)):
        raise ValueError("initial must be a nonempty finite vector")
    if not np.isclose(np.vdot(initial, initial).real, 1., atol=1e-12, rtol=0):
        raise ValueError("initial must be normalized")
    seen = set()
    for pair in pairs:
        if len(pair) != 2:
            raise ValueError("each instrument needs two permutations")
        for perm in pair:
            if id(perm) in seen:
                continue
            perm = np.asarray(perm)
            if (perm.shape != initial.shape or not np.issubdtype(perm.dtype, np.integer)
                    or not np.array_equal(np.sort(perm), np.arange(initial.size))):
                raise ValueError("branch map is not a permutation of the work space")
            seen.add(id(perm))
    return initial


def conditional_children(psi, pair, phase):
    """Two UNNORMALIZED postmeasurement work states. Caller validates maps."""
    left, right = np.empty_like(psi), np.empty_like(psi)
    left[pair[0]] = psi
    right[pair[1]] = psi
    right *= np.exp(1j * phase)
    return (left + right) / 2, (left - right) / 2


def distribution(pairs, initial, *, feedback=True, max_width=10):
    """Enumerate the iterative instrument, subject to the module's scope."""
    if len(pairs) > max_width:
        raise ValueError("full distribution exceeds the explicit enumeration limit")
    initial = validate_inputs(pairs, initial)
    frontier = [(0, 1., initial)]
    profile = []
    for step, pair in enumerate(reversed(pairs)):
        next_frontier = []
        max_support, max_normalization_error = 0, 0.
        for prefix, weight, psi in frontier:
            phase = -np.pi * prefix / (1 << step) if feedback else 0.
            children = conditional_children(psi, pair, phase)
            probs = [float(np.vdot(v, v).real) for v in children]
            max_normalization_error = max(max_normalization_error, abs(sum(probs)-1))
            for bit, (child, prob) in enumerate(zip(children, probs)):
                if prob == 0:  # no numerical amplitude/probability truncation
                    continue
                child = child / np.sqrt(prob)
                max_support = max(max_support, int(np.count_nonzero(np.abs(child) > 1e-10)))
                next_frontier.append((prefix | (bit << step), weight * prob, child))
        frontier = next_frontier
        profile.append(dict(step=step, branches=len(frontier), max_active_amplitudes=max_support,
                            normalization_error=max_normalization_error))
    probabilities = np.zeros(1 << len(pairs))
    for output, weight, _ in frontier:
        probabilities[output] = weight
    return probabilities, profile


def sample(pairs, initial, rng, *, validated=False):
    psi = np.asarray(initial, dtype=complex).copy() if validated else validate_inputs(pairs, initial).copy()
    prefix, path_probability, peak_support = 0, 1., int(np.count_nonzero(psi))
    for step, pair in enumerate(reversed(pairs)):
        children = conditional_children(psi, pair, -np.pi * prefix / (1 << step))
        weights = np.array([np.vdot(v, v).real for v in children])
        probs = weights / weights.sum()
        bit = int(rng.random() >= probs[0])
        psi = children[bit] / np.sqrt(weights[bit])
        prefix |= bit << step
        path_probability *= float(probs[bit])
        peak_support = max(peak_support, int(np.count_nonzero(np.abs(psi) > 1e-10)))
    return dict(output=prefix, path_probability=path_probability,
                peak_active_amplitudes=peak_support, stored_amplitudes=psi.size,
                steps=len(pairs))


def order_finding_probability(output, width, period):
    """Independent geometric-series reference for ONE ideal output probability.

    Period is supplied only to the reference, never to the conditional sampler.
    Integer modular reduction keeps large phase arguments out of sin().
    """
    if width < 0 or period < 1 or not 0 <= output < (1 << width):
        raise ValueError("invalid output, width or period")
    Q = 1 << width
    quotient, extra = divmod(Q, period)
    remainder = (output * period) % Q

    def geometric_norm_squared(length):
        if remainder == 0:
            return float(length * length)
        numerator = (remainder * length) % Q
        num = np.sin(np.pi * min(numerator, Q-numerator) / Q)
        den = np.sin(np.pi * min(remainder, Q-remainder) / Q)
        return float((num / den)**2)

    return ((period-extra)*geometric_norm_squared(quotient)
            + extra*geometric_norm_squared(quotient+1)) / Q**2


def eigenphase_path(period, width, eigenphase, *, rng=None, output=None,
                    input_phases=None):
    """Order-INFORMED classical latent-eigenphase baseline, not order discovery.

    Choose eigenphase k uniformly in range(period) OUTSIDE this function.
    The eigenvalue is exp(2*pi*i*k/period). Conditional on it, the iterative
    measurement is a sequence of scalar Bernoulli updates, with no work vector.
    The returned probability is CONDITIONAL on k, NOT the marginal output
    probability. Averaging over k gives ideal order-finding statistics.
    Integer phase reduction avoids huge floating-point angles. Complexity is
    O(width) scalar updates, not a constant bit-complexity or precision claim.
    Optional input_phases[i] applies diag(1, exp(i*phase)) to exponent bit i
    BEFORE the inverse QFT. This covers product exponent phase kicks, not
    arbitrary noncommuting work operations. A work-to-exponent replacement
    requires its own proof on the actual input. The default is ideal arithmetic.
    """
    if (not all(isinstance(x, (int, np.integer)) for x in (period, width, eigenphase))
            or period < 1 or width < 0 or not 0 <= eigenphase < period):
        raise ValueError("invalid period, width or eigenphase")
    if (rng is None) == (output is None):
        raise ValueError("supply exactly one of rng or forced output")
    if output is not None and (not isinstance(output, (int, np.integer))
                              or not 0 <= output < 1 << width):
        raise ValueError("invalid forced output")
    period, width, eigenphase = int(period), int(width), int(eigenphase)
    phases = None
    if input_phases is not None:
        if np.iscomplexobj(input_phases):
            raise ValueError("input phases must be real")
        phases = np.asarray(input_phases, dtype=float)
        if phases.shape != (width,) or not np.all(np.isfinite(phases)):
            raise ValueError("input phases must be a finite length-width vector")
    prefix, probability = 0, 1.
    for step in range(width):
        exponent_bit = width-1-step
        angle = (2*np.pi*((eigenphase * (1 << exponent_bit)) % period)/period
                 - np.pi*prefix/(1 << step)
                 + (0. if phases is None else phases[exponent_bit]))
        weights = np.array([np.cos(angle/2)**2, np.sin(angle/2)**2])
        probs = weights / weights.sum()
        bit = int(rng.random() >= probs[0]) if output is None else (int(output) >> step) & 1
        prefix |= bit << step
        probability *= float(probs[bit])
    return dict(output=prefix, eigenphase=eigenphase,
                conditional_path_probability=probability, steps=width)


def _qft_branch_operators(B0, B1, phase):
    """Shared inverse-QFT branch ordering for float and verified instruments."""
    return (B0+phase*B1)/2, (B0-phase*B1)/2


def _density_forward_step(B0, B1, rho, *, multiply, adjoint):
    """Density form of sequential_path's forward factor update, without QR.

    Arithmetic callbacks preserve the same finite-work channel; verified
    execution retains enclosures instead of using float QR/renormalization.
    """
    return (multiply(multiply(B0,rho),adjoint(B0))
            + multiply(multiply(B1,rho),adjoint(B1)))/2


def sequential_path(pairs, initial, *, rng=None, output=None,
                    max_payload_bytes=32 << 20):
    """Terminating inverse-QFT sampling for a sequential finite-work circuit.

    Each initially |+> control i acts ONCE, in ascending i order, via supplied
    unitary work matrices (B_i0,B_i1). Work-only gates following that control
    are folded into BOTH branches on the left. Work is traced at the end.
    No commutation/reversal of these gates is assumed. This is a finite-work
    instrument contraction, not a new gate parser or circuit propagator.

    Forward factors R_i represent the unmeasured work state rho_i=R_i R_i^dag.
    Backward factors C represent an effect E=C^dag C. Child weights are
    ||C K_bit R_i||_F^2, K_bit=(B_i0+(-1)^bit exp(i phase)B_i1)/2.
    QR compresses a b-by-2b factor to b-by-b without truncating singular values.
    Positive squared norms avoid negative probabilities from density/effect
    roundoff. All math identities are exact; the code remains complex128.

    Caps: width<=63, work dimension<=64, conservative dense payload estimate
    checked BEFORE conversion/copying. O(width*b^3) arithmetic and O(width*b^2)
    scalar storage; integer/precision costs and supplied-matrix construction
    are additional. The return probability is conditional on these inputs.
    """
    from lab.fourier_sampling import _ints, unit_phase
    if (not hasattr(pairs, "__len__") or not 0 <= len(pairs) <= 63
            or not isinstance(initial, np.ndarray) or initial.ndim != 1
            or not 1 <= initial.size <= 64
            or not _ints(max_payload_bytes) or max_payload_bytes < 1):
        raise ValueError("invalid sequential dimensions or allocation budget")
    width, dim = len(pairs), int(initial.size)
    estimate = 16*(5*width+16)*dim*dim
    if estimate > max_payload_bytes:
        raise ValueError("sequential instrument exceeds conservative payload budget")
    if (rng is None) == (output is None):
        raise ValueError("supply exactly one of rng or forced output")
    if output is not None and (not _ints(output) or not 0 <= output < 1 << width):
        raise ValueError("forced output outside exponent register")
    for pair in pairs:
        if (not hasattr(pair, "__len__") or len(pair) != 2
                or any(not isinstance(B, np.ndarray) or B.shape != (dim, dim) for B in pair)):
            raise ValueError("each branch must already be a dimension-matched ndarray")
    psi = np.asarray(initial, dtype=complex)
    if not np.all(np.isfinite(psi)) or abs(float(np.vdot(psi, psi).real)-1) > 1e-12:
        raise ValueError("initial work state must be finite and normalized")
    identity = np.eye(dim, dtype=complex)
    branches = []
    for pair in pairs:
        converted = tuple(np.asarray(B, dtype=complex) for B in pair)
        for B in converted:
            if (not np.all(np.isfinite(B))
                    or not np.allclose(B.conj().T @ B, identity, rtol=0, atol=1e-12)):
                raise ValueError("sequential branches must be finite unitaries")
        branches.append(converted)
    R = np.zeros((dim, dim), dtype=complex)
    R[:, 0] = psi
    roots, normalization_error = [R], 0.
    for B0, B1 in branches:
        combined = np.concatenate((B0 @ R, B1 @ R), axis=1)/np.sqrt(2.)
        _, triangular = np.linalg.qr(combined.conj().T, mode="reduced")
        R = triangular.conj().T
        mass = float(np.vdot(R, R).real)
        normalization_error = max(normalization_error, abs(mass-1))
        if not np.isfinite(mass) or abs(mass-1) > 1e-10:
            raise ArithmeticError("invalid forward sequential normalization")
        R = R/np.sqrt(mass)
        roots.append(R)
    C, prefix, probability = identity/np.sqrt(dim), 0, 1.
    measured = 0
    for step, i in enumerate(range(width-1, -1, -1)):
        B0, B1 = branches[i]
        phase = unit_phase(-prefix, 1 << (step+1))
        children = tuple(C @ K for K in _qft_branch_operators(B0,B1,phase))
        weights = np.array([float(np.vdot(M, M).real)
                            for M in (child @ roots[i] for child in children)])
        total = float(weights.sum())
        if not np.all(np.isfinite(weights)) or total <= 0:
            raise ArithmeticError("zero/nonfinite sequential conditional weight")
        probs = weights/total
        bit = int(rng.random() >= probs[0]) if output is None else (int(output) >> step) & 1
        prefix |= bit << step
        probability *= float(probs[bit])
        measured += 1
        if weights[bit] == 0:
            prefix = int(output)
            break
        # Effect scale cancels from conditional probabilities. Normalize by
        # its own norm, not a tiny branch probability, to avoid blow-up.
        child = children[bit]
        C = child/np.sqrt(float(np.vdot(child, child).real))
    return dict(output=prefix, conditional_path_probability=probability,
                steps=measured, forward_steps=width, matrix_dimension=dim,
                forward_factor_payload_bytes=16*(width+1)*dim*dim,
                branch_matrix_payload_bytes=32*width*dim*dim,
                conservative_payload_estimate_bytes=estimate,
                forward_normalization_error=normalization_error)
