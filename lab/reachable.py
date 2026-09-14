"""Reachable-state conditional instruments; no dense scratch-space tables.

The instrument/commutation scope is the same as lab.semiclassical. Branch
actions are vectorized callables on integer basis labels. Only reached labels
are stored. No tolerance truncation is used: only exactly zero amplitudes are
removed. Floating-point evaluation is not symbolically exact.

This is a sparse work-state baseline, NOT a sub-orbit compression algorithm.
For clean order finding it uses O(r) amplitudes in the worst case; an optional
transition cache can use O(d*r) entries for d distinct block multipliers.
"""
from __future__ import annotations
import math
import numpy as np
import walsh
from toffoli_arith import ToffoliModExp


class CompiledBranch:
    """Memoized selected-input replay of one existing controlled circuit.

    Both branches are replayed; the inactive branch is NOT assumed to be
    identity on arbitrary scratch states. Caches contain only queried labels.
    """
    def __init__(self, circuit, control, bit):
        if circuit.n != control + 1 or bit not in (0, 1):
            raise ValueError("requires a highest-bit control and bit 0 or 1")
        self.circuit, self.control, self.bit = circuit, control, bit
        self.cache = {}
        self.replayed_states = 0
        self.replay_batches = 0

    def __call__(self, ids):
        ids = np.asarray(ids)
        if (ids.ndim != 1 or not np.issubdtype(ids.dtype, np.integer)
                or (ids.size and (int(ids.min()) < 0
                                  or int(ids.max()) >= 1 << self.control))):
            raise ValueError("invalid work basis labels")
        missing = np.array(sorted(set(map(int, ids)) - self.cache.keys()), dtype=np.int64)
        if missing.size:
            images = walsh.classical_images(
                self.circuit, missing | (self.bit << self.control))
            if not np.all((images >> self.control) == self.bit):
                raise ValueError("compiled branch changes its control")
            images &= (1 << self.control) - 1
            self.cache.update(zip(map(int, missing), map(int, images)))
            self.replayed_states += missing.size
            self.replay_batches += 1
        return np.array([self.cache[int(x)] for x in ids], dtype=np.int64)


def compiled_pairs(N, a, width):
    """Compile distinct existing multipliers, without order/orbit enumeration.

    Constant generation uses repeated modular squaring. Setup is proportional
    to the emitted circuits, not the Hilbert-space dimension. No branch images
    are evaluated until the returned actions are queried.
    """
    if not isinstance(width, int) or width < 0:
        raise ValueError("width must be a nonnegative integer")
    me = ToffoliModExp(N, a, n_exp=1)
    if me.n_qubits > 63:
        raise ValueError("compiled replay currently uses int64 basis labels")
    cache, pairs = {}, []
    multiplier = a % N
    for _ in range(width):
        if multiplier not in cache:
            qc = me.u_a(me.exp[0], multiplier)
            cache[multiplier] = (CompiledBranch(qc, me.exp[0], 0),
                                 CompiledBranch(qc, me.exp[0], 1))
        pairs.append(cache[multiplier])
        multiplier = multiplier * multiplier % N
    return me, pairs


def action_stats(pairs):
    actions = {id(action): action for pair in pairs for action in pair}.values()
    actions = list(actions)
    circuits = {id(action.circuit): action.circuit for action in actions}
    return dict(cached_transitions=sum(len(a.cache) for a in actions),
                replayed_states=int(sum(a.replayed_states for a in actions)),
                replay_batches=sum(a.replay_batches for a in actions),
                distinct_blocks=len(circuits),
                logical_gates=sum(len(qc.logical) for qc in circuits.values()))


def sparse_children(ids, amps, pair, phase):
    """Return common reached labels and two unnormalized amplitude arrays.

    Caller supplies distinct labels and finite amplitudes. Injectivity is
    checked on the queried set; global bijectivity is a branch contract.
    """
    images = [np.asarray(action(ids)) for action in pair]
    if len(images) != 2 or any(
        img.shape != ids.shape or not np.issubdtype(img.dtype, np.integer)
        or (img.size and (int(img.min()) < 0 or int(img.max()) > np.iinfo(np.int64).max))
        or np.unique(img).size != ids.size for img in images
    ):
        raise ValueError("branch must be an injective integer basis action")
    labels = np.union1d(*images).astype(np.int64)
    left = np.zeros(labels.size, dtype=complex)
    right = np.zeros(labels.size, dtype=complex)
    left[np.searchsorted(labels, images[0])] = amps
    right[np.searchsorted(labels, images[1])] = amps * np.exp(1j * phase)
    return labels, ((left + right) / 2, (left - right) / 2)


def sparse_path(pairs, initial, *, rng=None, output=None, feedback=True):
    """Sample one path OR evaluate a forced output (including rare outcomes).

    initial maps work basis labels to amplitudes. Return per-step diagnostics,
    including the UNION support before selecting a measurement outcome. Counts
    are not process peak bytes: arrays, gate lists and caches coexist.
    """
    if (rng is None) == (output is None):
        raise ValueError("supply exactly one of rng or forced output")
    if output is not None and (not isinstance(output, (int, np.integer))
                              or not 0 <= output < 1 << len(pairs)):
        raise ValueError("output outside the exponent register")
    if not initial or any(not isinstance(x, (int, np.integer))
                          or not 0 <= x <= np.iinfo(np.int64).max for x in initial):
        raise ValueError("invalid initial basis labels")
    ids = np.array(sorted(initial), dtype=np.int64)
    amps = np.array([initial[int(x)] for x in ids], dtype=complex)
    if (not np.all(np.isfinite(amps))
            or not np.isclose(np.vdot(amps, amps).real, 1., atol=1e-12, rtol=0)):
        raise ValueError("initial state must be finite and normalized")
    nonzero = amps != 0
    ids, amps = ids[nonzero], amps[nonzero]
    prefix, path_probability, log_probability = 0, 1., 0.
    peak_stored, peak_union, profile = ids.size, ids.size, []
    for step, pair in enumerate(reversed(pairs)):
        phase = -np.pi * prefix / (1 << step) if feedback else 0.
        labels, children = sparse_children(ids, amps, pair, phase)
        weights = np.array([np.vdot(v, v).real for v in children])
        total = float(weights.sum())
        if not math.isfinite(total) or total <= 0:
            raise ArithmeticError("invalid instrument normalization")
        probs = weights / total
        bit = int(rng.random() >= probs[0]) if output is None else (int(output) >> step) & 1
        probability = float(probs[bit])
        prefix |= bit << step
        path_probability *= probability
        peak_union = max(peak_union, labels.size)
        if probability == 0:
            return dict(output=int(output) if output is not None else prefix,
                        path_probability=0., log_probability=None,
                        peak_stored_amplitudes=int(peak_stored),
                        peak_union_labels=int(peak_union), profile=profile,
                        zero_probability_step=step, final_ids=[], final_amplitudes=[])
        log_probability += math.log(probability)
        child = children[bit] / np.sqrt(weights[bit])
        nonzero = child != 0  # no tolerance cutoff
        ids, amps = labels[nonzero], child[nonzero]
        peak_stored = max(peak_stored, ids.size)
        profile.append(dict(step=step, exponent_bit=len(pairs)-1-step,
                            stored=int(ids.size), union=int(labels.size),
                            occupied_1e10=int(np.count_nonzero(np.abs(amps) > 1e-10)),
                            probability=probability, normalization_error=abs(total-1)))
    return dict(output=prefix, path_probability=path_probability,
                log_probability=log_probability, peak_stored_amplitudes=int(peak_stored),
                peak_union_labels=int(peak_union), profile=profile,
                final_ids=ids.tolist(),
                final_amplitudes=[[float(z.real), float(z.imag)] for z in amps])
