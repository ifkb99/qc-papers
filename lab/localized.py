"""Order/index-informed sampling for ONE finite-basis-support work defect.

V is identity outside d supplied orbit indices. This is not the promise that
V acts on d physical qubits, and order discovery/index lookup are not supplied.
No arrays of length r, 2^split or 2^width are created. Mathematical sampling
identities are exact; the implementation uses complex128/float64 and bounded
rejection attempts, raising on exhaustion rather than returning a fallback.
The theorem concerns uncapped rejection: discarding failed calls and keeping
only successes can reweight latent phases, so exhaustion must not be ignored.
"""
from __future__ import annotations
import math
import numpy as np

from lab.fourier_sampling import (_ints, unit_phase, geometric_sum,
                                  interval_path, progression_sample)
from lab.semiclassical import eigenphase_path


class LocalizedOrbitDefect:
    """A unitary block on at most 64 known orbit basis indices.

    unitary must already be an ndarray of shape (d,d), output rows / input
    columns. Its shape and the support cap are checked BEFORE numeric copying.
    The complement is the identity. Empty support represents the identity.
    """

    def __init__(self, period, indices, unitary):
        if not _ints(period) or not 1 <= period <= np.iinfo(np.int64).max:
            raise ValueError("period must be a positive int64-range integer")
        if not hasattr(indices, "__len__") or len(indices) > 64:
            raise ValueError("supply at most 64 known orbit indices")
        d = len(indices)
        if (not isinstance(unitary, np.ndarray) or unitary.shape != (d, d)
                or not all(_ints(p) and 0 <= p < period for p in indices)
                or len(set(indices)) != d):
            raise ValueError("invalid support or unitary block shape")
        self.period, self.indices = int(period), tuple(int(p) for p in indices)
        block = np.array(unitary, dtype=complex, copy=True)
        if (not np.all(np.isfinite(block))
                or not np.allclose(block.conj().T @ block, np.eye(d), rtol=0, atol=1e-12)):
            raise ValueError("defect block must be finite and unitary")
        block.flags.writeable = False
        self.unitary, self.d = block, d
        self.D, self.m = max(1, d), d+1

    def _phase_data(self, split, k):
        if not _ints(split, k) or not 0 <= split <= 63 or not 0 <= k < self.period:
            raise ValueError("invalid split or final eigenphase")
        split, k = int(split), int(k)
        L = 1 << split
        counts = tuple(0 if p >= L else 1+(L-1-p)//self.period for p in self.indices)
        phases = np.array([unit_phase(k*p, self.period) for p in self.indices])
        active = phases @ self.unitary
        delta = active - phases
        # Positive disjoint-support formula avoids subtracting O(L) quantities.
        norm = (L-sum(counts))/L + sum((n/L)*abs(h)**2 for n, h in zip(counts, active))
        weights = np.array([1.] + [(n/L)*abs(a)**2 for n, a in zip(counts, delta)])
        return dict(split=split, k=k, L=L, counts=counts, delta=delta,
                    norm=float(norm), weights=weights, T=float(weights.sum()))

    def phase_probability(self, split, eigenphase):
        return self._phase_data(split, eigenphase)["norm"] / self.period

    def _components(self, width, data, output):
        Q, L, k = 1 << width, data["L"], data["k"]
        b = [geometric_sum(L, k*Q-output*self.period, self.period*Q)/L]
        b.extend(a*unit_phase(-output*p, Q)*geometric_sum(n, -output*self.period, Q)/L
                 for p, n, a in zip(self.indices, data["counts"], data["delta"]))
        return np.array(b, dtype=complex)

    def _validate_path(self, width, split, output=None):
        if not _ints(width, split) or not 0 <= split <= width <= 63:
            raise ValueError("require 0<=split<=width<=63")
        if output is not None and (not _ints(output) or not 0 <= output < 1 << int(width)):
            raise ValueError("forced output outside exponent register")

    def _result(self, width, data, output, late_probability):
        mag = float(abs(self._components(width, data, output).sum())**2)
        norm = data["norm"]
        return dict(output=output, eigenphase=data["k"],
                    phase_probability=norm/self.period,
                    conditional_path_probability=late_probability*mag/norm if norm > 0 else None,
                    joint_latent_output_probability=late_probability*mag/self.period if norm > 0 else 0.,
                    expected_early_attempts=self.m*data["T"]/norm if norm > 0 else None,
                    component_count=self.m, late_steps=width-data["split"])

    def forced_joint(self, width, split, eigenphase, output):
        """Evaluate joint p(k,y), NOT the marginal p(y) (which sums over k)."""
        self._validate_path(width, split, output)
        if output is None:
            raise ValueError("forced output is required")
        width, split, output = int(width), int(split), int(output)
        data = self._phase_data(split, eigenphase)
        H = 1 << (width-split)
        late = eigenphase_path(self.period, width-split, (int(eigenphase)*data["L"]) % self.period,
                               output=output % H)
        return self._result(width, data, output, late["conditional_path_probability"])

    def sample(self, width, split, rng, *, max_attempts=100_000):
        """Sample k, late q, then early z by coherent-component rejection.

        Phase proposals average D. Early proposals average at most 5(d+1)
        OVER THE SAMPLED k distribution; this is NOT a pointwise-in-k bound.
        No asymptotic bit/precision guarantee follows from float64 timings.
        """
        self._validate_path(width, split)
        if not _ints(max_attempts) or max_attempts < 1:
            raise ValueError("rejection cap must be a positive integer")
        width, split = int(width), int(split)
        for phase_attempt in range(1, int(max_attempts)+1):
            k = int(rng.integers(self.period))
            data = self._phase_data(split, k)
            acceptance = data["norm"]/self.D
            if not math.isfinite(acceptance) or not 0 <= acceptance <= 1+1e-12:
                raise ArithmeticError("invalid phase rejection envelope")
            if rng.random() < min(1., acceptance):
                break
        else:
            raise RuntimeError("final-eigenphase rejection cap exhausted")
        L, H, Q = data["L"], 1 << (width-split), 1 << width
        late = eigenphase_path(self.period, width-split, (k*L) % self.period, rng=rng)
        q = late["output"]
        queries = 0
        for early_attempt in range(1, int(max_attempts)+1):
            component = int(rng.choice(self.m, p=data["weights"]/data["T"]))
            if component == 0:
                proposal = interval_path(split, L, q*self.period-k*Q, Q*self.period, rng=rng)
            else:
                ip = component-1
                proposal = progression_sample(split, self.indices[ip], self.period,
                                               data["counts"][ip], q, Q, rng)
            queries += proposal["marginal_queries"]
            output = q+H*proposal["output"]
            b = self._components(width, data, output)
            denom = self.m*float(np.vdot(b, b).real)
            if not math.isfinite(denom) or denom <= 0:
                raise ArithmeticError("zero/nonfinite weight at a proposed output")
            acceptance = float(abs(b.sum())**2)/denom
            if not math.isfinite(acceptance) or not 0 <= acceptance <= 1+1e-12:
                raise ArithmeticError("invalid coherent-component rejection envelope")
            if rng.random() < min(1., acceptance):
                result = self._result(width, data, output, late["conditional_path_probability"])
                result.update(phase_attempts=phase_attempt, early_attempts=early_attempt,
                              marginal_queries=queries)
                return result
        raise RuntimeError("early-component rejection cap exhausted")

    def stats(self):
        """Owned matrix payload, not process peak; scalar/list overhead excluded."""
        return dict(support_size=self.d, component_count=self.m,
                    matrix_payload_bytes=self.unitary.nbytes, orbit_table_entries=0,
                    early_vector_entries=0, output_table_entries=0,
                    expected_phase_attempts=self.D,
                    mean_early_attempt_bound=5*self.m,
                    scalar_workspace="O(d) beyond supplied/copied d-by-d matrix",
                    order_and_index_discovery_included=False)
