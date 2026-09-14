"""Final-eigenphase sampling with a bounded coherent early-control prefix.

This is an order- and orbit-index-INFORMED sampler, not order discovery or a
generic gate propagator. A caller supplies sparse normalized columns for
V|a^l>, expressed as known powers a^j. Arbitrary computational labels are NOT
orbit indices; finding their discrete logarithms is not performed here.

For at most D entries per column, h_k(l)=sum_j alpha_lj exp(2*pi*i*k*j/r)
gives final eigenphase weight ||h_k||²/(L*r), L=2^split. Uniform-k rejection
accepts with probability ||h_k||²/(L*D). Expected attempts are D in exact
arithmetic. After accepting k, retain only L early-control amplitudes, sample
the later controls with the existing scalar instrument, then take an L-point
FFT with their phase feedback. No orbit-sized array is created internally.
"""
from __future__ import annotations
import math
import numpy as np
from lab.semiclassical import eigenphase_path


class SparseOrbitPrefix:
    """Stream sparse columns; allocate O(L+D) temporary numerical data.

    column(l) returns (distinct orbit_indices, complex amplitudes), including
    all nonzero amplitudes, for 0<=l<L. Every call is checked. Columns must be
    deterministic and normalized; their being images of one physical unitary
    is the caller's independently validated contract. No columns are cached.
    Oracle-owned storage/setup are additional resources and must be reported.
    """

    def __init__(self, period, split, column, max_terms, *, max_prefix=4096):
        if (not all(isinstance(v, (int, np.integer)) for v in (period, split, max_terms))
                or not 1 <= period <= np.iinfo(np.int64).max
                or not 0 <= split <= 12 or max_terms < 1 or not callable(column)):
            raise ValueError("invalid period, split, sparse bound or column oracle")
        self.period, self.split, self.max_terms = int(period), int(split), int(max_terms)
        self.size = 1 << self.split
        if self.size > max_prefix:
            raise ValueError("early prefix exceeds its allocation budget")
        self.column = column
        self.column_calls = 0
        self.phase_rows = 0
        self.peak_column_entries = 0
        for l in range(self.size):
            self._column(l)

    def _column(self, l):
        ids, amps = self.column(l)
        ids, amps = np.asarray(ids), np.asarray(amps, dtype=complex)
        self.column_calls += 1
        if (ids.ndim != 1 or not np.issubdtype(ids.dtype, np.integer)
                or amps.shape != ids.shape or not 1 <= ids.size <= self.max_terms
                or int(ids.min()) < 0 or int(ids.max()) >= self.period
                or np.unique(ids).size != ids.size or not np.all(np.isfinite(amps))
                or abs(float(np.vdot(amps, amps).real)-1) > 1e-12):
            raise ValueError("column must have bounded distinct indices and normalized finite amplitudes")
        self.peak_column_entries = max(self.peak_column_entries, int(ids.size))
        return ids, amps

    def phase_row(self, eigenphase):
        """Unnormalized early amplitudes with the common 1/sqrt(r) omitted."""
        if (not isinstance(eigenphase, (int, np.integer))
                or not 0 <= eigenphase < self.period):
            raise ValueError("eigenphase outside orbit")
        k = int(eigenphase)
        row = np.empty(self.size, dtype=complex)
        for l in range(self.size):
            ids, amps = self._column(l)
            # Python-integer modular multiplication before conversion avoids
            # overflow and huge floating-point phase arguments.
            fractions = np.array([((k*int(j)) % self.period)/self.period for j in ids])
            row[l] = np.dot(amps, np.exp(2j*np.pi*fractions))
        self.phase_rows += 1
        return row

    def phase_probability(self, eigenphase):
        row = self.phase_row(eigenphase)
        return float(np.vdot(row, row).real) / self.size / self.period

    def draw_phase(self, rng, *, max_attempts=100_000):
        """Rejection sample k and its row; cap exhaustion raises, never biases."""
        if not isinstance(max_attempts, int) or max_attempts < 1:
            raise ValueError("max_attempts must be positive")
        for attempt in range(1, max_attempts+1):
            k = int(rng.integers(self.period))
            row = self.phase_row(k)
            acceptance = float(np.vdot(row, row).real) / self.size / self.max_terms
            if not math.isfinite(acceptance) or not 0 <= acceptance <= 1+1e-12:
                raise ArithmeticError("invalid sparse Fourier rejection envelope")
            # Min corrects only a possible rounding overshoot of a unit bound.
            if rng.random() < min(1., acceptance):
                return k, row, attempt
        raise RuntimeError("final-eigenphase rejection cap exhausted")

    def _conditional_path(self, width, k, row, *, rng=None, output=None):
        if (not isinstance(width, (int, np.integer)) or not self.split <= width <= 63
                or (rng is None) == (output is None)):
            raise ValueError("require split<=width<=63 and exactly one of rng/output")
        width = int(width)
        Q, H = 1 << width, 1 << (width-self.split)
        if output is not None and (not isinstance(output, (int, np.integer))
                                   or not 0 <= output < Q):
            raise ValueError("forced output outside exponent register")
        norm2 = float(np.vdot(row, row).real)
        phase_probability = norm2 / self.size / self.period
        if norm2 == 0:
            return dict(output=int(output) if output is not None else None, eigenphase=k,
                        phase_probability=0., conditional_path_probability=None,
                        joint_latent_output_probability=0., early_amplitudes=self.size)
        late_args = dict(rng=rng) if output is None else dict(output=int(output) % H)
        late = eigenphase_path(self.period, width-self.split,
                               (k*self.size) % self.period, **late_args)
        q = late["output"]
        feedback = np.exp(-2j*np.pi*np.array([((q*l) % Q)/Q for l in range(self.size)]))
        transformed = np.fft.fft(row*feedback/np.sqrt(norm2))/np.sqrt(self.size)
        weights = np.abs(transformed)**2
        total = float(weights.sum())
        if not math.isfinite(total) or abs(total-1) > 1e-10:
            raise ArithmeticError("invalid early-prefix normalization")
        probabilities = weights / total
        z = int(rng.choice(self.size, p=probabilities)) if output is None else int(output) // H
        conditional = late["conditional_path_probability"] * float(probabilities[z])
        return dict(output=q+H*z, eigenphase=k, phase_probability=phase_probability,
                    conditional_path_probability=conditional,
                    joint_latent_output_probability=phase_probability*conditional,
                    early_amplitudes=self.size, late_steps=width-self.split,
                    normalization_error=abs(total-1))

    def forced_joint(self, width, eigenphase, output):
        """A joint p(k,y), NOT the marginal output p(y); sum over k to validate."""
        row = self.phase_row(eigenphase)  # validate before coercing an invalid label
        return self._conditional_path(width, int(eigenphase), row, output=output)

    def sample(self, width, rng, *, max_attempts=100_000):
        """One marginal output sample via a sampled final latent eigenphase."""
        # Validate before any potentially expensive rejection work.
        if not isinstance(width, (int, np.integer)) or not self.split <= width <= 63:
            raise ValueError("require split<=width<=63")
        k, row, attempts = self.draw_phase(rng, max_attempts=max_attempts)
        result = self._conditional_path(width, k, row, rng=rng)
        result["rejection_attempts"] = attempts
        return result

    def stats(self):
        return dict(column_calls=self.column_calls, phase_rows=self.phase_rows,
                    early_amplitudes=self.size, phase_row_payload_bytes=16*self.size,
                    peak_column_entries=self.peak_column_entries,
                    orbit_table_entries=0, supplied_support_bound=self.max_terms,
                    oracle_storage_included=False)
