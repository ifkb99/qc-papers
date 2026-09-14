"""Does a short-prefix eigenphase sampler survive one orbit-preserving mixer?

DERIVED BEFORE MEASUREMENT.  After the first s ascending-power controls, put
one orbit-preserving work unitary V, and then apply the remaining controls and
the inverse QFT.  For a final work eigenphase k, with U|phi_k>=lambda_k|phi_k>,

    c[k,l] = <phi_k|V|a^l>,
    w[k] = (1/L) sum_l |c[k,l]|^2,

and the conditional control amplitudes at e=l+Lh are
    c[k,l]/sqrt(sum_l |c[k,l]|^2) * lambda_k**(L*h)/sqrt(H).
The k mixture is therefore weighted by w[k], not uniformly weighted in
general.  This experiment checks that formula against the existing dense
spectral effects reference, while using only the existing control-only
state-vector inverse QFT for the independent route.

P1: for r=3 and r=6, s=1 and s=2, seeded arbitrary complex orbit unitaries
    give normalized weights and the weighted conditional sampler exactly
    reproduces lab.spectral.single_defect_effects.
P2: the physical N=7,a=3 mixer of |1> and |5> at theta=pi/2 obeys the same
    weighted formula, including the nonuniform-or-uniform status of its k
    weights (checked before using the uniform-k control).
P3: for every tested case, each conditional control state is normalized and
    the full conditional output factors as y=q+H*z with the stated feedback.
P4: for generic normalized complex D-sparse orbit columns, the rejection
    acceptance a_k=||h_k||^2/(L*D) is at most one and averages exactly 1/D
    over uniform k (Parseval audit, not a statistical estimate).
P5: the independently implemented SparseOrbitPrefix API agrees with the
    coherent control-only route for every forced joint output and its reported
    phase probabilities sum to the same marginal distribution.
C1: replacing w[k] by uniform 1/r must fail whenever the measured weights are
    nonuniform.  This is the must-fail guard against the old initial-eigenphase
    mixture.
C2: erasing the complex phases of c[k,l] must fail for a seeded mixer.  This
    guards against treating the early prefix as a classical probability table.

Optional analytic audit: if every V|a^l> has at most D known orbit-index
amplitudes, then w[k] <= D/r and rejection from uniform k has expected cost at
most D.  Finding those orbit indices can itself hide a discrete-log/orbit
enumeration cost, so this is not a speedup claim.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_prefix_probe
"""
from __future__ import annotations

import time

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.prefix import SparseOrbitPrefix
from lab.spectral import single_defect_effects
import statevec


WIDTH = 4
THETA = np.pi / 2


def orbit_fourier(period: int) -> np.ndarray:
    """Columns are U-forward eigenvectors in the orbit-label basis."""
    j = np.arange(period)
    return np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)


def seeded_unitary(period: int, seed: int) -> np.ndarray:
    """Deterministic dense complex unitary, with no circuit propagation."""
    rng = np.random.default_rng(seed)
    z = rng.normal(size=(period, period)) + 1j * rng.normal(size=(period, period))
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    q *= np.where(np.abs(diagonal) > 0, np.conj(diagonal) / np.abs(diagonal), 1.)
    return q


def physical_mixer(orbit: np.ndarray, theta: float) -> np.ndarray:
    """exp(-i theta X_{|1>,|5>}/2) in the cyclic orbit basis."""
    lookup = {int(label): i for i, label in enumerate(orbit)}
    i, j = lookup[1], lookup[5]
    # Rotate only the selected two-dimensional subspace; leave all other
    # orbit labels untouched (cos(theta/2) times the full identity would not
    # be unitary when the generator is a partial swap).
    result = np.eye(len(orbit), dtype=complex)
    result[i, i] = result[j, j] = np.cos(theta / 2)
    result[i, j] = result[j, i] = -1j * np.sin(theta / 2)
    return result


def exact_outputs(v_eigen: np.ndarray, period: int, split: int) -> np.ndarray:
    effects = single_defect_effects(period, WIDTH, split, v_eigen)
    # Initial |1> has equal amplitudes in the Fourier eigenbasis.
    return effects.sum(axis=(1, 2)).real / period


def prefix_data(v_orbit: np.ndarray, period: int, split: int):
    """Return weights, conditional control distributions and early c arrays."""
    fourier = orbit_fourier(period)
    v_eigen = fourier.conj().T @ v_orbit @ fourier
    length = 1 << split
    height = 1 << (WIDTH - split)
    weights = np.empty(period)
    conditionals = np.empty((period, 1 << WIDTH))
    early = np.empty((period, length), dtype=complex)
    qft = Circuit(WIDTH).qft(list(range(WIDTH)), inverse=True)
    for k in range(period):
        # V|a^l> is column l in orbit-label coordinates; phi_k is column k.
        # Early powers wrap around the orbit when L exceeds r.
        c = fourier[:, k].conj() @ v_orbit[:, np.arange(length) % period]
        early[k] = c
        norm2 = float(np.vdot(c, c).real)
        weights[k] = norm2 / length
        if norm2 <= 1e-15:
            raise AssertionError("zero final-eigenphase weight is unsupported in this probe")
        control = np.zeros(1 << WIDTH, dtype=complex)
        alpha = c / np.sqrt(norm2)
        eigenvalue = np.exp(2j * np.pi * k / period)
        for l in range(length):
            for h in range(height):
                control[l + length * h] = alpha[l] * eigenvalue ** (length * h) / np.sqrt(height)
        transformed = statevec.run(qft, control)
        conditionals[k] = np.abs(transformed) ** 2
        if not np.isclose(np.sum(conditionals[k]), 1., atol=2e-12, rtol=0):
            raise AssertionError("conditional control state is not normalized")
    return v_eigen, weights, conditionals, early


def weighted_distribution(weights: np.ndarray, conditionals: np.ndarray) -> np.ndarray:
    return weights @ conditionals


def split_distribution(early: np.ndarray, period: int, split: int, k: int) -> np.ndarray:
    """Independent scalar late sum plus existing statevec early inverse-QFT."""
    length = 1 << split
    height = 1 << (WIDTH - split)
    qft_early = Circuit(split).qft(list(range(split)), inverse=True)
    alpha = early[k] / np.sqrt(np.vdot(early[k], early[k]).real)
    eigenvalue = np.exp(2j * np.pi * k / period)
    late = np.empty(height)
    for q in range(height):
        late[q] = abs(sum(eigenvalue ** (length * h)
                           * np.exp(-2j * np.pi * h * q / height)
                           for h in range(height)) / height) ** 2
    result = np.zeros(1 << WIDTH)
    for q in range(height):
        early_state = np.array([alpha[l] * np.exp(-2j * np.pi * q * l / (length * height))
                                for l in range(length)], dtype=complex)
        early_probs = np.abs(statevec.run(qft_early, early_state)) ** 2
        for z in range(length):
            result[q + height * z] = late[q] * early_probs[z]
    return result


def erased_distribution(weights: np.ndarray, early: np.ndarray, period: int, split: int) -> np.ndarray:
    """Same sampler after incorrectly replacing c[k,l] by |c[k,l]|."""
    length = 1 << split
    height = 1 << (WIDTH - split)
    qft = Circuit(WIDTH).qft(list(range(WIDTH)), inverse=True)
    result = np.zeros(1 << WIDTH)
    for k in range(period):
        norm2 = float(np.vdot(early[k], early[k]).real)
        alpha = np.abs(early[k]) / np.sqrt(norm2)
        control = np.zeros(1 << WIDTH, dtype=complex)
        eigenvalue = np.exp(2j * np.pi * k / period)
        for l in range(length):
            for h in range(height):
                control[l + length * h] = alpha[l] * eigenvalue ** (length * h) / np.sqrt(height)
        result += weights[k] * np.abs(statevec.run(qft, control)) ** 2
    return result


def sparse_rejection_audit(period: int, length: int, max_d: int, seed: int):
    """Parseval/envelope check for arbitrary normalized complex sparse columns."""
    rng = np.random.default_rng(seed)
    phases = np.exp(2j * np.pi * np.arange(period)[:, None]
                    * np.arange(period)[None, :] / period)
    details = []
    for d in range(1, max_d + 1):
        columns = np.zeros((length, period), dtype=complex)
        for l in range(length):
            support = rng.choice(period, size=d, replace=False)
            values = rng.normal(size=d) + 1j * rng.normal(size=d)
            values /= np.linalg.norm(values)
            columns[l, support] = values
        # h[k,l] = sum_j alpha[l,j] exp(+2*pi*i*k*j/r).
        h = phases @ columns.T
        acceptance = np.sum(np.abs(h) ** 2, axis=1) / (length * d)
        envelope_error = float(max(0., np.max(acceptance) - 1.))
        mean_error = float(abs(np.mean(acceptance) - 1. / d))
        details.append(dict(period=period, length=length, D=d,
                            max_acceptance=float(np.max(acceptance)),
                            mean_acceptance=float(np.mean(acceptance)),
                            envelope_error=envelope_error, mean_error=mean_error))
    return details


def sparse_oracle(matrix: np.ndarray, period: int, length: int):
    """Expose orbit-index columns, including periodic early powers."""
    indices = np.arange(length) % period

    def column(l):
        values = matrix[:, indices[l]]
        keep = np.flatnonzero(values != 0)  # no tolerance truncation
        return keep.astype(np.int64), values[keep]

    return column


def main() -> None:
    exp = Experiment("prefix_probe", doc=__doc__)
    exp.predict("P1", "weighted conditional eigenphase sampler matches spectral effects")
    exp.predict("P2", "N=7 physical |1>-|5> mixer obeys the weighted formula")
    exp.predict("P3", "conditional controls and y=q+H*z split factorization are correct")
    exp.predict("P4", "generic complex D-sparse rejection envelope and Parseval mean hold")
    exp.predict("P5", "SparseOrbitPrefix forced joints and phase weights match independent route")
    exp.must_fail("C1", "uniform final eigenphase weights are wrong when weights are nonuniform")
    exp.must_fail("C2", "erasing early conditional phases changes the output")

    rows = []
    start = time.perf_counter()
    cases = []
    for period in (3, 6):
        for split in (1, 2):
            cases.append((f"random_r{period}_s{split}", period, split,
                          seeded_unitary(period, 1200 + 100 * period + split), None))
    orbit = np.array([1, 3, 2, 6, 4, 5])
    cases.append(("N7a3_mixer_s2", 6, 2, physical_mixer(orbit, THETA), orbit))

    nonuniform_instances = 0
    phase_sensitive_instances = 0
    max_formula_error = 0.
    max_norm_error = 0.
    max_split_error = 0.
    uniform_errors = []
    erased_errors = []
    for name, period, split, v_orbit, physical_orbit in cases:
        v_eigen, weights, conditionals, early = prefix_data(v_orbit, period, split)
        exact = exact_outputs(v_eigen, period, split)
        predicted = weighted_distribution(weights, conditionals)
        uniform = np.mean(conditionals, axis=0)
        erased = erased_distribution(weights, early, period, split)
        split_errors = [float(np.max(np.abs(split_distribution(early, period, split, k)
                                     - conditionals[k]))) for k in range(period)]
        split_error = max(split_errors)
        formula_error = float(np.max(np.abs(predicted - exact)))
        normalization_error = float(abs(np.sum(weights) - 1.))
        uniform_error = float(np.max(np.abs(uniform - exact)))
        erased_error = float(np.max(np.abs(erased - exact)))
        max_terms = period if physical_orbit is None else 2
        prefix = SparseOrbitPrefix(period, split, sparse_oracle(v_orbit, period, 1 << split), max_terms)
        prefix_phase_error = float(max(abs(prefix.phase_probability(k) - weights[k])
                                      for k in range(period)))
        prefix_marginal = np.zeros(1 << WIDTH)
        for k in range(period):
            for y in range(1 << WIDTH):
                forced = prefix.forced_joint(WIDTH, k, y)
                prefix_marginal[y] += forced["joint_latent_output_probability"]
        prefix_joint_error = float(np.max(np.abs(prefix_marginal - exact)))
        max_formula_error = max(max_formula_error, formula_error)
        max_norm_error = max(max_norm_error, normalization_error)
        max_split_error = max(max_split_error, split_error)
        uniform_errors.append(uniform_error)
        erased_errors.append(erased_error)
        nonuniform = float(np.max(weights) - np.min(weights)) > 1e-10
        phase_sensitive = erased_error > 1e-10
        nonuniform_instances += int(nonuniform)
        phase_sensitive_instances += int(phase_sensitive)
        unitary_error = float(np.max(np.abs(v_orbit.conj().T @ v_orbit-np.eye(period))))
        exp.check("P1" if physical_orbit is None else "P2", formula_error < 2e-10
                  and unitary_error < 1e-12,
                  f"{name}: max weighted-vs-effects error={formula_error:.2e}")
        exp.check("P3", normalization_error < 2e-12 and split_error < 2e-10,
                  f"{name}: weight norm={normalization_error:.2e}, split={split_error:.2e}")
        exp.check("P5", prefix_phase_error < 2e-12 and prefix_joint_error < 2e-10,
                  f"{name}: prefix phase={prefix_phase_error:.2e}, forced marginal={prefix_joint_error:.2e}")
        rows.append(dict(series=name, period=period, split=split, width=WIDTH,
                         weights=weights.tolist(), exact=exact.tolist(),
                         weighted=predicted.tolist(), uniform=uniform.tolist(),
                         erased=erased.tolist(), formula_error=formula_error,
                         normalization_error=normalization_error,
                         split_factorization_error=split_error,
                         prefix_phase_error=prefix_phase_error,
                         prefix_forced_marginal_error=prefix_joint_error,
                         prefix_stats=prefix.stats(),
                         uniform_error=uniform_error, erased_phase_error=erased_error,
                         unitary_error=unitary_error,
                         mixer_theta_over_pi=None if physical_orbit is None else .5))
        exp.log(f"{name}: w range=[{weights.min():.6g},{weights.max():.6g}], "
                f"uniform error={uniform_error:.3g}, phase-erased error={erased_error:.3g}")

    # These controls are evaluated only after checking that they are not vacuous.
    exp.fail_check("C1", nonuniform_instances > 0 and max(uniform_errors) > 1e-8,
                   f"nonuniform cases={nonuniform_instances}, largest uniform-k error={max(uniform_errors):.6g}")
    exp.fail_check("C2", phase_sensitive_instances > 0 and max(erased_errors) > 1e-8,
                   f"phase-sensitive cases={phase_sensitive_instances}, largest erased-phase error={max(erased_errors):.6g}")
    sparse_rows = sparse_rejection_audit(6, 4, 3, 917)
    sparse_ok = all(row["envelope_error"] < 2e-12 and row["mean_error"] < 2e-12
                    for row in sparse_rows)
    exp.check("P4", sparse_ok, f"sparse rejection rows={sparse_rows}")
    elapsed = time.perf_counter() - start
    exp.finish(report_path="out/prefix_probe.json", rows=rows,
               metadata=dict(width=WIDTH, theta_over_pi=.5, cases=len(cases),
                             max_formula_error=max_formula_error,
                             max_weight_normalization_error=max_norm_error,
                             max_split_factorization_error=max_split_error,
                             sparse_rejection_audit=sparse_rows,
                             setup="dense orbit Fourier and single_defect_effects are validation references",
                             independent_route="control-only statevec inverse QFT; no work propagator",
                             rejection_audit="w_k <= D/r if each V|a^l> has at most D known orbit-index amplitudes; locating indices may hide discrete-log cost",
                             dense_payload_cap_bytes=1_000_000 * 16,
                             elapsed_seconds=elapsed, numpy=np.__version__))


if __name__ == "__main__":
    main()
