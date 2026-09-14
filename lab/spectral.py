"""Small, order-INFORMED spectral references for one intervening work unitary.

These are finite Fourier sums/effects, not a new circuit propagator or scalable
sampler. Preserve ascending arithmetic order: U^low, V, U^high, inverse QFT.
All arrays are dense and explicitly capped. The work unitary must preserve the
chosen cyclic orbit. No claim is made about arbitrary dirty scratch states.
"""
from __future__ import annotations
import numpy as np


def split_phase_filters(period, width, split, *, dtype=np.float64,
                        max_entries=1_000_000):
    """Return Q-by-r early/late normalized sums for all Fourier outputs.

    L=2^split, H=Q/L; eigenvalues are exp(2*pi*i*k/r).
    A[y,k] = mean_l exp(2*pi*i*(k/r-y/Q)*l), 0<=l<L.
    B[y,k] = mean_h exp(2*pi*i*(k/r-y/Q)*L*h), 0<=h<H.
    The cap budgets the eventual Q*r*r effects as well as these filters.
    """
    if (not all(isinstance(v, (int, np.integer)) for v in (period, width, split))
            or period < 1 or not 0 <= split <= width or width > 12):
        raise ValueError("invalid period, width or split (width cap is 12)")
    period, width, split = int(period), int(width), int(split)
    Q, L = 1 << width, 1 << split
    if Q * period * period > max_entries:
        raise ValueError("dense spectral reference exceeds its entry budget")
    dtype = np.dtype(dtype)
    if dtype.kind != "f" or dtype.itemsize < 8:
        raise ValueError("use float64 or extended floating precision")
    real = dtype.type
    pi = np.arccos(real(-1))
    frequencies = (np.arange(period, dtype=dtype)[None, :] / period
                   - np.arange(Q, dtype=dtype)[:, None] / Q)
    # Accumulate 2-D arrays: do not allocate a Q*r*L intermediate.
    def average(length, stride):
        result = np.zeros(frequencies.shape, dtype=np.result_type(dtype, 1j))
        for j in range(length):
            result += np.exp(2j * pi * frequencies * (stride*j))
        return result / length
    return average(L, 1), average(Q // L, L)


def single_defect_effects(period, width, split, defect, *, dtype=np.float64,
                          max_entries=1_000_000):
    """Return E_y=K_y^dagger K_y in the U eigenbasis.

    K_y = diag(B_y) V diag(A_y). For an arbitrary orbit density matrix rho,
    p(y)=trace(E_y rho). With initial |1>, rho[k,k']=1/r. Replacing rho by
    I/r is the initial-eigenphase-dephasing baseline, NOT the same circuit.
    Final eigenbasis coherences never contribute to the traced-work output.
    """
    A, B = split_phase_filters(period, width, split, dtype=dtype,
                               max_entries=max_entries)
    V = np.asarray(defect, dtype=A.dtype)
    if (V.shape != (period, period) or not np.all(np.isfinite(V))
            or np.max(np.abs(V.conj().T @ V - np.eye(period))) > 1e-10):
        raise ValueError("defect must be a finite unitary on the chosen orbit")
    effects = np.empty((1 << width, period, period), dtype=A.dtype)
    for y, (early, late) in enumerate(zip(A, B)):
        K = late[:, None] * V * early[None, :]
        effects[y] = K.conj().T @ K
    return effects


def coherence_response(effects, *, tolerance=1e-10):
    """Numerical response rank for off-diagonal Hermitian input coordinates.

    Coordinates are Re(rho[k,l]), Im(rho[k,l]) for k<l. Rank counts independent
    real linear statistics sufficient for THIS fixed measurement and arbitrary
    inputs. It is not a count of required stored matrix entries or a sampling
    memory lower bound. SVD diagnostics use float64, even for extended inputs.
    """
    effects = np.asarray(effects)
    if (effects.ndim != 3 or min(effects.shape) < 1
            or effects.shape[1] != effects.shape[2]
            or not np.all(np.isfinite(effects)) or tolerance <= 0
            or not np.isfinite(tolerance)):
        raise ValueError("invalid effects or tolerance")
    pairs = np.triu_indices(effects.shape[1], 1)
    upper = effects[:, pairs[0], pairs[1]]
    response = np.concatenate((2*upper.real, 2*upper.imag), axis=1).astype(float)
    singular = np.linalg.svd(response, compute_uv=False)
    return dict(rank=int(np.count_nonzero(singular > tolerance)),
                tolerance=tolerance, singular_values=singular.tolist(),
                detectable_pairs=int(np.count_nonzero(
                    np.max(np.abs(upper), axis=0) > tolerance)) if upper.size else 0)
