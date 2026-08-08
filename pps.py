"""Pauli Path Simulation with coefficient-threshold truncation.

Heisenberg picture: propagate the observable backwards through the circuit,
expanding it in the Pauli basis and discarding terms with |c| < delta.
Follows Gharibyan et al. (arXiv:2507.10771) Eqs. 8-9.
"""
from __future__ import annotations
import numpy as np
from pauli import commutes, i_sigma_p


class PPSResult:
    def __init__(self):
        self.n_terms: list[int] = []      # unique Pauli terms after each gate
        self.norm: list[float] = []       # sqrt(sum c^2) after each gate
        self.final_coeffs: np.ndarray | None = None
        self.expectation: float = 0.0
        self.n_max: int = 0
        self.truncated_weight: float = 0.0
        self.hit_cap: bool = False
        self.final_terms: dict[tuple[int, int], float] = {}


def propagate(circuit, observable, delta=0.0, max_terms=2_000_000,
              snapshot_at=None) -> PPSResult:
    """Evolve `observable` backwards through `circuit`.

    observable : dict {(x, z): coeff}
    delta      : coefficient truncation threshold (0 = exact)
    Returns a PPSResult; `snapshot_at` optionally records the coefficient
    multiset after that many gates.
    """
    terms = dict(observable)
    res = PPSResult()
    snapshot = None

    for gi, (sigma, theta) in enumerate(reversed(circuit.gates)):
        c, s = np.cos(theta), np.sin(theta)
        new: dict[tuple[int, int], float] = {}

        # U^dag P U = P                          if [P, sigma] = 0
        #           = cos(t) P + sin(t) (i sigma P)   if {P, sigma} = 0
        #
        # The only true no-op is theta = 0 mod 2pi. At theta = pi mod 2pi we
        # have sin = 0 but cos = -1, so anticommuting terms must be NEGATED --
        # that is exactly what X, Y and Z gates do, and skipping them here
        # silently turned every Pauli gate into the identity.
        if abs(s) < 1e-12 and c > 0:            # theta = 0 mod 2pi: identity
            new = terms
        else:
            for P, coeff in terms.items():
                if commutes(P, sigma):
                    new[P] = new.get(P, 0.0) + coeff
                    continue
                if abs(c) > 1e-12:
                    new[P] = new.get(P, 0.0) + c * coeff
                if abs(s) > 1e-12:
                    Q, sign = i_sigma_p(sigma, P)
                    new[Q] = new.get(Q, 0.0) + s * sign * coeff

        if delta > 0:
            kept = {}
            dropped = 0.0
            for P, v in new.items():
                if abs(v) >= delta:
                    kept[P] = v
                else:
                    dropped += v * v
            res.truncated_weight += dropped
            new = kept
        else:
            # floor at fp noise: exact cancellations leave ~1e-16 residue that
            # would otherwise be counted as real terms
            new = {P: v for P, v in new.items() if abs(v) > 1e-13}

        terms = new
        res.n_terms.append(len(terms))
        res.norm.append(float(np.sqrt(sum(v * v for v in terms.values()))))
        if snapshot_at is not None and gi + 1 == snapshot_at:
            snapshot = np.array([abs(v) for v in terms.values()])

        if len(terms) > max_terms:
            res.hit_cap = True
            break

    res.final_terms = terms
    res.n_max = max(res.n_terms) if res.n_terms else 0
    res.final_coeffs = snapshot if snapshot is not None else np.array(
        [abs(v) for v in terms.values()])
    # <0|P|0> is 1 for pure-Z strings, 0 otherwise
    res.expectation = float(sum(v for (x, _), v in terms.items() if x == 0))
    return res


def exact_expectation(circuit, observable) -> float:
    """Dense <0|U^dag O U|0> for cross-checking. Exponential in n."""
    from pauli import to_matrix
    n = circuit.n
    U = circuit.to_unitary()
    O = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for P, c in observable.items():
        O += c * to_matrix(P, n)
    psi = np.zeros(2 ** n, dtype=complex); psi[0] = 1
    out = U @ psi
    return float(np.real(out.conj() @ O @ out))


def fit_power_law(coeffs: np.ndarray, delta: float):
    """Fit rho(t) ~ A / t^(m+1) for |t| > delta via MLE on the Pareto tail.

    Returns (m, n_samples). MLE for a Pareto with exponent alpha = m+1:
        alpha_hat = 1 + N / sum(log(t_i / delta))
    """
    t = coeffs[coeffs >= delta]
    t = t[t > 0]
    if len(t) < 50:
        return None, len(t)
    alpha = 1.0 + len(t) / np.sum(np.log(t / delta))
    return alpha - 1.0, len(t)
