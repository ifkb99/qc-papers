"""Bounded tests of interval-prefix Fourier probabilities and progression reduction.

DERIVED BEFORE MEASUREMENT.  For a normalized interval state

    |psi> = n^(-1/2) sum_(m=0)^(n-1) exp(-2*pi*i*eta*m) |m>,  n <= T=2^b,

the negative-sign Fourier output has the following LOW-j-bit marginal.  Put
d=T/2^j, q=floor(n/d), and e=n mod d.  For a in [0,2^j),

  P(a) = ((d-e)|G_q(a/2^j + eta*d)|^2
          + e|G_(q+1)(a/2^j + eta*d)|^2)/(n*2^j),

where G_h(x)=sum_(m=0)^(h-1) exp(-2*pi*i*x*m).  This experiment compares
that finite formula with a direct numpy negative-sign FFT and with the
existing control-only inverse-QFT state-vector circuit at b<=8.

For a progression l=p+r*m in [0,L), feedback exp(-2*pi*i*q_out*l/Q), let
g=gcd(r,L), T=L/g and beta=r/g.  Since beta is odd for T>1,
the full L-output distribution should be the interval distribution with
eta=r*q_out/Q, evaluated at w=beta*z mod T and divided by g.  The offset p
is a global phase and therefore cannot change probabilities, although the
coherent sum over m must be retained.

P1: Every tested low-bit prefix formula agrees with direct FFT, including
    n=1, n=T, non-powers, signed/zero phases, and near resonances.
P2: Existing inverse-QFT state-vector probabilities agree with the FFT route
    on all b<=8 interval rows selected below.
P3: Progression reduction agrees with direct full distributions for several
    powers of two, odd/even strides, gcd lifts, offsets, feedback values, and
    the stride>L singleton case; p is probability-invariant.
C1: Wrong high-bit/endian marginal must fail against the LOW-bit formula.
C2: Ignoring gcd/lift (using an L-dimensional interval) must fail.
C3: Using the wrong progression count must fail.
C4: Dropping inverse-QFT feedback must fail.
C5: Replacing coherent progression interference by an incoherent uniform
    distribution must fail.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_progression_prefix
"""
from __future__ import annotations

import math
import time
from fractions import Fraction

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.fourier_sampling import interval_prefix_probability, interval_path
import statevec


MAX_VECTOR_ENTRIES = 1 << 8
MAX_REFERENCE_BYTES = 1 << 20
TOL = 3e-11


def _guard_entries(size: int, *, itemsize: int = 16, label: str) -> None:
    if size < 0 or size > MAX_VECTOR_ENTRIES:
        raise ValueError(f"{label} allocation {size} exceeds {MAX_VECTOR_ENTRIES}")
    if size * itemsize > MAX_REFERENCE_BYTES:
        raise ValueError(f"{label} allocation exceeds {MAX_REFERENCE_BYTES} bytes")


def geometric_sum(h: int, x: float) -> complex:
    """Finite sum, deliberately direct so the reference has no singular branch."""
    if h <= 0:
        return 0j
    _guard_entries(h, label="direct geometric sum")
    m = np.arange(h, dtype=float)
    return complex(np.sum(np.exp(-2j * np.pi * x * m)))


def interval_state(T: int, n: int, eta: float) -> np.ndarray:
    if T <= 0 or T & (T - 1) or not (1 <= n <= T):
        raise ValueError("interval requires power-of-two T and 1 <= n <= T")
    _guard_entries(T, label="interval state")
    result = np.zeros(T, dtype=complex)
    m = np.arange(n, dtype=float)
    result[:n] = np.exp(-2j * np.pi * eta * m) / np.sqrt(n)
    return result


def fft_distribution(T: int, n: int, eta: float) -> np.ndarray:
    state = interval_state(T, n, eta)
    return np.abs(np.fft.fft(state) / np.sqrt(T)) ** 2


def statevec_distribution(T: int, n: int, eta: float) -> np.ndarray:
    """Existing inverse-QFT engine; compare probabilities, not global phase."""
    b = T.bit_length() - 1
    state = interval_state(T, n, eta)
    qft = Circuit(b).qft(list(range(b)), inverse=True)
    return np.abs(statevec.run(qft, state)) ** 2


def prefix_formula(T: int, n: int, eta: float, j: int) -> np.ndarray:
    if not (0 <= j <= T.bit_length() - 1):
        raise ValueError("prefix width out of range")
    count = 1 << j
    _guard_entries(count, label="prefix formula")
    d = T // count
    q, e = divmod(n, d)
    result = np.empty(count, dtype=float)
    for a in range(count):
        x = a / count + eta * d
        result[a] = ((d - e) * abs(geometric_sum(q, x)) ** 2
                     + (e * abs(geometric_sum(q + 1, x)) ** 2 if e else 0.)) / (n * count)
    return result


def low_prefix(p: np.ndarray, j: int) -> np.ndarray:
    count = 1 << j
    return np.array([np.sum(p[a::count]) for a in range(count)])


def high_prefix(p: np.ndarray, j: int) -> np.ndarray:
    count = 1 << j
    block = len(p) // count
    return np.array([np.sum(p[a * block:(a + 1) * block]) for a in range(count)])


def progression_distribution(L: int, p: int, r: int, n: int,
                             q_out: int, Q: int) -> np.ndarray:
    if L <= 0 or L & (L - 1) or n < 1 or p < 0:
        raise ValueError("bad progression parameters")
    _guard_entries(L, label="progression state")
    _guard_entries(n, label="progression support")
    support = p + r * np.arange(n, dtype=np.int64)
    if np.any(support < 0) or np.any(support >= L):
        raise ValueError("progression is not contained in [0,L)")
    _guard_entries(L, label="progression state")
    state = np.zeros(L, dtype=complex)
    state[support] = np.exp(-2j * np.pi * q_out * support / Q) / np.sqrt(n)
    return np.abs(np.fft.fft(state) / np.sqrt(L)) ** 2


def reduced_progression(L: int, r: int, n: int, q_out: int, Q: int) -> np.ndarray:
    _guard_entries(L, label="reduced progression")
    g = math.gcd(r, L)
    T = L // g
    beta = r // g
    eta = r * q_out / Q
    reduced = fft_distribution(T, n, eta)
    result = np.empty(L, dtype=float)
    for z in range(L):
        result[z] = reduced[(beta * z) % T] / g
    return result


def wrong_count_progression(L: int, r: int, n: int, q_out: int, Q: int) -> np.ndarray:
    """Control: the same reduction with one extra term in the interval."""
    g = math.gcd(r, L)
    T = L // g
    beta = r // g
    reduced = fft_distribution(T, n + 1, r * q_out / Q)
    return np.array([reduced[(beta * z) % T] / g for z in range(L)])


def main() -> None:
    exp = Experiment("progression_prefix", doc=__doc__)
    # All predictions and controls are declared before any numerical route runs.
    exp.predict("P1", "finite geometric-sum low-bit marginals equal direct FFT")
    exp.predict("P2", "existing inverse-QFT state-vector probabilities equal FFT")
    exp.predict("P3", "gcd/odd-stride progression reduction equals direct full distributions")
    exp.predict("P4", "production scalar marginals and every forced bit path agree with FFT")
    exp.must_fail("C1", "wrong high-bit/endian marginal disagrees with LOW-bit formula")
    exp.must_fail("C2", "ignoring gcd/lift disagrees for a nontrivial gcd")
    exp.must_fail("C3", "wrong progression length changes the distribution")
    exp.must_fail("C4", "dropping inverse-QFT feedback changes the distribution")
    exp.must_fail("C5", "dropping coherent interference is not a valid progression distribution")

    rows = []
    started = time.perf_counter()

    # Interval formula sweep: T is the only changing size parameter inside a
    # row family; eta includes signed, zero, rational, and near-resonant cases.
    exp.section("P1/P2 interval prefixes")
    max_fft_error = max_statevec_error = max_prefix_error = 0.0
    max_production_error = 0.
    interval_rows = 0
    for b in range(1, 9):
        T = 1 << b
        n_values = sorted({1, T, 2, 3, 5, T - 1} & set(range(1, T + 1)))
        near = Fraction(1, T) + Fraction(1, 10 * T * T)
        eta_values = (Fraction(0), Fraction(1, 7), Fraction(-2, 9), near, -near)
        for n in n_values:
            for eta_fraction in eta_values:
                eta = float(eta_fraction)
                p_fft = fft_distribution(T, n, eta)
                p_sv = statevec_distribution(T, n, eta)
                fft_error = float(abs(np.sum(p_fft) - 1.0))
                sv_error = float(np.max(np.abs(p_sv - p_fft)))
                prefix_errors = []
                production_error = 0.
                for j in range(b + 1):
                    measured = low_prefix(p_fft, j)
                    predicted = prefix_formula(T, n, eta, j)
                    prefix_errors.append(float(np.max(np.abs(measured - predicted))))
                    production = np.array([interval_prefix_probability(
                        b, n, eta_fraction.numerator, eta_fraction.denominator, j, a)
                        for a in range(1 << j)])
                    production_error = max(production_error, float(np.max(np.abs(measured-production))))
                    exp.check("P1", prefix_errors[-1] < TOL,
                              f"T={T}, n={n}, eta={eta_fraction}, j={j}: "
                              f"prefix error={prefix_errors[-1]:.2e}")
                exp.check("P2", sv_error < TOL,
                          f"T={T}, n={n}, eta={eta_fraction}: statevec-vs-FFT={sv_error:.2e}")
                paths = np.array([interval_path(b, n, eta_fraction.numerator,
                                               eta_fraction.denominator, output=y)["path_probability"]
                                  for y in range(T)])
                production_error = max(production_error, float(np.max(np.abs(paths-p_fft))))
                exp.check("P4", production_error < TOL,
                          f"T={T}, n={n}, eta={eta_fraction}: production={production_error:.2e}")
                max_production_error = max(max_production_error, production_error)
                max_fft_error = max(max_fft_error, fft_error)
                max_statevec_error = max(max_statevec_error, sv_error)
                max_prefix_error = max(max_prefix_error, max(prefix_errors))
                interval_rows += 1
                rows.append(dict(kind="interval", T=T, b=b, n=n,
                                 eta=str(eta_fraction), fft_normalization_error=fft_error,
                                 statevec_error=sv_error, max_prefix_error=max(prefix_errors)))

    # Progression family.  Include every 2-adic gcd at L=32 and selected odd
    # residues; larger rows remain tiny (only the direct Fourier output array).
    exp.section("P3 progression reduction")
    progression_rows = 0
    max_reduction_error = max_offset_error = 0.0
    for L in (4, 8, 16, 32, 64, 128, 256):
        strides = sorted({1, 2, 3, 4, 6, 8, 16, 32, L, L + 1} & set(range(1, L + 2)))
        for r in strides:
            max_n = 1 if r > L else min(4, (L - 1) // r + 1)
            n_values = sorted({1, max_n, min(2, max_n), min(3, max_n)})
            for n in n_values:
                max_p = L - 1 - r * (n - 1) if r <= L else 0
                if max_p < 0:
                    continue
                offsets = sorted({0, min(1, max_p), max_p})
                for q_out in (0, 1, 3):
                    Q = 4 * L
                    prediction = reduced_progression(L, r, n, q_out, Q)
                    for p in offsets:
                        direct = progression_distribution(L, p, r, n, q_out, Q)
                        reduction_error = float(np.max(np.abs(direct - prediction)))
                        max_reduction_error = max(max_reduction_error, reduction_error)
                        exp.check("P3", reduction_error < TOL,
                                  f"L={L}, r={r}, n={n}, p={p}, q={q_out}: "
                                  f"reduction error={reduction_error:.2e}")
                        if p != 0:
                            offset_error = float(np.max(np.abs(direct - progression_distribution(
                                L, 0, r, n, q_out, Q))))
                            max_offset_error = max(max_offset_error, offset_error)
                        progression_rows += 1
                        rows.append(dict(kind="progression", L=L, r=r, gcd=math.gcd(r, L),
                                         T=L // math.gcd(r, L), beta=r // math.gcd(r, L),
                                         n=n, p=p, q_out=q_out, Q=Q,
                                         max_error=reduction_error))

    # Must-fail controls use fixed, non-vacuous rows with visible interference.
    exp.section("must-fail controls")
    T, n, eta, j = 16, 5, 1 / 7, 2
    p = fft_distribution(T, n, eta)
    endian_error = float(np.max(np.abs(high_prefix(p, j) - prefix_formula(T, n, eta, j))))
    exp.fail_check("C1", endian_error > 1e-5, f"wrong-endian error={endian_error:.6g}")

    L, r, n, q_out, Q = 16, 4, 3, 1, 4 * 16
    direct = progression_distribution(L, 0, r, n, q_out, Q)
    ignored_gcd = fft_distribution(L, n, r * q_out / Q)
    gcd_error = float(np.max(np.abs(direct - ignored_gcd)))
    exp.fail_check("C2", gcd_error > 1e-5, f"ignored-gcd/lift error={gcd_error:.6g}")

    L, r, n, q_out, Q = 16, 1, 3, 1, 4 * 16
    direct = progression_distribution(L, 0, r, n, q_out, Q)
    count_error = float(np.max(np.abs(direct - wrong_count_progression(L, r, n, q_out, Q))))
    exp.fail_check("C3", count_error > 1e-5, f"wrong-count error={count_error:.6g}")

    L, r, n, q_out, Q = 16, 3, 3, 1, 4 * 16
    direct = progression_distribution(L, 0, r, n, q_out, Q)
    no_feedback = progression_distribution(L, 0, r, n, 0, Q)
    feedback_error = float(np.max(np.abs(direct - no_feedback)))
    exp.fail_check("C4", feedback_error > 1e-5, f"dropped-feedback error={feedback_error:.6g}")

    L, r, n, Q = 16, 1, 3, 4 * 16
    coherent = progression_distribution(L, 0, r, n, 0, Q)
    incoherent = np.full(L, 1 / L)
    interference_error = float(np.max(np.abs(coherent - incoherent)))
    exp.fail_check("C5", interference_error > 1e-5,
                   f"deleted-interference error={interference_error:.6g}")

    exp.finish(report_path="out/progression_prefix_production_v2.json", rows=rows,
               metadata=dict(interval_rows=interval_rows, progression_rows=progression_rows,
                             max_fft_normalization_error=max_fft_error,
                             max_statevec_error=max_statevec_error,
                             max_prefix_error=max_prefix_error,
                             max_production_error=max_production_error,
                             max_reduction_error=max_reduction_error,
                             max_offset_error=max_offset_error,
                             controls=dict(wrong_endian=endian_error,
                                           ignored_gcd=gcd_error, wrong_count=count_error,
                                           dropped_feedback=feedback_error,
                                           deleted_interference=interference_error),
                             allocation_cap_entries=MAX_VECTOR_ENTRIES,
                             numpy=np.__version__, elapsed_seconds=time.perf_counter() - started,
                             interpretation="finite probability formula and bounded verifier; not a generic sampler"))


if __name__ == "__main__":
    main()
