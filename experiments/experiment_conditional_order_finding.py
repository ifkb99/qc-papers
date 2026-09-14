"""Actual order-finding distributions via conditional branches (TODO 15c).

DERIVATION BEFORE MEASUREMENT. A measured |+> exponent bit induces
(P0+(-1)^b exp(i theta)P1)/2 on the work state. Inverse-QFT feedback is
theta=-pi*prefix/2^step; largest input exponent power is processed first.
Keeping the outcome preserves the cross terms lost by averaging the control.

P1: entire small distributions equal independent FFT and existing Fourier-
arithmetic state-vector circuits. P2: individual wide-output probabilities
equal the independent geometric-series formula. P3: dense work storage is
independent of exponent width, but remains exponential in work qubits.
C1: omitting feedback must fail for an odd-order example. C2: reversing output
bits must fail. C3: replacing the instrument by uninformative fair bits fails.
C4: noncommuting controlled work gates must defeat naive order reversal.

The descending-power schedule uses commutation of modular multiplications on
the valid, clean reachable subspace. No commutation of the arbitrary full-
scratch maps is assumed. C4 was added on auditing this implementation assumption.

This is a known semiclassical-QFT route, not a new factoring algorithm. The
sampler never receives the period or an enumerated orbit. Actual branch maps
are extracted by full work-space replay of existing arithmetic circuits;
this setup is explicitly exponential in work width. Small full distributions
enumerate all outcomes solely for validation. Wide runs sample paths only.

Run: OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_conditional_order_finding
"""
from __future__ import annotations
import time
import numpy as np
from lab import Experiment
from circuits import Circuit
from lab.semiclassical import distribution, sample, validate_inputs, order_finding_probability
from toffoli_arith import ToffoliModExp
from modexp import ModExp
import accel
import walsh
import statevec


def arithmetic_pairs(N, a, width):
    me = ToffoliModExp(N=N, a=a, n_exp=1)
    dim = 1 << me.exp[0]
    pairs, cache = [], {}
    start = time.perf_counter()
    for i in range(width):
        multiplier = pow(a, 1 << i, N)
        if multiplier not in cache:
            qc = me.u_a(me.exp[0], multiplier)
            perm = (accel.classical_permutation(qc) if accel.enabled()
                    else walsh.classical_permutation(qc))
            ids = np.arange(2*dim)
            assert np.array_equal(perm // dim, ids // dim), "control changed"
            pair = (perm[:dim], perm[dim:] - dim)
            # Gate all legitimate work inputs, including scratch cleanup.
            clean = np.arange(N) << me.x[0]
            assert np.array_equal(pair[0][clean], clean)
            assert np.array_equal(pair[1][clean], ((np.arange(N)*multiplier) % N) << me.x[0])
            cache[multiplier] = pair
        pairs.append(cache[multiplier])
    psi = np.zeros(dim, complex)
    psi[1 << me.x[0]] = 1.
    validate_inputs(pairs, psi)
    return me, pairs, psi, time.perf_counter()-start


def fft_reference(N, a, width):
    Q = 1 << width
    amplitudes = np.zeros((Q, N), complex)
    amplitudes[np.arange(Q), [pow(a, e, N) for e in range(Q)]] = 1 / np.sqrt(Q)
    transformed = np.fft.fft(amplitudes, axis=0) / np.sqrt(Q)
    return np.sum(np.abs(transformed)**2, axis=1)


def circuit_reference(N, a, width):
    me = ModExp(N=N, a=a, n_exp=width)
    psi = statevec.run(me.build_shor())
    return np.sum(np.abs(psi.reshape(1 << width, -1))**2, axis=1)


def main():
    exp = Experiment("conditional_order_finding", doc=__doc__)
    exp.predict("P1", "full small distributions match FFT and independent compiled circuits")
    exp.predict("P2", "wide sampled path probabilities match the geometric-series reference")
    exp.predict("P3", "work storage does not grow with exponent width")
    exp.must_fail("C1", "omitting phase feedback fails for odd-order arithmetic")
    exp.must_fail("C2", "wrong output bit order fails")
    exp.must_fail("C3", "discarding interference and returning fair bits fails")
    exp.must_fail("C4", "noncommuting gates defeat the coherent-to-iterative reordering")
    rows = []
    ids = np.arange(4)
    v0, v1 = ids ^ ((ids & 1) << 1), ids ^ ((ids >> 1) & 1)
    initial = np.array([0., 1., 0., 0.], complex)
    iterative, _ = distribution([(ids, v0), (ids, v1)], initial)
    control = Circuit(4).x(0).h(2).h(3).toffoli(2, 0, 1).toffoli(3, 1, 0)
    control.qft([2, 3], inverse=True)
    coherent = np.sum(np.abs(statevec.run(control).reshape(4, 4))**2, axis=1)
    exp.fail_check("C4", np.max(np.abs(iterative-coherent)) > .1,
                   f"noncommuting reordering error={np.max(np.abs(iterative-coherent)):.6g}")
    for N, a in [(7, 6), (7, 3), (15, 7)]:
        me, all_pairs, psi, setup = arithmetic_pairs(N, a, 32)
        period = me.order()  # only validation/metadata; not passed to sampler
        for width in (3, 4, 6, 8):
            start = time.perf_counter()
            probabilities, profile = distribution(all_pairs[:width], psi)
            seconds = time.perf_counter()-start
            reference = fft_reference(N, a, width)
            error = float(np.max(np.abs(probabilities-reference)))
            formula = np.array([order_finding_probability(y, width, period)
                                for y in range(1 << width)])
            formula_error = float(np.max(np.abs(probabilities-formula)))
            exp.check("P1", error < 1e-10 and formula_error < 1e-10
                      and abs(probabilities.sum()-1) < 1e-10,
                      f"N={N} a={a} t={width}: FFT={error:.2e}, formula={formula_error:.2e}")
            if width == 3:
                compiled = circuit_reference(N, a, width)
                compiled_error = float(np.max(np.abs(probabilities-compiled)))
                exp.check("P1", compiled_error < 1e-9,
                          f"independent Fourier compilation: error={compiled_error:.2e}")
            if (N, a, width) == (7, 3, 4):
                wrong, _ = distribution(all_pairs[:width], psi, feedback=False)
                exp.fail_check("C1", np.max(np.abs(wrong-reference)) > 1e-3,
                               f"missing-feedback error={np.max(np.abs(wrong-reference)):.6g}")
            if (N, a, width) == (7, 6, 4):
                reverse = [int(f"{y:0{width}b}"[::-1], 2) for y in range(1 << width)]
                exp.fail_check("C2", np.max(np.abs(probabilities[reverse]-reference)) > .1,
                               "output endianness is not optional")
                exp.fail_check("C3", np.max(np.abs(reference-1/(1 << width))) > .1,
                               "averaged work channel cannot replace the measured instrument")
            rows.append(dict(kind="distribution", N=N, a=a, period=period, width=width,
                             max_error=error, formula_error=formula_error, seconds=seconds,
                             profile=profile, probabilities=probabilities.tolist(),
                             setup_seconds=setup, stored_amplitudes=psi.size))
        rng = np.random.default_rng(20260909)
        for width in (4, 8, 16, 32):
            start = time.perf_counter()
            paths = [sample(all_pairs[:width], psi, rng, validated=True) for _ in range(16)]
            elapsed = time.perf_counter()-start
            errors = [abs(p["path_probability"]-order_finding_probability(p["output"], width, period))
                      for p in paths]
            exp.check("P2", max(errors) < 1e-9,
                      f"N={N} a={a} t={width}: sampled probability error={max(errors):.2e}")
            exp.check("P3", all(p["stored_amplitudes"] == psi.size
                                and p["peak_active_amplitudes"] <= period for p in paths),
                      f"allocated={psi.size}, occupied peak={max(p['peak_active_amplitudes'] for p in paths)}, period={period}")
            rows.append(dict(kind="samples", N=N, a=a, period=period, width=width,
                             max_error=max(errors), seconds=elapsed, paths=paths,
                             setup_seconds=setup, stored_amplitudes=psi.size))
    exp.finish(report_path="out/conditional_order_finding.json", rows=rows,
               metadata=dict(numpy=np.__version__, gpu=accel.enabled(),
                             sampler="dense work vector; no exponent vector, no period input"))


if __name__ == "__main__":
    main()
