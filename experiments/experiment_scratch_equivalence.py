"""Same logical isometry, different scratch extensions (TODO 15d).

DERIVATION BEFORE MEASUREMENT. Append W_k = product_{j<k}
CCX(b_j,t_j,x_0) after modular exponentiation. On the clean-output code all
b and t bits are zero, hence every W_k is identity, including on coherent
superpositions. The same-layout full-space permutations can nevertheless
have different Walsh support. The exact valid-subspace pullback and actual
order-finding output distribution must be unchanged.

P1: all valid e,x<N images and clean-subspace observables agree.
P2: coherent gate-level guards preserve the code up to one global phase.
P3: the entire order-finding output distribution agrees after inverse QFT.
H1: at least one guard changes full Walsh support (measured, not guaranteed).
C1: an exponent-controlled output flip must change legitimate arithmetic.
C2: for N=7,a=3 it must also change order-finding output statistics.
P4: that invalid CNOT nevertheless preserves full Walsh sparsity (it translates
the frequency mask by the unchanged exponent bit). This tests why counts alone
cannot certify physical equivalence.

Uses existing permutation replay, state-vector gate evolution and inverse QFT.
No layout, modulus, base or observable changes within a guard-count series.
Run: OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_scratch_equivalence
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from lab import Experiment
from toffoli_arith import ToffoliModExp
import accel
import walsh
import statevec


def replay(qc):
    return accel.classical_permutation(qc) if accel.enabled() else walsh.classical_permutation(qc)


def main():
    exp = Experiment("scratch_equivalence", doc=__doc__)
    exp.predict("P1", "same logical isometry on all valid e,x<N; reduced observables identical")
    exp.predict("P2", "guard gates preserve coherent code states up to a global phase")
    exp.predict("P3", "entire post-inverse-QFT distribution unchanged")
    exp.predict("H1", "measure a full-spectrum difference despite logical equivalence")
    exp.predict("P4", "the invalid CNOT translates Walsh masks and preserves support count")
    exp.must_fail("C1", "the invalid guard changes legitimate arithmetic")
    exp.must_fail("C2", "the invalid guard changes the odd-order output distribution")
    rows, changed = [], False
    rng = np.random.default_rng(20260909)
    for base in (6, 3):
        me = ToffoliModExp(N=7, a=base, n_exp=4)
        qc = me.build(init_x=False)
        perm = replay(qc)
        es = np.arange(1 << me.n_exp)
        clean = ((es[:, None] << me.exp[0]) | (np.arange(me.N)[None, :] << me.x[0])).reshape(-1)
        expected = np.array([(int(e) << me.exp[0]) | ((x*pow(base, int(e), me.N)) % me.N << me.x[0])
                             for e in es for x in range(me.N)])
        assert np.array_equal(perm[clean], expected), "baseline arithmetic/cleanup failure"
        Q, dim = 1 << me.n_exp, 1 << me.n_qubits
        chi = (1 - 2*((perm >> me.x[0]) & 1)).astype(float)
        baseline_coeffs = walsh.wht(chi) / dim
        baseline_support = int(np.count_nonzero(baseline_coeffs))
        initial_indices = (es << me.exp[0]) | (1 << me.x[0])
        baseline_scalar = chi[initial_indices]
        baseline_restricted = chi[clean]
        arithmetic_state = np.zeros(dim, complex)
        arithmetic_state[perm[initial_indices]] = 1/np.sqrt(Q)
        qft = Circuit(me.n_qubits).qft(me.exp, inverse=True)
        reference = np.sum(np.abs(statevec.run(qft, arithmetic_state).reshape(Q, -1))**2, axis=1)
        coherent = np.zeros(dim, complex)
        amplitudes = rng.normal(size=clean.size) + 1j*rng.normal(size=clean.size)
        amplitudes /= np.linalg.norm(amplitudes)
        coherent[clean] = amplitudes
        for count in range(me.m + 1):
            guard = Circuit(me.n_qubits)
            for j in range(count):
                guard.toffoli(me.b[j], me.t[j], me.x[0])
            wp = replay(guard)
            modified = wp[perm]
            sign = (1 - 2*((modified >> me.x[0]) & 1)).astype(float)
            coeffs = walsh.wht(sign) / dim
            full_support = int(np.count_nonzero(coeffs))
            restricted_error = float(np.max(np.abs(sign[clean]-baseline_restricted)))
            scalar_error = float(np.max(np.abs(sign[initial_indices]-baseline_scalar)))
            exp.check("P1", np.array_equal(modified[clean], expected)
                      and restricted_error == 0 and scalar_error == 0,
                      f"a={base} guards={count}: all {clean.size} logical basis inputs agree")
            evolved = statevec.run(guard, coherent)
            overlap = np.vdot(coherent, evolved)
            phase = overlap / abs(overlap)
            coherent_error = float(np.max(np.abs(evolved-phase*coherent)))
            exp.check("P2", coherent_error < 1e-11,
                      f"gate-level coherent code error={coherent_error:.2e}")
            output = statevec.run(qft, statevec.run(guard, arithmetic_state))
            probabilities = np.sum(np.abs(output.reshape(Q, -1))**2, axis=1)
            output_error = float(np.max(np.abs(probabilities-reference)))
            exp.check("P3", output_error < 1e-10,
                      f"a={base} guards={count}: distribution error={output_error:.2e}")
            changed |= full_support != baseline_support
            rows.append(dict(a=base, N=me.N, width=me.n_exp, guards=count,
                             qubits=me.n_qubits, full_support=full_support,
                             full_operator_changed=int(np.count_nonzero(sign-chi)),
                             exponent_support=int(np.count_nonzero(walsh.wht(sign[initial_indices]))),
                             restricted_error=restricted_error, coherent_error=coherent_error,
                             output_error=output_error))
            exp.log(f"full support={full_support}, legitimate exponent support={rows[-1]['exponent_support']}")
        bad = Circuit(me.n_qubits).cnot(me.exp[0], me.x[0])
        bad_perm = replay(bad)[perm]
        bad_chi = (1 - 2*((bad_perm >> me.x[0]) & 1)).astype(float)
        bad_coeffs = walsh.wht(bad_chi) / dim
        exp.fail_check("C1", not np.array_equal(bad_perm[clean], expected), f"a={base}")
        translated = baseline_coeffs[np.arange(dim) ^ (1 << me.exp[0])]
        exp.check("P4", np.array_equal(bad_coeffs, translated)
                  and np.count_nonzero(bad_coeffs) == baseline_support,
                  f"a={base}: exactly translated Walsh vector, same support count")
        output = statevec.run(qft, statevec.run(bad, arithmetic_state))
        bad_prob = np.sum(np.abs(output.reshape(Q, -1))**2, axis=1)
        bad_error = float(np.max(np.abs(bad_prob-reference)))
        if base == 3:
            exp.fail_check("C2", bad_error > 1e-3, f"same-count invalid circuit output error={bad_error}")
        rows.append(dict(kind="invalid_control", N=me.N, a=base, width=me.n_exp,
                         full_support=int(np.count_nonzero(bad_coeffs)), output_error=bad_error))
    exp.check("H1", changed, "at least one clean-equivalent guard changes full support")
    exp.finish(report_path="out/scratch_equivalence.json", rows=rows,
               metadata=dict(numpy=np.__version__, gpu=accel.enabled(),
                             valid_code="ancillas zero, 0<=x<N, all exponent values"))


if __name__ == "__main__":
    main()
