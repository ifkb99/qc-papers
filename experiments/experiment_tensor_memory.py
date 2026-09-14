"""Walsh sparsity versus tensor memory, 2026-09-09 (TODO 15b).

DERIVATION BEFORE MEASUREMENT. The normalized Walsh transform is a product
of local orthogonal matrices, so every fixed bipartition has identical
singular values before and after transformation. An ideal period-r scalar
function has cut rank <= r via its running-residue realization. Neither
statement bounds the full scratch-space circuit function by r.

P1: singular values agree; P2: ideal periodic functions obey rank<=r and
the automaton reproduces them; P3: measured ranks are stable across tolerances.
H1/H2: full-space and intermediate ranks may be small OR may grow; report
either outcome without promoting a finite scan into a compression theorem.
C1: a seeded random sign table must defeat the small (<=6) rank hypothesis.

Only selected intermediate suffixes are scanned: no claim of measuring the
maximum rank over all gates. Uses existing replay and Walsh routines.
Run: OPENBLAS_NUM_THREADS=1 LAB_GPU=1 uv run python -m experiments.experiment_tensor_memory
"""
from __future__ import annotations
import numpy as np
from circuits import Circuit
from lab import Experiment, bit_table, order, random_boolean
from lab.tensors import cut_spectrum, rank_profile, residue_automaton
from toffoli_arith import ToffoliModExp
import accel
import walsh


def suffix(circuit, length):
    qc = Circuit(circuit.n)
    for op in circuit.logical[len(circuit.logical)-length:] if length else []:
        getattr(qc, op[0])(*op[1:])
    return qc


def replay(circuit):
    return (accel.classical_permutation(circuit) if accel.enabled()
            else walsh.classical_permutation(circuit))


def diagnostics(exp, chi, label, order_bits=None):
    coeffs = walsh.wht(chi) / chi.size
    profile = rank_profile(chi, order_bits)
    n = chi.size.bit_length()-1
    order_bits = list(range(n)) if order_bits is None else order_bits
    # Representative small, middle and large cuts. Equal singular values
    # also imply exact mathematical rank invariance; numerical ranks are
    # reported separately with a tolerance sweep.
    errors = []
    for k in sorted({1, n//2, n-1} - {0}):
        before = cut_spectrum(chi / np.sqrt(chi.size), order_bits[:k])
        after = cut_spectrum(coeffs, order_bits[:k])
        errors.append(float(np.max(np.abs(before-after))))
    exp.check("P1", max(errors, default=0.) < 1e-11,
              f"{label}: normalized cut-spectrum error={max(errors, default=0.):.2e}")
    stable = all(r["rank"] == r["rank_loose"] == r["rank_tight"] for r in profile)
    exp.check("P3", stable, f"{label}: ranks stable over relative tolerances 1e-8..1e-12")
    return dict(label=label, qubits=n, support=int(np.count_nonzero(coeffs)),
                max_sampled_cut_rank=max((r["rank"] for r in profile), default=1),
                profile=profile, max_spectrum_error=max(errors, default=0.))


def main():
    exp = Experiment("tensor_memory", doc=__doc__)
    exp.predict("P1", "local Walsh transforms preserve normalized cut singular values")
    exp.predict("P2", "period-r scalar functions have rank<=r and a residue realization")
    exp.predict("P3", "reported numerical ranks are tolerance stable")
    exp.predict("H1/H2", "measure small OR growing full-space/intermediate ranks; no committed outcome")
    exp.must_fail("C1", "random signs must defeat rank<=6")
    rows = []
    for N, a in [(7, 3), (13, 4), (15, 7)]:
        r = order(a, N)
        table = bit_table(N, a, 0)
        for width in (4, 8, 12, 16):
            g = table[np.arange(1 << width) % r]
            chi = (1 - 2*g).astype(float)
            row = diagnostics(exp, chi, f"ideal N={N} a={a} t={width}")
            samples = range(min(1 << width, 257))
            matches = all(residue_automaton(table, [(e >> b) & 1
                               for b in range(width-1, -1, -1)]) == g[e] for e in samples)
            exp.check("P2", row["max_sampled_cut_rank"] <= r and matches,
                      f"period={r}, rank={row['max_sampled_cut_rank']}, automaton agrees")
            row.update(kind="ideal", N=N, a=a, width=width, period=r)
            rows.append(row)
    random_chi = (1 - 2*random_boolean(12, rng=20260909)).astype(float)
    row = diagnostics(exp, random_chi, "fixed random 12-bit function")
    exp.fail_check("C1", row["max_sampled_cut_rank"] > 6,
                   f"random maximum cut rank={row['max_sampled_cut_rank']}")
    row["kind"] = "random_control"
    rows.append(row)

    for base in (6, 3):
        for width in (2, 3, 4):
            me = ToffoliModExp(N=7, a=base, n_exp=width)
            qc = me.build()
            perm = replay(qc)
            # This is the existing, independently state-vector-gated builder,
            # not a new variant. Check every legitimate input against pow()
            # using the full permutation already needed for this diagnostic.
            e = np.arange(1 << width)
            inputs = e << me.exp[0]  # build() initializes x from zero to one
            expected = inputs | np.array([pow(base, int(j), 7) << me.x[0] for j in e])
            assert np.array_equal(perm[inputs], expected), "arithmetic or scratch cleanup failed"
            chi = (1 - 2*((perm >> me.x[0]) & 1)).astype(float)
            row = diagnostics(exp, chi, f"full N=7 a={base} t={width}")
            row.update(kind="full", N=7, a=base, width=width)
            rows.append(row)
            exp.log(f"full support={row['support']}, max cut rank={row['max_sampled_cut_rank']}")
            if width == 3:
                block_len = len(me.u_a(me.exp[-1], pow(base, 1 << (width-1), 7)).logical)
                for length in sorted({0, block_len//4, block_len//2, 3*block_len//4, block_len}):
                    sqc = suffix(qc, length)
                    schi = (1 - 2*((replay(sqc) >> me.x[0]) & 1)).astype(float)
                    sr = diagnostics(exp, schi, f"suffix a={base} gates={length}")
                    sr.update(kind="suffix", N=7, a=base, width=width,
                              reverse_gates=length, total_gates=len(qc.logical))
                    rows.append(sr)
                    exp.log(f"sampled suffix support={sr['support']}, max cut rank={sr['max_sampled_cut_rank']}")
    exp.check("H1/H2", any(r["kind"] == "suffix" for r in rows),
              "full and selected intermediate profiles recorded; magnitudes are measurements")
    exp.finish(report_path="out/tensor_memory.json", rows=rows,
               metadata=dict(numpy=np.__version__, gpu=accel.enabled(),
                             ranks="numerical over R, not GF(2)", order="little-endian natural"))


if __name__ == "__main__":
    main()
