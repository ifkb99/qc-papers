"""Balanced controls versus full Walsh complexity, 2026-09-09 (TODO 15a).

DERIVATION BEFORE MEASUREMENT. For independent uniform windows controlling the
same work involution V through phi, let p=Pr(phi=1), psi=(-1)^phi, and
O+/-=(O +/- V^dag O V)/2. After K windows the full diagonal function is
O+ + O- product psi; its window average is O+ + (1-2p)^K O-.
Balanced phi gives the same idempotent reduced map regardless of Walsh support.
The supremum-norm error from replacing the tail by O+ is exactly
abs(1-2p)^K * norm(O-, infinity). Constant controls need not converge.

P1: all reduced entries and the full support-factorization formula agree.
P2: parity and nonlinear majority have identical nonzero reduced operators.
P3: biased OR has the predicted residual, including alternating sign.
C1: OR must fail one-step projection; C2: a constant control must not mix.

Toy uses CNOT on two work bits. Arithmetic uses the actual u_a(control,1)
permutation on all work states, extracted with the existing Walsh replay.
Composition is at the function/block level, NOT a new low-level compiler or
propagator. Gate-level PPS independently checks the one-control arithmetic
projection. The earlier toy was explored in conversation; this saved version
adds arithmetic checks and explicit support/error formulas.

Run: uv run --no-project --python 3.12 --with 'numpy<2.5' python -m
     experiments.experiment_balanced_controls
"""
from __future__ import annotations
import numpy as np
from lab import Experiment
from perm_pps import propagate_perm
from toffoli_arith import ToffoliModExp
import walsh


def control_table(name, width=3):
    ids = np.arange(1 << width)
    wt = np.array([int(i).bit_count() for i in ids])
    if name == "parity":
        return wt % 2
    if name == "majority":
        return (wt > width // 2).astype(int)
    if name == "or":
        return (wt > 0).astype(int)
    if name == "constant":
        return np.ones_like(ids)
    raise ValueError(name)


def identity_branches():
    me = ToffoliModExp(N=7, a=6, n_exp=1)
    qc = me.u_a(me.exp[0], 1)
    perm = walsh.classical_permutation(qc)
    dim = 1 << me.exp[0]
    ids = np.arange(dim)
    assert np.array_equal(perm[:dim], ids), "inactive control is not identity"
    v = perm[dim:] - dim
    assert np.array_equal(v[v], ids), "active work map is not an involution"
    f0 = (1 - 2 * ((ids >> me.x[0]) & 1)).astype(float)
    f1 = f0[v]
    return me, qc, f0, f1


def main():
    exp = Experiment("balanced_controls", doc=__doc__)
    exp.predict("P1", "full support factorization and every reduced entry agree")
    exp.predict("P2", "balanced controls give identical nonvacuous projections")
    exp.predict("P3", "OR residual is exactly (-3/4)^K O-")
    exp.must_fail("C1", "OR must fail one-step idempotent projection")
    exp.must_fail("C2", "constant control must retain a nonzero residual")
    me, qc, arith0, arith1 = identity_branches()
    rp = propagate_perm(qc, 1 << me.x[0], trace_plus=me.exp)
    expected = walsh.wht((arith0 + arith1) / 2) / len(arith0)
    actual = np.zeros_like(expected)
    for z, value in rp.final_terms.items():
        actual[z] = value
    exp.check("P1", not rp.hit_cap and np.array_equal(expected, actual),
              "arithmetic one-control gate PPS equals independent projection")

    y = np.arange(4)
    f0 = (1 - 2*((y >> 1) & 1)).astype(float)
    f1 = (1 - 2*(((y >> 1) & 1) ^ (y & 1))).astype(float)
    rows = []
    for label, f0, f1, max_k in [("CNOT", f0, f1, 4),
                                 ("arithmetic_identity", arith0, arith1, 3)]:
        plus, minus = (f0+f1)/2, (f0-f1)/2
        assert np.any(plus) and np.any(minus), "vacuous invariant or observable"
        minus_support = int(np.count_nonzero(walsh.wht(minus)))
        for name in ("parity", "majority", "or", "constant"):
            phi = control_table(name)
            psi = 1 - 2*phi
            bias = float(psi.mean())
            phi_support = int(np.count_nonzero(walsh.wht(psi)))
            for k in range(1, max_k + 1):
                activation = np.array([0])
                for _ in range(k):
                    activation = (phi[:, None] ^ activation[None, :]).reshape(-1)
                full = np.where(activation[:, None], f1[None, :], f0[None, :])
                reduced = full.mean(axis=0)
                predicted = plus + bias**k * minus
                reduced_coeffs = walsh.wht(reduced) / reduced.size
                full_support = int(np.count_nonzero(walsh.wht(full.reshape(-1))))
                reduced_support = int(np.count_nonzero(reduced_coeffs))
                # The zero-window Fourier sector contains only the reduced
                # operator; every other sector contains a scaled O-.
                nonzero_sectors = phi_support**k - int(bias != 0)
                support_prediction = reduced_support + minus_support * nonzero_sectors
                error = float(np.max(np.abs(reduced-predicted)))
                exp.check("P1", error == 0 and full_support == support_prediction,
                          f"{label} {name} K={k}: support={full_support}, reduced={reduced_support}, error={error}")
                if name in ("parity", "majority"):
                    exp.check("P2", np.array_equal(reduced, plus), f"{label} {name} K={k}")
                residual = float(np.max(np.abs(reduced-plus)))
                bound = abs(bias)**k * float(np.max(np.abs(minus)))
                if name == "or":
                    exp.check("P3", residual == bound, f"K={k}: exact norm error={residual}")
                    if k == 1:
                        exp.fail_check("C1", residual > 0, label)
                elif name == "constant" and k == max_k:
                    exp.fail_check("C2", residual > 0, label)
                rows.append(dict(family=label, control=name, windows=k,
                                 p=float(phi.mean()), control_walsh=phi_support,
                                 full_support=full_support, reduced_support=reduced_support,
                                 vector_error=error, projection_error=residual,
                                 predicted_error=bound))
    exp.finish(report_path="out/balanced_controls.json", rows=rows,
               metadata=dict(numpy=np.__version__, window_width=3,
                             task="full diagonal block function and reduced operator"))


if __name__ == "__main__":
    main()
