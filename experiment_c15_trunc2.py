"""P2 continuation: incremental truncation for the remaining two cases."""
from __future__ import annotations
from toffoli_arith import ToffoliModExp
from perm_pps import propagate_perm

for N, a, label in ((5, 2, "N=5 a=2 (r=4, beta=1)"),
                    (7, 3, "N=7 a=3 (r=6, beta=3) CONTROL")):
    print(f"\n  {label}")
    print(f"    {'n_exp':>6} {'delta':>9} {'N_max':>9} {'N_final':>9} {'<O>':>11}")
    per = {}
    for ne in (3, 4, 5):
        me = ToffoliModExp(N=N, a=a, n_exp=ne)
        qc = me.build()
        for d in (0.0, 1e-4, 1e-2, 1e-1):
            r = propagate_perm(qc, 1 << me.x[0], delta=d, max_terms=8_000_000)
            per.setdefault(d, []).append((r.n_max, len(r.final_terms),
                                          round(r.expectation, 9)))
            print(f"    {ne:6d} {d:9.0e} {r.n_max:9d} {len(r.final_terms):9d} "
                  f"{r.expectation:+11.6f}", flush=True)
    print("    constancy across n_exp, per delta:")
    for d, v in per.items():
        f = lambda i: "SAME" if all(x[i] == v[0][i] for x in v) else "differs"
        print(f"      delta={d:8.0e}  N_max {f(0):8s} N_final {f(1):8s} <O> {f(2)}")
