"""Does the DD-native PPS turn point replicate within a fixed (N, beta) class?

experiment_dd_ordering refuted C103's crossover as the mechanism for the turn
in the DD-native PPS peak-node ratio.  It did not merely miss the constant; it
got the ORDER backwards:

    family      C103 t* = alpha+mu+2     MEASURED turn
    N=7, a=2            3                     5
    N=7, a=3            4                     4
    N=5, a=2            5                  (degenerate: c_k = 1 for k >= 2,
                                            dictionary peak frozen at 48972;
                                            control C1 caught it)

Two valid points remain, and they are consistent with turn = 5 - alpha.  That
is ONE free parameter fitted to TWO points and is therefore worth nothing until
it predicts something it has not seen.  At 3-bit moduli the only non-degenerate
modexp families are N=7 with beta=3, and exactly two bases are unused:

    N=7, a=4: ord=3, beta=3, alpha=0, c_k = 4,2,4,2,...   (alpha-partner of a=2)
    N=7, a=5: ord=6, beta=3, alpha=1, c_k = 5,4,2,4,2,... (alpha-partner of a=3)

Both are non-degenerate.  Every other base is either already measured or has
c_k -> 1 (a=6), and at N=5 every base is degenerate.  So these two are the
complete remaining supply of cheap valid evidence, and they put BOTH live
hypotheses on the line at once:

    instance      C103 predicts    alpha-replicate predicts
    N=7, a=4            3                    5             <-- DISCRIMINATING
    N=7, a=5            4                    4                 (agree)

a=4 is the discriminating instance: the two hypotheses differ by two full
steps there.  a=5 is a consistency point that both must clear.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  R1  COMMITTED.  ALPHA-REPLICATE.  argmax_t r_plain(t) is 5 for N=7 a=4 and
      4 for N=7 a=5, both INTERIOR to the measured range.

  R2  COMMITTED.  C103, SECOND CHANCE.  argmax_t r_plain(t) is 3 for N=7 a=4
      and 4 for N=7 a=5, both interior.  R1 and R2 cannot both pass.  Stating
      the refuted hypothesis as a live prediction rather than assuming it dead
      is the point: it was refuted on two instances, not twenty, and it gets
      the same shot at a=4 that the replacement gets.

  R3  COMMITTED.  ORDERING ONLY, weaker than either.
      turn(N=7, a=4) > turn(N=7, a=5).  True under R1 (5 > 4), false under R2
      (3 < 4).  This survives a wrong constant in whichever is right.

  R4  NO COMMITMENT.  Whether the alpha-partners agree point for point:
      r_plain for a=4 close to a=2, and a=5 close to a=3.  Reported, not
      predicted; alpha equality is not circuit equality (see control D1).

CONTROLS (must fail).

  D1  NOT THE SAME CIRCUIT.  a=4 and a=2 share (N, ord, beta, alpha), and a=5
      and a=3 likewise.  If a replicate reproduced its partner's r_plain series
      EXACTLY it would be the same computation relabelled and the replication
      would be vacuous.  The control asserts exact series equality with the
      partner.  It must fail.

  D2  DENOMINATOR ARTIFACT.  The dictionary peak must not have an interior
      maximum in either replicate family.  It must fail.  This is the control
      that caught the degenerate N=5 family in the ordering experiment, so it
      is carried forward unchanged rather than retired after one success.

  D3  THE PIPELINE CAN REPORT DISAGREEMENT.  The norm-mutated Toffoli rule must
      DISAGREE with perm_pps.propagate_perm at N=7, a=4, n_exp=1.  It must
      fail.  Carried forward for the same reason.

Run (Python 3.12 per CLAUDE.md; 3.14 dies intermittently on long perm_pps jobs):
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      experiments/experiment_dd_alpha_replicate.py run 7 4 6
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      experiments/experiment_dd_alpha_replicate.py run 7 5 6
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      experiments/experiment_dd_alpha_replicate.py verdict

Limits, stated before the result is known:
  * This is a REPLICATION within one (N, beta) class, not an extension of the
    range of arithmetic.  mu(beta) = 1 throughout, N = 7 throughout.  Passing
    R1 shows the turn tracks alpha WITHIN beta=3 at one modulus; it says
    nothing about mu, and nothing about any other modulus.
  * turn = 5 - alpha is a FITTED form.  Even if R1 passes, the "5" is fitted
    to N=7 beta=3 and has no derivation.  Two passing instances would make it
    worth deriving, not derived.
  * Plain-convention node counts at the propagation peak.  Not bytes.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out/dd-ordering"

_spec = importlib.util.spec_from_file_location(
    "_ddord", ROOT / "experiments/experiment_dd_ordering.py")
_o = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_o)                      # verifies the frozen manifest

from lab.harness import Experiment                # noqa: E402
from perm_pps import propagate_perm               # noqa: E402
from toffoli_arith import ToffoliModExp           # noqa: E402

REPLICATES = [(7, 4), (7, 5)]
PARTNER = {(7, 4): (7, 2), (7, 5): (7, 3)}
R1_TURN = {(7, 4): 5, (7, 5): 4}                  # alpha-replicate, turn = 5-alpha
R2_TURN = {(7, 4): 3, (7, 5): 4}                  # C103, turn = alpha+mu+2


def rows_path(N: int, a: int) -> Path:
    return OUT / f"rows_N{N}_a{a}.json"


def load(N: int, a: int):
    p = rows_path(N, a)
    return json.loads(p.read_text()) if p.is_file() else None


def series(d, field="ratio_plain"):
    return {int(t): r[field] for t, r in d["points"].items()}


def declare(exp: Experiment) -> None:
    exp.predict("R1", "alpha-replicate: turn is 5 for a=4 and 4 for a=5, "
                      "both interior")
    exp.predict("R2", "C103 second chance: turn is 3 for a=4 and 4 for a=5, "
                      "both interior")
    exp.predict("R3", "ordering only: turn(a=4) > turn(a=5)")
    exp.predict("R4", "no commitment: alpha-partner agreement is reported")
    exp.must_fail("D1", "a replicate's r_plain series equals its alpha-"
                        "partner's exactly (it would be the same circuit)")
    exp.must_fail("D2", "a replicate's dictionary peak has an interior maximum")
    exp.must_fail("D3", "the norm-mutated Toffoli rule still agrees with "
                        "perm_pps.propagate_perm")


def stage_run(N: int, a: int, tmax: int) -> int:
    return _o.stage_run(N, a, tmax)


def stage_verdict() -> int:
    exp = Experiment("dd_alpha_replicate", doc=__doc__, exit_on_fail=False)
    declare(exp)
    exp.log(f"frozen ddprop.py sha256 {_o.MANIFEST['ddprop.py']}")
    exp.log(f"python {sys.version.split()[0]}")

    data = {f: load(*f) for f in REPLICATES}
    missing = [f for f, d in data.items() if d is None]
    if missing:
        exp.log(f"MISSING rows for {missing}: run the `run` stage first")

    exp.section("measured rows")
    turns = {}
    for f, d in data.items():
        if d is None:
            continue
        p = d["params"]
        exp.log(f"N={f[0]} a={f[1]}: ord={p['order']} beta={p['beta']} "
                f"alpha={p['alpha']} mu={p['mu']}")
        for t in sorted(d["points"], key=int):
            r = d["points"][t]
            exp.log(f"    t={t} q={r['qubits']} dict={r['dict_peak']} "
                    f"plain={r['dd_plain_peak']} "
                    f"r_plain={r['ratio_plain']:.4f} xcheck={r['cross_check']}")
        turns[f] = _o.argmax_interior(series(d))
        exp.log(f"    argmax t={turns[f][0]} "
                f"({'interior' if turns[f][1] else 'AT THE EDGE, not a turn'})")

    xbad = [f"N={f[0]} a={f[1]} t={t}" for f, d in data.items() if d
            for t, r in d["points"].items() if not r["cross_check"]]
    exp.log(f"cross-check against perm_pps.propagate_perm: "
            f"{'all points agree' if not xbad else 'DISAGREE ' + ', '.join(xbad)}")

    def resolve(pid, want, label):
        got, ok = {}, []
        for f in REPLICATES:
            best, interior = turns.get(f, (None, False))
            got[f] = best
            ok.append(interior and best == want[f])
        exp.check(pid, bool(ok) and all(ok) and len(got) == len(REPLICATES),
                  f"{label}: predicted {[want[f] for f in REPLICATES]}, "
                  f"measured {[got[f] for f in REPLICATES]}")

    exp.section("R1  alpha-replicate (turn = 5 - alpha)")
    resolve("R1", R1_TURN, "alpha-replicate")
    exp.section("R2  C103 second chance (turn = alpha + mu + 2)")
    resolve("R2", R2_TURN, "C103")

    exp.section("R3  ordering only")
    seq = [turns.get(f, (None,))[0] for f in REPLICATES]
    exp.check("R3", all(s is not None for s in seq) and seq[0] > seq[1],
              f"turn(a=4)={seq[0]}, turn(a=5)={seq[1]}")

    exp.section("R4  alpha-partner agreement, reported not predicted")
    same_exact = []
    for f in REPLICATES:
        d, pd = data.get(f), load(*PARTNER[f])
        if d is None or pd is None:
            continue
        s, ps = series(d), series(pd)
        shared = sorted(set(s) & set(ps))
        same_exact.append(bool(shared) and all(s[t] == ps[t] for t in shared))
        exp.log(f"N={f[0]} a={f[1]} vs partner a={PARTNER[f][1]}: "
                + ", ".join(f"t={t} {s[t]:.4f} vs {ps[t]:.4f}" for t in shared))
    exp.check("R4", True, "recorded; no band was preregistered")

    exp.section("D1  control: a replicate must not be its partner relabelled")
    exp.fail_check("D1", not (same_exact and all(same_exact)),
                   "replicate series differ from their alpha-partners, so "
                   "these are distinct circuits")

    exp.section("D2  control: the dictionary peak must not turn")
    d2 = []
    for f, d in data.items():
        if d is None:
            continue
        _, interior = _o.argmax_interior(series(d, "dict_peak"))
        d2.append(interior)
        exp.log(f"N={f[0]} a={f[1]}: dict_peak "
                f"{[d['points'][t]['dict_peak'] for t in sorted(d['points'], key=int)]}"
                f" -> {'TURNS' if interior else 'monotone'}")
    exp.fail_check("D2", not any(d2),
                   "no replicate's dictionary peak has an interior maximum")

    exp.section("D3  control: a broken Toffoli rule must be caught")
    me = ToffoliModExp(N=7, a=4, n_exp=1)
    qc = me.build()
    zmask = 1 << me.x[0]
    ref = propagate_perm(qc, zmask)
    res = _o.dd_propagate(qc, zmask, mutant="norm")
    mut_ok = _o.agrees({res["frame"].physical(k): v
                        for k, v in res["dd"].to_dict(res["root"]).items()},
                       ref.final_terms)
    exp.fail_check("D3", not mut_ok,
                   "norm-mutated Toffoli disagrees with propagate_perm"
                   if not mut_ok else "MUTANT AGREED: the cross-check is vacuous")

    rows = {f"N{f[0]}_a{f[1]}": d for f, d in data.items() if d}
    ok = exp.finish(report_path=OUT / "report_dd_alpha_replicate.json",
                    rows=rows,
                    metadata=dict(frozen_manifest=_o.MANIFEST,
                                  python=sys.version,
                                  source_sha256=_o.sha256(__file__)))
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["run", "verdict"])
    ap.add_argument("N", nargs="?", type=int)
    ap.add_argument("a", nargs="?", type=int)
    ap.add_argument("tmax", nargs="?", type=int)
    ns = ap.parse_args()
    if ns.stage == "run":
        if None in (ns.N, ns.a, ns.tmax):
            ap.error("run needs N a TMAX")
        return stage_run(ns.N, ns.a, ns.tmax)
    return stage_verdict()


if __name__ == "__main__":
    sys.exit(main())
