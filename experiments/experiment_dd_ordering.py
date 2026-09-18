"""Does the DD-native PPS peak-node turn point track C103's alpha + mu(beta)?

C103 proves the exponent-first ROBDD of the modexp pullback has widths at most
2^alpha D_mu(k - alpha), and that this bound is the TRIVIAL 2^k for every
k <= alpha + mu(beta).  So exactly min(t, alpha+mu(beta)+1) exponent levels are
trivial, and the first t at which ANY level is non-trivial is

    t* = alpha + mu(beta) + 2,     r = ord_N(a) = beta * 2^alpha, beta odd.

Separately, the DD-native PPS route (frozen pilot dd-pilot-20260915, promoted
under board task T163d291c2cdd4cd0, asymptotics accepted as S18ba76c1588945c1)
shows the peak-node ratio r_plain(t) = DD peak nodes / dictionary peak terms
RISING then FALLING, with a maximum at n_exp = 4 for N=7, a=3 over six points:
0.2808, 0.5185, 0.6343, 0.6496, 0.5797, 0.4729.

Nobody has a mechanism for that turn.  C103's crossover is a hypothesis with a
numerical coincidence attached: for N=7, a=3 we have r=6, beta=3, alpha=1,
mu(3)=1, so t* = 4, which is where the ratio turns.  ONE agreeing instance is
not evidence.  This experiment is the cheapest available discriminator: an
ORDERING across three instances, each a one-parameter change, whose predicted
turns are DIFFERENT and are fixed in advance by arithmetic alone.

    family        r   beta alpha mu   predicted turn t* = alpha+mu+2
    N=7,  a=2     3     3    0    1        3
    N=7,  a=3     6     3    1    1        4   (ALREADY MEASURED = 4)
    N=5,  a=2     4     1    2    1        5

An ordering across three instances is far harder to satisfy by accident than a
single turn.  A coincidence has to be right three times, in a specified order,
with the spacing fixed.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  COMMITTED.  ABSOLUTE.  For each family, argmax_t r_plain(t) equals
      t* = alpha + mu(beta) + 2 exactly: 3 for N=7 a=2, 4 for N=7 a=3, 5 for
      N=5 a=2.  The argmax must be INTERIOR to the measured range (a maximum
      at the last measured t is not a turn, it is an unfinished curve) or P1
      fails for that family.

  P2  COMMITTED.  ORDERING ONLY, weaker and separately resolved.
      turn(N=7,a=2) < turn(N=7,a=3) < turn(N=5,a=2).
      P2 can hold while P1 fails: that is the case where the turn tracks
      alpha + mu(beta) up to an additive constant other than 2.  Splitting
      them means a wrong constant does not masquerade as a refutation of the
      mechanism, and a right ordering does not masquerade as confirmation of
      the formula.

  P3  COMMITTED.  INSTRUMENT REPRODUCTION.  This file's driver, re-measuring
      N=7, a=3 at n_exp = 1..4, reproduces the accepted rows of
      S18ba76c1588945c1 EXACTLY (integer node counts and dictionary peaks,
      not approximately): dict_peak 8194, 24412, 48855, 98018 and plain DD
      peak 2301, 12657, 30989, 63673.  Same frozen propagator, same rules, so
      anything other than exact equality means this driver is not the accepted
      instrument and every other number here is void.

  P4  NO COMMITMENT.  The WEIGHTED (QMDD) ratio.  The accepted weighted series
      0.223, 0.417, 0.390, 0.440, 0.373, 0.290 has its maximum at n_exp=4 as
      well, but it is not monotone below the turn, so its argmax is a less
      stable statistic.  Recorded for every family, resolved as a check only
      in the sense of being reported.  No band is preregistered.

CONTROLS (must fail).

  C1  DENOMINATOR ARTIFACT.  The turn must be a property of the NUMERATOR.
      The control asserts that the dictionary peak itself turns, i.e. that
      dict_peak(t) has an interior maximum in some family.  It must fail:
      dict_peak must be strictly increasing in t everywhere.  If the baseline
      turned where the ratio turns, the effect would be an artifact of the
      denominator and nothing about decision diagrams.

  C2  GENERIC-DENSITY ARTIFACT.  The turn must be a property of the ARITHMETIC,
      not of support size.  For every measured point the pilot draws a MATCHED
      NULL: a random coefficient function with the same qubit count, the same
      support cardinality and the same value multiset as the real peak
      diagram.  The control asserts that the null ratio null_median/dict_peak
      has its interior maximum at the same t as the real ratio, in any family
      that turns.  It must fail.  If a structureless function with identical
      size statistics turned in the same place, the turn would be about
      density, and C103 would be irrelevant to it.

  C3  THE PIPELINE CAN REPORT DISAGREEMENT.  Propagate N=7, a=2, n_exp=1 with
      the Toffoli branch rule's normalisation deliberately broken (v0+v1+v2-v3
      instead of (v0+v1+v2-v3)/2, the frozen mutants.py 'norm' mutant) and
      require the cross-check against perm_pps.propagate_perm to FAIL.  It
      must fail.  Three separate measuring instruments failed in this line of
      work and every must-fail control passed while they were broken, so a
      control that exercises the actual comparison path is not optional.

Run (Python 3.12 per CLAUDE.md; 3.14 dies intermittently on long perm_pps jobs):
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      experiments/experiment_dd_ordering.py run N A TMAX
  ... once per family, then:
  uv run --no-project --python 3.12 --with 'numpy<2.5' python -u \
      experiments/experiment_dd_ordering.py verdict

Scope and limits, stated before the result is known:
  * All three families have mu(beta) = 1, so the ordering varies ALPHA only.
    This test cannot distinguish alpha + mu + 2 from alpha + f(mu) + 2 for any
    f with f(1) = 1.  A mu-varying point needs beta = 5 (mu = 2), whose
    smallest instance is N=11, a=3 at 4-bit work width: q = 17..21 over the
    range needed, roughly an order of magnitude more compute than all of this.
  * N=5, a=2 has beta = 1, which is C103's DEGENERATE family (the tail is one
    repeated involution and the ROBDD is exactly affine in t).  mu(1) = 1 is a
    convention adopted here, not a value C103 states.  At 3-bit moduli alpha=2
    forces beta=1, so this is unavoidable at this cost; it is a real weakness
    of the third point and not a detail.
  * Plain and weighted NODE counts at the propagation peak.  Not bytes.
  * One variable order (index-ascending), as in the accepted asymptotics.
  * A coincidence surviving three points is still a coincidence, not a
    mechanism.  Agreement here licenses the mu-varying test, nothing more.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "out/agent-board/artifacts/dd-pilot-20260915"
OUT = ROOT / "out/dd-ordering"
ACCEPTED = ROOT / "out/agent-board/workers/A4201728b909c4d58/rows_v6.json"

sys.setrecursionlimit(100000)
sys.path.insert(0, str(ROOT))
sys.path.append(str(FROZEN))

DDPROP_SHA = "d00cc1233a4191753017c86bf9d0457d331363f3ca8873f601205907b8720105"

# The accepted N=7 a=3 rows this driver must reproduce (S18ba76c1588945c1).
ACCEPTED_A3 = {1: (8194, 2301), 2: (24412, 12657),
               3: (48855, 30989), 4: (98018, 63673)}

FAMILIES = [(7, 2), (7, 3), (5, 2)]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_frozen() -> dict:
    """Refuse to run against anything but the sealed pilot bytes."""
    got = sha256(FROZEN / "ddprop.py")
    if got != DDPROP_SHA:
        raise SystemExit(f"frozen ddprop.py hash mismatch: {got}")
    manifest = {}
    for line in (FROZEN / "MANIFEST.sha256").read_text().splitlines():
        h, name = line.split()
        manifest[name] = h
        actual = sha256(FROZEN / name)
        if actual != h:
            raise SystemExit(f"frozen {name} hash mismatch: {actual} != {h}")
    return manifest


MANIFEST = verify_frozen()

from lab.harness import Experiment                 # noqa: E402
from perm_pps import propagate_perm, _BinaryFrame  # noqa: E402
from toffoli_arith import ToffoliModExp            # noqa: E402
from ddprop import MTBDD                           # noqa: E402  (FROZEN)


# ---------------------------------------------------------------- arithmetic --
def mult_order(a: int, N: int) -> int:
    r, x = 1, a % N
    while x != 1:
        x = x * a % N
        r += 1
    return r


def beta_alpha(r: int) -> tuple[int, int]:
    alpha = 0
    while r % 2 == 0:
        r //= 2
        alpha += 1
    return r, alpha


def mu(beta: int) -> int:
    """Order of 2 in (Z/beta)^x / {+-1}.  mu(1) = 1: see the header's limits."""
    if beta == 1:
        return 1
    s, x = 1, 2 % beta
    while x % beta not in (1, beta - 1):
        x = x * 2 % beta
        s += 1
    return s


def predicted_turn(N: int, a: int) -> dict:
    r = mult_order(a, N)
    b, al = beta_alpha(r)
    m = mu(b)
    return dict(N=N, a=a, order=r, beta=b, alpha=al, mu=m, t_star=al + m + 2)


# ------------------------------------------------------- counting and driver --
def weighted_count(dd: MTBDD, root: int) -> int:
    """QMDD-style class count: nodes equal up to a nonzero scalar are one
    class, the scalar moving to the incoming edge.  Plus one terminal."""
    classes: dict = {}
    memo: dict = {}
    order, stack = [], [(root, False)]
    while stack:                       # iterative post-order, no deep recursion
        nid, expanded = stack.pop()
        if nid in memo:
            continue
        if expanded:
            order.append(nid)
            continue
        stack.append((nid, True))
        if not dd.is_terminal(nid):
            _, lo, hi = dd._nodes[nid]
            stack.append((hi, False))
            stack.append((lo, False))
    for nid in order:
        if nid in memo:
            continue
        if dd.is_terminal(nid):
            v = dd.value(nid)
            memo[nid] = (("Z",), Fraction(0)) if v == 0 else (("T",), v)
            continue
        lv, lo, hi = dd._nodes[nid]
        (clo, slo), (chi, shi) = memo[lo], memo[hi]
        scale = slo if slo != 0 else shi
        key = (lv, clo, slo / scale, chi, shi / scale)
        cid = classes.setdefault(key, len(classes))
        memo[nid] = ((cid,), scale)
    return len(classes) + 1


def dd_propagate(qc, zmask, mutant: str | None = None):
    """Heisenberg back-propagation on the FROZEN MTBDD, tracking the peak.

    CNOT is absorbed by the binary Walsh frame and never touches the diagram;
    X is a parity sign; Toffoli is the frozen scatter rule.  `mutant` selects a
    deliberately wrong Toffoli terminal rule for control C3.
    """
    dd = MTBDD(qc.n, None)
    frame = _BinaryFrame(qc.n)
    root = dd.single_key(zmask, 1)
    peak_plain, peak_root, peak_step = dd.node_count(root), root, 0
    peak_w = weighted_count(dd, root)
    step = 0
    t0 = time.perf_counter()
    for fs in range(len(qc.logical) - 1, -1, -1):
        op = qc.logical[fs]
        if op[0] == "cnot":
            frame.cnot(op[1], op[2])
            continue
        if op[0] == "x":
            root = dd.parity_sign(root, frame.rows[op[1]])
        else:
            a, b, t = op[1], op[2], op[3]
            A, B, pred = (frame.inverse_columns[a], frame.inverse_columns[b],
                          frame.rows[t])
            root = (dd.toffoli(root, A, B, pred) if mutant is None
                    else toffoli_mutant(dd, root, A, B, pred, mutant))
        step += 1
        c = dd.node_count(root)
        if c > peak_plain:
            peak_plain, peak_root, peak_step = c, root, step
        peak_w = max(peak_w, weighted_count(dd, root))
    return dict(dd=dd, frame=frame, root=root, peak_plain=peak_plain,
                peak_weighted=peak_w, peak_root=peak_root,
                peak_step=peak_step, steps=step,
                nodes_ever_created=len(dd._nodes),
                seconds=time.perf_counter() - t0)


def toffoli_mutant(dd, root, A, B, pred, mutant):
    """The frozen mutants.py rule with the stated defect (control C3)."""
    memo = {}

    def rec(d0, d1, d2, d3, i, par):
        if i == dd.n:
            if not par:
                return dd.terminal(dd.value(d0))
            v0, v1, v2, v3 = (dd.value(d0), dd.value(d1),
                              dd.value(d2), dd.value(d3))
            v = (v0 + v1 + v2 - v3) if mutant == "norm" else (v0 + v1 + v2 - v3) / 2
            return dd.terminal(v)
        key = (d0, d1, d2, d3, i, par)
        if key in memo:
            return memo[key]
        v = dd.order[i]
        av, bv, pv = (A >> v) & 1, (B >> v) & 1, (pred >> v) & 1
        out = [rec(dd.child(d0, i, x), dd.child(d1, i, x ^ av),
                   dd.child(d2, i, x ^ bv), dd.child(d3, i, x ^ av ^ bv),
                   i + 1, par ^ (x & pv)) for x in (0, 1)]
        memo[key] = dd.mk(i, out[0], out[1])
        return memo[key]

    return rec(root, root, root, root, 0, 0)


def agrees(got: dict, ref) -> bool:
    ref = dict(ref)
    if got.keys() != ref.keys():
        return False
    return all(abs(float(got[k]) - ref[k]) < 1e-9 for k in ref)


# ------------------------------------------------------------ matched null --
def null_ratio(peak_terms: dict, n: int, draws: int = 3, seed: int = 20260916):
    """Median plain node count of random functions with the SAME qubit count,
    support cardinality and value multiset as the real peak diagram."""
    import random
    rng = random.Random(seed ^ (n << 8) ^ len(peak_terms))
    values = list(peak_terms.values())
    counts = []
    for _ in range(draws):
        keys = rng.sample(range(1 << n), len(values))
        vals = values[:]
        rng.shuffle(vals)
        dd = MTBDD(n, None)
        root = dd.terminal(Fraction(0))
        # build bottom-up from the table: cheaper and independent of the
        # propagator's node table
        table = [Fraction(0)] * (1 << n)
        for k, v in zip(keys, vals):
            table[k] = v
        ids = [dd.terminal(v) for v in table]
        for i in range(n - 1, -1, -1):
            ids = [dd.mk(i, ids[j], ids[j + 1]) for j in range(0, len(ids), 2)]
        root = ids[0]
        counts.append(dd.node_count(root))
    counts.sort()
    return counts[len(counts) // 2]


def measure_point(N: int, a: int, n_exp: int, log=print) -> dict:
    me = ToffoliModExp(N=N, a=a, n_exp=n_exp)
    qc = me.build()
    zmask = 1 << me.x[0]
    t0 = time.perf_counter()
    ref = propagate_perm(qc, zmask)
    t1 = time.perf_counter()
    res = dd_propagate(qc, zmask)
    terms = res["dd"].to_dict(res["root"])
    peak_terms = res["dd"].to_dict(res["peak_root"])
    row = dict(N=N, a=a, n_exp=n_exp, qubits=qc.n,
               gates=len(qc.logical),
               toffolis=sum(1 for op in qc.logical if op[0] == "toffoli"),
               dict_peak=ref.n_max, dict_final=len(ref.final_terms),
               dd_plain_peak=res["peak_plain"],
               dd_weighted_peak=res["peak_weighted"],
               nodes_ever_created=res["nodes_ever_created"],
               peak_step=res["peak_step"], nonlinear_steps=res["steps"],
               peak_support=len(peak_terms),
               dict_seconds=t1 - t0, dd_seconds=res["seconds"])
    row["ratio_plain"] = res["peak_plain"] / ref.n_max
    row["ratio_weighted"] = res["peak_weighted"] / ref.n_max
    row["dict_density"] = ref.n_max / (1 << qc.n)
    row["cross_check"] = agrees(
        {res["frame"].physical(k): v for k, v in terms.items()},
        ref.final_terms)
    nm = null_ratio(peak_terms, qc.n)
    row["null_plain_median"] = nm
    row["ratio_null_over_dict"] = nm / ref.n_max
    log(f"  N={N} a={a} t={n_exp}: {qc.n}q toff={row['toffolis']} "
        f"dict_peak={ref.n_max} plain={res['peak_plain']} "
        f"r={row['ratio_plain']:.4f} w={res['peak_weighted']} "
        f"rw={row['ratio_weighted']:.4f} null/dict="
        f"{row['ratio_null_over_dict']:.4f} xcheck={row['cross_check']} "
        f"[{row['dict_seconds']:.0f}s + {row['dd_seconds']:.0f}s]")
    return row


# ------------------------------------------------------------------ stages --
def rows_path(N: int, a: int) -> Path:
    return OUT / f"rows_N{N}_a{a}.json"


def declare(exp: Experiment) -> None:
    exp.predict("P1", "for each family argmax_t r_plain(t) equals "
                      "alpha+mu(beta)+2 exactly (3, 4, 5) and is interior")
    exp.predict("P2", "turn(N=7,a=2) < turn(N=7,a=3) < turn(N=5,a=2)")
    exp.predict("P3", "this driver reproduces the accepted N=7 a=3 rows "
                      "exactly at n_exp=1..4")
    exp.predict("P4", "no commitment: the weighted argmax is reported, not "
                      "predicted")
    exp.must_fail("C1", "dict_peak itself has an interior maximum "
                        "(the turn would be a denominator artifact)")
    exp.must_fail("C2", "the matched-null ratio turns at the same t as the "
                        "real ratio (the turn would be about density)")
    exp.must_fail("C3", "the norm-mutated Toffoli rule still agrees with "
                        "perm_pps.propagate_perm")


def stage_run(N: int, a: int, tmax: int) -> int:
    par = predicted_turn(N, a)
    print(f"=== dd_ordering run N={N} a={a} tmax={tmax} ===")
    print(f"    r={par['order']} beta={par['beta']} alpha={par['alpha']} "
          f"mu={par['mu']} predicted turn t*={par['t_star']}")
    print(f"    frozen ddprop.py sha256 {MANIFEST['ddprop.py']}")
    print(f"    python {sys.version.split()[0]}")
    OUT.mkdir(parents=True, exist_ok=True)
    path = rows_path(N, a)
    rows = json.loads(path.read_text()) if path.is_file() else {}
    rows["params"] = par
    for t in range(1, tmax + 1):
        if str(t) in rows.get("points", {}):
            print(f"  N={N} a={a} t={t}: already measured, skipping")
            continue
        row = measure_point(N, a, t)
        rows.setdefault("points", {})[str(t)] = row
        path.write_text(json.dumps(rows, indent=1) + "\n")
    print(f"  wrote {path}")
    return 0


def argmax_interior(series: dict[int, float]) -> tuple[int | None, bool]:
    """(argmax over measured t, whether it is interior to the measured range)."""
    if not series:
        return None, False
    ts = sorted(series)
    best = max(ts, key=lambda t: series[t])
    return best, ts[0] < best < ts[-1]


def stage_verdict() -> int:
    exp = Experiment("dd_ordering", doc=__doc__, exit_on_fail=False)
    declare(exp)
    exp.log(f"frozen ddprop.py sha256 {MANIFEST['ddprop.py']}")
    exp.log(f"python {sys.version.split()[0]}")

    data = {}
    for N, a in FAMILIES:
        path = rows_path(N, a)
        if not path.is_file():
            exp.log(f"MISSING {path.name}: run stage `run {N} {a} TMAX` first")
            continue
        data[(N, a)] = json.loads(path.read_text())

    exp.section("measured rows")
    for (N, a), d in data.items():
        p = d["params"]
        exp.log(f"N={N} a={a}: r={p['order']} beta={p['beta']} "
                f"alpha={p['alpha']} mu={p['mu']} -> predicted t*={p['t_star']}")
        for t in sorted(d["points"], key=int):
            r = d["points"][t]
            exp.log(f"    t={t} q={r['qubits']} dict={r['dict_peak']} "
                    f"plain={r['dd_plain_peak']} r_plain={r['ratio_plain']:.4f} "
                    f"r_weighted={r['ratio_weighted']:.4f} "
                    f"null/dict={r['ratio_null_over_dict']:.4f} "
                    f"xcheck={r['cross_check']}")

    # --- P3 instrument reproduction (resolved first: everything else rests on it)
    exp.section("P3  this driver reproduces the accepted N=7 a=3 rows")
    got, bad = {}, []
    d = data.get((7, 3), {})
    for t, (dp, pl) in sorted(ACCEPTED_A3.items()):
        r = d.get("points", {}).get(str(t))
        if r is None:
            bad.append(f"t={t} not measured")
            continue
        got[t] = (r["dict_peak"], r["dd_plain_peak"])
        if got[t] != (dp, pl):
            bad.append(f"t={t} got {got[t]} want {(dp, pl)}")
    exp.check("P3", not bad and len(got) == len(ACCEPTED_A3),
              "exact match on all four rows" if not bad else "; ".join(bad))

    # --- cross-checks (not a declared prediction; a gate on the rows)
    xbad = [f"N={N} a={a} t={t}" for (N, a), d in data.items()
            for t, r in d.get("points", {}).items() if not r["cross_check"]]
    exp.log(f"cross-check against perm_pps.propagate_perm: "
            f"{'all points agree' if not xbad else 'DISAGREE ' + ', '.join(xbad)}")

    # --- P1 / P2 turns
    exp.section("P1  absolute turn point = alpha + mu(beta) + 2")
    turns = {}
    p1 = []
    for (N, a), d in data.items():
        series = {int(t): r["ratio_plain"] for t, r in d["points"].items()}
        best, interior = argmax_interior(series)
        turns[(N, a)] = (best, interior)
        star = d["params"]["t_star"]
        ok = interior and best == star
        p1.append(ok)
        exp.log(f"N={N} a={a}: argmax at t={best} "
                f"({'interior' if interior else 'AT THE EDGE, not a turn'}), "
                f"predicted {star} -> {'match' if ok else 'MISS'}")
    exp.check("P1", bool(p1) and all(p1) and len(p1) == len(FAMILIES),
              f"{sum(p1)}/{len(FAMILIES)} families turn exactly at alpha+mu+2")

    exp.section("P2  ordering only")
    seq = [turns.get(f, (None, False))[0] for f in FAMILIES]
    ordered = (all(t is not None for t in seq)
               and all(x < y for x, y in zip(seq, seq[1:])))
    exp.check("P2", ordered, f"turns in family order {FAMILIES} are {seq}")

    exp.section("P4  weighted convention, reported not predicted")
    for (N, a), d in data.items():
        series = {int(t): r["ratio_weighted"] for t, r in d["points"].items()}
        best, interior = argmax_interior(series)
        exp.log(f"N={N} a={a}: weighted argmax at t={best} "
                f"({'interior' if interior else 'at the edge'}), "
                f"plain argmax at t={turns.get((N, a), (None,))[0]}")
    exp.check("P4", True, "recorded; no band was preregistered")

    # --- C1 denominator
    exp.section("C1  control: the dictionary peak must not turn")
    c1 = []
    for (N, a), d in data.items():
        series = {int(t): r["dict_peak"] for t, r in d["points"].items()}
        _, interior = argmax_interior(series)
        c1.append(interior)
        exp.log(f"N={N} a={a}: dict_peak "
                f"{[d['points'][t]['dict_peak'] for t in sorted(d['points'], key=int)]}"
                f" -> {'TURNS' if interior else 'monotone'}")
    exp.fail_check("C1", not any(c1),
                   "no family's dictionary peak has an interior maximum"
                   if not any(c1) else "a dictionary peak turned: the effect "
                                       "is a denominator artifact")

    # --- C2 matched null
    exp.section("C2  control: the matched null must not turn in the same place")
    same = []
    for (N, a), d in data.items():
        real = turns.get((N, a), (None, False))
        if not real[1]:
            continue
        series = {int(t): r["ratio_null_over_dict"] for t, r in d["points"].items()}
        best, interior = argmax_interior(series)
        same.append(interior and best == real[0])
        exp.log(f"N={N} a={a}: null argmax t={best} "
                f"({'interior' if interior else 'at the edge'}), "
                f"real argmax t={real[0]}")
    exp.fail_check("C2", not (same and all(same)),
                   "the matched null does not turn where the real ratio does")

    # --- C3 mutant
    exp.section("C3  control: a broken Toffoli rule must be caught")
    me = ToffoliModExp(N=7, a=2, n_exp=1)
    qc = me.build()
    zmask = 1 << me.x[0]
    ref = propagate_perm(qc, zmask)
    res = dd_propagate(qc, zmask, mutant="norm")
    mut_ok = agrees({res["frame"].physical(k): v
                     for k, v in res["dd"].to_dict(res["root"]).items()},
                    ref.final_terms)
    exp.fail_check("C3", not mut_ok,
                   "norm-mutated Toffoli disagrees with propagate_perm, so the "
                   "cross-check path can report a disagreement"
                   if not mut_ok else "MUTANT AGREED: the cross-check is vacuous")

    rows = {f"N{N}_a{a}": d for (N, a), d in data.items()}
    ok = exp.finish(report_path=OUT / "report_dd_ordering.json", rows=rows,
                    metadata=dict(frozen_manifest=MANIFEST, python=sys.version,
                                  source_sha256=sha256(__file__)))
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
