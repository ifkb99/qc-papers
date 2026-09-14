"""Independent bounded audit of perm_pps affine-frame physical coordinates.

The Walsh route below replays each classical circuit and applies an explicit
Walsh transform.  It is kept separate from the PPS dictionary walk so that a
frame bug cannot agree with a second copy of the same propagation logic.
"""
from __future__ import annotations

import json
import random
from pathlib import Path
import sys

RESEARCH_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESEARCH_ROOT))

import numpy as np

from circuits import Circuit
from lab.harness import Experiment
from perm_pps import propagate_perm
from walsh import classical_permutation, pullback_coefficients, wht


OUT = RESEARCH_ROOT / "out"
REPORT = OUT / "cnot_frame_laws_main.json"


def general_walsh(circuit: Circuit, zmask: int) -> np.ndarray:
    """Independent truth-table Walsh coefficients for any Z mask."""
    perm = classical_permutation(circuit)
    # numpy arrays do not expose int.bit_count; use parity of popcount.
    phase = np.where(np.array([int(v).bit_count() & 1 for v in (perm & zmask)]),
                     -1.0, 1.0)
    return wht(phase) / (1 << circuit.n)


def mapping_dict(view):
    return {int(k): float(v) for k, v in view.items()}


def close_maps(actual, expected, tol=2e-12):
    keys = set(actual) | {int(i) for i, c in enumerate(expected)
                           if abs(float(c)) > tol}
    bad = []
    for k in sorted(keys):
        a = float(actual.get(k, 0.0))
        e = float(expected[k]) if 0 <= k < len(expected) else 0.0
        if abs(a - e) > tol:
            bad.append((k, a, e))
    return bad


def make_fixture(rng: random.Random, n: int, depth: int,
                 force_first_cnot: bool = False) -> Circuit:
    c = Circuit(n)
    if force_first_cnot:
        c.cnot(0, 1)
    ntof = 0
    for _ in range(depth - int(force_first_cnot)):
        choices = ["x", "cnot"]
        if ntof < 3 and n >= 3:
            choices.append("toffoli")
        kind = rng.choice(choices)
        if kind == "x":
            c.x(rng.randrange(n))
        elif kind == "cnot":
            a, b = rng.sample(range(n), 2)
            c.cnot(a, b)
        else:
            a, b, d = rng.sample(range(n), 3)
            c.toffoli(a, b, d)
            ntof += 1
    return c


def run():
    exp = Experiment(
        "experiment_cnot_frame",
        doc=__doc__,
    )
    exp.predict("P1", "Every retained physical coefficient from affine_frame=True matches the independent truth-table Walsh law and default PPS on bounded mixed circuits.")
    exp.predict("P2", "Positive delta, arbitrary Z masks, identity, Mapping boundaries, and trace/cap ordering preserve the documented semantics.")
    exp.must_fail("C1", "A deliberately unconverted internal-key view must disagree with physical Walsh keys on a nontrivial CNOT circuit.")

    rng = random.Random(0xC0FFEE)
    rows = []
    all_ok = True
    fixture_count = 0

    exp.section("P1 mixed fixed-seed full-coefficient fixtures")
    # Vary depth while keeping each circuit small enough that three Toffolis
    # cannot exceed the n<=7 / <100 fixture budget.
    for n in range(2, 8):
        for j in range(12):
            depth = 2 + ((j + n) % 8)
            c = make_fixture(rng, n, depth)
            target = rng.randrange(n)
            zmask = 1 << target
            affine = propagate_perm(c, zmask, delta=0.0, affine_frame=True)
            default = propagate_perm(c, zmask, delta=0.0)
            expected = pullback_coefficients(c, target)
            got = mapping_dict(affine.final_terms)
            bad = close_maps(got, expected)
            default_bad = close_maps(dict(default.final_terms), expected)
            same = close_maps(got, np.array([default.final_terms.get(i, 0.0)
                                               for i in range(1 << n)]))
            ok = not bad and not default_bad and not same
            all_ok &= ok
            fixture_count += 1
            rows.append({"kind": "mixed", "n": n, "depth": depth,
                         "target": target, "ops": len(c.logical),
                         "terms": len(got), "bad": bad[:3],
                         "default_bad": default_bad[:3],
                         "same_bad": same[:3], "ok": ok})
    exp.check("P1", all_ok and fixture_count == 72,
              f"{fixture_count} fixtures; first mismatches={next((r for r in rows if not r['ok']), None)}")

    exp.section("P2 boundary and truncation fixtures")
    p2 = True
    details = []

    # Identity/no gates has a single coefficient and exposes the read-only
    # Mapping contract without any frame updates.
    identity = Circuit(4)
    ri = propagate_perm(identity, 1 << 2, affine_frame=True)
    identity_ok = (mapping_dict(ri.final_terms) == {4: 1.0}
                   and ri.expectation == 1.0
                   and ri.frame_cnot_updates == 0)
    p2 &= identity_ok
    details.append(("identity", identity_ok))

    # Arbitrary multi-bit observable mask uses the same independent replay
    # and Walsh transform (single-bit targets are covered by pullback above).
    arb = make_fixture(rng, 6, 8)
    arbitrary_mask = 1 | (1 << 3) | (1 << 5)
    ra = propagate_perm(arb, arbitrary_mask, affine_frame=True)
    ga = mapping_dict(ra.final_terms)
    ea = general_walsh(arb, arbitrary_mask)
    arbitrary_ok = not close_maps(ga, ea)
    p2 &= arbitrary_ok
    details.append(("arbitrary_zmask", arbitrary_ok))

    # First forward gate is CNOT.  delta>1 is a deliberately sharp empty
    # support boundary, and must agree across both coordinate paths.
    first_empty = Circuit(4).cnot(0, 1)
    re_f = propagate_perm(first_empty, 1 << 1, delta=1.1,
                          affine_frame=True)
    re_d = propagate_perm(first_empty, 1 << 1, delta=1.1)
    empty_ok = (not re_f.final_terms and not re_d.final_terms
                and re_f.n_terms == re_d.n_terms == [0]
                and re_f.frame_cnot_updates == 1)
    p2 &= empty_ok
    details.append(("first_cnot_delta_gt_one_empty", empty_ok,
                    {"n_terms": re_f.n_terms, "terms": len(re_f.final_terms)}))

    # A nonempty positive truncation through CNOT and Toffoli checks that the
    # first reverse CNOT does not alter the threshold semantics.  max_weight
    # filters physical keys after each update, so compare to that explicit
    # physical-key filtered Walsh law FOR THIS FIXTURE: the first reverse
    # CNOT leaves Z3 unchanged, the sole Toffoli does all branching, and the
    # last reverse CNOT changes no produced key. Final filtering is NOT a
    # general independent reference for sequential truncation.
    truncated = Circuit(4).cnot(0, 1).toffoli(0, 2, 3).cnot(1, 2)
    trunc_delta, trunc_weight = 0.5, 2
    tr_f = propagate_perm(truncated, 1 << 3, delta=trunc_delta,
                          max_weight=trunc_weight, affine_frame=True)
    tr_d = propagate_perm(truncated, 1 << 3, delta=trunc_delta,
                          max_weight=trunc_weight)
    tr_w = pullback_coefficients(truncated, 3)
    tr_expected = np.array([float(v) if abs(float(v)) >= trunc_delta
                            and int(i).bit_count() <= trunc_weight else 0.0
                            for i, v in enumerate(tr_w)])
    trunc_ok = (len(tr_f.final_terms) > 0
                and not close_maps(mapping_dict(tr_f.final_terms), tr_expected)
                and not close_maps(dict(tr_d.final_terms), tr_expected)
                and tr_f.n_terms == tr_d.n_terms
                and tr_f.hit_cap == tr_d.hit_cap)
    p2 &= trunc_ok
    details.append(("first_cnot_nonempty_delta_weight", trunc_ok,
                    {"n_terms": tr_f.n_terms, "terms": len(tr_f.final_terms)}))

    # Cap is checked BEFORE trace.  A Toffoli with a traced control yields four
    # pre-contraction terms; max_terms=2 must report the peak and stop before
    # emitting a trace event, so contraction cannot hide the cap hit.
    split = Circuit(3).toffoli(0, 1, 2)
    rt = propagate_perm(split, 1 << 2, trace_plus=[0], max_terms=2,
                        affine_frame=True)
    rt0 = propagate_perm(split, 1 << 2, trace_plus=[0], max_terms=2)
    trace_cap_ok = (rt.hit_cap and rt0.hit_cap and rt.n_terms == [4]
                    and rt0.n_terms == [4] and not rt.trace_events
                    and not rt0.trace_events)
    p2 &= trace_cap_ok
    details.append(("cap_before_trace", trace_cap_ok,
                    {"n_terms": rt.n_terms, "events": rt.trace_events,
                     "hit_cap": rt.hit_cap}))

    # Uncapped tracing should equal the independent full Walsh law filtered in
    # the final physical input coordinate.  This fixture keeps nonzero terms
    # after tracing, testing transformed trace constraints through CNOTs and
    # Toffolis rather than only checking an empty result.
    traced = (Circuit(5).toffoli(0, 1, 2).cnot(3, 4).cnot(2, 3)
              .toffoli(0, 3, 4))
    qtrace, qtarget = 0, 4
    tv = propagate_perm(traced, 1 << qtarget, trace_plus=[qtrace],
                        affine_frame=True)
    tw = pullback_coefficients(traced, qtarget)
    t_expected = np.array([float(v) if not (i & (1 << qtrace)) else 0.0
                           for i, v in enumerate(tw)])
    trace_ok = (len(tv.final_terms) > 0
                and not close_maps(mapping_dict(tv.final_terms), t_expected)
                and tv.trace_events and tv.trace_events[-1][3] > 0)
    p2 &= trace_ok
    details.append(("trace_uncapped_physical_filter", trace_ok,
                    {"events": tv.trace_events, "terms": len(tv.final_terms)}))

    # Mapping is read-only and physical-key lookup round-trips exactly.
    lookup_ok = True
    for k, v in ri.final_terms.items():
        lookup_ok &= (ri.final_terms[k] == v)
    try:
        _ = ri.final_terms["bad"]
        lookup_ok = False
    except KeyError:
        pass
    try:
        _ = ri.final_terms[-1]
        lookup_ok = False
    except KeyError:
        pass
    try:
        ri.final_terms[4] = 2.0
        lookup_ok = False
    except (TypeError, AttributeError):
        pass
    p2 &= lookup_ok
    details.append(("mapping_boundaries", lookup_ok))

    exp.check("P2", p2, json.dumps(details, sort_keys=True))

    exp.section("C1 deliberately broken frame conversion")
    broken_c = Circuit(3).cnot(0, 1).toffoli(0, 1, 2)
    broken = propagate_perm(broken_c, 1 << 2, affine_frame=True)
    # This is intentionally wrong: treat logical backing keys as physical
    # keys, omitting FramedTerms' deferred physical conversion.
    broken_internal = dict(broken.final_terms._terms)
    broken_expected = pullback_coefficients(broken_c, 2)
    broken_bad = close_maps(broken_internal, broken_expected)
    control_failed = bool(broken_bad)
    exp.fail_check("C1", control_failed,
                   f"broken internal-key view mismatches at {broken_bad[:3]}")

    exp.finish(
        report_path=REPORT,
        rows=rows + [{"kind": "boundary", "details": details}],
        metadata={"seed": "0xC0FFEE", "fixtures": fixture_count,
                  "max_n": 7, "source": "walsh.pullback_coefficients and truth-table WHT",
                  "scope": "bounded correctness; no timing sweep"},
    )


if __name__ == "__main__":
    run()
