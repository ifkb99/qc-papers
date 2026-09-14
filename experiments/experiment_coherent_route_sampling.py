"""Independent audit of the sparse-block coherent-route sampler.

The reference here is deliberately finite and tiny: it constructs the indexed
r-by-r work matrices only for r<=14, and compares prefix amplitudes and the
complete joint (coarse-sector, QFT-output) law against those matrices.  A
second check enumerates every state of the proposed classical transition law
(initial sector/exponent, work label, route sector and output prefix), using
the production prefix oracle only as the amplitude oracle.  This is the
gate-by-gate induction test, not a new generic propagator.

Predictions (written before measurement):
P1. prefix_vector is sqrt(M) times the direct projected work vector at every
    arithmetic/background/reflection boundary and at every partial QFT prefix.
P2. deterministic sparse-block transition enumeration equals the independent
    full-r joint law and its marginal equals sequential_path.
P3. k=0, insertion-0/end/self-loop and zero-background controls remain exact;
    wide samples allocate no orbit/sector/output tables and charge histories.
C1. freezing the initial coarse sector must change a nontrivial output law.
C2. incoherent coarse-sector dephasing must change that same coherent law.
C3. deleting route cross terms (incoherent history mixture) must fail.

All dense allocations are guarded at 16 MiB.  Failed controls are retained in
timestamped JSON artifacts; no tolerance-selected amplitude is deleted.
"""
from __future__ import annotations

import itertools
import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.fourier_sampling import unit_phase
from lab.semiclassical import sequential_path


MAX_DENSE_BYTES = 16 << 20
MAX_REFERENCE_OPERATIONS = 1 << 28
# This is a scalar-visit budget for the complete float transition law, not a
# dense array allocation. The bounded t=8,r=9 audit is below it; larger calls are
# rejected before the frontier is constructed.
MAX_TRANSITION_WORK = 1_000_000_000

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "prefix amplitudes equal the independent full-r projected vectors at all boundaries and partial QFT prefixes")
exp.predict("P2", "deterministic sparse-block transition enumeration equals the full-r joint law and sequential_path marginal")
exp.predict("P3", "edge cases normalize, and wide samples remain table-free with charged history/prefix work")
exp.must_fail("C1", "freezing the initial coarse sector changes a nontrivial output law")
exp.must_fail("C2", "incoherent coarse-sector dephasing changes the coherent output law")
exp.must_fail("C3", "deleting coherent route-history cross terms changes the output law")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def reference_guard(circuit, *, max_width=4):
    """Preallocation/operation gate for the independent tiny full-r reference.

    The default width cap remains four.  Wider diagnostics must opt in with an
    explicit integer cap no larger than eight; this does not enable legacy
    history enumeration.
    """
    if type(max_width) is not int or not 0 <= max_width <= 8:
        raise ValueError("max_width must be an integer in [0,8]")
    if circuit.period > 14 or circuit.width > max_width:
        raise ValueError(f"full-r reference is capped at r<=14,width<={max_width}")
    r, q = int(circuit.period), 1 << int(circuit.width)
    guard((r, r), label="full-r reference matrix")
    guard((r, q), label="full-r output workspace")
    if r * q > (1 << 15):
        raise MemoryError("full-r reference product cap")
    t = int(circuit.width)
    # Charge matrix entries touched by setup, all Q^2 exponent/output paths,
    # sector projections, and a conservative constant for vector arithmetic.
    # This is not a CPU/bit-runtime count or an oracle cost estimate.
    estimated_ops = 4*r*r*(q*q*(t+2)+(t+1)*r)
    if estimated_ops > MAX_REFERENCE_OPERATIONS:
        raise MemoryError("estimated direct full-r operations exceed opt-in budget")
    payload = 16*((2*t+16)*r*r+8*r*q)
    if payload > MAX_DENSE_BYTES:
        raise MemoryError("full-r simultaneous arithmetic payload exceeds 16 MiB")
    return dict(estimated_direct_operations=estimated_ops,
                max_width=max_width, dense_payload_bytes_upper_bound=payload)


def stamp(prefix):
    now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{now}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"{prefix}_{now}_{serial}.json"
    return path


def rx(theta):
    return np.array([[np.cos(theta / 2), -1j*np.sin(theta / 2)],
                     [-1j*np.sin(theta / 2), np.cos(theta / 2)]], complex)


def sector_basis(period, block, alpha):
    M = period // block
    guard((period, block), label="sector basis")
    q = np.zeros((period, block), complex)
    for m in range(M):
        for p in range(block):
            q[block*m+p, p] = unit_phase(-alpha*m, M) / np.sqrt(M)
    return q


def shift_matrix(period, power):
    guard((period, period), label="orbit shift")
    out = np.zeros((period, period), complex)
    for j in range(period):
        out[(j + power) % period, j] = 1.
    return out


def repeated(period, block, W):
    guard((period, period), label="repeated defect")
    return np.kron(np.eye(period // block, dtype=complex), np.asarray(W, complex))


def route_matrix(period, block, q):
    M = period // block
    guard((period, period), label="route matrix")
    out = np.zeros((period, period), complex)
    for m in range(M):
        for p in range(block):
            out[block*((-m) % M)+p, block*m+p] = unit_phase(-q*m, M)
    return out


def fixture(period=10, width=4, route_count=2, background=True):
    defects = {1: rx(np.pi/2), 3: rx(np.pi/2)} if background else {}
    reflections = {}
    if route_count >= 1:
        reflections[2] = (0, np.pi/4)
    if route_count >= 2:
        reflections[3] = (1, np.pi/4)
    if route_count >= 3:
        reflections[4] = (0, np.pi/4)
    return CoherentReflectionCircuit(period, 2, width, defects, reflections)


def full_gates(circuit, *, max_width=4):
    """Independent full-r pair matrices and insertion-0 initial vector."""
    r, b, t = circuit.period, circuit.b, circuit.width
    reference_guard(circuit, max_width=max_width)
    identity = np.eye(r, dtype=complex)
    U0 = repeated(r, b, circuit.background.defects.get(0, circuit.background.identity))
    if 0 in circuit.reflections:
        q, theta = circuit.reflections[0]
        U0 = (np.cos(theta/2)*identity - 1j*np.sin(theta/2)*route_matrix(r, b, q)) @ U0
    pairs = []
    for i in range(t):
        W = repeated(r, b, circuit.background.defects.get(i+1, circuit.background.identity))
        if i+1 in circuit.reflections:
            q, theta = circuit.reflections[i+1]
            J = route_matrix(r, b, q)
            K = np.cos(theta/2)*identity - 1j*np.sin(theta/2)*J
            W = K @ W
        S = shift_matrix(r, 1 << i)
        pairs.append((W, W @ S))
    initial = U0[:, 0]
    return pairs, initial


def direct_prefix(circuit, final, exponent, stop, boundary, measured=0, output=0, *, max_width=4):
    """Direct full-r coherent prefix contraction from physical orbit label 0.

    Unlike the production route expansion, this applies each controlled gate
    as a full-r matrix to the coherent physical initial state, then projects
    onto the requested sector.  The sqrt(M) factor converts the normalized
    physical amplitude to the production prefix convention.
    """
    r, b, t = circuit.period, circuit.b, circuit.width
    reference_guard(circuit, max_width=max_width)
    state = np.zeros(r, complex)
    state[0] = 1.
    if stop > 0 or boundary != "arithmetic":
        state = repeated(r, b, circuit.background.defects.get(0,
                      circuit.background.identity)) @ state
        if (stop > 0 or boundary == "reflection") and 0 in circuit.reflections:
            q, theta = circuit.reflections[0]
            state = (np.cos(theta/2)*np.eye(r) -
                     1j*np.sin(theta/2)*route_matrix(r,b,q)) @ state
    for i in range(stop):
        S = shift_matrix(r, 1 << i)
        if i >= t-measured:
            phase = unit_phase(-int(output)*(1 << i), 1 << t)
            state = (state + phase*S @ state) / 2
        elif int(exponent) & (1 << i):
            state = S @ state / np.sqrt(2)
        else:
            state = state / np.sqrt(2)
        s = i+1
        if s < stop or boundary != "arithmetic":
            state = repeated(r, b, circuit.background.defects.get(s,
                          circuit.background.identity)) @ state
        if s in circuit.reflections and (s < stop or boundary == "reflection"):
            q, theta = circuit.reflections[s]
            state = (np.cos(theta/2)*np.eye(r) -
                     1j*np.sin(theta/2)*route_matrix(r,b,q)) @ state
    state *= math.exp2(-(t-stop)/2)
    return np.sqrt(r//b) * sector_basis(r,b,int(final)).conj().T @ state


def direct_joint(circuit, *, max_width=4):
    """Full-r joint p(gamma,y), from explicit controlled matrices + inverse QFT."""
    reference_guard(circuit, max_width=max_width)
    pairs, initial = full_gates(circuit, max_width=max_width)
    r, b, t = circuit.period, circuit.b, circuit.width
    M, Q = r//b, 1 << t
    out = np.zeros((M, Q), float)
    for y in range(Q):
        projected = [np.zeros(b, complex) for _ in range(M)]
        for e in range(Q):
            state = initial.copy()
            for i, (B0, B1) in enumerate(pairs):
                state = (B1 if e & (1 << i) else B0) @ state
            state *= unit_phase(-y*e, Q) / Q
            for gamma in range(M):
                projected[gamma] += sector_basis(r, b, gamma).conj().T @ state
        for gamma in range(M):
            out[gamma, y] = float(np.vdot(projected[gamma], projected[gamma]).real)
    return out


def direct_marginal_for_alpha(circuit, alpha):
    reference_guard(circuit)
    pairs, _ = full_gates(circuit)
    r, b, t = circuit.period, circuit.b, circuit.width
    Q = 1 << t
    initial = sector_basis(r, b, alpha)[:, 0]
    if 0 in circuit.background.defects:
        initial = repeated(r, b, circuit.background.defects[0]) @ initial
    if 0 in circuit.reflections:
        q, theta = circuit.reflections[0]
        I = np.eye(r, dtype=complex)
        initial = (np.cos(theta/2)*I - 1j*np.sin(theta/2)*route_matrix(r,b,q)) @ initial
    result = np.zeros(Q)
    for y in range(Q):
        state = np.zeros(r, complex)
        for e in range(Q):
            v = initial.copy()
            for i, (B0, B1) in enumerate(pairs):
                v = (B1 if e & (1 << i) else B0) @ v
            state += unit_phase(-y*e, Q) * v / Q
        result[y] = float(np.vdot(state, state).real)
    return result


def enumerate_transitions(circuit, *, max_width=4):
    """Same enumeration with explicit (sector, work-label, output-prefix)."""
    reference_guard(circuit, max_width=max_width)
    M, Q, b, t = circuit.sectors, 1 << circuit.width, circuit.b, circuit.width
    # For each of M*Q initial (sector, exponent) pairs, the pre-QFT frontier
    # has at most M*b states.  Across the partial-QFT refinement it visits
    # fewer than 2Q child labels.  Charge a conservative constant for the
    # two-child and b-coordinate weight evaluations; this is an estimate of
    # scalar visits, not a claim about simultaneous storage.
    estimated_frontier_work = (4 * M * M * b * b * Q
                               * max(1, t + 1 + 2 * Q))
    if estimated_frontier_work > MAX_TRANSITION_WORK:
        raise MemoryError("estimated transition/frontier work exceeds opt-in budget")
    # Both old/new frontiers can be live together. This explicit bookkeeping
    # allowance is not a measured Python allocator/RSS bound.
    if 2*M*b*Q*1024 > MAX_DENSE_BYTES:
        raise MemoryError("simultaneous transition frontier allowance exceeds 16 MiB")
    joint = np.zeros((M, Q), float)
    for initial_sector in range(M):
        for exponent in range(Q):
            frontier = {(initial_sector, 0, 0): 1./(M*Q)}
            for stop in range(t+1):
                if stop and (exponent & (1 << (stop-1))):
                    frontier = { (a,(p+(1 << (stop-1))) % b,y): w
                                 for (a,p,y),w in frontier.items() }
                if stop in circuit.background.defects:
                    nxt = {}
                    for (a,p,y), w in frontier.items():
                        vec = circuit.prefix_vector(a, exponent, stop, boundary="background")
                        weights = np.abs(vec)**2
                        total = float(weights.sum())
                        if not np.isfinite(total) or total <= 0:
                            raise ArithmeticError("zero/nonfinite background transition mass")
                        for pp,z in enumerate(weights):
                            if z != 0:
                                nxt[(a,pp,y)] = nxt.get((a,pp,y),0.) + w*float(z/total)
                    frontier = nxt
                if stop in circuit.reflections:
                    q,_ = circuit.reflections[stop]
                    nxt = {}
                    for (a,p,y),w in frontier.items():
                        target=(-a-q)%M
                        if target == a:
                            nxt[(a,p,y)] = nxt.get((a,p,y),0.)+w
                            continue
                        weights=np.array([abs(circuit.prefix_vector(x, exponent, stop)[p])**2
                                          for x in (a,target)])
                        total=float(weights.sum())
                        if not np.isfinite(total) or total <= 0:
                            raise ArithmeticError("zero/nonfinite reflection transition mass")
                        for x,z in zip((a,target),weights):
                            if z != 0:
                                nxt[(x,p,y)] = nxt.get((x,p,y),0.)+w*float(z/total)
                    frontier=nxt
            for measured in range(1,t+1):
                nxt={}
                for (a,p,y0),w in frontier.items():
                    children=(y0,y0 | (1 << (measured-1)))
                    weights=np.array([abs(circuit.prefix_vector(a,exponent,t,
                                      measured=measured,output=y)[p])**2 for y in children])
                    total=float(weights.sum())
                    if not np.isfinite(total) or total <= 0:
                        raise ArithmeticError("zero/nonfinite QFT transition mass")
                    for y,z in zip(children,weights):
                        if z != 0:
                            nxt[(a,p,y)] = nxt.get((a,p,y),0.)+w*float(z/total)
                frontier=nxt
            for (a,p,y),w in frontier.items():
                joint[a,y]+=w
    return joint


def maxerr(a,b):
    return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))


def prefix_audit(circuit):
    errors=[]
    for stop in range(circuit.width+1):
        boundaries = ("arithmetic", "background", "reflection")
        for boundary in boundaries:
            if stop == circuit.width and boundary == "reflection":
                measured_values=range(circuit.width+1)
            else:
                measured_values=(0,)
            for measured in measured_values:
                outs=range(1 << measured) if measured else (0,)
                for gamma in range(circuit.sectors):
                    for exponent in (0, 1, (1<<circuit.width)-1):
                        for output in outs:
                            got=circuit.prefix_vector(gamma,exponent,stop,
                                      boundary=boundary,measured=measured,output=output)
                            ref=direct_prefix(circuit,gamma,exponent,stop,boundary,measured,output)
                            errors.append(maxerr(got,ref))
    return max(errors) if errors else 0.


def audit_one(circuit):
    pref=prefix_audit(circuit)
    target=direct_joint(circuit)
    enum=enumerate_transitions(circuit)
    prod=np.array([[circuit.joint_probability(g,y) for y in range(target.shape[1])]
                   for g in range(target.shape[0])])
    pairs,initial=full_gates(circuit)
    seq=np.array([sequential_path(pairs,initial,output=y)["conditional_path_probability"]
                  for y in range(target.shape[1])])
    return dict(prefix_max_abs_error=pref, joint_max_abs_error=maxerr(prod,target),
                transition_max_abs_error=maxerr(enum,target),
                transition_mass=float(enum.sum()), direct_mass=float(target.sum()),
                marginal_vs_sequential=maxerr(enum.sum(axis=0),seq),
                sequential_mass=float(seq.sum()))


def controls(circuit, target):
    if 0 in circuit.background.defects or 0 in circuit.reflections:
        raise AssertionError("control fixture must not use insertion-0 gates")
    actual=target.sum(axis=0)
    frozen=direct_marginal_for_alpha(circuit,0)
    incoherent=np.mean([direct_marginal_for_alpha(circuit,a)
                        for a in range(circuit.sectors)],axis=0)
    # Incoherent history mixture: replace the coherent final norm by a sum of
    # history norms, retaining the correct final sector and normalization.
    incoherent_hist=np.zeros_like(target)
    r,b,t=circuit.period,circuit.b,circuit.width
    for gamma in range(circuit.sectors):
        for y in range(1<<t):
            total=0.
            for coeff,selected in circuit._histories(t,"reflection"):
                # Reconstruct each selected route branch independently; this
                # intentionally removes cross terms while retaining the exact
                # branch amplitudes and route inverse.
                current=gamma
                for q in reversed(tuple(selected.values())):
                    current=(-current-q)%(r//b)
                state=sector_basis(r,b,current)[:,0]
                for i in range(t):
                    S=shift_matrix(r,1<<i); phase=unit_phase(-y*(1<<i),1<<t)
                    state=(state+phase*S@state)/2
                    state=repeated(r,b,circuit.background.defects.get(i+1,circuit.background.identity))@state
                    if i+1 in selected:
                        state=route_matrix(r,b,selected[i+1])@state
                amp=sector_basis(r,b,gamma).conj().T@state
                total += abs(coeff)**2*float(np.vdot(amp,amp).real)
            incoherent_hist[gamma,y]=total/(r//b)
    return dict(frozen_sector_tv=float(.5*np.abs(frozen-actual).sum()),
                incoherent_sector_tv=float(.5*np.abs(incoherent-actual).sum()),
                deleted_history_cross_terms_tv=float(.5*np.abs(incoherent_hist.sum(axis=0)-actual).sum()),
                actual_mass=float(actual.sum()), frozen_mass=float(frozen.sum()),
                incoherent_sector_mass=float(incoherent.sum()),
                deleted_history_mass=float(incoherent_hist.sum()))


def edge_audit():
    rows=[]
    cases=[("k0",fixture(route_count=0)),
           ("W0",CoherentReflectionCircuit(10,2,4,{0:rx(np.pi/2)},{})),
           ("K0",CoherentReflectionCircuit(10,2,4,{}, {0:(0,np.pi/4)})),
           ("end",CoherentReflectionCircuit(10,2,4,{}, {4:(0,np.pi/4)})),
           ("selfloop",CoherentReflectionCircuit(10,2,4,{}, {2:(0,np.pi/4)}))]
    for name,c in cases:
        audited = audit_one(c)
        rows.append(dict(case=name, **audited,
                         prefix_evaluations=c.sample(np.random.default_rng(19))["prefix_vector_evaluations"]))
    return rows


def wide_audit():
    rows=[]
    for k in (2,3):
        defects={1:rx(np.pi/2),3:rx(np.pi/2),4:rx(np.pi/2)}
        refs={2:(0,np.pi/4),3:(1,np.pi/4)}
        if k==3: refs[4]=(0,np.pi/4)
        c=CoherentReflectionCircuit(2*1000000007,2,63,defects,refs)
        rng=np.random.default_rng(100+k)
        sample_start=time.perf_counter()
        samples=[c.sample(rng) for _ in range(2)]
        sample_seconds=time.perf_counter()-sample_start
        rejection_cap=10000
        rejection_start=time.perf_counter()
        rej=[c.sample_rejection(rng,max_proposals=rejection_cap) for _ in range(2)]
        rejection_seconds=time.perf_counter()-rejection_start
        rows.append(dict(k=k, stats=c.stats(), sample_wall_seconds=sample_seconds,
                         rejection_wall_seconds=rejection_seconds,
                         rejection_max_proposals_cap=rejection_cap,
                         sample_prefix_evaluations=[x["prefix_vector_evaluations"] for x in samples],
                         sample_histories=c.stats()["history_count_upper_bound"],
                         rejection_proposals=[x["rejection_proposals"] for x in rej]))
    return rows


def main():
    started=time.time()
    report={"prediction":"P1/P2/P3 with C1/C2/C3 must-fail controls",
            "algorithm2_audit":"BGL gate-by-gate block update: unchanged-coordinate marginal is invariant under each local unitary; this experiment uses the specialized sparse blocks.",
            "fixtures":[], "controls":{}, "edges":[], "wide":[], "status":"PASS"}
    try:
        tiny=fixture()
        audit=audit_one(tiny)
        report["fixtures"].append({"name":"r10_t4_k2",**audit})
        second=fixture(period=14,width=4,route_count=2)
        audit2=audit_one(second)
        report["fixtures"].append({"name":"r14_t4_k2",**audit2})
        report["controls"]=controls(tiny,direct_joint(tiny))
        report["edges"]=edge_audit()
        report["wide"]=wide_audit()
        all_fixtures=(audit, audit2)
        p1_error=max(row[k] for row in all_fixtures
                     for k in ("prefix_max_abs_error", "joint_max_abs_error"))
        p2_error=max(row[k] for row in all_fixtures
                     for k in ("transition_max_abs_error", "marginal_vs_sequential"))
        mass_error=max(abs(row[k]-1) for row in all_fixtures
                       for k in ("transition_mass", "direct_mass", "sequential_mass"))
        edge_ok=all(max(row[key] for key in
                        ("prefix_max_abs_error", "joint_max_abs_error",
                         "transition_max_abs_error", "marginal_vs_sequential")) < 5e-11
                    and max(abs(row[key]-1) for key in
                            ("transition_mass", "direct_mass", "sequential_mass")) < 5e-11
                    for row in report["edges"])
        control_masses = [report["controls"][key] for key in
                          ("actual_mass", "frozen_mass", "incoherent_sector_mass",
                           "deleted_history_mass")]
        controls_normalize = all(abs(m-1) < 5e-11 for m in control_masses)
        wide_ok=all(row["sample_histories"] == 1 << row["k"]
                    and row["stats"]["orbit_table_entries"] == 0
                    and row["stats"]["sector_table_entries"] == 0
                    and row["stats"]["output_table_entries"] == 0
                    and all(x > 0 for x in row["sample_prefix_evaluations"])
                    for row in report["wide"])
        exp.check("P1", p1_error < 5e-11,
                  f"tiny prefix/joint max error {p1_error:.3g}")
        exp.check("P2", p2_error < 5e-11 and mass_error < 5e-11,
                  f"all tiny transition/sequential max error {p2_error:.3g}, mass error {mass_error:.3g}")
        exp.check("P3", edge_ok and wide_ok and controls_normalize,
                  f"edge_ok={edge_ok}, wide_table_free={wide_ok}, control_masses={control_masses}")
        # Both controls are required to fail, and route interference must be
        # visible; a vacuous control is recorded as failure of the experiment.
        exp.fail_check("C1", report["controls"]["frozen_sector_tv"] > 1e-5,
                       f"frozen-sector TV {report['controls']['frozen_sector_tv']:.6g}")
        exp.fail_check("C2", report["controls"]["incoherent_sector_tv"] > 1e-5,
                       f"incoherent-sector TV {report['controls']['incoherent_sector_tv']:.6g}")
        exp.fail_check("C3", report["controls"]["deleted_history_cross_terms_tv"] > 1e-5,
                       f"deleted-history TV {report['controls']['deleted_history_cross_terms_tv']:.6g}")
        if not exp.finish():
            raise AssertionError("lab.Experiment harness reported a failed check")
    except Exception as exc:
        report["status"]="FAIL"
        report["error"]={"type":type(exc).__name__,"message":str(exc),"traceback":traceback.format_exc()}
    report["elapsed_seconds"]=time.time()-started
    path=stamp("coherent_route_sampling")
    path.write_text(json.dumps(report,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"status":report["status"],"report":str(path),
                      "elapsed_seconds":report["elapsed_seconds"]},sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
