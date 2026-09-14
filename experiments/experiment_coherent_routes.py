"""Charged coherent-route formula for a few sector reflections.

Use a known indexed orbit with b=2, M=5, r=10 and t=5.  Repeated Rx(pi/2)
blocks act at s=1,3,4.  At s=2 insert exp(-i theta J_0/2), and at s=3
insert exp(-i theta J_1/2) AFTER its repeated background block, where

    R_cell |b*m+p> = |b*(-m mod M)+p>,
    D_q |b*m+p> = exp(+2*pi*i*q*m/M)|b*m+p>,
    J_q = D_q R_cell.

Each J_q is Hermitian/unitary and routes coarse alpha to -alpha-q.  Expanding
the two inserted gates gives histories h with coefficient c_h.  For a fixed
FINAL sector gamma, history h starts from alpha_h=Pi_h^{-1}(gamma), not from
one common initial alpha.  The complete-output formula tested here is

    p(gamma,y) = ||sum_h c_h a_h(gamma,y)||^2 / M,
    q(gamma,y) = sum_h |c_h| ||a_h(gamma,y)||^2 / (M*B),
    B = sum_h |c_h|.

All a_h retain their original Fourier-filtered complex phases.  Cauchy-Schwarz
predicts p<=B^2 q, normalized accepted law p, and exact mean proposal count
B^2.  This is a fixed-output formula audit, not a production sampler.

PREDICTIONS, WRITTEN BEFORE MEASUREMENT.

  P1  The final-gamma coherent-history formula agrees with the existing full-r
      sequential_path output law for theta=0, pi/8, pi/4, pi/2 and for the
      separate k=0,1,2 route-count rows.
  P2  J_q is Hermitian/unitary and its projected routing is exactly
      alpha -> -alpha-q; all history output laws normalize.
  P3  The B^2 q envelope holds, the normalized accepted law is exact, and
      the measured mean proposal count is B^2.  Final gamma weights need not be
      uniform.

  C1  Deleting coherent cross terms must fail.
  C2  Starting every history from alpha=0 instead of alpha_h must fail.
  C3  Treating uniform final gamma as the target marginal must fail on a
      nontrivial route row; if it is vacuous, preserve that fact explicitly.

All dense allocations are capped before allocation (r<=18,t<=6,k<=3,16 MiB).
The finite history sum charges H and M; no second generic propagator or
large-sector/output array is introduced.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_coherent_routes
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


MAX_DENSE_BYTES = 16 * 1024 * 1024
BLOCK, M, PERIOD, WIDTH = 2, 5, 10, 5
THETA_SWEEP = np.pi * np.array([-1/4, 0., 1/8, 1/4, 1/2])
BACKGROUND_THETA = np.pi / 2
ROUTE_SCHEDULE = ((2, 0), (3, 1), (4, 0))


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def report_path(prefix: str = "coherent_routes_formula") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def rx(theta: float) -> np.ndarray:
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return np.array([[c, -1j*s], [-1j*s, c]], dtype=complex)


def background_block(theta: float = BACKGROUND_THETA) -> np.ndarray:
    return rx(theta)


def repeated_matrix(block: np.ndarray) -> np.ndarray:
    guard((PERIOD, PERIOD), np.complex128, "repeated orbit block")
    result = np.kron(np.eye(M, dtype=complex), np.asarray(block, dtype=complex))
    return result


def cell_reflection() -> np.ndarray:
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for m in range(M):
        for p in range(BLOCK):
            result[BLOCK*((-m) % M) + p, BLOCK*m + p] = 1.
    return result


def route_gate(q: int) -> np.ndarray:
    if not isinstance(q, int) or not 0 <= q < M:
        raise ValueError("q outside cell character range")
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for m in range(M):
        for p in range(BLOCK):
            target = BLOCK*((-m) % M) + p
            # D_q is evaluated on the reflected target cell -m.
            result[target, BLOCK*m+p] = unit_phase(-q*m, M)
    return result


def sector_basis(alpha: int) -> np.ndarray:
    """Columns |alpha,p> in the indexed orbit basis."""
    if not 0 <= alpha < M:
        raise ValueError("alpha outside sector range")
    result = np.zeros((PERIOD, BLOCK), dtype=complex)
    for m in range(M):
        for p in range(BLOCK):
            result[BLOCK*m+p, p] = unit_phase(-alpha*m, M) / np.sqrt(M)
    return result


def route_alpha(alpha: int, q: int) -> int:
    return (-alpha - q) % M


def shift(power: int) -> np.ndarray:
    result = np.zeros((PERIOD, PERIOD), dtype=complex)
    for j in range(PERIOD):
        result[(j + power) % PERIOD, j] = 1.
    return result


def history_coefficients(theta: float, route_count: int, schedule=ROUTE_SCHEDULE) -> list[tuple[tuple[int, ...], complex]]:
    if not 0 <= route_count <= len(schedule):
        raise ValueError("route count exceeds bounded experiment cap")
    c, s = np.cos(theta / 2), -1j*np.sin(theta / 2)
    result = []
    for choices in itertools.product((0, 1), repeat=route_count):
        coeff = 1.+0j
        for bit in choices:
            coeff *= c if bit == 0 else s
        result.append((choices, coeff))
    return result


def route_insertions(route_count: int, schedule=ROUTE_SCHEDULE) -> tuple[tuple[int, int], ...]:
    if not 0 <= route_count <= len(schedule):
        raise ValueError("route count exceeds supplied route schedule")
    if route_count == 0:
        return ()
    return tuple(schedule[:route_count])


def history_matrices(theta: float, choices: tuple[int, ...], route_count: int,
                     schedule=ROUTE_SCHEDULE) -> dict[int, np.ndarray]:
    """Full-r work matrices after each control for one normalized history."""
    identity = np.eye(PERIOD, dtype=complex)
    result = {1: repeated_matrix(background_block()),
              3: repeated_matrix(background_block()),
              4: repeated_matrix(background_block())}
    for (s, q), bit in zip(route_insertions(route_count, schedule), choices):
        if bit and 0 < s < WIDTH:
            result[s] = route_gate(q) @ result.get(s, identity)
    return result


def history_route(alpha: int, choices: tuple[int, ...], route_count: int,
                  schedule=ROUTE_SCHEDULE) -> int:
    current = alpha
    for (s, q), bit in zip(route_insertions(route_count, schedule), choices):
        if bit:
            current = route_alpha(current, q)
    return current


def inverse_initial_sector(gamma: int, choices: tuple[int, ...], route_count: int,
                           schedule=ROUTE_SCHEDULE) -> int:
    candidates = [alpha for alpha in range(M)
                  if history_route(alpha, choices, route_count, schedule) == gamma]
    if len(candidates) != 1:
        raise AssertionError("history route must be a sector permutation")
    return candidates[0]


def history_amplitude(theta: float, choices: tuple[int, ...], route_count: int,
                      gamma: int, output: int, *, force_initial: int | None = None,
                      schedule=ROUTE_SCHEDULE) -> np.ndarray:
    alpha = (inverse_initial_sector(gamma, choices, route_count, schedule)
             if force_initial is None else force_initial)
    sector = sector_basis(gamma)
    state = sector_basis(alpha)[:, 0]
    matrices = history_matrices(theta, choices, route_count, schedule)
    selected = {s: q for (s, q), bit in zip(route_insertions(route_count, schedule), choices)
                if bit}
    if 0 in selected:
        state = route_gate(selected[0]) @ state
    accumulator = np.zeros(BLOCK, dtype=complex)
    for e in range(1 << WIDTH):
        work = state.copy()
        for i in range(WIDTH):
            if (e >> i) & 1:
                work = shift(1 << i) @ work
            w = matrices.get(i + 1)
            if w is not None:
                work = w @ work
        if WIDTH in selected:
            work = route_gate(selected[WIDTH]) @ work
        accumulator += np.exp(-2j*np.pi*output*e/(1 << WIDTH)) * (sector.conj().T @ work)
    # Initial |+> amplitude and inverse-QFT kernel each contribute 1/sqrt(Q).
    return accumulator / (1 << WIDTH)


def sequential_output(theta: float, route_count: int, *, initial_alpha: int | None = None,
                      choices_override: tuple[int, ...] | None = None) -> np.ndarray:
    """Existing sequential_path output law for the coherent inserted gates."""
    identity = np.eye(PERIOD, dtype=complex)
    # Actual gate at each insertion is W_s times exp(-i theta J_q/2) after it.
    matrices = {1: repeated_matrix(background_block()),
                3: repeated_matrix(background_block()),
                4: repeated_matrix(background_block())}
    for s, q in ROUTE_SCHEDULE[:route_count]:
        matrices[s] = (np.cos(theta/2)*identity - 1j*np.sin(theta/2)*route_gate(q)) @ matrices.get(s, identity)
    pairs = []
    for i in range(WIDTH):
        w = matrices.get(i + 1, identity)
        pairs.append((w, w @ shift(1 << i)))
    initial = np.zeros(PERIOD, dtype=complex)
    if initial_alpha is None:
        initial[0] = 1.
    else:
        initial = sector_basis(initial_alpha)[:, 0]
    result = np.zeros(1 << WIDTH, dtype=float)
    for y in range(1 << WIDTH):
        result[y] = sequential_path(pairs, initial, output=y)["conditional_path_probability"]
    return result


def coherent_formula(theta: float, route_count: int, *, wrong_initial: bool = False,
                     delete_cross_terms: bool = False,
                     schedule=ROUTE_SCHEDULE) -> dict:
    histories = history_coefficients(theta, route_count, schedule)
    qsize = 1 << WIDTH
    if PERIOD > 18 or WIDTH > 6 or route_count > len(schedule):
        raise ValueError("formula row exceeds bounded r/t/k scope")
    guard((M, qsize), np.float64, "sector-output formula table")
    p_joint = np.zeros((M, qsize), dtype=float)
    q_joint = np.zeros((M, qsize), dtype=float)
    incoherent = np.zeros((M, qsize), dtype=float)
    bnorm = float(sum(abs(coeff) for _, coeff in histories))
    max_envelope = 0.
    for gamma in range(M):
        for y in range(qsize):
            coherent = np.zeros(BLOCK, dtype=complex)
            proposal_norm = 0.
            incoherent_mass = 0.
            for choices, coeff in histories:
                amp = history_amplitude(theta, choices, route_count, gamma, y,
                                       force_initial=0 if wrong_initial else None,
                                       schedule=schedule)
                coherent += coeff * amp
                proposal_norm += abs(coeff) * float(np.vdot(amp, amp).real)
                incoherent_mass += abs(coeff)**2 * float(np.vdot(amp, amp).real)
            p_joint[gamma, y] = float(np.vdot(coherent, coherent).real) / M
            q_joint[gamma, y] = proposal_norm / (M*bnorm) if bnorm else 0.
            incoherent[gamma, y] = incoherent_mass / M
            if q_joint[gamma, y] > 0:
                max_envelope = max(max_envelope,
                                   p_joint[gamma, y] / (bnorm*bnorm*q_joint[gamma, y]))
    if delete_cross_terms:
        p_joint = incoherent
    accepted = np.zeros_like(p_joint)
    for gamma in range(M):
        for y in range(qsize):
            if q_joint[gamma, y] > 0:
                accepted[gamma, y] = q_joint[gamma, y] * (
                    p_joint[gamma, y] / (bnorm*bnorm*q_joint[gamma, y]))
    return dict(p_joint=p_joint, q_joint=q_joint, accepted=accepted,
                incoherent=incoherent, B=bnorm, history_count=len(histories),
                max_envelope=max_envelope, mean_proposals=bnorm*bnorm,
                target_output=p_joint.sum(axis=0), final_gamma=p_joint.sum(axis=1),
                accepted_total=float(accepted.sum()), proposal_total=float(q_joint.sum()))


def tv(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum(np.abs(np.asarray(a) - np.asarray(b))) / 2)


def production_law(theta: float, route_count: int,
                   schedule=ROUTE_SCHEDULE) -> dict:
    """Stream the new production target/proposal tables for one bounded row."""
    if route_count > len(schedule):
        raise ValueError("route count exceeds production schedule")
    defects = {s: background_block() for s in (1, 3, 4)}
    reflections = {s: (q, float(theta)) for s, q in schedule[:route_count]}
    circuit = CoherentReflectionCircuit(PERIOD, BLOCK, WIDTH, defects, reflections)
    formula = coherent_formula(theta, route_count, schedule=schedule)
    shape = (M, 1 << WIDTH)
    guard(shape, np.float64, "production audit table")
    target = np.zeros(shape, dtype=float)
    proposal = np.zeros(shape, dtype=float)
    acceptance = np.zeros(shape, dtype=float)
    for gamma in range(M):
        for output in range(1 << WIDTH):
            weights = circuit.rejection_weights(gamma, output)
            target[gamma, output] = weights["joint_probability"]
            proposal[gamma, output] = weights["proposal_joint_probability"]
            acceptance[gamma, output] = weights["acceptance_probability"]
    accepted = proposal * acceptance
    b2 = float(circuit.coefficient_l1 ** 2)
    accepted_expected = target / b2
    target_formula = formula["p_joint"]
    proposal_formula = formula["q_joint"]
    acceptance_formula = np.divide(
        target_formula, b2 * proposal_formula,
        out=np.zeros_like(target_formula), where=proposal_formula > 0)
    # Compare the production acceptance on its exact positive support.  The
    # independent formula retains roundoff-only nonzero amplitudes at outputs
    # that production arithmetic cancels exactly; those are reported below,
    # never thresholded or deleted from any law.
    acceptance_consistency = np.divide(
        target, b2 * proposal, out=np.zeros_like(target), where=proposal > 0)
    support = proposal > 0
    support_error = (float(np.max(np.abs(acceptance[support]-acceptance_consistency[support])))
                     if np.any(support) else 0.)
    return dict(circuit=circuit, formula=formula, target=target, proposal=proposal,
                acceptance=acceptance, accepted=accepted, b2=b2,
                accepted_expected=accepted_expected,
                target_error=float(np.max(np.abs(target-target_formula))),
                proposal_error=float(np.max(np.abs(proposal-proposal_formula))),
                acceptance_error=float(np.max(np.abs(acceptance-acceptance_formula))),
                acceptance_support_error=support_error,
                accepted_error=float(np.max(np.abs(accepted-accepted_expected))),
                target_mass=float(target.sum()), proposal_mass=float(proposal.sum()),
                accepted_mass=float(accepted.sum()),
                target_normalized_error=tv(target/target.sum(), target_formula/target_formula.sum()),
                accepted_normalized_error=tv(accepted/accepted.sum(), target/target.sum()),
                stats=circuit.stats())


def normalized_control_tv(control: dict, reference: np.ndarray) -> tuple[float, float]:
    mass = float(control["p_joint"].sum())
    if not math.isfinite(mass) or mass <= 0:
        return float("inf"), mass
    return tv(control["target_output"] / mass, reference / reference.sum()), mass


def main() -> None:
    exp = Experiment("coherent_routes_formula", doc=__doc__)
    exp.predict("P1", "coherent final-gamma histories match sequential_path output laws")
    exp.predict("P2", "route gates are Hermitian/unitary and routed histories normalize")
    exp.predict("P3", "B^2 envelope and normalized accepted law hold; gamma need not be uniform")
    exp.must_fail("C1", "delete coherent history cross terms")
    exp.must_fail("C2", "start every history from alpha=0")
    exp.must_fail("C3", "treat uniform final gamma as target marginal")
    start = time.perf_counter()
    rows = []
    exp.check("P2", PERIOD <= 18 and WIDTH <= 6 and len(ROUTE_SCHEDULE) <= 3,
              f"preflight r={PERIOD}, t={WIDTH}, k_cap=3, dense_budget={MAX_DENSE_BYTES}")

    # Verify J structure and projection routing before output laws.
    routing_error = hermitian_error = unitary_error = 0.
    for q in (0, 1):
        j = route_gate(q)
        hermitian_error = max(hermitian_error, float(np.max(np.abs(j - j.conj().T))))
        unitary_error = max(unitary_error, float(np.max(np.abs(j.conj().T @ j - np.eye(PERIOD)))))
        for alpha in range(M):
            expected = sector_basis(route_alpha(alpha, q))
            actual = j @ sector_basis(alpha)
            routing_error = max(routing_error, float(np.max(np.abs(actual - expected))))
    exp.check("P2", hermitian_error < 3e-14 and unitary_error < 3e-14
              and routing_error < 3e-14,
              f"J structure: hermitian={hermitian_error:.2e}, unitary={unitary_error:.2e}, "
              f"routing={routing_error:.2e}")

    max_deleted_tv = max_wrong_tv = 0.
    max_deleted_normalized_tv = max_wrong_normalized_tv = 0.
    control_masses = []
    max_formula_error = max_reference_norm = 0.
    for theta in THETA_SWEEP:
        theta = float(theta)
        formula = coherent_formula(theta, 2)
        reference = sequential_output(theta, 2)
        target = formula["target_output"]
        formula_error = float(np.max(np.abs(target - reference)))
        norm_error = max(abs(float(formula["p_joint"].sum()) - 1),
                         abs(float(reference.sum()) - 1))
        max_formula_error = max(max_formula_error, formula_error)
        max_reference_norm = max(max_reference_norm, norm_error)
        exp.check("P1", formula_error < 3e-9,
                  f"theta/pi={theta/np.pi:g}: formula/sequential={formula_error:.2e}")
        exp.check("P2", norm_error < 3e-9,
                  f"theta/pi={theta/np.pi:g}: norm={norm_error:.2e}")
        accepted_error = float(np.max(np.abs(
            formula["accepted"] / formula["accepted_total"] - formula["p_joint"])))
        envelope_ok = formula["max_envelope"] <= 1 + 3e-12
        exp.check("P3", envelope_ok and accepted_error < 3e-10
                  and abs(formula["proposal_total"] - 1) < 3e-10,
                  f"theta/pi={theta/np.pi:g}: envelope={formula['max_envelope']:.6g}, "
                  f"accepted={accepted_error:.2e}, qmass={formula['proposal_total']:.2e}")
        wrong = coherent_formula(theta, 2, wrong_initial=True)
        wrong_tv = tv(wrong["target_output"], reference)
        wrong_normalized_tv, wrong_mass = normalized_control_tv(wrong, reference)
        deleted = coherent_formula(theta, 2, delete_cross_terms=True)
        deleted_tv = tv(deleted["target_output"], reference)
        deleted_normalized_tv, deleted_mass = normalized_control_tv(deleted, reference)
        max_wrong_tv = max(max_wrong_tv, wrong_tv)
        max_deleted_tv = max(max_deleted_tv, deleted_tv)
        max_wrong_normalized_tv = max(max_wrong_normalized_tv, wrong_normalized_tv)
        max_deleted_normalized_tv = max(max_deleted_normalized_tv, deleted_normalized_tv)
        control_masses.extend((wrong_mass, deleted_mass))
        exp.check("P2", math.isfinite(wrong_mass) and wrong_mass > 0
                  and math.isfinite(deleted_mass) and deleted_mass > 0,
                  f"theta/pi={theta/np.pi:g}: control masses wrong={wrong_mass:.6g}, "
                  f"deleted={deleted_mass:.6g} (wrong-sector control is intentionally "
                  "not normalized before its failure comparison)")
        rows.append(dict(series="theta", theta=theta, theta_over_pi=theta/np.pi,
                         output=target.tolist(), reference=reference.tolist(),
                         final_gamma=formula["final_gamma"].tolist(),
                         wrong_initial=wrong["target_output"].tolist(),
                         deleted_cross=deleted["target_output"].tolist(),
                         formula_error=formula_error, norm_error=norm_error,
                         accepted_error=accepted_error,
                         wrong_initial_raw_half_l1=wrong_tv,
                         wrong_initial_normalized_tv=wrong_normalized_tv,
                         deleted_cross_raw_half_l1=deleted_tv,
                         deleted_cross_normalized_tv=deleted_normalized_tv,
                         B=formula["B"], history_count=formula["history_count"],
                         mean_proposals=formula["mean_proposals"],
                         max_envelope=formula["max_envelope"],
                         proposal_mass=formula["proposal_total"],
                         accepted_mass=formula["accepted_total"]))

    # Audit the production target, proposal and accepted law, including exact
    # zero-angle histories and the two boundary insertion positions.
    production_rows = []
    for theta in (float(-np.pi/4), 0., float(np.pi/4)):
        for route_count in (0, 1, 2, 3):
            audit = production_law(theta, route_count)
            stats = audit["stats"]
            exp.check("P1", audit["target_error"] < 3e-12,
                      f"production theta/pi={theta/np.pi:g}, k={route_count}: "
                      f"target={audit['target_error']:.2e}")
            exp.check("P3", audit["proposal_error"] < 3e-12
                      and audit["acceptance_support_error"] < 3e-12
                      and audit["accepted_error"] < 3e-12
                      and abs(audit["target_mass"]-1) < 3e-12
                      and abs(audit["proposal_mass"]-1) < 3e-12
                      and abs(audit["accepted_mass"]-1/audit["b2"]) < 3e-12
                      and audit["accepted_normalized_error"] < 3e-12
                      and abs(stats["mathematical_mean_rejection_proposals"]-audit["b2"]) < 3e-12,
                      f"production theta/pi={theta/np.pi:g}, k={route_count}: "
                      f"q={audit['proposal_error']:.2e}, acc-support={audit['acceptance_support_error']:.2e}, "
                      f"acc-independent-raw={audit['acceptance_error']:.2e}, "
                      f"accepted={audit['accepted_error']:.2e}, masses="
                      f"({audit['target_mass']:.6g},{audit['proposal_mass']:.6g},"
                      f"{audit['accepted_mass']:.6g}), B2={audit['b2']:.6g}")
            production_rows.append(dict(theta=theta, theta_over_pi=theta/np.pi,
                                        route_count=route_count,
                                        target_error=audit["target_error"],
                                        proposal_error=audit["proposal_error"],
                                        acceptance_error=audit["acceptance_error"],
                                        acceptance_support_error=audit["acceptance_support_error"],
                                        accepted_error=audit["accepted_error"],
                                        target_mass=audit["target_mass"],
                                        proposal_mass=audit["proposal_mass"],
                                        accepted_mass=audit["accepted_mass"],
                                        accepted_normalized_error=audit["accepted_normalized_error"],
                                        B2=audit["b2"], stats=stats))

    for boundary_name, schedule in (("initial", ((0, 0),)),
                                    ("end", ((WIDTH, 0),))):
        for theta in (float(-np.pi/4), 0., float(np.pi/4)):
            audit = production_law(theta, 1, schedule=schedule)
            exp.check("P1", audit["target_error"] < 3e-12,
                      f"production {boundary_name} theta/pi={theta/np.pi:g}: "
                      f"target={audit['target_error']:.2e}")
            exp.check("P3", audit["proposal_error"] < 3e-12
                      and audit["acceptance_support_error"] < 3e-12
                      and audit["accepted_error"] < 3e-12
                      and abs(audit["target_mass"]-1) < 3e-12
                      and abs(audit["proposal_mass"]-1) < 3e-12,
                      f"production {boundary_name} theta/pi={theta/np.pi:g}: "
                      f"q={audit['proposal_error']:.2e}, accepted={audit['accepted_error']:.2e}, "
                      f"masses=({audit['target_mass']:.6g},{audit['proposal_mass']:.6g})")
            production_rows.append(dict(boundary=boundary_name, theta=theta,
                                        theta_over_pi=theta/np.pi,
                                        target_error=audit["target_error"],
                                        proposal_error=audit["proposal_error"],
                                        acceptance_support_error=audit["acceptance_support_error"],
                                        accepted_error=audit["accepted_error"],
                                        target_mass=audit["target_mass"],
                                        proposal_mass=audit["proposal_mass"],
                                        accepted_mass=audit["accepted_mass"],
                                        B2=audit["b2"]))

    # Vary k separately at a fixed nonzero angle; k=3 is the explicit cap.
    for route_count in (0, 1, 2, 3):
        # The third row uses q=0 at s=4 only to exercise the stated k<=3 cap.
        formula = coherent_formula(np.pi/4, route_count)
        reference = sequential_output(np.pi/4, route_count)
        formula_error = float(np.max(np.abs(formula["target_output"] - reference)))
        accepted_error = float(np.max(np.abs(
            formula["accepted"] / formula["accepted_total"] - formula["p_joint"])))
        exp.check("P1", formula_error < 3e-9,
                  f"k={route_count}: formula/sequential={formula_error:.2e}")
        exp.check("P3", formula["max_envelope"] <= 1 + 3e-12 and accepted_error < 3e-10,
                  f"k={route_count}: B={formula['B']:.6g}, H={formula['history_count']}, "
                  f"envelope={formula['max_envelope']:.6g}, accepted={accepted_error:.2e}")
        rows.append(dict(series="k", route_count=route_count,
                         output=formula["target_output"].tolist(), reference=reference.tolist(),
                         final_gamma=formula["final_gamma"].tolist(), formula_error=formula_error,
                         accepted_error=accepted_error, B=formula["B"],
                         history_count=formula["history_count"],
                         mean_proposals=formula["mean_proposals"],
                         max_envelope=formula["max_envelope"]))

    exp.fail_check("C1", max_deleted_normalized_tv > 1e-5,
                   f"max deleted-cross-term raw_half_l1={max_deleted_tv:.6g}, "
                   f"normalized_TV={max_deleted_normalized_tv:.6g}, "
                   f"masses={control_masses[1::2]}")
    exp.fail_check("C2", max_wrong_normalized_tv > 1e-5,
                   f"max wrong-initial-sector raw_half_l1={max_wrong_tv:.6g}, "
                   f"normalized_TV={max_wrong_normalized_tv:.6g}, "
                   f"masses={control_masses[0::2]}")
    representative = coherent_formula(np.pi/4, 2)
    uniform = np.full(M, 1/M)
    uniform_tv = tv(representative["final_gamma"], uniform)
    exp.fail_check("C3", uniform_tv > 1e-5,
                   f"uniform-final-gamma TV={uniform_tv:.6g}")
    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(period=PERIOD, block_size=BLOCK, sectors=M, width=WIDTH,
                             background="Rx(pi/2) repeated at s=1,3,4",
                             route_schedule_initial=((2, 0), (3, 1)),
                             theta_sweep=[float(x) for x in THETA_SWEEP],
                             production_rows=production_rows,
                             control_metrics=dict(
                                 deleted_raw_half_l1=max_deleted_tv,
                                 deleted_normalized_tv=max_deleted_normalized_tv,
                                 wrong_initial_raw_half_l1=max_wrong_tv,
                                 wrong_initial_normalized_tv=max_wrong_normalized_tv,
                                 masses=control_masses),
                             control_mass_range=(min(control_masses), max(control_masses)),
                             prior_failed_artifacts=[
                                 {"artifact": "out/coherent_routes_formula_run.log",
                                  "cause": "independent fixture used the wrong D_q phase sign, "
                                           "applied arithmetic shifts on every exponent bit, "
                                           "and used 1/sqrt(Q) rather than 1/Q normalization"},
                                 {"artifact": "out/coherent_routes_formula_run_2.log",
                                  "cause": "phase/bit-order fixes were present but the inverse-QFT "
                                           "normalization was still 1/sqrt(Q)"},
                                 {"artifact": "out/coherent_routes_formula_run_3.log",
                                  "cause": "all formula checks passed before the preallocation guard "
                                           "was moved ahead of np.kron"}],
                             k_cap=3, dense_budget_bytes=MAX_DENSE_BYTES,
                             structure_errors=dict(hermitian=hermitian_error,
                                                    unitary=unitary_error,
                                                    routing=routing_error),
                             cost_scope="complete-output formula; charge M and H; no sampler",
                             no_amplitude_cutoff=True, elapsed_seconds=elapsed,
                             numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("coherent_routes_formula_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
