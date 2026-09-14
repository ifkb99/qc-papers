"""Audit deterministic coarse-sector routing by orbit-character phases.

DERIVED BEFORE MEASUREMENT.  For r=b*M, let W be repeated on every
consecutive b-label block and let

  D_q |b*m+p> = exp(+2*pi*i*q*m/M) |b*m+p>.

Since D_q maps |alpha,p> to |alpha-q mod M,p>, the gate D_q W routes a
coarse sector bijectively while applying the same b-by-b W.  Several such
gates therefore route an initially uniform alpha deterministically; the final
work trace still removes cross-sector coherences.  The route changes the
sector used by later T_alpha^(2^i), so this must be checked in original
ascending control order.

P1: explicit projections of D_q and repeated W map source sectors to the
    predicted target/sign and preserve the supplied b-block.
P2: an independent full-r sequential_path reference agrees with both an
    independently assembled routed sequential contraction and RoutedOrbitCircuit
    for multiple noncommuting blocks and q values.
P3: q=0, q divisible by M, and negative q modulo M give the expected equal
    probabilities; a small r=4 Circuit/statevec fixture agrees with full-r.
P4: width-63 known-r routed samples execute with bounded finite-work payload.

C1: wrong target/update sign must fail at the matrix level.  A frozen-label
    probability comparison is diagnostic only; its uniform mixture can be
    vacuous through relabeling (M=3 avoids that matrix-level ambiguity).
C2: freezing alpha despite a nonzero D_q must fail at probability level.
C3: replacing an integer route by a fractional character phase must fail at
    probability level; coherent sector mixing is not deterministic routing.

This is standard block-permutation/symmetry routing under supplied known order
and indexed labels, not a factoring result or a physical single-qubit claim.
All dense diagnostics are capped before allocation at 16 MiB, r<=18 and
width<=7 except the width-63 finite-work sample.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
     'numpy<2.5' python -m experiments.experiment_sector_routing
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from circuits import Circuit
from lab import Experiment
from lab.periodic import RoutedOrbitCircuit
from lab.semiclassical import sequential_path
import statevec


MAX_BYTES = 16 * 1024 * 1024
MAX_R = 18
MAX_T = 7
TOL = 3e-10


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def report_path(prefix: str = "sector_routing") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    return root / f"{prefix}_{stamp}.json"


def dft_unitary(b: int, phase: float = 0.0) -> np.ndarray:
    p = np.arange(b)
    f = np.exp(2j * np.pi * p[:, None] * p[None, :] / b) / np.sqrt(b)
    return f @ np.diag(np.exp(1j * phase * np.arange(b)))


def sector_basis(r: int, b: int, alpha: int) -> np.ndarray:
    if r > MAX_R or r % b or not 0 <= alpha < r // b:
        raise ValueError("invalid sector dimensions")
    M = r // b
    guard((r, b), np.complex128, "sector basis")
    q = np.zeros((r, b), dtype=complex)
    for m in range(M):
        q[b * m:b * m + b, :] = np.eye(b) * np.exp(
            -2j * np.pi * alpha * m / M) / np.sqrt(M)
    return q


def periodic_block(r: int, b: int, block: np.ndarray) -> np.ndarray:
    if not 1 <= r <= MAX_R or r % b:
        raise ValueError("invalid periodic dimensions")
    guard((r, r), np.complex128, "periodic block")
    result = np.zeros((r, r), dtype=complex)
    for m in range(r // b):
        result[m * b:(m + 1) * b, m * b:(m + 1) * b] = block
    return result


def character_diagonal(r: int, b: int, q: float) -> np.ndarray:
    """D_q; fractional q is deliberate only for the C3 control."""
    if not 1 <= r <= MAX_R or not 1 <= b <= r or r % b:
        raise ValueError("b must divide r")
    # Guard the actual resulting matrix before arange/diag temporaries.
    guard((r, r), np.complex128, "character diagonal")
    M = r // b
    phases = np.exp(2j * np.pi * q * (np.arange(r) // b) / M)
    return np.diag(phases)


def shift_matrix(r: int, power: int) -> np.ndarray:
    guard((r, r), np.complex128, "orbit shift")
    result = np.zeros((r, r), dtype=complex)
    for j in range(r):
        result[(j + power) % r, j] = 1.0
    return result


def sector_shift_matrix(b: int, M: int, current: int, power: int) -> np.ndarray:
    """T_current^power with integer phase reduction before trig conversion."""
    if not 1 <= b <= 64 or not 1 <= M <= np.iinfo(np.int64).max:
        raise ValueError("invalid sector shift dimensions")
    guard((b, b), np.complex128, "sector shift")
    result = np.zeros((b, b), dtype=complex)
    for col in range(b):
        wraps, target = divmod(col + int(power), b)
        phase_index = (int(current) * wraps) % M
        result[target, col] = np.exp(2j * np.pi * phase_index / M)
    return result


def full_distribution(r: int, b: int, width: int,
                      defects: dict[int, np.ndarray],
                      shifts: dict[int, float]) -> np.ndarray:
    """Independent full-r original-order sequential_path marginal."""
    if width > MAX_T or r > MAX_R:
        raise ValueError("full reference cap exceeded")
    guard((r, r), np.complex128, "full reference identity")
    identity = np.eye(r, dtype=complex)
    blocks = {s: periodic_block(r, b, w) for s, w in defects.items()}
    initial = np.zeros(r, dtype=complex)
    initial[0] = 1.0
    if 0 in blocks or 0 in shifts:
        initial = (character_diagonal(r, b, shifts.get(0, 0.0))
                   @ blocks.get(0, identity) @ initial)
    pairs = []
    for i in range(width):
        gate = identity
        insertion = i + 1
        if insertion in blocks or insertion in shifts:
            gate = (character_diagonal(r, b, shifts.get(insertion, 0.0))
                    @ blocks.get(insertion, identity))
        pairs.append((gate, gate @ shift_matrix(r, 1 << i)))
    guard((1 << width,), np.float64, "full marginal")
    return np.array([
        sequential_path(pairs, initial, output=y)["conditional_path_probability"]
        for y in range(1 << width)
    ])


def local_routed_distribution(r: int, b: int, width: int,
                              defects: dict[int, np.ndarray],
                              shifts: dict[int, int], *,
                              initial_sector: int | None = None,
                              update_sign: int = -1,
                              frozen: bool = False) -> np.ndarray:
    """Independent b-dimensional route assembled around sequential_path."""
    if not 1 <= r <= MAX_R or not 1 <= b <= r or r % b or not 0 <= width <= MAX_T:
        raise ValueError("local route exceeds bounded diagnostic dimensions")
    M = r // b
    guard((b, b), np.complex128, "local route identity")
    guard((1 << width,), np.float64, "local route marginal")
    identity = np.eye(b, dtype=complex)
    sectors = tuple(range(M)) if initial_sector is None else (int(initial_sector),)
    if not 0 <= sectors[0] < M:
        raise ValueError("initial sector outside range")
    marginal = np.zeros(1 << width, dtype=float)
    for alpha0 in sectors:
        current = (alpha0 + update_sign * shifts.get(0, 0)) % M
        initial = defects.get(0, identity)[:, 0]
        pairs = []
        for i in range(width):
            insertion = i + 1
            W = defects.get(insertion, identity)
            p = sector_shift_matrix(b, M, current, 1 << i)
            pairs.append((W, W @ p))
            if not frozen:
                current = (current + update_sign * shifts.get(insertion, 0)) % M
        marginal += np.array([
            sequential_path(pairs, initial, output=y)[
                "conditional_path_probability"]
            for y in range(1 << width)
        ])
    return marginal / len(sectors)


def local_routed_sample(r: int, b: int, width: int,
                        defects: dict[int, np.ndarray], shifts: dict[int, int],
                        rng: np.random.Generator) -> dict:
    if not 1 <= b <= 64 or not 1 <= r or r % b or not 0 <= width <= 63:
        raise ValueError("invalid routed sample dimensions")
    M = r // b
    guard((b, b), np.complex128, "routed sample identity")
    alpha0 = int(rng.integers(M))
    current = (alpha0 - shifts.get(0, 0)) % M
    identity = np.eye(b, dtype=complex)
    initial = defects.get(0, identity)[:, 0]
    pairs = []
    for i in range(width):
        insertion = i + 1
        W = defects.get(insertion, identity)
        # Integer-reduce current*wraps modulo M before trig conversion, as in
        # production's shift() path; this remains valid for width 63.
        p = sector_shift_matrix(b, M, current, 1 << i)
        pairs.append((W, W @ p))
        current = (current - shifts.get(insertion, 0)) % M
    result = sequential_path(pairs, initial, rng=rng)
    return dict(initial_sector=alpha0, final_sector=current, **result)


def abstract_r4_circuit(theta0: float, theta1: float) -> np.ndarray:
    """r=4,b=2 D_1 fixture using existing Circuit/statevec gates."""
    qc = Circuit(4)
    qc.h(2).h(3)
    qc.toffoli(2, 0, 1).cnot(2, 0)  # controlled U
    qc.rz(0, theta0)                 # repeated W after one control
    qc.rz(1, np.pi)                  # D_1: cell-character phase, up to global
    qc.cnot(3, 1)                    # controlled U^2
    qc.rx(0, theta1)                 # second noncommuting W
    qc.qft([2, 3], inverse=True)
    psi = statevec.run(qc)
    probs = np.zeros(4, dtype=float)
    for index, amplitude in enumerate(psi):
        probs[index >> 2] += abs(amplitude) ** 2
    return probs


def main() -> None:
    exp = Experiment("sector_routing", doc=__doc__)
    exp.predict("P1", "D_q routes alpha to alpha-q and W preserves the b-state")
    exp.predict("P2", "full-r and routed sequential marginals agree for multiple q/W")
    exp.predict("P3", "identity/modulo/negative routes and r=4 Circuit agree")
    exp.predict("P4", "width-63 known-r routed sample stays finite-work")
    exp.must_fail("C1", "wrong routing sign/target fails the explicit sector matrix map")
    exp.must_fail("C2", "freezing alpha despite D_q fails probabilities")
    exp.must_fail("C3", "fractional character phase is not deterministic q routing")

    started = time.perf_counter()
    rows: list[dict] = []
    max_projection_error = 0.0
    max_full_route_error = 0.0
    max_production_error = 0.0

    # P1: source/target projections for b=2,3 and M=2,3,4, including q=0,
    # negative q and multiples of M.  W is independently varied by b only.
    exp.section("P1 source-to-target projections")
    for b, block in ((2, dft_unitary(2, 0.23)), (3, dft_unitary(3, 0.41))):
        for M in (2, 3, 4):
            r = b * M
            Wfull = periodic_block(r, b, block)
            for q in (0, 1, -1, M, 2 * M + 1):
                qmod = q % M
                for alpha in range(M):
                    source = sector_basis(r, b, alpha)
                    target = sector_basis(r, b, (alpha - q) % M)
                    d_err = float(np.max(np.abs(
                        target.conj().T @ character_diagonal(r, b, q) @ source
                        - np.eye(b))))
                    w_err = float(np.max(np.abs(
                        source.conj().T @ Wfull @ source - block)))
                    max_projection_error = max(max_projection_error, d_err, w_err)
                    exp.check("P1", d_err < 3e-11 and w_err < 3e-11,
                              f"b={b},M={M},q={q},alpha={alpha}: "
                              f"D={d_err:.2e}, W={w_err:.2e}")
                rows.append(dict(series="projection", b=b, M=M, r=r, q=q,
                                 q_mod=qmod, max_D_error=max(
                                     float(np.max(np.abs(
                                         sector_basis(r, b, (a - q) % M).conj().T
                                         @ character_diagonal(r, b, q)
                                         @ sector_basis(r, b, a) - np.eye(b))))
                                     for a in range(M)),
                                 block_error=w_err))

    # P2: multiple noncommuting W_s and several routed charges.  The full-r
    # path is independent of the sector decomposition; production is audited
    # against the separately assembled local route as well.
    exp.section("P2 routed/full original-order marginals")
    for b, M in ((2, 3), (3, 2), (3, 3)):
        r, width = b * M, 4
        w0 = dft_unitary(b, 0.17)
        w1 = dft_unitary(b, 0.53)
        w2 = dft_unitary(b, -0.31)
        defects = {1: w0, 2: w1, 3: w2}
        shifts = {0: 1, 2: -1, 3: 2 * M + 1}
        full = full_distribution(r, b, width, defects, shifts)
        local = local_routed_distribution(r, b, width, defects,
                                          {s: int(q) for s, q in shifts.items()})
        production = RoutedOrbitCircuit(r, b, width, defects, shifts)
        prod = np.array([sum(production.forced_joint(alpha, y)[
            "conditional_path_probability"] / M for alpha in range(M))
                          for y in range(1 << width)])
        local_error = float(np.max(np.abs(full - local)))
        production_error = float(np.max(np.abs(full - prod)))
        max_full_route_error = max(max_full_route_error, local_error)
        max_production_error = max(max_production_error, production_error)
        exp.check("P2", local_error < 3e-10 and production_error < 3e-10
                  and abs(full.sum() - 1) < 3e-10,
                  f"b={b},M={M},r={r}: local={local_error:.2e}, "
                  f"production={production_error:.2e}, norm={full.sum():.12g}")
        rows.append(dict(series="multiple_routes", b=b, M=M, r=r, width=width,
                         shifts=shifts, full=full.tolist(), local=local.tolist(),
                         production=prod.tolist(), local_error=local_error,
                         production_error=production_error))

    # P3 identity, q modulo M, and negative q controls at probability level.
    exp.section("P3 modulo/negative routes and Circuit fixture")
    b, M, r, width = 2, 3, 6, 4
    defects = {1: dft_unitary(2, 0.27), 2: dft_unitary(2, -0.44)}
    no_route = local_routed_distribution(r, b, width, defects, {})
    q0 = local_routed_distribution(r, b, width, defects, {1: 0, 2: 0})
    qM = local_routed_distribution(r, b, width, defects, {1: M, 2: -M})
    qneg = local_routed_distribution(r, b, width, defects, {1: -1})
    qmod = local_routed_distribution(r, b, width, defects, {1: M - 1})
    exp.check("P3", max(np.max(np.abs(no_route - q0)),
                         np.max(np.abs(no_route - qM)),
                         np.max(np.abs(qneg - qmod))) < 3e-10,
              f"q0/qM/negative-mod errors={max(np.max(np.abs(no_route-q0)), np.max(np.abs(no_route-qM)), np.max(np.abs(qneg-qmod))):.2e}")
    theta0, theta1 = np.pi / 3, -np.pi / 5
    b, M, r, width = 2, 2, 4, 2
    circuit_defects = {1: np.diag([np.exp(-0.5j * theta0), np.exp(0.5j * theta0)]),
                       2: np.array([[np.cos(theta1 / 2), -1j * np.sin(theta1 / 2)],
                                    [-1j * np.sin(theta1 / 2), np.cos(theta1 / 2)]])}
    circuit_full = full_distribution(r, b, width, circuit_defects, {1: 1})
    circuit_gate = abstract_r4_circuit(theta0, theta1)
    circuit_error = float(np.max(np.abs(circuit_full - circuit_gate)))
    exp.check("P3", circuit_error < 3e-10,
              f"r=4 Circuit/statevec vs full-r={circuit_error:.2e}")
    rows.append(dict(series="modulo_and_circuit", q0_error=float(np.max(abs(no_route-q0))),
                     qM_error=float(np.max(abs(no_route-qM))),
                     negative_mod_error=float(np.max(abs(qneg-qmod))),
                     circuit_error=circuit_error))

    # P4: wide known-r supplied route.  No output table is allocated; only a
    # finite-work sample is drawn from both independent local and production
    # implementations, and payload metadata is retained.
    exp.section("P4 width-63 supplied-period sample")
    b, M, r, width = 3, 1_000_000_007, 3_000_000_021, 63
    wide_defects = {16: dft_unitary(3, 0.19), 32: dft_unitary(3, -0.37),
                    48: dft_unitary(3, 0.61)}
    wide_shifts = {0: -1, 16: 1, 32: 2, 48: M + 1}
    # The production class accepts int64-range known periods; no dense r-array.
    production_wide = RoutedOrbitCircuit(r, b, width, wide_defects, wide_shifts)
    local_rng = np.random.default_rng(20260910)
    prod_rng = np.random.default_rng(20260910)
    local_sample = local_routed_sample(r, b, width, wide_defects,
                                       {s: int(q) for s, q in wide_shifts.items()},
                                       local_rng)
    prod_sample = production_wide.sample(prod_rng)
    stats = production_wide.stats()
    exp.check("P4", 0 <= local_sample["output"] < 1 << width
              and 0 <= prod_sample["output"] < 1 << width
              and stats["sampled_work_dimension"] == b
              and stats["orbit_table_entries"] == 0
              and stats["output_table_entries"] == 0,
              f"outputs in range; local/prod={local_sample['output']}/{prod_sample['output']}, "
              f"payload={prod_sample['conservative_payload_estimate_bytes']}")
    rows.append(dict(series="wide_sample", period=r, block_size=b, width=width,
                     local_sample=local_sample, production_sample=prod_sample,
                     stats=stats))

    # Must-fail controls.  Matrix sign is explicit because a uniform alpha
    # mixture can hide sign relabeling.  Probability controls use a fixed
    # nonuniform alpha=0 route and a genuinely fractional character phase.
    exp.section("must-fail controls")
    r, b, M, width = 6, 2, 3, 3
    q = 1
    source = sector_basis(r, b, 0)
    wrong_target = sector_basis(r, b, (0 + q) % M)
    wrong_matrix_error = float(np.max(np.abs(
        wrong_target.conj().T @ character_diagonal(r, b, q) @ source - np.eye(b))))
    defects = {1: dft_unitary(2, 0.21), 2: dft_unitary(2, -0.39)}
    shifts = {1: 1, 2: -1}
    correct_frozen_case = local_routed_distribution(
        r, b, width, defects, shifts, initial_sector=0)
    wrong_update_case = local_routed_distribution(
        r, b, width, defects, shifts, initial_sector=0, update_sign=1)
    frozen_case = local_routed_distribution(
        r, b, width, defects, shifts, initial_sector=0, frozen=True)
    wrong_probability_error = float(np.max(np.abs(correct_frozen_case - wrong_update_case)))
    frozen_probability_error = float(np.max(np.abs(correct_frozen_case - frozen_case)))
    exp.fail_check("C1", wrong_matrix_error > 1e-4 or wrong_probability_error > 1e-4,
                   f"wrong target matrix={wrong_matrix_error:.6g}, "
                   f"wrong update probability={wrong_probability_error:.6g}")
    exp.fail_check("C2", frozen_probability_error > 1e-4,
                   f"frozen-alpha probability error={frozen_probability_error:.6g}")
    # Put one route at s=0 and one at s=1.  A fractional character at the
    # earlier insertion then remains observable; the previous s=1/s=2 choice
    # happened to cancel in this small output row.
    fractional_defects = {0: dft_unitary(2, 0.21), 1: dft_unitary(2, -0.39)}
    fractional = full_distribution(6, 2, 2, fractional_defects,
                                   {1: 0.37, 0: -1.0})
    integer_route = full_distribution(6, 2, 2, fractional_defects,
                                      {1: 1, 0: -1})
    fractional_error = float(np.max(np.abs(fractional - integer_route)))
    exp.fail_check("C3", fractional_error > 1e-4,
                   f"fractional-vs-integer route probability error={fractional_error:.6g}")
    rows.append(dict(series="controls", wrong_matrix_error=wrong_matrix_error,
                     wrong_probability_error=wrong_probability_error,
                     frozen_probability_error=frozen_probability_error,
                     fractional_error=fractional_error))

    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(
                   max_dense_bytes=MAX_BYTES, max_r=MAX_R, max_t=MAX_T,
                   max_projection_error=max_projection_error,
                   max_full_route_error=max_full_route_error,
                   max_production_error=max_production_error,
                   references="full-r sequential_path, independently assembled routed sequential_path, RoutedOrbitCircuit, and r=4 Circuit/statevec",
                   assumption="known order, indexed orbit labels, exact integer character charges; D_q physical construction is an input cost",
                   routing_note="uniform initial alpha remains valid because routes are deterministic bijections; individual alpha is not conserved",
                   wrong_sign_probability_control="vacuous at fixed alpha in this row; C1 is matrix-level",
                   elapsed_seconds=time.perf_counter() - started,
               ))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("sector_routing_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
