"""Orbit-phase structure and coarse-sector routing for physical work phases.

DERIVED BEFORE MEASUREMENT.  For the clean modular orbit, a physical work
Rz on bit q has diagonal phase sign z_q(j)=(-1)**bit_q(a**j mod N).  This
experiment first audits the exact integer orbit and sign sequences.  It then
projects the phase gate into explicit U^b coarse sectors and compares small
full-orbit sequential_path marginals with sector regroupings.

The N=7,a=3,r=6,b=3,bit-1 sequence is
    z=[+,-,-,-,+,+],  z(j+3)=-z(j).
At theta=pi, Rz(pi)=-i Z therefore routes the two U^3 sectors bijectively;
initial coarse-sector dephasing is valid again even though generic angles
have both within-sector and cross-sector blocks.  For N=13,a=2,r=12,bit 1
has minimal period 6, so the natural regrouping is b'=6, not b=3.

P1: exact integer orbit, explicit pre-derived sign-period/complement relations,
    and projected coarse blocks agree for N=7 bit 0/1/2 and N=13 bits 0/1.
P2: for N=7, theta=pi has a pure sector permutation, while generic theta has
    nonzero within- and cross-sector blocks and differs from initial coarse
    dephasing on the fixed b=3 background.
P3: for N=13 bit 1, the full r=12 sequential marginal agrees with the
    regrouped b'=6 PeriodicOrbitCircuit using the same b=3 background blocks.
P4: N=13 bit-1 differs as an operator from a forced b=3 phase, while the
    forced b=3 and true operators are input-specifically invisible at s=2.

C1: forcing the old b=3 periodic phase on the selected N=13 bit-0 control
    must fail (bit 1 agrees on the first four reached labels at s=2).
C2: deleting cross-sector coherence by initial coarse dephasing must fail for
    the selected generic N=7 angle.
C3: treating the generic rotation as deterministic sector routing must fail.

All exact orbit checks precede dense matrices.  This is an indexed, known-order
diagnostic with r<=24 and t=5; numerical DFTs are diagnostics only, not a
production promise or hardness claim.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
    'numpy<2.5' python -m experiments.experiment_phase_structure
"""
from __future__ import annotations

import json
import math
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import sequential_path


MAX_DENSE_BYTES = 16 * 1024 * 1024
MAX_R = 24
WIDTH = 5
BACKGROUND_THETA = np.pi / 2
N7, A7 = 7, 3
N13, A13 = 13, 2


def guard(shape, dtype, label: str) -> int:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")
    return payload


def report_path(prefix: str = "phase_structure") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def orbit(N: int, base: int) -> np.ndarray:
    if N < 2 or math.gcd(base, N) != 1:
        raise ValueError("orbit input must have gcd(base,N)=1")
    values = [1]
    value = base % N
    while value != 1:
        if len(values) >= MAX_R:
            raise ValueError("orbit exceeds experiment cap")
        values.append(value)
        value = value * base % N
    return np.asarray(values, dtype=np.int64)


def sign_sequence(values: np.ndarray, bit: int) -> np.ndarray:
    return np.where((values >> bit) & 1, -1, 1).astype(np.int8)


def cyclic_min_period(seq: np.ndarray) -> int:
    r = len(seq)
    for d in range(1, r + 1):
        if r % d == 0 and np.array_equal(seq, np.roll(seq, -d)):
            return d
    raise AssertionError("period search failed")


def translation(period: int, step: int) -> np.ndarray:
    guard((period, period), np.complex128, "translation")
    result = np.zeros((period, period), dtype=complex)
    for j in range(period):
        result[(j + step) % period, j] = 1.
    return result


def rx3(kind: str, theta: float = BACKGROUND_THETA) -> np.ndarray:
    if kind not in ("W01", "W12"):
        raise ValueError("unknown 3-block")
    i, j = (0, 1) if kind == "W01" else (1, 2)
    out = np.eye(3, dtype=complex)
    c, s = np.cos(theta / 2), -1j * np.sin(theta / 2)
    out[i, i] = out[j, j] = c
    out[i, j] = out[j, i] = s
    return out


def repeated_block(period: int, block: np.ndarray) -> np.ndarray:
    b = block.shape[0]
    if period % b:
        raise ValueError("period must be a block multiple")
    guard((period, period), np.complex128, "repeated block")
    out = np.zeros((period, period), dtype=complex)
    for start in range(0, period, b):
        out[start:start+b, start:start+b] = block
    return out


def phase_gate(seq: np.ndarray, theta: float) -> np.ndarray:
    seq = np.asarray(seq)
    if seq.ndim != 1 or len(seq) < 1 or len(seq) > MAX_R or not np.isfinite(theta):
        raise ValueError("phase gate exceeds finite orbit cap")
    guard((len(seq), len(seq)), np.complex128, "phase gate")
    return np.diag(np.exp(-0.5j * theta * seq.astype(float))).astype(complex)


def sector_basis(period: int, block: int, alpha: int) -> np.ndarray:
    if period % block or not 0 <= alpha < period // block:
        raise ValueError("invalid coarse sector")
    M = period // block
    guard((period, block), np.complex128, "sector basis")
    q = np.zeros((period, block), dtype=complex)
    for p in range(block):
        for m in range(M):
            q[block*m+p, p] = np.exp(-2j*np.pi*alpha*m/M) / np.sqrt(M)
    return q


def projected_blocks(operator: np.ndarray, period: int, block: int) -> list[list[np.ndarray]]:
    M = period // block
    bases = [sector_basis(period, block, alpha) for alpha in range(M)]
    return [[bases[beta].conj().T @ operator @ bases[alpha]
             for alpha in range(M)] for beta in range(M)]


def branch_pairs(period: int, block: int, phase: np.ndarray | None,
                 *, phase_angle: float = 0.) -> list[tuple[np.ndarray, np.ndarray]]:
    """Original ascending order: W after controls 1,3,4; phase after 2."""
    if period < 1 or period > MAX_R or block != 3 or not np.isfinite(phase_angle):
        raise ValueError("branch dimensions or block promise exceed audit scope")
    if phase is not None and (np.asarray(phase).shape != (period,)
                              or not np.all(np.isfinite(phase))):
        raise ValueError("phase sequence shape/finite check failed")
    guard((period, period), np.complex128, "branch identity")
    identity = np.eye(period, dtype=complex)
    blocks = {1: repeated_block(period, rx3("W01")),
              3: repeated_block(period, rx3("W12")),
              4: repeated_block(period, rx3("W01"))}
    if phase is not None:
        blocks[2] = phase_gate(phase, phase_angle)
    pairs = []
    for i in range(WIDTH):
        post = blocks.get(i + 1, identity)
        shift = translation(period, 1 << i)
        pairs.append((post, post @ shift))
    return pairs


def sequential_marginal(pairs: list[tuple[np.ndarray, np.ndarray]], initial: np.ndarray) -> np.ndarray:
    result = np.zeros(1 << WIDTH, dtype=float)
    for y in range(1 << WIDTH):
        result[y] = sequential_path(pairs, initial, output=y)["conditional_path_probability"]
    return result


def coarse_dephased_marginal(pairs: list[tuple[np.ndarray, np.ndarray]],
                             period: int, block: int) -> np.ndarray:
    M = period // block
    result = np.zeros(1 << WIDTH, dtype=float)
    for alpha in range(M):
        initial = sector_basis(period, block, alpha)[:, 0]
        result += sequential_marginal(pairs, initial) / M
    return result


def tv(left: np.ndarray, right: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(left - right)))


def complex_payload(matrix: np.ndarray) -> dict:
    """JSON-safe complex diagnostic payload without discarding phases."""
    matrix = np.asarray(matrix)
    return dict(real=matrix.real.tolist(), imag=matrix.imag.tolist())


def main() -> None:
    exp = Experiment("phase_structure", doc=__doc__)
    exp.predict("P1", "integer orbit/sign/period and projected sector identities hold")
    exp.predict("P2", "N=7 pi routing is exact while generic angle needs cross-sector coherence")
    exp.predict("P3", "N=13 period-6 phase regrouping matches full r=12 output")
    exp.predict("P4", "N=13 bit-1 b=3 mismatch is operator-level but input-specific invisible at s=2")
    exp.must_fail("C1", "forcing old b=3 periodicity on selected N=13 bit-0 phase")
    exp.must_fail("C2", "initial coarse dephasing at generic N=7 angle")
    exp.must_fail("C3", "deterministic routing approximation for a generic rotation")

    # Exact integer orbit and sign audit, before any Fourier/matrix work.
    rows = []
    integer_rows = []
    expected = {
        (N7, 0): (6, True), (N7, 1): (6, True), (N7, 2): (6, True),
        (N13, 0): (12, True), (N13, 1): (6, False),
    }
    orbit_enumeration_steps = 0
    values_by_N = {}
    for N, base, bits in ((N7, A7, (0, 1, 2)), (N13, A13, (0, 1, 2))):
        values = orbit(N, base)
        values_by_N[N] = values
        r = len(values)
        orbit_enumeration_steps += r
        for bit in bits:
            seq = sign_sequence(values, bit)
            d = cyclic_min_period(seq)
            complement = None
            if r % 2 == 0:
                complement = bool(np.array_equal(np.roll(seq, -(r // 2)), -seq))
            integer_rows.append(dict(N=N, base=base, r=r, bit=bit,
                                     orbit=values.tolist(), signs=seq.tolist(),
                                     minimal_period=d, half_complement=complement))
            if (N, bit) in expected:
                expected_period, expected_half = expected[(N, bit)]
                exp.check("P1", d == expected_period and complement is expected_half,
                          f"N={N}, bit={bit}: r={r}, minimal period={d} "
                          f"(expected {expected_period}), half-complement={complement} "
                          f"(expected {expected_half})")
            else:
                # N=13 bit 2 is retained as a diagnostic control only; no
                # universal period/complement claim is made for it.
                exp.log(f"diagnostic N={N}, bit={bit}: period={d}, half-complement={complement}")
        if N == N7:
            exp.check("P1", np.array_equal(values[(np.arange(r)+3) % r],
                                             (-values) % N),
                      "N=7: a^(j+3) == -a^j mod N")

    n7b1 = next(row for row in integer_rows if row["N"] == N7 and row["bit"] == 1)
    n13b1 = next(row for row in integer_rows if row["N"] == N13 and row["bit"] == 1)
    exp.check("P1", n7b1["signs"] == [1, -1, -1, -1, 1, 1]
              and n7b1["minimal_period"] == 6
              and n13b1["minimal_period"] == 6,
              f"selected signs: N7={n7b1['signs']}, N13 period={n13b1['minimal_period']}")

    # Explicit coarse blocks for N=7 bit 1.  At pi, only off-diagonal sector
    # blocks survive; at a generic angle both diagonal and off-diagonal parts
    # are present.
    values7 = values_by_N[N7]
    period7 = len(values7)
    z7 = sign_sequence(values7, 1)
    phase_pi = phase_gate(z7, np.pi)
    phase_generic = phase_gate(z7, np.pi / 4)
    blocks_pi = projected_blocks(phase_pi, period7, 3)
    blocks_generic = projected_blocks(phase_generic, period7, 3)
    diag_pi = max(float(np.max(np.abs(blocks_pi[a][a]))) for a in range(2))
    cross_pi = max(float(np.max(np.abs(blocks_pi[1-a][a]))) for a in range(2))
    diag_generic = max(float(np.max(np.abs(blocks_generic[a][a]))) for a in range(2))
    cross_generic = max(float(np.max(np.abs(blocks_generic[1-a][a]))) for a in range(2))
    exp.check("P2", diag_pi < 3e-12 and abs(cross_pi - 1.) < 3e-12
              and diag_generic > 1e-4 and cross_generic > 1e-4,
              f"N7 bit1 blocks: pi diag={diag_pi:.2e}, cross={cross_pi:.6g}; "
              f"generic diag={diag_generic:.6g}, cross={cross_generic:.6g}")

    # Full N=7 sequential output with the fixed periodic background.
    initial7 = np.eye(period7, dtype=complex)[:, 0]
    pairs_generic7 = branch_pairs(period7, 3, z7, phase_angle=np.pi / 4)
    pairs_pi7 = branch_pairs(period7, 3, z7, phase_angle=np.pi)
    full_generic7 = sequential_marginal(pairs_generic7, initial7)
    full_pi7 = sequential_marginal(pairs_pi7, initial7)
    dephased_generic7 = coarse_dephased_marginal(pairs_generic7, period7, 3)
    dephased_pi7 = coarse_dephased_marginal(pairs_pi7, period7, 3)
    generic_dephase_tv = tv(full_generic7, dephased_generic7)
    pi_dephase_tv = tv(full_pi7, dephased_pi7)
    exp.check("P2", generic_dephase_tv > 1e-5 and pi_dephase_tv < 3e-10,
              f"N7 output: generic dephase TV={generic_dephase_tv:.6g}, "
              f"pi dephase TV={pi_dephase_tv:.2e}")
    rows.append(dict(series="N7_bit1", orbit=values7.tolist(), signs=z7.tolist(),
                     generic_output=full_generic7.tolist(), generic_dephased=dephased_generic7.tolist(),
                     pi_output=full_pi7.tolist(), pi_dephased=dephased_pi7.tolist(),
                     generic_dephase_tv=generic_dephase_tv, pi_dephase_tv=pi_dephase_tv,
                     projected_pi=[[complex_payload(x) for x in row] for row in blocks_pi],
                     projected_generic=[[complex_payload(x) for x in row]
                                        for row in blocks_generic]))

    # N=13 bit1 has period six.  Group two b=3 blocks into b'=6 and compare
    # the full r=12 instrument with the existing PeriodicOrbitCircuit.
    values13 = values_by_N[N13]
    z13 = sign_sequence(values13, 1)
    phase13 = phase_gate(z13, np.pi / 4)
    full_pairs13 = branch_pairs(len(values13), 3, z13, phase_angle=np.pi / 4)
    full13 = sequential_marginal(full_pairs13, np.eye(len(values13), dtype=complex)[:, 0])
    w01_6 = np.kron(np.eye(2, dtype=complex), rx3("W01"))
    w12_6 = np.kron(np.eye(2, dtype=complex), rx3("W12"))
    # The phase sequence has period six, so this same 6x6 diagonal block is
    # repeated on both b'=6 orbit blocks.
    k6 = np.diag(np.exp(-0.5j * (z13[:6] * np.pi / 4)))
    grouped = PeriodicOrbitCircuit(12, 6, WIDTH,
                                   {1: w01_6, 2: k6, 3: w12_6, 4: w01_6})
    grouped13 = np.array([sum(grouped.forced_joint(alpha, y)[
        "joint_latent_output_probability"] for alpha in range(grouped.sectors))
                          for y in range(1 << WIDTH)])
    grouped_error = float(np.max(np.abs(full13 - grouped13)))
    grouped_tv = tv(full13, grouped13)
    exp.check("P3", grouped_error < 3e-10 and grouped_tv < 3e-10,
              f"N13 b'=6 regrouping: max={grouped_error:.2e}, TV={grouped_tv:.2e}")
    rows.append(dict(series="N13_bit1_grouped_b6", orbit=values13.tolist(), signs=z13.tolist(),
                     minimal_period=cyclic_min_period(z13), full=full13.tolist(),
                     grouped=grouped13.tolist(), max_error=grouped_error, tv=grouped_tv))

    # The old b=3 phase is a different operator for N=13 bit 1, but at s=2
    # only labels 0..3 are reached and those signs happen to agree.  Record
    # this as an input-specific positive result rather than treating it as an
    # operator identity.
    wrong_b3_bit1 = np.tile(z13[:3], 4)
    phase_operator_error = float(np.max(np.abs(
        phase_gate(z13, np.pi / 4) - phase_gate(wrong_b3_bit1, np.pi / 4))))
    bit1_true_s2 = sequential_marginal(
        branch_pairs(12, 3, z13, phase_angle=np.pi / 4), np.eye(12, dtype=complex)[:, 0])
    bit1_wrong_s2 = sequential_marginal(
        branch_pairs(12, 3, wrong_b3_bit1, phase_angle=np.pi / 4), np.eye(12, dtype=complex)[:, 0])
    phase_input_vacuity_tv = tv(bit1_true_s2, bit1_wrong_s2)
    exp.check("P4", phase_operator_error > 1e-4 and phase_input_vacuity_tv < 3e-10,
              f"N13 bit1 operator error={phase_operator_error:.6g}, "
              f"s=2 input TV={phase_input_vacuity_tv:.2e}")

    # C1: use a deliberately wrong b=3 repeated phase on the N=13 bit-0
    # control.  Bit 1 is intentionally not used here: at s=2 its first four
    # reached signs happen to agree with the tiled b=3 sequence, making that
    # particular old-period control vacuous.
    z13_bit0 = sign_sequence(values13, 0)
    wrong_z13 = np.tile(z13_bit0[:3], 4)
    bit0_pairs13 = branch_pairs(12, 3, z13_bit0, phase_angle=np.pi / 4)
    wrong_pairs13 = branch_pairs(12, 3, wrong_z13, phase_angle=np.pi / 4)
    bit0_full13 = sequential_marginal(bit0_pairs13, np.eye(12, dtype=complex)[:, 0])
    wrong13 = sequential_marginal(wrong_pairs13, np.eye(12, dtype=complex)[:, 0])
    wrong_period_error = tv(bit0_full13, wrong13)
    exp.fail_check("C1", wrong_period_error > 1e-5,
                   f"N13 forced b=3 phase TV={wrong_period_error:.9g}")

    # C2 is the generic-angle initial coarse-sector dephasing already
    # measured above; resolve the must-fail control explicitly.
    exp.fail_check("C2", generic_dephase_tv > 1e-5,
                   f"N7 generic initial coarse dephasing TV={generic_dephase_tv:.9g}")

    # C3: replace generic theta=pi/4 by deterministic pi-style sector routing
    # (same phase signs, but no within-sector component).
    deterministic_pairs7 = branch_pairs(period7, 3, z7, phase_angle=np.pi)
    deterministic7 = sequential_marginal(deterministic_pairs7, initial7)
    deterministic_error = tv(full_generic7, deterministic7)
    exp.fail_check("C3", deterministic_error > 1e-5,
                   f"generic vs deterministic-routing TV={deterministic_error:.9g}")

    path = report_path()
    exp.finish(report_path=path, rows=rows, metadata=dict(
        width=WIDTH, max_r=MAX_R, dense_budget_bytes=MAX_DENSE_BYTES,
        integer_rows=integer_rows, orbit_enumeration_steps=orbit_enumeration_steps,
        selected_N7_bit1=n7b1,
        selected_N13_bit1=n13b1, generic_angle=float(np.pi / 4), pi_angle=float(np.pi),
        background="W01@s1,W12@s3,W01@s4, theta=pi/2",
        N7_generic_dephase_tv=generic_dephase_tv, N7_pi_dephase_tv=pi_dephase_tv,
        N13_grouped_b6_error=grouped_error, N13_wrong_b3_tv=wrong_period_error,
        N13_bit1_b3_operator_error=phase_operator_error,
        N13_bit1_s2_input_vacuity_tv=phase_input_vacuity_tv,
        deterministic_routing_tv=deterministic_error,
        reference="lab.sequential_path and PeriodicOrbitCircuit; exact integer orbit setup charged",
        numerical_DFT_role="diagnostic projected blocks only; no tolerance-selected production promise"))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("phase_structure_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
