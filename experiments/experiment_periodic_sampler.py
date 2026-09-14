"""Production and physical validation for the periodic-orbit sampler.

This is the first bounded test of ``lab.periodic.PeriodicOrbitCircuit`` against
an actual original-time-order modular-exponentiation circuit.  The physical
instance is N=7, a=3, with clean orbit [1,3,2,6,4,5], block size b=3, and two
noncommuting repeated blocks: W01 mixes orbit positions (0,1), while W12
mixes positions (1,2).  Both are implemented in the existing Circuit engine
as work-register Pauli rotations; no circuit propagator is added here.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  For fixed t=5 and insertions s=1 (W01), s=3 (W12), optionally s=4
      (W01), the sum of all coarse-sector forced joints equals the physical
      exponent marginal.  The one-defect row also agrees with lab.spectral.
  P2  Identity W recovers the known ideal order-finding distribution; a W at
      s=t is output-invisible because it is a traced work-only operation after
      all arithmetic powers.
  P3  The production sampler remains bounded at t=63 for r=6 and for the
      indexed abstract period r=3*1_000_000_007, with no orbit/prefix/output
      arrays and only b-dimensional sampled work factors.
  P4  The two physical W constructions are unitary, match their intended
      8x8 work matrices, and leave invalid labels 0 and 7 unchanged.

  C1  Omitting inverse-QFT feedback must fail against the physical marginal.
  C2  Replacing coarse-sector dephasing by a single alpha=0 sector must fail.
  C3  Reversing the arithmetic/defect schedule while preserving initialization
      and final inverse QFT must fail; the defect-free reverse remains valid.

The sampler's coarse latent variable is alpha in [0,r/b), not the final fine
orbit eigenphase k.  Reports preserve physical state-vector setup costs,
sampler payloads, allocation guards, and all failures.  This is a validation
of a known-order indexed-orbit construction, not order discovery, factoring,
or a physical locality claim.

Run from research/:
  OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with
  'numpy<2.5' python -m experiments.experiment_periodic_sampler
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
from lab.periodic import PeriodicOrbitCircuit
from lab.semiclassical import order_finding_probability
from lab.spectral import single_defect_effects
from modexp import ModExp
import statevec


MAX_DENSE_BYTES = 16 * 1024 * 1024
N, BASE, WIDTH, PERIOD, BLOCK = 7, 3, 5, 6, 3
ORBIT = np.array([1, 3, 2, 6, 4, 5], dtype=np.int64)
THETA = np.pi / 2


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def report_path(prefix: str = "periodic_sampler") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    root = Path("out")
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"{prefix}_{stamp}.json"
    serial = 0
    while path.exists():
        serial += 1
        path = root / f"{prefix}_{stamp}_{serial}.json"
    return path


def complex_payload(matrix: np.ndarray) -> dict:
    return dict(real=np.asarray(matrix).real.tolist(), imag=np.asarray(matrix).imag.tolist())


def block_rotation(kind: str, theta: float) -> np.ndarray:
    """3x3 orbit-index Rx block, with the untouched local basis explicit."""
    if kind not in ("W01", "W12"):
        raise ValueError("unknown block kind")
    target, other = ((0, 1) if kind == "W01" else (1, 2))
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    w = np.eye(3, dtype=complex)
    w[target, target] = w[other, other] = c
    w[target, other] = w[other, target] = -1j * s
    return w


def work_mixer(kind: str, theta: float, work_qubits: list[int]) -> Circuit:
    """Physical repeated block using the specified work-bit parity projector."""
    if kind == "W01":
        target = work_qubits[1]
        spectators = (1 << work_qubits[0]) | (1 << work_qubits[2])
    elif kind == "W12":
        target = work_qubits[0]
        spectators = (1 << work_qubits[1]) | (1 << work_qubits[2])
    else:
        raise ValueError("unknown block kind")
    qc = Circuit(max(work_qubits) + 1)
    # X_target (I-Z_spectators)/2.  The two Pauli rotations commute.
    qc.rot((1 << target, 0), theta / 2)
    qc.rot((1 << target, spectators), -theta / 2)
    return qc


def expected_work8(kind: str, theta: float) -> np.ndarray:
    """Expected computational-basis 8x8 matrix for the physical block."""
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    result = np.eye(8, dtype=complex)
    pairs = ((1, 3), (6, 4)) if kind == "W01" else ((3, 2), (4, 5))
    for left, right in pairs:
        result[left, left] = result[right, right] = c
        result[left, right] = result[right, left] = -1j * s
    return result


def validate_physical_blocks(exp: Experiment) -> dict:
    rows = []
    for kind in ("W01", "W12"):
        local = work_mixer(kind, THETA, [0, 1, 2])
        guard((8, 8), np.complex128, f"{kind} local work unitary")
        actual = local.to_unitary()
        expected = expected_work8(kind, THETA)
        unitary_error = float(np.max(np.abs(actual.conj().T @ actual - np.eye(8))))
        matrix_error = float(np.max(np.abs(actual - expected)))
        invalid_error = 0.
        statevec_error = 0.
        for label in range(8):
            basis = statevec.basis(3, label)
            evolved = statevec.run(local, basis)
            statevec_error = max(statevec_error,
                                 float(np.max(np.abs(evolved - actual @ basis))))
            if label in (0, 7):
                invalid_error = max(invalid_error,
                                    float(np.max(np.abs(evolved - basis))))
        exp.check("P4", unitary_error < 2e-14 and matrix_error < 2e-14
                  and statevec_error < 2e-14 and invalid_error < 2e-14,
                  f"{kind}: unitary={unitary_error:.2e}, matrix={matrix_error:.2e}, "
                  f"statevec={statevec_error:.2e}, invalid={invalid_error:.2e}")
        rows.append(dict(kind=kind, unitary_error=unitary_error,
                         matrix_error=matrix_error, statevec_error=statevec_error,
                         invalid_label_error=invalid_error,
                         actual=complex_payload(actual), expected=complex_payload(expected)))
    return dict(rows=rows)


def ideal_output(period: int, width: int) -> np.ndarray:
    if not 0 <= width <= 8 or not 1 <= period <= 18:
        raise ValueError("ideal enumeration cap")
    return np.array([order_finding_probability(y, width, period)
                     for y in range(1 << width)], dtype=float)


def sampler_marginal(circuit: PeriodicOrbitCircuit) -> tuple[np.ndarray, float]:
    if circuit.width > 8 or circuit.period > 18:
        raise ValueError("marginal enumeration cap (not a large-period sampler)")
    result = np.zeros(1 << circuit.width, dtype=float)
    max_joint_mass_error = 0.
    for alpha in range(circuit.sectors):
        mass = 0.
        for y in range(1 << circuit.width):
            row = circuit.forced_joint(alpha, y)
            mass += row["joint_latent_output_probability"]
            result[y] += row["joint_latent_output_probability"]
        max_joint_mass_error = max(max_joint_mass_error, abs(mass - 1/circuit.sectors))
    return result, max_joint_mass_error


def repeated_orbit_matrix(block: np.ndarray, period: int = PERIOD) -> np.ndarray:
    if not 1 <= period <= 18 or period % BLOCK:
        raise ValueError("period must be divisible by block size")
    guard((period, period), np.complex128, "repeated orbit matrix")
    return np.kron(np.eye(period // BLOCK, dtype=complex), block)


def orbit_fourier(period: int) -> np.ndarray:
    guard((period, period), np.complex128, "orbit Fourier")
    j = np.arange(period, dtype=float)
    return np.exp(-2j * np.pi * j[:, None] * j[None, :] / period) / np.sqrt(period)


def spectral_marginal(v_orbit: np.ndarray, period: int, width: int, split: int) -> np.ndarray:
    f = orbit_fourier(period)
    eigen_defect = f.conj().T @ v_orbit @ f
    effects = single_defect_effects(period, width, split, eigen_defect,
                                    max_entries=1_000_000)
    return effects.sum(axis=(1, 2)).real / period


def physical_circuit(defects: dict[int, str], width: int = WIDTH,
                     *, reverse_all: bool = False, feedback: bool = True) -> tuple[Circuit, ModExp]:
    if not 1 <= width <= 8:
        raise ValueError("physical reference width cap")
    me = ModExp(N, BASE, n_exp=width)
    qc = Circuit(me.n_qubits).x(me.x[0])
    for q in me.exp:
        qc.h(q)
    stages = []
    for i, q in enumerate(me.exp):
        stages.append(me.u_a(q, pow(BASE, 1 << i, N)))
        kind = defects.get(i + 1)
        if kind is not None:
            stages.append(work_mixer(kind, THETA, me.x))
    for stage in reversed(stages) if reverse_all else stages:
        qc.extend(stage)
    if feedback:
        qc.qft(me.exp, inverse=True)
    else:
        # Deleting all adaptive/controlled Fourier phases leaves Hadamards
        # plus the SAME bit reversal; arithmetic and initialization unchanged.
        for q in me.exp:
            qc.h(q)
        for i in range(width//2):
            qc.swap(me.exp[i], me.exp[-1-i])
    return qc, me


def physical_output(defects: dict[int, str], width: int = WIDTH,
                    *, reverse_all: bool = False, feedback: bool = True) -> tuple[np.ndarray, dict]:
    qc, me = physical_circuit(defects, width, reverse_all=reverse_all, feedback=feedback)
    payload = (1 << me.n_qubits) * np.dtype(np.complex128).itemsize
    if payload > MAX_DENSE_BYTES:
        raise MemoryError(f"physical state allocation {payload} bytes exceeds 16 MiB")
    psi = np.zeros(1 << me.n_qubits, dtype=complex)
    psi[0] = 1.
    psi = statevec.run(qc, psi)
    probabilities = np.sum(np.abs(psi.reshape(1 << width, -1)) ** 2, axis=1)
    return probabilities, dict(qubits=me.n_qubits, vector_payload_bytes=psi.nbytes,
                               norm_error=abs(float(np.vdot(psi, psi).real) - 1),
                               gates=len(qc.gates), rotations=qc.n_nonclifford())


def main() -> None:
    exp = Experiment("periodic_sampler", doc=__doc__)
    exp.predict("P1", "multi-defect coarse forced joints match original-order physical output")
    exp.predict("P2", "identity and s=t controls recover ideal output")
    exp.predict("P3", "t=63 physical-dimension sampler remains bounded for small and huge indexed r")
    exp.predict("P4", "physical repeated blocks match unitary/statevec constructions")
    exp.must_fail("C1", "removing inverse-QFT feedback")
    exp.must_fail("C2", "wrong single-alpha coarse dephasing")
    exp.must_fail("C3", "reversing arithmetic/defect stages with initialization and QFT unchanged")
    start = time.perf_counter()
    rows = []
    setup = validate_physical_blocks(exp)

    # Main physical multi-defect row: W01 after s=1, W12 after s=3, W01 after s=4.
    defects = {1: "W01", 3: "W12", 4: "W01"}
    blocks = {1: block_rotation("W01", THETA), 3: block_rotation("W12", THETA),
              4: block_rotation("W01", THETA)}
    sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, WIDTH,
                                   {s: w for s, w in blocks.items()})
    physical, physical_info = physical_output(defects)
    sampled, joint_mass_error = sampler_marginal(sampler)
    multi_error = float(np.max(np.abs(sampled - physical)))
    exp.check("P1", multi_error < 2e-9 and joint_mass_error < 2e-11
              and physical_info["norm_error"] < 2e-12,
              f"multi defects: sampler/physical={multi_error:.2e}, "
              f"joint mass={joint_mass_error:.2e}, norm={physical_info['norm_error']:.2e}")
    rows.append(dict(name="physical_multi_defect", defects=defects,
                     physical=physical.tolist(), sampler=sampled.tolist(),
                     max_error=multi_error, joint_mass_error=joint_mass_error,
                     physical_info=physical_info, stats=sampler.stats()))

    # One-defect independent spectral route at width 8, as a separate reference.
    one_width, one_split = 8, 1
    one_w = block_rotation("W01", THETA)
    one_sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, one_width, {1: one_w})
    one_sampled, one_mass_error = sampler_marginal(one_sampler)
    one_spectral = spectral_marginal(repeated_orbit_matrix(one_w), PERIOD, one_width, one_split)
    one_spectral_error = float(np.max(np.abs(one_sampled - one_spectral)))
    exp.check("P1", one_spectral_error < 3e-9 and one_mass_error < 2e-11,
              f"one W01 spectral: sampler/spectral={one_spectral_error:.2e}, "
              f"joint mass={one_mass_error:.2e}")
    rows.append(dict(name="one_defect_spectral", width=one_width, split=one_split,
                     sampler=one_sampled.tolist(), spectral=one_spectral.tolist(),
                     max_error=one_spectral_error, joint_mass_error=one_mass_error))

    ideal = ideal_output(PERIOD, WIDTH)
    identity_sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, WIDTH, {})
    identity_output_values, identity_mass_error = sampler_marginal(identity_sampler)
    identity_error = float(np.max(np.abs(identity_output_values - ideal)))
    exp.check("P2", identity_error < 3e-9 and identity_mass_error < 2e-11,
              f"identity: sampler/ideal={identity_error:.2e}, mass={identity_mass_error:.2e}")

    end_sampler = PeriodicOrbitCircuit(PERIOD, BLOCK, WIDTH, {WIDTH: blocks[1]})
    end_output, end_mass_error = sampler_marginal(end_sampler)
    end_error = float(np.max(np.abs(end_output - ideal)))
    exp.check("P2", end_error < 3e-9 and end_mass_error < 2e-11,
              f"s=t: sampler/ideal={end_error:.2e}, mass={end_mass_error:.2e}")
    rows.append(dict(name="identity_and_end", identity_error=identity_error,
                     end_error=end_error, identity_mass_error=identity_mass_error,
                     end_mass_error=end_mass_error))

    # Intentionally wrong controls against the genuine physical marginal.
    wrong_feedback, _ = physical_output(defects, feedback=False)
    no_feedback_tv = float(np.sum(np.abs(wrong_feedback - physical)) / 2)
    exp.fail_check("C1", no_feedback_tv > 1e-5,
                   f"no-feedback TV from physical={no_feedback_tv:.6g}")

    alpha_zero = np.array([sampler.forced_joint(0, y)["conditional_path_probability"]
                           for y in range(1 << WIDTH)])
    wrong_coarse_tv = float(np.sum(np.abs(alpha_zero - physical)) / 2)
    exp.fail_check("C2", wrong_coarse_tv > 1e-5,
                   f"single-alpha=0 TV from physical={wrong_coarse_tv:.6g}")

    reversed_physical, reversed_info = physical_output(defects, reverse_all=True)
    reversed_tv = float(np.sum(np.abs(reversed_physical - physical)) / 2)
    exp.fail_check("C3", reversed_tv > 1e-5,
                   f"reversed-arithmetic/defect-schedule TV={reversed_tv:.6g}")
    clean_reversed, _ = physical_output({}, reverse_all=True)
    clean_reverse_error = float(np.max(np.abs(clean_reversed-ideal)))
    exp.check("P2", clean_reverse_error < 2e-11,
              f"defect-free reversal is valid: max error={clean_reverse_error:.2e}")
    rows.append(dict(name="controls", no_feedback_tv=no_feedback_tv,
                     wrong_coarse_tv=wrong_coarse_tv, reversed_tv=reversed_tv,
                     reversed_info=reversed_info, clean_reverse_error=clean_reverse_error))

    # Low-cost production sampler rows: no orbit, prefix, output or dense-work arrays.
    large_rows = []
    for period, label in ((6, "r6"), (3 * 1_000_000_007, "r3000000021")):
        large = PeriodicOrbitCircuit(period, BLOCK, 63,
                                     {16: blocks[1], 32: blocks[3], 48: blocks[1]})
        rng = np.random.default_rng(20260910 + period % 1000)
        draws = []
        for _ in range(8):
            draw = large.sample(rng)
            draws.append(dict(output=int(draw["output"]), sector=int(draw["coarse_eigenphase"]),
                              probability=float(draw["conditional_path_probability"]),
                              joint_probability=float(draw["joint_latent_output_probability"]),
                              forward_factor_payload_bytes=draw["forward_factor_payload_bytes"],
                              branch_matrix_payload_bytes=draw["branch_matrix_payload_bytes"],
                              conservative_payload_estimate_bytes=draw["conservative_payload_estimate_bytes"],
                              forward_normalization_error=draw["forward_normalization_error"]))
        stats = large.stats()
        bounded = (stats["sampled_work_dimension"] == BLOCK
                   and stats["orbit_table_entries"] == 0
                   and stats["early_vector_entries"] == 0
                   and stats["output_table_entries"] == 0
                   and stats["owned_block_payload_bytes"] <= MAX_DENSE_BYTES
                   and all(0 <= d["output"] < 1 << 63 and 0 < d["probability"] <= 1
                           and 0 < d["joint_probability"] <= 1
                           and d["forward_factor_payload_bytes"] == 64*9*16
                           for d in draws))
        exp.check("P3", bounded,
                  f"{label}: stats={stats}, draws={len(draws)}")
        large_rows.append(dict(name=label, draws=draws, stats=stats, bounded=bounded))
    rows.extend(large_rows)

    elapsed = time.perf_counter() - start
    path = report_path()
    exp.finish(report_path=path, rows=rows,
               metadata=dict(N=N, base=BASE, orbit=ORBIT.tolist(), period=PERIOD,
                             block_size=BLOCK, width=WIDTH, theta=float(THETA),
                             defects=defects, max_dense_bytes=MAX_DENSE_BYTES,
                             physical_setup=setup,
                             production_contract="coarse alpha sectors; known order/orbit indices",
                             no_amplitude_cutoff=True,
                             elapsed_seconds=elapsed, numpy=np.__version__))
    print(f"report: {path}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("periodic_sampler_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc(),
                                           python=__import__("sys").version), indent=2) + "\n")
        raise
