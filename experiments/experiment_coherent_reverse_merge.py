"""Bounded law audit of the coherent sparse reverse merge.

The candidate reverse sampler keeps complex fine-work vectors for the reached
coarse sectors and merges histories by coherent amplitude addition.  This
experiment checks the complete accepted (coarse-sector, QFT-output) law, not
only a proposal marginal.  The primary fixture is r=9, b=3, t=4 with an
insertion-0 mixer, noncommuting W1/W3 mixers, and four q=0/1 reflections.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  accepted_joint_submass sums to 1/(b*S), where S is the reported global
      peak-sector support bound, and its globally normalized law agrees with
      both the existing coherent joint_probability and independent direct_joint.
  P2  every forced proposal law normalizes, the coherent acceptance equals the
      supplied terminal coherent numerator, and no history table/enumeration is
      used; sampled counters charge all proposals.
  P3  endpoint, fixed-point, modular-wrap, collision, t=0/1, and b=2 cases
      obey the same accepted-law identity.

  C1  omitting terminal acceptance changes the target law.
  C2  freezing gamma through rejection and then mixing normalized gamma laws
      changes the target law.
  C3  replacing the coherent terminal numerator by the incoherent sum of
      sector norms changes the target law.
  C4  merging colliding sector vectors coherently differs from merging their
      probabilities, including cancellation at a two-cycle fixture.

This is a bounded floating-point implementation/reference diagnostic, not a
precision certificate.  It does not introduce a generic propagator: the
independent reference is the existing tiny full-r direct_joint helper.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit


MAX_DENSE_BYTES = 16 << 20
TOL = 2e-10

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "coherent accepted submass and normalized joint law agree with both references")
exp.predict("P2", "proposal/terminal identities and table-free charged counters hold")
exp.predict("P3", "small endpoint, fixed-point, wrap, collision, and b=2 cases pass")
exp.must_fail("C1", "omitting terminal acceptance changes the target law")
exp.must_fail("C2", "freezing gamma through rejection changes the target law")
exp.must_fail("C3", "incoherent terminal sector norms change the target law")
exp.must_fail("C4", "probability-only merging misses coherent collision cancellation")


def guard(shape, dtype=np.float64, label="array"):
    entries = math.prod(int(x) for x in shape)
    payload = entries * np.dtype(dtype).itemsize
    if entries < 0 or payload > MAX_DENSE_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def stamp(prefix="coherent_reverse_merge"):
    now = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{now}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    serial = 0
    while path.exists():
        serial += 1
        path = Path("out") / f"{prefix}_{now}_{serial}.json"
    return path


def rx(theta, block):
    out = np.eye(block, dtype=complex)
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    out[:2, :2] = ((c, -1j * s), (-1j * s, c))
    return out


def rz(theta, block):
    out = np.eye(block, dtype=complex)
    out[0, 0] = np.exp(-1j * theta / 2)
    out[1, 1] = np.exp(1j * theta / 2)
    return out


def fixture_primary():
    """M=3: fixed points and repeated q=0/1 words collide by construction."""
    return CoherentReflectionCircuit(
        9, 3, 4,
        {0: rx(np.pi / 4, 3),
         1: rx(np.pi / 7, 3),
         3: rz(np.pi / 5, 3)},
        {0: (0, np.pi / 5),
         2: (1, np.pi / 3),
         3: (0, np.pi / 7),
         4: (1, np.pi / 4)})


def fixture_t0():
    return CoherentReflectionCircuit(
        9, 3, 0, {0: rx(np.pi / 6, 3)}, {0: (0, np.pi / 3)})


def fixture_t1():
    return CoherentReflectionCircuit(
        9, 3, 1, {0: rx(np.pi / 6, 3), 1: rz(np.pi / 5, 3)},
        {0: (1, np.pi / 3), 1: (0, np.pi / 4)})


def fixture_b2():
    """M=5: wrap and fixed-point routing with a noncommuting W schedule."""
    return CoherentReflectionCircuit(
        10, 2, 4,
        {0: rx(np.pi / 5, 2), 1: rz(np.pi / 7, 2), 3: rx(np.pi / 6, 2)},
        {0: (4, np.pi / 5), 2: (0, np.pi / 3),
         3: (1, np.pi / 7), 4: (4, np.pi / 4)})


def law_array(circuit):
    shape = (circuit.sectors, 1 << circuit.width)
    guard(shape, np.float64, "joint law")
    return np.zeros(shape, dtype=float)


def direct_reference(circuit):
    # Existing helper owns the independent full-r construction and its cap.
    from experiments.experiment_coherent_route_sampling import direct_joint
    return np.asarray(direct_joint(circuit), dtype=float)


def target_reference(circuit):
    target = law_array(circuit)
    for gamma in range(circuit.sectors):
        for output in range(1 << circuit.width):
            target[gamma, output] = float(circuit.joint_probability(gamma, output))
    return target


def total_variation(left, right):
    return float(np.sum(np.abs(np.asarray(left) - np.asarray(right))) / 2)


def finite(v):
    return np.isfinite(np.asarray(v, dtype=float)).all()


def expected_reverse_cost(circuit, final_sector):
    """Exact counter recurrence from the reached-sector support sets."""
    current = {int(final_sector)}
    work_matvecs = 0
    qft_pairs = 0
    reflection_contributions = 0
    peak = 1
    for depth in range(circuit.width):
        insertion = circuit.width - depth
        before = len(current)
        if insertion in circuit.reflections:
            q, _ = circuit.reflections[insertion]
            current = current | {(-alpha - q) % circuit.sectors
                                 for alpha in current}
            reflection_contributions += 2 * before
            peak = max(peak, len(current))
        work_matvecs += 3 * len(current)
        qft_pairs += len(current)
    before = len(current)
    if 0 in circuit.reflections:
        q, _ = circuit.reflections[0]
        current = current | {(-alpha - q) % circuit.sectors
                             for alpha in current}
        reflection_contributions += 2 * before
        peak = max(peak, len(current))
    work_matvecs += len(current)
    return dict(work_matvecs=work_matvecs,
                qft_matrix_pair_constructions=qft_pairs,
                reflection_vector_contributions=reflection_contributions,
                peak_sector_count=peak)


def collision_control(SparseCoherentReverse):
    """A two-source collision cancels one label and doubles the other."""
    circuit = CoherentReflectionCircuit(
        10, 2, 0, {}, {0: (4, np.pi / 2)})
    worker = SparseCoherentReverse(circuit)
    e0 = np.array([1.0 + 0j, 0j])
    vector = {0: e0, 1: 1j * e0}
    cost = dict(peak_sector_count=1, reflection_vector_contributions=0)
    merged = worker._reflection_adjoint(vector, 0, cost)
    norms = np.array([float(np.vdot(merged.get(i, np.zeros(2, complex)),
                                  merged.get(i, np.zeros(2, complex))).real)
                      for i in (0, 1)])
    coherent = norms / norms.sum()
    probability_only = np.array([0.5, 0.5])
    return dict(norms=norms.tolist(),
                coherent_vs_probability_tv=total_variation(coherent,
                                                            probability_only),
                cancellation=bool(abs(norms[0]) < TOL
                                  and abs(norms[1] - 2.) < TOL),
                cost=cost)


def run_fixture(name, circuit, SparseCoherentReverse, *, sample_seed):
    # Guard every retained law before any production/reference allocation.
    guard((circuit.sectors, 1 << circuit.width), np.float64,
          f"{name} law")
    worker = SparseCoherentReverse(circuit)
    target = target_reference(circuit)
    direct = direct_reference(circuit)
    if not finite(target) or not finite(direct):
        raise AssertionError(f"{name}: nonfinite reference law")
    target_mass = float(target.sum())
    direct_mass = float(direct.sum())

    proposal = law_array(circuit)
    accepted = law_array(circuit)
    incoherent_accepted = law_array(circuit)
    coherent_identity_error = 0.0
    accepted_identity_error = 0.0
    proposal_mass_error = 0.0
    history_counts = []
    support_bounds = []
    peak_counts = []
    mathematical_means = []
    terminal_amplitude_error = 0.0
    attempt_cost_ok = True
    forced_rows = 0
    for gamma in range(circuit.sectors):
        for boundary in range(circuit.b):
            per_boundary = law_array(circuit)
            for output in range(1 << circuit.width):
                result = worker.attempt(gamma, boundary, output=output)
                required = ("proposal_path_probability", "acceptance_probability",
                            "accepted_joint_submass", "support_bound",
                            "peak_sector_count", "mathematical_mean_attempts",
                            "terminal_incoherent_numerator", "terminal_norm",
                            "terminal_coherent_numerator", "terminal_coherent_amplitude",
                            "history_enumerations", "reverse_steps",
                            "reflection_vector_contributions", "work_matvecs",
                            "qft_matrix_pair_constructions",
                            "max_relative_completeness_error")
                missing = [key for key in required if key not in result]
                if missing:
                    raise AssertionError(f"{name}: API result missing {missing}")
                q = float(result["proposal_path_probability"])
                a = float(result["acceptance_probability"])
                accepted_mass = float(result["accepted_joint_submass"])
                if not all(np.isfinite(x) for x in (q, a, accepted_mass)):
                    raise AssertionError(f"{name}: nonfinite forced result")
                proposal[gamma, output] += q / (circuit.sectors * circuit.b)
                accepted[gamma, output] += accepted_mass
                per_boundary[gamma, output] += q
                expected_accepted = q * a / (circuit.sectors * circuit.b)
                accepted_identity_error = max(
                    accepted_identity_error, abs(accepted_mass - expected_accepted))
                coherent_acceptance = float(result["terminal_coherent_numerator"])
                norm = float(result["terminal_norm"])
                support = float(result["support_bound"])
                if norm <= 0 or support <= 0:
                    raise AssertionError(f"{name}: invalid terminal support/norm")
                coherent_acceptance /= support * norm
                coherent_identity_error = max(
                    coherent_identity_error, abs(a - coherent_acceptance))
                expected_amplitude = complex(circuit.prefix_vector(
                    gamma, 0, circuit.width, measured=circuit.width,
                    output=output)[boundary])
                terminal_amplitude_error = max(
                    terminal_amplitude_error,
                    abs(complex(result["terminal_coherent_amplitude"]).conjugate()
                        - expected_amplitude))
                incoherent_acceptance = (
                    float(result["terminal_incoherent_numerator"])
                    / (support * norm))
                incoherent_accepted[gamma, output] += (
                    q * incoherent_acceptance / (circuit.sectors * circuit.b))
                history_counts.append(int(result["history_enumerations"]))
                support_bounds.append(int(result["support_bound"]))
                peak_counts.append(int(result["peak_sector_count"]))
                mathematical_means.append(float(result["mathematical_mean_attempts"]))
                expected_cost = expected_reverse_cost(circuit, gamma)
                attempt_cost_ok &= (
                    int(result["reverse_steps"]) == circuit.width
                    and int(result["reflection_vector_contributions"])
                        == expected_cost["reflection_vector_contributions"]
                    and int(result["work_matvecs"])
                        == expected_cost["work_matvecs"]
                    and int(result["qft_matrix_pair_constructions"])
                        == expected_cost["qft_matrix_pair_constructions"]
                    and int(result["peak_sector_count"])
                        == expected_cost["peak_sector_count"]
                    and float(result["max_relative_completeness_error"]) < 1e-10)
                forced_rows += 1
            proposal_mass_error = max(
                proposal_mass_error, abs(float(per_boundary[gamma].sum()) - 1.0))

    if target_mass <= 0 or direct_mass <= 0:
        raise AssertionError(f"{name}: reference law has no mass")
    accepted_total = float(accepted.sum())
    support = max(support_bounds)
    peak = max(peak_counts)
    if support < peak or support <= 0:
        raise AssertionError(f"{name}: invalid support bound {support} < peak {peak}")
    expected_total = 1.0 / (circuit.b * support)
    accepted_normalized = accepted / accepted_total
    incoherent_total = float(incoherent_accepted.sum())
    incoherent_normalized = incoherent_accepted / incoherent_total
    proposal_normalized = proposal / float(proposal.sum())
    frozen = np.zeros_like(target)
    for gamma in range(circuit.sectors):
        row_mass = float(accepted[gamma].sum())
        if row_mass <= 0:
            raise AssertionError(f"{name}: zero accepted mass for gamma={gamma}")
        frozen[gamma] = accepted[gamma] / row_mass / circuit.sectors

    # Charge every sampled attempt independently from its starting label;
    # nonnegative aggregate counters alone would let missing work pass.
    sampled_cost = dict(work_matvecs=0, qft_matrix_pair_constructions=0,
                        reflection_vector_contributions=0, peak_sector_count=0)
    sampled_attempts = 0
    original_attempt = worker.attempt
    def counted_attempt(gamma, boundary, **kwargs):
        nonlocal sampled_attempts
        sampled_attempts += 1
        expected = expected_reverse_cost(circuit, gamma)
        for key in sampled_cost:
            if key == "peak_sector_count":
                sampled_cost[key] = max(sampled_cost[key], expected[key])
            else:
                sampled_cost[key] += expected[key]
        return original_attempt(gamma, boundary, **kwargs)
    worker.attempt = counted_attempt
    sample = worker.sample(np.random.default_rng(sample_seed), max_proposals=100_000)
    sample_required = ("output", "final_coarse_sector", "rejection_proposals",
                       "history_enumerations", "mathematical_mean_attempts",
                       "peak_sector_count", "support_bound", "reverse_steps",
                       "reflection_vector_contributions", "work_matvecs",
                       "qft_matrix_pair_constructions",
                       "max_relative_completeness_error")
    missing = [key for key in sample_required if key not in sample]
    if missing:
        raise AssertionError(f"{name}: sample result missing {missing}")
    sample_ok = (
        0 <= int(sample["output"]) < (1 << circuit.width)
        and 0 <= int(sample["final_coarse_sector"]) < circuit.sectors
        and 1 <= int(sample["rejection_proposals"]) <= 100_000
        and int(sample["history_enumerations"]) == 0
        and int(sample["peak_sector_count"]) <= int(sample["support_bound"])
        and int(sample["reverse_steps"]) ==
            int(sample["rejection_proposals"]) * circuit.width
        and int(sample["rejection_proposals"]) == sampled_attempts
        and all(int(sample[key]) == expected
                for key, expected in sampled_cost.items())
        and float(sample["max_relative_completeness_error"]) < 1e-10
        and abs(float(sample["mathematical_mean_attempts"])
                - circuit.b * int(sample["support_bound"])) < TOL)
    return {
        "name": name,
        "period": circuit.period,
        "block_size": circuit.b,
        "width": circuit.width,
        "forced_rows": forced_rows,
        "target_mass": target_mass,
        "direct_mass": direct_mass,
        "target_mass_error": abs(target_mass - 1.0),
        "direct_mass_error": abs(direct_mass - 1.0),
        "accepted_total": accepted_total,
        "expected_accepted_total": expected_total,
        "accepted_mass_error": abs(accepted_total - expected_total),
        "proposal_mass_error": proposal_mass_error,
        "coherent_acceptance_identity_error": coherent_identity_error,
        "terminal_amplitude_error": terminal_amplitude_error,
        "accepted_submass_identity_error": accepted_identity_error,
        "max_history_enumerations": max(history_counts),
        "attempt_cost_ok": attempt_cost_ok,
        "support_bound": support,
        "peak_sector_count": peak,
        "mathematical_mean_attempts": sorted(set(mathematical_means)),
        "target_vs_direct_tv": total_variation(target, direct),
        "accepted_vs_target_tv": total_variation(accepted_normalized, target),
        "accepted_vs_direct_tv": total_variation(accepted_normalized, direct),
        "proposal_vs_target_tv": total_variation(proposal_normalized, target),
        "frozen_gamma_vs_target_tv": total_variation(frozen, target),
        "incoherent_vs_target_tv": total_variation(incoherent_normalized, target),
        "sample_ok": sample_ok,
        "sample_rejection_proposals": int(sample["rejection_proposals"]),
        "sample_history_enumerations": int(sample["history_enumerations"]),
        "target": target.tolist(),
        "direct": direct.tolist(),
        "accepted_normalized": accepted_normalized.tolist(),
        "proposal_normalized": proposal_normalized.tolist(),
        "incoherent_normalized": incoherent_normalized.tolist(),
    }


def main():
    started = time.perf_counter()
    report = {"status": "PASS", "rows": [], "controls": {}}
    try:
        from lab.coherent_reverse import SparseCoherentReverse

        fixtures = (("primary_r9_b3_t4", fixture_primary()),
                    ("endpoint_t0", fixture_t0()),
                    ("endpoint_t1", fixture_t1()),
                    ("wrap_fixed_collision_r10_b2", fixture_b2()))
        for index, (name, circuit) in enumerate(fixtures):
            report["rows"].append(run_fixture(
                name, circuit, SparseCoherentReverse, sample_seed=991 + index))

        primary = report["rows"][0]
        p1 = all(row["accepted_mass_error"] <= TOL
                 and row["target_mass_error"] <= TOL
                 and row["direct_mass_error"] <= TOL
                 and row["accepted_vs_target_tv"] <= TOL
                 and row["accepted_vs_direct_tv"] <= TOL
                 and row["target_vs_direct_tv"] <= TOL
                 for row in report["rows"])
        p2 = all(row["proposal_mass_error"] <= TOL
                 and row["coherent_acceptance_identity_error"] <= TOL
                 and row["terminal_amplitude_error"] <= TOL
                 and row["attempt_cost_ok"]
                 and row["accepted_submass_identity_error"] <= TOL
                 and row["max_history_enumerations"] == 0
                 and row["sample_ok"] for row in report["rows"])
        p3 = all(row["accepted_vs_target_tv"] <= TOL for row in report["rows"])
        c1 = primary["proposal_vs_target_tv"] > 1e-5
        c2 = primary["frozen_gamma_vs_target_tv"] > 1e-5
        c3 = primary["incoherent_vs_target_tv"] > 1e-5
        collision = collision_control(SparseCoherentReverse)
        c4 = (collision["cancellation"]
              and collision["coherent_vs_probability_tv"] > 0.4)
        report["controls"] = {
            "omit_acceptance_tv": primary["proposal_vs_target_tv"],
            "freeze_gamma_tv": primary["frozen_gamma_vs_target_tv"],
            "incoherent_terminal_tv": primary["incoherent_vs_target_tv"],
            "collision": collision,
            "all_controls_fail": bool(c1 and c2 and c3 and c4),
        }
        exp.check("P1", p1, "accepted submass and independent/full-r laws")
        exp.check("P2", p2, "proposal, terminal identities, and charged counters")
        exp.check("P3", p3, "endpoint, fixed-point, wrap, collision, and b=2 rows")
        exp.fail_check("C1", c1, f"omit-acceptance TV={primary['proposal_vs_target_tv']:.6g}")
        exp.fail_check("C2", c2, f"freeze-gamma TV={primary['frozen_gamma_vs_target_tv']:.6g}")
        exp.fail_check("C3", c3, f"incoherent-terminal TV={primary['incoherent_vs_target_tv']:.6g}")
        exp.fail_check("C4", c4, f"collision TV={collision['coherent_vs_probability_tv']:.6g}")
        report["status"] = "PASS" if all((p1, p2, p3, c1, c2, c3, c4)) else "FAIL"
        path = stamp()
        if not exp.finish(report_path=path, rows=report["rows"], metadata={
                "allocation_bound_bytes": MAX_DENSE_BYTES,
                "tolerance": TOL,
                "method": "SparseCoherentReverse vs joint_probability and direct_joint",
                "controls": report["controls"],
                "counter_recurrence": "support expansion at each descending insertion; 3 matvecs and 1 QFT pair per reached sector, then one final W0 matvec",
                "elapsed_seconds": time.perf_counter() - started}):
            raise AssertionError("experiment harness failed")
    except BaseException as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
        path = stamp("coherent_reverse_merge_failure")
        path.write_text(json.dumps(report, indent=2, default=str) + "\n")
        print(json.dumps({"status": "FAIL", "report": str(path)}))
        raise
    print(json.dumps({"status": report["status"], "report": str(path)}))


if __name__ == "__main__":
    main()
