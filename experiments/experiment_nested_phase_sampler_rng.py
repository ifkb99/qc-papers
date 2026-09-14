"""Complete weighted traversal of actual nested-phase sampler decisions.

PREDICTIONS, WRITTEN BEFORE MEASURING.
P1: For r3,b1,Q8,s2 and phases at0,1, literal columns are
    (-1)^e |e mod3>. Actual sample() paths at cuts0,2, with strides3,6,
    reconstruct their FFT joint law under fixed-work retries in both modes.
P2: A forced rejected attempt followed by acceptance draws work exactly once;
    finite exhaustion retains the work of every executed attempt.
P3: An independently specified unequal-mass conditional row exercises actual
    mass/root_mass component, interval-bit, lift and accept/reject paths.
C1: Restarting work after rejection biases its known (3,3,2)/8 marginal.

This enumerates weighted decisions, not uniform finite-RNG frequencies or a
precision certificate. The unequal-mass hook is NOT a physical joint circuit.
All bounds are frozen before execution; failed reports remain in out/.
"""
from __future__ import annotations
import math
import traceback
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from lab import Experiment
from lab.work_first import NestedPhaseProgressions

Q = 8
MAX_WORK = 1_000_000
MAX_BYTES = 1 << 20


class ScriptedRNG:
    def __init__(self, integers, reals):
        self.ints, self.reals = list(integers), list(reals)
        self.integer_bounds, self.real_calls = [], 0
    def integers(self, high):
        self.integer_bounds.append(int(high))
        value = self.ints.pop(0)
        assert 0 <= value < high
        return value
    def random(self):
        self.real_calls += 1
        value = self.reals.pop(0)
        assert 0 <= value < 1
        return value
    def exhausted(self):
        assert not self.ints and not self.reals


class Budget:
    def __init__(self):
        self.reserved = self.observed_checksum = 0
        self.observed = defaultdict(int)
        self.reference = defaultdict(int)
        self.calls = 0
    def reserve(self, amount):
        assert self.reserved+amount <= MAX_WORK
        self.reserved += amount
    def collect(self, counters, cap):
        amount = sum(counters.values())
        assert amount <= cap, (amount, cap, counters)
        self.observed_checksum += amount
        for key, value in counters.items():
            self.observed[key] += value
        self.calls += 1


def run_path(model, rng, budget, *, mode, accept, attempts=1, row=None):
    cap = 512 if attempts == 1 else 1024
    budget.reserve(cap)
    budget.reference["scripted_rng_values"] += len(rng.ints)+len(rng.reals)
    counters = model._new_counters() if row is not None else None
    result = None
    try:
        if row is None:
            result = model.sample(rng, max_attempts=attempts, proposal=mode)
        else:
            result = model._sample_row(row, rng, counters,
                                       max_attempts=attempts, proposal=mode)
        assert accept, "a prescribed rejection unexpectedly returned"
    except RuntimeError as exc:
        assert not accept and "cap exhausted" in str(exc)
    finally:
        counters = model.last_counters if counters is None else counters
        if counters is not None:
            budget.collect(counters, cap)
    rng.exhausted()
    assert counters["acceptance_draws"] == attempts
    assert counters["work_draws"] == (1 if row is None else 0)
    return result, counters


def physical(budget):
    budget.reference["literal_assignments"] += Q
    budget.reference["fft_input_entries"] += Q*3
    literal = np.zeros((Q, 3), complex)
    for e in range(Q):
        literal[e, e % 3] = (-1)**e
    target = abs(np.fft.fft(literal, axis=0).T)**2/Q**2
    work_mass = np.array([3., 3., 2.])/Q
    results = []
    for cut in (0, 2):
        for mode in ("mass", "root_mass"):
            model = NestedPhaseProgressions(3, 1, 3, 2, np.ones((1, 1)),
                np.ones((1, 1)), lambda j: 1,
                early_phases=((0, lambda j: 1j**j), (1, lambda j: (-1)**j)), cut=cut)
            budget.reference["setup_unitary_product_terms"] += 2
            budget.reference["setup_geometry_phase_terms"] += model.stats()["geometry_setup_phase_terms"]
            stride = 3 if cut == 0 else 6
            g = math.gcd(stride, Q)
            T, beta = Q//g, stride//g
            accepted, rejected, paths, retry_done = np.zeros((3, Q)), 0., 0, False
            for e in range(Q):
                j, m = e % 3, (3, 3, 2)[e % 3]
                for component in range(m):
                    for k in range(T):
                        bits = [.25 if not (k >> bit)&1 else .75 for bit in range(T.bit_length()-1)]
                        for lift in range(g):
                            y = (pow(beta, -1, T)*k) % T + T*lift
                            a = Q**2*target[j, y]/m**2
                            assert 0 <= a <= 1
                            weight = 1/(Q*m*Q)  # uniform exponent, component, bits, lift
                            accepted[j, y] += weight*a
                            rejected += weight*(1-a)
                            paths += 1
                            budget.reference["weighted_path_terms"] += 1
                            for branch in (True, False):
                                if (branch and a == 0) or (not branch and a == 1):
                                    continue
                                ints = [e]+([lift] if g > 1 else [])
                                rng = ScriptedRNG(ints, [.5, (component+.5)/m, *bits,
                                    a/2 if branch else (1+a)/2])
                                result, counters = run_path(model, rng, budget, mode=mode, accept=branch)
                                assert rng.integer_bounds == [Q]+([g] if g > 1 else [])
                                assert counters["marginal_queries"] == 2*(T.bit_length()-1)
                                if result:
                                    assert (result["work"], result["output"]) == (j, y)
                            if not retry_done and 0 < a < 1:
                                rng = ScriptedRNG([e]+([lift, lift] if g > 1 else []),
                                    [.5, (component+.5)/m, *bits, (1+a)/2,
                                     (component+.5)/m, *bits, a/2])
                                result, counters = run_path(model, rng, budget, mode=mode,
                                                            accept=True, attempts=2)
                                assert result["attempts"] == 2 and result["output"] == y
                                retry_done = True
            accepted_work = accepted.sum(axis=1)
            restored = accepted/accepted_work[:, None]*work_mass[:, None]
            error = float(np.max(abs(restored-target)))
            budget.reference["comparison_probability_entries"] += Q*3
            assert error < 2e-15 and retry_done
            assert np.max(abs(accepted_work-work_mass/np.array([3., 3., 2.]))) < 2e-15
            assert abs(accepted.sum()+rejected-1) < 2e-14
            bias = float(abs(accepted_work/accepted_work.sum()-work_mass).sum()/2)
            results.append(dict(cut=cut, mode=mode, stride=stride, gcd=g,
                weighted_paths=paths, joint_error=error, accepted_mass=float(accepted.sum()),
                rejected_mass=rejected, restart_work_tv=bias, actual_retry_checked=retry_done))
    return results


def unequal_row(budget):
    budget.reference["literal_assignments"] += 3
    budget.reference["fft_input_entries"] += Q
    budget.reference["setup_unitary_product_terms"] += 2
    model = NestedPhaseProgressions(3, 1, 3, 2, np.ones((1, 1)), np.ones((1, 1)),
                                   lambda j: 1, early_phases=(), cut=2)
    row = dict(work=0, stride=2, norm=1.,
               components=((0, 1, math.sqrt(.9)), (2, 2, math.sqrt(.05))))
    vector = np.zeros(Q, complex)
    vector[0], vector[2], vector[4] = math.sqrt(.9), math.sqrt(.05), math.sqrt(.05)
    target = abs(np.fft.fft(vector))**2/Q
    results = []
    for mode, weights, envelope in (("mass", (.9, .1), 2.), ("root_mass", (.75, .25), 1.6)):
        accepted, rejected, paths = np.zeros(Q), 0., 0
        for component in range(2):
            selection = weights[0]/2 if component == 0 else weights[0]+weights[1]/2
            reduced = (.25, .25, .25, .25) if component == 0 else (.5, .25, 0., .25)
            for k, probability in enumerate(reduced):
                if probability == 0:
                    continue
                bits = [.25 if k % 2 == 0 else .75,
                        .5 if component == 1 and k == 0 else (.25 if k < 2 else .75)]
                for lift in range(2):
                    y = k+4*lift
                    q = weights[0]/Q+weights[1]*(.5, .25, 0., .25)[k]/2
                    a = target[y]/(envelope*q)
                    assert 0 < a < 1
                    weight = weights[component]*probability/2
                    accepted[y] += weight*a
                    rejected += weight*(1-a)
                    paths += 1
                    budget.reference["weighted_path_terms"] += 1
                    for branch in (True, False):
                        rng = ScriptedRNG([lift], [selection, *bits, a/2 if branch else (1+a)/2])
                        result, counters = run_path(model, rng, budget, mode=mode,
                                                    accept=branch, row=row)
                        assert rng.integer_bounds == [2] and rng.real_calls == 4
                        assert counters["weighted_ratio_divisions"] == (2 if mode == "root_mass" else 0)
                        if result:
                            assert result["output"] == y
                            assert abs(result["expected_attempts"]-envelope) < 2e-15
        error = float(np.max(abs(envelope*accepted-target)))
        budget.reference["comparison_probability_entries"] += Q
        assert error < 3e-16 and abs(accepted.sum()+rejected-1) < 2e-15
        results.append(dict(mode=mode, weighted_paths=paths, conditional_error=error,
                            envelope=envelope, accepted_mass=float(accepted.sum())))
    return results


def main():
    exp = Experiment("nested_phase_sampler_rng", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "complete physical weighted RNG paths reconstruct the independent FFT joint law")
    exp.predict("P2", "actual sample retries hold work fixed; exhausted paths retain counters")
    exp.predict("P3", "actual unequal-mass conditional decisions reconstruct both proposal laws")
    exp.must_fail("C1", "restarting work after rejection fails to preserve its marginal")
    budget, result = Budget(), {}
    flags = [False]*4
    # At most 2 cuts * 2 modes * 8 seeds * 3 components * 8 outcomes * 2
    # branches, plus four retry calls and 56 conditional calls. These are
    # reserves, distinct from the subsequently counted actual operations.
    preflight = dict(max_work=MAX_WORK, max_bytes=MAX_BYTES,
        reserved_work=1536*512+4*1024+56*512+25_000,
        numeric_payload_bytes=512 << 10)
    assert preflight["reserved_work"] <= MAX_WORK and preflight["numeric_payload_bytes"] <= MAX_BYTES
    budget.reserve(25_000)  # independent reference/path/setup units, counted below
    try:
        result["physical"] = physical(budget)
        flags[0] = True
        flags[1] = all(row["actual_retry_checked"] for row in result["physical"])
        result["conditional_hook"] = unequal_row(budget)
        flags[2] = True
        flags[3] = all(row["restart_work_tv"] > .08 for row in result["physical"])
        assert sum(budget.reference.values()) <= 25_000
    except Exception:
        result["failure"] = traceback.format_exc()
    exp.check("P1", flags[0], "both strides, both proposal modes, all positive decision branches")
    exp.check("P2", flags[1], "full sample calls exercise rejection then success and finite exhaustion")
    exp.check("P3", flags[2], "separate unequal-mass conditional-row hook, compared coordinatewise to FFT")
    exp.fail_check("C1", flags[3], "restart-work TV exceeds .08 in the explicit physical fixture")
    result.update(preflight=preflight, reserved_used=budget.reserved,
                  observed_sampler_checksum=budget.observed_checksum,
                  observed_sampler_counters=dict(budget.observed), executed_sampler_calls=budget.calls)
    result["observed_reference_counters"] = dict(budget.reference)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out")/f"nested_phase_sampler_rng_{stamp}.json"
    ok = exp.finish(report_path=path, metadata=result)
    print("Report:", path)
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
