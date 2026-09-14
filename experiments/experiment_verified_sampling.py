"""Does the verified exact-input sampler realize the C60 law-level budget?

PREDICTIONS BEFORE MEASUREMENT (TEMPLATE.py workflow).
P1. On fixed r=10,b=2,t=4, exact enumeration of the actual dyadic-weight,
    finite-bit transition kernels normalizes and its full final law stays
    inside the rational planned TV budget at each requested accuracy.
P2. The verified final-amplitude balls agree with the existing independent
    full-r matrix-product reference (a float diagnostic, not the certificate).
P3. A supplied wide indexed instance needs no orbit/output arrays and reports
    its observed working precision, prefix count and bound <= requested TV.
C1. A deliberately capped evaluator cannot furnish the requested oracle;
    a successful ordinary-float sample alone must not carry our certificate.
C2. Coarse all-zero approximate blocks must retain probability through the
    declared fallback; dropping them fails normalization.

Only target TV varies in the primary sweep. Exact Fraction frontiers contain
at most M*2^t*b=160 unique labels after consumed exponent bits are removed.
Reference arrays are <32 KiB total payload and existing full-r guards apply.
The frontier/cache reference is capped at r<=10,t<=4 and 5000 queries; an
8 MiB conservative object allowance is checked before enumerating, plus
32 MiB peak traced-Python allocation cap after execution (not total RSS).
The wide run has no reference frontier/cache. Seeded Random supplies testing
bits only; ideal independent unbiased bits remain a mathematical assumption.
"""
from __future__ import annotations

import json
import time
import traceback
import tracemalloc
from fractions import Fraction
from random import Random

import numpy as np
from flint import arb, ctx, fmpq

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit
from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction, integer_cdf_counts
from experiments.experiment_coherent_route_sampling import direct_joint, rx, stamp


def fixtures():
    backgrounds = {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))}
    reflections = {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))}
    verified = VerifiedReflectionCircuit(10, 4, backgrounds, reflections)
    float_circuit = CoherentReflectionCircuit(10, 2, 4,
        {1: rx(np.pi/7), 3: np.diag([np.exp(-1j*np.pi/10), np.exp(1j*np.pi/10)])},
        {s: (q, float(angle)*np.pi) for s, (q, angle) in reflections.items()})
    return verified, float_circuit


def exact_transition_law(circuit, accuracy_bits, random_bits, *, discard_zeros=False):
    """Reference-only exact Markov enumeration for b=2/3; no sampled histogram."""
    M, t = circuit.sectors, circuit.width
    if circuit.period > 10 or t > 4 or circuit.b not in (2, 3):
        raise ValueError("reference enumeration exceeds fixed tiny budget")
    Q = 1 << t
    # Each prefix has M*Q unique (sector, exponent, measured-output) labels.
    capacity = M * Q * circuit.b
    query_capacity = (len(circuit.backgrounds) + len(circuit.reflections) + t) * M * Q
    if (query_capacity + 2*capacity) * 4096 > 8 << 20:
        raise MemoryError("preflight Python cache/frontier allowance exceeded")
    frontier = {(a, e, 0, 0): Fraction(1, M * Q) for a in range(M) for e in range(Q)}
    cache = {}
    zero_blocks, max_working, query_count, peak_states = 0, 0, 0, len(frontier)

    def query(a, e, stop, boundary="reflection", measured=0, output=0):
        nonlocal max_working, query_count
        e &= (1 << (t - measured)) - 1
        key = (a, e, stop, boundary, measured, output)
        if key not in cache:
            if len(cache) >= min(5000, query_capacity):
                raise MemoryError("reference query cache cap exceeded")
            result = circuit.prefix_dyadic(a, e, stop, boundary=boundary,
                measured=measured, output=output, accuracy_bits=accuracy_bits)
            cache[key] = result.weights
            max_working = max(max_working, result.working_precision)
            query_count += 1
        return cache[key]

    def add(out, key, mass):
        if mass:
            out[key] = out.get(key, Fraction(0)) + mass

    def update(kind, index):
        nonlocal frontier, zero_blocks, peak_states
        nxt = {}
        for key, mass in frontier.items():
            a, e, p, y = key
            if kind == "background":
                children = [(a, e, pp, y) for pp in range(circuit.b)]
                weights = query(a, e, index, "background")
            elif kind == "reflection":
                target = (-a - circuit.reflections[index][0]) % M
                if target == a:
                    add(nxt, key, mass)
                    continue
                children = [(aa, e, p, y) for aa in (a, target)]
                weights = [query(aa, e, index)[p] for aa in (a, target)]
            else:
                children = [(a, e & ((1 << (t-index))-1), p, yy)
                            for yy in (y, y | (1 << (index-1)))]
                weights = [query(a, e, t, measured=index, output=yy)[p]
                           for _, _, _, yy in children]
            if not any(weights):
                zero_blocks += 1
                if discard_zeros:
                    continue
            for child, count in zip(children, integer_cdf_counts(weights, random_bits)):
                add(nxt, child, mass * Fraction(count, 1 << random_bits))
        frontier = nxt
        peak_states = max(peak_states, len(frontier))
        if len(frontier) > capacity:
            raise AssertionError("unique-state frontier count exceeded preflight")
        if not discard_zeros and sum(frontier.values()) != 1:
            raise AssertionError("exact transition law lost probability")

    for stop in range(t + 1):
        if stop:
            bit, nxt = 1 << (stop - 1), {}
            for (a, e, p, y), mass in frontier.items():
                add(nxt, (a, e, (p + bit) % circuit.b if e & bit else p, y), mass)
            frontier = nxt
        if stop in circuit.backgrounds:
            update("background", stop)
        if stop in circuit.reflections:
            update("reflection", stop)
    for measured in range(1, t + 1):
        update("qft", measured)
    joint = {}
    for (a, _, p, y), mass in frontier.items():
        add(joint, (a, p, y), mass)
    return joint, dict(zero_blocks=zero_blocks, max_working_precision=max_working,
                      unique_prefix_queries=query_count, peak_frontier_states=peak_states)


def enclosed_tv(circuit, law):
    """Exact rational outward bound from ideal final balls, including work-p."""
    with ctx.workprec(160):
        total = arb(0)
        reference_marginal = np.zeros((circuit.sectors, 1 << circuit.width))
        for a in range(circuit.sectors):
            for y in range(1 << circuit.width):
                vector = circuit.prefix_enclosure(a, 0, circuit.width,
                    measured=circuit.width, output=y, working_precision=160)
                for p, z in enumerate(vector):
                    ideal = (z.real*z.real + z.imag*z.imag) / circuit.sectors
                    mass = law.get((a, p, y), Fraction(0))
                    total += abs(ideal - arb(fmpq(mass.numerator, mass.denominator)))
                    reference_marginal[a, y] += float(ideal.mid())
        return binary_fraction((total/2).upper()), reference_marginal


def main():
    exp = Experiment("verified_sampling", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "exact finite-bit full transition laws normalize and obey rational budgets")
    exp.predict("P2", "final verified balls match the independent float full-r reference")
    exp.predict("P3", "supplied-wide execution is table-free with charged accuracy and history work")
    exp.must_fail("C1", "a diagnostic precision cap/ordinary float result is not certification")
    exp.must_fail("C2", "discarding zero approximate blocks loses probability")
    started = time.perf_counter()
    circuit, comparator = fixtures()
    independent = direct_joint(comparator)
    rows = []
    tracemalloc.start()
    for target in (Fraction(1, 4), Fraction(1, 64), Fraction(1, 4096), Fraction(1, 1_000_000)):
        plan = circuit.accuracy_plan(target)
        law, stats = exact_transition_law(circuit, plan["absolute_accuracy_bits"], plan["random_bits"])
        tv_upper, ball_midpoints = enclosed_tv(circuit, law)
        discrepancy = float(np.max(np.abs(ball_midpoints - independent)))
        rows.append(dict(target_tv=str(target), planned_tv_upper=str(plan["total_tv_upper_bound"]),
                         observed_tv_upper=str(tv_upper), observed_tv_upper_float=float(tv_upper),
                         mass=str(sum(law.values())), accuracy_bits=plan["absolute_accuracy_bits"],
                         categorical_bits=plan["random_bits"], independent_max_error=discrepancy, **stats))
        exp.check("P1", sum(law.values()) == 1 and tv_upper <= plan["total_tv_upper_bound"] <= target,
                  f"target={target}, outward TV<={float(tv_upper):.6g}, planned<={plan['total_tv_upper_bound']}")
    _, traced_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    exp.check("P2", max(r["independent_max_error"] for r in rows) < 2e-14
              and abs(float(independent.sum())-1) < 2e-14,
              f"full-r diagnostic max difference={max(r['independent_max_error'] for r in rows):.3g}")

    wide = VerifiedReflectionCircuit((1 << 61)-2, 63,
        {1: ("x", Fraction(1, 7)), 31: ("z", Fraction(1, 5)), 63: ("x", Fraction(1, 7))},
        {2: (0, Fraction(1, 5)), 32: (1, Fraction(1, 5)), 63: (2, Fraction(1, 5))})
    wide_start = time.perf_counter()
    sample = wide.sample(Random(624), target_tv=Fraction(1, 1_000_000))
    sample["elapsed_seconds"] = time.perf_counter()-wide_start
    exp.check("P3", sample["total_tv_upper_bound"] <= Fraction(1, 1_000_000)
              and sample["prefix_vector_evaluations"] <= 3 + 2*3 + 2*63
              and not sample["orbit_or_output_tables"] and sample["precision_cap"] is None
              and traced_peak < 32 << 20,
              f"wide t63: {sample['prefix_vector_evaluations']} queries, "
              f"P={sample['max_working_precision']}, {sample['elapsed_seconds']:.3g}s; traced tiny peak={traced_peak}")

    capped = False
    try:
        circuit.prefix_dyadic(1, 0, 4, accuracy_bits=100, max_precision=64)
    except ArithmeticError:
        capped = True
    baseline = comparator.sample(np.random.default_rng(625))
    exp.fail_check("C1", capped and "total_tv_upper_bound" not in baseline,
                   f"cap raised={capped}; float sample has certificate={'total_tv_upper_bound' in baseline}")
    coarse, coarse_stats = exact_transition_law(circuit, 0, 16)
    dropped, dropped_stats = exact_transition_law(circuit, 0, 16, discard_zeros=True)
    exp.fail_check("C2", sum(coarse.values()) == 1 and coarse_stats["zero_blocks"] > 0
                   and sum(dropped.values()) < 1 and dropped_stats["zero_blocks"] > 0,
                   f"coarse fallback mass={sum(coarse.values())}, deleted mass={sum(dropped.values())}")
    path = stamp("verified_sampling")
    ok = exp.finish(report_path=path, rows=rows, metadata=dict(
        package="python-flint==0.9.0", fixture="r10,t4,W1 Rx(pi/7),W3 Rz(pi/5),K2 q0 pi/5,K3 q1 pi/5",
        wide={k: str(v) if isinstance(v, Fraction) else v for k,v in sample.items()},
        elapsed_seconds=time.perf_counter()-started, tiny_peak_traced_python_bytes=traced_peak,
        traced_peak_excludes_native_allocations=True, seeded_rng_is_not_unbiasedness_proof=True,
        no_same_accuracy_certified_rejection_comparator=True,
        coarse_fallback=coarse_stats, discarded_mass=str(sum(dropped.values()))))
    print(json.dumps(dict(report=str(path), ok=ok)))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    try:
        main()
    except BaseException as error:
        path = stamp("verified_sampling_failure")
        path.write_text(json.dumps(dict(error=repr(error), traceback=traceback.format_exc()), indent=2)+"\n")
        raise
