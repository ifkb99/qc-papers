"""Does reachable compiled arithmetic buy more than removing scratch space?

DERIVATION BEFORE MEASUREMENT (TODO 14).
After s descending-power controls, every work label is a^(2^(t-s)*j) mod N
for 0 <= j < 2^s. Thus even the premeasurement union contains at most
min(2^s, r/gcd(r,2^(t-s))) distinct labels. Measurement only changes their
amplitudes. All intermediate classical gates preserve the number of labels
within a branch. No order is required to discover these labels on demand.
This eliminates exponential scratch-space tables but is O(r), not sub-orbit
compression. The coherent pre-QFT state has exponent/work Schmidt rank
min(r,2^t), since its distinct orbit labels correspond to disjoint, nonempty
exponent classes. A scalar output's smaller tensor rank is not that rank.

P1: forced full small distributions, sampled wide paths, and forced rare
outputs match independent FFT/geometric-series references. Relative error
on positive rare outputs at t<=32 must be <1e-6 (absolute <1e-12 elsewhere).
P2: all per-step union sizes obey the derived reachable-set bound; compiled
setup evaluates ZERO basis images, and every reached block boundary is clean.
P3: the exact coherent Schmidt rank is min(r,Q), independent of scalar rank.
P4 (added before the second run): once r has been discovered, a uniform latent
eigenphase k/r followed by scalar Bernoulli updates gives the SAME marginal
distribution with no work vector. Charge O(r) return-to-one order discovery;
this uses constant many growing integers, not a stored orbit. Check by exact
mixture enumeration at small width. Conditional-on-k probabilities must NOT
be compared to the marginal geometric-series formula. This strengthens the
baseline after the first run found ordinary orbit enumeration faster.
H1: runtime advantage beyond a scratch-dense baseline is OPEN. Charge compiled
setup, cold replay, warm sampling, transition cache, and ordinary orbit setup
and sampling separately. There is no predicted win over orbit enumeration.
C1: omitted feedback fails on odd-factor arithmetic. C2: scalar rank does not
bound joint rank. C3: growing order defeats any universal constant-work-memory
interpretation of this sparse representation. C4: noncommuting actions still
defeat naive coherent-to-iterative reordering.

Two one-parameter series: fixed N=21,a=2 with t=2..8,16,24,32; then fixed
a=2,t=32 and N=15,21,35,77,143,221,323,437,667. Order varies with N in the
second series and is recorded as a consequence, not independently controlled.
Small validation uses N=7,a=3. Timings are local engineering measurements,
not competitive MPS or decision-diagram benchmarks.

Run: OPENBLAS_NUM_THREADS=1 uv run --no-project --python 3.12 --with 'numpy<2.5' python -m experiments.experiment_reachable_order_finding
"""
from __future__ import annotations
import math
import time
import numpy as np
from lab import Experiment
from lab.reachable import compiled_pairs, sparse_path, action_stats
from lab.semiclassical import sample, order_finding_probability, eigenphase_path
from experiments.experiment_conditional_order_finding import (
    arithmetic_pairs, fft_reference, circuit_reference)
from lab.tensors import rank_profile
from circuits import Circuit
import statevec
import accel


def orbit_baseline(N, a, width):
    """Ordinary orbit enumeration, charged as baseline setup, not free input."""
    orbit, value = [], 1
    while not orbit or value != 1:
        orbit.append(value)
        value = value * a % N
    lookup = {x: i for i, x in enumerate(orbit)}
    ids = np.arange(len(orbit))
    multiplier, cache, pairs = a % N, {}, []
    for _ in range(width):
        if multiplier not in cache:
            cache[multiplier] = np.array([lookup[x*multiplier % N] for x in orbit])
        pairs.append((ids, cache[multiplier]))
        multiplier = multiplier * multiplier % N
    initial = np.zeros(len(orbit), complex)
    initial[0] = 1
    return orbit, pairs, initial


def modular_oracle_pairs(N, a, width):
    """Cheap high-level arithmetic comparator; not the compiled replay route."""
    pairs, multiplier = [], a % N
    for _ in range(width):
        pairs.append((lambda ids: ids.copy(), lambda ids, c=multiplier: ids*c % N))
        multiplier = multiplier * multiplier % N
    return pairs


def order_by_return(N, a):
    """Ordinary sequential order discovery, storing no orbit table."""
    if math.gcd(N, a) != 1:
        raise ValueError("order requires a coprime base")
    value, period = a % N, 1
    while value != 1:
        value = value * a % N
        period += 1
    return period


def validate_blocks(me, pairs):
    # Fresh replay wrappers for the correctness gate: do NOT warm benchmark
    # action caches with all N valid inputs.
    from lab.reachable import CompiledBranch
    clean = np.arange(me.N, dtype=np.int64) << me.x[0]
    checked = set()
    multiplier = me.a
    for pair in pairs:
        if id(pair[0]) not in checked:
            for bit, action in enumerate(pair):
                gated = CompiledBranch(action.circuit, action.control, bit)(clean)
                expected = clean if bit == 0 else ((np.arange(me.N)*multiplier % me.N) << me.x[0])
                if not np.array_equal(gated, expected):
                    raise AssertionError("compiled modular block failed clean-code verification")
            checked.add(id(pair[0]))
        multiplier = multiplier * multiplier % me.N


def slim(path):
    return {k: v for k, v in path.items() if k not in ("final_ids", "final_amplitudes")}


def benchmark(exp, N, a, width, series):
    start = time.perf_counter()
    me, pairs = compiled_pairs(N, a, width)
    compile_seconds = time.perf_counter()-start
    exp.check("P2", action_stats(pairs)["cached_transitions"] == 0,
              f"N={N},t={width}: setup compiled circuits without basis enumeration")
    validate_blocks(me, pairs)
    initial = {1 << me.x[0]: 1.}
    start = time.perf_counter()
    cold = sparse_path(pairs, initial, rng=np.random.default_rng(11))
    cold_seconds = time.perf_counter()-start
    cold_stats = action_stats(pairs)
    rng = np.random.default_rng(42)
    start = time.perf_counter()
    paths = [sparse_path(pairs, initial, rng=rng) for _ in range(8)]
    warm_seconds = time.perf_counter()-start
    sample_stats = action_stats(pairs)

    start = time.perf_counter()
    orbit, orbit_pairs, psi = orbit_baseline(N, a, width)
    orbit_setup = time.perf_counter()-start
    period = len(orbit)  # reference/baseline only; NEVER given to compiled sampler
    rng = np.random.default_rng(42)
    start = time.perf_counter()
    orbit_paths = [sample(orbit_pairs, psi, rng, validated=True) for _ in range(8)]
    orbit_seconds = time.perf_counter()-start

    start = time.perf_counter()
    oracle_pairs = modular_oracle_pairs(N, a, width)
    oracle_setup = time.perf_counter()-start
    rng = np.random.default_rng(42)
    start = time.perf_counter()
    oracle_paths = [sparse_path(oracle_pairs, {1: 1.}, rng=rng) for _ in range(8)]
    oracle_seconds = time.perf_counter()-start
    start = time.perf_counter()
    discovered_period = order_by_return(N, a)
    phase_setup = time.perf_counter()-start
    assert discovered_period == period
    rng = np.random.default_rng(42)
    start = time.perf_counter()
    phase_paths = [eigenphase_path(discovered_period, width,
                                  int(rng.integers(discovered_period)), rng=rng)
                   for _ in range(8)]
    phase_seconds = time.perf_counter()-start
    errors = [abs(p["path_probability"]-order_finding_probability(p["output"], width, period))
              for p in [cold] + paths]
    exp.check("P1", max(errors) < 1e-12 and all(
        p["output"] == q["output"] == v["output"]
        and abs(p["path_probability"]-q["path_probability"]) < 1e-12
        and abs(p["path_probability"]-v["path_probability"]) < 1e-12
        for p, q, v in zip(paths, orbit_paths, oracle_paths)),
        f"N={N},t={width}: compiled/orbit/oracle paths agree, error={max(errors):.2e}")

    Q = 1 << width
    rare = []
    for output in sorted({1, Q//5 + 1, Q//2 - 1}):
        path = sparse_path(pairs, initial, output=output)
        reference = order_finding_probability(output, width, period)
        error = abs(path["path_probability"]-reference)
        relative = error/reference if reference > 0 else None
        exp.check("P1", relative < 1e-6 if relative is not None else error < 1e-25,
                  f"N={N},t={width},forced y={output}: p={reference:.3e}, rel={relative}")
        rare.append(dict(output=output, reference=reference, absolute_error=error,
                         relative_error=relative, path=slim(path)))
    all_paths = [cold] + paths
    bounds_ok = all(row["union"] <= min(1 << (row["step"]+1),
                         period//math.gcd(period, 1 << row["exponent_bit"]))
                    for path in all_paths for row in path["profile"])
    # Every memoized transition was reached by a sample or forced path. Check
    # both input and output clean labels, not only final state support.
    clean_mask = ((1 << me.n) - 1) << me.x[0]
    actions = {id(action): action for pair in pairs for action in pair}.values()
    clean_ok = all((x & ~clean_mask) == 0 and (x >> me.x[0]) < N
                   for action in actions for kv in action.cache.items() for x in kv)
    exp.check("P2", bounds_ok and clean_ok,
              f"N={N},t={width}: reachable subgroup bound and clean boundaries hold")
    return dict(series=series, N=N, a=a, width=width, period=period,
                work_qubits=me.exp[0], dense_work_amplitudes=1 << me.exp[0],
                compile_seconds=compile_seconds, cold_sample_seconds=cold_seconds,
                warm_eight_samples_seconds=warm_seconds, cold_stats=cold_stats,
                sample_stats=sample_stats, final_stats=action_stats(pairs),
                orbit_setup_seconds=orbit_setup, orbit_eight_samples_seconds=orbit_seconds,
                orbit_stored_amplitudes=len(orbit), oracle_setup_seconds=oracle_setup,
                oracle_eight_samples_seconds=oracle_seconds, max_error=max(errors),
                phase_order_discovery_seconds=phase_setup,
                phase_eight_samples_seconds=phase_seconds, phase_paths=phase_paths,
                cold_path=slim(cold), paths=[slim(p) for p in paths], rare=rare)


def main():
    exp = Experiment("reachable_order_finding", doc=__doc__)
    exp.predict("P1", "small distributions, sampled paths and rare probabilities agree")
    exp.predict("P2", "no full-space setup; clean reachable unions obey subgroup bound")
    exp.predict("P3", "coherent joint Schmidt rank is min(r,Q), not scalar output rank")
    exp.predict("P4", "order-informed latent-eigenphase mixture reproduces marginal outputs")
    exp.predict("H1", "open timing comparison, with setup and sampling both charged")
    exp.must_fail("C1", "missing phase feedback changes odd-factor output statistics")
    exp.must_fail("C2", "small scalar rank does not bound joint-state rank")
    exp.must_fail("C3", "growing order defeats a constant-memory interpretation")
    exp.must_fail("C4", "noncommuting arithmetic invalidates naive reordering")
    rows = []
    me, pairs = compiled_pairs(7, 3, 8)
    validate_blocks(me, pairs)
    for width in (3, 4, 6, 8):
        probabilities = np.array([sparse_path(pairs[:width], {1 << me.x[0]: 1.}, output=y)
                                  ["path_probability"] for y in range(1 << width)])
        reference = fft_reference(7, 3, width)
        error = float(np.max(np.abs(probabilities-reference)))
        exp.check("P1", error < 1e-12 and abs(probabilities.sum()-1) < 1e-12,
                  f"all {1 << width} outputs at t={width}: FFT error={error:.2e}")
        if width == 3:
            exp.check("P1", np.max(np.abs(probabilities-circuit_reference(7, 3, width))) < 1e-9,
                      "independent Fourier-compiled full-state reference")
        if width == 4:
            wrong = np.array([sparse_path(pairs[:width], {1 << me.x[0]: 1.}, output=y,
                                          feedback=False)["path_probability"] for y in range(16)])
            exp.fail_check("C1", np.max(np.abs(wrong-reference)) > .01,
                           f"missing feedback error={np.max(np.abs(wrong-reference)):.6g}")
        rows.append(dict(series="full_distribution", N=7, a=3, width=width, error=error))

    # Same existing dense coherent simulator for the noncommuting gate control.
    ids = np.arange(4)
    actions = [(lambda x: x.copy(), lambda x: x ^ ((x & 1) << 1)),
               (lambda x: x.copy(), lambda x: x ^ ((x >> 1) & 1))]
    iterative = np.array([sparse_path(actions, {1: 1.}, output=y)["path_probability"]
                          for y in range(4)])
    qc = Circuit(4).x(0).h(2).h(3).toffoli(2, 0, 1).toffoli(3, 1, 0)
    qc.qft([2, 3], inverse=True)
    coherent = np.sum(np.abs(statevec.run(qc).reshape(4, 4))**2, axis=1)
    exp.fail_check("C4", np.max(np.abs(iterative-coherent)) > .1,
                   f"noncommuting reordering error={np.max(np.abs(iterative-coherent)):.6g}")

    for N, a in [(7, 3), (13, 4), (15, 7)]:
        orbit, _, _ = orbit_baseline(N, a, 1)
        period = len(orbit)
        for width in (3, 4):
            mixture = np.array([sum(eigenphase_path(period, width, k, output=y)
                                    ["conditional_path_probability"] for k in range(period))/period
                                for y in range(1 << width)])
            mixture_error = float(np.max(np.abs(mixture-fft_reference(N, a, width))))
            exp.check("P4", mixture_error < 1e-12,
                      f"N={N},a={a},t={width}: latent eigenphase mixture error={mixture_error:.2e}")
            rows.append(dict(series="eigenphase_mixture", N=N, a=a, width=width,
                             period=period, max_error=mixture_error))
        for width in (1, 2, 3, 6):
            Q = 1 << width
            joint = np.zeros((Q, N))
            joint[np.arange(Q), [pow(a, e, N) for e in range(Q)]] = 1/np.sqrt(Q)
            singular = np.linalg.svd(joint, compute_uv=False)
            rank = int(np.count_nonzero(singular > 1e-12))
            exp.check("P3", rank == min(period, Q), f"N={N},t={width}: joint rank={rank}")
        if N == 13:
            values = np.array([(-1.)**(pow(a, e, N) & 1) for e in range(64)])
            scalar = rank_profile(values)
            scalar_rank = max(row["rank"] for row in scalar)
            exp.fail_check("C2", scalar_rank == 1 and rank == 6,
                           f"same arithmetic: scalar rank={scalar_rank}, joint rank={rank}")

    for width in (*range(2, 9), 16, 24, 32):
        rows.append(benchmark(exp, 21, 2, width, "width"))
    for N in (15, 21, 35, 77, 143, 221, 323, 437, 667):
        rows.append(benchmark(exp, N, 2, 32, "modulus"))

    modulus_rows = [row for row in rows if row["series"] == "modulus"]
    peaks = [max(p["peak_stored_amplitudes"] for p in row["paths"]) for row in modulus_rows]
    exp.fail_check("C3", peaks[-1] > 10*peaks[0], f"stored peaks as order varies: {peaks}")
    timed = [row for row in rows if row["series"] in ("width", "modulus")]
    exp.check("H1", all(row["compile_seconds"] > 0 and row["orbit_setup_seconds"] > 0
                        and row["phase_order_discovery_seconds"] > 0
                        for row in timed), "setup, cold replay, warm sampling and orbit costs recorded; no speedup prediction")

    # Dense full-scratch reference only where it is cheap; identical task/seed.
    start = time.perf_counter()
    dense_me, dense_pairs, psi, setup = arithmetic_pairs(7, 3, 8)
    rng = np.random.default_rng(42)
    dense_paths = [sample(dense_pairs, psi, rng, validated=True) for _ in range(8)]
    rows.append(dict(series="dense_scratch", N=7, a=3, width=8,
                     setup_seconds=setup, total_seconds=time.perf_counter()-start,
                     stored_amplitudes=psi.size, paths=dense_paths))
    rng = np.random.default_rng(42)
    sparse_paths = [sparse_path(pairs, {1 << me.x[0]: 1.}, rng=rng) for _ in range(8)]
    exp.check("P1", all(p["output"] == q["output"]
                         and abs(p["path_probability"]-q["path_probability"]) < 1e-12
                         for p, q in zip(sparse_paths, dense_paths)),
              "matched samples agree with full-scratch dense compiled baseline")
    exp.finish(report_path="out/reachable_order_finding.json", rows=rows,
               metadata=dict(numpy=np.__version__, gpu=accel.enabled(), samples_per_row=8,
                             amplitude_cutoff=None, precision="complex128",
                             timing="single local run; not a competitive simulator benchmark"))


if __name__ == "__main__":
    main()
