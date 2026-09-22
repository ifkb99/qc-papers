"""Research reproduction of the bounded TODO66 clean-block memory pilot.

P1 candidate and existing references give the same conditional probabilities.
P2 normalization error stays below 1e-10; probabilities agree within 1e-10.
C1 replacing coherent interference by a mixture changes the two-label witness.
C2 reversing feedback sign changes that witness with an imaginary overlap.

The accepted pilot validated the listed cases but found no complete-call memory
or runtime benefit over the logical dense baseline: certification dominated.
This is a research experiment, not a new default simulator or public lab API.
Details and measured scope: notes/MC-memory-structure-pilots.md.

Run validation first, then one backend/case per fresh capped process. Keep
timing and allocation modes separate. See the note for serialized commands.
Run: uv run python -m experiments.experiment_clean_block --help
The --reviewed switch retains the original manual design/reference gate.
"""
from __future__ import annotations
import argparse
import json
import math
import os
from pathlib import Path
import resource
import sys
import threading
import time
import tracemalloc

ROOT = next(p for p in Path(__file__).resolve().parents if (p/"CLAUDE.md").exists())
sys.path.insert(0, str(ROOT))
import numpy as np
from lab import Experiment
from lab.reachable import compiled_pairs, sparse_path
from lab.semiclassical import sample, order_finding_probability, conditional_children
from toffoli_arith import ToffoliModExp
import statevec
from experiments.clean_block import certify, arithmetic_pairs, multipliers, two_buffer_path


def controls(exp):
    # Exact witness <psi|Ppsi>=1/2 or -i/2, proved and gate-instantiated
    # in diagnostic_v1. This does not claim arbitrary phase mutants differ
    # on the plain |1> Shor input, whose distribution can be sign symmetric.
    for a in (2,3):
        N = 31
        x = np.arange(N, dtype=np.int64)
        pair = (x, a*x % N)
        psi = np.zeros(N, complex); psi[1] = psi[a] = 1/np.sqrt(2)
        plus, _ = conditional_children(psi, pair, 0.)
        coherent = float(np.vdot(plus, plus).real)
        mixture = float(np.vdot(psi,psi).real)/2
        exp.fail_check("C1", abs(coherent-.75)<1e-12 and abs(coherent-mixture)>.2,
                       f"a={a}: existing reference detects dephasing mutant")
        psi[a] *= 1j
        correct, _ = conditional_children(psi, pair, np.pi/2)
        mutant, _ = conditional_children(psi, pair, -np.pi/2)
        pc, pm = (float(np.vdot(v,v).real) for v in (correct,mutant))
        exp.fail_check("C2", abs(pc-.75)<1e-12 and abs(pm-.25)<1e-12,
                       f"a={a}: existing reference detects phase-sign mutant")


def validate(exp):
    rows = []
    # Actual Pauli-rotation interpreter is independent of classical_images.
    N, a, width = 7, 2, 3
    me = ToffoliModExp(N, a, n_exp=width)
    final = statevec.run(me.build_shor())
    probabilities = np.sum(np.abs(final.reshape(1 << width, -1))**2, axis=1)
    pairs, _ = certify(N, a, width)
    initial = np.zeros(N, complex); initial[1] = 1
    got = np.array([two_buffer_path(pairs, initial, output=y)["path_probability"]
                    for y in range(1 << width)])
    exp.check("P1", np.max(np.abs(got-probabilities)) < 1e-10,
              "N=7,a=2,width=3: actual full gate statevector")
    exp.check("P2", abs(float(got.sum())-1) < 1e-10, "candidate distribution norm")
    rows.append(dict(N=N, a=a, width=width,
                     reference="statevec.run(actual build_shor)",
                     max_absolute_error=float(np.max(np.abs(got-probabilities)))))
    for a in (2, 3):
        N = 31
        # Reference-only classical order discovery is charged to validation;
        # neither candidate nor measured baselines receive the period.
        period, value = 1, a
        while value != 1:
            value = value*a % N; period += 1
        for width in range(3, 9):
            pairs, _ = certify(N, a, width)
            actions = arithmetic_pairs(N, a, width)
            initial = np.zeros(N, complex); initial[1] = 1
            err = normerr = stateerr = 0.
            for y in range(1 << width):
                r = two_buffer_path(pairs, initial, output=y, return_state=True)
                s = sparse_path(actions, {1: 1.+0j}, output=y)
                ref = order_finding_probability(y, width, period)
                err = max(err, abs(r["path_probability"]-s["path_probability"]),
                          abs(r["path_probability"]-ref))
                normerr = max(normerr, r["max_normalization_error"])
                if r["state"] is not None and s["path_probability"] > 1e-12:
                    v = np.zeros(N, complex)
                    v[s["final_ids"]] = [complex(*z) for z in s["final_amplitudes"]]
                    stateerr = max(stateerr, float(np.linalg.norm(r["state"]-v)))
            exp.check("P1", err < 1e-10 and stateerr < 1e-8,
                      f"N={N},a={a},width={width}: probability/state reference")
            exp.check("P2", normerr < 1e-10, "normalization")
            rows.append(dict(N=N,a=a,width=width,max_absolute_error=err,
                             max_conditional_state_error=stateerr,
                             max_normalization_error=normerr))
    return rows


def measure(args, exp):
    N, a, width = 31, args.a, args.width
    if args.instrument == "allocation":
        tracemalloc.start()
    start = time.perf_counter()
    initial = None
    if args.backend in ("candidate", "logical-dense"):
        initial = np.zeros(N, complex); initial[1] = 1
    if args.backend == "candidate":
        pairs, cert = certify(N, a, width)
    elif args.backend == "logical-dense":
        x = np.arange(N, dtype=np.int64)
        maps = {m: m*x % N for m in set(multipliers(N, a, width))}
        pairs = [(x, maps[m]) for m in multipliers(N, a, width)]
        del maps
        cert = []
    elif args.backend == "logical-sparse":
        pairs, cert = arithmetic_pairs(N, a, width), []
    else:
        me, pairs = compiled_pairs(N, a, width)
        cert = []
    setup_seconds = time.perf_counter()-start
    setup_peak = None
    if args.instrument == "allocation":
        _, setup_peak = tracemalloc.get_traced_memory()
        tracemalloc.reset_peak()
    rng = np.random.default_rng(args.seed)
    outputs, errors = [], []
    sample_start = time.perf_counter()
    for _ in range(args.samples):
        if args.backend == "candidate":
            result = two_buffer_path(pairs, initial, rng=rng)
        elif args.backend == "logical-dense":
            result = sample(pairs, initial, rng, validated=True)
        elif args.backend == "logical-sparse":
            result = sparse_path(pairs, {1: 1.+0j}, rng=rng)
        else:
            result = sparse_path(pairs, {1 << me.x[0]: 1.+0j}, rng=rng)
        outputs.append([result["output"], result["path_probability"]])
        errors.append(result.get("max_normalization_error", 0.))
    sampling_seconds = time.perf_counter()-sample_start
    hot_peak = None
    if args.instrument == "allocation":
        _, hot_peak = tracemalloc.get_traced_memory()
    total_seconds = time.perf_counter()-start
    if args.instrument == "allocation":
        tracemalloc.stop()
    # Period is used only AFTER the timed interval, never to choose transitions.
    period, value = 1, a
    while value != 1:
        value = value*a % N; period += 1
    err = max(abs(p-order_finding_probability(y, width, period)) for y,p in outputs)
    exp.check("P1", err < 1e-10, "all sampled paths match independent scalar reference")
    exp.check("P2", max(errors) < 1e-10, "candidate normalization when available")
    return [dict(N=N,a=a,width=width,backend=args.backend,samples=args.samples,
                 seed=args.seed,instrument=args.instrument,setup_seconds=setup_seconds,
                 sampling_seconds=sampling_seconds,total_seconds=total_seconds,
                 setup_tracemalloc_peak=setup_peak,hot_tracemalloc_peak=hot_peak,
                 total_tracemalloc_peak=(max(setup_peak,hot_peak)
                                         if setup_peak is not None else None),
                 rss_max_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 outputs=outputs,certificates=cert,max_probability_error=err,
                 note="Only timing-mode wall times decide slowdown; allocation-mode times are diagnostic. RSS includes imports; traced peaks exclude imports and include setup/input/output.")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("validate", "measure"))
    parser.add_argument("--reviewed", action="store_true")
    parser.add_argument("--backend", choices=("candidate","logical-dense","logical-sparse","compiled-sparse"), default="candidate")
    parser.add_argument("--a", type=int, choices=(2,3), default=3)
    parser.add_argument("--width", type=int, choices=range(3,9), default=6)
    parser.add_argument("--samples", type=int, choices=(1,128), default=128)
    parser.add_argument("--instrument", choices=("allocation","timing"), default="allocation")
    parser.add_argument("--seed", type=int, default=20260921)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if not args.reviewed:
        parser.error("candidate execution requires accepted design/reference dependencies")
    if args.report.exists():
        parser.error("refuse to overwrite evidence")
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    # CPU imports can map >2 GiB of virtual addresses while RSS is small.
    # RSS watchdog, not RLIMIT_AS: the check interval can permit overshoot.
    # Planned arrays at N=31 are independently tiny; this is a backstop.
    def rss_watchdog():
        while True:
            if 1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 2 << 30:
                print("RSS budget exceeded (2 GiB)", flush=True)
                os._exit(99)
            time.sleep(.02)
    threading.Thread(target=rss_watchdog, daemon=True).start()
    exp = Experiment("memory_block_pilot_v1", doc=__doc__)
    exp.predict("P1", "candidate path probabilities match independent references")
    exp.predict("P2", "normalization error below 1e-10")
    exp.must_fail("C1", "incoherent mixture disagrees on real-overlap witness")
    exp.must_fail("C2", "wrong phase sign disagrees on imaginary-overlap witness")
    controls(exp)
    rows = validate(exp) if args.mode == "validate" else measure(args,exp)
    exp.finish(report_path=args.report, rows=rows,
               metadata=dict(draft_version=1,mode=args.mode,
                             resource_hypothesis="open",dtype="complex128"))


if __name__ == "__main__":
    main()
