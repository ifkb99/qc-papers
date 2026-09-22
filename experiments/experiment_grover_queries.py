"""Structured Grover coherent-query validation and isolated resource pilot.

Run from research/ with python -m and single-thread BLAS environment.
Each invocation writes one new --report. Compare only successfully completed,
hash-verified case reports. See notes/MC-memory-structure-pilots.md.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import resource
import signal
import statistics
import time
import tracemalloc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("validate", "case", "compare"), required=True)
    parser.add_argument("--reports-dir", type=Path)
    parser.add_argument("--method", choices=("compact", "streaming", "dense"))
    parser.add_argument("--mode", choices=("time", "allocation"))
    parser.add_argument("--n", type=int)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if os.environ.get(key) != "1":
            raise ValueError(f"set {key}=1 before launching")
    if args.report.exists():
        raise FileExistsError(args.report)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    # Resource runs need an external RSS guard; virtual address size is not RSS.
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    signal.alarm(60)
    import numpy as np  # common warmed imports; excluded from timed region
    import walsh
    from experiments.grover_queries import compact, dense, streaming, oracle, query_masks
    from lab.harness import Experiment

    if args.stage == "compare":
        from fractions import Fraction
        if args.reports_dir is None:
            raise ValueError("compare requires --reports-dir")
        reference = json.loads((args.reports_dir / "reference_v1.json").read_text())
        if not (reference.get("ok") is True and reference.get("name") == "grover_reference_validation_v1"
                and not reference.get("warnings") and reference.get("checks")
                and all(row.get("ok") is True for row in reference["checks"])):
            raise ValueError("missing or unsuccessful validated reference report")
        exp = Experiment("grover_case_comparison_v1", doc=__doc__)
        exp.predict("P1", "All isolated methods/modes agree on the same six probabilities; rational routes agree exactly.")
        exp.predict("P2", "Every case respects the registered domain, finite-output and process-resource contract.")
        exp.must_fail("C1", "The omitted-oracle output is incompatible with actual returned p(1).")
        exp.must_fail("C2", "The dephased output is incompatible with actual returned p(0).")
        rows = []
        for n in range(12, 21):
            reports = {}
            protocol_ok = True
            for method in ("compact", "streaming", "dense"):
                for mode in ("time", "allocation"):
                    path = args.reports_dir / f"{method}_n{n}_{mode}_v1.json"
                    report = json.loads(path.read_text())
                    reports[method, mode] = report
                    protocol_ok &= (report.get("method") == method and report.get("mode") == mode
                                    and report.get("n") == n and report.get("steps") == 3
                                    and report.get("masks") == list(query_masks(n))
                                    and report.get("preconditions_valid") is True
                                    and report.get("comparison_pending") is True
                                    and report["raw_process_peak_rss_bytes"] <= (1 << 30)
                                    and len(report["probabilities"]) == len(query_masks(n))
                                    and all(math.isfinite(x) and 0 <= x <= 1+1e-11 for x in report["probabilities"]))
                    if mode == "time":
                        protocol_ok &= len(report["wall_seconds"]) == 5 and all(0 < t < 60 for t in report["wall_seconds"])
                    else:
                        protocol_ok &= isinstance(report["traced_peak_bytes"], int) and report["traced_peak_bytes"] >= 0
            exact = tuple(Fraction(x) for x in reports["compact", "time"]["exact_probabilities"])
            agrees = True
            for (method, mode), report in reports.items():
                agrees &= max(abs(float(x)-y) for x, y in zip(exact, report["probabilities"])) <= 1e-11
                if method != "dense":
                    agrees &= tuple(Fraction(x) for x in report["exact_probabilities"]) == exact
            exp.check("P1", agrees, f"n={n}")
            exp.check("P2", protocol_ok, f"n={n}")
            exp.fail_check("C1", exact[query_masks(n).index(1)] != 0, f"n={n}")
            exp.fail_check("C2", exact[0] != Fraction(1, 1 << n), f"n={n}")
            # Cost is tabulated only for cases whose same-output checks passed.
            if agrees and protocol_ok:
                for baseline in ("dense", "streaming"):
                    ratio = reports["compact", "time"]["median_wall_seconds"] / reports[baseline, "time"]["median_wall_seconds"]
                    smaller = reports["compact", "allocation"]["traced_peak_bytes"] < reports[baseline, "allocation"]["traced_peak_bytes"]
                    rows.append(dict(n=n, baseline=baseline, wall_ratio=ratio,
                                     smaller_traced_peak=smaller, provisional_useful=bool(smaller and ratio <= 2),
                                     strongest_structured_baseline="same dynamic program; no superiority claim"))
        exp.finish(report_path=args.report, rows=rows, metadata={
            "comparison": "same output and accuracy gates precede cost conclusions",
            "metrics": "untraced median callable wall time and separate traced live allocation; raw RSS retained in input reports",
            "reference": str(args.reports_dir / "reference_v1.json")})
        return

    if args.stage == "validate":
        exp = Experiment("grover_reference_validation_v1", doc=__doc__)
        exp.predict("P1", "All exact compact and streaming probabilities match independent vector updates within 1e-11.")
        exp.predict("P2", "Each complete tiny distribution sums to one within 1e-11, including zero/eight iterations.")
        exp.must_fail("C1", "Omitting the oracle destroys predicted nonzero y=1 probability at all planned widths.")
        exp.must_fail("C2", "Dephasing before the last H layer produces uniform output and disagrees at y=0.")
        exp.must_fail("C3", "Omitting the adjacency clause (0,1) is caught by input x=3.")
        rows = []
        for n in range(2, 9):
            for steps in (0, 1, 3, 8):
                masks = tuple(range(1 << n))
                c, d, s = compact(n, steps, masks), dense(n, steps, masks), streaming(n, steps, masks, block=7)
                error = max(abs(float(x)-y) for x, y in zip(c, d))
                exp.check("P1", c == s and error <= 1e-11, f"n={n}; steps={steps}; max_abs_error={error}")
                exp.check("P2", sum(c) == 1 and abs(sum(d)-1) <= 1e-11, f"n={n}; steps={steps}")
                rows.append(dict(n=n, steps=steps, max_absolute_error=error))
        for n in range(12, 21):
            c = compact(n, 3, (0, 1))
            exp.fail_check("C1", c[1] != 0, f"n={n}; omission is exactly zero at mask 1")
            exp.fail_check("C2", c[0] != 1/(1 << n), f"n={n}; dephased oracle output is exactly uniform")
            exp.fail_check("C3", not oracle(n, 3), f"n={n}; defective clause-free predicate accepts 3")
        exp.finish(report_path=args.report, rows=rows, metadata={
            "scope": "independent vector route shares NumPy and walsh.wht only; author's reference is not blind",
            "tolerance": 1e-11, "scalar_domain": "n=2..22, steps=0..8, in-range integer masks"})
        return

    if args.n not in range(12, 21) or args.method is None or args.mode is None:
        raise ValueError("case requires n=12..20, method and mode")
    funcs = dict(compact=compact, streaming=streaming, dense=dense)
    fn = funcs[args.method]
    # Common output specification construction included in each measured call.
    def invoke():
        return fn(args.n, 3, query_masks(args.n))
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    walls = []
    traced = None
    probabilities = None
    if args.mode == "allocation":
        tracemalloc.start()
        probabilities = invoke()
        _, traced = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    else:
        for _ in range(5):
            probabilities = None  # do not retain prior output inside next trial
            start = time.perf_counter_ns()
            probabilities = invoke()
            walls.append((time.perf_counter_ns() - start) * 1e-9)
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    values = [float(x) for x in probabilities]
    # This is a single-arm resource report, NOT a correctness verdict.
    report = dict(schema_version=1, kind="single-arm resource observation",
                  method=args.method, n=args.n, steps=3, mode=args.mode,
                  masks=list(query_masks(args.n)), probabilities=values,
                  exact_probabilities=([str(x) for x in probabilities]
                                       if args.method != "dense" else None),
                  wall_seconds=walls, median_wall_seconds=(statistics.median(walls) if walls else None),
                  traced_peak_bytes=traced, raw_process_peak_rss_bytes=rss_after,
                  raw_import_peak_rss_bytes=rss_before,
                  rss_definition="Linux ru_maxrss; whole-process high-water mark, includes imports; difference is not total allocation",
                  allocation_definition="tracemalloc live-allocation peak during one complete call; separate process from untraced time",
                  comparison_pending=True, preconditions_valid=(rss_after <= (1 << 30) and all(math.isfinite(p) and 0 <= p <= 1+1e-11 for p in values)))
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, allow_nan=False))


if __name__ == "__main__":
    main()
