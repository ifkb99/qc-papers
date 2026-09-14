"""Boundary test for late-background cancellation beyond b=2.

The fixed comparison uses M=3,t=4 and paired indexed orbits r=6,b=2 and
r=9,b=3.  Each has W1 an embedded Rx(pi/7), W3 an embedded Rz(pi/5), and
coherent reflections K2(q=0,pi/5), K3(q=1,pi/5).  Only the joint change caused
by removing W1 and W3 together is measured; changing b also changes r, so this
is a boundary comparison rather than a width-scaling sweep.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The b=2 late-background removal leaves the complete (gamma,y) law
      unchanged, while the b=3 case is an input-specific visibility test and
      may differ; no universal odd-b claim is assumed.
  P2  The existing coherent route contraction agrees with the independent
      full-r product and all compared laws normalize.
  P3  A background placed only at the final insertion s=t is invisible for
      both b=2 and b=3 because work is traced immediately afterward.
  C1  If b=3 is visibly changed, the false assertion that all late backgrounds
      are removable must fail.

No new propagator or production orbit/random-word table is introduced.
The diagnostic DOES allocate dense full-r matrices, guarded at r<=14,t<=4.
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
from experiments.experiment_coherent_route_sampling import direct_joint, rx
from lab.coherent_routes import CoherentReflectionCircuit


MAX_BYTES = 16 << 20
WIDTH = 4
M = 3

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "b=2 late-background removal preserves the joint law; b=3 visibility is tested separately by C1")
exp.predict("P2", "existing route contraction agrees with independent full-r laws and normalizes")
exp.predict("P3", "a final-only background is invisible for both paired block sizes")
exp.must_fail("C1", "the universal claim that every late background can be removed")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def report_path():
    return Path("out") / f"odd_block_tail_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def embedded_rx(block, theta):
    out = np.eye(block, dtype=complex)
    out[:2, :2] = rx(theta)
    return out


def embedded_rz(block, theta):
    out = np.eye(block, dtype=complex)
    out[0, 0] = np.exp(-1j * theta / 2)
    out[1, 1] = np.exp(1j * theta / 2)
    return out


def fixture(block, *, late=True, final_only=False):
    period = M * block
    guard((period, period), label="full-r work matrix")
    guard((period, 1 << WIDTH), label="full-r output workspace")
    backgrounds = {}
    if late:
        backgrounds.update({1: embedded_rx(block, np.pi / 7),
                            3: embedded_rz(block, np.pi / 5)})
    if final_only:
        backgrounds = {4: embedded_rx(block, np.pi / 7)}
    return CoherentReflectionCircuit(
        period, block, WIDTH, backgrounds,
        {2: (0, np.pi / 5), 3: (1, np.pi / 5)})


def joint_matrix(circuit):
    return np.array([[circuit.joint_probability(gamma, output)
                      for output in range(1 << circuit.width)]
                     for gamma in range(circuit.sectors)], dtype=float)


def compare_pair(block):
    with_late = fixture(block, late=True)
    without_late = fixture(block, late=False)
    full_with = direct_joint(with_late)
    full_without = direct_joint(without_late)
    route_with = joint_matrix(with_late)
    route_without = joint_matrix(without_late)
    return dict(block_size=block, period=block * M,
                with_late_mass=float(full_with.sum()),
                without_late_mass=float(full_without.sum()),
                route_with_mass=float(route_with.sum()),
                route_without_mass=float(route_without.sum()),
                full_with_route_max_error=float(np.max(np.abs(full_with-route_with))),
                full_without_route_max_error=float(np.max(np.abs(full_without-route_without))),
                late_removal_max_abs_error=float(np.max(np.abs(full_with-full_without))),
                late_removal_tv=float(np.abs(full_with-full_without).sum() / 2),
                full_with=full_with.tolist(), full_without=full_without.tolist())


def final_only_compare(block):
    no_background = fixture(block, late=False)
    final_background = fixture(block, late=False, final_only=True)
    a = direct_joint(no_background)
    b = direct_joint(final_background)
    return dict(block_size=block, period=block * M,
                no_background_mass=float(a.sum()), final_background_mass=float(b.sum()),
                max_abs_error=float(np.max(np.abs(a-b))),
                tv=float(np.abs(a-b).sum() / 2))


def main():
    started = time.time()
    report = {"status": "PASS", "pairs": {}, "final_only": {}, "controls": {}}
    try:
        for block in (2, 3):
            report["pairs"][f"b{block}"] = compare_pair(block)
            report["final_only"][f"b{block}"] = final_only_compare(block)
        b2 = report["pairs"]["b2"]
        b3 = report["pairs"]["b3"]
        p1 = (abs(b2["with_late_mass"] - 1) < 1e-13
              and abs(b2["without_late_mass"] - 1) < 1e-13
              and b2["late_removal_max_abs_error"] < 2e-14)
        p2 = all(abs(row[key] - 1) < 1e-13
                  for row in report["pairs"].values()
                  for key in ("with_late_mass", "without_late_mass",
                              "route_with_mass", "route_without_mass")) and all(
            row[key] < 2e-14 for row in report["pairs"].values()
            for key in ("full_with_route_max_error", "full_without_route_max_error"))
        p3 = all(abs(row["no_background_mass"] - 1) < 1e-13
                 and abs(row["final_background_mass"] - 1) < 1e-13
                 and row["max_abs_error"] < 2e-14
                 for row in report["final_only"].values())
        # This is deliberately the universal shortcut control: it must fail
        # exactly when the b=3 fixture supplies a visible late-background
        # counterexample.  If b=3 is also invisible, retain that failed
        # prediction rather than silently changing the predicate.
        report["controls"]["universal_late_removal_tv_b3"] = b3["late_removal_tv"]
        report["controls"]["universal_claim_fails"] = b3["late_removal_tv"] > 1e-10
        exp.check("P1", p1, f"b2 late-removal max={b2['late_removal_max_abs_error']:.3g}; b3 TV={b3['late_removal_tv']:.8g}")
        exp.check("P2", p2, "full-r and existing route contraction masses/errors checked")
        exp.check("P3", p3, "final-only background removal checked for b=2 and b=3")
        exp.fail_check("C1", report["controls"]["universal_claim_fails"],
                       f"b3 late-removal TV={b3['late_removal_tv']:.8g}")
        path = report_path()
        if not exp.finish(report_path=path,
                          rows=list(report["pairs"].values()),
                          metadata=dict(final_only=report["final_only"],
                                        controls=report["controls"],
                                        reference="existing full-r complex128 product",
                                        max_reference_period=14, max_reference_width=4,
                                        per_array_payload_ceiling_bytes=MAX_BYTES)):
            raise AssertionError("Experiment harness failed")
        print(json.dumps(dict(ok=True, report=str(path))))
        return
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time() - started
    path = report_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(report), indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": str(path),
                      "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
