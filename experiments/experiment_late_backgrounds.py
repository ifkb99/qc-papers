"""Full-r audit of the b=2 late-background invisibility boundary.

For the frozen r=10, b=2, t=4 coherent fixture, the two backgrounds at
insertions s=1 and s=3 are removed together while reflections and gate order
are held fixed.  The existing independent full-r product computes the joint
(coarse-sector, QFT-output) law.  A separate full-r projected work vector
checks a within-block basis observable, since equality of the requested joint
law is not equality of the retained work state.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  Removing all W_s with s>=1 leaves the complete joint (gamma,y) law
      unchanged for this b=2 fixture, to full-r floating reference precision.
  P2  The final within-block basis probabilities can still change under the
      same removal, so the invariance is only for the traced output contract.
  C1  Adding an insertion-0 Rx(pi/4) to the same fixed reflection fixture
      changes the joint law; the late-background rule must not be extended to
      W_0.

This is a bounded algebra/reference audit, not a new propagator or a claim
about arbitrary block sizes.  Dense reference allocations are preflighted
under the existing r<=14,t<=4 cap.
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
from experiments.experiment_coherent_route_sampling import direct_joint, direct_prefix, rx


MAX_BYTES = 16 << 20
PERIOD = 10
BLOCK = 2
WIDTH = 4

exp = Experiment("late_backgrounds", doc=__doc__, exit_on_fail=False)
exp.predict("P1", "all late b=2 backgrounds are invisible to the complete joint output law")
exp.predict("P2", "a within-block work-basis observable can still distinguish the states")
exp.must_fail("C1", "an insertion-0 mixer changes the joint law in the same fixture")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def report_path():
    return Path("out") / f"late_backgrounds_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def frozen_backgrounds():
    return {
        1: rx(np.pi / 7),
        3: np.diag([np.exp(-1j * np.pi / 10), np.exp(1j * np.pi / 10)]),
    }


def frozen_reflections():
    return {2: (0, np.pi / 5), 3: (1, np.pi / 5)}


def make_circuit(backgrounds):
    guard((PERIOD, PERIOD), label="full-r work matrix budget")
    guard((PERIOD, 1 << WIDTH), label="full-r output workspace budget")
    return CoherentReflectionCircuit(PERIOD, BLOCK, WIDTH, backgrounds,
                                     frozen_reflections())


def full_joint_metrics(original, late_free):
    original_joint = direct_joint(original)
    late_free_joint = direct_joint(late_free)
    difference = np.abs(original_joint - late_free_joint)
    return dict(original_mass=float(original_joint.sum()),
                late_free_mass=float(late_free_joint.sum()),
                max_abs_error=float(difference.max()),
                tv=float(difference.sum() / 2),
                original_joint=original_joint.tolist(),
                late_free_joint=late_free_joint.tolist())


def within_block_metrics(original, late_free):
    max_vector_error = 0.0
    max_basis_probability_error = 0.0
    vector_location = None
    probability_location = None
    for gamma in range(original.sectors):
        for output in range(1 << original.width):
            # Existing direct_prefix is the independent full-r projected
            # vector at the terminating inverse-QFT output, not a new state
            # propagator.  Its two entries are the retained p=0,1 amplitudes.
            vector = direct_prefix(original, gamma, 0, original.width,
                                  "reflection", measured=original.width,
                                  output=output)
            late_vector = direct_prefix(late_free, gamma, 0, late_free.width,
                                        "reflection", measured=late_free.width,
                                        output=output)
            vector_error = float(np.max(np.abs(vector - late_vector)))
            # direct_prefix returns the sqrt(M)-scaled p-vector.  Divide its
            # squared-coordinate difference by M before calling this an
            # actual joint (gamma,y,p) probability difference.
            probability_error = float(np.max(np.abs(np.abs(vector) ** 2
                                                    - np.abs(late_vector) ** 2))
                                      / original.sectors)
            if vector_error > max_vector_error:
                max_vector_error = vector_error
                vector_location = (gamma, output)
            if probability_error > max_basis_probability_error:
                max_basis_probability_error = probability_error
                probability_location = (gamma, output)
    return dict(max_projected_vector_abs_error=max_vector_error,
                vector_location=vector_location,
                max_joint_within_block_basis_probability_error=max_basis_probability_error,
                scaled_vector_squared_error=max_basis_probability_error * original.sectors,
                basis_probability_location=probability_location,
                basis_observable_changes=max_basis_probability_error > 1e-10)


def insertion_zero_control(late_backgrounds):
    without = make_circuit(late_backgrounds)
    with_initial = CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH,
        {0: rx(np.pi / 4), **late_backgrounds}, frozen_reflections())
    a = direct_joint(without)
    b = direct_joint(with_initial)
    difference = np.abs(a - b)
    return dict(theta_over_pi=1 / 4,
                baseline_mass=float(a.sum()), initial_mix_mass=float(b.sum()),
                max_abs_error=float(difference.max()),
                tv=float(difference.sum() / 2),
                changes_joint=float(difference.sum() / 2) > 1e-8)


def main():
    started = time.time()
    report = {"status": "PASS", "predictions": {}, "control": {}}
    try:
        late_backgrounds = frozen_backgrounds()
        original = make_circuit(late_backgrounds)
        late_free = make_circuit({})
        report["predictions"]["P1_joint"] = full_joint_metrics(original, late_free)
        report["predictions"]["P2_work_basis"] = within_block_metrics(original, late_free)
        report["control"]["C1_insertion_zero"] = insertion_zero_control(late_backgrounds)
        p1 = (abs(report["predictions"]["P1_joint"]["original_mass"] - 1) < 1e-13
              and abs(report["predictions"]["P1_joint"]["late_free_mass"] - 1) < 1e-13
              and report["predictions"]["P1_joint"]["max_abs_error"] < 2e-14)
        p2 = report["predictions"]["P2_work_basis"]["basis_observable_changes"]
        c1 = report["control"]["C1_insertion_zero"]["changes_joint"]
        exp.check("P1", p1, f"joint max error={report['predictions']['P1_joint']['max_abs_error']:.3g}")
        exp.check("P2", p2, f"within-block joint-basis max error={report['predictions']['P2_work_basis']['max_joint_within_block_basis_probability_error']:.3g}")
        c1 = (c1 and abs(report["control"]["C1_insertion_zero"]["baseline_mass"] - 1) < 1e-13
              and abs(report["control"]["C1_insertion_zero"]["initial_mix_mass"] - 1) < 1e-13)
        exp.fail_check("C1", c1, f"insertion-0 TV={report['control']['C1_insertion_zero']['tv']:.8g}")
        path = report_path()
        if not exp.finish(report_path=path,rows=[json_safe(report["predictions"])],metadata=dict(
                controls=json_safe(report["control"]),reference="existing full-r complex128 matrix product")):
            raise AssertionError("Experiment harness failed")
        print(json.dumps(dict(ok=True,report=str(path))))
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
