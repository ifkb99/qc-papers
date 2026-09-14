"""Complete-output comparison for the TODO39 earlier-phase schedule.

Frozen fixture: N=61, a=2, r=60, b=3, t=8, split s=7, with W0=Rx/phase
work_block(pi/4), W1=work_block(-pi/10), and G1(j)=exp(2*pi*i*2^j/61).
An identical G1 is inserted after v low controls, v=0,...,7, before the
remaining low controls and the unchanged W1/G1/128-shift suffix.

The target is assembled from literal indexed vector shifts, pointwise phases,
and repeated 3-by-3 blocks.  Its output law is evaluated by a zero-padded
256-point FFT.  A second baseline dephases only the single high input bit but
retains the early phase and the exact z=y mod 2 feedback.  This is a finite
complete-law comparison, not a sampler or a general propagator.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  ER's retained literal target columns give finite, normalized
      zero-padded FFT laws, with each column preserving norm.
  P2  The dephased law agrees with its independent 128-point feedback-aware
      conditional construction and normalizes for every v.
  P3  The first output bit is uniform for every v under the [-2,4] support
      and dist_60(128,0)=8>6 certificate; full-law TV outcomes remain open.
  C1  Omitting the nonzero feedback phase in a synthetic Fourier-basis
      boundary gives a measurably wrong law.
"""
from __future__ import annotations

import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from experiments.experiment_clean_orbit_output import work_block
from experiments.experiment_earlier_phase_cycles import direct_column as er_direct_column
from lab import Experiment
from lab.fourier_sampling import geometric_sum


N, BASE, PERIOD, BLOCK = 61, 2, 60, 3
WIDTH, SPLIT, Q, L, H = 8, 7, 1 << 8, 1 << 7, 2
MAX_BYTES = 16 << 20
# The first preflight (preserved as a failure artifact) totaled 1,122,304
# units after replacing duplicated omitted columns by one build.  Keep a
# modest explicit 2M bound for the frozen eight-position run.
MAX_SCALAR_TERMS = 2_000_000
TOL = 3e-10
V_VALUES = tuple(range(8))


def report_path():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"earlier_phase_feedback_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard_bytes(payload, label):
    payload = int(payload)
    if payload < 0 or payload > MAX_BYTES:
        raise MemoryError(f"{label} payload {payload} exceeds 16 MiB")
    return payload


def phase(index, counters, bounds):
    if counters["phase_queries"] + 1 > bounds["phase_queries"]:
        raise MemoryError("phase-query bound before modular power")
    counters["phase_queries"] += 1
    counters["modular_power_queries"] += 1
    return complex(np.exp(2j * np.pi * pow(BASE, int(index) % PERIOD, N) / N))


def preflight():
    # Conservative simultaneous numeric payload.  This includes the old and
    # new target matrices, the retained omitted matrix, complex FFT input/
    # output, abs/square temporaries, feedback temporaries, converted report
    # laws, and work blocks. Python object headers/RSS are not claimed.
    matrix = Q * PERIOD * np.dtype(np.complex128).itemsize
    chi = L * PERIOD * np.dtype(np.complex128).itemsize
    payload_components = {
        "target_old_matrix": matrix,
        "target_new_matrix": matrix,
        "omitted_matrix": matrix,
        "dephased_chi_matrix": chi,
        "fft_complex_input_output": 2 * matrix,
        "fft_abs_square_temporaries": 2 * Q * PERIOD * np.dtype(np.float64).itemsize,
        "feedback_complex_temporaries": 2 * chi,
        "feedback_abs_square_temporaries": 2 * L * PERIOD * np.dtype(np.float64).itemsize,
        "retained_converted_law_lists": 5 * len(V_VALUES) * Q * np.dtype(np.float64).itemsize,
        "work_blocks": 2 * BLOCK * BLOCK * np.dtype(np.complex128).itemsize,
        "synthetic_control_laws": 2 * L * np.dtype(np.float64).itemsize,
        "work_vectors_and_live_law_temporaries": 32 * PERIOD * 16 + 32 * Q * 8,
    }
    payload = sum(payload_components.values())
    guard_bytes(payload, "aggregate earlier-phase output payload")
    operation_bounds = {
        "direct_column_calls": len(V_VALUES) * Q,
        "omitted_column_calls": Q,
        "block_product_entries": (len(V_VALUES) + 1) * Q
                                * (PERIOD // BLOCK) * BLOCK * BLOCK,
        # Target: 3 initial + at most 6 late; omitted: at most 6; chi: 3.
        "phase_queries": (len(V_VALUES) * Q * 3 * BLOCK
                           + Q * 2 * BLOCK
                           + len(V_VALUES) * L * BLOCK),
        "shift_operations": len(V_VALUES) * Q * 3 + Q * 2,
        "chi_entries": len(V_VALUES) * L * BLOCK,
        "fft_law_calls": 44,  # 8 target + 1 omitted + 8*(1+2+1) + 3 control
        "fft_input_entries": (9 * Q * PERIOD + 8 * Q * PERIOD
                               + 8 * 3 * L * PERIOD + 3 * L),
        "synthetic_geometric_terms": L,
    }
    scalar = sum(operation_bounds.values())
    if scalar > MAX_SCALAR_TERMS:
        raise MemoryError("earlier-phase output scalar preflight exceeds cap")
    return {
        "payload_bytes": payload,
        "payload_components": payload_components,
        "operation_bounds": operation_bounds,
        "scalar_work_bound_total": scalar,
        "max_scalar_terms": MAX_SCALAR_TERMS,
    }


def ensure_capacity(counters, bounds, key, increment):
    if counters.get(key, 0) + int(increment) > bounds[key]:
        raise MemoryError(f"{key} bound before operation")


def omitted_column(exponent, W0, W1, counters, bounds):
    """The C78 schedule after omitting the added early phase.

    This baseline is independent of the nominal insertion v, so it is built
    once rather than redundantly eight times.
    """
    high, low = divmod(int(exponent), L)
    ensure_capacity(counters, bounds, "block_product_entries",
                    (PERIOD // BLOCK) * BLOCK * BLOCK)
    ensure_capacity(counters, bounds, "phase_queries", 2 * BLOCK)
    ensure_capacity(counters, bounds, "shift_operations", 2)
    vector = np.zeros(PERIOD, dtype=complex)
    vector[:BLOCK] = W0[:, 0]
    counters["shift_operations"] += 1
    vector = np.roll(vector, low)
    for cell in range(PERIOD // BLOCK):
        counters["block_product_entries"] += BLOCK * BLOCK
        start = cell * BLOCK
        vector[start:start+BLOCK] = W1 @ vector[start:start+BLOCK]
    for index in np.flatnonzero(vector):
        vector[index] *= phase(index, counters, bounds)
    counters["shift_operations"] += 1
    vector = np.roll(vector, L * high)
    return vector


def build_columns(insertion, W0, W1, counters, bounds, er_counters, *, omit_early=False):
    columns = np.zeros((Q, PERIOD), dtype=complex)
    for exponent in range(Q):
        if omit_early:
            ensure_capacity(counters, bounds, "omitted_column_calls", 1)
            counters["omitted_column_calls"] += 1
            columns[exponent] = omitted_column(exponent, W0, W1, counters, bounds)
        else:
            # Use ER's retained literal direct-column oracle for the target.
            ensure_capacity(counters, bounds, "direct_column_calls", 1)
            ensure_capacity(counters, bounds, "phase_queries", 3 * BLOCK)
            ensure_capacity(counters, bounds, "block_product_entries",
                            (PERIOD // BLOCK) * BLOCK * BLOCK)
            ensure_capacity(counters, bounds, "shift_operations", 3)
            counters["direct_column_calls"] += 1
            before = dict(er_counters)
            columns[exponent] = er_direct_column(
                exponent, insertion, W0, W1, er_counters)
            deltas = {key: er_counters[key] - before[key] for key in er_counters}
            if deltas["direct_phase_queries"] > 3 * BLOCK:
                raise AssertionError("ER phase count exceeded per-column bound")
            if deltas["direct_block_products"] > (PERIOD // BLOCK) * BLOCK * BLOCK:
                raise AssertionError("ER block count exceeded per-column bound")
            if deltas["direct_shift_ops"] > 3:
                raise AssertionError("ER shift count exceeded per-column bound")
            counters["phase_queries"] += deltas["direct_phase_queries"]
            counters["modular_power_queries"] += deltas["modular_pow_queries"]
            counters["block_product_entries"] += deltas["direct_block_products"]
            counters["shift_operations"] += deltas["direct_shift_ops"]
    return columns


def fft_law(columns, counters, bounds):
    ensure_capacity(counters, bounds, "fft_law_calls", 1)
    ensure_capacity(counters, bounds, "fft_input_entries",
                    columns.shape[0] * columns.shape[1])
    counters["fft_law_calls"] += 1
    counters["fft_input_entries"] += int(columns.shape[0] * columns.shape[1])
    amplitudes = np.fft.fft(columns, axis=0) / Q
    return np.sum(np.abs(amplitudes)**2, axis=1)


def dephased_chi(insertion, W0, counters, bounds):
    chi = np.zeros((L, PERIOD), dtype=complex)
    for low in range(L):
        a = low % (1 << insertion)
        ensure_capacity(counters, bounds, "chi_entries", BLOCK)
        ensure_capacity(counters, bounds, "phase_queries", BLOCK)
        for u in range(BLOCK):
            index = (u + low) % PERIOD
            phase_index = (u + a) % PERIOD
            counters["chi_entries"] += 1
            chi[low, index] += W0[u, 0] * phase(
                phase_index, counters, bounds)
    return chi


def dephased_law(chi, counters, bounds):
    ensure_capacity(counters, bounds, "fft_law_calls", 1)
    ensure_capacity(counters, bounds, "fft_input_entries", Q * chi.shape[1])
    counters["fft_law_calls"] += 1
    counters["fft_input_entries"] += int(Q * chi.shape[1])
    return np.sum(np.abs(np.fft.fft(chi, n=Q, axis=0))**2, axis=1) / (Q*L)


def feedback_law(chi, counters, bounds):
    result = np.zeros(Q, dtype=float)
    for z in range(H):
        coefficient = np.exp(-2j*np.pi*z*np.arange(L)/Q)[:, None]
        ensure_capacity(counters, bounds, "fft_law_calls", 1)
        ensure_capacity(counters, bounds, "fft_input_entries", L * chi.shape[1])
        counters["fft_law_calls"] += 1
        counters["fft_input_entries"] += int(L * chi.shape[1])
        low_amplitudes = np.fft.fft(chi*coefficient, axis=0) / L
        conditional = np.sum(np.abs(low_amplitudes)**2, axis=1)
        for w in range(L):
            result[z+H*w] = conditional[w] / H
    return result


def no_feedback_law(chi, counters, bounds):
    ensure_capacity(counters, bounds, "fft_law_calls", 1)
    ensure_capacity(counters, bounds, "fft_input_entries", L * chi.shape[1])
    counters["fft_law_calls"] += 1
    counters["fft_input_entries"] += int(L * chi.shape[1])
    low_amplitudes = np.fft.fft(chi, axis=0) / L
    conditional = np.sum(np.abs(low_amplitudes)**2, axis=1)
    result = np.zeros(Q, dtype=float)
    for z in range(H):
        result[z+H*np.arange(L)] = conditional / H
    return result


def synthetic_feedback_control(counters, bounds):
    # A low-register Fourier basis state with frequency k=0.  Nonzero z
    # shifts it away from the integer Fourier bin; omission leaves a delta.
    coeff = np.ones((L, 1), dtype=complex)
    z = 1
    # Exercise the SAME candidate implementations, not a second standalone
    # FFT expression that would leave feedback_law itself uncontrolled.
    proper = H * feedback_law(coeff, counters, bounds)[z::H]
    missing = H * no_feedback_law(coeff, counters, bounds)[z::H]
    ensure_capacity(counters, bounds, "synthetic_geometric_terms", L)
    counters["synthetic_geometric_terms"] += L
    expected = np.array([abs(geometric_sum(L, -z-H*w, Q))**2 / L**2
                         for w in range(L)])
    if (not np.all(np.isfinite(proper)) or not np.all(np.isfinite(missing))
            or abs(float(proper.sum())-1.) >= TOL
            or abs(float(missing.sum())-1.) >= TOL
            or np.max(np.abs(proper-expected)) >= TOL
            or abs(float(missing[0])-1.) >= TOL
            or np.max(np.abs(missing[1:])) >= TOL):
        raise AssertionError("synthetic feedback law or geometric control invalid")
    return float(np.sum(np.abs(proper-missing))/2), proper, missing


def law_diagnostics(values):
    values = np.asarray(values, dtype=float)
    even_mass = float(values[::2].sum()) if values.shape == (Q,) else float("nan")
    odd_mass = float(values[1::2].sum()) if values.shape == (Q,) else float("nan")
    return {
        "shape": list(values.shape),
        "finite": bool(np.all(np.isfinite(values))),
        "nonnegative": bool(values.size and np.min(values) >= -TOL),
        "normalization_error": abs(float(values.sum()) - 1.)
        if values.size else float("inf"),
        "even_mass": even_mass,
        "odd_mass": odd_mass,
        "even_mass_error": abs(even_mass - .5),
        "odd_mass_error": abs(odd_mass - .5),
        "valid": bool(values.shape == (Q,)
                      and np.all(np.isfinite(values))
                      and np.min(values) >= -TOL
                      and abs(float(values.sum()) - 1.) < TOL
                      and abs(even_mass - .5) < TOL
                      and abs(odd_mass - .5) < TOL),
    }


def main():
    exp = Experiment("earlier_phase_feedback", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "ER literal columns produce finite normalized FFT laws")
    exp.predict("P2", "dephased law equals feedback-aware 128-point conditionals")
    exp.predict("P3", "first output bit is uniform for every earlier insertion")
    exp.must_fail("C1", "synthetic omission of nonzero Fourier feedback is wrong")
    started = time.perf_counter()
    report = {"status": "FAIL", "rows": []}
    p1 = p2 = p3 = c1 = False
    try:
        report["preflight"] = preflight()
        bounds = report["preflight"]["operation_bounds"]
        counters = {"direct_column_calls": 0, "omitted_column_calls": 0,
                    "block_product_entries": 0,
                    "phase_queries": 0, "modular_power_queries": 0,
                    "shift_operations": 0, "chi_entries": 0,
                    "fft_law_calls": 0,
                    "fft_input_entries": 0,
                    "synthetic_geometric_terms": 0,
                    "dephased_calls": 0, "feedback_calls": 0,
                    "omitted_calls": 0}
        er_counters = {"direct_phase_queries": 0,
                       "modular_pow_queries": 0,
                       "direct_block_products": 0,
                       "direct_shift_ops": 0}
        W0 = work_block(np.pi/4)
        W1 = work_block(-np.pi/10)
        law_rows = []
        deph_rows = []
        omit_rows = []
        target_column_norm_errors = {}
        target_first_bit_errors = {}
        target_vs_omit = {}
        target_vs_deph = {}
        deph_feedback_errors = {}
        deph_no_feedback_tvs = {}
        feedback_rows = []
        no_feedback_rows = []
        diagnostics = {"target": [], "omitted": [], "dephased": [],
                       "feedback": [], "no_feedback": []}
        # Omitting the added phase gives the same schedule for all nominal v.
        omit = build_columns(0, W0, W1, counters, bounds, er_counters,
                             omit_early=True)
        counters["omitted_calls"] += 1
        omit_law = fft_law(omit, counters, bounds)
        diagnostics["omitted"] = [law_diagnostics(omit_law)
                                   for _ in V_VALUES]
        for insertion in V_VALUES:
            target = build_columns(insertion, W0, W1, counters, bounds,
                                   er_counters)
            target_law = fft_law(target, counters, bounds)
            chi = dephased_chi(insertion, W0, counters, bounds)
            counters["dephased_calls"] += 1
            deph_law = dephased_law(chi, counters, bounds)
            counters["feedback_calls"] += 1
            feedback = feedback_law(chi, counters, bounds)
            no_feedback = no_feedback_law(chi, counters, bounds)
            law_rows.append(target_law.tolist())
            omit_rows.append(omit_law.tolist())
            deph_rows.append(deph_law.tolist())
            feedback_rows.append(feedback.tolist())
            no_feedback_rows.append(no_feedback.tolist())
            diagnostics["target"].append(law_diagnostics(target_law))
            diagnostics["dephased"].append(law_diagnostics(deph_law))
            diagnostics["feedback"].append(law_diagnostics(feedback))
            diagnostics["no_feedback"].append(law_diagnostics(no_feedback))
            target_first_bit_errors[str(insertion)] = abs(float(target_law[::2].sum())-.5)
            target_vs_omit[str(insertion)] = float(np.sum(abs(target_law-omit_law))/2)
            target_vs_deph[str(insertion)] = float(np.sum(abs(target_law-deph_law))/2)
            deph_feedback_errors[str(insertion)] = float(np.max(abs(deph_law-feedback)))
            deph_no_feedback_tvs[str(insertion)] = float(np.sum(abs(deph_law-no_feedback))/2)
            target_column_norm_errors[str(insertion)] = float(
                np.max(np.abs(np.sum(np.abs(target)**2, axis=1)-1.)))
        synthetic_tv, synthetic_proper, synthetic_missing = synthetic_feedback_control(
            counters, bounds)
        counters_ok = all(counters[key] <= value for key, value in bounds.items())
        er_reconciled = {
            "direct_phase_queries": er_counters["direct_phase_queries"],
            "modular_pow_queries": er_counters["modular_pow_queries"],
            "direct_block_products": er_counters["direct_block_products"],
            "direct_shift_ops": er_counters["direct_shift_ops"],
        }
        report.update({
            "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                        "t": WIDTH, "s": SPLIT, "Q": Q, "L": L, "H": H},
            "target_laws": law_rows,
            "omit_early_laws": omit_rows,
            "dephased_laws": deph_rows,
            "feedback_laws": feedback_rows,
            "no_feedback_laws": no_feedback_rows,
            "law_diagnostics": diagnostics,
            "target_first_bit_errors": target_first_bit_errors,
            "target_vs_omit_tv": target_vs_omit,
            "target_vs_dephased_tv": target_vs_deph,
            "dephased_feedback_max_errors": deph_feedback_errors,
            "dephased_missing_feedback_tv": deph_no_feedback_tvs,
            "target_column_norm_errors": target_column_norm_errors,
            "synthetic_feedback_tv": synthetic_tv,
            "synthetic_feedback_proper": synthetic_proper.tolist(),
            "synthetic_feedback_missing": synthetic_missing.tolist(),
            "counters": counters,
            "er_counters": er_reconciled,
            "counter_bounds_ok": counters_ok,
            "actual_named_work_total": sum(counters[key] for key in bounds),
            "python_version": sys.version,
            "numpy_version": np.__version__,
            "platform": platform.platform(),
        })
        p1 = (counters_ok
              and all(item["valid"] for item in diagnostics["target"])
              and diagnostics["omitted"][0]["valid"]
              and max(target_column_norm_errors.values()) < TOL
              and counters["direct_column_calls"] == len(V_VALUES)*Q
              and counters["omitted_column_calls"] == Q
              and all(counters[key] == bounds[key]
                      for key in ("block_product_entries", "shift_operations",
                                  "chi_entries", "fft_law_calls", "fft_input_entries",
                                  "synthetic_geometric_terms"))
              and counters["phase_queries"] == counters["modular_power_queries"])
        p2 = (all(item["valid"] for item in diagnostics["dephased"])
              and all(item["valid"] for item in diagnostics["feedback"])
              and all(item["valid"] for item in diagnostics["no_feedback"])
              and max(deph_feedback_errors.values()) < TOL)
        p3 = (max(target_first_bit_errors.values()) < TOL
              and diagnostics["omitted"][0]["valid"])
        c1 = synthetic_tv > 1e-3
        report["checks"] = {"P1": p1, "P2": p2, "P3": p3, "C1": c1}
        report["full_law_tv_classification"] = {
            str(v): {"omit": ("below1e-3" if target_vs_omit[str(v)] < 1e-3
                               else "between1e-3_and1e-2" if target_vs_omit[str(v)] < 1e-2
                               else "above1e-2"),
                     "dephased": ("below1e-3" if target_vs_deph[str(v)] < 1e-3
                                  else "between1e-3_and1e-2" if target_vs_deph[str(v)] < 1e-2
                                  else "above1e-2")}
            for v in V_VALUES}
        report["threshold_booleans"] = {
            str(v): {
                "omit_gt_1e-3": target_vs_omit[str(v)] > 1e-3,
                "omit_gt_1e-2": target_vs_omit[str(v)] > 1e-2,
                "dephased_gt_1e-3": target_vs_deph[str(v)] > 1e-3,
                "dephased_gt_1e-2": target_vs_deph[str(v)] > 1e-2,
            } for v in V_VALUES}
        report["status"] = "PASS" if p1 and p2 and p3 and c1 else "FAIL"
        exp.check("P1", p1, "ER direct FFT laws are normalized")
        exp.check("P2", p2, "dephased law matches feedback-aware conditional construction")
        exp.check("P3", p3, "first output bit uniform for every insertion")
        exp.fail_check("C1", c1, "missing-feedback synthetic law differs")
    except Exception as exc:
        report["status"] = "FAIL"
        report["exception"] = repr(exc)
        report["traceback"] = __import__("traceback").format_exc()
        exp.log("EXCEPTION", repr(exc))
        for name in ("P1", "P2", "P3"):
            exp.check(name, False, "exception before completion")
        exp.fail_check("C1", False, "exception before completion")
    report["elapsed_seconds"] = time.perf_counter() - started
    path = report_path()
    ok = exp.finish(report_path=path, rows=[report], metadata={
        "fixture": {"N": N, "a": BASE, "r": PERIOD, "b": BLOCK,
                    "t": WIDTH, "s": SPLIT},
        "max_numeric_payload_bytes": MAX_BYTES,
        "max_scalar_work": MAX_SCALAR_TERMS,
        "reference": "ER direct-column oracle plus zero-padded FFT; independent dephased feedback formulas",
        "no_new_propagator": True, "no_sampler": True,
    })
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
