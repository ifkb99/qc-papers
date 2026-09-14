"""Bounded precision diagnostic for the C59 coherent-route contraction.

This is a numerical gap audit, not an arbitrary-accuracy certificate.  The
fixed fixture is r=10=b*M, b=2, t=4, with Rx(pi/7) repeated blocks at
insertions 1 and 3 and two coherent reflections K_2(q=0,pi/5),
K_3(q=1,pi/5).  The same finite full-r product formula is evaluated with
complex64, complex128, long-double complex, and a stored-complex128 matrix
set merely widened to complex256.  The latter distinguishes input
coefficient rounding from later arithmetic precision.

P1: the complex128 full-r law agrees with the existing production
    CoherentReflectionCircuit marginal, while mantissa/phase dtypes are
    actually retained through the finite matrix products.
P2: complex64, complex128, and exact-rational-pi long-double inputs have
    finite, normalized laws whose TV gaps to the high-precision reference are
    recorded rather than treated as a certificate.
P3: widening stored complex128 coefficients is distinguishable from
    reconstructing the rational-pi inputs at long-double precision.
C1: stored-coefficient widening must fail to restore exact unitary structure;
    its matrix norm drift remains above the high-precision reconstructed input.

The independent calculation is the existing bounded full-r controlled-product
formula, not a second generic propagator.  Dense arrays are guarded before
allocation and remain far below 16 MiB.  Exact-zero and underflow claims are
outside this test; global TV is computed only for this 16-output fixture.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from fractions import Fraction
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from lab import Experiment
from lab.coherent_routes import CoherentReflectionCircuit


PERIOD, BLOCK, SECTORS, WIDTH = 10, 2, 5, 4
MAX_BYTES = 16 << 20
Q = 1 << WIDTH
BACKGROUND_INSERTIONS = (1, 3)
REFLECTION_INSERTIONS = (2, 3)


def guard(shape, dtype, label: str) -> None:
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} bytes exceeds 16 MiB")


def report_path(prefix: str = "precision_inputs") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def dtype_pair(name: str) -> tuple[np.dtype, np.dtype]:
    if name == "complex64":
        return np.dtype(np.float32), np.dtype(np.complex64)
    if name == "complex128":
        return np.dtype(np.float64), np.dtype(np.complex128)
    if name in ("clongdouble", "longdouble"):
        return np.dtype(np.longdouble), np.dtype(np.clongdouble)
    raise ValueError(f"unknown precision mode {name}")


def pi_value(real_dtype: np.dtype):
    real = real_dtype.type
    return real(np.arctan(real(1)) * real(4))


def rational_angle(num: int, den: int, real_dtype: np.dtype):
    return real_dtype.type(num) * pi_value(real_dtype) / real_dtype.type(den)


def unit_phase(numerator: int, denominator: int,
               real_dtype: np.dtype, complex_dtype: np.dtype):
    rem = int(numerator) % int(denominator)
    real = real_dtype.type
    angle = real(2) * pi_value(real_dtype) * real(rem) / real(denominator)
    return complex_dtype.type(np.cos(angle) + 1j * np.sin(angle))


def rx_matrix(theta, real_dtype: np.dtype, complex_dtype: np.dtype) -> np.ndarray:
    real = real_dtype.type
    c, s = real(np.cos(theta / real(2))), real(np.sin(theta / real(2)))
    return np.array([[c, -1j * s], [-1j * s, c]], dtype=complex_dtype)


def repeated_matrix(block: np.ndarray, real_dtype: np.dtype,
                    complex_dtype: np.dtype) -> np.ndarray:
    del real_dtype
    guard((PERIOD, PERIOD), complex_dtype, "repeated precision block")
    out = np.zeros((PERIOD, PERIOD), dtype=complex_dtype)
    for m in range(SECTORS):
        out[BLOCK*m:BLOCK*(m+1), BLOCK*m:BLOCK*(m+1)] = block
    return out


def shift_matrix(power: int, complex_dtype: np.dtype) -> np.ndarray:
    guard((PERIOD, PERIOD), complex_dtype, "precision shift")
    out = np.zeros((PERIOD, PERIOD), dtype=complex_dtype)
    for j in range(PERIOD):
        out[(j + int(power)) % PERIOD, j] = complex_dtype.type(1)
    return out


def route_matrix(q: int, real_dtype: np.dtype,
                 complex_dtype: np.dtype) -> np.ndarray:
    guard((PERIOD, PERIOD), complex_dtype, "precision route")
    out = np.zeros((PERIOD, PERIOD), dtype=complex_dtype)
    for m in range(SECTORS):
        for p in range(BLOCK):
            source = BLOCK*m + p
            target = BLOCK*((-m) % SECTORS) + p
            out[target, source] = unit_phase(-q*m, SECTORS,
                                             real_dtype, complex_dtype)
    return out


def fixture(mode: str) -> tuple[list[tuple[np.ndarray, np.ndarray]], np.ndarray,
                                 list[np.ndarray], list[np.ndarray], dict]:
    """Build the same finite product inputs at one explicitly selected dtype."""
    if mode == "stored_wide":
        source_real, source_complex = dtype_pair("complex128")
        target_real, target_complex = dtype_pair("longdouble")
        # Inputs are first constructed and rounded as complex128, then widened.
        cast_wide = True
    elif mode == "longdouble":
        source_real, source_complex = dtype_pair(mode)
        target_real, target_complex = source_real, source_complex
        cast_wide = False
    else:
        source_real, source_complex = dtype_pair(mode)
        target_real, target_complex = source_real, source_complex
        cast_wide = False

    def cast_matrix(matrix: np.ndarray) -> np.ndarray:
        return matrix.astype(target_complex, copy=False) if cast_wide else matrix

    real_dtype = source_real
    complex_dtype = source_complex
    background_block = rx_matrix(rational_angle(1, 7, real_dtype),
                                  real_dtype, complex_dtype)
    guard((PERIOD, PERIOD), target_complex, "precision identity")
    identity = np.eye(PERIOD, dtype=target_complex)
    repeated = cast_matrix(repeated_matrix(background_block, real_dtype,
                                           complex_dtype))
    blocks = [repeated if i + 1 in BACKGROUND_INSERTIONS else identity.copy()
              for i in range(WIDTH)]
    reflections = []
    for q in (0, 1):
        theta = rational_angle(1, 5, real_dtype)
        j = cast_matrix(route_matrix(q, real_dtype, complex_dtype))
        c = source_real.type(np.cos(theta / source_real.type(2)))
        s = source_real.type(np.sin(theta / source_real.type(2)))
        k = source_complex.type(c) * np.eye(PERIOD, dtype=source_complex)
        k = k - source_complex.type(1j * s) * j.astype(source_complex, copy=False)
        reflections.append(cast_matrix(k))

    pairs: list[tuple[np.ndarray, np.ndarray]] = []
    for i in range(WIDTH):
        gate = blocks[i]
        if i + 1 == REFLECTION_INSERTIONS[0]:
            gate = reflections[0] @ gate
        if i + 1 == REFLECTION_INSERTIONS[1]:
            gate = reflections[1] @ gate
        shift = shift_matrix(1 << i, target_complex)
        pairs.append((gate, gate @ shift))
    initial = np.zeros(PERIOD, dtype=target_complex)
    initial[0] = target_complex.type(1)
    return pairs, initial, blocks, reflections, dict(
        source_real=str(source_real), source_complex=str(source_complex),
        target_complex=str(target_complex), input_pi_precision=source_real.name,
        input_coefficients=("complex128 then widened" if cast_wide else "native"),
        _phase_real_dtype=source_real,
        qft_phase_policy=("target longdouble" if cast_wide else "native"),
    )


def full_distribution(mode: str) -> tuple[np.ndarray, dict, float, float]:
    pairs, initial, blocks, reflections, details = fixture(mode)
    phase_real_dtype = (np.dtype(np.longdouble) if mode == "stored_wide"
                        else details.pop("_phase_real_dtype"))
    details.pop("_phase_real_dtype", None)
    ctype = pairs[0][0].dtype
    if any(a.dtype != ctype or b.dtype != ctype for a,b in pairs):
        raise AssertionError("matrix-product dtype changed")
    guard((Q, PERIOD), ctype, f"{mode} amplitude workspace")
    amplitudes = np.zeros((Q, PERIOD), dtype=ctype)
    for exponent in range(Q):
        state = initial.copy()
        for i, (b0, b1) in enumerate(pairs):
            state = (b1 if exponent & (1 << i) else b0) @ state
        amplitudes[exponent] = state
    probs = np.zeros(Q, dtype=np.longdouble if ctype.itemsize > 16 else np.float64)
    for y in range(Q):
        state = np.zeros(PERIOD, dtype=ctype)
        for exponent in range(Q):
            phase = unit_phase(-y*exponent, Q,
                               phase_real_dtype, ctype)
            if np.asarray(phase).dtype != ctype:
                raise AssertionError("phase dtype changed")
            state += phase * amplitudes[exponent] / ctype.type(Q)
        probs[y] = np.vdot(state, state).real
    norm_errors = []
    eye = np.eye(PERIOD, dtype=ctype)
    for matrix in (*blocks, *reflections):
        norm_errors.append(float(np.max(np.abs(matrix.conj().T @ matrix - eye))))
    details.update(
        mantissa_bits=int(np.finfo(np.dtype(np.float32 if ctype.itemsize == 8
                                             else np.float64 if ctype.itemsize == 16
                                             else np.longdouble)).nmant),
        probability_sum=float(probs.sum()),
        max_input_unitary_norm_defect=max(norm_errors),
    )
    return np.asarray(probs, dtype=np.longdouble), details, max(norm_errors), float(probs.sum())


def assert_fixture_schedule() -> None:
    """Check the intended production-matched insertion schedule pre-measurement."""
    _, _, blocks, reflections, _ = fixture("complex128")
    identity = np.eye(PERIOD, dtype=np.complex128)
    expected_block = repeated_matrix(
        rx_matrix(rational_angle(1, 7, np.dtype(np.float64)),
                  np.dtype(np.float64), np.dtype(np.complex128)),
        np.dtype(np.float64), np.dtype(np.complex128))
    for i, block in enumerate(blocks):
        expected = expected_block if i + 1 in BACKGROUND_INSERTIONS else identity
        if not np.allclose(block, expected, rtol=0, atol=0):
            raise AssertionError(f"background schedule mismatch at insertion {i+1}")
    if len(reflections) != len(REFLECTION_INSERTIONS):
        raise AssertionError("reflection insertion schedule mismatch")


def production_marginal() -> np.ndarray:
    real, complex_dtype = dtype_pair("complex128")
    block = rx_matrix(rational_angle(1, 7, real), real, complex_dtype)
    theta = float(rational_angle(1, 5, real))
    circuit = CoherentReflectionCircuit(
        PERIOD, BLOCK, WIDTH, {1: block, 3: block},
        {REFLECTION_INSERTIONS[0]: (0, theta),
         REFLECTION_INSERTIONS[1]: (1, theta)})
    return np.array([
        sum(circuit.joint_probability(sector, y) for sector in range(SECTORS))
        for y in range(Q)
    ], dtype=np.float64)


def main() -> None:
    exp = Experiment("precision_inputs", doc=__doc__)
    exp.predict("P1", "complex128 full-r law agrees with production and dtype precision is retained")
    exp.predict("P2", "lower precision has finite normalized law with measurable TV gap to long-double input")
    exp.predict("P3", "widened stored coefficients differ from reconstructed rational-pi high precision")
    exp.must_fail("C1", "widening stored complex128 coefficients restores exact unitary structure")
    started = time.perf_counter()
    rows: list[dict] = []
    # Schedule assertion occurs before any dtype sweep measurement.
    assert_fixture_schedule()
    exp.log(f"schedule asserted: background insertions={BACKGROUND_INSERTIONS}, "
            f"reflection insertions={REFLECTION_INSERTIONS}")

    exp.section("P1/P2 dtype-parametric full-r contraction")
    results = {}
    for mode in ("complex64", "complex128", "longdouble", "stored_wide"):
        distribution, details, norm_defect, probability_sum = full_distribution(mode)
        results[mode] = distribution
        rows.append(dict(series="dtype", mode=mode, **details,
                         norm_defect=norm_defect,
                         probabilities=[float(x) for x in distribution]))

    raw_results = results
    native_masses = {mode: raw_results[mode].sum() for mode in raw_results}
    raw_masses = {mode: float(value) for mode,value in native_masses.items()}
    normalized = {mode: raw_results[mode] / native_masses[mode]
                  for mode in raw_results}
    if any(not np.isfinite(v).all() or np.min(v)<0
           or abs(v.sum()-np.longdouble(1)) > np.longdouble('1e-18')
           for v in normalized.values()):
        raise AssertionError("comparison laws did not normalize at longdouble precision")
    high = normalized["longdouble"]
    production = production_marginal()
    production_normalized = production / production.sum()
    c128_raw_error = float(np.max(np.abs(np.asarray(raw_results["complex128"], dtype=np.float64)
                                         - production)))
    c128_error = float(np.max(np.abs(np.asarray(normalized["complex128"], dtype=np.float64)
                                     - production_normalized)))
    raw_half_l1 = {mode: float(np.sum(np.abs(raw_results[mode] - raw_results["longdouble"])) / 2)
                   for mode in raw_results}
    tv = {mode: float(np.sum(np.abs(normalized[mode] - high)) / 2)
          for mode in results}
    max_pointwise = {mode: float(np.max(np.abs(normalized[mode] - high)))
                     for mode in results}
    rows.append(dict(series="law_comparison",
                     production_complex128_raw_max_abs_error=c128_raw_error,
                     production_complex128_normalized_max_abs_error=c128_error,
                     raw_masses=raw_masses, raw_half_l1_vs_longdouble=raw_half_l1,
                     normalized_total_variation=tv,
                     normalized_max_pointwise_error=max_pointwise,
                     longdouble_raw_sum=raw_masses["longdouble"],
                     production_raw_sum=float(production.sum()),
                     production_normalized_sum=float(production_normalized.sum())))
    exp.check("P1", c128_error < 3e-12
              and abs(float(normalized["complex128"].sum()) - 1) < 3e-12
              and all(np.isfinite(normalized[mode]).all() for mode in normalized),
              f"production normalized max abs error={c128_error:.2e}, "
              f"complex128 raw sum={raw_masses['complex128']:.16g}")
    exp.check("P2", abs(float(high.sum()) - 1) < 3e-18
              and abs(float(normalized["complex64"].sum()) - 1) < 3e-18
              and tv["complex64"] >= tv["complex128"]
              and tv["complex64"] > 0,
              f"TV complex64={tv['complex64']:.3e}, complex128={tv['complex128']:.3e}, "
              f"longdouble={tv['longdouble']:.3e}, complex64 raw sum={raw_masses['complex64']:.12g}")
    exp.check("P3", tv["stored_wide"] > 0
              and max_pointwise["stored_wide"] > 0,
              f"stored-wide TV={tv['stored_wide']:.3e}, "
              f"max pointwise={max_pointwise['stored_wide']:.3e}")

    high_defect = next(row["norm_defect"] for row in rows
                       if row.get("mode") == "longdouble")
    stored_defect = next(row["norm_defect"] for row in rows
                         if row.get("mode") == "stored_wide")
    widened_drift = stored_defect - high_defect
    _, _, stored_blocks, _, _ = fixture("complex128")
    c_exact = Fraction.from_float(float(stored_blocks[0][0,0].real))
    s_exact = Fraction.from_float(float(stored_blocks[0][0,1].imag))
    exact_rational_norm_residual = c_exact*c_exact+s_exact*s_exact-1
    exp.fail_check("C1", widened_drift > 1e-18 and stored_defect > high_defect
                   and exact_rational_norm_residual != 0,
                   f"stored-wide norm defect={stored_defect:.3e}, "
                   f"reconstructed-longdouble={high_defect:.3e}, "
                   f"gap={widened_drift:.3e}")
    rows.append(dict(series="unitary_gap_control", stored_norm_defect=stored_defect,
                     longdouble_norm_defect=high_defect,
                     stored_Rx_exact_rational_norm_residual=str(exact_rational_norm_residual),
                     widened_input_gap=widened_drift))

    report = report_path()
    exp.finish(report_path=report, rows=rows,
               metadata=dict(period=PERIOD, block_size=BLOCK, sectors=SECTORS,
                             width=WIDTH, max_dense_bytes=MAX_BYTES,
                             mantissa_float32=int(np.finfo(np.float32).nmant),
                             mantissa_float64=int(np.finfo(np.float64).nmant),
                             mantissa_longdouble=int(np.finfo(np.longdouble).nmant),
                             precision_scope="diagnostic only; no sampling certificate",
                             high_reference="long-double rational-pi finite product",
                             background_insertions=BACKGROUND_INSERTIONS,
                             reflection_insertions=REFLECTION_INSERTIONS,
                             stored_wide_qft_phase_policy="target longdouble; only work-gate coefficients are widened stored inputs",
                             elapsed_seconds=time.perf_counter() - started))
    print(f"report: {report}")


if __name__ == "__main__":
    try:
        main()
    except BaseException as exc:
        failure = report_path("precision_inputs_failure")
        failure.write_text(json.dumps(
            dict(ok=False, error=repr(exc), traceback=traceback.format_exc()),
            indent=2) + "\n")
        raise
