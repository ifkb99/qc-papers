"""Initial contraction-regime probe for b=3 forward work channels.

This is an exploratory diagnostic for TODO26, not a sampler or a numerical
certificate.  A fixed r=9,b=3 exact-input circuit repeats W=Rx(pi/7) after
each of 16 controls, with no reflections.  Existing VerifiedFiniteWork._build
supplies the branch matrices; the channel representation below is only the
traceless Hermitian coordinate representation of those existing density
updates.  Channel-product length is the sole varying parameter.

PREDICTIONS, WRITTEN BEFORE MEASURING.

  P1  The repeated-Rx schedule can have unit one-step traceless norm while a
      multi-step product contracts for every initial coarse sector; report the
      SVD as an explicitly uncertified diagnostic only.
  P2  The no-background b=3 null retains a nonzero conserved traceless mode
      through every product length, so a maximally mixed fixed point alone is
      not evidence of mixing.
  C1  The must-fail shortcut "unit one-step norm means no product contraction"
      is contradicted if the measured repeated schedule contracts.

No new propagator is introduced: all branches come from _build and all
channel applications use the existing density-channel formula.  Width 16,
work dimension 3, all three sectors, and 8-by-8 traceless representations
are preflighted explicitly.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from fractions import Fraction

import numpy as np

from lab import Experiment
from lab.verified_finite_work import VerifiedFiniteWork
from lab.verified_prefix import VerifiedReflectionCircuit
from lab.semiclassical import _density_forward_step


PERIOD = 9
BLOCK = 3
WIDTH = 16
MAX_BYTES = 16 << 20
DIM = BLOCK * BLOCK - 1

exp = Experiment(__name__.rsplit(".", 1)[-1], doc=__doc__, exit_on_fail=False)
exp.predict("P1", "repeated exact b=3 channels can contract in products despite unit one-step diagnostics")
exp.predict("P2", "the no-background conserved traceless mode survives every product length")
exp.must_fail("C1", "unit one-step contraction coefficient one rules out all multi-step product contraction")


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def report_path():
    return Path("out") / f"forward_mixing_probe_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    return value


def acb_matrix_to_complex(matrix):
    return np.array([[complex(matrix[i, j]) for j in range(matrix.ncols())]
                     for i in range(matrix.nrows())], dtype=complex)


def traceless_hermitian_basis():
    """Orthonormal Hilbert-Schmidt Gell-Mann basis for Herm_0(C^3)."""
    basis = [np.diag([1, -1, 0]).astype(complex) / np.sqrt(2),
             np.diag([1, 1, -2]).astype(complex) / np.sqrt(6)]
    for i in range(BLOCK):
        for j in range(i + 1, BLOCK):
            h = np.zeros((BLOCK, BLOCK), complex)
            h[i, j] = h[j, i] = 1 / np.sqrt(2)
            basis.append(h)
            h = np.zeros((BLOCK, BLOCK), complex)
            h[i, j] = -1j / np.sqrt(2)
            h[j, i] = 1j / np.sqrt(2)
            basis.append(h)
    if len(basis) != DIM:
        raise AssertionError("traceless basis dimension mismatch")
    gram = np.array([[np.trace(a.conj().T @ b) for b in basis] for a in basis])
    if np.max(np.abs(gram - np.eye(DIM))) > 1e-12:
        raise AssertionError("traceless basis is not Hilbert-Schmidt orthonormal")
    return tuple(basis)


def channel_matrix(branches, basis):
    """Real Herm_0 coordinate matrix of the existing density channel."""
    b0, b1 = branches
    out = np.empty((DIM, DIM), dtype=float)
    max_coefficient_imag = 0.0
    max_hermitian_residual = 0.0
    for col, h in enumerate(basis):
        image = _density_forward_step(b0,b1,h,multiply=lambda A,B:A@B,
                                      adjoint=lambda A:A.conj().T)
        max_hermitian_residual = max(max_hermitian_residual,
                                     float(np.max(np.abs(image - image.conj().T))))
        max_trace = abs(np.trace(image))
        if max_trace > 1e-10:
            raise AssertionError("density channel did not preserve traceless subspace")
        for row, g in enumerate(basis):
            coefficient = np.trace(g.conj().T @ image)
            max_coefficient_imag = max(max_coefficient_imag, float(abs(coefficient.imag)))
            out[row, col] = float(coefficient.real)
    return out, max_coefficient_imag, max_hermitian_residual


def channel_products(circuit, initial_sector):
    worker = VerifiedFiniteWork(circuit, 0)
    from flint import ctx
    with ctx.workprec(192):
        branches, _states, _final = worker._build(initial_sector)
    basis = traceless_hermitian_basis()
    guard((DIM, DIM), dtype=np.float64, label="traceless channel")
    product = np.eye(DIM)
    rows = []
    one_step_max = []
    trace_imag_max = 0.0
    for length, pair in enumerate(branches, 1):
        pair = tuple(acb_matrix_to_complex(x) for x in pair)
        if any(np.max(np.abs(B.conj().T@B-np.eye(BLOCK))) > 1e-12 for B in pair):
            raise AssertionError("branch unitarity diagnostic failed")
        matrix, imag, hermitian = channel_matrix(pair, basis)
        trace_imag_max = max(trace_imag_max, imag)
        one_step = float(np.linalg.svd(matrix, compute_uv=False)[0])
        product = matrix @ product
        product_norm = float(np.linalg.svd(product, compute_uv=False)[0])
        one_step_max.append(one_step)
        rows.append(dict(length=length, one_step_svd_norm=one_step,
                         product_svd_norm=product_norm,
                         product_frobenius_norm=float(np.linalg.norm(product, "fro")),
                         representation_dimension=DIM,
                         hermitian_residual=hermitian,
                         svd_is_uncertified_diagnostic=True))
    return rows, trace_imag_max


def conserved_mode_null(circuit, initial_sector):
    """Track an explicit Fourier-eigenprojector difference in the no-background null."""
    worker = VerifiedFiniteWork(circuit, 0)
    from flint import ctx
    with ctx.workprec(192):
        branches, _states, _final = worker._build(initial_sector)
    first_b0 = acb_matrix_to_complex(branches[0][0])
    first_b1 = acb_matrix_to_complex(branches[0][1])
    shift = first_b0.conj().T @ first_b1
    eigenvalues, eigenvectors = np.linalg.eig(shift)
    if np.max(np.abs(np.abs(eigenvalues) - 1)) > 1e-10:
        raise AssertionError("null shift is not unitary")
    v0 = eigenvectors[:, 0]
    v1 = eigenvectors[:, 1]
    mode = v0[:, None] @ v0.conj()[None, :] - v1[:, None] @ v1.conj()[None, :]
    mode = mode / np.linalg.norm(mode, "fro")
    residuals = []
    current = mode
    for b0_raw, b1_raw in branches:
        b0, b1 = acb_matrix_to_complex(b0_raw), acb_matrix_to_complex(b1_raw)
        current = _density_forward_step(b0,b1,current,multiply=lambda A,B:A@B,
                                        adjoint=lambda A:A.conj().T)
        residuals.append(float(np.linalg.norm(current - mode, "fro")))
    return dict(max_residual=max(residuals), residuals=residuals,
                mode_frobenius_norm=float(np.linalg.norm(mode, "fro")),
                mode_is_traceless=abs(np.trace(mode)) < 1e-12)


def build_circuit(backgrounds):
    return VerifiedReflectionCircuit(PERIOD, WIDTH, backgrounds, {}, block_size=BLOCK)


def main():
    started = time.time()
    report = {"status": "PASS", "repeated": {}, "null": {}, "controls": {}}
    try:
        # Aggregate allowances: all _build matrices, basis and temporaries,
        # both live coordinate matrices, and every retained report row.
        if ((3*WIDTH+40)*BLOCK**2*512 + 6*WIDTH*1024 + 4*DIM**2*16) > MAX_BYTES:
            raise MemoryError("aggregate retained allocation allowance exceeded")
        repeated_backgrounds = {s: ("x", Fraction(1, 7)) for s in range(1, WIDTH + 1)}
        repeated = build_circuit(repeated_backgrounds)
        null = build_circuit({})
        for sector in range(repeated.sectors):
            rows, imag = channel_products(repeated, sector)
            report["repeated"][str(sector)] = dict(rows=rows,
                                                    max_channel_imaginary_coordinate=imag)
        for sector in range(null.sectors):
            rows, imag = channel_products(null, sector)
            report["null"][str(sector)] = dict(rows=rows,
                                               max_channel_imaginary_coordinate=imag,
                                               conserved_mode=conserved_mode_null(null, sector))
        repeated_rows = [r for d in report["repeated"].values() for r in d["rows"]]
        all_one_step = max(abs(r["one_step_svd_norm"] - 1) for r in repeated_rows) < 2e-12
        contracted = all(report["repeated"][str(sector)]["rows"][-1]["product_svd_norm"] < 0.9
                         for sector in range(repeated.sectors))
        null_mode = all(d["conserved_mode"]["max_residual"] < 2e-12
                        and d["conserved_mode"]["mode_is_traceless"]
                        and abs(d["rows"][-1]["product_svd_norm"] - 1) < 2e-12
                        for d in report["null"].values())
        report["controls"] = dict(one_step_all_unit=all_one_step,
                                   repeated_product_contracts=contracted,
                                   single_step_shortcut_fails=all_one_step and contracted,
                                   repeated_final_product_norms={s: report["repeated"][str(s)]["rows"][-1]["product_svd_norm"]
                                                                for s in range(repeated.sectors)},
                                   null_final_product_norms={s: report["null"][str(s)]["rows"][-1]["product_svd_norm"]
                                                             for s in range(null.sectors)})
        p1 = (all_one_step and contracted
              and all(d["max_channel_imaginary_coordinate"] < 1e-12
                      and max(r["hermitian_residual"] for r in d["rows"]) < 1e-12
                      for d in report["repeated"].values()))
        p2 = (null_mode and all(d["max_channel_imaginary_coordinate"] < 1e-12
                                and max(r["hermitian_residual"] for r in d["rows"]) < 1e-12
                                for d in report["null"].values()))
        exp.check("P1", p1, "repeated-Rx b=3 product contracts despite unit one-step SVD diagnostics")
        exp.check("P2", p2, "no-background conserved traceless mode retained")
        exp.fail_check("C1", report["controls"]["single_step_shortcut_fails"],
                       f"final product norms={report['controls']['repeated_final_product_norms']}")
    except Exception as exc:
        report["status"] = "FAIL"
        report["error"] = {"type": type(exc).__name__, "message": str(exc),
                            "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.time() - started
    path = report_path()
    if report["status"] == "FAIL":
        exp.check("P1",False,str(report.get("error")))
    ok = exp.finish(report_path=path,metadata=json_safe(report))
    report["status"] = "PASS" if ok else "FAIL"
    print(json.dumps({"status": report["status"], "report": str(path),
                      "elapsed_seconds": report["elapsed_seconds"]}, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
