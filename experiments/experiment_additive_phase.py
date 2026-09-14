"""Physical additive residue phases for the frozen Gauss-sum/output probes.

TODO 35, derived before measuring: exp(2*pi*i*k*x/N) is a product of
binary one-qubit phases with an explicitly retained global phase. It does
not require a multiplicative discrete logarithm. The main experiment lives
in experiment_gauss_kernels and experiment_additive_phase_output, with
predictions, controls, raw reports and resource counters.

This builder reuses the capped clean-orbit mixer/multiplier fixtures and
Circuit/statevec. The additive phase itself takes n+1 Pauli rotations; the
rest of this tiny arithmetic compiler still enumerates truth tables.
No second propagator, scalable arithmetic compiler or sampler is added.
"""
from __future__ import annotations

import math
import numpy as np

from circuits import Circuit
from experiments.experiment_clean_orbit_gates import (
    Fixture, append, checked_run, controlled_multiplier, mixer)
import statevec

MAX_BYTES = 16 * 1024 * 1024


def additive_gate(N, k, qubits):
    """Exactly G_k on the listed binary register, including its global phase.

    Extends to every binary label, including x>=N, by the SAME displayed
    phase. Floats are not an exact synthesis/precision certificate.
    """
    if (type(N) is not int or type(k) is not int or not 3 <= N <= 17
            or not -2 <= k <= 2 or not isinstance(qubits, (tuple, list))
            or len(qubits) != (N-1).bit_length()
            or any(type(q) is not int or not 0 <= q < 10 for q in qubits)
            or len(set(qubits)) != len(qubits)):
        raise ValueError("outside frozen additive phase fixture")
    qc = Circuit(max(qubits)+1)
    angle_sum = 0.
    for bit, q in enumerate(qubits):
        angle = 2*np.pi*k*(1 << bit)/N
        qc.rz(q, angle)
        angle_sum += angle
    qc.rot((0, 0), -angle_sum)
    return qc


def physical_circuit(k, *, location="interior"):
    if type(k) is not int or k not in (0, 1, 2) or location not in ("interior", "end"):
        raise ValueError("frozen output uses k=0,1,2 and interior/end")
    f, width = Fixture(13, 2, 3, 12), 3
    # Retain the previous ten-qubit layout; flag q=6 is now unused and clean.
    offset = f.nb+1
    qc = Circuit(offset+width).x(0)
    exponent = list(range(offset, qc.n))
    for control in exponent:
        qc.h(control)
    append(qc, mixer(f, np.pi/4))
    for i, control in enumerate(exponent):
        controlled_multiplier(qc, f, 1 << i, control)
        append(qc, mixer(f, (-1)**i*np.pi/(5+i)))
        if i == 1 and location == "interior":
            append(qc, additive_gate(f.N, k, f.work))
    if location == "end":
        append(qc, additive_gate(f.N, k, f.work))
    qc.qft(exponent, inverse=True)
    return qc, f, offset


def physical_output(k, *, location="interior", budget=None):
    qc, f, offset = physical_circuit(k, location=location)
    if 16*24*(1 << qc.n) > MAX_BYTES:
        raise MemoryError("additive full-output pre-allocation cap")
    if budget is None:
        budget = [0]
    psi, costs = checked_run(qc, statevec.basis(qc.n, 0), budget)
    probabilities = np.sum(abs(psi.reshape(8, 1 << offset))**2, axis=1)
    idx = np.arange(len(psi))
    dirty = ((idx >> f.nw) & ((1 << (f.np+1))-1)) != 0
    leakage = float(np.sum(abs(psi[dirty])**2))
    norm = float(np.vdot(psi, psi).real)
    if not np.all(np.isfinite(probabilities)) or not math.isfinite(norm+leakage):
        raise ArithmeticError("nonfinite physical output")
    return dict(probabilities=probabilities.tolist(), leakage=leakage, norm=norm,
                qubits=qc.n, cumulative_gate_entry_updates=budget[0],
                phase_rotation_count=f.nw+1, **costs)
