"""Research-only clean arithmetic blocks and a two-buffer conditional instrument.

Specialized finite-work instrument, not a general circuit propagator. Exact
identities are evaluated in complex128; no amplitude truncation or known order.
Valid scope: clean scratch, x<N, commuting modular maps, terminating inverse QFT.
The bounded TODO66 pilot found no complete-call memory/runtime advantage over
the existing logical dense baseline when construction/certification was charged.
See notes/MC-memory-structure-pilots.md for scope, results and provenance.
"""
from __future__ import annotations
import math
import numpy as np
import walsh
from toffoli_arith import ToffoliModExp


def multipliers(N, a, width):
    if not isinstance(width, int) or not 1 <= width <= 24:
        raise ValueError("pilot requires integer width in 1..24")
    if not isinstance(N, int) or not 3 <= N <= 127 or math.gcd(a, N) != 1:
        raise ValueError("pilot requires 3<=N<=127 and gcd(a,N)=1")
    out, m = [], a % N
    for _ in range(width):
        out.append(m)
        m = m*m % N
    return out


def certify(N, a, width):
    """Exact selected-basis replay; shared existing interpreter, no new kernel.

    Retain only compact permutations, discarding each gate list before making
    the next. Both control branches and every valid x are checked. Dense full
    scratch tables are never created. Construction/replay must be charged.
    """
    ms = multipliers(N, a, width)
    me = ToffoliModExp(N, a, n_exp=1)
    if me.n_qubits > 63:
        raise ValueError("selected-input replay uses signed int64 labels")
    x = np.arange(N, dtype=np.int64)
    identity = x.copy()
    table, records = {}, []
    for m in dict.fromkeys(ms):
        qc = me.u_a(me.exp[0], m)
        expected = (m*x) % N
        for bit in (0, 1):
            labels = (x << me.x[0]) | (bit << me.exp[0])
            actual = walsh.classical_images(qc, labels)
            target = ((expected if bit else x) << me.x[0]) | (bit << me.exp[0])
            if not np.array_equal(actual, target):
                raise ArithmeticError("actual block violates clean valid-code map")
        table[m] = expected
        records.append(dict(multiplier=m, valid_basis_cases=2*N,
                            logical_gates=len(qc.logical)))
        del qc
    pairs = [(identity, table[m]) for m in ms]
    return pairs, records


def arithmetic_pairs(N, a, width):
    """Strong sparse logical baseline: no order, gate replay, table or cache."""
    def identity(ids):
        return ids.copy()
    return [(identity, lambda ids, m=m: (ids*m) % N)
            for m in multipliers(N, a, width)]


def two_buffer_path(pairs, initial, *, rng=None, output=None, return_state=False):
    """M_b=(I+(-1)^b exp(i*phase)P)/2, selected branch only.

    Caller supplies maps certified by certify(), or independently validated
    identity/valid modular permutations. No global assumption on dirty inputs.
    Pair order and feedback match lab.semiclassical.sample/sparse_path.
    Persistent ndarray amplitude payload: two vectors of 16*N bytes each.
    np.vdot returns scalar; ufuncs use out= or in-place updates throughout.
    Python/NumPy headers, maps, input and output ownership are additional.
    """
    if (rng is None) == (output is None):
        raise ValueError("supply RNG xor forced output")
    if output is not None and (not isinstance(output, (int, np.integer))
                               or not 0 <= output < 1 << len(pairs)):
        raise ValueError("forced output outside register")
    psi = np.array(initial, dtype=complex, copy=True)
    if psi.ndim != 1 or not psi.size or not np.all(np.isfinite(psi)):
        raise ValueError("initial must be a finite vector")
    if abs(float(np.vdot(psi, psi).real)-1) > 1e-12:
        raise ValueError("initial must have norm one")
    tmp = np.empty_like(psi)
    prefix, probability, log_probability = 0, 1., 0.
    max_norm_error = 0.
    for step, pair in enumerate(reversed(pairs)):
        # Source-to-destination scatter, not an inverse-map gather.
        tmp[pair[1]] = psi
        tmp *= np.exp(-1j*np.pi*prefix/(1 << step))
        tmp += psi
        w0 = float(np.vdot(tmp, tmp).real)/4
        # Now tmp is -(psi-phase*Ppsi), preserving destructive interference.
        tmp -= psi
        tmp -= psi
        w1 = float(np.vdot(tmp, tmp).real)/4
        total = w0+w1
        max_norm_error = max(max_norm_error, abs(total-1))
        if not math.isfinite(total) or abs(total-1) > 1e-10:
            raise ArithmeticError("instrument lost normalization")
        bit = int(rng.random() >= w0/total) if output is None else (int(output) >> step) & 1
        prefix |= bit << step
        weight = w1 if bit else w0
        p = weight/total
        probability *= p
        if p == 0:
            result = dict(output=int(output) if output is not None else prefix,
                          path_probability=0., log_probability=None,
                          max_normalization_error=max_norm_error,
                          zero_probability_step=step)
            if return_state:
                result["state"] = None
            return result
        log_probability += math.log(p)
        if bit:
            tmp *= -1
        else:
            tmp += psi
            tmp += psi
        tmp /= 2*math.sqrt(weight)
        psi, tmp = tmp, psi
    result = dict(output=prefix, path_probability=probability,
                  log_probability=log_probability,
                  max_normalization_error=max_norm_error)
    if return_state:
        result["state"] = psi
    return result
