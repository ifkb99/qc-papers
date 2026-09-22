"""Permutation-native Pauli propagation.

Standard PPS decomposes a Toffoli into Clifford+T and propagates through the
rotations. Those rotations include Hadamards, so the operator leaves the
diagonal mid-circuit even though the Toffoli as a whole is a permutation --
the peak Pauli count is therefore much larger than the final one, and is not
a Walsh quantity.

Treating X / CNOT / Toffoli as ATOMIC permutation gates fixes both problems.
Conjugating a Z-string by a permutation gives a diagonal operator, so the
expansion stays purely Z-type at every step, and the count after k gates is
exactly the Walsh sparsity of the k-gate suffix. Peak retained support is then

    N_max = max over suffixes of (Walsh sparsity of that suffix)

which is exact and computable, closing the gap Fable identified.

Conjugation rules for Z^z (all derived from (-1)^{popcount(perm(y) & z)}):

  X(q)            : Z^z -> (-1)^{z_q} Z^z                       1 term
  CNOT(c,t)       : Z^z -> Z^{z XOR (z_t << c)}                 1 term
  Toffoli(a,b,c)  : z_c = 0 -> Z^z                              1 term
                    z_c = 1 -> 1/2 [ Z^z + Z^{z^a} + Z^{z^b}
                                     - Z^{z^a^b} ]              4 terms
  (last line uses (-1)^{y_a y_b} = 1/2(1 + (-1)^{y_a} + (-1)^{y_b}
                                        - (-1)^{y_a + y_b}))
"""
from __future__ import annotations
from collections.abc import Iterable, Mapping
from operator import index


class _BinaryFrame:
    """Physical Walsh key z = F k; store F rows and F^-1 columns."""

    def __init__(self, n: int):
        self.rows = [1 << q for q in range(n)]
        self.inverse_columns = self.rows.copy()

    def cnot(self, control: int, target: int) -> None:
        self.rows[control] ^= self.rows[target]
        self.inverse_columns[target] ^= self.inverse_columns[control]

    def physical(self, key: int) -> int:
        return sum(((row & key).bit_count() & 1) << q
                   for q, row in enumerate(self.rows))

    def logical(self, physical: int) -> int:
        key = 0
        while physical:
            bit = physical & -physical
            key ^= self.inverse_columns[bit.bit_length() - 1]
            physical ^= bit
        return key

    def constraints(self, mask: int) -> tuple[int, ...]:
        return tuple(row for q, row in enumerate(self.rows) if (mask >> q) & 1)


class FramedTerms(Mapping[int, float]):
    """Read-only physical-key view of one dictionary in binary coordinates.

    Iteration decodes one key at a time; lookup applies the inverse frame.
    ``dict(view)`` explicitly allocates a second, physical-key dictionary.
    Values/items stream without building an intermediate key collection.
    The propagation owns the backing dictionary and stops mutating it before
    exposing this view. This interface defers key conversion, not coefficients.
    """

    def __init__(self, terms: dict[int, float], frame: _BinaryFrame):
        self._terms = terms
        self._frame = frame

    def __len__(self):
        return len(self._terms)

    def __iter__(self):
        for key in self._terms:
            yield self._frame.physical(key)

    def __getitem__(self, key):
        try:
            physical = index(key)
        except TypeError:
            raise KeyError(key) from None
        if physical < 0 or physical.bit_length() > len(self._frame.rows):
            raise KeyError(key)
        return self._terms[self._frame.logical(physical)]

    def values(self):
        return self._terms.values()

    def items(self):
        for key, value in self._terms.items():
            yield self._frame.physical(key), value


class PermPPSResult:
    def __init__(self):
        self.n_terms: list[int] = []
        self.n_max: int = 0
        self.final_terms: Mapping[int, float] = {}
        self.expectation: float = 0.0
        self.hit_cap: bool = False
        self.frame_cnot_updates: int = 0
        # (reverse steps completed, qubit mask, terms before, terms after).
        self.trace_events: list[tuple[int, int, int, int]] = []
        # Same layout, for `reset_before` clears.
        self.reset_events: list[tuple[int, int, int, int]] = []


def propagate_perm(circuit, zmask: int, delta: float = 0.0,
                   max_terms: int = 4_000_000,
                   max_weight: int | None = None, *,
                   trace_plus: Iterable[int] = (),
                   affine_frame: bool = False,
                   reset_before: Mapping[int, int] | None = None) -> PermPPSResult:
    """Back-propagate the Z-type observable Z^zmask through a classical
    reversible circuit, keeping only Z-strings (keys are z bitmasks).

    `delta`      : coefficient truncation (drop |c| < delta)
    `max_weight` : weight truncation (drop Pauli weight popcount(z) > k), the
                   other standard PPS knob. For permutation circuits the Pauli
                   weight of Z^z is popcount(z), which is exactly the Fourier
                   DEGREE of that Walsh coefficient -- so weight truncation here
                   is literally low-degree Fourier truncation of the pulled-back
                   Boolean function.

    `trace_plus`: qubits initially in an independent |+> product state.
                  Contract each immediately after its final use in reverse
                  propagation (its FIRST occurrence in forward order). All
                  remaining gates then act on other qubits, so the contraction
                  commutes with their propagation. Z-terms containing that
                  qubit vanish; the result is a REDUCED observable, not the full
                  Heisenberg expansion. Unspecified input qubits are |0> for
                  `expectation`, as before. Exact at delta=0, max_weight=None
                  up to the existing numerical floor; truncation need not
                  commute with tracing.

    `n_terms` and `n_max` count retained terms BEFORE each contraction so that
    the temporary pre-contraction expansion is included in the reported peak.
    They do not count simultaneous old/new dictionary allocations or bytes.

    `affine_frame`: opt in to binary Walsh-key coordinates. Each CNOT updates
                  two n-bit frame vectors instead of rebuilding the support
                  dictionary. Toffoli branches and |+> contractions use that
                  frame directly. The final_terms result is then a read-only
                  FramedTerms Mapping of PHYSICAL keys; iteration converts keys
                  on demand. Call dict(result.final_terms) only if a concrete
                  dictionary is required and budget that additional allocation.
                  The frame stores O(n^2) bits plus Python containers. It does
                  not reduce support cardinality or nonlinear-gate allocations.
                  max_weight still scans physical keys after each CNOT.

    `reset_before`: {k: mask} inserts a reset to |0> of the qubits in `mask`
                  immediately before forward op k (0 <= k <= len(logical); k = 0
                  is the input boundary, k = len(logical) the output). Its
                  adjoint on a diagonal observable is C105's clear rule:
                  Z^z -> Z^(z & ~mask), colliding keys summed, then the usual
                  coefficient filter. A unitary circuit whose qubits in `mask`
                  are |0> at that point on some input set computes, on that set,
                  the same expectation with or without the reset (a "virtual
                  reset", TODO 66); off that set it computes a different
                  function. Counts are recorded in n_terms before each clear.
                  Not supported with affine_frame; masks may not meet trace_plus.
    """
    if not circuit.is_classical():
        bad = next(op[0] for op in circuit.logical
                   if op[0] not in ("x", "cnot", "toffoli"))
        raise ValueError(f"circuit is not a permutation (found {bad!r})")

    frame = _BinaryFrame(circuit.n) if affine_frame else None
    if frame is not None:
        zmask = index(zmask)
        if zmask < 0 or zmask.bit_length() > circuit.n:
            raise ValueError("zmask contains a qubit outside the circuit")
        for op in circuit.logical:
            if len(set(op[1:])) != len(op) - 1 or any(
                    q < 0 or q >= circuit.n for q in op[1:]):
                raise ValueError("affine frame requires distinct in-range gate qubits")

    resets = dict(reset_before or {})
    if resets:
        if frame is not None:
            raise ValueError("reset_before is not supported with affine_frame")
        full = (1 << circuit.n) - 1
        for k, mask in resets.items():
            if not 0 <= k <= len(circuit.logical) or mask < 0 or mask & ~full:
                raise ValueError("reset_before step or mask out of range")

    def clear(terms: dict[int, float], mask: int, step: int) -> dict[int, float]:
        new: dict[int, float] = {}
        for z, c in terms.items():
            zz = z & ~mask
            new[zz] = new.get(zz, 0.0) + c
        floor_ok = ((lambda v: abs(v) >= delta) if delta > 0
                    else (lambda v: abs(v) > 1e-13))
        new = {z: v for z, v in new.items() if floor_ok(v)}
        res.reset_events.append((step, mask, len(terms), len(new)))
        return new

    terms: dict[int, float] = {zmask: 1.0}
    res = PermPPSResult()
    coefficients_filtered = False
    if resets.get(len(circuit.logical)):
        terms = clear(terms, resets[len(circuit.logical)], 0)

    trace_qubits = {index(q) for q in trace_plus}
    if any(q < 0 or q >= circuit.n for q in trace_qubits):
        raise ValueError("trace_plus contains a qubit outside the circuit")
    trace_mask_all = sum(1 << q for q in trace_qubits)
    if any(mask & trace_mask_all for mask in resets.values()):
        raise ValueError("reset_before masks may not contain trace_plus qubits")
    trace_at: dict[int, int] = {}
    unseen = set(trace_qubits)
    if unseen:
        for forward_step, op in enumerate(circuit.logical):
            for q in op[1:]:
                if q in unseen:
                    trace_at[forward_step] = trace_at.get(forward_step, 0) | (1 << q)
                    unseen.remove(q)
            if not unseen:
                break
    unused_mask = sum(1 << q for q in unseen)
    if unused_mask:
        before = len(terms)
        terms = {z: v for z, v in terms.items() if not (z & unused_mask)}
        res.trace_events.append((0, unused_mask, before, len(terms)))

    for forward_step in range(len(circuit.logical) - 1, -1, -1):
        op = circuit.logical[forward_step]
        deferred_cnot = frame is not None and op[0] == "cnot"
        new: dict[int, float] = terms if deferred_cnot else {}

        if op[0] == "x":
            q = op[1]
            predicate = frame.rows[q] if frame is not None else 1 << q
            for z, c in terms.items():
                new[z] = new.get(z, 0.0) + (
                    -c if (z & predicate).bit_count() & 1 else c)

        elif op[0] == "cnot":
            cq, tq = op[1], op[2]
            if frame is not None:
                frame.cnot(cq, tq)
                res.frame_cnot_updates += 1
            else:
                for z, c in terms.items():
                    zz = z ^ (((z >> tq) & 1) << cq)
                    new[zz] = new.get(zz, 0.0) + c

        else:                                   # toffoli(a, b, c)
            a, b, cq = op[1], op[2], op[3]
            ba = frame.inverse_columns[a] if frame is not None else 1 << a
            bb = frame.inverse_columns[b] if frame is not None else 1 << b
            predicate = frame.rows[cq] if frame is not None else 1 << cq
            for z, c in terms.items():
                if not ((z & predicate).bit_count() & 1):
                    new[z] = new.get(z, 0.0) + c
                    continue
                h = 0.5 * c
                for zz, s in ((z, h), (z ^ ba, h), (z ^ bb, h), (z ^ ba ^ bb, -h)):
                    new[zz] = new.get(zz, 0.0) + s

        if not deferred_cnot or not coefficients_filtered:
            if delta > 0:
                new = {z: v for z, v in new.items() if abs(v) >= delta}
            else:
                new = {z: v for z, v in new.items() if abs(v) > 1e-13}
            coefficients_filtered = True
        if max_weight is not None:
            new = {z: v for z, v in new.items()
                   if (frame.physical(z) if frame is not None else z).bit_count()
                   <= max_weight}
        terms = new

        res.n_terms.append(len(terms))
        if len(terms) > max_terms:
            res.hit_cap = True
            break
        mask = trace_at.get(forward_step, 0)
        if mask:
            before = len(terms)
            if frame is None:
                terms = {z: v for z, v in terms.items() if not (z & mask)}
            else:
                constraints = frame.constraints(mask)
                terms = {z: v for z, v in terms.items()
                         if not any((z & row).bit_count() & 1 for row in constraints)}
            res.trace_events.append((len(res.n_terms), mask, before, len(terms)))
        rmask = resets.get(forward_step, 0)
        if rmask:
            terms = clear(terms, rmask, len(res.n_terms))

    res.n_max = max(res.n_terms) if res.n_terms else 0
    res.final_terms = FramedTerms(terms, frame) if frame is not None else terms
    # Traced |+> qubits have been eliminated; <0|Z^z|0> = 1 on the rest.
    res.expectation = float(sum(terms.values()))
    return res
