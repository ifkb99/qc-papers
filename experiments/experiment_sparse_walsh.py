"""Exact prefix-Walsh contraction with mask-free bit blocks.

Opt-in implementation of accepted derivation S023deae7d4044fda. The original
C89, C94 and packed comparator remain unchanged. A separate versioned harness
owns predictions and finite checks. At most 2K+1 transfers for K active mask
positions; big-integer products, scans and workspace still cost resources.
No measured performance, outside-average or novelty claim is made here.
"""
from __future__ import annotations


MAX_BITS = 1_048_576


def _mask_free_block(state, width, shift, limit):
    """Joint carry/borrow counts; local words lie in [0,2^width)."""
    size = 1 << width
    row = [0, 0, 0, 0]
    for index, weight in enumerate(state):
        if not weight:
            continue
        carry_cut = size-shift-(index >> 1)
        borrow_cut = limit+(index & 1)
        row[0] += weight*max(0, carry_cut-borrow_cut)
        row[1] += weight*min(carry_cut, borrow_cut)
        row[2] += weight*(size-max(carry_cut, borrow_cut))
        row[3] += weight*max(0, borrow_cut-carry_cut)
    return row


def _signed_bit(state, shift, limit, alpha, beta):
    row = [0, 0, 0, 0]
    for index, weight in enumerate(state):
        if not weight:
            continue
        carry, borrow = index >> 1, index & 1
        for x in (0, 1):
            total = x+shift+carry
            after = 2*(total >> 1)+int(x-limit-borrow < 0)
            phase = (alpha*x) ^ (beta*(total & 1))
            row[after] += -weight if phase else weight
    return row


def prefix_walsh(bits, limit, shift, alpha, beta, *, stats=None):
    """Return sum(x<limit) chi_alpha(x) chi_beta(x+shift mod2^bits).

    Integer parameters only, with 0<=bits<=MAX_BITS, valid masks/limit, and
    arbitrary signed shift. State order is (carry,borrow)=00,01,10,11.
    No list of active positions or block matrices is retained. Normalizing a
    much longer supplied shift is charged separately from the transfer bound.
    """
    if any(type(value) is not int for value in (bits, limit, shift, alpha, beta)):
        raise TypeError("integer width, threshold, shift and masks required")
    if not 0 <= bits <= MAX_BITS:
        raise ValueError("width outside the experimental bound")
    size = 1 << bits
    if not 0 <= limit <= size or not 0 <= alpha < size or not 0 <= beta < size:
        raise ValueError("threshold/mask out of range")
    if stats is not None:
        stats.update(signed_bits=0, free_blocks=0)
    if not limit:
        return 0
    # For a power-of-two modulus this also normalizes negative shifts exactly.
    shift &= size-1
    full = limit == size
    low_limit = 0 if full else limit
    state = [1, 0, 0, 0]
    active, position = alpha | beta, 0
    signed_bits = free_blocks = 0
    while active:
        next_bit = (active & -active).bit_length()-1
        gap = next_bit-position
        if gap:
            mask = (1 << gap)-1
            state = _mask_free_block(state, gap, (shift >> position) & mask,
                                     (low_limit >> position) & mask)
            free_blocks += 1
        state = _signed_bit(state, (shift >> next_bit) & 1,
                            (low_limit >> next_bit) & 1,
                            (alpha >> next_bit) & 1, (beta >> next_bit) & 1)
        signed_bits += 1
        position = next_bit+1
        active &= active-1
    if position < bits:
        gap = bits-position
        mask = (1 << gap)-1
        state = _mask_free_block(state, gap, (shift >> position) & mask,
                                 (low_limit >> position) & mask)
        free_blocks += 1
    if stats is not None:
        stats.update(signed_bits=signed_bits, free_blocks=free_blocks)
    return sum(state) if full else state[1]+state[3]


def contract_pieces(width, pieces, queries):
    """Exact integer numerators for either packed or streamed intervals.

    The caller owns interval validity and schedule/representation guards.
    Prefix subtraction and every query use the same new kernel explicitly.
    """
    sums = [0]*len(queries)
    for lo, hi, shift in pieces:
        for index, (alpha, beta) in enumerate(queries):
            sums[index] += (prefix_walsh(width+1, hi, shift, alpha, beta)
                            - prefix_walsh(width+1, lo, shift, alpha, beta))
    return tuple(sums)
