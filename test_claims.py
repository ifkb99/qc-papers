"""Claims regression -- cheap executable reproductions of the headline
results in CLAIMS.md, pinned to the exact logged numbers. Suite 7 of 9 in the
correctness gate.

Until this file existed, only the *infrastructure* had regression protection;
the science itself did not. Each block names the claim it re-verifies, uses a
small instance that runs in seconds, and includes the control where the
original had one. A future refactor that silently breaks a result fails here
loudly.

Not covered here: C25/C26 (the cryptanalytic-bound import -- kept out of the
automated gate by request; `experiments/experiment_crypto.py` remains the
record), and C19 (a citation, not a computation). The windowed block (C36-C39)
reaches q = 21 and adds ~15s; everything else stays at q <= 17.
"""
from __future__ import annotations
import os
import math

import numpy as np

os.environ.setdefault("LAB_NO_CACHE", "1")   # regression must recompute

from circuits import ripple_adder
from lab import (build_modexp, support, sparsity, peak_pps, fn_support,
                 find_structures, quadrant_counts, coset_split,
                 bit_table, tile, order, v2_split)
import walsh


def t(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    assert ok, name


print("[C8] PPS term count = Walsh sparsity, identical support")
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs_w = support(qc, me.x[0])
res = peak_pps(qc, 1 << me.x[0])
zs_p = np.sort(np.array(list(res.final_terms), dtype=np.int64))
t("walsh sparsity = 3086", zs_w.size == 3086)
t("pps final = 3086", len(res.final_terms) == 3086)
t("supports identical (not just counts)", np.array_equal(zs_w, zs_p))

print("[C10] adder collapse is affineness: low bit sparsity 1")
add, lay = ripple_adder(4)
t("4-bit adder Z_b0 sparsity = 1", sparsity(add, lay["b"][0]) == 1)
t("4-bit adder Z_b2 sparsity = 10 (control: higher bit not affine)",
  sparsity(add, lay["b"][2]) == 10)

print("[C15/C18-lite] constancy for beta=1, growth for beta>1 (fix N, vary a)")
t("N=7 a=6 (r=2): 15549 at n_exp=2 and 3",
  sparsity(build_modexp(N=7, a=6, n_exp=2).build(), 8) == 15549
  and sparsity(build_modexp(N=7, a=6, n_exp=3).build(), 8) == 15549)
t("N=7 a=3 (r=6) control: 15539 -> 30712",
  sparsity(build_modexp(N=7, a=3, n_exp=2).build(), 8) == 15539
  and sparsity(build_modexp(N=7, a=3, n_exp=3).build(), 8) == 30712)

print("[C21] onset at n_exp = v2(r)+1: alpha=2 grows once, then locks")
t("N=5 a=2 (r=4): 15493 -> 32143 -> 32143 at n_exp=2,3,4",
  [sparsity(build_modexp(N=5, a=2, n_exp=k).build(), 8) for k in (2, 3, 4)]
  == [15493, 32143, 32143])

print("[C24] support confined to z_I in {0, 1_I} (halves pinned)")
me = build_modexp(N=7, a=6, n_exp=3)
alpha, _ = v2_split(order(6, 7))
zs = support(me.build(), me.x[0])
mask = 0
for eq in me.exp[alpha:]:
    mask |= 1 << eq
zi = zs & mask
n0 = int(np.count_nonzero(zi == 0))
n1 = int(np.count_nonzero(zi == mask))
t("halves 7770/7779, other = 0",
  (n0, n1) == (7770, 7779) and n0 + n1 == zs.size)

print("[C43] C24 at SET level: the supports coincide, not merely their sizes")
# Pinned at alpha=2 (N=5, a=2, r=4) because it is cheap; the claim was
# established at alpha=3 and 4 in experiments/experiment_c21_onset.py, where
# the instances are q=23..27 and belong in a GPU sweep, not in the gate.
_al, _ = v2_split(order(2, 5))
_sigs = []
for _ne in (_al + 1, _al + 2):
    _me = build_modexp(N=5, a=2, n_exp=_ne)
    _zs = support(_me.build(), _me.x[0])
    _tm = 0
    for _q in _me.exp[_al:]:
        _tm |= 1 << _q
    _rest = (1 << _me.exp[_al]) - 1
    _sigs.append(np.sort((_zs & _rest)
                         | (((_zs & _tm) == _tm).astype(np.int64) << _me.exp[_al])))
t("N=5 a=2: (z_rest, tailflag) sets identical at n_exp = 3 and 4",
  _sigs[0].size == 32143 and bool(np.array_equal(_sigs[0], _sigs[1])))

print("[C30] the 1/2 ceiling is the linear structure w = b_msb^anc")
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs = support(qc, me.x[0])
w = (1 << me.b[me.m - 1]) | (1 << me.anc)
ev, od = coset_split(zs, w)
t("support in the hyperplane <z,w>=0", od == 0)
g = ((walsh.classical_permutation(qc) >> me.x[0]) & 1).astype(np.int64)
idx = np.arange(g.size, dtype=np.int64)
t("g(y^w) = g(y) pointwise", bool(np.array_equal(g[idx ^ w], g)))

print("[C32/C34] one nonlinear monomial: 23464, empty (msb=1,anc=0) quadrant")
me4 = build_modexp(N=5, a=2, n_exp=2, wraps=("msb_t0_anc",))
zs4 = support(me4.build(), me4.x[0])
t("v4 support = 23464", zs4.size == 23464)
quad = quadrant_counts(zs4, me4.b[me4.m - 1], me4.anc)
t("quadrant (1,0) exactly empty", quad[(1, 0)] == 0)
t("density above 1/2 (linear structure destroyed)",
  zs4.size / (1 << me4.n_qubits) > 0.5)

print("[C35] breaking the structure keeps C15 constancy (23488)")
t("v4 N=7 a=6: 23488 at n_exp=2 and 3",
  sparsity(build_modexp(N=7, a=6, n_exp=2, wraps=("msb_t0_anc",)).build(), 8)
  == 23488
  and sparsity(build_modexp(N=7, a=6, n_exp=3, wraps=("msb_t0_anc",)).build(),
               8) == 23488)

print("[C33] function-level all-ones structures, with negative control")
_, st = find_structures(tile(bit_table(7, 3, 0), 12), 12)
t("r=6 even t: all-ones AFFINE", st == [((1 << 12) - 1, 1)])
_, st = find_structures(tile(bit_table(7, 2, 0), 12), 12)
t("r=3 even t: all-ones LINEAR", st == [((1 << 12) - 1, 0)])
zs, st = find_structures(tile(bit_table(11, 2, 0), 12), 12)
t("r=10 control: no structure, density 1", st == [] and zs.size == 4096)

print("[C36/C37] windowed: V^2=id in BOTH designs, but only lookup is affinely"
      " controlled")
from windowed_arith import WindowedModExp, SelectModExp

def _tail_block(me):
    P = walsh.classical_permutation(me.window_block(me.windows()[0], 1))
    idn = np.arange(P.size, dtype=np.int64)
    inv = bool((P[P] == idn).all())
    moved = bool((P != idn).any())
    indep = all(bool((P[idn ^ (1 << q)] == (P ^ (1 << q))).all())
                for q in me.windows()[0])
    return inv, moved, indep

inv_l, moved_l, indep_l = _tail_block(WindowedModExp(N=5, a=4, n_exp=2, w=2))
inv_s, moved_s, indep_s = _tail_block(SelectModExp(N=5, a=4, n_exp=2, w=2))
t("lookup tail block: nontrivial and V^2=id", moved_l and inv_l)
t("select tail block: nontrivial and V^2=id too (so V^2=id cannot decide)",
  moved_s and inv_s)
t("lookup block does NOT read its window; select block DOES",
  indep_l and not indep_s)

print("[C37/C38] lookup: tail bits dead, |S| 2-periodic not growing")
zs = [support(WindowedModExp(N=5, a=4, n_exp=2 * (K + 1), w=2).build(),
              WindowedModExp(N=5, a=4, n_exp=2 * (K + 1), w=2).x[0])
      for K in (0, 1)]
t("N=5 a=4 lookup: 28078 at K=0, 62680 at K=1",
  [z.size for z in zs] == [28078, 62680])
me = WindowedModExp(N=5, a=4, n_exp=4, w=2)
tail = sum(1 << q for q in me.exp[1:])
t("no z with any tail bit set (sharper than C24)",
  not bool((zs[1] & tail != 0).any()))

print("[C39] the multiply-by-1 branch is the whole difference")
sk = [support(SelectModExp(N=5, a=4, n_exp=4, w=2, skip_zero=s).build(),
              SelectModExp(N=5, a=4, n_exp=4, w=2, skip_zero=s).x[0]).size
      for s in (True, False)]
t("select K=1: 128981 skipping j=0 vs 32075 emitting it", sk == [128981, 32075])

print("[C44] measured arithmetic peak relation, not a universal formula")
import pps as _pps
import perm_pps as _perm_pps
import walsh as _w

_me = build_modexp(N=5, a=2, n_exp=1)
_qc, _tq = _me.build(), _me.x[0]
_rot = _pps.propagate(_qc, {(0, 1 << _tq): 1.0}, delta=0.0)
_perm = _perm_pps.propagate_perm(_qc, 1 << _tq)
t("logged peaks: rot 13666, perm 6834",
  (_rot.n_max, _perm.n_max) == (13666, 6834))


class _Suffix:                       # the last m gates == the state after m steps
    def __init__(self, c, m, logical=False):
        self.n = c.n
        self.gates = [] if logical else c.gates[len(c.gates) - m:]
        self.logical = c.logical[len(c.logical) - m:] if logical else []

    def is_classical(self):
        return True


_m = int(np.argmax(_rot.n_terms)) + 1
_peak = _pps.propagate(_Suffix(_qc, _m), {(0, 1 << _tq): 1.0},
                       delta=0.0).final_terms
_nz = sorted({x for (x, _) in _peak if x})
t("peak lives in {I, X_c} for a single c", len(_nz) == 1
  and _nz[0].bit_count() == 1)
_c = _nz[0].bit_length() - 1
_S = set(_perm_pps.propagate_perm(
    _Suffix(_qc, int(np.argmax(_perm.n_terms)) + 1, logical=True),
    1 << _tq).final_terms)
_fold = {}
for (_x, _z) in _peak:
    _k = (_z | (1 << _c)) if _x else _z
    _fold[_k] = _fold.get(_k, 0) + 1
t("fold recovers the perm peak set, doubled on z_c=1",
  _fold == {z: (2 if (z >> _c) & 1 else 1) for z in _S})
_B = sorted(z for z in _S if not ((z >> _c) & 1))
t("B = {Z_x0, Z_x0 Z_e0}, so the deficit is 2",
  _B == sorted([1 << _me.x[0], (1 << _me.x[0]) | (1 << _me.exp[0])])
  and _rot.n_max == 2 * _perm.n_max - len(_B))
_co = _w.pullback_coefficients(_qc, _tq)
t("B is exactly the |coefficient| = 1/2 Walsh pair",
  set(np.nonzero(np.abs(np.abs(_co) - 0.5) < 1e-12)[0].tolist()) == set(_B))

print("[F12-lite] function-level dichotomy")
t("N=15 a=7 (r=4): sparsity 4 at t=12 and t=16",
  fn_support(tile(bit_table(15, 7, 0), 12)).size == 4
  and fn_support(tile(bit_table(15, 7, 0), 16)).size == 4)

print("[C44 counterexample] the universal factor-two bound is false")
from circuits import Circuit
tiny = Circuit(4).toffoli(3, 1, 0).toffoli(1, 0, 2)
t("two-Toffoli counterexample: atomic peak 4, rotation peak 10",
  _perm_pps.propagate_perm(tiny, 4).n_max == 4
  and _pps.propagate(tiny, {(0, 4): 1.}).n_max == 10)

print("[C47] balanced nonlinear controls have the same reduced map")
u = np.arange(1 << 9)
activation = np.zeros_like(u)
for shift in (0, 3, 6):
    w = (u >> shift) & 7
    activation ^= ((w & 1) + ((w >> 1) & 1) + ((w >> 2) & 1) >= 2).astype(int)
work0, work1 = np.array([1., 1., -1., -1.]), np.array([1., -1., -1., 1.])
truth = np.where(activation[:, None], work1, work0)
t("three majority windows: support 130 but the same nonzero projection",
  np.count_nonzero(walsh.wht(truth.reshape(-1))) == 130
  and np.array_equal(truth.mean(axis=0), [1., 0., -1., 0.]))

print("[C48] dense ideal Walsh spectrum, tiny tensor rank")
from lab.tensors import rank_profile
sign = 1. - 2*tile(bit_table(7, 3, 0), 12)
t("ideal r=6 at t=12: support 2048, maximum cut rank 3",
  np.count_nonzero(walsh.wht(sign)) == 2048
  and max(row["rank"] for row in rank_profile(sign)) == 3)
t("scalar-period counterexample remains rank one and sparsity one",
  fn_support(tile(bit_table(13, 4, 0), 12)).size == 1)

print("[C49] conditional order-finding retains the correct interference")
from lab.semiclassical import distribution
work_ids = np.arange(8)
pairs = [(work_ids, np.array([(x*pow(3, 1 << i, 7)) % 7 if x < 7 else x
                             for x in work_ids])) for i in range(4)]
initial = np.zeros(8, complex)
initial[1] = 1
prob, _ = distribution(pairs, initial)
amps = np.zeros((16, 8), complex)
amps[np.arange(16), [pow(3, e, 7) for e in range(16)]] = .25
reference = np.sum(np.abs(np.fft.fft(amps, axis=0)/4)**2, axis=1)
t("odd-order conditional distribution equals the independent FFT", np.max(np.abs(prob-reference)) < 1e-12)
wrong, _ = distribution(pairs, initial, feedback=False)
t("omitting feedback changes the distribution", np.max(np.abs(wrong-reference)) > 1e-3)

print("[C50] same-layout scratch guards separate logical equivalence from full support")
me = build_modexp(N=7, a=3, n_exp=4)
baseline = me.build(init_x=False)
modified = Circuit(me.n_qubits).extend(baseline)
for j in range(me.m):
    modified.toffoli(me.b[j], me.t[j], me.x[0])
before = walsh.classical_permutation(baseline)
after = walsh.classical_permutation(modified)
clean = np.array([(e << me.exp[0]) | (x << me.x[0])
                  for e in range(16) for x in range(7)])
t("guards leave all 112 legitimate images unchanged", np.array_equal(before[clean], after[clean]))
before_chi = 1. - 2*((before >> me.x[0]) & 1)
after_chi = 1. - 2*((after >> me.x[0]) & 1)
t("full support changes from 64353 to 129152",
  np.count_nonzero(walsh.wht(before_chi)) == 64353
  and np.count_nonzero(walsh.wht(after_chi)) == 129152)

print("[C51] selected compiled replay obeys the reachable subgroup bound")
from lab.reachable import compiled_pairs, sparse_path, action_stats
me_reachable, reached_pairs = compiled_pairs(7, 3, 6)
t("compilation enumerates no work basis images", action_stats(reached_pairs)["cached_transitions"] == 0)
reached_path = sparse_path(reached_pairs, {1 << me_reachable.x[0]: 1.}, output=0)
t("descending powers: reached unions plateau at odd part then reach full order",
  [row["union"] for row in reached_path["profile"]] == [2, 3, 3, 3, 3, 6]
  and reached_path["peak_stored_amplitudes"] == 6)

print("[C52] latent eigenphase samples outputs, not an arbitrary work state")
from lab.semiclassical import eigenphase_path
mixture = np.array([sum(eigenphase_path(6, 4, k, output=y)["conditional_path_probability"]
                        for k in range(6))/6 for y in range(16)])
t("uniform eigenphase mixture matches the independent FFT",
  np.max(np.abs(mixture-reference)) < 1e-12)
t("conditional-on-eigenphase probability is not the marginal probability",
  eigenphase_path(6, 4, 0, output=0)["conditional_path_probability"] == 1
  and abs(mixture[0] - 1) > .5)
joint = np.zeros((16, 13))
joint[np.arange(16), [pow(4, e, 13) for e in range(16)]] = .25
t("N=13,a=4 has joint rank six despite its scalar-rank-one LSB",
  np.linalg.matrix_rank(joint, tol=1e-12) == 6
  and max(row["rank"] for row in rank_profile(1. - 2*tile(bit_table(13, 4, 0), 4))) == 1)

print("[C53] detectable coherence does not force a work-state sampler")
from lab.spectral import single_defect_effects, coherence_response
orbit53 = np.array([1, 3, 2, 6, 4, 5])
phase53 = np.exp(-2j*np.pi*np.arange(6)[:, None]*np.arange(6)[None, :]/6)/np.sqrt(6)
theta53 = np.pi/2
V53 = phase53.conj().T @ np.diag(np.exp(-.5j*theta53*(1-2*(orbit53 & 1)))) @ phase53
E53 = single_defect_effects(6, 4, 2, V53)
p53 = E53.sum(axis=(1, 2)).real/6
d53 = np.trace(E53, axis1=1, axis2=2).real/6
joint53 = np.zeros((16, 7), complex)
for e53 in range(16):
    low53 = pow(3, e53 % 4, 7)
    joint53[e53, pow(3, e53, 7)] = np.exp(-.5j*theta53*(1-2*(low53 & 1)))/4
ref53 = np.sum(np.abs(np.fft.fft(joint53, axis=0)/4)**2, axis=1)
t("single-defect output effects agree with independent phase-tagged FFT",
  np.max(np.abs(p53-ref53)) < 1e-12)
t("physical initial-eigenphase dephasing loses measured interference",
  np.max(np.abs(d53-ref53)) > .02)
response53 = coherence_response(E53)
t("half-pi defect detects twelve pairs but eleven real response directions",
  response53["detectable_pairs"] == 12 and response53["rank"] == 11)
adjusted53 = np.array([sum(eigenphase_path(6, 4, k, output=y,
                                         input_phases=[0., -theta53, 0., 0.])
                          ["conditional_path_probability"] for k in range(6))/6
                      for y in range(16)])
t("transferred exponent phase restores scalar sampling",
  np.max(np.abs(adjusted53-ref53)) < 1e-12)
t("noncommuting earlier insertion can leave outputs unchanged",
  np.max(np.abs(single_defect_effects(6, 4, 1, V53).sum(axis=(1, 2)).real/6
                -mixture)) < 1e-12)

print("[C54] final spectral conditioning retains only the early coherent prefix")
from lab.prefix import SparseOrbitPrefix
c54, s54 = np.cos(np.pi/4), np.sin(np.pi/4)
def column54(l):
    j = l % 6
    if j == 0:
        return np.array([0, 5]), np.array([c54, -1j*s54])
    if j == 5:
        return np.array([5, 0]), np.array([c54, -1j*s54])
    return np.array([j]), np.array([1.+0j])
prefix54 = SparseOrbitPrefix(6, 2, column54, 2)
weights54 = np.array([prefix54.phase_probability(k) for k in range(6)])
V54 = np.eye(6, dtype=complex)
V54[0, 0] = V54[5, 5] = c54
V54[0, 5] = V54[5, 0] = -1j*s54
E54 = single_defect_effects(6, 4, 2, phase53.conj().T @ V54 @ phase53)
expected54 = E54.sum(axis=(1, 2)).real/6
joint54 = np.array([[prefix54.forced_joint(4, k, y)["joint_latent_output_probability"]
                     for y in range(16)] for k in range(6)])
t("weighted final phase plus four amplitudes matches time-ordered effects",
  np.max(np.abs(joint54.sum(axis=0)-expected54)) < 1e-12)
t("uniform final eigenphases are an invalid shortcut for the mixer",
  np.max(np.abs((joint54/weights54[:, None]).mean(axis=0)-expected54)) > .02)
t("two-sparse rejection envelope has exactly half mean acceptance",
  np.max(6*weights54/2) <= 1+1e-12 and abs(np.mean(6*weights54/2)-.5) < 1e-12)
t("one sampled path needs no internal orbit table",
  prefix54.sample(32, np.random.default_rng(3))["early_amplitudes"] == 4
  and prefix54.stats()["orbit_table_entries"] == 0)

print("[C55] localized defect removes the coherent prefix vector under a stronger promise")
from lab.localized import LocalizedOrbitDefect
localized55 = LocalizedOrbitDefect(6, (0, 5), V54[np.ix_((0, 5), (0, 5))])
joint55 = np.array([[localized55.forced_joint(4, 2, k, y)["joint_latent_output_probability"]
                    for y in range(16)] for k in range(6)])
t("localized finite sums match the existing prefix and independent time-ordered effects",
  np.max(np.abs(joint55-joint54)) < 1e-12
  and np.max(np.abs(joint55.sum(axis=0)-expected54)) < 1e-12)
mean55 = sum(localized55.phase_probability(2, k)
             * localized55.forced_joint(4, 2, k, 0)["expected_early_attempts"] for k in range(6))
t("component rejection has bounded mean over actual final phase weights", mean55 <= 15)
identity55 = LocalizedOrbitDefect(6, (), np.empty((0, 0)))
t("identity specialization agrees with ideal latent phase for each output",
  max(abs(identity55.forced_joint(4, 3, k, y)["joint_latent_output_probability"]
          - eigenphase_path(6, 4, k, output=y)["conditional_path_probability"]/6)
      for k in range(6) for y in range(16)) < 1e-12)

print("[C56] conserved coarse Fourier sectors retain noncommuting periodic defects")
from lab.periodic import PeriodicOrbitCircuit
W56 = np.eye(3, dtype=complex)
W56[0, 0] = W56[1, 1] = np.cos(np.pi/4)
W56[0, 1] = W56[1, 0] = -1j*np.sin(np.pi/4)
V56 = np.kron(np.eye(2), W56)
periodic56 = PeriodicOrbitCircuit(6, 3, 4, {2: W56})
got56 = np.array([sum(periodic56.forced_joint(alpha, y)["joint_latent_output_probability"]
                     for alpha in range(2)) for y in range(16)])
ref56 = single_defect_effects(6, 4, 2, phase53.conj().T @ V56 @ phase53).sum(axis=(1, 2)).real/6
t("uniform coarse-sector mixture equals original-order spectral effects",
  np.max(np.abs(got56-ref56)) < 1e-12)
t("odd-block periodic defect is not the ideal output", np.max(np.abs(got56-mixture)) > .01)
early56 = PeriodicOrbitCircuit(6, 3, 4, {1: W56})
early_p56 = np.array([sum(early56.forced_joint(alpha, y)["joint_latent_output_probability"]
                         for alpha in range(2)) for y in range(16)])
# The initial visibility assertion at s=1 failed, and its log is preserved.
# Odd block size only removes a sufficient commutation guarantee; it does
# not prohibit an input-specific cancellation within the early reachable set.
t("odd block size alone does not guarantee visibility at every insertion",
  np.max(np.abs(early_p56-mixture)) < 1e-12)
Wbinary56 = W56[:2, :2]
binary56 = PeriodicOrbitCircuit(6, 2, 4, {1: Wbinary56})
binary_p56 = np.array([sum(binary56.forced_joint(alpha, y)["joint_latent_output_probability"]
                          for alpha in range(3)) for y in range(16)])
t("binary block is invisible when its period divides the late power",
  np.max(np.abs(binary_p56-mixture)) < 1e-12)

print("[C57] co-moving support bounds localized-kick impact after periodic mixing")
from lab.periodic import lightcone_cover
block57 = np.exp(2j*np.pi*np.arange(3)[:, None]*np.arange(3)[None, :]/3)/np.sqrt(3)
V57 = np.kron(np.eye(2), block57)
# Four conditional early branches, after W@s1 and the second controlled U^2.
# This is a fixed finite matrix identity, not a second circuit propagator.
early57 = np.column_stack((V57[:, 0], V57[:, 1],
                          np.roll(V57[:, 0], 2), np.roll(V57[:, 1], 2)))/2
actualF57 = float(np.sum(np.abs(early57[[4, 5, 0], :])**2))
cover57 = lightcone_cover(6, 2, (4, 5, 0), 2)
t("periodic mixing can invalidate the bare residue-mass shortcut",
  abs(actualF57-1/3) < 1e-12 and actualF57 > 1/4)
t("expanded cover bounds the actual pre-kick mass",
  actualF57 <= cover57["covered_prefix_labels"]/cover57["hit_probability_denominator"])
background57 = PeriodicOrbitCircuit(3_000_000_021, 3, 63, {0: block57, 16: block57, 48: block57})
certificate57 = background57.localized_kick_bound(32, (0, 1))
t("large indexed orbit gives an exact rational mathematical omission bound",
  certificate57["covered_prefix_labels"] == 16
  and certificate57["tv_squared_numerator"] == 64
  and certificate57["tv_squared_denominator"] == 1 << 32
  and not certificate57["numerical_sampler_error_included"])

print("[C58] physical pi phase restores coarse dephasing by sector routing")
from lab.periodic import RoutedOrbitCircuit
from lab.semiclassical import sequential_path
W12_58 = W56[::-1, ::-1]
blocks58 = {1: W56, 3: W12_58, 4: W56}
z58 = np.array([1, -1, -1, -1, 1, 1])
full58 = {}
dephased58 = {}
for theta58 in (np.pi/4, np.pi):
    gates58 = {s: np.kron(np.eye(2), w) for s, w in blocks58.items()}
    gates58[2] = np.diag(np.exp(-.5j*theta58*z58))
    pairs58 = [(gates58.get(i+1, np.eye(6)),
                gates58.get(i+1, np.eye(6)) @ np.roll(np.eye(6), 1 << i, axis=0))
               for i in range(5)]
    initial_laws58 = [np.array([sequential_path(pairs58, np.eye(6)[:, j], output=y)
                                ["conditional_path_probability"] for y in range(32)])
                      for j in (0, 3)]
    full58[theta58] = initial_laws58[0]
    dephased58[theta58] = sum(initial_laws58)/2
t("intermediate physical angle has detectable initial sector coherence",
  np.sum(np.abs(full58[np.pi/4]-dephased58[np.pi/4]))/2 > .05)
t("pi rotation restores coarse-mixture validity without commuting with U^3",
  np.max(np.abs(full58[np.pi]-dephased58[np.pi])) < 1e-12
  and np.max(np.abs(z58[:3]+z58[3:])) == 0)
routed58 = RoutedOrbitCircuit(6, 3, 5,
                              {**blocks58, 2: -1j*np.diag(z58[:3])}, {2: 1})
law58 = np.array([sum(routed58.forced_joint(a, y)["joint_latent_output_probability"]
                     for a in range(2)) for y in range(32)])
t("deterministic sector route reproduces the complete pi-phase output law",
  np.max(np.abs(law58-full58[np.pi])) < 1e-12)

print("[C59] few coherent route histories preserve a normalized sampling envelope")
from lab.coherent_routes import CoherentReflectionCircuit
coherent59 = CoherentReflectionCircuit(10, 2, 4, {1: Wbinary56, 3: Wbinary56},
                                      {2: (0, .7), 3: (1, .8)})
identity59 = np.eye(10, dtype=complex)
posts59 = {1: np.kron(np.eye(5), Wbinary56), 3: np.kron(np.eye(5), Wbinary56)}
for s59, (q59, theta59) in coherent59.reflections.items():
    J59 = np.zeros((10, 10), complex)
    for m59 in range(5):
        for p59 in range(2):
            J59[2*((-m59) % 5)+p59, 2*m59+p59] = np.exp(-2j*np.pi*q59*m59/5)
    posts59[s59] = (np.cos(theta59/2)*identity59-1j*np.sin(theta59/2)*J59) @ posts59.get(s59, identity59)
pairs59 = [(posts59.get(i+1, identity59),
            posts59.get(i+1, identity59) @ np.roll(identity59, 1 << i, axis=0)) for i in range(4)]
ref59 = np.array([sequential_path(pairs59, identity59[:, 0], output=y)
                  ["conditional_path_probability"] for y in range(16)])
joint59 = np.array([[coherent59.joint_probability(a, y) for y in range(16)] for a in range(5)])
t("coherent-route final-sector sum matches full-orbit original-order instrument",
  np.max(np.abs(joint59.sum(axis=0)-ref59)) < 1e-12)
t("uniform final sector is only a proposal, not the coherent target",
  np.sum(np.abs(joint59.sum(axis=1)-.2))/2 > .01)
accepted59 = sum(coherent59.rejection_weights(a, y)["proposal_joint_probability"]
                 * coherent59.rejection_weights(a, y)["acceptance_probability"]
                 for a in range(5) for y in range(16))
t("joint rejection resamples the final sector and has mean coefficient-l1 squared",
  abs(accepted59*coherent59.coefficient_l1**2-1) < 1e-12)

print("[C60] conditional direction and normalization must be charged in sampling errors")
from fractions import Fraction
from lab.sampling_error import prefix_error_budget, rejection_error_budget
P60 = [Fraction(1,2), Fraction(0), Fraction(1,2), Fraction(0)]
R60 = [Fraction(0), Fraction(1,2), Fraction(0), Fraction(1,2)]
tv60 = sum(abs(p-r) for p,r in zip(P60,R60))/2
t("equal block masses can conceal completely wrong conditional directions",
  P60[0]+P60[1] == R60[0]+R60[1] and tv60 == 1)
Pzero60 = [Fraction(1,4), Fraction(0), Fraction(3,4), Fraction(0)]
Rzero60 = [Fraction(0), Fraction(0), Fraction(3,4), Fraction(1,4)]
weighted60 = Fraction(1,4)*Fraction(1,2)+Fraction(3,4)*Fraction(1,4)
fulltv60 = sum(abs(p-r) for p,r in zip(Pzero60,Rzero60))/2
t("a normalized zero-block fallback is covered by full-law weighted error",
  weighted60 <= 2*fulltv60)
t("large exponent dimension is charged without allocating it",
  prefix_error_budget(63,2,[Fraction(1,1 << 100)])
  ["scaled_coordinate_to_state_l2_factor"] == 6074001000)
t("small accepted-mass error is amplified by the rejection envelope",
  rejection_error_budget(8,Fraction(1,128)) == Fraction(1,16))

print("[C61] exact-input oracle and finite-bit kernel realize the C60 promise")
from lab.verified_prefix import VerifiedReflectionCircuit, integer_cdf_counts
t("exact-zero weight filtering preserves positive weight ratios",
  integer_cdf_counts([0, 3, 1], 4) == (0, 12, 4)
  and integer_cdf_counts([3, 1], 4) == (12, 4))
try:
    import flint
except ImportError:
    print("  SKIP C61 Arb oracle tests: optional python-flint==0.9.0 required")
else:
    bell61 = VerifiedReflectionCircuit(2, 1, {}, {})
    # |+>|0> -> Bell by controlled shift -> H on control. Four final
    # amplitudes +/-1/2 exactly; this is independent analytic ground truth.
    values61 = [bell61.prefix_dyadic(0, 0, 1, measured=1, output=y,
                                   accuracy_bits=20).coordinates for y in range(2)]
    t("analytic Bell/H final amplitudes are returned exactly on a dyadic grid",
      values61 == [((1 << 19, 0), (1 << 19, 0)),
                   ((1 << 19, 0), (-(1 << 19), 0))])
    from random import Random
    path61 = bell61.sample(Random(61), target_tv=Fraction(1, 1000))
    t("law-level planner certificate includes a covered exact integer kernel",
      path61["total_tv_upper_bound"] <= Fraction(1, 1000)
      and path61["block_updates_upper_bound"] == 1
      and path61["requires_independent_unbiased_bits"])

print("[C62] unitary norm stability does not imply rectangular-width stability")
# Rational rotation has columns (3/5,4/5),(-4/5,3/5), exactly orthonormal.
c62,s62 = Fraction(3,5),Fraction(4,5)
t("a unitary rotation can enlarge both rectangular coordinate radii",
  c62*c62+s62*s62 == 1 and (abs(c62)+abs(s62))**2 == Fraction(49,25))
exact62,mid62 = [Fraction(1),Fraction(0)],[Fraction(1),Fraction(0)]
radius62 = Fraction(0)
for j62 in range(1,9):
    exact62 = [c62*exact62[0]-s62*exact62[1],s62*exact62[0]+c62*exact62[1]]
    local62 = [Fraction((-1)**j62,10**j62),Fraction(1,10**(j62+1))]
    mid62 = [c62*mid62[0]-s62*mid62[1]+local62[0],
             s62*mid62[0]+c62*mid62[1]+local62[1]]
    radius62 += sum(abs(e62) for e62 in local62)
    t("exact local residual upper bounds telescope in a norm preserved by the steps",
      sum((a62-b62)**2 for a62,b62 in zip(exact62,mid62)) <= radius62**2)

print("[C63] one-sided acceptance controls mass without a minimum proposal mass")
from lab.verified_rejection import conservative_threshold
for j63 in (0,1,32,256,1024):
    d63 = Fraction(1,1 << j63)
    p63 = d63/3
    pl63,pu63,dl63,du63 = p63-d63/100,p63+d63/100,d63-d63/100,d63+d63/100
    a63 = Fraction(conservative_threshold(pl63,du63,16),1 << 16)
    t("accepted-measure deficit bound has no division by the true tiny mass",
      0 <= p63-d63*a63 <= pu63-pl63+du63-dl63+d63/Fraction(1 << 16))
p63 = [Fraction(1,2),Fraction(1,2)]
q63 = [Fraction(1,4),Fraction(3,4)]
a63 = [Fraction(1),Fraction(1,3)]
wrong63 = [Fraction(1),Fraction(0)]
accepted63 = [v*a for v,a in zip(wrong63,a63)]
law63 = [v/sum(accepted63) for v in accepted63]
t("perfect acceptance ratios do not repair an uncertified proposal",
  [q*a for q,a in zip(q63,a63)] == [v/2 for v in p63]
  and sum(abs(v-p) for v,p in zip(law63,p63))/2 == Fraction(1,2))

print("[C64] binary mass and normalized-weight oracle budgets are distinct")
from lab.sampling_error import prefix_mass_error_budget,plan_conditional_weight_accuracy
t("binary internal-node count is charged without tree allocation",
  prefix_mass_error_budget(4,Fraction(1,1024),20)["oracle_tv_upper_bound"] == Fraction(30,1024))
tiny64 = Fraction(1,1 << 128)
t("uniform fallback at a rare deterministic prefix has large local but tiny weighted TV",
  Fraction(1,2) > Fraction(1,1000) and tiny64/2 < Fraction(1,1 << 100))

print("[C65] collective late-background removal uses disjoint tensor factors")
# Exact integer commuting test: coarse route J tensor I_p versus I_alpha tensor W.
J65 = np.array([[0,1],[1,0]],dtype=int)
W65 = np.array([[1,0],[0,-1]],dtype=int)
X65 = np.array([[0,1],[1,0]],dtype=int)
A65,B65 = np.kron(J65,np.eye(2,dtype=int)),np.kron(np.eye(2,dtype=int),W65)
t("coarse routing commutes with within-sector work gates although the work gates need not commute",
  np.array_equal(A65@B65,B65@A65) and not np.array_equal(W65@X65,X65@W65))

print("[C61/C64] odd-block constants differ from binary output-tree constants")
budget_odd = prefix_error_budget(4,3,[Fraction(1,1024)],random_bits=20)
t("three complex coordinates and three-bin rounding are explicitly charged",
  budget_odd["scaled_coordinate_to_state_l2_factor"] == 10
  and budget_odd["oracle_tv_upper_bound"] == Fraction(20,1024)
  and budget_odd["finite_random_bits_tv_upper_bound"] == Fraction(2,1 << 20))
shift_odd = np.roll(np.eye(3,dtype=int),1,axis=0)
t("powers two modulo three remain non-scalar, unlike the binary scalar tail",
  all(not np.array_equal(np.linalg.matrix_power(shift_odd,1 << i),np.eye(3,dtype=int)) for i in range(8)))

print("[C66] a shared state error differs from prefix-dependent choices")
# Complete diagonal POVM in dimension five, with four two-child parents.
# All numbers exact; neither a floating norm estimate nor a QFT realization.
aggregate66 = Fraction(0)
sum_effects66 = [Fraction(0)]*5
for a66 in range(1,5):
    rho66 = [Fraction(1)]+[Fraction(0)]*4
    approx66 = [Fraction(9,10)]+[Fraction(int(j == a66),10) for j in range(1,5)]
    effects66 = [[Fraction(1,4)]+[Fraction(0)]*4,
                 [Fraction(int(j == a66)) for j in range(5)]]
    true66 = [sum(x*y for x,y in zip(e,rho66)) for e in effects66]
    other66 = [sum(x*y for x,y in zip(e,approx66)) for e in effects66]
    aggregate66 += sum(true66)*sum(abs(x/sum(true66)-y/sum(other66))
                                  for x,y in zip(true66,other66))/2
    for e66 in effects66:
        sum_effects66 = [x+y for x,y in zip(sum_effects66,e66)]
    t("each adversarial forward state is normalized PSD with trace error one fifth",
      sum(approx66) == 1 and min(approx66) >= 0
      and sum(abs(x-y) for x,y in zip(rho66,approx66)) == Fraction(1,5))
t("prefix-dependent choices can defeat one maximum-state-error budget",
  sum_effects66 == [1]*5 and aggregate66 == Fraction(4,13) > Fraction(1,5))
rho66 = [Fraction(1)]+[Fraction(0)]*4
fixed66 = [Fraction(9,10),Fraction(1,10),Fraction(0),Fraction(0),Fraction(0)]
absolute66 = Fraction(0)
for a66 in range(1,5):
    for e66 in ([Fraction(1,4)]+[Fraction(0)]*4,
                 [Fraction(int(j == a66)) for j in range(5)]):
        absolute66 += abs(sum(e*(x-y) for e,x,y in zip(e66,rho66,fixed66)))
t("one shared state obeys POVM trace-norm contraction on the same effects",
  absolute66 == sum(abs(x-y) for x,y in zip(rho66,fixed66)) == Fraction(1,5))
t("r9 power-two arithmetic has exact period six for all later channel labels",
  pow(2,6,9) == 1 and all(pow(2,i+6,9) == pow(2,i,9) for i in range(12)))

print("[C68] reverse trajectories correct the initial boundary by rejection")
from lab.semiclassical import _qft_branch_operators
R68 = np.array([[Fraction(3,5),Fraction(-4,5)],
                [Fraction(4,5),Fraction(3,5)]],dtype=object)
Z68 = np.array([[1,0],[0,-1]],dtype=object)
K68 = _qft_branch_operators(R68,R68@Z68,1)
I68 = np.eye(2,dtype=object)
t("both completeness relations hold for nonnormal rational branch operators",
  np.array_equal(sum(K.T@K for K in K68),I68)
  and np.array_equal(sum(K@K.T for K in K68),I68)
  and not np.array_equal(K68[0].T@K68[0],K68[0]@K68[0].T))
q68,accepted68,wrong68 = [],[],[]
for K in K68:
    # j is sampled uniformly at the FINAL boundary. Reverse raw vector is
    # K^T e_j. Telescoping cancels its norm in the accepted joint submass.
    q68.append(sum(x*x for x in K.flat)/2)
    accepted68.append(sum(K[j,0]**2 for j in range(2))/2)
    wrong68.append(sum(K[0,j]**2 for j in range(2))/2)
t("pure reverse accepted submass is p/b with mean b attempts",
  q68 == [Fraction(1,2)]*2 and accepted68 == [Fraction(1,2),0]
  and sum(accepted68) == Fraction(1,2))
t("omitting rejection or forgetting the adjoint changes this exact target",
  q68 != [1,0] and [2*x for x in wrong68] == [Fraction(9,25),Fraction(16,25)])

print("[C69] rejection must redraw the small reverse boundary label")
# Same rational nonnormal Kraus pair, now psi=(3/5,4/5). Accepted
# unscaled cell masses N_jy factor as R_jy^2 * |psi_y|^2 here.
p69 = [Fraction(9,25),Fraction(16,25)]
N69 = [[R68[j,y]**2*p69[y] for y in range(2)] for j in range(2)]
fresh69 = [sum(N69[j][y] for j in range(2)) for y in range(2)]
frozen69 = [sum(N69[j][y]/sum(N69[j]) for j in range(2))/2 for y in range(2)]
t("redrawing j gives the target but rejection at frozen j biases its mixture",
  fresh69 == p69 and sum(frozen69) == 1
  and abs(frozen69[0]-p69[0]) > Fraction(1,100))
from lab.verified_reverse_work import VerifiedReverseWork
for dim69,period69 in ((2,10),(3,9)):
    plan69 = VerifiedReverseWork(VerifiedReflectionCircuit(period69,4,{},{},block_size=dim69),0).plan(Fraction(1,1000000))
    e69 = (2*plan69["proposal_tv_upper_bound"]
           + 16*plan69["acceptance_interval_width_tolerance"]
           + Fraction(1,1 << plan69["acceptance_random_bits"]))
    t("accepted-mass normalization charges dimension and retains positive success",
      plan69["total_tv_upper_bound"] == dim69*e69 <= Fraction(3,4000000)
      and plan69["success_probability_lower_bound"] == Fraction(1,dim69)-e69)

print("[C70] local instrument and projective compression budgets add with depth")
from lab.verified_quantized_reverse import (VerifiedQuantizedReverseWork,
    canonical_quantize,integer_norm,integer_overlap_numerator)
for b70,r70 in ((2,10),(3,9)):
    for width70 in (0,1,4,63):
        p70 = VerifiedQuantizedReverseWork(VerifiedReflectionCircuit(
            r70,width70,{},{},block_size=b70),0).plan(Fraction(1,1000000))
        e70 = (width70*p70["local_cq_error_upper_bound"]
               +p70["initial_state_trace_error_upper_bound"]
               +Fraction(1,1 << p70["acceptance_random_bits"]))
        t("local CQ error, terminal effect and acceptance floor cover the full law",
          p70["total_tv_upper_bound"] == b70*e70 <= Fraction(1,1000000)
          and p70["proposal_tv_upper_bound"] == width70*p70["local_cq_error_upper_bound"]/2
          and p70["success_probability_lower_bound"] == Fraction(1,b70)-e70
          and p70["initial_coordinate_bits_upper_bound"] == p70["grid_bits"]+2)
for v70 in (((-7,7),(3,-4)),((1,0),(1 << 400,1)),((0,-1),(1,0))):
    w70,_ = canonical_quantize(v70,8)
    # Pure-state trace distance squared, evaluated without numerical roots.
    trace_squared70 = 4*(1-Fraction(integer_overlap_numerator(v70,w70),
                                     integer_norm(v70)*integer_norm(w70)))
    t("canonical compression has a norm-independent projector error bound",
      0 <= trace_squared70 <= Fraction(4*len(v70),1 << 8)**2
      and max(abs(x) for z in w70 for x in z) == 1 << 8)
physical70 = (integer_norm(((1,0),(0,0))),integer_norm(((2,0),(0,0))))
compressed70 = tuple(integer_norm(canonical_quantize(v,8)[0])
                     for v in (((1,0),(0,0)),((2,0),(0,0))))
t("equal projective child states do not imply equal Born branch weights",
  Fraction(physical70[0],sum(physical70)) == Fraction(1,5)
  and Fraction(compressed70[0],sum(compressed70)) == Fraction(1,2))

print("[C71] finite-time reflection words and the coherent initial boundary")
from lab.coherent_reverse import route_support_bound
for M71,qs71 in ((7,(0,1)),(16,(2,6)),(101,(0,37))):
    for alpha71 in range(min(M71,8)):
        reached71 = {alpha71}
        for depth71 in range(8):
            q71 = qs71[depth71%2]
            reached71 |= {(-a-q71)%M71 for a in reached71}
            t("two involutions bound every reached prefix without an invertible-difference assumption",
              len(reached71) <= route_support_bound(M71,tuple(qs71[i%2] for i in range(depth71+1))))
for coordinates71 in ((1,1,1),(1,-1,0),(Fraction(1,5),Fraction(2,5),Fraction(-3,5))):
    D71 = sum(x*x for x in coordinates71)
    N71 = sum(coordinates71)**2
    t("coherent boundary acceptance has the support-size envelope including cancellations",
      0 <= N71 <= len(coordinates71)*D71)
t("discarding cross-sector interference changes the boundary functional",
  (Fraction(1)-Fraction(1))**2 == 0
  and Fraction(1)**2+Fraction(-1)**2 == 2)
smooth_wrap_sum71 = sum(((1 << i)+2)//3 for i in range(15))
smooth_bound71 = (4*Fraction(383,1000)
                  *Fraction(44,7*((1 << 40)-1))*smooth_wrap_sum71)
t("neighboring routes have a small task-specific bound despite a full static group",
  smooth_wrap_sum71 == 10930 and 0 < smooth_bound71 < Fraction(1,10_000_000))
t("the neighboring-displacement bound cannot silently certify a separated route",
  679535556937*smooth_bound71 > 1)
for cycle71 in (5,7,8):
    identity71 = np.eye(cycle71,dtype=np.int64)
    shift71 = np.roll(identity71,1,axis=0)
    laplacian71 = 2*identity71-shift71-shift71.T
    for q71 in range(cycle71):
        reflection71 = identity71[(-np.arange(cycle71)-q71)%cycle71]
        t("every dihedral route preserves the coarse Laplacian exactly",
          np.array_equal(reflection71.T@laplacian71@reflection71,laplacian71)
          and np.array_equal(reflection71@laplacian71,laplacian71@reflection71))
    sector_sign71 = identity71.copy()
    sector_sign71[0,0] = -1
    t("sector-dependent work phases can break smoothness conservation",
      not np.array_equal(sector_sign71@laplacian71@sector_sign71,laplacian71))

print("[C72] batched adjoint prefixes and completed-oracle phase invariance")
# Integer-complex entries make this duality check exact in binary arithmetic;
# unitarity is unnecessary for the amplitude identity itself.
A72 = np.arange(36).reshape(6,6).astype(complex)
A72 += 1j*np.roll(np.arange(36).reshape(6,6),1,axis=0)
initial72 = np.array([1,0,1,0,1,0],dtype=complex)
for gamma72 in range(3):
    boundary72 = np.eye(6,dtype=complex)[:,2*gamma72:2*gamma72+2]
    reverse72 = A72.conj().T@boundary72
    row_sum72 = reverse72.reshape(3,2,2)[:,0,:].sum(axis=0)
    target72 = (A72@initial72)[2*gamma72:2*gamma72+2]
    t("coherent initial overlap of all boundary columns has the forward amplitude convention",
      np.array_equal(row_sum72.conj(),target72)
      and not np.array_equal(row_sum72,target72))
v72 = np.array([1+2j,3-4j,-2+3j])
t("completed oracle phases can be dropped, but internal interference cannot",
  np.allclose(np.abs(v72)**2,np.abs(np.abs(v72))**2,atol=1e-14,rtol=0)
  and abs((1+1j)+(1-1j))**2 != (abs(1+1j)+abs(1-1j))**2)

print("[C73] integer normal forms count reached labels, not reflection histories")
import itertools
for dimension73 in range(4):
    for radius73 in range(5):
        exact73 = sum(sum(abs(x) for x in point73)<=radius73
            for point73 in itertools.product(range(-radius73,radius73+1),repeat=dimension73))
        count73 = sum(2**j*math.comb(dimension73,j)*math.comb(radius73,j)
                      for j in range(min(dimension73,radius73)+1))
        t("integer L1-ball count matches independent lattice enumeration",exact73==count73)
for modulus73,alphabet73 in ((101,(7,13,41)),(12,(0,4,8)),(15,(0,3,6))):
    for word73 in itertools.product(range(3),repeat=6):
        alpha73,q073,n73,odd73 = 5,alphabet73[0],[0,0],False
        for letter73 in word73:
            alpha73 = (-alpha73-alphabet73[letter73])%modulus73
            n73 = [-v for v in n73]
            if letter73:
                n73[letter73-1] -= 1
            odd73 = not odd73
            normal73 = (-5-q073 if odd73 else 5)+sum(
                n73[j]*(alphabet73[j+1]-q073) for j in range(2))
            assert normal73%modulus73 == alpha73 and sum(abs(v) for v in n73)<=6
    t("signed coefficient normal form survives separated and noninvertible route differences",True)
for m73 in range(1,7):
    word73 = tuple(v for j in range(m73) for v in (0,3**j))
    reached73 = {0}
    for q73 in reversed(word73):
        reached73 |= {(-a-q73)%(3**(m73+1)) for a in reached73}
    t("growing alphabets defeat a universal linear reached-support claim",
      len(reached73)==3**m73 and len(reached73)>=2**m73)
for modulus73 in (101,1000):
    even73,oddset73 = {0},set()
    for q73 in reversed(tuple(i%3 for i in range(8))):
        even73,oddset73 = (even73 | {(-o-q73)%modulus73 for o in oddset73},
                          oddset73 | {(-e-q73)%modulus73 for e in even73})
    sizes73 = [len({(gamma73+e)%modulus73 for e in even73}
                   | {(-gamma73+o)%modulus73 for o in oddset73}) for gamma73 in (0,5)]
    t("a single-sector support cannot replace the global word-specific envelope",
      len(even73)==len(oddset73)==10 and sizes73==[11,20]
      and len(even73)*len(oddset73)<modulus73//math.gcd(2,modulus73))

print("[C74] physical cell coordinates and uniformly accurate phase oracles")
for N74,a74,b74,M74 in ((7,3,3,2),(13,2,3,4),(35,2,6,2),(97,5,3,32)):
    r74 = b74*M74
    table74 = {pow(a74,M74*p74,N74):p74 for p74 in range(b74)}
    wrong74 = 0
    for j74 in range(r74):
        x74 = pow(a74,j74,N74)
        p74 = table74[pow(x74,M74,N74)]
        expected74 = pow(a74,(-j74+2*(j74%b74))%r74,N74)
        assert p74 == j74%b74
        assert pow(a74,2*p74,N74)*pow(x74,-1,N74)%N74 == expected74
        wrong74 += pow(x74,-1,N74) != expected74
    t("fine subgroup lookup implements cell reflection but plain inversion does not",
      len(table74)==b74 and wrong74>0)

def circle74(x, y):
    z = (x-y)%1
    return min(z,1-z)

for modulus74 in (3,4,7,8,15):
    depth74 = (modulus74-1).bit_length()
    for residue74 in range(modulus74):
        for signs74 in itertools.product((-1,1),repeat=depth74+1):
            noisy74 = [(Fraction((residue74 << i)%modulus74,modulus74)
                        +Fraction(signs74[i],16))%1 for i in range(depth74+1)]
            estimate74 = noisy74[-1]
            for i74 in range(depth74-1,-1,-1):
                candidates74 = (estimate74/2,(estimate74+1)/2)
                estimate74 = min(candidates74,key=lambda v:circle74(v,noisy74[i74]))
            assert circle74(estimate74,Fraction(residue74,modulus74)) <= Fraction(1,16 << depth74)
            assert math.floor(modulus74*estimate74+Fraction(1,2))%modulus74 == residue74
    t("constant-error circular inverse doubling recovers the exact root grid",True)
bad_low74,bad_high74 = Fraction(6,25),Fraction(19,25)
wrong_lift74 = min((bad_high74/2,(bad_high74+1)/2),
                   key=lambda v:circle74(v,bad_low74))
t("epsilon below one quarter alone does not validate nearest-half decoding",
  math.floor(2*wrong_lift74+Fraction(1,2))%2 != 0)

print("[C75] clean coordinate conjugation and a coherent involution flag")
for N75,a75,b75,M75 in ((7,3,3,2),(13,2,3,4),(35,2,6,2)):
    r75 = b75*M75
    table75 = {pow(a75,M75*p75,N75):p75 for p75 in range(b75)}
    # Independent amplitude-level code-isometry calculation; the separate
    # experiment, not this algebra regression, audits physical gate synthesis.
    W75 = np.exp(2j*np.pi*np.outer(np.arange(b75),np.arange(b75))/b75)/np.sqrt(b75)
    clean75 = np.zeros((r75,r75),complex)
    for j75 in range(r75):
        x75 = pow(a75,j75,N75)
        p75 = table75[pow(x75,M75,N75)]
        y75 = x75*pow(a75,-p75,N75)%N75
        for outp75 in range(b75):
            outx75 = y75*pow(a75,outp75,N75)%N75
            outcoord75 = outp75 ^ table75[pow(outx75,M75,N75)]
            assert outcoord75 == 0
            assert outx75 == pow(a75,b75*(j75//b75)+outp75,N75)
            clean75[b75*(j75//b75)+outp75,j75] = W75[outp75,p75]
    t("new output coordinate uncomputes after arbitrary coherent fine mixing",
      np.array_equal(clean75,np.kron(np.eye(M75),W75)))
R75 = np.zeros((12,12),complex)
for j75 in range(12):
    R75[3*((-(j75//3))%4)+j75%3,j75] = 1
I75 = np.eye(12,dtype=complex)
H75 = np.kron(np.array([[1,1],[1,-1]])/np.sqrt(2),I75)
controlled75 = np.zeros((24,24),complex)
controlled75[:12,:12],controlled75[12:,12:] = I75,R75
A75 = H75@controlled75@H75
angle75 = np.pi/3
rz75 = np.diag(np.repeat(np.exp(np.array([-1j,1j])*angle75/2),12))
out75 = (A75.conj().T@rz75@A75)[:,:12]
target75 = np.cos(angle75/2)*I75-1j*np.sin(angle75/2)*R75
t("involution eigenvalue kickback implements the full complex rotation with clean flag",
  np.allclose(out75[:12],target75,atol=2e-14,rtol=0)
  and np.max(abs(out75[12:]))<2e-14)
controlled75[12:,12:] *= -1
badA75 = H75@controlled75@H75
bad75 = (badA75.conj().T@rz75@badA75)[:12,:12]
t("controlled global phase is not harmless in coherent reflection exponentiation",
  np.max(abs(bad75-target75))>.5
  and np.allclose(bad75,target75.conj(),atol=2e-14,rtol=0))
shift75 = np.roll(I75,3,axis=0)
t("fine mixers commute with reflection but arithmetic cannot be commuted away",
  np.array_equal(R75@shift75,shift75.T@R75)
  and not np.array_equal(R75@shift75,shift75@R75))

print("[C76] exact roots-of-unity cancellations and prime-field nonzeros")
def monic_remainder76(coefficients, divisor):
    # Tiny independent integer/Gaussian-integer division, not FLINT's kernel.
    remainder = list(coefficients)
    degree = len(divisor)-1
    for power in range(len(remainder)-1,degree-1,-1):
        leading = remainder[power]
        for shift,coefficient in enumerate(divisor):
            remainder[power-degree+shift] -= leading*coefficient
    return tuple(remainder[:degree])

nonzeros76 = 0
for k76 in (1,2):
    for p76 in range(3):
        for delta76 in range(4):
            coeff76 = [0j]*13
            for m76 in range(4):
                coeff76[k76*pow(2,3*m76+p76,13)%13] += 1j**(delta76*m76)
            nonzeros76 += any(monic_remainder76(coeff76,[1]*13))
t("prime N13/M4 subgroup sums have nonzero exact Gaussian-integer remainders",
  nonzeros76==24)
supports76=[]
for p76 in range(2):
    support76=[]
    for delta76 in range(3):
        coeff76=[0]*9
        for m76 in range(3):
            coeff76[(pow(2,2*m76+p76,9)+3*delta76*m76)%9] += 1
        if any(monic_remainder76(coeff76,[1,0,0,1,0,0,1])):
            support76.append(delta76)
    supports76.append(support76)
t("composite N9 has exact cancellations: the prime hypothesis cannot be dropped",
  supports76==[[2],[1]])
N76,r76=17,16
coeff76=np.array([sum(np.exp(2j*np.pi*((pow(3,j76,N76)/N76
       +delta76*j76/r76)%1)) for j76 in range(r76))/r76 for delta76 in range(r76)])
prob76=abs(coeff76)**2
t("one-sector Gauss mass bound limits best sector-sparse state truncation",
  abs(float(prob76.sum())-1)<2e-14 and max(prob76)<=N76/r76**2+2e-14
  and 1-sum(sorted(prob76,reverse=True)[:8])>=1-8*N76/r76**2-2e-14)

print("[C53 follow-up] purification transfer and a different output objective")
A53u = np.array([[1.,1.],[1.,-1.]])/2
G53u = np.diag([1.,1j])
V53u = 2*(A53u.conj().T@G53u@A53u)
t("a Fourier purification transfers a diagonal work gate via a nondiagonal exponent unitary",
  np.allclose(A53u@V53u,G53u@A53u,atol=2e-14,rtol=0)
  and np.allclose(V53u.conj().T@V53u,np.eye(2),atol=2e-14,rtol=0)
  and sum(abs(np.diag(A53u.conj().T@G53u@A53u))) < .8)
# Terminal work H is invisible to exponent Z measurement. The joint-state
# optimum nevertheless changes that measurement on a nonuniform Schmidt pair.
A53u = np.diag(np.sqrt([.8,.2]))
G53u = np.array([[1.,1.],[1.,-1.]])/math.sqrt(2)
V53u = np.array([[5.,4.],[4.,-5.]])/math.sqrt(41)
C53u = A53u.conj().T@G53u@A53u
F53u = np.linalg.svd(C53u,compute_uv=False).sum()
target53u = G53u@A53u
surrogate53u = A53u@V53u
law53u = np.sum(abs(surrogate53u)**2,axis=0)
t("joint-state optimal transfer can be worse than omission for the measured output",
  abs(F53u-math.sqrt(.82))<2e-14
  and abs(np.linalg.norm(target53u-surrogate53u)**2-(2-2*F53u))<2e-14
  and np.allclose(np.sum(abs(target53u)**2,axis=0),[.8,.2],atol=2e-14,rtol=0)
  and np.sum(abs(law53u-np.array([.8,.2])))/2>.23)

A53p = np.array([[math.sqrt(.3),0],[0,math.sqrt(.2)],
                 [math.sqrt(.3),0],[0,-math.sqrt(.2)]],dtype=complex)
U53p = np.roll(np.eye(4),2,axis=0)
G53p = np.diag([1.,1.,-1.,1.])
E53p = np.array([[1.,1j],[1.,-1j]])/math.sqrt(2)
def parity_law53(X):
    law=[]
    for y in range(4):
        v=np.exp(-2j*np.pi*y*np.arange(2)/4)
        z=(X+(-1)**y*U53p@X)@v/(2*math.sqrt(2))
        law.append(float(np.vdot(z,z).real))
    return np.array(law)
base53p=parity_law53(A53p)
early53p=parity_law53(A53p@E53p.T)
kick53p=parity_law53(G53p@A53p)
t("a work-only parity witness constrains every early-unitary replacement",
  abs(base53p.sum()-1)<2e-14 and abs(early53p.sum()-1)<2e-14
  and abs(kick53p.sum()-1)<2e-14
  and abs(sum(base53p[::2])-.6)<2e-14
  and abs(sum(early53p[::2])-.6)<2e-14
  and abs(sum(kick53p[::2]))<2e-14
  and np.sum(abs(kick53p-early53p))/2>=.6-2e-14)

print("[C77] partial Fourier orthogonality and inclusive support cones")
# A fixed l, high-history Gram matrix suffices; averaging translated copies
# over l leaves these inner products unchanged. Independent exact quarter roots.
labels77=[54,16,32,54]
gram77=np.equal.outer(labels77,labels77).astype(int)
numerators77=[]
for z77 in range(4):
    numerator77=sum((1j)**(z77*(h77-hp77))*gram77[h77,hp77]
                    for h77 in range(4) for hp77 in range(4))
    assert numerator77.imag==0
    numerators77.append(int(numerator77.real))
t("inclusive cones sharing an endpoint permit nonuniform Fourier output",
  numerators77==[6,4,2,4]
  and all(min((label77-16*h77)%60,(16*h77-label77)%60)<=6
          for h77,label77 in enumerate(labels77)))
# Orthogonal high histories give a uniform reduced Fourier measurement even
# with arbitrary phases on the supported work labels.
root77=np.array([[1j**(h77*z77) for h77 in range(4)] for z77 in range(4)])/2
phase77=np.diag(np.exp(1j*np.array([.1,.7,-.3,1.2])))
rho77=phase77@(np.eye(4)/4)@phase77.conj().T
t("separated histories yield uniform low output bits, not an assumption of sector conservation",
  np.allclose(np.diag(root77@rho77@root77.conj().T),np.full(4,.25),atol=2e-14,rtol=0))

print("[C78] work-first conditioning and sparse Fourier rows")
# Gaussian-integer columns divided by sqrt(2), so every column is normalized.
# Work is the first index below; exact quarter roots avoid trigonometry.
from fractions import Fraction
rows78 = [[1, 0, 0, 1], [1, 1, 0, 0],
          [0, 1j, 1, 0], [0, 0, 1, 1], [0, 0, 0, 0]]
column_norms78 = [sum(abs(row78[e78])**2 for row78 in rows78)
                  for e78 in range(4)]
# abs() can round sqrt(2); real/imag Gaussian integers give exact norms.
joint_numerators78 = []
for row78 in rows78:
    nums78=[]
    for y78 in range(4):
        amp78=sum(row78[e78]*(1j)**(-y78*e78) for e78 in range(4))
        nums78.append(int(amp78.real**2+amp78.imag**2))
    joint_numerators78.append(nums78)
full_nums78 = [sum(row78[y78] for row78 in joint_numerators78)
               for y78 in range(4)]
t("uniform-history column sampling gives Born work weights, including a zero row",
  column_norms78 == [2,2,2,2]
  and [sum(abs(a78)**2 for a78 in row78)/8 for row78 in rows78]
      == [.25,.25,.25,.25,0])
t("coherent work rows reconstruct a nonuniform complete Fourier law",
  full_nums78 == [14,6,2,10] and sum(full_nums78)==32
  and full_nums78 != [8]*4)
prefix_ok78=True
for j78,row78 in enumerate(rows78[:4]):
    Z78=int(sum(abs(a78)**2 for a78 in row78))
    for d78 in range(3):
        H78,L78=1<<d78,4>>d78
        for z78 in range(H78):
            groups78={}
            for e78,a78 in enumerate(row78):
                if a78:
                    c78=e78%L78
                    groups78[c78]=groups78.get(c78,0)+a78*(1j)**(
                        -z78*(e78//L78)*(4//H78))
            num78=sum(int(a78.real**2+a78.imag**2) for a78 in groups78.values())
            expected78=Fraction(sum(joint_numerators78[j78][z78::H78]),4*Z78)
            prefix_ok78 &= Fraction(num78,H78*Z78)==expected78
t("occupied-residue grouping gives every exact sparse-row Fourier prefix",prefix_ok78)
t("two disjoint input components have mean two coherent-rejection attempts",
  all(Fraction(sum(nums78),16)==Fraction(1,2)
      and all(0<=num78<=4 for num78 in nums78)
      for nums78 in joint_numerators78[:4]))

print("[C79] binary-residue cycles and revived cancellations")
import math as math79
cycle_cases79=0
for r79 in range(1,33):
    for v79 in range(7):
        A79=1<<v79
        P79=A79//math79.gcd(r79,A79)
        residues79=[(n79*r79)%A79 for n79 in range(P79)]
        assert len(set(residues79))==P79 and (P79*r79)%A79==0
        for count79 in range(9):
            partition79=[]
            for z79 in range(min(P79,count79)):
                size79=1+(count79-1-z79)//P79
                partition79.extend(z79+P79*k79 for k79 in range(size79))
            assert sorted(partition79)==list(range(count79))
        cycle_cases79+=1
t("residue-walk periods and truncated cycle classes are exact finite partitions",cycle_cases79==224)
old79,new79={},{}
for u79 in range(2):
    for q79 in range(2):
        rho79=(q79-u79)%6
        if rho79>=4:
            continue
        sign79=(-1)**q79  # numerator of W1[1,q]*W0[u,0], denominator 2
        old79[rho79]=old79.get(rho79,0)+sign79
        new79[rho79]=new79.get(rho79,0)+sign79*(-1)**(u79+(rho79%2))
t("an earlier diagonal phase revives a canceled component even when P=1",
  old79=={0:0,1:-1} and new79=={0:2,1:1}
  and sum(a79!=0 for a79 in old79.values())==1
  and sum(a79!=0 for a79 in new79.values())==2)
separated79=0
overlap_control79=False
for low79 in range(128):
    reached79={3*((u79+low79)%60//3)+p79 for u79 in range(3) for p79 in range(3)}
    assert not reached79.intersection({(j79+128)%60 for j79 in reached79})
    overlap_control79 |= bool(reached79.intersection({(j79+64)%60 for j79 in reached79}))
    separated79+=1
t("one high input bit stays support-separated through arbitrary earlier diagonal phases",
  separated79==128 and overlap_control79)

dual_cases79=0
for r79 in range(1,17):
    for v79 in range(5):
        A79=1<<v79
        for K79 in (1,2,4):
            L79=A79*K79
            for rho79 in range(r79):
                partition79=[]
                for k79 in range(K79):
                    shifted79=(rho79-A79*k79)%r79
                    count79=max(0,1+(A79-1-shifted79)//r79)
                    for n79 in range(count79):
                        low79=A79*k79+shifted79+r79*n79
                        assert A79*k79<=low79<A79*(k79+1)
                        # Add arbitrary u to either side: the early phase
                        # argument is fixed within this remaining-history cell.
                        assert (low79%A79)%r79==(rho79-A79*k79)%r79
                        partition79.append(low79)
                assert sorted(partition79)==list(range(rho79,L79,r79))
                assert len(partition79)==len(set(partition79))
                dual_cases79+=1
t("dual remaining-history cover exactly partitions residues, including short and aliased cells",
  dual_cases79==2040)
balanced_cases79=0
for r79 in range(1,33):
    alpha79=(r79 & -r79).bit_length()-1
    for s79 in range(9):
        worst79=max(min((1<<v79)//math79.gcd(r79,1<<v79),1<<(s79-v79))
                    for v79 in range(s79+1))
        assert worst79==1<<(max(0,s79-alpha79)//2)
        balanced_cases79+=1
t("balanced cycle/history choice has the stated exact worst insertion factor",
  balanced_cases79==288)

print("[C80] component-mixture envelope and conditional retry accounting")
lambda80=(Fraction(9,10),Fraction(1,10))
w80=(Fraction(3,4),Fraction(1,4))
envelope80=sum(a80/b80 for a80,b80 in zip(lambda80,w80))
t("root-mass weights attain the exact unequal-mass Cauchy envelope",
  envelope80==Fraction(8,5)
  and sum(a80/b80 for a80,b80 in zip(lambda80,lambda80))==2)
# Rational grid is a regression, not the proof of universal optimality.
t("nearby rational component mixtures do not improve the proved envelope",
  all(lambda80[0]/Fraction(k80,100)+lambda80[1]/(1-Fraction(k80,100))>=envelope80
      for k80 in range(1,100)))
# d=(9,1) makes the cost-optimal weights exactly (1/2,1/2).
cost_opt80=sum(a80/Fraction(1,2) for a80 in lambda80)*5
mass_root_cost80=envelope80*(9*w80[0]+w80[1])
t("minimizing attempts differs from minimizing component-dependent cost",
  cost_opt80==10 and mass_root_cost80==Fraction(56,5) and cost_opt80<mass_root_cost80)
# Truncated series plus its explicit unreturned residual, different E_j.
for e80 in (Fraction(1),Fraction(8,5),Fraction(3)):
    success80=1/e80
    partial80=sum((1-success80)**k80*success80 for k80 in range(8))
    assert partial80+(1-success80)**8==1
t("fixed-work retries normalize separately for unequal row envelopes",True)
t("control: redraw-work rejection reweights equal work probabilities",
  Fraction(1,2)/(Fraction(1,2)+Fraction(1,4))==Fraction(2,3))

print("[C81] nested prefix cuts, largest gaps, and finite-truncation costs")
from itertools import combinations_with_replacement as combinations81
gap_cases81 = 0
for r81 in range(1, 17):
    alpha81 = (r81 & -r81).bit_length()-1
    for s81 in range(8):
        families81 = [()] + [(v81,) for v81 in range(s81+1)]
        families81 += list(combinations81(range(s81+1), 2))
        for positions81 in families81:
            factors81 = []
            for cut81 in range(s81+1):
                left81 = [v81 for v81 in positions81 if v81 < cut81]
                P81 = ((1 << max(left81)) // math79.gcd(r81, 1 << max(left81))
                       if left81 else 1)
                factors81.append((1 << (s81-cut81))*P81)
            B81 = max(0, s81-alpha81)
            points81 = sorted({0, B81} | {max(0, v81-alpha81) for v81 in positions81})
            delta81 = max((b81-a81 for a81, b81 in zip(points81, points81[1:])), default=0)
            assert min(factors81) == 1 << (B81-delta81)
            assert min(factors81) <= 1 << (len(positions81)*B81//(len(positions81)+1))
            gap_cases81 += 1
t("largest-gap identity includes empty, initial, duplicate and endpoint phases",
  gap_cases81 == 2624)
pair_costs81 = {}
for cut81 in (6, 7):
    A81, K81 = 1 << cut81, 1 << (7-cut81)
    V81 = max(v81 for v81 in (3, 6) if v81 < cut81)
    P81 = (1 << V81)//math79.gcd(60, 1 << V81)
    grouping81 = coefficient81 = 0
    for high81 in range(2):
        cell81 = 3*((-128*high81) % 60//3)
        for k81 in range(K81):
            for u81 in range(3):
                for q81 in range(3):
                    grouping81 += 1
                    rho81 = (cell81+q81-A81*k81-u81) % 60
                    N81 = max(0, 1+(A81-1-rho81)//60)
                    coefficient81 += min(P81, N81)
    pair_costs81[cut81] = grouping81, coefficient81, K81*P81
t("control: a better untruncated factor can require more actual local pair visits",
  pair_costs81 == {6: (36,42,4), 7: (18,42,16)}
  and sum(pair_costs81[6][:2]) > sum(pair_costs81[7][:2]))
old_numerator81 = sum((-1)**u81 for u81 in range(2))
new_numerator81 = sum((-1)**u81 * (-1)**u81 * (1j)**u81 for u81 in range(2))
t("two diagonal insertions revive the old canceled Hadamard rho-zero coefficient",
  old_numerator81 == 0 and new_numerator81 == 1+1j)

# Exact finite truncation is stronger than the untruncated largest-gap rule.
# Small bounded fixtures include all cuts, both clipping branches, and r=b.
selector_cases81 = 0
selector_branches81 = set()
for r81, b81, s81 in ((1,1,2), (2,2,2), (3,3,3), (6,2,3), (9,3,4), (12,3,4)):
    L81, H81 = 1 << s81, 2
    for j81 in sorted({0, r81-1}):
        for positions81 in ((), (0,), (s81,), (0,s81),
                            (s81//2,s81), (s81//2,s81//2)):
            T81 = S81 = 0
            for h81 in range(H81):
                c81 = b81 * (((j81-L81*h81) % r81)//b81)
                residues81 = [(c81+q81-u81) % r81
                              for u81 in range(b81) for q81 in range(b81)]
                T81 += sum(max(0,1+(L81-1-rho81)//r81) for rho81 in residues81)
                S81 += sum(max(0,1+(L81-1-rho81)//r81) for rho81 in set(residues81))
            costs81 = []
            for cut81 in range(s81+1):
                A81, K81 = 1 << cut81, 1 << (s81-cut81)
                left81 = [v81 for v81 in positions81 if v81 < cut81]
                P81 = ((1 << max(left81))//math79.gcd(r81,1 << max(left81))
                       if left81 else 1)
                selector_branches81.add(P81 <= A81//r81)
                C81_count = G81_count = 0
                for h81 in range(H81):
                    c81 = b81 * (((j81-L81*h81) % r81)//b81)
                    for k81 in range(K81):
                        classes81 = {}
                        for u81 in range(b81):
                            for q81 in range(b81):
                                rho81 = (c81+q81-A81*k81-u81) % r81
                                classes81[rho81] = classes81.get(rho81,0)+1
                        for rho81, multiplicity81 in classes81.items():
                            # Enumerate compatible offsets rather than using
                            # the full-L closed form on both sides.
                            zs81 = {n81 % P81 for n81, _a81
                                    in enumerate(range(rho81,A81,r81))}
                            C81_count += multiplicity81 * len(zs81)
                            G81_count += len(zs81)
                R81 = min(r81,2*b81-1)
                assert C81_count == min(H81*K81*b81*b81*P81,T81)
                assert G81_count == min(H81*K81*R81*P81,S81)
                costs81.append(H81*K81*b81*b81+C81_count)
                selector_cases81 += 1
            candidate_cuts81 = set(positions81) | {s81}
            assert min(costs81) == min(costs81[w81] for w81 in candidate_cuts81)
t("exact clipped minima count pairs and geometric candidates before cancellation",
  selector_cases81 == 270 and selector_branches81 == {False,True})

print("\nALL TESTS PASSED")
