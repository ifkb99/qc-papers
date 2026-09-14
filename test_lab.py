"""Engine tests for lab/ -- every extracted helper is pinned to a number
already logged in NOTES.md/CLAIMS.md before the extraction, so a refactor
that changes behaviour fails against the historical record, not against
itself. Part of the correctness gate (suite 6 of 9).
"""
from __future__ import annotations
import os

import numpy as np

os.environ["LAB_NO_CACHE"] = "1"          # engine tests never touch the cache

from lab import (Experiment, rank_kernel, structures, find_structures,
                 check_structure, coset_split, quadrant_counts, slice_fn,
                 order, v2_split, ord2, carmichael, bit_table, tile,
                 build_modexp, verify_correctness, WRAPS,
                 random_boolean, random_table, planted_structure,
                 fn_spectrum, fn_support, support, sparsity)
import lab.measure as measure


def t(name, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    assert ok, name


print("[1] modarith vs known values")
t("order(6,7)=2", order(6, 7) == 2)
t("order(3,7)=6", order(3, 7) == 6)
t("order(7,15)=4", order(7, 15) == 4)
t("v2_split(12)=(2,3)", v2_split(12) == (2, 3))
t("v2_split(16)=(4,1)", v2_split(16) == (4, 1))
t("ord2(9)=6", ord2(9) == 6)
t("ord2(1)=0", ord2(1) == 0)
t("carmichael: 15->4, 21->6, 51->16, 143->60",
  [carmichael(x) for x in (15, 21, 51, 143)] == [4, 6, 16, 60])
t("bit_table(7,3,0) = [1,1,0,0,0,1]",
  bit_table(7, 3, 0).tolist() == [1, 1, 0, 0, 0, 1])
t("tile length and content", tile(np.array([1, 0, 1]), 4).tolist()
  == [1, 0, 1] * 5 + [1])

print("[2] gf2: rank/kernel against experiment_linstruct logged values")
# modexp N=5 a=2 n_exp=1 (q=14): |supp|=3086, rank 13, w = {b3, anc12}
me = build_modexp(N=5, a=2, n_exp=1)
qc = me.build()
zs_modexp = support(qc, me.x[0])
t("N=5 modexp |supp| = 3086", zs_modexp.size == 3086)
rank, kern = rank_kernel(zs_modexp.tolist(), qc.n)
t("rank 13, defect 1", (rank, len(kern)) == (13, 1))
w0 = (1 << me.b[me.m - 1]) | (1 << me.anc)
t("kernel element = b_msb^anc", kern[0] == w0)
ev, od = coset_split(zs_modexp, w0)
t("support entirely in <z,w>=0", od == 0 and ev == 3086)

print("[3] gf2: affine finder on planted structures (C33 machinery)")
g, w, eps = planted_structure(12, rng=1, eps=1)
zs, st = find_structures(g, 12)
t("planted affine found", st == [(w, 1)])
g, w, eps = planted_structure(12, rng=2, eps=0)
zs, st = find_structures(g, 12)
t("planted linear (dead var) found", st == [(w, 0)])
g = random_boolean(12, rng=3)
zs, st = find_structures(g, 12)
t("random function: no structure (must-fail control)", st == [])

print("[4] gf2: function-level all-ones structures (C33 pinned)")
g = tile(bit_table(7, 3, 0), 12)
zs, st = find_structures(g, 12)
t("r=6 t=12: density exactly 1/2", zs.size == 2048)
t("r=6 t=12: w=all-ones AFFINE", st == [((1 << 12) - 1, 1)])
g13 = tile(bit_table(7, 3, 0), 13)
_, st13 = find_structures(g13, 13)
t("r=6 t=13: none (odd t)", st13 == [])
g3 = tile(bit_table(7, 2, 0), 12)
_, st3 = find_structures(g3, 12)
t("r=3 t=12: w=all-ones LINEAR", st3 == [((1 << 12) - 1, 0)])

print("[5] slices and quadrants against experiment_resid2 logged values")
me4 = build_modexp(N=5, a=2, n_exp=2, wraps=("msb_t0_anc",))
qc4 = me4.build()
zs4 = support(qc4, me4.x[0])
t("v4 |supp| = 23464", zs4.size == 23464)
q = quadrant_counts(zs4, me4.b[me4.m - 1], me4.anc)
t("v4 quadrant (msb=1,anc=0) EMPTY", q[(1, 0)] == 0)
t("v4 quadrants sum to |supp|", sum(q.values()) == 23464)
import walsh as _walsh
g4 = ((_walsh.classical_permutation(qc4) >> me4.x[0]) & 1).astype(np.int64)
s1 = slice_fn(g4, qc4.n, me4.t[0], 1)
w14 = 1 << (me4.b[me4.m - 1])          # msb < t0, position unchanged in slice
t("v4 t0=1 slice structure = msb alone (pinned)",
  check_structure(s1, w14, 0))
t("slice_fn shape", s1.size == 1 << (qc4.n - 1))

print("[6] variants: registry, correctness gate, C35 constancy")
t("known wraps registered",
  set(WRAPS) >= {"msb_t0_anc", "msb_t0_t1", "msb_t1_anc"})
t("v4 computes a^e mod N", verify_correctness(me4))
t("baseline computes a^e mod N",
  verify_correctness(build_modexp(N=5, a=2, n_exp=2)))
try:
    build_modexp(N=5, a=2, n_exp=1, wraps=("nope",))
    t("unknown wrap rejected", False)
except ValueError:
    t("unknown wrap rejected", True)
s2 = sparsity(build_modexp(N=7, a=6, n_exp=2, wraps=("msb_t0_anc",)).build(),
              build_modexp(N=7, a=6, n_exp=2).x[0])
s3 = sparsity(build_modexp(N=7, a=6, n_exp=3, wraps=("msb_t0_anc",)).build(),
              build_modexp(N=7, a=6, n_exp=3).x[0])
t("C35: v4 constancy 23488 at n_exp=2,3", s2 == 23488 and s3 == 23488)

print("[7] measure: cache round-trip (content-addressed)")
os.environ["LAB_NO_CACHE"] = "0"
key_path = measure.CACHE_DIR / f"supp_{measure._key(qc, me.x[0])}.npy"
if key_path.exists():
    key_path.unlink()
zs_a = support(qc, me.x[0])
t("cache file created", key_path.exists())
zs_b = support(qc, me.x[0])
t("cache round-trip identical", np.array_equal(zs_a, zs_b))
t("cached == uncached", np.array_equal(zs_a, zs_modexp))
os.environ["LAB_NO_CACHE"] = "1"

print("[8] nulls: trap defaults")
h = random_table(3, rng=0)
t("random_table rejects constants (r=3)",
  all(0 < int(random_table(3, rng=s).sum()) < 3 for s in range(20)))
t("fn_spectrum normalised (Parseval)",
  abs((fn_spectrum(random_boolean(10, rng=4)) ** 2).sum() - 1) < 1e-9)

print("[9] harness: protocol enforcement")
exp = Experiment("selftest", exit_on_fail=False)
exp.predict("P1", "trivially true")
exp.check("P1", True)
exp.predict("LATE", "declared after measuring")     # must warn
exp.check("LATE", True)
ok = exp.finish()
t("passing experiment reports ok", ok)
t("derive-then-test violation warned",
  any("AFTER measurement" in w for w in exp.warnings))
t("missing must-fail control warned",
  any("must-fail" in w for w in exp.warnings))
exp2 = Experiment("selftest2", exit_on_fail=False)
exp2.must_fail("C1", "control")
exp2.fail_check("C1", False)                        # control did NOT fail
t("vacuous control detected", not exp2.finish())

print("[10] structured experiment reports retain failures and declarations")
import json
from pathlib import Path
import tempfile
with tempfile.TemporaryDirectory(prefix="qsim-report-test-") as tmp:
    report = Path(tmp) / "failed.json"
    failed = Experiment("report-test", exit_on_fail=True)
    failed.predict("P1", "deliberately false")
    failed.must_fail("C1", "known control")
    failed.check("P1", False, "must be recorded before exit")
    failed.fail_check("C1", True)
    try:
        failed.finish(report_path=report, rows=[dict(value=3)], metadata=dict(seed=1))
        t("failing report exits nonzero", False)
    except SystemExit as exc:
        t("failing report exits nonzero", exc.code == 1)
    saved = json.loads(report.read_text())
    t("failure report exists and is complete", not saved["ok"]
      and saved["checks"][0]["ok"] is False and saved["rows"] == [dict(value=3)]
      and saved["metadata"]["seed"] == 1 and "P1" in saved["predictions"])

print("[11] tensor-cut diagnostics and residue automaton")
from lab.tensors import cut_matrix, cut_spectrum, numerical_rank, rank_profile, residue_automaton
bell = np.array([1., 0., 0., 1.]) / np.sqrt(2)
t("Bell vector has cut rank two", numerical_rank(cut_spectrum(bell, [0])) == 2)
t("product sign tensor has rank one", rank_profile(np.array([1., -1., 1., -1.]))[0]["rank"] == 1)
t("complex tensor supported", numerical_rank(cut_spectrum(bell * 1j, [1])) == 2)
t("little-endian cut indices", cut_matrix(np.arange(8), [2, 0]).tolist()
  == [[0, 2], [4, 6], [1, 3], [5, 7]])
t("zero tensor rank zero", numerical_rank(cut_spectrum(np.zeros(8), [1])) == 0)
t("residue automaton MSB convention", residue_automaton(np.arange(6), [1, 0, 1, 1]) == 5)
for vector, left in [(np.arange(6), [0]), (np.arange(8), [0, 0]), (np.arange(8), [3])]:
    try:
        cut_matrix(vector, left)
        t("invalid tensor cut rejected", False)
    except ValueError:
        t("invalid tensor cut rejected", True)

print("[12] conditional measurement instruments preserve interference")
from lab.semiclassical import conditional_children, distribution, sample, validate_inputs
identity = np.arange(2)
flip = np.array([1, 0])
zero = np.array([1., 0.], dtype=complex)
plus, minus = conditional_children(zero, (identity, flip), 0.)
t("conditional X branches retain opposite coherent signs",
  np.array_equal(plus, [.5, .5]) and np.array_equal(minus, [.5, -.5]))
t("instrument probabilities sum to one", abs(np.vdot(plus, plus)+np.vdot(minus, minus)-1) < 1e-12)
out, _ = distribution([(identity, identity)] * 3, zero)
t("identity arithmetic produces only Fourier output zero", np.allclose(out, [1, 0, 0, 0, 0, 0, 0, 0]))
path = sample([(identity, identity)] * 12, zero, np.random.default_rng(1))
t("sampling follows one normalized deterministic path", path["output"] == 0
  and path["path_probability"] == 1 and path["steps"] == 12)
try:
    validate_inputs([(identity, np.array([0, 0]))], zero)
    t("non-permutation branch rejected", False)
except ValueError:
    t("non-permutation branch rejected", True)
try:
    distribution([(identity, identity)] * 11, zero)
    t("unintended exponential enumeration rejected", False)
except ValueError:
    t("unintended exponential enumeration rejected", True)

print("\n[13] reachable-state conditional instruments")
from lab.reachable import sparse_children, sparse_path, compiled_pairs, action_stats
identity_action = lambda ids: ids.copy()
flip_action = lambda ids: ids ^ 1
labels, children = sparse_children(np.array([0]), np.array([1.+0j]),
                                   (identity_action, flip_action), 0.)
t("sparse instrument preserves signs and interference", labels.tolist() == [0, 1]
  and np.array_equal(children[0], [.5, .5]) and np.array_equal(children[1], [.5, -.5]))
sp = sparse_path([(identity_action, identity_action)] * 12, {0: 1.}, output=0)
t("sparse deterministic identity needs one amplitude", sp["path_probability"] == 1
  and sp["peak_stored_amplitudes"] == 1 and len(sp["profile"]) == 12)
sp = sparse_path([(identity_action, identity_action)], {0: 1.}, output=1)
t("forced impossible branch does not divide by zero", sp["path_probability"] == 0
  and sp["zero_probability_step"] == 0)
me, pairs = compiled_pairs(7, 3, 4)
t("compiled setup performs no work-space enumeration", action_stats(pairs)["cached_transitions"] == 0)
small_ids = np.arange(8)
dense_pairs = [(small_ids, np.where(small_ids < 7, small_ids*pow(3, 1 << i, 7) % 7, small_ids))
               for i in range(4)]
dense_probs, _ = distribution(dense_pairs, np.eye(8, dtype=complex)[1])
sparse_probs = [sparse_path(pairs, {1 << me.x[0]: 1.}, output=y)["path_probability"]
                for y in range(16)]
t("compiled sparse full distribution agrees with dense ideal work maps",
  np.max(np.abs(sparse_probs-dense_probs)) < 1e-12)
try:
    sparse_children(np.array([0, 1]), np.array([1., 0.]),
                    (identity_action, lambda ids: np.zeros_like(ids)), 0.)
    t("noninjective sparse branch rejected", False)
except ValueError:
    t("noninjective sparse branch rejected", True)
try:
    sparse_path([], {0: 2.}, output=0)
    t("unnormalized sparse input rejected", False)
except ValueError:
    t("unnormalized sparse input rejected", True)

from lab.semiclassical import eigenphase_path
mixture = np.array([sum(eigenphase_path(6, 4, k, output=y)
                        ["conditional_path_probability"] for k in range(6))/6
                    for y in range(16)])
t("latent-eigenphase mixture agrees with independent work-state instruments",
  np.max(np.abs(mixture-dense_probs)) < 1e-12)
phase_sample = eigenphase_path(4, 8, 1, rng=np.random.default_rng(3))
t("exact representable eigenphase has deterministic Fourier output",
  phase_sample["output"] == 64 and abs(phase_sample["conditional_path_probability"]-1) < 1e-12)
try:
    eigenphase_path(6, 4, 6, output=0)
    t("invalid eigenphase index rejected", False)
except ValueError:
    t("invalid eigenphase index rejected", True)

print("[14] phase-adjusted scalar sampling and bounded single-defect effects")
from lab.spectral import single_defect_effects, split_phase_filters, coherence_response
from circuits import Circuit
import statevec
# Independent control-register circuit: eigenphase kickback plus arbitrary
# product phases, with no arithmetic or spectral-effects helper involved.
phase_offsets = np.array([.2, -.7, .4])
for eigen_k in range(3):
    qc_phase = Circuit(3)
    for bit in range(3):
        qc_phase.h(bit).rz(bit, 2*np.pi*eigen_k*(1 << bit)/3 + phase_offsets[bit])
    qc_phase.qft([0, 1, 2], inverse=True)
    phase_reference = np.abs(statevec.run(qc_phase))**2
    phase_prob = np.array([eigenphase_path(3, 3, eigen_k, output=y,
                                          input_phases=phase_offsets)
                           ["conditional_path_probability"] for y in range(8)])
    t("arbitrary product input phases agree with coherent inverse QFT",
      np.max(np.abs(phase_reference-phase_prob)) < 1e-12)
    sampled_phase = eigenphase_path(3, 3, eigen_k, rng=np.random.default_rng(9),
                                    input_phases=phase_offsets)
    t("sampled conditional probability equals forced-path evaluation",
      abs(sampled_phase["conditional_path_probability"]
          - phase_prob[sampled_phase["output"]]) < 1e-12)

empty = eigenphase_path(1, 0, 0, output=0, input_phases=[])
t("zero-width scalar path remains normalized", empty["conditional_path_probability"] == 1)
for bad_phases in ([0.], [0., float("nan")], [[0., 0.]], [0., 1j]):
    try:
        eigenphase_path(2, 2, 0, output=0, input_phases=bad_phases)
        t("invalid input phase vector rejected", False)
    except ValueError:
        t("invalid input phase vector rejected", True)

# Non-diagonal unitary on the eigenbasis; compare the closed-form factorization
# against its defining sum over every exponent, retaining the stated ordering.
V_test = np.array([[1., -1j], [-1j, 1.]]) / np.sqrt(2)
effects_test = single_defect_effects(2, 3, 1, V_test)
eigenvalues = np.array([1., -1.])
direct_effects = []
for output_y in range(8):
    K_direct = np.zeros((2, 2), complex)
    for exponent in range(8):
        low = exponent % 2
        high = exponent-low
        K_direct += (np.exp(-2j*np.pi*output_y*exponent/8)/8
                     * eigenvalues[:, None]**high * V_test * eigenvalues[None, :]**low)
    direct_effects.append(K_direct.conj().T @ K_direct)
t("single-defect filters preserve ascending time order",
  np.max(np.abs(effects_test-direct_effects)) < 1e-12)
t("single-defect effects sum to identity",
  np.max(np.abs(effects_test.sum(axis=0)-np.eye(2))) < 1e-12)
t("defect after all arithmetic cannot change output effects",
  np.max(np.abs(single_defect_effects(2, 3, 3, V_test)
                - single_defect_effects(2, 3, 3, np.eye(2)))) < 1e-12)
for split in range(4):
    no_defect = single_defect_effects(3, 3, split, np.eye(3))
    t("identity defect has zero coherence-response rank at every split",
      coherence_response(no_defect)["rank"] == 0)
# Effects of a Y-basis measurement have imaginary off-diagonals: this also
# checks that the diagnostic does not silently ignore their real coordinates.
y_effects = np.array([[[.5, -.5j], [.5j, .5]], [[.5, .5j], [-.5j, .5]]])
t("imaginary coherence is detected", coherence_response(y_effects)["rank"] == 1
  and coherence_response(y_effects)["detectable_pairs"] == 1)
one_effect = np.ones((1, 1, 1))
t("one-dimensional orbit has no off-diagonal response",
  coherence_response(one_effect)["rank"] == 0)
for params in ((2, 13, 1), (2, 3, -1), (2, 3, 4), (0, 3, 1)):
    try:
        split_phase_filters(*params)
        t("invalid spectral dimensions rejected", False)
    except ValueError:
        t("invalid spectral dimensions rejected", True)
try:
    single_defect_effects(2, 3, 1, np.zeros((2, 2)))
    t("nonunitary orbit defect rejected", False)
except ValueError:
    t("nonunitary orbit defect rejected", True)
try:
    split_phase_filters(32, 10, 2)
    t("dense spectral allocation cap enforced", False)
except ValueError:
    t("dense spectral allocation cap enforced", True)
try:
    split_phase_filters(np.int64(1 << 62), 2, 1)
    t("allocation budget cannot overflow int64", False)
except ValueError:
    t("allocation budget cannot overflow int64", True)

print("[15] final-eigenphase rejection and a coherent early prefix")
from lab.prefix import SparseOrbitPrefix
identity_prefix = SparseOrbitPrefix(3, 2,
                                    lambda l: (np.array([l % 3]), np.array([1.])), 1)
for output in range(16):
    for k in range(3):
        joint_prefix = identity_prefix.forced_joint(4, k, output)
        expected_joint = eigenphase_path(3, 4, k, output=output)["conditional_path_probability"]/3
        t("identity prefix equals ideal conditional eigenphase sampling",
          abs(joint_prefix["joint_latent_output_probability"]-expected_joint) < 1e-12)
prefix_sample = identity_prefix.sample(32, np.random.default_rng(5))
t("unit sparse envelope always accepts first attempt",
  prefix_sample["rejection_attempts"] == 1 and prefix_sample["early_amplitudes"] == 4)
t("streamed prefix owns no orbit table", identity_prefix.stats()["orbit_table_entries"] == 0
  and identity_prefix.stats()["phase_row_payload_bytes"] == 64)
try:
    identity_prefix.forced_joint(4, None, 0)
    t("missing phase label has a validation error", False)
except ValueError:
    t("missing phase label has a validation error", True)
cancelled = SparseOrbitPrefix(2, 0,
                             lambda l: (np.array([0, 1]), np.array([1., -1.])/np.sqrt(2)), 2)
zero_joint = cancelled.forced_joint(2, 0, 0)
t("impossible latent phase avoids division by zero",
  zero_joint["phase_probability"] == 0
  and zero_joint["conditional_path_probability"] is None
  and zero_joint["joint_latent_output_probability"] == 0)
class RejectedPhaseRNG:
    def integers(self, high):
        return 0
    def random(self):
        return 0.
try:
    cancelled.draw_phase(RejectedPhaseRNG(), max_attempts=2)
    t("rejection-cap exhaustion raises instead of returning a biased sample", False)
except RuntimeError:
    t("rejection-cap exhaustion raises instead of returning a biased sample", True)
for bad_column in (lambda l: ([3], [1.]), lambda l: ([0, 0], [.6, .8]),
                   lambda l: ([0], [2.]), lambda l: ([0.5], [1.]),
                   lambda l: ([0, 1], [1., 1j])):
    try:
        SparseOrbitPrefix(3, 1, bad_column, 2)
        t("malformed sparse orbit column rejected", False)
    except ValueError:
        t("malformed sparse orbit column rejected", True)
for bad_width in (1, 64):
    before_calls = identity_prefix.stats()["column_calls"]
    try:
        identity_prefix.sample(bad_width, np.random.default_rng(1))
        t("invalid width rejected before rejection work", False)
    except ValueError:
        t("invalid width rejected before rejection work",
          identity_prefix.stats()["column_calls"] == before_calls)

print("[16] finite Fourier marginals and localized defects without prefix arrays")
from lab.fourier_sampling import interval_path, interval_prefix_probability, geometric_sum
from lab.localized import LocalizedOrbitDefect
for width16, count16 in ((0, 1), (3, 3), (5, 21)):
    state16 = np.zeros(1 << width16, complex)
    state16[:count16] = np.exp(-2j*np.pi*2*np.arange(count16)/7)/np.sqrt(count16)
    expected16 = np.abs(np.fft.fft(state16))**2/(1 << width16)
    got16 = np.array([interval_path(width16, count16, 2, 7, output=y)["path_probability"]
                      for y in range(1 << width16)])
    t("finite Fourier bit paths equal independent FFT", np.max(np.abs(got16-expected16)) < 1e-12)
t("geometric sum has exact root-of-unity zero", geometric_sum(8, 1, 8) == 0)
t("singleton low-bit marginal is uniform at width 63",
  interval_prefix_probability(63, 1, 7, 13, 63, (1 << 63)-1) == 1/(1 << 63))
empty16 = LocalizedOrbitDefect(3, (), np.empty((0, 0)))
sample16 = empty16.sample(63, 62, np.random.default_rng(16))
t("identity localized sampler accepts both proposals immediately",
  sample16["phase_attempts"] == sample16["early_attempts"] == 1)
t("no prefix or output array in wide localized sampler",
  empty16.stats()["early_vector_entries"] == empty16.stats()["output_table_entries"] == 0)
cancel16 = LocalizedOrbitDefect(3, (0, 1), np.array([[1., 1.], [-1., 1.]])/np.sqrt(2))
zero16 = cancel16.forced_joint(3, 0, 0, 0)
t("localized impossible phase is explicitly undefined conditionally",
  zero16["phase_probability"] == 0 and zero16["conditional_path_probability"] is None
  and zero16["joint_latent_output_probability"] == 0)
try:
    cancel16.sample(3, 0, RejectedPhaseRNG(), max_attempts=2)
    t("localized phase cap raises, no biased fallback", False)
except RuntimeError:
    t("localized phase cap raises, no biased fallback", True)
for bad16 in (lambda: LocalizedOrbitDefect(100, range(65), np.empty((0, 0))),
              lambda: LocalizedOrbitDefect(3, (0, 0), np.eye(2)),
              lambda: LocalizedOrbitDefect(3, (0,), np.zeros((1, 1))),
              lambda: empty16.forced_joint(3, 1, None, 0),
              lambda: empty16.forced_joint(3, 1, 0, None),
              lambda: empty16.sample(64, 1, None),
              lambda: interval_path(64, 1, 0, 1, output=0)):
    try:
        bad16()
        t("localized dimensions/inputs rejected before expensive work", False)
    except ValueError:
        t("localized dimensions/inputs rejected before expensive work", True)

print("[17] sequential forward/backward instruments preserve noncommuting time order")
from lab.semiclassical import sequential_path
from lab.periodic import PeriodicOrbitCircuit
theta17, phi17 = .7, .4
X17 = np.array([[0., 1.], [1., 0.]])
Z17 = np.diag(np.exp(np.array([-.5j, .5j])*theta17))
R17 = np.cos(phi17/2)*np.eye(2)-1j*np.sin(phi17/2)*X17
pairs17 = [(Z17, Z17 @ X17), (R17, R17 @ X17)]
qc17 = Circuit(3).h(1).h(2).cnot(1, 0).rz(0, theta17).cnot(2, 0).rx(0, phi17)
qc17.qft([1, 2], inverse=True)
psi17 = np.zeros(8, complex)
psi17[0] = 1.
ref17 = np.sum(np.abs(statevec.run(qc17, psi17).reshape(4, 2))**2, axis=1)
got17 = np.array([sequential_path(pairs17, np.array([1., 0.]), output=y)
                  ["conditional_path_probability"] for y in range(4)])
t("square-root forward/backward sampler equals existing coherent circuit", np.max(np.abs(got17-ref17)) < 1e-12)
wrong17 = np.array([sequential_path(list(reversed(pairs17)), np.array([1., 0.]), output=y)
                    ["conditional_path_probability"] for y in range(4)])
t("reversing noncommuting sequential blocks changes the output", np.max(np.abs(wrong17-ref17)) > .01)
zero17 = sequential_path([(np.eye(1), np.eye(1))]*3, np.ones(1), output=1)
t("zero-weight forced sequential branch terminates without division", zero17["conditional_path_probability"] == 0)
periodic17 = PeriodicOrbitCircuit(6, 1, 5, {})
for k17 in range(6):
    for y17 in (0, 1, 9, 31):
        t("one-dimensional conserved sectors reduce to ideal scalar sampling",
          abs(periodic17.forced_joint(k17, y17)["conditional_path_probability"]
              - eigenphase_path(6, 5, k17, output=y17)["conditional_path_probability"]) < 1e-12)
wide17 = PeriodicOrbitCircuit(3_000_000_021, 3, 63, {16: np.eye(3), 32: np.eye(3)})
sample17 = wide17.sample(np.random.default_rng(17))
t("wide coarse-sector path uses a three-dimensional work register",
  sample17["matrix_dimension"] == 3 and sample17["forward_factor_payload_bytes"] == 64*9*16)
for bad17 in (lambda: sequential_path(pairs17, np.array([1., 0.]), output=0, max_payload_bytes=1),
              lambda: sequential_path([(np.eye(2), np.zeros((2, 2)))], np.array([1., 0.]), output=0),
              lambda: PeriodicOrbitCircuit(7, 3, 4, {}),
              lambda: PeriodicOrbitCircuit(6, 3, 4, {5: np.eye(3)}),
              lambda: PeriodicOrbitCircuit(6, 3, 4, {1: np.zeros((3, 3))})):
    try:
        bad17()
        t("invalid sequential/periodic inputs rejected before propagation", False)
    except ValueError:
        t("invalid sequential/periodic inputs rejected before propagation", True)

print("[18] exact integer light-cone covers certify localized-kick omission")
from lab.periodic import lightcone_cover
cover18 = lightcone_cover(9, 4, (0, 8), 1)
t("wrapped intervals merge and count the truncated orbit traversal exactly",
  cover18["intervals"] == [(0, 2), (7, 9)]
  and cover18["covered_prefix_labels"] == 6
  and cover18["hit_probability_denominator"] == 16)
t("empty support has exactly zero squared TV bound",
  lightcone_cover(9, 63, (), 100)["tv_squared_numerator"] == 0)
t("a radius spanning the whole orbit gives full coverage",
  lightcone_cover(9, 4, (4,), 4)["covered_prefix_labels"] == 16)
wide18 = PeriodicOrbitCircuit(3_000_000_021, 3, 63,
                            {0: np.eye(3), 16: np.eye(3), 48: np.eye(3)})
bound18 = wide18.localized_kick_bound(32, (0, 1))
t("only preceding background mixers determine the cover radius",
  bound18["pre_kick_mixer_count"] == 2 and bound18["radius"] == 4)
t("large-order certificate is nontrivial and uses exact rational integers",
  0 < bound18["tv_squared_numerator"] < bound18["tv_squared_denominator"]
  and bound18["orbit_table_entries"] == bound18["prefix_table_entries"] == 0)
for bad18 in (lambda: lightcone_cover(0, 3, (0,), 0),
              lambda: lightcone_cover(6, 64, (0,), 1),
              lambda: lightcone_cover(6, 3, (0, 0), 1),
              lambda: lightcone_cover(6, 3, (6,), 1),
              lambda: lightcone_cover(6, 3, (0,), -1),
              lambda: lightcone_cover(6, 3, [0]*65, 1),
              lambda: lightcone_cover(6, 3, np.zeros((1, 1)), 1),
              lambda: wide18.localized_kick_bound(64, (0,))):
    try:
        bad18()
        t("invalid cover inputs rejected before interval construction", False)
    except ValueError:
        t("invalid cover inputs rejected before interval construction", True)

print("[19] deterministic sector routes preserve a finite-work sampler")
from lab.periodic import RoutedOrbitCircuit
qc19 = Circuit(4).h(2).h(3).toffoli(2, 0, 1).cnot(2, 0)
qc19.rz(0, theta17).rz(1, np.pi).cnot(3, 1).rx(0, phi17)
qc19.qft([2, 3], inverse=True)
ref19 = np.sum(np.abs(statevec.run(qc19).reshape(4, 4))**2, axis=1)
route19 = RoutedOrbitCircuit(4, 2, 2, {1: Z17, 2: R17}, {1: 1})
got19 = np.array([sum(route19.forced_joint(a, y)["joint_latent_output_probability"]
                     for a in range(2)) for y in range(4)])
t("routed mixture equals an independent controlled-cycle physical circuit",
  np.max(np.abs(got19-ref19)) < 1e-12)
same19 = RoutedOrbitCircuit(4, 2, 2, {1: Z17, 2: R17}, {0: 2, 1: -1, 2: 4})
t("initial/final and negative charges reduce modulo sector count",
  max(abs(route19.forced_joint(a, y)["joint_latent_output_probability"]
          - same19.forced_joint(a, y)["joint_latent_output_probability"])
      for a in range(2) for y in range(4)) < 1e-12)
empty19 = RoutedOrbitCircuit(6, 2, 0, {0: R17}, {0: -1})
draw19 = empty19.sample(np.random.default_rng(19))
t("zero-control path still applies the initial sector route",
  draw19["output"] == 0 and draw19["conditional_path_probability"] == 1
  and draw19["final_coarse_eigenphase"] == (draw19["coarse_eigenphase"]+1) % 3)
wide19 = RoutedOrbitCircuit(3_000_000_021, 3, 63,
                           {0: np.eye(3), 16: np.eye(3), 48: np.eye(3)},
                           {0: -1, 16: 1, 32: 2, 48: 1_000_000_008})
draw19 = wide19.sample(np.random.default_rng(1919))
t("wide route uses bounded work with the predicted final sector",
  draw19["matrix_dimension"] == 3
  and draw19["final_coarse_eigenphase"] == (draw19["coarse_eigenphase"]-3) % wide19.sectors
  and wide19.stats()["orbit_table_entries"] == 0)
bound19 = wide19.localized_kick_bound(32, (0, 1))
t("diagonal routing phases add zero displacement to the omission bound",
  bound19["covered_prefix_labels"] == bound18["covered_prefix_labels"]
  and bound19["diagonal_routing_phases_add_zero_displacement"])
for bad19 in (lambda: RoutedOrbitCircuit(6, 2, 2, {}, {1: .5}),
              lambda: RoutedOrbitCircuit(6, 2, 2, {}, {3: 1}),
              lambda: RoutedOrbitCircuit(6, 2, 2, {}, {1.0: 1}),
              lambda: RoutedOrbitCircuit(6, 2, 64, {}, {}),
              lambda: route19.forced_joint(2, 0),
              lambda: route19.forced_joint(0, 4),
              lambda: route19.sector_path(0),
              lambda: route19.sector_path(0, rng=np.random.default_rng(1), output=0)):
    try:
        bad19()
        t("invalid routed dimensions/charges/paths rejected", False)
    except ValueError:
        t("invalid routed dimensions/charges/paths rejected", True)

print("[20] coherent route histories retain prefix amplitudes and final interference")
from lab.coherent_routes import CoherentReflectionCircuit
theta20 = .9
coherent20 = CoherentReflectionCircuit(4, 2, 2, {1: Z17, 2: R17}, {1: (1, theta20)})
qc20 = Circuit(4).h(2).h(3)
prefix_error20 = 0.
for stage20 in ("arithmetic", "background", "reflection"):
    if stage20 == "arithmetic":
        qc20.toffoli(2, 0, 1).cnot(2, 0)
    elif stage20 == "background":
        qc20.rz(0, theta17)
    else:
        qc20.rz(1, theta20)
    psi20 = statevec.run(qc20).reshape(4, 2, 2)
    expected_all20, actual_all20 = [], []
    for a20 in range(2):
        for x20 in range(4):
            expected20 = (psi20[x20, 0]+(-1)**a20*psi20[x20, 1])/np.sqrt(2)
            actual20 = coherent20.prefix_vector(a20, x20, 1, boundary=stage20)/np.sqrt(2)
            expected_all20.append(expected20)
            actual_all20.append(actual20)
    # Circuit's H/CNOT decomposition is explicitly up to global phase.
    # Align ONCE for the entire state, never separately by sector or history.
    expected_all20, actual_all20 = np.array(expected_all20), np.array(actual_all20)
    overlap20 = np.vdot(actual_all20, expected_all20)
    common_phase20 = overlap20/abs(overlap20)
    prefix_error20 = max(prefix_error20,
                         np.max(np.abs(expected_all20-common_phase20*actual_all20)))
t("all three work-gate boundaries match an independent physical prefix state",
  prefix_error20 < 1e-12)
qc20.cnot(3, 1).rx(0, phi17).qft([2, 3], inverse=True)
psi20 = statevec.run(qc20).reshape(4, 2, 2)
joint20 = np.array([[np.sum(np.abs((psi20[y, 0]+(-1)**a*psi20[y, 1])/np.sqrt(2))**2)
                     for y in range(4)] for a in range(2)])
got20 = np.array([[coherent20.joint_probability(a, y) for y in range(4)] for a in range(2)])
t("complete final-sector/output joint law equals the independently compiled circuit",
  np.max(np.abs(joint20-got20)) < 1e-12 and abs(got20.sum()-1) < 1e-12)
accepted20 = np.array([[coherent20.rejection_weights(a, y)["proposal_joint_probability"]
                         * coherent20.rejection_weights(a, y)["acceptance_probability"]
                         for y in range(4)] for a in range(2)])
t("coherent rejection has the expected normalization and accepted law",
  abs(accepted20.sum()*coherent20.coefficient_l1**2-1) < 1e-12
  and np.max(np.abs(accepted20/accepted20.sum()-got20)) < 1e-12)
wide20 = CoherentReflectionCircuit(2_000_000_014, 2, 63,
                                  {16: R17, 48: Z17}, {17: (0, .6), 49: (1, .4)})
draw20 = wide20.sample(np.random.default_rng(20))
t("wide coherent sample uses prefix amplitudes without a rejection loop or sector table",
  0 <= draw20["output"] < 1 << 63 and draw20["history_count_upper_bound"] == 4
  and draw20["rejection_proposals"] == draw20["sector_table_entries"] == 0)
empty20 = CoherentReflectionCircuit(6, 2, 0, {0: R17}, {0: (1, .7)})
t("zero-control coherent input rotations give a normalized final-sector law",
  abs(sum(empty20.joint_probability(a, 0) for a in range(3))-1) < 1e-12
  and empty20.sample(np.random.default_rng(201))["output"] == 0)
for bad20 in (lambda: CoherentReflectionCircuit(6, 2, 2, {}, {1: (.5, .2)}),
              lambda: CoherentReflectionCircuit(6, 2, 2, {}, {1: (1, np.nan)}),
              lambda: CoherentReflectionCircuit(6, 2, 2, {}, {1: (1, 1j)}),
              lambda: CoherentReflectionCircuit(6, 2, 2, {}, {3: (1, .2)}),
              lambda: CoherentReflectionCircuit(6, 2, 63, {}, {s: (1, .2) for s in range(9)}),
              lambda: coherent20.prefix_vector(0, 0, 1, measured=1),
              lambda: coherent20.prefix_vector(0, 0, 2, measured=1, output=2),
              lambda: coherent20.prefix_vector(0, 0, 2, boundary="unknown"),
              lambda: coherent20.joint_probability(2, 0),
              lambda: coherent20.sample_rejection(None, max_proposals=0)):
    try:
        bad20()
        t("invalid coherent routes and prefix requests rejected before contraction", False)
    except ValueError:
        t("invalid coherent routes and prefix requests rejected before contraction", True)

print("[21] exact-rational sampling budgets state their oracle contract")
from fractions import Fraction
from lab.sampling_error import prefix_error_budget, plan_prefix_accuracy, rejection_error_budget
budget21 = prefix_error_budget(4, 2, [Fraction(1, 1024)]*8, random_bits=16)
t("scaled amplitude budget includes real/imag, global dimension and categorical bits",
  budget21["scaled_coordinate_to_state_l2_factor"] == 8
  and budget21["total_tv_upper_bound"] == Fraction(1, 8)+Fraction(8, 1 << 16)
  and budget21["requires_uniform_oracle_error_bound"]
  and not budget21["certifies_existing_float_sampler"])
plan21 = plan_prefix_accuracy(63, 2, 135, Fraction(1, 1_000_000))
t("wide accuracy plan uses integer arithmetic without claiming mantissa sufficiency",
  plan21["total_tv_upper_bound"] <= Fraction(1, 1_000_000)
  and plan21["absolute_accuracy_bits"] < 70
  and plan21["machine_mantissa_bits"] is None)
t("empty stochastic circuit has zero conditional budget",
  prefix_error_budget(0, 1, [], random_bits=1)["total_tv_upper_bound"] == 0)
t("rejection budget charges amplification by inverse acceptance mass",
  rejection_error_budget(4, Fraction(1, 100)) == Fraction(1, 25))
for bad21 in (lambda: prefix_error_budget(64, 2, []),
              lambda: prefix_error_budget(4, 2, [.01]),
              lambda: prefix_error_budget(4, 2, [-1]),
              lambda: prefix_error_budget(4, 2, [], random_bits=0),
              lambda: plan_prefix_accuracy(4, 2, 0, Fraction(1, 10)),
              lambda: plan_prefix_accuracy(4, 2, 1, 1),
              lambda: rejection_error_budget(Fraction(1, 2), 0)):
    try:
        bad21()
        t("invalid or inexact error-budget inputs rejected", False)
    except ValueError:
        t("invalid or inexact error-budget inputs rejected", True)

print("[22] verified exact-input prefix primitives (optional Arb backend)")
from lab.verified_prefix import (VerifiedReflectionCircuit, binary_fraction,
                                 integer_cdf_index, integer_cdf_counts,
                                 round_dyadic, uniform_integer)
from random import Random
t("signed ties-to-even rounding never passes through floats",
  [round_dyadic(Fraction(v, 8), 2) for v in (-3, -1, 1, 3)] == [-2, 0, 0, 2])
t("integer CDF includes zero weights and exact endpoint ties",
  [integer_cdf_index([1, 0, 3], w, 2) for w in range(4)] == [0, 2, 2, 2]
  and integer_cdf_counts([1, 0, 3], 2) == (1, 0, 3)
  and integer_cdf_counts([0, 0, 0], 2) == (2, 1, 1))
class ScriptedBits22:
    def __init__(self):
        self.words = iter((3, 2))
        self.calls = 0
    def getrandbits(self, bits):
        assert bits == 2
        self.calls += 1
        return next(self.words)
bits22 = ScriptedBits22()
t("initial integer rejects out-of-range words instead of modulo bias",
  uniform_integer(3, bits22) == 2 and bits22.calls == 2
  and uniform_integer(1, bits22) == 0 and bits22.calls == 2)
verified22 = VerifiedReflectionCircuit(10, 4,
    {1: ("x", Fraction(1, 7)), 3: ("z", Fraction(1, 5))},
    {2: (0, Fraction(1, 5)), 3: (1, Fraction(1, 5))})
for bad22 in (lambda: VerifiedReflectionCircuit(9, 4, {}, {}),
              lambda: VerifiedReflectionCircuit(10, 64, {}, {}),
              lambda: VerifiedReflectionCircuit(10, 4, {1: ("x", .1)}, {}),
              lambda: VerifiedReflectionCircuit(10, 4, {1: ("y", 1)}, {}),
              lambda: VerifiedReflectionCircuit(10, 4, {}, {1: (.1, 1)}),
              lambda: VerifiedReflectionCircuit(10, 4, {1: ("x", 1 << 2048)}, {}),
              lambda: verified22.prefix_dyadic(0, 0, 2, accuracy_bits=4, measured=1),
              lambda: verified22.accuracy_plan(.001),
              lambda: integer_cdf_index([1, -1], 0, 4),
              lambda: integer_cdf_index([1], 16, 4),
              lambda: round_dyadic(.1, 4)):
    try:
        bad22()
        t("invalid or inexact verified-input requests rejected", False)
    except ValueError:
        t("invalid or inexact verified-input requests rejected", True)
try:
    from flint import acb, arb, ctx, fmpq
except ImportError:
    print("  SKIP Arb oracle tests: install optional python-flint==0.9.0 to exercise these")
else:
    old22 = ctx.prec
    item22 = verified22.prefix_dyadic(1, 3, 4, accuracy_bits=40)
    with ctx.workprec(192):
        ball22 = verified22.prefix_enclosure(1, 3, 4, working_precision=192)
        errors22 = []
        for coordinates22, z22 in zip(item22.coordinates, ball22):
            for integer22, scalar22 in zip(coordinates22, (z22.real, z22.imag)):
                d22 = arb(fmpq(integer22, 1 << 40)) - scalar22
                errors22.append(binary_fraction(abs(d22).upper()))
    t("returned dyadics obey coordinate bound including midpoint rounding",
      max(errors22) <= item22.coordinate_error_bound
      and item22.max_ball_radius <= item22.coordinate_error_bound/2 and ctx.prec == old22)
    empty22 = VerifiedReflectionCircuit(6, 0, {}, {})
    t("zero-width and exact-zero coordinates terminate without relative error",
      empty22.prefix_dyadic(0, 0, 0, accuracy_bits=80).coordinates == ((1 << 80, 0), (0, 0))
      and empty22.sample(Random(22))["total_tv_upper_bound"] == 0)
    class InflatedEnclosure22(VerifiedReflectionCircuit):
        def _enclosure(self, *args):
            vector = super()._enclosure(*args)
            # Adding a zero-containing ball is a valid deliberate loss of
            # accuracy, forcing the real refinement branch without a new API.
            inflation = acb(arb(0, (1, -ctx.prec//2)), arb(0, (1, -ctx.prec//2)))
            return tuple(z + inflation for z in vector)
    retry22 = InflatedEnclosure22(2, 0, {}, {})
    retry_item22 = retry22.prefix_dyadic(0, 0, 0, accuracy_bits=40)
    t("refinement retries a valid loose enclosure and restores working context",
      retry_item22.refinements > 0 and retry_item22.working_precision > 64
      and ctx.prec == old22)
    try:
        retry22.prefix_dyadic(0, 0, 0, accuracy_bits=40, max_precision=64)
        t("diagnostic precision exhaustion raises instead of fallback sampling", False)
    except ArithmeticError:
        t("diagnostic precision exhaustion raises instead of fallback sampling", ctx.prec == old22)
    wide22 = VerifiedReflectionCircuit((1 << 61)-2, 63,
        {1: ("x", Fraction(1, 7)), 31: ("z", Fraction(1, 5)), 63: ("x", Fraction(1, 7))},
        {2: (0, Fraction(1, 5)), 32: (1, Fraction(1, 5)), 63: (2, Fraction(1, 5))})
    sample22 = wide22.sample(Random(624), target_tv=Fraction(1, 1_000_000))
    t("wide verified path charges all updates and finite bits without orbit tables",
      sample22["total_tv_upper_bound"] <= Fraction(1, 1_000_000)
      and sample22["prefix_vector_evaluations"] <= 135
      and sample22["requires_independent_unbiased_bits"]
      and not sample22["certifies_existing_float_sampler"]
      and sample22["precision_cap"] is None and not sample22["orbit_or_output_tables"])

print("[23] contraction-aware norm enclosures preserve absolute accuracy")
try:
    from flint import arb, ctx, fmpq
except ImportError:
    print("  SKIP norm enclosure tests: optional python-flint==0.9.0 required")
else:
    norm23 = VerifiedReflectionCircuit(10, 4, dict(verified22.backgrounds),
                                       dict(verified22.reflections), enclosure_mode="norm")
    for label23 in ((0,0,0,"arithmetic",0,0), (1,3,1,"background",0,0),
                    (4,7,3,"reflection",0,0), (1,0,4,"reflection",2,3),
                    (4,0,4,"reflection",4,9)):
        a23,e23,s23,bnd23,m23,y23 = label23
        dy23 = norm23.prefix_dyadic(a23,e23,s23,accuracy_bits=40,
                                   boundary=bnd23,measured=m23,output=y23)
        with ctx.workprec(192):
            reference23 = verified22.prefix_enclosure(a23,e23,s23,boundary=bnd23,
                                    measured=m23,output=y23,working_precision=192)
            upper23 = max(binary_fraction(abs(arb(fmpq(i23,1 << 40))-z23).upper())
                          for pair23,v23 in zip(dy23.coordinates,reference23)
                          for i23,z23 in zip(pair23,(v23.real,v23.imag)))
        t("norm-mode dyadic coordinates meet outward high-precision reference bound",
          upper23 <= dy23.coordinate_error_bound)
    norm_wide23 = VerifiedReflectionCircuit(wide22.period,wide22.width,
        dict(wide22.backgrounds),dict(wide22.reflections),enclosure_mode="norm")
    sample23 = norm_wide23.sample(Random(624),target_tv=Fraction(1,1_000_000))
    t("both enclosure geometries expose the same law-level contract",
      sample23["total_tv_upper_bound"] == sample22["total_tv_upper_bound"]
      and sample23["enclosure_mode"] == "norm"
      and sample22["enclosure_mode"] == "rectangular")
    t("norm bookkeeping removes retries on the fixed wide regression, not universally",
      sample23["refinement_retries"] == 0 and sample22["refinement_retries"] > 0)
try:
    VerifiedReflectionCircuit(2,0,{}, {},enclosure_mode="unsafe")
    t("unknown enclosure mode rejected",False)
except ValueError:
    t("unknown enclosure mode rejected",True)

print("[24] verified rejection charges proposal and accepted-measure errors")
from lab.verified_rejection import VerifiedRejectionSampler, conservative_threshold
t("one-sided exact integer acceptance handles boundaries without zero division",
  conservative_threshold(0, 0, 16) == 0
  and conservative_threshold(-1, 3, 16) == 0
  and conservative_threshold(1, 1, 16) == 1 << 16
  and conservative_threshold(Fraction(1,3), 1, 4) == 5)
rejection24 = VerifiedRejectionSampler(verified22)
for delta24 in (Fraction(1,4), Fraction(1,1000), Fraction(1,10**12)):
    plan24 = rejection24.plan(delta24)
    t("rejection planner charges history/component proposal, interval and random bits",
      plan24["total_tv_upper_bound"] <= delta24
      and plan24["proposal_tv_upper_bound"] >= plan24["component_target_tv"]
      and plan24["success_probability_lower_bound"] > 0
      and plan24["expected_attempts_upper_bound"] * plan24["success_probability_lower_bound"] == 1)
for invalid24 in (lambda: rejection24.plan(.001), lambda: rejection24.plan(1),
                  lambda: rejection24.component(4),
                  lambda: conservative_threshold(1, 0, 8),
                  lambda: conservative_threshold(0, -1, 8)):
    try:
        invalid24()
        t("unsupported rejection input rejected", False)
    except ValueError:
        t("unsupported rejection input rejected", True)
try:
    from flint import arb, ctx, fmpq
except ImportError:
    print("  SKIP verified rejection arithmetic tests: optional python-flint==0.9.0 required")
else:
    exactcoins24 = VerifiedRejectionSampler(VerifiedReflectionCircuit(6,3,{},
                                            {0:(0,0),1:(1,1),2:(2,2),3:(1,3)}))
    t("exact rational-pi zero and unit history probabilities retain endpoints",
      exactcoins24.history_thresholds(20)["thresholds"] == (0,1 << 20,0,1 << 20))
    component24 = rejection24.component(3)
    physical24 = VerifiedReflectionCircuit(10,4,dict(verified22.backgrounds),
                                          {2:(0,1),3:(1,1)})
    with ctx.workprec(192):
        for a24,e24,s24,m24,y24 in ((0,3,2,0,0),(2,7,3,0,0),(4,0,4,4,9)):
            v24 = component24.prefix_enclosure(a24,e24,s24,measured=m24,output=y24,
                                               working_precision=192)
            w24 = physical24.prefix_enclosure(a24,e24,s24,measured=m24,output=y24,
                                              working_precision=192)
            # Omitted K(pi) factors are global: compare Born weights, NOT
            # complex amplitudes that legitimately differ by powers of -i.
            diff24 = max(binary_fraction(abs(abs(v)**2-abs(w)**2).upper()) for v,w in zip(v24,w24))
            t("single-history prefix has the K(pi) circuit's Born weights", diff24 < Fraction(1,1 << 170))
    zero24 = VerifiedRejectionSampler(VerifiedReflectionCircuit(6,2,{}, {2:(1,0)}))
    decision24 = zero24.acceptance(0,1,interval_tolerance=Fraction(1,10**9),random_bits=30)
    t("exact-zero target/proposal rejection decision terminates with zero acceptance",
      decision24.probability == 0 and decision24.width_sum <= Fraction(1,10**9))
    sampled24 = rejection24.sample(Random(630),target_tv=Fraction(1,1_000_000))
    t("rejection returns only certified joint labels and charges all proposals",
      sampled24["total_tv_upper_bound"] <= Fraction(1,1_000_000)
      and "final_within_sector" not in sampled24
      and sampled24["acceptance_enclosure_evaluations"] >= sampled24["attempts"] >= 1
      and sampled24["attempt_cap"] is None and sampled24["component_history_count"] == 1)
print("[25] verified finite-work and scalar component proposals")
from lab.sampling_error import plan_prefix_mass_accuracy,plan_conditional_weight_accuracy
from lab.verified_finite_work import VerifiedFiniteWork,VerifiedScalarWork
mass25 = plan_prefix_mass_accuracy(63,Fraction(1,1_000_000))
scalar25 = plan_conditional_weight_accuracy(63,Fraction(1,1_000_000))
t("normalized scalar oracle needs fewer requested weight bits than prefix masses",
  mass25["total_tv_upper_bound"] <= Fraction(1,1_000_000)
  and scalar25["total_tv_upper_bound"] <= Fraction(1,1_000_000)
  and scalar25["weight_accuracy_bits"] < mass25["weight_accuracy_bits"])
for bad25 in (lambda: VerifiedFiniteWork(verified22,True),
              lambda: VerifiedScalarWork(verified22,4),
              lambda: plan_prefix_mass_accuracy(4,.01),
              lambda: plan_conditional_weight_accuracy(4,1),
              lambda: VerifiedRejectionSampler(verified22,component_method="float")):
    try:
        bad25(); t("invalid finite-work/scalar request rejected",False)
    except ValueError:
        t("invalid finite-work/scalar request rejected",True)
try:
    from flint import arb,ctx,acb_mat
except ImportError:
    print("  SKIP verified finite-work/scalar arithmetic: optional backend required")
else:
    edge25 = VerifiedReflectionCircuit(6,0,{0:("x",Fraction(1,4))},{0:(1,Fraction(1,3))})
    for factory25 in (VerifiedFiniteWork,VerifiedScalarWork):
        result25 = factory25(edge25,1).path(0,output=0)
        t("zero-width deterministic route preserves normalized output",
          result25["conditional_path_probability"] == 1 and result25["final_coarse_sector"] == 2
          and result25["total_tv_upper_bound"] == 0)
    # Include the observable early mixer; late-only fixtures are insufficient.
    full25 = VerifiedReflectionCircuit(10,4,{0:("x",Fraction(1,4)),**dict(verified22.backgrounds)},dict(verified22.reflections))
    for factory25 in (VerifiedFiniteWork,VerifiedScalarWork):
        for mask25 in (0,3):
            sampler25 = factory25(full25,mask25)
            laws25 = [sampler25.path(1,output=y,target_tv=Fraction(1,1000))["conditional_path_probability"] for y in range(16)]
            t("early-mixer component finite-bit law normalizes exactly",sum(laws25)==1)
    cursor25 = VerifiedFiniteWork(full25,2).cursor(1,accuracy_bits=30)
    cursor25.weights(); cursor25.advance(0)
    with ctx.workprec(cursor25.precision):
        cursor25.effect += acb_mat([[arb(0,1),0],[0,arb(0,1)]])
    weights25 = cursor25.weights()
    t("valid loose effect triggers precision rebuild with chosen prefix replay",
      cursor25.builds > 1 and cursor25.replayed >= 1 and cursor25.measured == 1 and cursor25.output == 0)
    s25 = VerifiedRejectionSampler(full25,component_method="scalar").sample(Random(640))
    t("scalar rejection preserves whole-law certificate and performs no prefix queries",
      s25["total_tv_upper_bound"] <= Fraction(1,1_000_000)
      and s25["component_scalar_evaluations"] >= 4*s25["attempts"]
      and s25["prefix_vector_evaluations"] == 0 and s25["component_forward_steps"] == 0)

print("[26] opt-in exact three-state work blocks retain dimension-aware certification")
odd26 = VerifiedReflectionCircuit(9,4,{1:("x",Fraction(1,7)),3:("z",Fraction(1,5))},
    {2:(0,Fraction(1,5)),3:(1,Fraction(1,5))},block_size=3)
plan26 = odd26.accuracy_plan(Fraction(1,1000))
from lab.verified_prefix import _apply_background
t("embedded two-pi rotation is not a global phase on the three-state block",
  _apply_background([1,2,3],"x",-1,0) == [-1,-2,3])
t("three-state input/planner charges six coordinates and three-bin CDFs",
  odd26.b == odd26.sectors == 3
  and plan26["scaled_coordinate_to_state_l2_factor"] == 10
  and plan26["finite_random_bits_tv_upper_bound"] == Fraction(2*plan26["updates"],1 << plan26["random_bits"]))
for invalid26 in (lambda: VerifiedReflectionCircuit(9,4,{},{}),
                  lambda: VerifiedReflectionCircuit(10,4,{},{},block_size=3),
                  lambda: VerifiedReflectionCircuit(9,4,{},{},block_size=True),
                  lambda: VerifiedReflectionCircuit(12,4,{},{},block_size=4),
                  lambda: VerifiedScalarWork(odd26,0),
                  lambda: VerifiedRejectionSampler(odd26,component_method="scalar")):
    try:
        invalid26(); t("invalid block/scalar specialization rejected",False)
    except ValueError:
        t("invalid block/scalar specialization rejected",True)
snapshot26 = VerifiedRejectionSampler(odd26)
t("rejection snapshots and deterministic components preserve block dimension",
  snapshot26.circuit.b == snapshot26.component(3).b == 3)
try:
    from flint import ctx,arb,fmpq
except ImportError:
    print("  SKIP verified odd-block arithmetic: optional backend required")
else:
    norm26 = VerifiedReflectionCircuit(9,4,dict(odd26.backgrounds),dict(odd26.reflections),
                                      block_size=3,enclosure_mode="norm")
    for c26 in (odd26,norm26):
        dy26 = c26.prefix_dyadic(1,0,4,measured=4,output=3,accuracy_bits=40)
        with ctx.workprec(192):
            ref26 = odd26.prefix_enclosure(1,0,4,measured=4,output=3,working_precision=192)
            error26 = max(binary_fraction(abs(arb(fmpq(n,1 << 40))-z).upper())
                          for pair,v in zip(dy26.coordinates,ref26)
                          for n,z in zip(pair,(v.real,v.imag)))
        t("rectangular/norm odd-block dyadics cover all three complex coordinates",
          len(dy26.coordinates) == 3 and error26 <= dy26.coordinate_error_bound)
    worker26 = VerifiedFiniteWork(odd26,3)
    forced26 = [worker26.path(1,output=y,target_tv=Fraction(1,1000)) for y in range(16)]
    t("three-state forward/effect forced law normalizes with correct matrix count",
      sum(r["conditional_path_probability"] for r in forced26) == 1
      and all(r["work_block_size"] == 3 and r["stored_matrix_scalar_count_upper_bound"] == 9*r["stored_matrix_count_upper_bound"] for r in forced26))
    for method26 in ("prefix","finite_work"):
        sample26 = VerifiedRejectionSampler(odd26,component_method=method26).sample(Random(651))
        t("three-state accepted output carries the requested full-law bound",
          sample26["total_tv_upper_bound"] <= Fraction(1,1_000_000)
          and 0 <= sample26["final_coarse_sector"] < 3 and 0 <= sample26["output"] < 16)

print("\n[27] opt-in checkpoint/recompute preserves the finite-work oracle")
for spacing27 in (0,-1,True,1.0,64):
    try:
        VerifiedFiniteWork(odd26,0,checkpoint_spacing=spacing27)
        t("invalid checkpoint spacing rejected",False)
    except ValueError:
        t("invalid checkpoint spacing rejected",True)
try:
    from flint import ctx,arb,acb_mat
except ImportError:
    print("  SKIP verified checkpoint arithmetic: optional backend required")
else:
    from weakref import ref
    for b27,r27 in ((2,10),(3,9)):
        for width27 in (0,1,3,5):
            refs27 = {0:(1,Fraction(1,5))}
            if width27:
                refs27[width27] = (1,Fraction(1,5))
            circuit27 = VerifiedReflectionCircuit(r27,width27,
                {0:("x",Fraction(1,7)),width27:("z",Fraction(1,5))},refs27,
                block_size=b27)
            mask27 = (1 << len(refs27))-1
            full27 = VerifiedFiniteWork(circuit27,mask27)
            for spacing27 in (1,2,4,63):
                compact27 = VerifiedFiniteWork(circuit27,mask27,checkpoint_spacing=spacing27)
                for y27 in {0,(1 << width27)-1}:
                    a27 = full27.path(1,output=y27)
                    z27 = compact27.path(1,output=y27)
                    t("checkpoint endpoints and complete forced-path probabilities match the default",
                      a27["conditional_path_probability"] == z27["conditional_path_probability"]
                      and a27["final_coarse_sector"] == z27["final_coarse_sector"]
                      and z27["total_tv_upper_bound"] <= Fraction(1,1000000)
                      and z27["forward_steps"] <= 2*width27*z27["forward_builds"]
                      and z27["stored_branch_matrix_count"] == 0)
    compact27 = VerifiedFiniteWork(odd26,0,checkpoint_spacing=2)
    cursor27 = compact27.cursor(0,accuracy_bits=40,initial_precision=64)
    cursor27.weights(); cursor27.advance(0)
    with ctx.workprec(64):
        for i27 in range(3):
            cursor27.effect[i27,i27] += arb(0,1)
    refined27 = cursor27.weights()
    fresh27 = compact27.cursor(0,accuracy_bits=40,initial_precision=cursor27.precision)
    fresh27.weights(); fresh27.advance(0)
    t("checkpoint refinement reconstructs the selected effect and charges its replay",
      cursor27.builds >= 2 and cursor27.replayed >= 1
      and refined27 == fresh27.weights()
      and cursor27.forward_steps >= odd26.width*cursor27.builds)
    live27 = ref(cursor27)
    del cursor27
    t("checkpoint storage does not retain its cursor through an ownership cycle",live27() is None)

print("\n[28] verified reverse vectors charge proposal and boundary acceptance")
from lab.verified_reverse_work import VerifiedReverseWork,_norm_squared
reverse28 = VerifiedReverseWork(odd26,3)
for invalid28 in (lambda: VerifiedReverseWork(odd26,True),
                  lambda: reverse28.plan(0),lambda: reverse28.plan(1),
                  lambda: reverse28.plan(1e-6)):
    try:
        invalid28(); t("invalid reverse input or accuracy rejected",False)
    except ValueError:
        t("invalid reverse input or accuracy rejected",True)
for b28,r28 in ((2,10),(3,9)):
    for width28 in (0,4,63):
        worker28 = VerifiedReverseWork(VerifiedReflectionCircuit(r28,width28,{},{},block_size=b28),0)
        for delta28 in (Fraction(1,1000),Fraction(1,10**15)):
            plan28 = worker28.plan(delta28)
            t("reverse planner covers all errors and a positive per-attempt success",
              0 <= plan28["total_tv_upper_bound"] <= 3*delta28/4
              and plan28["expected_attempts_upper_bound"] == 1/plan28["success_probability_lower_bound"])
try:
    from flint import ctx,arb,acb,acb_mat
except ImportError:
    print("  SKIP verified reverse arithmetic: optional backend required")
else:
    t("polynomial norm remains finite for coordinate intervals crossing zero",
      _norm_squared(acb_mat([[acb(arb(0,1),arb(0,1))],[0]])).is_finite())
    # Both refinement sites must preserve the same measured prefix and j.
    c28 = reverse28.cursor(1,2,accuracy_bits=40,initial_precision=64)
    c28.weights(); c28.advance(1)
    with ctx.workprec(64):
        c28.vector[0,0] += arb(0,1)
    weights28 = c28.weights()
    fresh28 = reverse28.cursor(1,2,accuracy_bits=40,initial_precision=c28.precision)
    fresh28.weights(); fresh28.advance(1)
    t("reverse child refinement replays fixed boundary and selected bits",
      c28.builds >= 2 and c28.replayed >= 1 and c28.boundary == 2
      and c28.output == 1 and weights28 == fresh28.weights())
    c28.advance(0); fresh28.advance(0)
    for bit28 in (1,0):
        c28.weights(); c28.advance(bit28)
        fresh28.weights(); fresh28.advance(bit28)
    builds28,replayed28 = c28.builds,c28.replayed
    with ctx.workprec(c28.precision):
        c28.vector[0,0] += arb(0,1)
    a28 = c28.acceptance(interval_tolerance=Fraction(1,10**12),random_bits=40)
    again28 = reverse28.cursor(1,2,accuracy_bits=40,initial_precision=c28.precision)
    for bit28 in (1,0,1,0):
        again28.weights(); again28.advance(bit28)
    b28 = again28.acceptance(interval_tolerance=Fraction(1,10**12),random_bits=40)
    t("terminal acceptance refinement replays the full fixed proposal",
      c28.builds > builds28 and c28.replayed >= replayed28+4
      and c28.output == 5 and c28.boundary == 2 and a28.threshold == b28.threshold)
    zero_terminal_count28 = 0
    for block28,period28 in ((2,6),(3,9)):
        zero28 = VerifiedReverseWork(VerifiedReflectionCircuit(period28,0,{},
            {0:(1,Fraction(1,5))},block_size=block28),1)
        paths28 = [zero28.attempt(1,j28,output=0) for j28 in range(block28)]
        t("zero-width boundary rejection has exact ideal acceptance mass and route",
          all(p["proposal_path_probability"] == 1 and p["final_coarse_sector"] == 1 for p in paths28)
          and sum(p["acceptance_probability"] for p in paths28) == 1)
        null28 = VerifiedReverseWork(VerifiedReflectionCircuit(period28,2,{},{},block_size=block28),0)
        paths28 = [null28.attempt(0,j28,output=y28) for j28 in range(block28) for y28 in range(4)]
        zero_terminal_count28 += sum(p["acceptance_denominator_bounds"][1] == 0 for p in paths28)
        t("all exact-zero reverse terminal norms are rejected without a ratio limit",
          all(p["acceptance_probability"] == 0 for p in paths28 if p["acceptance_denominator_bounds"][1] == 0)
          and sum(p["proposal_path_probability"] for p in paths28) == block28)
    t("reverse zero-denominator regression is nonvacuous",zero_terminal_count28 > 0)
    from unittest.mock import patch
    with patch.object(VerifiedFiniteWork,"_build",side_effect=AssertionError("hidden forward history")):
        sample28 = reverse28.sample(Random(6828))
    t("reverse returned output charges all rejected trajectories and no forward states",
      sample28["total_tv_upper_bound"] <= Fraction(1,1000000)
      and sample28["reverse_child_evaluations"] >= 4*sample28["attempts"]
      and sample28["stored_forward_matrix_count"] == sample28["stored_branch_matrix_count"] == 0
      and not sample28["vector_normalization"] and "boundary_j" not in sample28)

print("\n[29] local instrument error supports bounded-grid integer trajectories")
from lab.verified_quantized_reverse import (VerifiedQuantizedReverseWork,canonical_quantize,
    integer_norm,integer_matvec,integer_overlap_numerator)
q29,z29 = canonical_quantize(((-7,7),(3,-4)),4)
t("canonical rounding preserves signed/tied maximal coordinates exactly",
  q29 == ((-16,16),(7,-9)) and not z29
  and canonical_quantize(((-700,700),(300,-400)),4)[0] == q29
  and max(abs(x) for pair in q29 for x in pair) == 16)
t("zero child has an explicit diagnostic-only canonical fallback",
  canonical_quantize(((0,0),(0,0)),4) == (((16,0),(0,0)),True))
t("integer norm and Hermitian overlap use the correct conjugation",
  integer_norm(((3,4),(0,1))) == 26
  and integer_overlap_numerator(((1,1),(0,0)),((1,1),(0,0))) == 4)
for bad29 in (lambda: canonical_quantize((),4),
              lambda: canonical_quantize(((1.0,0),),4),
              lambda: canonical_quantize(((1,0),),0)):
    try:
        bad29(); t("invalid canonical grid inputs rejected",False)
    except ValueError:
        t("invalid canonical grid inputs rejected",True)
for block29,period29 in ((2,10),(3,9)):
    plans29 = [VerifiedQuantizedReverseWork(VerifiedReflectionCircuit(period29,width29,{},{},block_size=block29),0).plan(Fraction(1,1000000)) for width29 in range(64)]
    t("all bounded widths have a valid local-error and integer-size plan",
      all(p["total_tv_upper_bound"] <= Fraction(1,1000000)
          and p["stacked_operator_error_upper_bound"] < 1
          and p["expected_attempts_upper_bound"] == 1/p["success_probability_lower_bound"]
          for p in plans29)
      and all(x["grid_bits"] <= y["grid_bits"] for x,y in zip(plans29,plans29[1:])))
try:
    from flint import ctx,arb
except ImportError:
    print("  SKIP quantized reverse oracle tests: optional backend required")
else:
    qworker29 = VerifiedQuantizedReverseWork(odd26,3)
    c29 = qworker29.cursor(1,2,target_tv=Fraction(1,1000000))
    c29.weights(); c29.advance(1)
    state29,measured29,output29 = c29.vector,c29.measured,c29.output
    original29 = qworker29.builder._branch
    calls29 = []
    def loose_once29(*args):
        pair29 = original29(*args)
        if not calls29:
            pair29[0][0,0] += arb(0,1)
        calls29.append(1)
        return pair29
    with patch.object(qworker29.builder,"_branch",side_effect=loose_once29):
        weights29 = c29.weights()
    t("local oracle refinement preserves the previous integer state and selected bits",
      len(calls29) >= 2 and c29.refinements >= 1
      and (c29.vector,c29.measured,c29.output) == (state29,measured29,output29)
      and sum(weights29) > 0)
    p29 = c29.plan
    t("integer operands obey explicit bounds after a nontrivial rounded prefix",
      c29.max_state_bits <= p29["state_coordinate_bits_upper_bound"]
      and c29.initial_coordinate_bits <= p29["initial_coordinate_bits_upper_bound"]
      and c29.max_child_bits <= p29["child_coordinate_bits_upper_bound"])
    # An exactly zero first branch is forceable but not selectable by the CDF.
    zero29 = VerifiedQuantizedReverseWork(VerifiedReflectionCircuit(6,2,{},{}),0)
    zero_cursor29 = zero29.cursor(0,0)
    zero_weights29 = zero_cursor29.weights()
    zero_counts29 = integer_cdf_counts(zero_weights29,12)
    forced29 = zero29.attempt(0,0,output=1)
    t("zero child cannot be sampled, while forced zero paths remain zero probability",
      zero_weights29[1] == zero_counts29[1] == 0
      and sum(zero_counts29) == 1 << 12
      and forced29["proposal_path_probability"] == 0
      and forced29["diagnostic_zero_child_fallbacks"] >= 1)
    rare29 = VerifiedQuantizedReverseWork(VerifiedReflectionCircuit(2*((1 << 20)-1),3,{},{}),0)
    rare_cursor29 = rare29.cursor(1,0)
    rare_weights29 = rare_cursor29.weights()
    t("a positive very rare physical child has bounded canonical work coordinates",
      0 < Fraction(min(rare_weights29),sum(rare_weights29)) < Fraction(1,10**8))
    rare_cursor29.advance(int(rare_weights29[1] < rare_weights29[0]))
    t("rare-child compression retains a maximal coordinate without relative-norm refinement",
      max(abs(x) for z in rare_cursor29.vector for x in z) == 1 << rare_cursor29.plan["grid_bits"]
      and rare_cursor29.refinements == 0 and rare_cursor29.zero_children == 0)
    with patch.object(VerifiedFiniteWork,"_build",side_effect=AssertionError("hidden forward history")):
        sample29 = qworker29.sample(Random(7029))
    t("quantized sample charges all attempts, no replay or forward history, bounded state integers",
      sample29["total_tv_upper_bound"] <= Fraction(1,1000000)
      and sample29["vector_compressions"] == odd26.width*sample29["attempts"]
      and sample29["max_state_coordinate_bits"] <= sample29["state_coordinate_bits_upper_bound"]
      and sample29["max_initial_coordinate_bits"] <= sample29["initial_coordinate_bits_upper_bound"]
      and sample29["stored_forward_matrix_count"] == sample29["stored_branch_matrix_count"] == 0
      and sample29["replayed_vector_steps"] == sample29["diagnostic_zero_child_fallbacks"] == 0)

print("\n[30] sparse coherent reverse paths merge amplitudes, not histories")
from lab.coherent_reverse import SparseCoherentReverse,route_support_bound
from unittest.mock import patch
t("route support bound distinguishes empty, one-label, two-label and generic inputs",
  route_support_bound(1009,()) == 1
  and route_support_bound(1009,(2,)*8) == 2
  and route_support_bound(1009,(0,1)*4) == 17
  and route_support_bound(1009,(0,1,0,3,0,9,0,27)) == 256
  and route_support_bound(1,(0,1)*4) == 1)
W30 = np.eye(3,dtype=complex)
W30[:2,:2] = [[np.cos(np.pi/7),-1j*np.sin(np.pi/7)],
               [-1j*np.sin(np.pi/7),np.cos(np.pi/7)]]
circuit30 = CoherentReflectionCircuit(9,3,3,{0:W30,1:W30,3:W30.conj().T},
    {0:(0,np.pi/5),2:(1,np.pi/3),3:(0,np.pi/7)})
worker30 = SparseCoherentReverse(circuit30)
mass30,amplitude_error30,proposal_error30 = 0.,0.,0.
for gamma30 in range(3):
    for j30 in range(3):
        prop30 = 0.
        for y30 in range(8):
            row30 = worker30.attempt(gamma30,j30,output=y30)
            target30 = circuit30.prefix_vector(gamma30,0,3,measured=3,output=y30)[j30]
            amplitude_error30 = max(amplitude_error30,
                abs(row30["terminal_coherent_amplitude"].conjugate()-target30))
            prop30 += row30["proposal_path_probability"]
            mass30 += row30["accepted_joint_submass"]
        proposal_error30 = max(proposal_error30,abs(prop30-1))
t("complex reverse amplitudes and each boundary's full proposal law match the history oracle",
  amplitude_error30 < 2e-13 and proposal_error30 < 2e-13
  and abs(mass30-1/(3*worker30.support_bound)) < 2e-13)
with patch.object(CoherentReflectionCircuit,"_histories",side_effect=AssertionError("hidden history expansion")):
    sample30 = worker30.sample(np.random.default_rng(7130))
t("returned sparse samples charge all attempts and never enumerate histories",
  sample30["reverse_steps"] == 3*sample30["rejection_proposals"]
  and sample30["peak_sector_count"] <= worker30.support_bound
  and sample30["history_enumerations"] == sample30["orbit_table_entries"] == 0
  and not sample30["precision_certificate"])
zero30 = SparseCoherentReverse(CoherentReflectionCircuit(6,2,2,{},{}))
forced30 = zero30.attempt(0,0,output=1)
t("forced exact-zero reverse prefix stays zero accepted mass without sampled fallback",
  forced30["proposal_path_probability"] == forced30["accepted_joint_submass"] == 0
  and forced30["forced_zero_prefixes"] >= 1)
rejected30 = dict(row30,acceptance_probability=0.)
try:
    with patch.object(worker30,"attempt",return_value=rejected30):
        worker30.sample(np.random.default_rng(7131),max_proposals=1)
    t("explicit float diagnostic cap raises instead of returning a fallback",False)
except RuntimeError:
    t("explicit float diagnostic cap raises instead of returning a fallback",True)
for bad30 in (lambda: route_support_bound(0,()),
              lambda: route_support_bound(3,(0,)*65),
              lambda: worker30.attempt(3,0,output=0),
              lambda: worker30.attempt(0,3,output=0),
              lambda: worker30.attempt(0,0,output=8)):
    try:
        bad30(); t("invalid sparse reverse inputs rejected",False)
    except ValueError:
        t("invalid sparse reverse inputs rejected",True)

print("\n[31] merged intermediate prefixes reuse the no-rejection sampler")
from lab.merged_prefix import MergedCoherentPrefixes
W31 = np.diag([np.exp(-1j*np.pi/5),np.exp(1j*np.pi/5),1.])
circuit31 = CoherentReflectionCircuit(9,3,3,{0:W30,1:W31,3:W30},
    {0:(0,np.pi/5),2:(1,np.pi/3),3:(0,np.pi/7)})
merged31 = MergedCoherentPrefixes(circuit31)
error31,cost_ok31,requests31 = 0.,True,0
for gamma31 in range(3):
    for exponent31 in range(8):
        requests31_list = [(s,z,0,0) for s in range(4)
                           for z in ("arithmetic","background","reflection")]
        requests31_list += [(3,"reflection",ell,y) for ell in range(1,4)
                            for y in range(1 << ell)]
        for stop31,boundary31,measured31,output31 in requests31_list:
            actual31 = merged31.prefix_vector(gamma31,exponent31,stop31,
                boundary=boundary31,measured=measured31,output=output31)
            reference31 = circuit31.prefix_vector(gamma31,exponent31,stop31,
                boundary=boundary31,measured=measured31,output=output31)
            error31 = max(error31,float(np.max(np.abs(actual31-reference31))))
            requests31 += 1
            reached31 = {gamma31}
            expected31 = dict(background_block_products=0,control_block_products=0,
                fixed_control_block_scales=0,qft_matrix_pair_constructions=0,
                reflection_block_contributions=0,peak_sector_count=1)
            for s31 in range(stop31,-1,-1):
                if s31 in circuit31.reflections and (s31 < stop31 or boundary31=="reflection"):
                    expected31["reflection_block_contributions"] += 2*len(reached31)
                    q31 = circuit31.reflections[s31][0]
                    reached31 |= {(-a-q31)%3 for a in reached31}
                    expected31["peak_sector_count"] = max(expected31["peak_sector_count"],len(reached31))
                if s31 < stop31 or boundary31!="arithmetic":
                    expected31["background_block_products"] += len(reached31)
                if s31:
                    if s31-1 >= 3-measured31:
                        expected31["control_block_products"] += len(reached31)
                        expected31["qft_matrix_pair_constructions"] += len(reached31)
                    else:
                        expected31["control_block_products"] += bool(exponent31 & (1 << (s31-1)))*len(reached31)
                        expected31["fixed_control_block_scales"] += len(reached31)
            local31 = merged31.last_prefix_stats
            cost_ok31 &= all(local31[key] == value for key,value in expected31.items())
            cost_ok31 &= local31["peak_sector_count"] <= local31["prefix_support_bound"]
t("every tiny intermediate and partial-QFT complex vector matches the history reference",
  requests31 == 624 and error31 < 3e-13)
t("batched operation counters match an independent reached-set recurrence",
  cost_ok31 and merged31.stats()["merged_prefix_queries"] == requests31)
for seed31 in (7231,7232):
    baseline31 = circuit31.sample(np.random.default_rng(seed31))
    with patch.object(CoherentReflectionCircuit,"_histories",side_effect=AssertionError("hidden history")):
        sampled31 = merged31.sample(np.random.default_rng(seed31))
    t("the reused sampler resets counts, keeps the seeded path and never rejects or enumerates histories",
      all(sampled31[key] == baseline31[key] for key in
          ("output","final_coarse_sector","final_within_sector","prefix_vector_evaluations"))
      and sampled31["merged_prefix_queries"] == sampled31["prefix_vector_evaluations"]
      and sampled31["history_enumerations"] == sampled31["rejection_proposals"] == 0
      and not sampled31["precision_certificate"])
zero31 = MergedCoherentPrefixes(CoherentReflectionCircuit(6,2,2,{},{}))
t("an exactly zero measured prefix stays zero without normalization",
  np.array_equal(zero31.prefix_vector(0,0,2,measured=2,output=1),np.zeros(2)))
snap31 = MergedCoherentPrefixes(circuit31)
before31 = snap31.prefix_vector(0,0,0,boundary="background")
circuit31.background.defects[0] = np.eye(3)
t("merged adapter snapshots input blocks instead of aliasing caller mutations",
  np.array_equal(before31,snap31.prefix_vector(0,0,0,boundary="background")))
for invalid31 in (lambda: merged31.prefix_vector(3,0,0),
                  lambda: merged31.prefix_vector(0,8,0),
                  lambda: merged31.prefix_vector(0,0,4),
                  lambda: merged31.prefix_vector(0,0,0,boundary="invalid"),
                  lambda: merged31.prefix_vector(0,0,2,measured=1),
                  lambda: merged31.prefix_vector(0,0,3,measured=1,output=2)):
    try:
        invalid31(); t("merged oracle shares the exact prefix input contract",False)
    except ValueError:
        t("merged oracle shares the exact prefix input contract",True)

print("\n[32] fixed alphabets and the explicit structural input boundary")
from lab.coherent_routes import CoherentReflectionInput
t("three-label L1 cover improves the history bound beyond the old cap",
  route_support_bound(1_000_003,tuple((0,1,101)[i%3] for i in range(9))) == 362
  and route_support_bound(1_000_003,tuple((0,1,101)[i%3] for i in range(32))) == 4226
  and route_support_bound(1009,(0,1)*32) == 129)
refs32 = {s:(s%3,np.pi/7) for s in range(9)}
refs32[8] = (2,0.)  # ninth insertion is identity; old reference omits it
input32 = CoherentReflectionInput(9,3,8,{0:W30,3:W31,8:W30},refs32)
old32 = CoherentReflectionCircuit(9,3,8,input32.background.defects,
    {s:r for s,r in refs32.items() if s!=8})
adapter32 = MergedCoherentPrefixes(input32)
error32 = 0.
for gamma32 in range(3):
    for exponent32 in (0,73,255):
        for s32,z32,m32,y32 in ((0,"arithmetic",0,0),(0,"reflection",0,0),
                (3,"background",0,0),(8,"arithmetic",0,0),
                (8,"reflection",0,0),(8,"reflection",4,7),(8,"reflection",8,197)):
            reference32 = old32.prefix_vector(gamma32,exponent32,s32,
                boundary=z32,measured=m32,output=y32)
            with patch.object(CoherentReflectionCircuit,"_histories",side_effect=AssertionError("hidden histories")):
                got32 = adapter32.prefix_vector(gamma32,exponent32,s32,
                    boundary=z32,measured=m32,output=y32)
            error32 = max(error32,float(np.max(np.abs(got32-reference32))))
t("explicit ninth identity insertion matches the legacy eight-insertion amplitudes",error32<3e-13)
t("structural inputs and their owned snapshots have no legacy history/sampler API",
  type(adapter32.circuit) is CoherentReflectionInput
  and not any(hasattr(adapter32.circuit,key) for key in ("_histories","sample","sample_rejection")))
with patch.object(CoherentReflectionCircuit,"_histories",side_effect=AssertionError("hidden histories")):
    sample32 = adapter32.sample(np.random.default_rng(7332))
    reverse32 = SparseCoherentReverse(input32).sample(np.random.default_rng(7333),max_proposals=256)
t("both existing merged samplers consume structural input without histories",
  sample32["history_enumerations"] == reverse32["history_enumerations"] == 0
  and sample32["rejection_proposals"] == 0
  and reverse32["reverse_steps"] == 8*reverse32["rejection_proposals"])
saved32 = adapter32.prefix_vector(1,73,8)
input32.reflections[1] = (2,1.3)
input32.background.defects[0] = np.eye(3)
t("structural snapshots own both routing and work inputs",
  np.array_equal(saved32,adapter32.prefix_vector(1,73,8)))
for invalid32 in (lambda: CoherentReflectionCircuit(9,3,8,{},refs32),
                  lambda: next(CoherentReflectionCircuit._histories(input32,8,"reflection")),
                  lambda: CoherentReflectionInput(9,3,64,{},{}),
                  lambda: CoherentReflectionInput(9,3,8,{}, {9:(0,1.)}),
                  lambda: CoherentReflectionInput(9,3,8,{}, {0:(0,float("nan"))})):
    try:
        invalid32(); t("legacy history caps and structural input validation remain enforced",False)
    except ValueError:
        t("legacy history caps and structural input validation remain enforced",True)
large32 = CoherentReflectionInput(3_000_009,3,63,{},
    {i:(i%9,np.pi/4) for i in range(64)})
for adapter_cls32 in (SparseCoherentReverse,MergedCoherentPrefixes):
    try:
        adapter_cls32(large32); t("excess structural covers reject before sparse propagation",False)
    except MemoryError:
        t("excess structural covers reject before sparse propagation",True)
endpoint32 = CoherentReflectionInput(9,3,63,{}, {i:(0,0.) for i in range(64)})
end32 = MergedCoherentPrefixes(endpoint32)
endref32 = MergedCoherentPrefixes(CoherentReflectionCircuit(9,3,63,{},{}))
t("64 explicit identity rotations retain all endpoint semantics at the width cap",
  np.max(np.abs(end32.prefix_vector(1,0,63,measured=63,output=13)
                -endref32.prefix_vector(1,0,63,measured=63,output=13))) < 3e-13)

print("\n[33] word-specific global envelopes preserve merged sampling")
from lab.coherent_reverse import word_support_plan
plan33 = word_support_plan(101,tuple(i%3 for i in range(8)))
t("word setup counts the two affine orientations, not support at sector zero",
  plan33["support_bound"]==20 and plan33["alphabet_support_bound"]==101
  and plan33["word_even_offset_count"]==plan33["word_odd_offset_count"]==10
  and plan33["word_setup_copied_offsets"]==plan33["word_setup_reflected_offsets"]==71
  and plan33["word_setup_peak_live_offset_entries"]==38
  and plan33["word_retained_offset_entries"]==0)
input33 = CoherentReflectionInput(14,2,3,{0:np.array([[1,0],[0,1]],complex),
    1:np.array([[1,-1j],[-1j,1]],complex)/np.sqrt(2)},
    {0:(0,.3),1:(1,.7),3:(0,.2)})
default33 = SparseCoherentReverse(input33)
word33 = SparseCoherentReverse(input33,support_mode="word")
merged_default33 = MergedCoherentPrefixes(input33)
merged_word33 = MergedCoherentPrefixes(input33,support_mode="word")
t("word mode improves the unsaturated global envelope without changing the default",
  default33.support_bound==7 and word33.support_bound==6
  and default33._stats()["word_setup_steps"]==0
  and word33._stats()["word_setup_steps"]==3)
mass33 = 0.
for gamma33 in range(7):
    for fine33 in range(2):
        for y33 in range(8):
            a33 = default33.attempt(gamma33,fine33,output=y33)
            b33 = word33.attempt(gamma33,fine33,output=y33)
            assert abs(a33["terminal_coherent_amplitude"]-b33["terminal_coherent_amplitude"])<1e-14
            assert abs(7*a33["accepted_joint_submass"]-6*b33["accepted_joint_submass"])<1e-14
            mass33 += b33["accepted_joint_submass"]
t("all accepted submasses rescale by the new global bound without changing the target",
  abs(mass33-1/12)<3e-13)
for gamma33 in range(7):
    for exponent33 in (0,7):
        for stop33 in range(4):
            for boundary33 in ("arithmetic","background","reflection"):
                assert np.array_equal(merged_default33.prefix_vector(gamma33,exponent33,stop33,
                    boundary=boundary33),merged_word33.prefix_vector(gamma33,exponent33,stop33,
                    boundary=boundary33))
with patch("lab.coherent_reverse.word_support_plan",side_effect=AssertionError("repeated word setup")):
    p33 = merged_word33.sample(np.random.default_rng(7533))
    q33 = word33.sample(np.random.default_rng(7534),max_proposals=256)
t("one-time envelope setup is not repeated by prefix queries or rejection attempts",
  p33["word_setup_steps"]==q33["word_setup_steps"]==3
  and p33["rejection_proposals"]==0 and q33["reverse_steps"]==3*q33["rejection_proposals"])
saved_ref33 = input33.reflections.pop(0)
input33.reflections[0] = saved_ref33
reordered33 = SparseCoherentReverse(input33,support_mode="word")
t("word preprocessing and snapshot both use sorted original-time insertions",
  reordered33.support_bound==6 and tuple(reordered33.circuit.reflections)==(0,1,3))
for bad33 in (lambda: word_support_plan(101,(),max_offset_entries=0),
              lambda: word_support_plan(101,(),max_offset_entries=65_537),
              lambda: word_support_plan(101,(0,)*65),
              lambda: SparseCoherentReverse(input33,support_mode="gamma0"),
              lambda: MergedCoherentPrefixes(input33,support_mode=None)):
    try:
        bad33(); t("invalid word-envelope requests reject explicitly",False)
    except ValueError:
        t("invalid word-envelope requests reject explicitly",True)
try:
    word_support_plan(101,tuple(i%3 for i in range(8)),max_offset_entries=8)
    t("word live-offset cap rejects before an oversized update",False)
except MemoryError:
    t("word live-offset cap rejects before an oversized update",True)

print("[uniform output prefixes] exact modular near-collision certificate")
from lab.periodic import uniform_prefix_certificate
checked34=0
for period34 in range(1,17):
    for width34 in range(7):
        for prefix34 in range(width34+1):
            for radius34 in (0,1,2,6):
                result34=uniform_prefix_certificate(period34,width34,prefix34,radius34)
                H34,L34=1<<prefix34,1<<(width34-prefix34)
                direct34=sum(min((L34*q)%period34,(-L34*q)%period34)<=2*radius34
                             for q in range(1,H34))
                assert result34["near_collision_count"]==direct34
                assert result34["certified_uniform"]==(direct34==0)
                assert result34["vacuous_prefix"]==(prefix34==0)
                assert result34["floor_sum_calls"] in (0,2)
                assert result34["euclidean_iterations"]<=4*period34.bit_length()+4
                checked34+=1
t("support certificates match small brute modular differences including empty/full circles",
  checked34==16*28*4)
touch34=uniform_prefix_certificate(60,6,2,6)
trim34=uniform_prefix_certificate(60,6,2,4)
t("touching inclusive cones are not separated; removing a terminal unitary can sharpen the promise",
  touch34["near_collision_count"]==1 and not touch34["certified_uniform"]
  and trim34["near_collision_count"]==0 and trim34["certified_uniform"])
wide34=uniform_prefix_certificate(15,63,63,0)
t("wide output count uses exact integers without enumerating exponentially many histories",
  wide34["near_collision_count"]==((1<<63)-1)//15
  and wide34["high_history_count"]==1<<63
  and wide34["orbit_table_entries"]==wide34["prefix_table_entries"]==0
  and 0<wide34["euclidean_iterations"]<=20)
for bad34 in ((0,4,1,0),(1,64,1,0),(5,4,5,0),(5,4,-1,0),
              (5,4,1,-1),(1<<63,4,1,0),(5,4,1,1<<63),
              (True,4,1,0),(5,4,True,0),(5.,4,1,0),(5,4,1,[])):
    try:
        uniform_prefix_certificate(*bad34)
        t("invalid prefix certificate request rejects before arithmetic",False)
    except ValueError:
        t("invalid prefix certificate request rejects before arithmetic",True)

print("[work-first progressions] coherent rows and actual bounded samples")
from lab.work_first import LateWorkProgressions
init78 = np.array([[1,1,0],[1,-1,0],[0,0,np.sqrt(2)]],complex)/np.sqrt(2)
mid78 = np.array([[np.sqrt(2),0,0],[0,1,-1j],[0,-1j,1]],complex)/np.sqrt(2)
phases78 = (1,1j,-1)
model78 = LateWorkProgressions(3,3,3,2,init78,mid78,lambda j:phases78[j])
# Independent tiny literal matrix schedule, then numpy's Fourier transform.
# 8x3 columns and a 3x8 joint law only; no new general propagator.
cols78 = np.array([np.roll(np.array(phases78)*(mid78 @ np.roll(init78[:,0],e%4)),
                          4*(e//4)) for e in range(8)])
reference78 = np.abs(np.fft.fft(cols78,axis=0).T)**2/64
got78 = np.array([[model78.forced_joint(j,y)["joint_probability"]
                  for y in range(8)] for j in range(3)])
t("work-first joint law matches literal matrices and an independent FFT",
  np.max(np.abs(reference78-got78))<3e-16 and abs(got78.sum()-1)<1e-14)
for j78 in range(3):
    row78=model78.row(j78)
    expanded78={}
    for start78,count78,amp78 in row78["components"]:
        for k78 in range(count78):
            e78=start78+3*k78
            assert e78 not in expanded78
            expanded78[e78]=amp78
    assert all(abs(expanded78.get(e78,0j)-cols78[e78,j78])<3e-16 for e78 in range(8))
    accepted78=sum(model78.forced_joint(j78,y78)["proposal_probability"]
                   *model78.forced_joint(j78,y78)["acceptance"] for y78 in range(8))
    assert abs(accepted78-1/len(row78["components"]))<3e-16
t("small-period residue aliases merge before disjoint progression rejection",True)
draws78=[model78.sample(np.random.default_rng(7800+i78),max_attempts=256) for i78 in range(12)]
t("actual draws retain one column/work selection across all output retries",
  all(0<=d["output"]<8 and 0<=d["work"]<3
      and d["counters"]["column_local_terms"]==9 and d["counters"]["row_local_terms"]==18
      and d["counters"]["work_draws"]==1
      and d["counters"]["progression_proposals"]==d["counters"]["acceptance_draws"]==d["attempts"]
      and d["counters"]["fourier_component_terms"]==d["attempts"]*d["component_count"]
      and d["counters"]["marginal_queries"]==6*d["attempts"] for d in draws78)
  and any(d["attempts"]>1 for d in draws78))
zero78=LateWorkProgressions(4,1,1,1,np.ones((1,1)),np.ones((1,1)),lambda j:1)
t("unreachable work label has zero joint mass and no invented conditional law",
  zero78.row(3)["components"]==() and zero78.forced_joint(3,0)["joint_probability"]==0
  and zero78.forced_joint(3,0)["conditional_probability"] is None)
point78=LateWorkProgressions(1,1,3,3,np.ones((1,1)),np.ones((1,1)),lambda j:1)
t("a mathematical zero proposal has zero acceptance without a threshold",
  point78.forced_joint(0,1)["proposal_probability"]==0
  and point78.forced_joint(0,1)["acceptance"]==0
  and point78.sample(np.random.default_rng(78))["output"]==0)
for bad78 in (lambda:LateWorkProgressions(True,1,1,1,np.ones((1,1)),np.ones((1,1)),lambda j:1),
              lambda:LateWorkProgressions(3,2,3,2,init78,mid78,lambda j:1),
              lambda:LateWorkProgressions(3,3,64,63,init78,mid78,lambda j:1),
              lambda:LateWorkProgressions(3,3,3,2,init78,2*mid78,lambda j:1),
              lambda:model78.row(True),lambda:model78.column(8),
              lambda:model78.forced_joint(0,8),
              lambda:model78.sample(np.random.default_rng(1),max_attempts=0)):
    try:
        bad78(); t("invalid work-first parameters reject",False)
    except ValueError:
        t("invalid work-first parameters reject",True)
for kwargs78 in ({"max_payload_bytes":1},{"max_components":1},{"max_local_terms":1}):
    try:
        LateWorkProgressions(3,3,3,2,init78,mid78,lambda j:1,**kwargs78)
        t("structural limits reject before work-first allocation/loops",False)
    except MemoryError:
        t("structural limits reject before work-first allocation/loops",True)
with patch("lab.work_first.progression_sample",return_value={"output":1,"marginal_queries":6}):
    try:
        point78.sample(np.random.default_rng(78),max_attempts=2)
        t("exhausted conditional rejection raises with no replacement output",False)
    except RuntimeError:
        t("exhausted conditional rejection raises with no replacement output",True)
t("supplied matrices are copied read-only and no hidden large tables are owned",
  not model78.initial.flags.writeable and not np.shares_memory(model78.initial,init78)
  and model78.stats()["orbit_table_entries"]==model78.stats()["expanded_progression_entries"]==0)

# Complete first-attempt RNG traversal, independently derived: W0=I, W1=H,
# r=2, Q=4, split=1. Work j has probability 1/2, four singleton components,
# a uniform component Fourier proposal, and acceptance exactly at y=2*j.
# This covers actual interval bits AND gcd lifts, not just forced_joint.
# 128 calls, each capped at ONE proposal; 4x4 numeric output table only.
decision78=LateWorkProgressions(2,2,2,1,np.eye(2),
    np.array([[1,1],[1,-1]],complex)/np.sqrt(2),lambda j:1)
class DecisionRNG78:
    def __init__(self,e,j,c,y):
        self.ints=iter((e,y//2))
        self.reals=iter(((j+.5)/2,(c+.5)/4,((y%2)+.5)/2,.5))
        self.integer_bounds=[]
        self.real_calls=0
    def integers(self,high):
        self.integer_bounds.append(high)
        return next(self.ints)
    def random(self):
        self.real_calls+=1
        return next(self.reals)
accepted78=np.zeros((2,4))
rejected78=0.
traversed78=0
for e78 in range(4):
    for j78 in range(2):
        for c78 in range(4):
            for y78 in range(4):
                rng78=DecisionRNG78(e78,j78,c78,y78)
                traversed78+=1
                try:
                    returned78=decision78.sample(rng78,max_attempts=1)
                    assert y78==2*j78
                    assert (returned78["work"],returned78["output"])==(j78,y78)
                    accepted78[j78,y78]+=1/128
                except RuntimeError:
                    assert y78!=2*j78
                    rejected78+=1/128
                assert rng78.integer_bounds==[4,2] and rng78.real_calls==4
t("complete actual one-attempt transition law includes component, Fourier bit and gcd lift",
  traversed78==128 and np.array_equal(accepted78,np.array([[.125,0,0,0],[0,0,.125,0]]))
  and rejected78==.75)
t("control: returning an unaccepted component proposal gives the wrong joint law",
  np.sum(np.abs(accepted78*4-np.full((2,4),.125)))/2==.75)

print("[earlier-phase progressions] shared sampler with two row strides")
from lab.work_first import EarlierPhaseProgressions
joint_cases79=0
for period79 in (3,6):
    early_values79=np.array([1,1j,-1,-1j,1,1j][:period79])
    late_values79=np.array([-1,1j,1,-1j,-1,1][:period79])
    initial_vector79=np.zeros(period79,complex)
    initial_vector79[:3]=init78[:,0]
    middle_matrix79=np.kron(np.eye(period79//3),mid78)
    for insertion79 in range(4):
        columns79=[]
        for exponent79 in range(16):
            high79,low79=divmod(exponent79,8)
            prefix79=low79%(1<<insertion79)
            state79=np.roll(initial_vector79,prefix79)*early_values79
            state79=np.roll(state79,low79-prefix79)
            columns79.append(np.roll(late_values79*(middle_matrix79@state79),8*high79))
        columns79=np.array(columns79)
        reference79=abs(np.fft.fft(columns79,axis=0).T)**2/256
        for cover79 in ("auto","cycle","dual"):
            model79=EarlierPhaseProgressions(period79,3,4,3,init78,mid78,
                lambda j:late_values79[j],early_split=insertion79,
                early_phase=lambda j:early_values79[j],cover=cover79)
            got79=np.array([[model79.forced_joint(j,y)["joint_probability"]
                            for y in range(16)] for j in range(period79)])
            assert np.all(np.isfinite(got79)) and np.min(got79)>=0
            assert abs(got79.sum()-1)<2e-14 and np.max(abs(got79-reference79))<5e-16
            for j79 in range(period79):
                row79=model79.row(j79)
                expanded79={}
                for start79,count79,amp79 in row79["components"]:
                    for k79 in range(count79):
                        exponent79=start79+row79["stride"]*k79
                        assert exponent79 not in expanded79
                        expanded79[exponent79]=amp79
                assert all(abs(expanded79.get(e79,0j)-columns79[e79,j79])<4e-16
                           for e79 in range(16))
                assert row79["counters"]["row_local_terms"]+9<=model79.local_term_bound
                assert row79["counters"]["phase_queries"]==(
                    row79["counters"]["early_phase_queries"]+row79["counters"]["late_phase_queries"])
            draw79=model79.sample(np.random.default_rng(7900+joint_cases79),max_attempts=512)
            assert draw79["counters"]["work_draws"]==1
            assert draw79["counters"]["fourier_component_terms"]==draw79["attempts"]*draw79["component_count"]
            assert draw79["counters"]["row_local_terms"]+9<=model79.local_term_bound
            assert draw79["counters"]["phase_queries"]<=model79.stats()["phase_queries_per_draw_bound"]
            joint_cases79+=1
t("both covers and auto agree with independently normalized literal FFT joint laws",joint_cases79==24)

H79=np.array([[1,1],[1,-1]],complex)/np.sqrt(2)
revival79=EarlierPhaseProgressions(6,2,2,2,H79,H79,lambda j:1,
    early_split=1,early_phase=lambda j:(-1)**j)
old_revival79=LateWorkProgressions(6,2,2,2,H79,H79,lambda j:1)
t("the implemented row retains a component revived by the early phase",
  len(old_revival79.row(1)["components"])==1
  and len(revival79.row(1)["components"])==2
  and abs(dict((s,a) for s,n,a in revival79.row(1)["components"])[0]-1)<3e-16)

for caps79 in ({"max_local_terms":1},{"max_components":1},{"max_payload_bytes":1}):
    with patch("lab.work_first.np.array",side_effect=AssertionError("copied before cap")):
        try:
            EarlierPhaseProgressions(3,3,4,3,init78,mid78,lambda j:1,
                early_split=2,early_phase=lambda j:1,**caps79)
            t("earlier-phase structural cap rejects before any matrix copy",False)
        except MemoryError:
            t("earlier-phase structural cap rejects before any matrix copy",True)
for kwargs79 in ({"early_split":True},{"early_split":4},{"early_split":-1},
                 {"early_phase":None},{"cover":"guess"},{"cover":False}):
    args79=dict(early_split=2,early_phase=lambda j:1)
    args79.update(kwargs79)
    try:
        EarlierPhaseProgressions(3,3,4,3,init78,mid78,lambda j:1,**args79)
        t("invalid earlier insertion/oracle/cover rejects",False)
    except ValueError:
        t("invalid earlier insertion/oracle/cover rejects",True)
for bad_phase79 in (0,2,float("nan"),complex(float("inf"),0)):
    bad_model79=EarlierPhaseProgressions(3,3,4,3,init78,mid78,lambda j:1,
        early_split=1,early_phase=lambda j:bad_phase79)
    try:
        bad_model79.column(0)
        t("queried early phase must be finite and unit modulus",False)
    except ValueError:
        t("queried early phase must be finite and unit modulus",True)

wide79=EarlierPhaseProgressions(6,2,63,62,H79,H79,lambda j:1,
    early_split=61,early_phase=lambda j:(-1)**j)
wide_row79=wide79.row(0)
t("wide input creates bounded long progressions, never exponent/orbit arrays",
  wide79.cover=="dual" and len(wide_row79["components"])<=12
  and any(count79>10**12 for _,count79,_ in wide_row79["components"])
  and wide79.stats()["output_table_entries"]==wide79.stats()["expanded_progression_entries"]==0
  and wide79.stats()["numeric_payload_bound_bytes"]<16384)

# C79 actual one-attempt transition: r=3,b=1,Q=8,s=2,v=1 and parity f.
# Independently, phi_e=(-1)^e |e mod 3>. Cycle rows have m=(3,3,2)
# singleton components, stride 6, two fair interval bits and two gcd lifts.
# A seed/component/k/lift path has probability 1/(64*m); its acceptance is
# 64*joint[j,y]/m^2. Enumerate both branches where they have positive mass.
# This tests the real sampler, not only an algebraic forced-law evaluator.
decision79=EarlierPhaseProgressions(3,1,3,2,np.ones((1,1)),np.ones((1,1)),
    lambda j:1,early_split=1,early_phase=lambda j:(-1)**j,cover="cycle")
literal79=np.zeros((8,3),complex)
for e79 in range(8):
    literal79[e79,e79%3]=(-1)**e79
decision_reference79=abs(np.fft.fft(literal79,axis=0).T)**2/64
class DecisionRNG79:
    def __init__(self,e,c,m,k,lift,acceptance_draw):
        self.ints=iter((e,lift))
        self.reals=iter((.5,(c+.5)/m,((k&1)+.5)/2,
                         (((k>>1)&1)+.5)/2,acceptance_draw))
        self.integer_bounds=[]
        self.real_calls=0
    def integers(self,high):
        self.integer_bounds.append(high)
        value=next(self.ints)
        assert 0<=value<high
        return value
    def random(self):
        self.real_calls+=1
        value=next(self.reals)
        assert 0<=value<1
        return value
accepted79=np.zeros((3,8))
rejected79=0.
traversed79=0
rejection_calls79=0
for e79 in range(8):
    j79=e79%3
    m79=(3,3,2)[j79]
    for c79 in range(m79):
        for k79 in range(4):
            for lift79 in range(2):
                y79=(3*k79)%4+4*lift79
                a79=64*decision_reference79[j79,y79]/m79**2
                assert 0<=a79<=1+1e-15
                a79=min(1.,a79)
                traversed79+=1
                accepted79[j79,y79]+=a79/(64*m79)
                rejected79+=(1-a79)/(64*m79)
                for accept_branch79 in (True,False):
                    if (accept_branch79 and a79==0) or (not accept_branch79 and a79>1-1e-14):
                        continue
                    rng79=DecisionRNG79(e79,c79,m79,k79,lift79,
                        a79/2 if accept_branch79 else (1+a79)/2)
                    try:
                        result79=decision79.sample(rng79,max_attempts=1)
                        assert accept_branch79
                        assert (result79["work"],result79["output"])==(j79,y79)
                        assert result79["counters"]["work_draws"]==1
                    except RuntimeError:
                        assert not accept_branch79
                        rejection_calls79+=1
                    assert rng79.integer_bounds==[8,2] and rng79.real_calls==5
work_weights79=np.array([3,3,2])/8
accepted_work79=accepted79.sum(axis=1)
restored79=accepted79/accepted_work79[:,None]*work_weights79[:,None]
t("earlier stride: complete weighted RNG law and both gcd lifts match an independent FFT",
  traversed79==176 and rejection_calls79>0
  and np.max(abs(restored79-decision_reference79))<3e-16
  and np.max(abs(accepted_work79-work_weights79/np.array([3,3,2])))<3e-16
  and abs(accepted79.sum()+rejected79-1)<2e-14)
t("control: restarting work after rejection changes its marginal",
  np.sum(abs(accepted_work79/accepted_work79.sum()-work_weights79))/2>0.08)

print("[component weighting] opt-in square-root proposal, same accepted law")
weighted_cases79=0
for cls79 in (LateWorkProgressions,EarlierPhaseProgressions):
    for r79 in (3,6):
        options79={} if cls79 is LateWorkProgressions else dict(
            early_split=1,early_phase=lambda j:(1,1j,-1)[j%3])
        weighted_model79=cls79(r79,3,4,3,init78,mid78,
            lambda j:(1,-1j,-1)[j%3],**options79)
        # Independent literal matrices, not forced_joint as the reference.
        input79=np.zeros(r79,complex)
        input79[:3]=init78[:,0]
        full_middle79=np.kron(np.eye(r79//3),mid78)
        literal79=[]
        for e79 in range(16):
            low79,high79=e79%8,e79//8
            prefix79=low79%2 if options79 else 0
            state79=np.roll(input79,prefix79)
            if options79:
                state79*=np.array([(1,1j,-1)[j%3] for j in range(r79)])
            state79=np.roll(state79,low79-prefix79)
            state79=full_middle79@state79
            state79*=np.array([(1,-1j,-1)[j%3] for j in range(r79)])
            literal79.append(np.roll(state79,8*high79))
        target79=abs(np.fft.fft(np.array(literal79),axis=0).T)**2/256
        for mode79 in ("mass","root_mass"):
            joint79=np.zeros((r79,16))
            for j79 in range(r79):
                row79=weighted_model79.row(j79)
                counters79=weighted_model79._new_counters()
                prepared79=weighted_model79._prepare_proposal(row79,mode79,counters79)
                m79=len(row79["components"])
                total_proposal79=total_accepted79=0.
                for y79 in range(16):
                    law79=weighted_model79._fourier(row79,y79,counters79,prepared=prepared79)
                    joint79[j79,y79]=row79["norm"]/16*(law79["conditional_probability"] or 0.)
                    total_proposal79+=law79["proposal_probability"]
                    total_accepted79+=law79["proposal_probability"]*law79["acceptance"]
                    if m79:
                        assert abs(law79["proposal_probability"]*law79["acceptance"]
                                   *prepared79["expected_attempts"]-law79["conditional_probability"])<5e-16
                if m79:
                    assert abs(total_proposal79-1)<2e-14
                    assert abs(total_accepted79*prepared79["expected_attempts"]-1)<2e-14
                    assert 1-1e-14<=prepared79["expected_attempts"]<=m79*(1+1e-14)
                assert counters79["component_weight_terms"]==m79
                assert counters79["component_sqrt_terms"]==(m79 if mode79=="root_mass" else 0)
                assert counters79["weighted_ratio_divisions"]==(16*m79 if mode79=="root_mass" else 0)
            assert np.max(abs(joint79-target79))<5e-16 and abs(joint79.sum()-1)<2e-14
            draw79=weighted_model79.sample(np.random.default_rng(8000+weighted_cases79),
                max_attempts=512,proposal=mode79)
            dc79=draw79["counters"]
            assert dc79["work_draws"]==1 and dc79["component_weight_terms"]==draw79["component_count"]
            assert dc79["component_sqrt_terms"]==(draw79["component_count"] if mode79=="root_mass" else 0)
            assert dc79["weighted_ratio_divisions"]==(dc79["fourier_component_terms"] if mode79=="root_mass" else 0)
            weighted_cases79+=1
t("C78/C79 mass and root-mass accepted full joint laws match independent FFTs",weighted_cases79==8)
for mode79 in (True,None,1,[],"root"):
    with patch.object(weighted_model79,"_column",side_effect=AssertionError("RNG/column before mode check")):
        try:
            weighted_model79.sample(None,proposal=mode79)
            t("invalid proposal mode rejects before any RNG/column work",False)
        except ValueError:
            t("invalid proposal mode rejects before any RNG/column work",True)
t("root-mass forced query keeps an empty work row undefined",
  zero78.forced_joint(3,0,proposal="root_mass")["conditional_probability"] is None
  and zero78.forced_joint(3,0,proposal="root_mass")["proposal_probability"]==0)

# Isolated conditional-row test of the SAME loop used by sample, not a
# claim that this injected row matches the model's physical work columns.
row80=dict(work=0,stride=2,norm=1.,components=((0,1,np.sqrt(.9)),(2,2,np.sqrt(.05))))
model80=LateWorkProgressions(3,3,3,2,np.eye(3),np.eye(3),lambda j:1)
vector80=np.zeros(8,complex)
vector80[0]=np.sqrt(.9)
vector80[[2,4]]=np.sqrt(.05)
target80=abs(np.fft.fft(vector80))**2/8
class RowRNG80:
    def __init__(self,values,lift):
        self.values=iter(values)
        self.lift=lift
        self.real_calls=0
        self.integer_bounds=[]
    def random(self):
        self.real_calls+=1
        value=next(self.values)
        assert 0<=value<1
        return value
    def integers(self,high):
        self.integer_bounds.append(high)
        assert 0<=self.lift<high
        return self.lift
for mode80,weights80,e80 in (("mass",(.9,.1),2.),("root_mass",(.75,.25),1.6)):
    accepted80=np.zeros(8)
    rejected80=0.
    for c80 in range(2):
        selection80=weights80[0]/2 if c80==0 else weights80[0]+weights80[1]/2
        reduced80=(.25,.25,.25,.25) if c80==0 else (.5,.25,0.,.25)
        for k80,pred80 in enumerate(reduced80):
            if pred80==0:  # algebraic zero, not a float cutoff
                continue
            bits80=[.25 if k80%2==0 else .75,
                    .5 if c80==1 and k80==0 else (.25 if k80<2 else .75)]
            for lift80 in range(2):
                y80=k80+4*lift80
                q80=weights80[0]/8+weights80[1]*(.5,.25,0.,.25)[k80]/2
                a80=target80[y80]/(e80*q80)
                assert 0<a80<1
                probability80=weights80[c80]*pred80/2
                accepted80[y80]+=probability80*a80
                rejected80+=probability80*(1-a80)
                for accept80 in (True,False):
                    rng80=RowRNG80([selection80,*bits80,a80/2 if accept80 else (1+a80)/2],lift80)
                    counters80=model80._new_counters()
                    try:
                        result80=model80._sample_row(row80,rng80,counters80,
                            max_attempts=1,proposal=mode80)
                        assert accept80 and result80["output"]==y80
                        assert abs(result80["expected_attempts"]-e80)<1e-14
                    except RuntimeError:
                        assert not accept80
                    assert rng80.real_calls==4 and rng80.integer_bounds==[2]
                    assert counters80["work_draws"]==0  # conditional hook only
                    assert counters80["fourier_component_terms"]==2
                    assert counters80["weighted_ratio_divisions"]==(2 if mode80=="root_mass" else 0)
    assert np.max(abs(e80*accepted80-target80))<2e-16
    assert abs(accepted80.sum()+rejected80-1)<2e-15
t("isolated row: actual weighted component/bit/lift/accept-reject paths give the FFT law",True)

print("[nested phases] C81 rows and row-dependent progression strides")
from lab.work_first import NestedPhaseProgressions

# Exact b=1 columns: phases at controls 0,1,2 act on the corresponding
# binary prefix modulo the supplied orbit period. No production row formula
# is used to form this literal reference.
phase81a=lambda j:(1j,-1,1)[j%3]
phase81b=lambda j:(1,-1j,-1)[j%3]
literal81=np.zeros((8,3),complex)
for e81 in range(8):
    literal81[e81,e81%3]=phase81a(0)*phase81b(e81%2)*phase81a((e81%4)%3)
target81=abs(np.fft.fft(literal81,axis=0).T)**2/64
for cache81 in (False,True):
    for cut81 in ("auto",0,1,2):
        model81=NestedPhaseProgressions(3,1,3,2,np.ones((1,1)),np.ones((1,1)),
            lambda j:1,early_phases=((0,phase81a),(1,phase81b),(2,phase81a)),
            cut=cut81,cache_right=cache81)
        for e81 in range(8):
            assert model81.column(e81)["amplitudes"]=={e81%3:literal81[e81,e81%3]}
        for j81 in range(3):
            row81=model81.row(j81)
            counters81=row81["counters"]
            assert (counters81["row_local_terms"]+counters81["selector_pair_terms"]
                    +counters81["selector_residue_terms"]+counters81["selector_candidate_terms"]
                    +counters81["phase_partition_terms"]+1<=model81.local_term_bound)
            assert counters81["phase_queries"]<=model81.stats()["phase_queries_per_draw_bound"]
            assert counters81["row_phase_products"]<=model81.stats()["row_phase_products_bound"]
            if row81["selection"] is not None:
                assert row81["selection"]["coefficients"]==counters81["coefficient_pair_visits"]
                assert row81["selection"]["geometric"]==counters81["geometric_components"]
            for mode81 in ("mass","root_mass"):
                for y81 in range(8):
                    forced81=model81.forced_joint(j81,y81,proposal=mode81)
                    assert abs(forced81["joint_probability"]-target81[j81,y81])<3e-15
t("multiple distinct phase oracles: all cuts/cache modes match exact-root literal columns and FFT",True)

# T changes with work; the exact selector uses cut3 for j0,1 and cut4 for
# j2,3,4. The associated strides are 5 and40, hence different gcd lifts.
variable81=NestedPhaseProgressions(5,1,5,4,np.ones((1,1)),np.ones((1,1)),
    lambda j:1,early_phases=((0,lambda j:1),(3,lambda j:(-1)**j)))
assert [variable81.row(j)["cut"] for j in range(5)]==[3,3,4,4,4]
for mode81 in ("mass","root_mass"):
    for seed81 in range(16):
        draw81=variable81.sample(np.random.default_rng(81000+seed81),
                               max_attempts=256,proposal=mode81)
        gcd81=1 if draw81["work"]<2 else 8
        expected_queries81=2*(5-(gcd81.bit_length()-1))*draw81["attempts"]
        assert draw81["counters"]["marginal_queries"]==expected_queries81
        assert draw81["counters"]["work_draws"]==1
t("actual shared sampler uses each selected row stride for its progression queries",True)

revived81=NestedPhaseProgressions(6,2,2,2,H79,H79,lambda j:1,
    early_phases=((1,lambda j:(-1)**j),(2,lambda j:1j**(j%2))),cut=2)
row81=revived81.row(1)
assert abs(next(a for start,count,a in row81["components"] if start==0)-(1+1j)/2)<3e-16
t("two distinct phases revive a previously canceled candidate",True)

for caps81 in ({"max_local_terms":1},{"max_components":1},{"max_payload_bytes":1}):
    with patch("lab.work_first.np.array",side_effect=AssertionError("copied before cap")):
        try:
            NestedPhaseProgressions(6,2,4,3,H79,H79,lambda j:1,
                early_phases=((1,lambda j:1),(2,lambda j:1)),**caps81)
            t("nested preflight rejects before matrix copies",False)
        except MemoryError:
            t("nested preflight rejects before matrix copies",True)
for kwargs81 in ({"early_phases":((True,lambda j:1),)},
                 {"early_phases":((4,lambda j:1),)},
                 {"early_phases":((1,None),)},
                 {"early_phases":((1,lambda j:1),)*65},
                 {"cut":True},{"cut":4},{"cache_right":1}):
    try:
        NestedPhaseProgressions(6,2,4,3,H79,H79,lambda j:1,**kwargs81)
        t("invalid nested schedule/cut/cache mode rejects",False)
    except ValueError:
        t("invalid nested schedule/cut/cache mode rejects",True)
for value81 in (0,2,float("nan"),complex(float("inf"),0)):
    invalid81=NestedPhaseProgressions(3,1,3,2,np.ones((1,1)),np.ones((1,1)),
        lambda j:1,early_phases=((1,lambda j:value81),))
    try:
        invalid81.column(0)
        t("invalid queried nested phase rejects with retained counters",False)
    except ValueError:
        t("invalid queried nested phase rejects with retained counters",
          invalid81.last_counters["phase_queries"]==1)

t("existing sampling/default methods are inherited unchanged",
  NestedPhaseProgressions.sample is LateWorkProgressions.sample
  and NestedPhaseProgressions._sample_row is LateWorkProgressions._sample_row
  and NestedPhaseProgressions.forced_joint is LateWorkProgressions.forced_joint)

print("\nALL TESTS PASSED")
