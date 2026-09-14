"""Bounded edge audit for the exact-input b=3 verified finite-work extension.

The production API under test is the opt-in generalized
``VerifiedReflectionCircuit(..., block_size=3)`` plus its existing
``VerifiedFiniteWork`` contraction.  This probe does not implement another
propagator: small full-r products use the existing ``CoherentReflectionCircuit``
and ``direct_joint`` reference only.

P1: malformed dimensions, divisibility, scalar/rounded gate requests, and
    invalid masks are rejected; width-zero and insertion-end routes remain
    valid and forced conditional laws normalize.
P2: exact-zero work coordinates and a 2^-128 near-zero prefix are retained
    while complete forced finite-work laws normalize; any all-zero fallback
    occurrence is reported rather than silently omitted.
P3: the C56 r=6,b=3 early span-preserving Rx(pi/2) at s=1 is invisible,
    while the same genuinely non-scalar block at s=2 changes a conditional
    output ratio and the independent full-r law.
P4: widening a valid effect AFTER a selected prefix forces production
    refinement and replay; resetting that effect at the SAME node changes
    its conditional ratio;
    rectangular and norm enclosures agree outwardly on all six real
    coordinates and contain a higher-precision midpoint.

C1 is a must-fail API control: invalid/scalar requests must not be accepted.
C2 rejects resetting the selected backward effect to identity at the same
depth/output prefix. Earlier agent reports compared different depths and did
not trigger production replay; they are retained, not interpreted as this test.
All checks are deterministic, dense references are guarded before allocation,
and no histogram is used as evidence.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path

import numpy as np

from lab import Experiment


BLOCK = 3
WIDTH = 4
PERIOD = 9
MAX_BYTES = 16 << 20
TARGET = Fraction(1, 1024)


def report_path(prefix: str = "odd_block_edges") -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    path = Path("out") / f"{prefix}_{stamp}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def guard(shape, dtype=np.complex128, label="array"):
    payload = math.prod(int(x) for x in shape) * np.dtype(dtype).itemsize
    if payload > MAX_BYTES:
        raise MemoryError(f"{label} allocation {payload} exceeds 16 MiB")


def guard_entries(entries, bytes_per_entry=4096, label="retained entries"):
    if type(entries) is not int or entries < 0 or entries * bytes_per_entry > MAX_BYTES:
        raise MemoryError(f"{label} exceeds bounded diagnostic storage")


def json_safe(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def embedded_rx(theta):
    """Existing full-r diagnostic: a 2-level rotation plus fixed third label."""
    result = np.eye(BLOCK, dtype=complex)
    half = float(theta) / 2
    result[:2, :2] = np.array([[np.cos(half), -1j * np.sin(half)],
                               [-1j * np.sin(half), np.cos(half)]])
    return result


def embedded_rz(theta):
    result = np.eye(BLOCK, dtype=complex)
    half = float(theta) / 2
    result[0, 0] = np.exp(-1j * half)
    result[1, 1] = np.exp(1j * half)
    return result


def make_verified(V, *, period=PERIOD, width=WIDTH, backgrounds=None,
                  reflections=None, mode="rectangular"):
    return V(period, width, backgrounds or {}, reflections or {},
             block_size=BLOCK, enclosure_mode=mode)


def make_reference(CoherentReflectionCircuit, direct_joint, *, period=PERIOD,
                   width=WIDTH, backgrounds=None, reflections=None):
    backgrounds = backgrounds or {}
    matrices = {}
    for insertion, (axis, theta) in backgrounds.items():
        matrices[insertion] = embedded_rx(float(theta) * np.pi) if axis == "x" \
            else embedded_rz(float(theta) * np.pi)
    return direct_joint(CoherentReflectionCircuit(
        period, BLOCK, width, matrices,
        {s: (int(q), float(theta) * np.pi)
         for s, (q, theta) in (reflections or {}).items()}))


def forced_component_law(verified, VerifiedFiniteWork, mask, *, target=TARGET):
    sectors, outputs = verified.sectors, 1 << verified.width
    guard_entries(sectors * outputs + sectors + 16,
                  label="forced finite-work law entries")
    worker = VerifiedFiniteWork(verified, mask)
    law = {}
    conditional = []
    fallbacks = 0
    max_precision = 0
    for initial in range(sectors):
        mass = Fraction(0)
        final_labels = set()
        for output in range(outputs):
            result = worker.path(initial, output=output, target_tv=target)
            probability = Fraction(result["conditional_path_probability"])
            if probability < 0:
                raise AssertionError("negative forced probability")
            mass += probability
            final_labels.add(result["final_coarse_sector"])
            key = (result["final_coarse_sector"], output)
            law[key] = law.get(key, Fraction(0)) + probability / sectors
            fallbacks += int(result.get("zero_approximate_block_fallbacks", 0))
            max_precision = max(max_precision, int(result.get("max_working_precision", 0)))
        conditional.append(dict(initial_sector=initial, mass=str(mass),
                                final_sectors=sorted(final_labels)))
        if mass != 1 or len(final_labels) != 1:
            raise AssertionError("forced conditional law lost mass or route determinism")
    if sum(law.values(), Fraction(0)) != 1:
        raise AssertionError("forced joint finite-work law did not normalize")
    return law, dict(conditional=conditional, fallback_count=fallbacks,
                     max_working_precision=max_precision,
                     mass=str(sum(law.values(), Fraction(0))))


def law_tv(left, right):
    keys = set(left) | set(right)
    return sum(abs(left.get(k, Fraction(0)) - right.get(k, Fraction(0)))
               for k in keys) / 2


def run_invalid_controls(V, VerifiedFiniteWork, VerifiedScalarWork, valid):
    attempts = [
        ("period_not_divisible", lambda: make_verified(V, period=10)),
        ("block_not_three", lambda: V(9, WIDTH, {}, {}, block_size=2)),
        ("scalar_background", lambda: V(9, WIDTH, {1: 1}, {}, block_size=3)),
        ("rounded_angle", lambda: make_verified(
            V, backgrounds={1: ("x", 0.5)})),
        ("bad_enclosure_mode", lambda: V(9, WIDTH, {}, {}, block_size=3,
                                           enclosure_mode="bad")),
        ("mask_out_of_range", lambda: VerifiedFiniteWork(valid, 4)),
        ("scalar_b3_rejected", lambda: VerifiedScalarWork(valid, 0)),
    ]
    rejected = {}
    for name, call in attempts:
        try:
            call()
        except (ArithmeticError, TypeError, ValueError, OverflowError):
            rejected[name] = True
        else:
            rejected[name] = False
    return rejected


def width_zero_and_end_routes(V, VerifiedFiniteWork):
    zero = make_verified(V, period=6, width=0,
                         backgrounds={0: ("x", Fraction(1, 2))},
                         reflections={0: (1, Fraction(1, 3))})
    zero_rows = []
    for mask in (0, 1):
        worker = VerifiedFiniteWork(zero, mask)
        for initial in range(zero.sectors):
            result = worker.path(initial, output=0, target_tv=TARGET)
            zero_rows.append(dict(mask=mask, initial=initial,
                                  probability=str(result["conditional_path_probability"]),
                                  final=result["final_coarse_sector"]))
    end = make_verified(V, reflections={WIDTH: (1, Fraction(1, 3))})
    end_worker = VerifiedFiniteWork(end, 1)
    end_rows = []
    for initial in range(end.sectors):
        mass = Fraction(0)
        final = set()
        for output in range(1 << WIDTH):
            result = end_worker.path(initial, output=output, target_tv=TARGET)
            mass += Fraction(result["conditional_path_probability"])
            final.add(result["final_coarse_sector"])
        end_rows.append(dict(initial=initial, mass=str(mass), finals=sorted(final)))
    return dict(width_zero=zero_rows, insertion_end=end_rows,
                width_zero_normalized=all(Fraction(row["probability"]) == 1
                                          and row["final"] == ((-row["initial"]-1)%zero.sectors
                                              if row["mask"] else row["initial"])
                                          for row in zero_rows),
                insertion_end_normalized=all(Fraction(row["mass"]) == 1
                                             and row["finals"] == [(-row["initial"]-1)%end.sectors]
                                             for row in end_rows))


def enclosure_audit(verified_rect, verified_norm, binary_fraction):
    label = (1, 0, WIDTH, "reflection", WIDTH, 3)
    from flint import ctx
    with ctx.workprec(128):
        rect = verified_rect.prefix_enclosure(*label[:3], boundary=label[3],
                                              measured=label[4], output=label[5],
                                              working_precision=128)
        norm = verified_norm.prefix_enclosure(*label[:3], boundary=label[3],
                                              measured=label[4], output=label[5],
                                              working_precision=128)
        rect_data = tuple(tuple((binary_fraction(getattr(z, component).lower()),
                                 binary_fraction(getattr(z, component).upper()),
                                 binary_fraction(getattr(z, component).mid()),
                                 binary_fraction(getattr(z, component).rad().upper()))
                                for component in ("real", "imag")) for z in rect)
        norm_data = tuple(tuple((binary_fraction(getattr(z, component).lower()),
                                 binary_fraction(getattr(z, component).upper()),
                                 binary_fraction(getattr(z, component).mid()),
                                 binary_fraction(getattr(z, component).rad().upper()))
                                for component in ("real", "imag")) for z in norm)
    with ctx.workprec(256):
        high = verified_rect.prefix_enclosure(*label[:3], boundary=label[3],
                                              measured=label[4], output=label[5],
                                              working_precision=256)
        high_data = tuple(tuple(binary_fraction(getattr(z, component).mid())
                                for component in ("real", "imag")) for z in high)
    rows = []
    for index, (a, b, h) in enumerate(zip(rect_data, norm_data, high_data)):
        for component in ("real", "imag"):
            position = 0 if component == "real" else 1
            a_lo, a_hi, a_mid, a_rad = a[position]
            b_lo, b_hi, b_mid, b_rad = b[position]
            h_mid = h[position]
            rows.append(dict(coordinate=2 * index + (component == "imag"),
                             rect_contains_high=a_lo <= h_mid <= a_hi,
                             norm_contains_high=b_lo <= h_mid <= b_hi,
                             midpoint_overlap=abs(a_mid - b_mid) <= a_rad + b_rad,
                             rect_radius=str(a_rad), norm_radius=str(b_rad)))
    return dict(coordinates=rows,
                six_coordinate_containment=len(rows) == 6 and all(row["rect_contains_high"]
                                               and row["norm_contains_high"]
                                               for row in rows),
                six_coordinate_overlap=all(row["midpoint_overlap"] for row in rows))


def zero_near_audit(VerifiedReflectionCircuit, binary_fraction):
    """Use exact zero coordinates plus a 2^-128 nonzero W0 amplitude."""
    zero = make_verified(VerifiedReflectionCircuit, backgrounds={}, reflections={})
    zero_prefix = zero.prefix_dyadic(0, 0, 0, boundary="arithmetic",
                                     accuracy_bits=32)
    zero_coordinates = sum(re == 0 and im == 0 for re, im in zero_prefix.coordinates)
    near = make_verified(VerifiedReflectionCircuit,
                         backgrounds={0: ("x", Fraction(1, 1 << 128))},
                         reflections={})
    from flint import ctx
    with ctx.workprec(256):
        near_ball = near.prefix_enclosure(0, 0, 0, boundary="background",
                                          working_precision=256)
        total_norm = sum(value.real * value.real + value.imag * value.imag
                         for value in near_ball)
        tiny_norm = near_ball[1].real * near_ball[1].real \
            + near_ball[1].imag * near_ball[1].imag
        total_bounds = (binary_fraction(total_norm.lower()),
                        binary_fraction(total_norm.upper()))
        tiny_bounds = (binary_fraction(tiny_norm.lower()),
                       binary_fraction(tiny_norm.upper()))
    return dict(zero_coordinate_count=zero_coordinates,
                zero_coordinates=[json_safe(item) for item in zero_prefix.coordinates],
                near_total_norm_lower=str(total_bounds[0]),
                near_total_norm_upper=str(total_bounds[1]),
                near_tiny_component_lower=str(tiny_bounds[0]),
                near_tiny_component_upper=str(tiny_bounds[1]),
                near_tiny_component_positive=tiny_bounds[0] > 0,
                near_angle=str(Fraction(1, 1 << 128)))


def forced_output_law(verified, VerifiedFiniteWork, initial, *, target=TARGET):
    worker = VerifiedFiniteWork(verified, 0)
    result = []
    for output in range(1 << verified.width):
        row = worker.path(initial, output=output, target_tv=target)
        result.append(Fraction(row["conditional_path_probability"]))
    if sum(result, Fraction(0)) != 1:
        raise AssertionError("fixed-initial output law did not normalize")
    return result


def c56_span_reference(CoherentReflectionCircuit, direct_joint):
    period = 6
    guard((period, period), label="r=6 full-r matrix")
    guard((period, 1 << WIDTH), label="r=6 full-r output law")
    baseline = make_reference(CoherentReflectionCircuit, direct_joint, period=period,
                              backgrounds={}, reflections={})
    s1 = make_reference(CoherentReflectionCircuit, direct_joint, period=period,
                        backgrounds={1: ("x", Fraction(1, 2))}, reflections={})
    s2 = make_reference(CoherentReflectionCircuit, direct_joint, period=period,
                        backgrounds={2: ("x", Fraction(1, 2))}, reflections={})
    baseline_y, s1_y, s2_y = (np.sum(value, axis=0) for value in (baseline, s1, s2))
    return dict(s1_y_tv=float(np.abs(baseline_y - s1_y).sum() / 2),
                s2_y_tv=float(np.abs(baseline_y - s2_y).sum() / 2),
                baseline_mass=float(baseline.sum()), s1_mass=float(s1.sum()),
                s2_mass=float(s2.sum()))


def main() -> None:
    exp = Experiment("odd_block_edges", doc=__doc__, exit_on_fail=False)
    exp.predict("P1", "edge constructors and width-zero/end forced laws are valid and normalized")
    exp.predict("P2", "exact-zero coordinates and a 2^-128 prefix are retained while forced laws normalize")
    exp.predict("P3", "C56 s=1 cancellation persists while nonscalar s=2 changes a conditional ratio")
    exp.predict("P4", "production precision refinement replays the same chosen prefix and enclosure checks pass")
    exp.must_fail("C1", "invalid dimensions and scalar/rounded requests must be rejected")
    exp.must_fail("C2", "resetting the effect at the same nonempty prefix preserves its conditional law")
    started = time.perf_counter()
    rows = []
    try:
        from lab.coherent_routes import CoherentReflectionCircuit
        from lab.verified_finite_work import VerifiedFiniteWork, VerifiedScalarWork
        from lab.verified_prefix import VerifiedReflectionCircuit, binary_fraction
        from experiments.experiment_coherent_route_sampling import direct_joint

        exp.section("P1 API edge validation and route endpoints")
        valid = make_verified(VerifiedReflectionCircuit,
                              backgrounds={1: ("x", Fraction(1, 2)),
                                            3: ("z", Fraction(1, 5))},
                              reflections={2: (1, Fraction(1, 5)),
                                           3: (0, Fraction(1, 5))})
        rejected = run_invalid_controls(VerifiedReflectionCircuit, VerifiedFiniteWork,
                                        VerifiedScalarWork, valid)
        route_edges = width_zero_and_end_routes(VerifiedReflectionCircuit, VerifiedFiniteWork)
        p1 = all(rejected.values()) and route_edges["width_zero_normalized"] \
            and route_edges["insertion_end_normalized"]
        exp.check("P1", p1, f"rejected={rejected}, route_edges={route_edges}")
        rows.append(dict(series="api_edges", rejected=rejected, route_edges=route_edges))

        exp.section("P2 exact-zero and near-zero complete forced laws")
        zero = make_verified(VerifiedReflectionCircuit, backgrounds={}, reflections={})
        zero_law, zero_stats = forced_component_law(zero, VerifiedFiniteWork, 0)
        near = make_verified(VerifiedReflectionCircuit,
                             backgrounds={0: ("x", Fraction(1, 1 << 128))}, reflections={})
        near_law, near_stats = forced_component_law(near, VerifiedFiniteWork, 0)
        zero_near = zero_near_audit(VerifiedReflectionCircuit, binary_fraction)
        p2 = sum(zero_law.values(), Fraction(0)) == 1 \
            and sum(near_law.values(), Fraction(0)) == 1 \
            and zero_near["zero_coordinate_count"] >= 2 \
            and zero_near["near_tiny_component_positive"] \
            and Fraction(zero_near["near_tiny_component_upper"]) < Fraction(1,1 << 250) \
            and Fraction(zero_near["near_total_norm_lower"]) <= Fraction(1,1 << WIDTH) \
                <= Fraction(zero_near["near_total_norm_upper"])
        # Coarse mass-oracle request, separate from TARGET: exact half-half
        # root weights round to zero at p=0, exercising the declared fallback.
        from lab.verified_prefix import integer_cdf_counts
        coarse_cursor = VerifiedFiniteWork(zero,0).cursor(0,accuracy_bits=0)
        coarse_weights = coarse_cursor.weights()
        coarse_counts = integer_cdf_counts(coarse_weights,8)
        p2 &= not any(coarse_weights) and coarse_counts == (128,128)
        exp.check("P2", p2,
                  f"zero={zero_stats}, near={near_stats}, zero_near={zero_near}, "
                  f"zero_law_keys={len(zero_law)}, near_law_keys={len(near_law)}")
        rows.append(dict(series="forced_zero_near", zero=zero_stats, near=near_stats,
                         zero_near=zero_near,
                         coarse_root_weights=coarse_weights,coarse_root_counts=coarse_counts,
                         coarse_weight_bits=0,coarse_test_not_at_TARGET=True,
                         zero_law=json_safe(zero_law), near_law=json_safe(near_law)))

        exp.section("P3 C56 cancellation and genuinely nonscalar visibility")
        span = c56_span_reference(CoherentReflectionCircuit, direct_joint)
        identity = make_verified(VerifiedReflectionCircuit, period=6,
                                 backgrounds={}, reflections={})
        s2 = make_verified(VerifiedReflectionCircuit, period=6,
                           backgrounds={2: ("x", Fraction(1, 2))}, reflections={})
        identity_law, _ = forced_component_law(identity, VerifiedFiniteWork, 0)
        s2_law, _ = forced_component_law(s2, VerifiedFiniteWork, 0)
        identity_output = forced_output_law(identity, VerifiedFiniteWork, 0)
        s2_output = forced_output_law(s2, VerifiedFiniteWork, 0)
        conditional_tv = sum(abs(a - b) for a, b in zip(identity_output, s2_output)) / 2
        p3 = span["s1_y_tv"] < 2e-12 and span["s2_y_tv"] > 1e-10 \
            and all(abs(span[key]-1) < 2e-13 for key in ("baseline_mass","s1_mass","s2_mass")) \
            and conditional_tv > Fraction(1, 10000)
        exp.check("P3", p3, f"span_y_marginal={span}, "
                  f"verified_s2_output_conditional_tv={conditional_tv}")
        rows.append(dict(series="span_visibility", reference=span,
                         verified_s2_output_conditional_tv=str(conditional_tv)))

        exp.section("P4 precision replay and rectangular/norm enclosure audit")
        replay_worker = VerifiedFiniteWork(valid, 0)
        cursor = replay_worker.cursor(0, accuracy_bits=40, initial_precision=64)
        first_weights = cursor.weights()
        cursor.advance(0)  # valid selected prefix; the effect is now non-identity
        from flint import acb_mat,arb,ctx
        before_builds = cursor.builds
        with ctx.workprec(cursor.precision):
            cursor.effect += acb_mat([[arb(0,1) if i == j else 0
                                      for j in range(BLOCK)] for i in range(BLOCK)])
        second_weights = cursor.weights()
        fresh = replay_worker.cursor(0, accuracy_bits=40, initial_precision=cursor.precision)
        fresh.weights()
        fresh.advance(0)
        fresh_second = fresh.weights()
        wrong_reset = replay_worker.cursor(0, accuracy_bits=40, initial_precision=cursor.precision)
        wrong_reset.weights()
        wrong_reset.advance(0)
        wrong_reset.effect = acb_mat([[int(i == j) for j in range(BLOCK)] for i in range(BLOCK)])
        wrong_weights = wrong_reset.weights()
        normalized = lambda values: [Fraction(value, sum(values)) for value in values]
        correct_ratio = normalized(fresh_second)
        wrong_ratio = normalized(wrong_weights)
        reset_ratio_tv = sum(abs(a - b) for a, b in zip(correct_ratio, wrong_ratio)) / 2
        replay = dict(production_builds=cursor.builds,
                      production_replayed=cursor.replayed,
                      first_weights=json_safe(first_weights),
                      second_weights=json_safe(second_weights),
                      fresh_high_precision_second=json_safe(fresh_second),
                      fresh_replay_consistent=second_weights == fresh_second,
                      refined_after_selected_prefix=cursor.builds > before_builds and cursor.replayed >= 1,
                      same_chosen_prefix=cursor.measured == wrong_reset.measured == 1
                                         and cursor.output == wrong_reset.output == 0,
                      wrong_reset_same_node=json_safe(wrong_weights),
                      reset_ratio_tv=str(reset_ratio_tv))
        rect = make_verified(VerifiedReflectionCircuit,
                             backgrounds={1: ("x", Fraction(1, 2)),
                                           3: ("z", Fraction(1, 5))},
                             reflections={2: (1, Fraction(1, 5))}, mode="rectangular")
        norm = make_verified(VerifiedReflectionCircuit,
                             backgrounds={1: ("x", Fraction(1, 2)),
                                           3: ("z", Fraction(1, 5))},
                             reflections={2: (1, Fraction(1, 5))}, mode="norm")
        enclosure = enclosure_audit(rect, norm, binary_fraction)
        p4 = replay["fresh_replay_consistent"] and replay["refined_after_selected_prefix"] \
            and replay["same_chosen_prefix"] \
            and enclosure["six_coordinate_containment"] \
            and enclosure["six_coordinate_overlap"]
        exp.check("P4", p4, f"replay={replay}, enclosure={enclosure}")
        rows.append(dict(series="precision_and_enclosures", replay=replay,
                         enclosure=enclosure))

        exp.fail_check("C1", all(rejected.values()), f"rejected={rejected}")
        exp.fail_check("C2", reset_ratio_tv > Fraction(1,10000),
                       f"same-node reset conditional TV={float(reset_ratio_tv):.8g}")
        path = report_path()
        if not exp.finish(report_path=path, rows=json_safe(rows),
                          metadata=dict(period=PERIOD, block_size=BLOCK, width=WIDTH,
                                        target_tv=str(TARGET),
                                        reference="existing full-r direct_joint only",
                                        allocation_bound_bytes=MAX_BYTES,
                                        arithmetic="exact rational inputs; Arb/Acb enclosures")):
            raise AssertionError("experiment harness failed")
        print(json.dumps(dict(ok=True, report=str(path))))
        return
    except BaseException as exc:
        failure = report_path("odd_block_edges_failure")
        failure.write_text(json.dumps(dict(ok=False, error=repr(exc),
                                           traceback=traceback.format_exc()),
                                         indent=2) + "\n")
        raise


if __name__ == "__main__":
    main()
