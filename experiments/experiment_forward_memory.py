"""Does one approximate forward state per depth control the complete QFT law?

Extends C64 on the frozen r9,b3,t4 fixture, all deterministic histories and
initial sectors. Only epsilon varies: rho_tilde_i=(1-epsilon)rho_i+epsilon I/3.
Exact backward effects are retained. This is a capped diagnostic, NOT a new
production sampler or an implementation that saves checkpoint memory.

PREDICTIONS, WRITTEN BEFORE MEASURING.
  P1 At each depth the adaptive child effects form a complete POVM. The
     independent full-r component law agrees with its density/effect law.
  P2 The exact finite-bit output law normalizes and its outward ideal TV is
     bounded by sum_i sqrt(b)||rho_i-rho_tilde_i||_F plus
     2(2^t-1)eta+t/2^L. The first term upper-bounds trace norms; numerical
     errors are charged DIRECTLY against the true tree, because arbitrary
     rho_tilde_i need not be channel-consistent.
  P3 With epsilon=0 the experiment's kernels agree with the production
     finite-work cursor at matched p,L; nonzero approximation is observable.
  P4 The C65 binary late-gate null can hide this state sensitivity; it is
     not used as the main approximation fixture.
  C1 Resetting the backward effect does not preserve the frozen b3 output law.
  C2 A maximally mixed work marginal does not imply output independence:
     for r=b=3,t=2 and mixed initial work, the exact y law is
     (3/8,1/4,1/8,1/4), whereas resetting E gives the uniform law, TV=1/8.
  C3 In a separate five-dimensional classical POVM, prefix-dependent state
     choices of trace error 1/5 produce aggregate conditional TV 4/13.
     This refutes a general max-error extension, NOT a b3 QFT-specific bound.

No amplitude cutoff. Arb/Acb enclosures use P192/P256; all endpoints and CDF
probabilities are exact Fractions. The dense independent comparison is only
a floating diagnostic. Retained tree/row allowances are capped before work.
"""
from __future__ import annotations

from fractions import Fraction
from datetime import datetime, timezone
from pathlib import Path
import time
import traceback
import numpy as np

from lab import Experiment
from lab.semiclassical import _qft_branch_operators
from lab.verified_finite_work import VerifiedFiniteWork, _adjoint
from lab.verified_prefix import (VerifiedReflectionCircuit, _rational_phase,
    binary_fraction, round_dyadic, integer_cdf_counts)
from lab.verified_rejection import VerifiedRejectionSampler
from experiments.experiment_verified_odd_block import verified_fixture, float_fixture
from experiments.experiment_verified_finite_work import deterministic_float_component
from experiments.experiment_verified_rejection import (
    target_interval, tv_from_intervals, tv_lower_from_intervals)

EPSILONS = (Fraction(0), Fraction(1,4096), Fraction(1,256), Fraction(1,16))
BITS = 80
CDF_BITS = 80
PRECISIONS = (192,256)
MAX_BYTES = 32 << 20


def identity(b):
    from flint import acb_mat
    return acb_mat([[int(i == j) for j in range(b)] for i in range(b)])


def check_zero_matrix(A):
    return all(A[i,j].contains(0) for i in range(A.nrows()) for j in range(A.ncols()))


def frobenius_trace_upper(A):
    """Outward ||A||_1 <= sqrt(b)||A||_F, without an eigensolver."""
    from flint import arb
    # Squaring a ball containing zero then taking sqrt can straddle the
    # nonnegative domain. Use exact nonnegative upper endpoints throughout.
    total = sum((abs(A[i,j]).upper()**2 for i in range(A.nrows())
                 for j in range(A.ncols())), arb(0))
    return binary_fraction((A.nrows()*total).upper().sqrt().upper())


def finite_tree(branches, states, epsilon, *, reset=False):
    """Enumerate existing two-sided contractions; fixed precision, tiny only.

    Unlike the production cursor, this diagnostic retains a whole tiny tree
    and checks its enclosures before rounding. It never claims adaptive
    convergence or certifies resource-conditioned execution on larger inputs.
    """
    from flint import fmpq
    t,b = len(branches),states[0].nrows()
    if t > 4 or b > 3 or t < 1 or b < 2:
        raise ValueError("bounded diagnostic requires 1<=t<=4, 2<=b<=3")
    if (2*(1 << t)*b*b + 4*(t+1)*b*b)*512 > MAX_BYTES:
        raise MemoryError("tree scalar allowance exceeds 32 MiB")
    I = identity(b)
    e = fmpq(epsilon.numerator,epsilon.denominator)
    approx = [(1-e)*rho+e*I/b for rho in states[:-1]]
    errors = [frobenius_trace_upper(rho-r) for rho,r in zip(states,approx)]
    for rho in states + approx:
        if not (rho.trace()-1).contains(0) or not check_zero_matrix(rho-_adjoint(rho)):
            raise AssertionError("forward state trace/Hermiticity check failed")
    # PSD follows from unitary channels and convex mixtures, not interval
    # eigenvalue thresholding. Input matrices were constructed by shared code.
    frontier = {0:(Fraction(1),I)}
    completeness, aggregate_errors, radii = [],[],[]
    for depth in range(t):
        i = t-1-depth
        children, sum_effects = {}, I*0
        summed_mass_error = Fraction(0)
        max_radius = Fraction(0)
        for prefix,(mass,E) in frontier.items():
            phase = _rational_phase(-prefix,1 << (depth+1))
            effects = [_adjoint(K)*(I if reset else E)*K
                       for K in _qft_branch_operators(*branches[i],phase)]
            weights = []
            for F in effects:
                w = (F*approx[i]).trace()
                if not w.imag.contains(0) or w.real.upper() < 0:
                    raise AssertionError("invalid approximate child mass")
                rad = binary_fraction(w.real.rad().upper())
                max_radius = max(max_radius,rad)
                if rad > Fraction(1,1 << (BITS+1)):
                    raise ArithmeticError("fixed diagnostic precision insufficient")
                weights.append(max(0,round_dyadic(binary_fraction(w.real.mid()),BITS)))
                sum_effects += F
                summed_mass_error += binary_fraction(abs((F*(states[i]-approx[i])).trace()).upper())
            counts = integer_cdf_counts(tuple(weights),CDF_BITS)
            for bit in (0,1):
                children[prefix | (bit << depth)] = (
                    mass*Fraction(counts[bit],1 << CDF_BITS),effects[bit])
        if sum(mass for mass,_ in children.values()) != 1:
            raise AssertionError("finite tree lost mass")
        completeness.append(check_zero_matrix(sum_effects-I))
        aggregate_errors.append(summed_mass_error)
        radii.append(max_radius)
        frontier = children
    law = {y:mass for y,(mass,_E) in frontier.items()}
    numerical = Fraction(2*((1 << t)-1),1 << BITS)+Fraction(t,1 << CDF_BITS)
    # Numerical widths in the norm estimates can slightly exceed zero; they
    # are retained, not silently chopped. Compare to a separately rounded-up
    # epsilon*||rho-I/b|| bound for the POVM audit.
    return law, dict(state_error_upper=errors, numerical_tv_upper=numerical,
        total_tv_upper=min(Fraction(1),sum(errors)+numerical),
        complete_povm=all(completeness), max_weight_radius=max(radii),
        aggregate_mass_error_upper=aggregate_errors)


def production_law(worker,initial):
    law = {}
    for y in range(1 << worker.circuit.width):
        cursor = worker.cursor(initial,accuracy_bits=BITS,initial_precision=256)
        mass = Fraction(1)
        for j in range(worker.circuit.width):
            bit = (y >> j)&1
            mass *= Fraction(integer_cdf_counts(cursor.weights(),CDF_BITS)[bit],1 << CDF_BITS)
            cursor.advance(bit)
        law[y] = mass
    return law


def fraction_tv(a,b):
    return sum(abs(a.get(k,0)-b.get(k,0)) for k in a.keys() | b.keys())/2


def main():
    from flint import ctx
    exp = Experiment("forward_memory",doc=__doc__,exit_on_fail=False)
    exp.predict("P1","adaptive effects are complete and independent references agree")
    exp.predict("P2","complete finite laws meet state-plus-numerical budgets")
    exp.predict("P3","zero epsilon matches production; nonzero epsilon is observable")
    exp.predict("P4","the binary scalar-tail null hides the chosen PSD-mixture sensitivity")
    exp.must_fail("C1","resetting effects preserves the original pure-input output law")
    exp.must_fail("C2","maximally mixed work implies uncorrelated measured outputs")
    exp.must_fail("C3","prefix-dependent states obey one max trace-norm budget for every POVM")
    rows,controls = [],{}
    start = time.monotonic()
    error = None
    try:
        # Includes all retained row scalars and temporary full-r diagnostics;
        # conservative allowance, not RSS or native allocation measurement.
        if len(EPSILONS)*4*3*len(PRECISIONS)*128*1024 > MAX_BYTES:
            raise MemoryError("retained report allowance exceeds cap")
        c,comparator = verified_fixture(),float_fixture()
        sampler = VerifiedRejectionSampler(c)
        independent_error = 0.0
        for mask in range(4):
            component = sampler.component(mask)
            independent = deterministic_float_component(comparator,c,mask)
            if abs(float(independent.sum())-1) > 2e-14:
                raise AssertionError("independent full-r law lost mass")
            worker = VerifiedFiniteWork(c,mask)
            for initial in range(3):
                reference_law = production_law(worker,initial)
                for P in PRECISIONS:
                    with ctx.workprec(P):
                        branches,states,final = worker._build(initial)
                        I = identity(3)
                        if not all(check_zero_matrix(_adjoint(B)*B-I)
                                   for pair in branches for B in pair):
                            raise AssertionError("branch unitarity failed")
                        target = {y:tuple(3*x for x in target_interval(component,final,y,P))
                                  for y in range(16)}
                        if not sum(lo for lo,_ in target.values()) <= 1 <= sum(hi for _,hi in target.values()):
                            raise AssertionError("ideal intervals lost mass")
                        independent_error = max(independent_error,max(
                            abs(float((lo+hi)/6)-independent[final,y])
                            for y,(lo,hi) in target.items()))
                        for epsilon in EPSILONS:
                            law,stats = finite_tree(branches,states,epsilon)
                            tv = tv_from_intervals(law,target)
                            rows.append(dict(mask=mask,initial=initial,precision=P,
                                epsilon=str(epsilon),law=[str(law[y]) for y in range(16)],
                                mass=str(sum(law.values())),tv_upper=str(tv),
                                tv_lower=str(tv_lower_from_intervals(law,target)),
                                **{k:[str(x) for x in v] if isinstance(v,list) else str(v)
                                   for k,v in stats.items()},
                                budget_pass=tv <= stats["total_tv_upper"],
                                povm_pass=stats["complete_povm"],
                                aggregate_pass=epsilon == 0 or all(
                                    v <= bound for v,bound in zip(
                                        stats["aggregate_mass_error_upper"],
                                        reversed(stats["state_error_upper"]))),
                                zero_matches=epsilon != 0 or law == reference_law))
                        if mask == 0 and initial == 0 and P == 256:
                            wrong,_ = finite_tree(branches,states,Fraction(0),reset=True)
                            controls["reset_original_tv_lower"] = str(tv_lower_from_intervals(wrong,target))
        with ctx.workprec(256):
            mixed = VerifiedFiniteWork(VerifiedReflectionCircuit(3,2,{}, {},block_size=3),0)
            branches,_states,_ = mixed._build(0)
            states = [identity(3)/3 for _ in range(3)]
            law,stats = finite_tree(branches,states,Fraction(0))
            wrong,_ = finite_tree(branches,states,Fraction(0),reset=True)
            exact = dict(enumerate((Fraction(3,8),Fraction(1,4),Fraction(1,8),Fraction(1,4))))
            controls.update(mixed_correct_tv=str(fraction_tv(law,exact)),
                mixed_reset_tv=str(fraction_tv(law,wrong)),
                mixed_reset_uniform=all(p == Fraction(1,4) for p in wrong.values()),
                mixed_povm=stats["complete_povm"],mixed_numerical=str(stats["numerical_tv_upper"]))
            binary = VerifiedFiniteWork(VerifiedReflectionCircuit(10,4,
                {1:("x",Fraction(1,7)),3:("z",Fraction(1,5))},{}),0)
            binary_errors = []
            for initial in range(5):
                branches,states,_ = binary._build(initial)
                exact_binary,stats = finite_tree(branches,states,Fraction(0))
                mixed_binary,_ = finite_tree(branches,states,Fraction(1,16))
                binary_errors.append(fraction_tv(exact_binary,mixed_binary))
            controls.update(binary_null_tv_max=str(max(binary_errors)),
                            binary_null_numerical=str(stats["numerical_tv_upper"]))
        # Exact diagonal five-state POVM, independently checked as vectors.
        # Four parents a=1..4 have F_a0=|0><0|/4, F_a1=|a><a|.
        # rho=|0><0|, rho_tilde_a=.9|0><0|+.1|a><a|.
        povm_total = [Fraction(0)]*5
        aggregate = Fraction(0)
        for a in range(1,5):
            F0 = [Fraction(1,4)]+[Fraction(0)]*4
            F1 = [Fraction(int(j == a)) for j in range(5)]
            rho = [Fraction(1)]+[Fraction(0)]*4
            approx = [Fraction(9,10)]+[Fraction(int(j == a),10) for j in range(1,5)]
            if sum(approx) != 1 or min(approx) < 0:
                raise AssertionError("prefix-dependent state invalid")
            distance = sum(abs(x-y) for x,y in zip(rho,approx))
            if distance != Fraction(1,5):
                raise AssertionError("counterexample trace norm wrong")
            w = [sum(x*y for x,y in zip(F,rho)) for F in (F0,F1)]
            v = [sum(x*y for x,y in zip(F,approx)) for F in (F0,F1)]
            aggregate += sum(w)*sum(abs(x/sum(w)-y/sum(v)) for x,y in zip(w,v))/2
            povm_total = [x+y+z for x,y,z in zip(povm_total,F0,F1)]
        if povm_total != [1]*5 or aggregate != Fraction(4,13):
            raise AssertionError("prefix-dependent POVM computation failed")
        controls.update(prefix_dependent_aggregate_tv=str(aggregate),
            prefix_dependent_max_trace_error="1/5",prefix_dependent_dimension=5,
            prefix_dependent_is_actual_qft_fixture=False)
        exp.check("P1",all(row["povm_pass"] for row in rows) and independent_error < 2e-14,
                  f"all-depth completeness; independent max error={independent_error:.3g}")
        exp.check("P2",all(row["budget_pass"] and row["aggregate_pass"]
                           and row["mass"] == "1" for row in rows),
                  f"{len(rows)} complete laws at two precisions")
        observable = max(Fraction(row["tv_lower"]) for row in rows if row["epsilon"] == "1/16")
        exp.check("P3",all(row["zero_matches"] for row in rows) and observable > Fraction(1,10000),
                  f"matched exact finite kernels; nonzero epsilon TV lower={float(observable):.6g}")
        exp.check("P4",Fraction(controls["binary_null_tv_max"]) <= 2*Fraction(controls["binary_null_numerical"]),
                  f"binary null complete-law change={controls['binary_null_tv_max']}")
        exp.fail_check("C1",Fraction(controls["reset_original_tv_lower"]) > Fraction(1,10000),
                       f"outward reset TV lower={float(Fraction(controls['reset_original_tv_lower'])):.6g}")
        exp.fail_check("C2",Fraction(controls["mixed_reset_tv"]) > Fraction(1,10000),
                       f"mixed marginal reset TV={controls['mixed_reset_tv']}")
        exp.fail_check("C3",aggregate > Fraction(1,5),
                       f"general POVM aggregate TV={aggregate} > max state error=1/5")
        if not (controls["mixed_povm"] and controls["mixed_reset_uniform"]
                and Fraction(controls["mixed_correct_tv"]) <= Fraction(controls["mixed_numerical"])):
            raise AssertionError("mixed-input exact reference check failed")
    except Exception:
        error = traceback.format_exc()
        exp.check("P2",False,error)
    path = Path("out")/f"forward_memory_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}.json"
    ok = exp.finish(report_path=path,rows=rows,metadata=dict(controls=controls,
        elapsed_seconds=time.monotonic()-start,error=error,bits=BITS,cdf_bits=CDF_BITS,
        precision_levels=PRECISIONS,allocation_allowance_bytes=MAX_BYTES,
        approximation_is_experiment_only=True,production_changes=False))
    print(path,flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
