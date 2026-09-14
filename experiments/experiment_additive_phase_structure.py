"""Exact dense sector coupling and the limits of an exponent-phase shortcut.

DERIVED BEFORE MEASURING (TODO 35). For prime N, r=N-1=bM and k!=0,
each subgroup kernel coefficient is nonzero: over Q(zeta_M), zeta_N still
has minimal polynomial Phi_N, whereas the putative vanishing sum has no
constant term and degree at most N-1. We independently check the sums by
integer-polynomial reduction modulo Phi_lcm(N,M), not a float zero cutoff.

P1: all subgroup coefficients are nonzero for primes (7,3),(13,2),(17,3),
    each divisor b<=6 of r, all p/delta, and k=1,2.
P2: the b-term full Gauss-sum decomposition matches the direct kernel;
    each one-sector output probability <=N/M^2, so its best L-sector
    truncation loses at least max(0,1-L*N/M^2).
P3: the known no-kick reachable-prefix vectors satisfy scalar phase transfer.
C1: composite N=9,a=2,b=2,r=6,k=1 has exact zeros, refuting extension of
    prime-field nonvanishing to arbitrary modular arithmetic.
C2: the frozen N13 b3 initial/post-control mixer prefix is not mapped to
    scalar multiples by interior G1 or G2; C53's particular transfer fails.

Exact checks cap cyclotomic order at 272, polynomial coefficient input slots
at 100k, and reduction calls at 500; floating direct sums cap 50k terms and
all retained numerical arrays at 16MiB. Exact polynomial costs are counted
as calls/input slots, NOT bit-operation/runtime/RSS measurements. Existing
full-r reference matrices are used, not a new propagator or sampler.
"""
from __future__ import annotations

import cmath
from datetime import datetime, timezone
import math
from pathlib import Path
import traceback

import numpy as np
from flint import fmpz_poly

from lab import Experiment

MAX_ROOT_ORDER = 272
MAX_POLY_SLOTS = 100_000
MAX_REDUCTIONS = 500
MAX_SUM_TERMS = 50_000
MAX_BYTES = 16*1024*1024
TOL = 3e-11


def phase(numerator, denominator):
    return cmath.exp(2j*math.pi*(numerator % denominator)/denominator)


def cases():
    return [(N,a,b,N-1,k) for N,a in ((7,3),(13,2),(17,3))
            for b in range(1,7) if (N-1) % b == 0 for k in (1,2)]


def preflight():
    rows = cases()+[(9,2,2,6,1)]
    slots = sum(r*math.lcm(N,r//b) for N,a,b,r,k in rows)
    reductions = sum(r for N,a,b,r,k in rows)
    sum_terms = sum(r*((r//b)+b*r) for N,a,b,r,k in cases())
    root_orders = {math.lcm(N,r//b) for N,a,b,r,k in rows}
    # Polynomial coefficient payload is separate from numerical ndarray
    # payload. Cache has one tiny Phi per root order and no kernel history.
    if (max(root_orders)>MAX_ROOT_ORDER or slots>MAX_POLY_SLOTS
            or reductions>MAX_REDUCTIONS or sum_terms>MAX_SUM_TERMS):
        raise MemoryError("frozen structure preflight exceeds budget")
    numeric_bytes = 16*100*16*16
    if numeric_bytes>MAX_BYTES:
        raise MemoryError("combined numerical payload budget")
    return dict(polynomial_input_slots=slots, reductions=reductions,
                float_summands=sum_terms, root_orders=sorted(root_orders),
                numeric_payload_bound_bytes=numeric_bytes,
                cached_polynomial_degree_sum=sum(int(fmpz_poly.cyclotomic(L).degree()) for L in root_orders))


def exact_numerator(N,a,b,r,k,p,delta,cache,cost):
    M = r//b
    L = math.lcm(N,M)
    if (L>MAX_ROOT_ORDER or cost["polynomial_input_slots"]+L>MAX_POLY_SLOTS
            or cost["reductions"]+1>MAX_REDUCTIONS):
        raise MemoryError("exact polynomial per-call cap")
    cost["polynomial_input_slots"] += L
    cost["reductions"] += 1
    coeff = [0]*L
    for m in range(M):
        x = pow(a,b*m+p,N)
        power = ((L//N)*k*x+(L//M)*delta*m)%L
        coeff[power] += 1
    if L not in cache:
        cache[L] = fmpz_poly.cyclotomic(L)
    remainder = fmpz_poly(coeff) % cache[L]
    return remainder


def kernels(N,a,b,r,k,cost):
    M = r//b
    count = r*(M+b*r)
    if cost["float_summands"]+count>MAX_SUM_TERMS:
        raise MemoryError("kernel sum preflight")
    # Every direct and reconstructed term is counted inside its actual loop.
    direct = np.zeros((b,M),complex)
    reconstructed = np.zeros_like(direct)
    values = [pow(a,j,N) for j in range(r)]
    for p in range(b):
        scale = pow(a,p,N)
        for delta in range(M):
            for m in range(M):
                direct[p,delta] += phase(k*values[b*m+p],N)*phase(delta*m,M)/M
                cost["float_summands"] += 1
            for t in range(b):
                s = delta+t*M
                for j,x in enumerate(values):
                    reconstructed[p,delta] += phase(k*scale*x,N)*phase(s*j,r)/r
                    cost["float_summands"] += 1
    return direct,reconstructed


def prefix_transfer_rows():
    # Ordinary products of the existing finite branch matrices for four
    # fixed histories; not a new circuit parser or state propagator.
    from experiments.experiment_additive_phase_output import full_pairs, additive_phase
    initial,pairs = full_pairs(0)
    rows=[]
    for k in (0,1,2):
        gate=additive_phase(k)
        for ell in range(4):
            state=initial[:,0]
            for i in range(2):
                state=pairs[i][(ell>>i)&1]@state
            norm=float(np.vdot(state,state).real)
            target=gate@state
            scalar=np.vdot(state,target)/norm
            residual=float(np.linalg.norm(target-scalar*state)/np.sqrt(norm))
            if not math.isfinite(residual) or abs(norm-1)>TOL:
                raise ArithmeticError("invalid reachable prefix vector")
            rows.append(dict(k=k,low_history=ell,norm=norm,residual=residual,
                             best_scalar=[float(scalar.real),float(scalar.imag)]))
    return rows


def main():
    exp=Experiment("additive_phase_structure",doc=__doc__,exit_on_fail=False)
    exp.predict("P1","prime subgroup kernels have no exact zeros")
    exp.predict("P2","Gauss decomposition, one-sector mass and truncation bounds hold")
    exp.predict("P3","zero kick has the reachable scalar-transfer property")
    exp.must_fail("C1","composite N9 violates the prime nonzero conclusion")
    exp.must_fail("C2","nonzero interior kick has a scalar exponent-phase replacement")
    rows=[]
    cost=dict(polynomial_input_slots=0,reductions=0,float_summands=0)
    metadata={}
    try:
        plan=preflight()
        metadata["preflight"]=plan
        cache={}
        for N,a,b,r,k in cases():
            if len({pow(a,j,N) for j in range(r)})!=r or pow(a,r,N)!=1:
                raise ValueError("frozen primitive-root order verification failed")
            M=r//b
            nonzero=[]
            for p in range(b):
                nonzero.append(sum(bool(exact_numerator(N,a,b,r,k,p,d,cache,cost))
                                   for d in range(M)))
            direct,reconstructed=kernels(N,a,b,r,k,cost)
            error=float(np.max(abs(direct-reconstructed)))
            mass=np.abs(direct)**2
            cap=N/(M*M)
            # Pure p states plus two normalized coherent fine vectors.
            # Their sector weights are a mixture because output p is orthogonal.
            weights=[np.eye(b)[p] for p in range(b)]
            weights += [np.full(b,1/b),np.arange(1,b+1)/sum(range(1,b+1))]
            max_truncation_violation=0.
            for weights_p in weights:
                probs=weights_p@mass
                ordered=np.sort(probs)[::-1]
                for L in range(M+1):
                    retained=float(np.sum(ordered[:L]))
                    max_truncation_violation=max(max_truncation_violation,retained-min(1,L*cap))
            row=dict(N=N,a=a,b=b,r=r,M=M,k=k,exact_nonzero_per_p=nonzero,
                     min_numeric_magnitude=float(np.min(abs(direct))),
                     decomposition_error=error,probability_cap=cap,
                     max_sector_probability=float(np.max(mass)),
                     norm_error=float(np.max(abs(np.sum(mass,axis=1)-1))),
                     max_truncation_violation=max_truncation_violation)
            rows.append(row)
        exp.check("P1",all(row["exact_nonzero_per_p"]==[row["M"]]*row["b"] for row in rows),
                  f"all {sum(row['r'] for row in rows)} prime coefficients nonzero by exact polynomial reduction")
        exp.check("P2",all(math.isfinite(row["decomposition_error"])
                  and row["decomposition_error"]<TOL and row["norm_error"]<TOL
                  and row["max_sector_probability"]<=row["probability_cap"]+TOL
                  and row["max_truncation_violation"]<TOL for row in rows),
                  "independent character-sum decomposition and every retained-sector size")
        composite=[]
        for p in range(2):
            remainders=[exact_numerator(9,2,2,6,1,p,d,cache,cost) for d in range(3)]
            composite.append(dict(p=p,nonzero=[bool(v) for v in remainders],
                                  remainders=[[int(c) for c in v] for v in remainders]))
        exp.fail_check("C1",[v["nonzero"] for v in composite]==[[False,False,True],[False,True,False]],
                       "N9 b2: two EXACT zeros per p; only one routed sector remains")
        rows.append(dict(series="composite_control",values=composite))
        transfer=prefix_transfer_rows()
        exp.check("P3",all(row["residual"]<TOL for row in transfer if row["k"]==0),
                  "all four no-kick histories are scalar-transfer nulls")
        exp.fail_check("C2",all(row["residual"]>.1 for row in transfer if row["k"]!=0),
                       "each frozen nonzero kick/history defeats a scalar phase replacement")
        rows.append(dict(series="transfer_control",values=transfer))
        exp.check("P1",all(cost[key]==plan[key] for key in cost),"all exact and floating-loop counters reconcile")
    except Exception as exc:
        exp.check("execution",False,repr(exc))
        rows.append(dict(error=repr(exc),traceback=traceback.format_exc()))
    metadata.update(cost=cost,limits=dict(root_order=MAX_ROOT_ORDER,polynomial_slots=MAX_POLY_SLOTS,
        reductions=MAX_REDUCTIONS,float_summands=MAX_SUM_TERMS,numeric_bytes=MAX_BYTES),
        scalar_transfer_scope="fixed clean-input low histories, not every possible classical rewriting",
        bound_scope="states initially supported in ONE coarse sector; not arbitrary coherent multi-sector inputs",
        exact_backend="python-flint fmpz_poly integer cyclotomic reduction; no floating zero cutoff")
    path=Path("out")/("additive_phase_structure_"+datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")+".json")
    ok=exp.finish(report_path=path,rows=rows,metadata=metadata)
    print(f"report: {path}")
    if not ok:
        raise SystemExit(1)


if __name__=="__main__":
    main()
