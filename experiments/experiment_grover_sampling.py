"""Exact structured Grover full-output sampling: accepted bounded TODO71 plan.

Stages are reference, validate, case and compare. Every output path is supplied;
existing reports are never overwritten. Resource cases import only the standard
library and sampler. Reference/validation additionally import the existing dense
helper, NumPy, Walsh and lab.harness. See the accepted design for decision rules.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import random
import resource
import signal
import statistics
import sys
import time
import tracemalloc


WIDTHS = (4, 8, 16, 32, 64, 128)
STEPS = (0, 1, 3, 8)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def successful(path):
    obj = json.loads(path.read_text())
    if obj.get('ok') is not True or obj.get('warnings') or not obj.get('checks') or not all(c['ok'] for c in obj['checks']):
        raise ValueError(f'unusable correctness report: {path}')
    return obj


def independent_mutant(n, steps, mode):
    """Direct dense operations, independent of candidate counts/environments."""
    import numpy as np
    from walsh import wht
    x = np.arange(1 << n, dtype=np.uint64)
    selected = x % 2 == 1
    for i in range(n-1):
        if mode == 'drop_edge' and i == 0:
            continue
        selected &= ~((((x >> i)&1) == 1) & (((x >> (i+1))&1) == 1))
    state = np.ones(1 << n)
    for _ in range(steps):
        if mode != 'omit_oracle':
            state[selected] *= -1
        state[:] = 2*state.mean()-state
    if mode == 'dephase':
        return tuple(float(np.dot(state,state)/(1 << (2*n))) for _ in state)
    return tuple(float(v*v) for v in wht(state)/(1 << n))


def reference(args):
    from lab.harness import Experiment
    from experiments.grover_queries import dense
    exp = Experiment('grover71_reference_v1', doc=__doc__)
    exp.predict('P1', 'Independent dense full distributions normalize and match exact hand fixtures.')
    for key in ('omit_oracle', 'dephase', 'drop_edge', 'label_shift'):
        exp.must_fail(key, 'Named dense/verifier mutation differs on asymmetric n4,t3.')
    rows = []
    for n in range(2,9):
        for t in STEPS:
            p = dense(n,t,tuple(range(1 << n)))
            exp.check('P1', all(0 <= z <= 1+1e-11 for z in p) and abs(sum(p)-1) <= 1e-11,
                      f'n={n}; t={t}; complete normalization')
            rows.append(dict(n=n, steps=t, probabilities=p))
    p = next(r['probabilities'] for r in rows if r['n']==4 and r['steps']==3)
    exp.check('P1', p[0] == float(Fraction(13225,16384)) and p[1] == float(Fraction(729,16384)), 'n4,t3 hand scalars')
    for t in (0,1):
        actual = next(r['probabilities'] for r in rows if r['n']==2 and r['steps']==t)
        exp.check('P1', actual == ((1.,0.,0.,0.) if t==0 else (.25,)*4), f'n2,t{t} exact boundary fixture')
    mutant_rows = []
    for mode in ('omit_oracle','dephase','drop_edge','label_shift'):
        q = p[1:]+p[:1] if mode=='label_shift' else independent_mutant(4,3,mode)
        error = max(abs(a-b) for a,b in zip(p,q))
        exp.fail_check(mode,error > 1e-11, f'n4,t3; max disagreement={error}')
        mutant_rows.append(dict(mutant=mode,n=4,steps=3,max_absolute_disagreement=error))
    exp.finish(report_path=args.report,rows=rows,metadata={
        'mutants':mutant_rows,'candidate_imported':False,
        'independence':'Existing dense phase/mean/Walsh route; author has read it, algorithmic independence only.',
        'shared_components':'Oracle contract; NumPy and Walsh used by dense only; no candidate recurrence.',
        'dense_source_sha256':sha('experiments/grover_queries.py'),
        'walsh_source_sha256':sha('walsh.py'),'tolerance':1e-11})


def sampler_module():
    if __package__ == 'experiments':
        from . import grover_sampling
    else:
        import grover_sampling
    return grover_sampling


def validate(args):
    from lab.harness import Experiment
    from experiments.grover_queries import compact, query_masks
    ref = successful(args.reference)
    sm = sampler_module()
    exp = Experiment('grover71_candidate_validation_v1', doc=__doc__)
    exp.predict('P1','Every tiny full distribution/prefix agrees with independent dense; integer identities exact.')
    exp.predict('P2','C114 selected exact scalars, n128 analytic rare branches/leaves, and conditioned outputs agree.')
    exp.predict('P3','Actual RNG rejection, all threshold boundaries, deterministic branches and invalid inputs obey contract.')
    exp.must_fail('C1','Deleting coherent cross terms changes physical mass on n4,t3.')
    exp.must_fail('C2','Zero-probability prefix conditioning is rejected and consumes no bits.')
    rows=[]
    for fixture in ref['rows']:
        n,t,p = fixture['n'],fixture['steps'],fixture['probabilities']
        roll,cache = (sm.StructuredGroverSampler(n,t,m) for m in ('rolling','cached'))
        ok,err,count = True,0.,0
        for length in range(n+1):
            for prefix in range(1 << length):
                wr,wc = roll.prefix_numerator(prefix,length),cache.prefix_numerator(prefix,length)
                expected=sum(p[y] for y in range(1 << n) if (y & ((1 << length)-1)) == prefix)
                err=max(err,abs(float(Fraction(wr,roll.denominator))-expected))
                ok &= wr==wc and wr>=0
                if length<n:
                    ok &= wr == roll.prefix_numerator(prefix,length+1)+roll.prefix_numerator(prefix+(1 << length),length+1)
                count+=1
        ok &= roll.prefix_numerator(0,0)==roll.denominator
        exp.check('P1',ok and err<=1e-11,f'n={n}; t={t}; prefixes={count}; max_error={err}')
        rows.append(dict(n=n,steps=t,prefixes=count,max_absolute_error=err))
    for n in range(2,23):
        for t in STEPS:
            masks=query_masks(n)
            exact=compact(n,t,masks)
            s=sm.StructuredGroverSampler(n,t)
            exp.check('P2',tuple(s.prefix_probability(y,n) for y in masks)==exact,f'C114 exact n={n}; t={t}')
    # Analytic fixture evaluated separately with reduced Fractions, not candidate coefficients.
    old=json.loads(args.preconditions.read_text())
    frozen128=next(r for r in old['rows'] if r['n']==128)
    f0,f1=0,1
    for _ in range(128):
        f0,f1=f1,f0+f1
    p=Fraction(f0,1 << 128)
    a,b=Fraction(1),Fraction(0)
    for _ in range(3):
        a,b=(1-4*p)*a-2*p*b,2*a+b
    rare=b*b*p/2
    analytic=dict(prefix1=rare,child01=rare/2,child11=rare/2,leaf1=(b*p)**2,leaf0=(a+b*p)**2)
    exp.check('P2',0<rare<Fraction(1,10**8) and rare==Fraction(frozen128['rare_prefix_mass'])
              and analytic['leaf0']==Fraction(frozen128['p0']) and analytic['leaf1']==Fraction(frozen128['p1']),
              'n128 independently evaluated Fraction fixture equals frozen preconditions')
    forced=[]
    for method in ('rolling','cached'):
        s=sm.StructuredGroverSampler(128,3,method)
        actual=dict(prefix1=s.prefix_probability(1,1),child01=s.prefix_probability(1,2),
                    child11=s.prefix_probability(3,2),leaf1=s.prefix_probability(1,128),leaf0=s.prefix_probability(0,128))
        exp.check('P2',actual==analytic,f'n128 {method}; exact analytic prefix/child/leaf probes')
        rng=random.Random(20260921)
        values=s.sample(16,rng.getrandbits,prefix=1,length=1)
        exp.check('P2',len(values)==16 and all(0<=x<(1 << 128) and x&1 for x in values),f'n128 {method};16 conditional outputs retain rare prefix')
        forced.append(dict(method=method,outputs=values))
    def forbidden(_):
        raise AssertionError('random bits consumed by deterministic branch')
    draws=iter((3,2))
    exp.check('P3',sm.exact_randbelow(3,lambda k:next(draws))==2,'rejection3 then acceptance2 at total3')
    exp.check('P3',sm.exact_randbelow(1,forbidden)==0 and sm.choose_bit(0,7,forbidden)==1
              and sm.choose_bit(7,0,forbidden)==0,'total1 and both zero branches consume no bits')
    for w0,w1 in ((5,7),(3,5)):
        bits=[sm.choose_bit(w0,w1,lambda k,v=v:v) for v in (0,w0-1,w0,w0+w1-1)]
        exp.check('P3',bits==[0,0,1,1],f'actual chooser all four boundaries; weights={w0,w1}')
    for bad in (-1,4,0.5,True):
        caught=False
        try:
            sm.exact_randbelow(3,lambda k:bad)
        except ValueError:
            caught=True
        exp.check('P3',caught,f'invalid random-bit output {bad!r} rejected')
    for n in WIDTHS:
        s=sm.StructuredGroverSampler(n,0)
        caught=False
        try:
            s.sample_one(forbidden,prefix=1,length=1)
        except ValueError:
            caught=True
        exp.fail_check('C2',caught and s.prefix_probability(1,1)==0,f'n={n}; impossible conditioning rejected')
        exp.check('P3',s.sample(2,forbidden)==[0,0],f'n={n}; t0 complete deterministic sampler')
    # Drive sample_one with lower endpoints; chooser thresholds are checked above.
    s=sm.StructuredGroverSampler(4,3)
    exp.check('P3',0 <= s.sample_one(lambda k:0)<16,'sample_one drives actual chooser with lower endpoints')
    exp.check('P3',s.sample(0,forbidden)==[] and s.sample(1,forbidden,prefix=1,length=4)==[1],
              'zero count and positive full-prefix conditioning')
    row=sm.branch_row(s.initial,0,0)
    E=next(s.environments())
    broken=tuple(tuple(0 if (i==0)!=(j==0) else E[i][j] for j in range(3)) for i in range(3))
    exp.fail_check('C1',sm.weight(row,E)!=sm.weight(row,broken),'actual environment mutant removes constant/predicate cross entries')
    exp.finish(report_path=args.report,rows=rows,metadata={
        'reference_sha256':sha(args.reference),'preconditions_sha256':sha(args.preconditions),
        'analytic128':{k:str(v) for k,v in analytic.items()},'forced_conditional_samples':forced,
        'fixed_seed_is_not_iid_guarantee':True,'helper_sha256':sha(sm.__file__)})


def case(args):
    # Keep reference libraries out of resource processes.
    sm=sampler_module()
    if args.n not in WIDTHS or args.method not in ('rolling','cached') or args.mode not in ('time','allocation'):
        raise ValueError('outside frozen resource cases')
    successful(args.validation)
    def invoke():
        return sm.seeded_samples(args.n,3,16,args.method,20260921)
    baseline_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    walls=[]
    traced=None
    outputs=None
    if args.mode=='time':
        for _ in range(3):
            outputs=None
            start=time.perf_counter_ns()
            outputs=invoke()
            walls.append((time.perf_counter_ns()-start)*1e-9)
    else:
        tracemalloc.start()
        outputs=invoke()
        _,traced=tracemalloc.get_traced_memory()
        tracemalloc.stop()
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    obj=dict(schema_version=1,kind='single-arm resource observation',method=args.method,n=args.n,
             steps=3,count=16,seed=20260921,mode=args.mode,outputs=outputs,
             wall_seconds=walls,median_wall_seconds=statistics.median(walls) if walls else None,
             traced_peak_bytes=traced,raw_process_peak_rss_bytes=rss,raw_import_peak_rss_bytes=baseline_rss,
             comparison_pending=True,preconditions_valid=rss<=(1 << 30) and len(outputs)==16 and all(0<=x<(1 << args.n) for x in outputs),
             validation_sha256=sha(args.validation),helper_sha256=sha(sm.__file__),
             imported_numpy='numpy' in sys.modules,imported_lab='lab' in sys.modules,
             measurement_scope='Complete sampler construction, coefficients/environments, RNG,16 integer outputs; imports/launch/locking/serialization excluded.',
             output_encoding='Python integers at known n; zero integers may be shared, not materialized n-character binary strings.',
             rss_definition='Linux ru_maxrss whole-process high-water bytes; includes imports; not incremental allocation.')
    args.report.write_text(json.dumps(obj,indent=2)+'\n')
    print(json.dumps(obj))


def compare(args):
    from lab.harness import Experiment
    successful(args.validation)
    exp=Experiment('grover71_resource_comparison_v1',doc=__doc__)
    exp.predict('P1','All24 cases conform, use frozen helper/validation and agree on seeded outputs.')
    exp.predict('P2','Rolling traced allocation is smaller at n>=32; time ratio recorded against2x aspiration.')
    exp.must_fail('C1','An all-zero histogram cannot validate n128; analytic rare mass stays positive.')
    val=successful(args.validation)
    exp.fail_check('C1',Fraction(val['metadata']['analytic128']['prefix1'])>0,'reuse exact rare-prefix gate; all-zero histogram is insufficient')
    rows=[]
    cases=[]
    for n in WIDTHS:
        reports={}
        for method in ('rolling','cached'):
            for mode in ('time','allocation'):
                path=args.reports_dir/f'{method}_n{n}_{mode}_v1.json'
                r=json.loads(path.read_text()); reports[method,mode]=r
                cases.append(dict(path=str(path),sha256=sha(path),report=r))
                ok=(r['n']==n and r['steps']==3 and r['count']==16 and r['method']==method
                    and r['mode']==mode and r['preconditions_valid'] is True
                    and r['validation_sha256']==sha(args.validation)
                    and r['helper_sha256']==val['metadata']['helper_sha256']
                    and r['raw_process_peak_rss_bytes']<=(1 << 30)
                    and not r['imported_numpy'] and not r['imported_lab'])
                ok &= (len(r['wall_seconds'])==3 and all(0<t<120 for t in r['wall_seconds'])) if mode=='time' else isinstance(r['traced_peak_bytes'],int)
                exp.check('P1',ok,f'n={n}; {method}/{mode}; provenance and resource contract')
        first=reports['rolling','time']['outputs']
        exp.check('P1',all(r['outputs']==first for r in reports.values()),f'n={n}; exact seeded outputs across both methods/modes')
        ratio=reports['rolling','time']['median_wall_seconds']/reports['cached','time']['median_wall_seconds']
        smaller=reports['rolling','allocation']['traced_peak_bytes']<reports['cached','allocation']['traced_peak_bytes']
        # This is an assessed prediction, never a reason to tune the fixture.
        if n>=32:
            exp.check('P2',smaller,f'n={n}; rolling has smaller measured traced peak')
        rows.append(dict(n=n,rolling_to_cached_time_ratio=ratio,rolling_smaller_traced_peak=smaller,
                         within_2x_runtime=ratio<=2,all_seeded_outputs_zero=all(x==0 for x in first),
                         rolling_time_seconds=reports['rolling','time']['median_wall_seconds'],
                         cached_time_seconds=reports['cached','time']['median_wall_seconds'],
                         rolling_peak_bytes=reports['rolling','allocation']['traced_peak_bytes'],
                         cached_peak_bytes=reports['cached','allocation']['traced_peak_bytes'],
                         rolling_rss_bytes=reports['rolling','allocation']['raw_process_peak_rss_bytes'],
                         cached_rss_bytes=reports['cached','allocation']['raw_process_peak_rss_bytes']))
    exp.finish(report_path=args.report,rows=rows,metadata={'case_artifacts':cases,
        'case_order':'n ascending, rolling then cached, time then allocation; timings exploratory, not randomized.',
        'baseline':'Same exact bond3 structured contraction; no advantage over algebraically identical best method asserted.'})


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=('reference','validate','case','compare'),required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--reference',type=Path)
    parser.add_argument('--preconditions',type=Path)
    parser.add_argument('--validation',type=Path)
    parser.add_argument('--reports-dir',type=Path)
    parser.add_argument('--n',type=int)
    parser.add_argument('--method')
    parser.add_argument('--mode')
    args=parser.parse_args()
    if args.report.exists():
        raise FileExistsError(args.report)
    args.report.parent.mkdir(parents=True,exist_ok=True)
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
        if os.environ.get(key)!='1':
            raise ValueError(f'set {key}=1')
    resource.setrlimit(resource.RLIMIT_CPU,(120,120))
    signal.alarm(120)
    sys.path.insert(0,str(Path.cwd()))
    dict(reference=reference,validate=validate,case=case,compare=compare)[args.stage](args)


if __name__=='__main__':
    main()
