"""Actual pinned-CADO cache lifetime validation and bounded RSS experiment.
Predictions: byte-identical caches and independent sparse-XOR products for U/R/B/RB;
terminal VSC raw and consumed staging capacities vanish only for builder patch.
Resource savings are an open hypothesis, not a success criterion.
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, resource, shutil, subprocess, sys, time
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser()
for name in ['source','build','bin-root','output','project']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--phase',choices=['validate','performance'],required=True)
a=p.parse_args()
for name in ['source','build','bin_root','output','project']:setattr(a,name,getattr(a,name).resolve())
a.output.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(a.project))
from lab.harness import Experiment
exp=Experiment('rsa73b_'+a.phase,doc=__doc__)
exp.predict('gate','All required exact-output, applicability and resource-mechanism gates pass; performance requires prior validation')
exp.must_fail('omission','Actual CADO product from F1 with its unique terminal coefficient omitted disagrees with original reference')
records=[];checks=[];t0=time.monotonic();parent0=time.process_time()
prior=json.loads((a.output/'validate_report_v1.json').read_text()) if a.phase=='performance' else None
prior_cpu=prior['metadata']['aggregate_cpu_seconds'] if prior else 0.0

def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while block:=f.read(1048576):h.update(block)
 return h.hexdigest()
def spent():return prior_cpu+time.process_time()-parent0+sum(x.get('user_seconds',0)+x.get('system_seconds',0) for x in records)
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,)*2)
 resource.setrlimit(resource.RLIMIT_CPU,(45,45))
def child(cmd,label,diagnostic=False):
 if spent()+45>600:raise RuntimeError('cumulative CPU budget reserve exhausted')
 log=a.output/(label+'.log');usage=a.output/(label+'.usage')
 env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 env.pop('QSIM_DIAGNOSTIC',None)
 if diagnostic:env['QSIM_DIAGNOSTIC']='1'
 start=time.monotonic()
 with log.open('x') as out:
  proc=subprocess.Popen(['/usr/bin/time','-f','%M %U %S','-o',str(usage),*map(str,cmd)],stdout=out,stderr=subprocess.STDOUT,env=env,preexec_fn=limits,start_new_session=True)
  try:rc=proc.wait(timeout=60)
  except subprocess.TimeoutExpired:
   import signal
   os.killpg(proc.pid,signal.SIGKILL);proc.wait();rc=124
  out.write(f'\n[exit {rc}]\n'.encode() if 'b' in out.mode else f'\n[exit {rc}]\n')
 rec={'label':label,'command':list(map(str,cmd)),'exit_code':rc,'wall_seconds':time.monotonic()-start,'log':str(log),'diagnostic':diagnostic}
 if usage.exists():
  values=usage.read_text().splitlines()[-1].split()
  if len(values)==3:
   rec.update(peak_rss_kib=int(values[0]),user_seconds=float(values[1]),system_seconds=float(values[2]))
 records.append(rec)
 (a.output/(a.phase+'_rows_v1.json')).write_text(json.dumps(records,indent=2)+'\n')
 if rc:raise RuntimeError(f'{label}: child exit {rc}')
 if spent()>600:raise RuntimeError('cumulative CPU budget exceeded')
 return rec,log.read_text()

def csr(path,nr,nc,rr,cc):
 rr=np.asarray(rr,dtype=np.int64);cc=np.asarray(cc,dtype=np.uint32)
 order=np.argsort(rr,kind='stable');rr=rr[order];cc=cc[order]
 rw=np.bincount(rr,minlength=nr).astype(np.uint32);cw=np.bincount(cc,minlength=nc).astype(np.uint32)
 starts=np.arange(nr,dtype=np.int64)+np.r_[0,np.cumsum(rw[:-1],dtype=np.int64)]
 words=np.empty(nr+len(cc),dtype=np.uint32);words[starts]=rw
 words[np.arange(len(cc),dtype=np.int64)+rr+1]=cc
 words.tofile(path);rw.tofile(path.with_suffix('.rw.bin'));cw.tofile(path.with_suffix('.cw.bin'))
 return {'file':str(path),'rows':nr,'cols':nc,'coefficients':len(cc),'sha256':digest(path),'bytes':path.stat().st_size}

def vec(n,indices):
 x=np.arange(n,dtype=np.uint64)*np.uint64(0x9e3779b97f4a7c15)
 x=(x & np.uint64(0x7fffffff00000000)) | np.uint64(1<<63)
 for k,i in enumerate(indices):x[i]|=np.uint64(1<<k)
 return x

def bench(arm,base,nr,nc,rr,cc,transpose,label,method='vsc',diagnostic=True,mutant=False):
 d=a.output/label;d.mkdir()
 m=d/'matrix.bin'
 if transpose:csr(m,nc,nr,cc,rr)
 else:csr(m,nr,nc,rr,cc)
 nin=nr if transpose else nc;nout=nc if transpose else nr
 support=sorted(set(map(int,rr if transpose else cc)))
 indices=support if len(support)<=24 else [0,16383,16384,32767,32768,nin-1]
 x=vec(nin,indices);v=d/'vector.bin';x.tofile(v)
 expected=np.zeros(nout,dtype=np.uint64)
 np.bitwise_xor.at(expected,cc if transpose else rr,x[rr if transpose else cc])
 cmd=[a.bin_root/arm/'bench_matcache','--nmax','0','--nchecks','0','-impl','bucket','srcvec='+str(v)]
 if method:cmd+=['matmul_bucket_methods='+method]
 if transpose:cmd+=['-t']
 cmd+=['--',m]
 rec,text=child(cmd,label,diagnostic)
 result=np.fromfile(str(v)+'.dst',dtype=np.uint64)
 if not np.array_equal(result,expected):raise RuntimeError(label+': direct sparse XOR mismatch')
 cache=Path(str(m)+'-bucket'+('T' if transpose else '')+'.bin')
 rec.update(cache_sha256=digest(cache),cache_bytes=cache.stat().st_size)
 if diagnostic:
  raw=re.findall(r'^QSIM RAW (\d+)$',text,re.M)
  cap=re.findall(r'^QSIM [XC] (\d+) (\d+)$',text,re.M)
  steps=[list(map(int,x)) for x in re.findall(r'^QSIM STEP (\d+) (\d+)$',text,re.M)]
  if method=='vsc':
   if not raw or not cap:raise RuntimeError(label+': missing VSC path')
   patched='B' in arm
   if any((int(v)==0)!=patched for v in raw):raise RuntimeError(label+': raw capacity control')
   if not any(int(x)>0 for x,y in cap):raise RuntimeError(label+': empty capacity control')
   if not all(int(y)==0 if patched else int(x)==int(y) for x,y in cap):raise RuntimeError(label+': staging capacity control')
  elif raw:raise RuntimeError(label+': negative applicability unexpectedly reached VSC')
  rec.update(raw_capacity_words=list(map(int,raw)),staging_capacity_pairs=[list(map(int,x)) for x in cap],steps=steps)
 # Reload in a fresh actual CADO process, preserving the original cache bytes.
 _,reload_text=child(cmd,label+'_reload',False)
 if 'Reusing cache file' not in reload_text or not np.array_equal(np.fromfile(str(v)+'.dst',dtype=np.uint64),expected):raise RuntimeError(label+': reload mismatch')
 return rec,result,indices

def generate_perf(name,nr,nc,degree):
 start=time.monotonic();cpu=time.process_time();d=a.output/'fixtures';d.mkdir(exist_ok=True)
 path=d/(name+'.bin');cw=np.zeros(nc,dtype=np.uint64)
 with path.open('xb') as f:
  for i in range(0,nr,4096):
   rows=np.arange(i,min(i+4096,nr),dtype=np.uint64)
   cols=((rows[:,None]*131+np.arange(degree,dtype=np.uint64)*16381)%nc).astype(np.uint32);cols.sort(axis=1)
   words=np.empty((len(rows),degree+1),dtype=np.uint32);words[:,0]=degree;words[:,1:]=cols;words.tofile(f)
   cw+=np.bincount(cols.ravel(),minlength=nc).astype(np.uint64)
 np.full(nr,degree,dtype=np.uint32).tofile(path.with_suffix('.rw.bin'));cw.astype(np.uint32).tofile(path.with_suffix('.cw.bin'))
 return path,{'fixture':name,'rows':nr,'cols':nc,'degree':degree,'bytes':path.stat().st_size,'sha256':digest(path),'generation_wall_seconds':time.monotonic()-start,'generation_cpu_seconds':time.process_time()-cpu,'generator_process_highwater_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}

fixtures=[];failure=None;mutant_detected=False
try:
 manifest=json.loads((a.bin_root.parent/'build_manifest_v1.json').read_text())
 for arm in manifest:
  if arm['exit_code']!=0:raise RuntimeError('unbuilt arm')
  for exe in arm['executables']:
   target=a.bin_root/arm['arm']/Path(exe['path']).name
   if digest(target)!=exe['sha256']:raise RuntimeError('executable identity changed')
   for lib in exe['libraries']:
    if digest(lib['path'])!=lib['sha256']:raise RuntimeError('runtime library identity changed')
 if a.phase=='validate':
  f1r=np.array([0,0,0,1,1,3,3,3],dtype=np.int64);f1c=np.array([0,43690,43691,87381,87382,131072,0,0],dtype=np.int64)
  deg=np.repeat([17,6,1],16384);f2r=np.repeat(np.arange(49152),deg)
  f2c=np.r_[np.tile(np.arange(17)*65536,16384),np.tile(np.arange(6)*65536,16384),np.full(16384,16*65536)]
  for name,nr,nc,rr,cc,method in [('F1',4,131073,f1r,f1c,'vsc'),('F2',49152,17*65536,f2r,f2c,'vsc'),('F3',16,16,np.arange(16),np.arange(16),'')]:
   for tr in ([False,True] if name!='F3' else [False]):
    baseline=None
    for arm in ['U','R','B','RB']:
     rec,_,_=bench(arm,None,nr,nc,rr,cc,tr,name+('_left_' if tr else '_right_')+arm,method)
     if baseline is None:baseline=rec['cache_sha256']
     if rec['cache_sha256']!=baseline:raise RuntimeError('arm cache byte mismatch')
     if name=='F2' and not tr and rec['steps']!=[[1,16384],[3,16384],[17,16384]]:raise RuntimeError('F2 defer-step precondition')
  keep=~((f1r==3)&(f1c==131072))
  _,got,indices=bench('RB',None,4,131073,f1r[keep],f1c[keep],False,'F1_omission_mutant')
  # Use the original packed input; basis lane for131072 is absent from the mutant's support.
  # Compare all-ones lane63, where deleting one coefficient necessarily flips row3.
  expected=np.zeros(4,dtype=np.uint64);origx=vec(131073,indices)
  np.bitwise_xor.at(expected,f1r,origx[f1c])
  mutant_detected=not np.array_equal(got,expected) and ((int(got[3])^int(expected[3]))>>63)==1
  if not mutant_detected:raise RuntimeError('wrong-result mutant survived')
 else:
  if not prior['ok'] or not prior['metadata']['mutant_detected']:raise RuntimeError('prior validation missing')
  mutant_detected=True
  for name,nr in [('P1_small',131072),('P1_large',262144)]:
   path,info=generate_perf(name,nr,1048576,64);fixtures.append(info)
   orders=[['U','R','B','RB']] if name=='P1_small' else [['U','R','B','RB'],['RB','B','R','U'],['R','U','RB','B']]
   reference=None
   for repetition,order in enumerate(orders):
    for arm in order:
     label=f'{name}_{repetition}_{arm}';d=a.output/label;d.mkdir()
     cmd=[a.bin_root/arm/'build_matcache','--matrix-file',path,'-impl','bucket','-direction','right','-tmpdir',d,'matmul_bucket_methods=vsc']
     rec,_=child(cmd,label)
     cache=d/(path.stem+'-bucket.bin');sha=digest(cache)
     if reference is None:reference=sha
     if sha!=reference:raise RuntimeError('timed cache byte mismatch')
     rec.update(arm=arm,fixture=name,repetition=repetition,cache_bytes=cache.stat().st_size,cache_sha256=sha,input_bytes=info['bytes'])
     print('MEASURE',json.dumps(rec),flush=True)
  # Prescribed optional default-policy contrast, one run per arm after path gate.
  path,info=generate_perf('P2_default',2097152,1048576,2);fixtures.append(info)
  d=a.output/'P2_applicability';d.mkdir()
  cmd=[a.bin_root/'RB'/'build_matcache','--matrix-file',path,'-impl','bucket','-direction','right','-tmpdir',d]
  _,text=child(cmd,'P2_applicability',True)
  if not re.search(r'^QSIM RAW 0$',text,re.M):raise RuntimeError('P2 default VSC applicability failed')
  reference=None
  for arm in ['U','R','B','RB']:
   label='P2_default_0_'+arm;d=a.output/label;d.mkdir()
   rec,_=child([a.bin_root/arm/'build_matcache','--matrix-file',path,'-impl','bucket','-direction','right','-tmpdir',d],label)
   cache=d/(path.stem+'-bucket.bin');sha=digest(cache)
   if reference is None:reference=sha
   if sha!=reference:raise RuntimeError('P2 cache byte mismatch')
   rec.update(arm=arm,fixture='P2_default',repetition=0,cache_bytes=cache.stat().st_size,cache_sha256=sha,input_bytes=info['bytes'])
   print('MEASURE',json.dumps(rec),flush=True)
except Exception as exc:
 failure=str(exc)
finally:
 exp.check('gate',failure is None,failure or 'Exact output, applicability and capacity gates satisfied')
 exp.fail_check('omission',mutant_detected,'actual CADO coefficient-omission result differs' if a.phase=='validate' else 'reused unchanged validated mutant from validate_report_v1.json')
 exp.finish(report_path=a.output/(a.phase+'_report_v1.json'),rows=records,metadata={'failure':failure,'fixtures':fixtures,'aggregate_cpu_seconds':spent(),'phase_wall_seconds':time.monotonic()-t0,'mutant_detected':mutant_detected,'driver_sha256':digest(__file__),'source':str(a.source),'build':str(a.build),'io':'warm/unspecified OS page cache; process RSS excludes global page cache','limits':{'address_space_bytes':4*1024**3,'child_cpu_seconds':45,'child_wall_seconds':60,'aggregate_cpu_seconds':600}})
