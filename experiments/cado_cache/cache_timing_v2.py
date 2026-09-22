"""Reviewed timing-only rerun: unchanged CADO binaries, inputs, and 20-row schedule.
Prediction: all caches match frozen v1 baselines. RSS/time improvement is open.
The passed arithmetic/capacity/omission gates are reused by exact evidence identity.
"""
from __future__ import annotations
import argparse,hashlib,json,os,platform,re,resource,signal,subprocess,sys,threading,time
from pathlib import Path
p=argparse.ArgumentParser()
for name in ['prior-root','inputs-manifest','output','project']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--inputs-sha256',required=True);p.add_argument('--validation-sha256',required=True)
a=p.parse_args();root=a.prior_root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(a.project.resolve()))
from lab.harness import Experiment
exp=Experiment('rsa73b_corrected_timing_v2',doc=__doc__)
exp.predict('gate','All pinned inputs, runtime identities, full usage fields, resource caps and 20 exact cache hashes pass')
exp.must_fail('omission','Reuse exact successful v1 detector of single-coefficient omission; no repeat science')
rows=[];failure=None;reused=False;manifest=None;phase_start=time.monotonic()
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def self_cpu():
 r=resource.getrusage(resource.RUSAGE_SELF);return r.ru_utime+r.ru_stime

def ledger():
 return 10+self_cpu()+sum(r['user_seconds']+r['system_seconds']+.02 for r in rows if 'user_seconds' in r)

def limits():
 resource.setrlimit(resource.RLIMIT_AS,(4294967296,4294967296));resource.setrlimit(resource.RLIMIT_CPU,(45,45))

def verify_entry(e):
 path=root/e['path']
 if path.stat().st_size!=e['bytes'] or sha(path)!=e['sha256']:raise RuntimeError('input identity mismatch: '+e['path'])
 return path

def measured_child(cmd,label):
 if ledger()+45>600:raise RuntimeError('CPU accounting reserve exhausted')
 log=out/(label+'.log');usage=out/(label+'.usage')
 env=os.environ.copy();env.pop('QSIM_DIAGNOSTIC',None)
 env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 timed_out=threading.Event();state_lock=threading.Lock();done=False;timer=None
 with log.open('x') as stream:
  start=time.monotonic()
  proc=subprocess.Popen(['/usr/bin/time','-f','%M %U %S %e','-o',str(usage),*map(str,cmd)],stdout=stream,stderr=subprocess.STDOUT,env=env,preexec_fn=limits,start_new_session=True)
  def expire():
   nonlocal done
   with state_lock:
    if done:return
    timed_out.set()
    try:os.killpg(proc.pid,signal.SIGKILL)
    except ProcessLookupError:pass
  try:
   timer=threading.Timer(max(0.0,start+60-time.monotonic()),expire);timer.start()
   rc=proc.wait()
   end=time.monotonic()  # Immediate timestamp, before watchdog cleanup or parsing.
  finally:
   with state_lock:done=True
   if timer is not None:timer.cancel();timer.join()
  if timed_out.is_set():rc=124
  stream.write(f'\n[exit {rc}]\n')
 rec={'label':label,'command':list(map(str,cmd)),'exit_code':rc,'timed_out':timed_out.is_set(),'launcher_elapsed_seconds':end-start,'log':log.name,'usage':usage.name}
 rows.append(rec)
 def save(): (out/'rows_v2.json').write_text(json.dumps(rows,indent=2)+'\n')
 save()
 if not usage.exists():raise RuntimeError(label+': missing GNU time usage')
 raw=usage.read_text();rec['raw_usage']=raw
 fields=raw.splitlines()[-1].split() if raw.splitlines() else []
 if len(fields)!=4 or not re.fullmatch(r'\d+',fields[0]) or any(not re.fullmatch(r'\d+(\.\d+)?',f) for f in fields[1:]):
  save();raise RuntimeError(label+': malformed GNU time usage')
 rec.update(peak_rss_kib=int(fields[0]),user_seconds=float(fields[1]),system_seconds=float(fields[2]),gnu_elapsed_seconds=float(fields[3]));save()
 if rc:raise RuntimeError(label+': child exit '+str(rc))
 if ledger()>600:raise RuntimeError('CPU accounting cap exceeded')
 return rec

try:
 if sha(a.inputs_manifest)!=a.inputs_sha256:raise RuntimeError('input manifest hash mismatch')
 manifest=json.loads(a.inputs_manifest.read_text())
 if manifest['prior_validation']['sha256']!=a.validation_sha256:raise RuntimeError('reviewed validation hash mismatch')
 for e in manifest['artifacts']:verify_entry(e)
 validation=json.loads(verify_entry(manifest['prior_validation']).read_text());prior=json.loads(verify_entry(manifest['prior_performance']).read_text())
 if not validation['ok'] or validation['metadata']['failure'] is not None or not validation['metadata']['mutant_detected'] or not all(c['ok'] for c in validation['checks']):raise RuntimeError('prior validation unsuccessful')
 if not prior['ok'] or not any(x['label']=='P2_applicability' and x['exit_code']==0 for x in prior['rows']):raise RuntimeError('prior P2 applicability absent')
 reused=True
 for key,wanted in manifest['allocator_environment'].items():
  if os.environ.get(key)!=wanted:raise RuntimeError('allocator environment changed: '+key)
 exes={}
 for exe in manifest['executables']:
  path=verify_entry(exe)
  for lib in exe['libraries']:
   if sha(lib['path'])!=lib['sha256']:raise RuntimeError('runtime library hash changed: '+lib['path'])
  if path.name=='build_matcache':exes[exe['arm']]=path
 if set(exes)!={'U','R','B','RB'}:raise RuntimeError('missing executable arm')
 fixtures={}
 for f in manifest['fixtures']:
  path=verify_entry(f)
  for suffix,dimension in [('rw',f['rows']),('cw',f['cols'])]:
   side=next(e for e in f['sidecars'] if e['path'].endswith('.'+suffix+'.bin'))
   if side['bytes']!=4*dimension:raise RuntimeError('sidecar dimension mismatch')
   verify_entry(side)
  fixtures[f['fixture']]=(f,path)
 identities={'interpreter':{'path':str(Path(sys.executable).resolve()),'sha256':sha(sys.executable),'version':sys.version},'gnu_time':{'path':'/usr/bin/time','sha256':sha('/usr/bin/time')},'platform':platform.platform(),'allocator_environment':manifest['allocator_environment'],'input_manifest_sha256':a.inputs_sha256,'driver_sha256':sha(__file__),'prior_validation_sha256':a.validation_sha256}
 (out/'runtime_identities_v2.json').write_text(json.dumps(identities,indent=2)+'\n')
 sequence=[('P1_small',0,['U','R','B','RB']),('P1_large',0,['U','R','B','RB']),('P1_large',1,['RB','B','R','U']),('P1_large',2,['R','U','RB','B']),('P2_default',0,['U','R','B','RB'])]
 for name,repetition,order in sequence:
  f,path=fixtures[name]
  for arm in order:
   label=f'{name}_{repetition}_{arm}';dest=out/label;dest.mkdir(exist_ok=False)
   cmd=[exes[arm],'--matrix-file',path,'-impl','bucket','-direction','right','-tmpdir',dest]
   if name.startswith('P1'):cmd+=['matmul_bucket_methods=vsc']
   rec=measured_child(cmd,label)
   cache=dest/(path.stem+'-bucket.bin');actual=sha(cache)
   rec.update(arm=arm,fixture=name,repetition=repetition,input_bytes=f['bytes'],cache_bytes=cache.stat().st_size,cache_sha256=actual)
   if actual!=f['cache_sha256'] or cache.stat().st_size!=f['cache_bytes']:raise RuntimeError(label+': frozen baseline cache mismatch')
   print('MEASURE',json.dumps(rec),flush=True)
 if len(rows)!=20:raise RuntimeError('incorrect row count')
except Exception as exc:failure=repr(exc)
finally:
 (out/'rows_v2.json').write_text(json.dumps(rows,indent=2)+'\n')
 exp.check('gate',failure is None,failure or 'All 20 corrected timing rows match frozen exact caches and required resource fields')
 exp.fail_check('omission',reused,'reused frozen successful v1 omission detector by required SHA256; not rerun')
 exp.finish(report_path=out/'performance_report_v2.json',rows=rows,metadata={'failure':failure,'driver_sha256':sha(__file__),'input_manifest_sha256':a.inputs_sha256,'validation_sha256':a.validation_sha256,'self_cpu_from_process_start_seconds':self_cpu(),'cpu_accounting_debit_seconds':ledger(),'prior_administrative_allowance_seconds':10,'prior_allowance_is_certified_bound':False,'cpu_reading_rounding_allowance_per_child_seconds':.02,'phase_wall_seconds':time.monotonic()-phase_start,'interval':'launcher-inclusive, blocking wait, timestamp before watchdog cleanup','gnu_time_resolution_seconds':.01,'cpu_resolution_seconds':.01,'io':'warm/unspecified page cache; no fsync; child RSS excludes global page cache','reuse':'unchanged v1 binaries, fixtures, direct sparse products, omission/resource/path validation; no regenerated fixtures or repeated diagnostics'})
