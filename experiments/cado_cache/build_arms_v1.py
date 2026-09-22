"""Build static CADO arm executables; engineering only, no science execution."""
import argparse,subprocess,time,json,hashlib,shutil,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--budget-seconds',type=float,default=420);a=p.parse_args()
t0=time.monotonic();records=[]
for arm in ['U','R','RB','B']:
 for f in ['matrix_u32.cpp','matmul-bucket.cpp']:
  shutil.copyfile(a.output/'arms'/('sources_'+arm)/f,a.source/'linalg/bwc'/f)
 cmd=['cmake','--build',str(a.build),'--target','build_matcache','bench_matcache','-j2']
 log=a.output/f'build_{arm}_v1.log'
 with log.open('w') as stream:
  start=time.monotonic()
  try:r=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,timeout=max(1,a.budget_seconds-(start-t0)));rc=r.returncode
  except subprocess.TimeoutExpired:rc=124
  stream.write(f'\n[exit {rc}]\n')
 rec={'arm':arm,'command':cmd,'exit_code':rc,'wall_seconds':time.monotonic()-start,'executables':[]}
 if rc==0:
  dest=a.output/'bin'/arm;dest.mkdir(parents=True,exist_ok=True)
  for name in ['build_matcache','bench_matcache']:
   src=a.build/'linalg/bwc'/name;target=dest/name;shutil.copy2(src,target)
   ldd=subprocess.check_output(['ldd',str(target)],text=True)
   (dest/(name+'.ldd.txt')).write_text(ldd)
   libs=[]
   for line in ldd.splitlines():
    for word in line.split():
     if word.startswith('/') and Path(word).is_file():libs.append({'path':str(Path(word).resolve()),'sha256':hashlib.sha256(Path(word).read_bytes()).hexdigest()})
   rec['executables'].append({'path':str(target),'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'libraries':libs})
 records.append(rec);(a.output/'build_manifest_v1.json').write_text(json.dumps(records,indent=2)+'\n')
 print(arm,rc,rec['wall_seconds'],flush=True)
 if rc:raise SystemExit(rc)
