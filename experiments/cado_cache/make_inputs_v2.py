"""Freeze explicit, relocatable inputs for the reviewed timing-only rerun.
This is provenance preparation: reads and hashes existing evidence; runs no science.
"""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--prior-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.prior_root.resolve()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def entry(path):return {'path':str(path),'sha256':sha(root/path),'bytes':(root/path).stat().st_size}
validation=Path('science_v1/validate_report_v1.json');performance=Path('science_v1/performance_report_v1.json')
v=json.loads((root/validation).read_text());r=json.loads((root/performance).read_text());b=json.loads((root/'build_manifest_v1.json').read_text());frozen=json.loads((root/'frozen_versions_v1.json').read_text());prov=json.loads((root/'provenance_v1.json').read_text())
assert v['ok'] and r['ok'] and v['metadata']['mutant_detected']
fixtures=[]
for f in r['metadata']['fixtures']:
 path=Path('science_v1/fixtures')/(f['fixture']+'.bin');e=entry(path)
 assert e['sha256']==f['sha256'] and e['bytes']==f['bytes']
 side=[]
 for suffix,n in [('rw',f['rows']),('cw',f['cols'])]:
  x=entry(path.with_suffix('.'+suffix+'.bin'));assert x['bytes']==4*n;side.append(x)
 baseline=next(x for x in r['rows'] if x.get('fixture')==f['fixture'] and x.get('arm')=='U')
 fixtures.append({**f,**e,'sidecars':side,'cache_sha256':baseline['cache_sha256'],'cache_bytes':baseline['cache_bytes']})
artifacts=[entry(validation),entry(performance),entry(Path('build_manifest_v1.json')),entry(Path('frozen_versions_v1.json')),entry(Path('provenance_v1.json'))]
# Bind the exact source identities already frozen before the v1 science run.
for f in frozen:
 path=Path(f['path']);parts=path.parts
 if 'arms' in parts:relative=Path(*parts[parts.index('arms'):])
 else:relative=Path(path.name)
 e=entry(relative);assert e['sha256']==f['sha256'];artifacts.append(e)
executables=[]
for arm in b:
 assert arm['exit_code']==0
 for old in arm['executables']:
  e=entry(Path('bin')/arm['arm']/Path(old['path']).name);assert e['sha256']==old['sha256']
  executables.append({**e,'arm':arm['arm'],'libraries':old['libraries']})
p2=next(x for x in r['rows'] if x['label']=='P2_applicability')
assert p2['exit_code']==0
log=entry(Path('science_v1/P2_applicability.log'));assert 'QSIM RAW 0' in (root/log['path']).read_text();artifacts.append(log)
obj={'schema_version':2,'prior_validation':entry(validation),'prior_performance':entry(performance),'artifacts':artifacts,'executables':executables,'fixtures':fixtures,'allocator_environment':prov['allocator_environment'],'prior_cpu_administrative_allowance_seconds':10,'prior_reported_cpu_seconds':prov['experiment_cpu_reported_seconds'],'prior_rounding_allowance_seconds':prov['cpu_rounding_upper_allowance_seconds'],'prior_cpu_limitation':'10s is an administrative allowance, not a certified bound; v1 interpreter setup was unmeasured','sidecar_provenance':'rw/cw hashes first recorded in this v2 manifest; v1 pinned reader dimensions are checked by byte lengths','v1_report_hashes':{'validation':sha(root/validation),'performance':sha(root/performance)}}
a.output.write_text(json.dumps(obj,indent=2)+'\n');print(sha(a.output))
