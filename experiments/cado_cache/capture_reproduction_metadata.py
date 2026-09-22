"""Administrative provenance for a fresh reproduction; never launches science."""
import argparse,hashlib,json,os
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--stage',choices=['freeze','finalize'],required=True);a=p.parse_args();root=a.root.resolve()
keys=['LD_PRELOAD','LD_LIBRARY_PATH','MALLOC_CONF','MALLOC_ARENA_MAX','TCMALLOC_RELEASE_RATE']
env={k:os.environ.get(k) for k in keys}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if a.stage=='freeze':
 files=[root/f for f in ['prepare_arms_v1.py','build_arms_v1.py','cache_experiment_v1.py','build_manifest_v1.json']]+list((root/'arms').rglob('*.cpp'))+list((root/'arms').glob('*.patch'))
 (root/'frozen_versions_v1.json').write_text(json.dumps([{'path':str(f),'sha256':sha(f)} for f in files],indent=2)+'\n')
 (root/'allocator_environment_v1.json').write_text(json.dumps(env,indent=2)+'\n')
else:
 assert env==json.loads((root/'allocator_environment_v1.json').read_text()),'allocator environment changed'
 v=json.loads((root/'science_v1/validate_report_v1.json').read_text());r=json.loads((root/'science_v1/performance_report_v1.json').read_text())
 assert v['ok'] and r['ok']
 (root/'provenance_v1.json').write_text(json.dumps({'allocator_environment':env,'experiment_cpu_reported_seconds':r['metadata']['aggregate_cpu_seconds'],'cpu_rounding_upper_allowance_seconds':.02*(len(v['rows'])+len(r['rows']))},indent=2)+'\n')
