"""Engineering setup: exact lifetime patches plus shared opt-in diagnostics."""
import argparse, difflib, json, hashlib
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
paths=['linalg/bwc/matrix_u32.cpp','linalg/bwc/matmul-bucket.cpp']
u={x:(a.source/x).read_text() for x in paths}
r=u[paths[0]].replace('#include <memory>','#include <memory>\n#include <utility>').replace('return { m, {','return { std::move(m), {')
b=u[paths[1]].replace('        mb.prepare_vsc_slices(V, main_i0, fence);','        mb.prepare_vsc_slices(V, main_i0, fence);\n        mb.rowhead = nullptr;\n        std::vector<uint32_t>().swap(mb.data);').replace('            S.x.clear();','            std::vector<uint16_t>().swap(S.x);').replace('            S.c.clear();','            std::vector<uint8_t>().swap(S.c);')
assert r.count('std::move(m)')==2 and b.count('swap(S.x)')==1 and b.count('swap(S.c)')==1
for name,changed in [('reader_v1.patch',{paths[0]:r}),('builder_v1.patch',{paths[1]:b})]:
 text=''.join(''.join(difflib.unified_diff(u[x].splitlines(True),y.splitlines(True),fromfile='a/'+x,tofile='b/'+x)) for x,y in changed.items())
 (a.output/name).write_text(text)
def diag(s):
 s=s.replace('        scratch3size = MAX(scratch3size, V->tbuf_space);', '''        if (std::getenv("QSIM_DIAGNOSTIC")) {
            fprintf(stderr, "QSIM RAW %zu\\n", mb.data.capacity());
            for (auto const & step : V->steps)
                fprintf(stderr, "QSIM STEP %u %u\\n", step.defer, step.nrows);
        }
        scratch3size = MAX(scratch3size, V->tbuf_space);''')
 for typ,key in [('uint16_t','x'),('uint8_t','c')]:
  for operation in [f'S.{key}.clear();',f'std::vector<{typ}>().swap(S.{key});']:
   if operation in s:
    s=s.replace('            '+operation, f'''            size_t const qsim_{key}_before = S.{key}.capacity();
            {operation}
            if (std::getenv("QSIM_DIAGNOSTIC"))
                fprintf(stderr, "QSIM {key.upper()} %zu %zu\\n", qsim_{key}_before, S.{key}.capacity());''')
 return s
for arm in ['U','R','B','RB']:
 d=a.output/('sources_'+arm);d.mkdir(exist_ok=True)
 for path,content in zip(paths,[r if 'R' in arm else u[paths[0]],b if 'B' in arm else u[paths[1]]]):
  if path==paths[1]:content=diag(content)
  (d/Path(path).name).write_text(content)
(a.output/'diagnostic_v1.patch').write_text(''.join(difflib.unified_diff(u[paths[1]].splitlines(True),diag(u[paths[1]]).splitlines(True),fromfile='a/'+paths[1],tofile='b/'+paths[1])))
# U is the first compiled arm.
for path in paths:(a.source/path).write_text((a.output/'sources_U'/Path(path).name).read_text())
