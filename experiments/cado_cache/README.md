# Reproduce the bounded CADO cache experiment

This package changes cache-construction allocation lifetimes at CADO revision
70354d7a8d54e985e46ca0fb6fb64d10716ba8bc. It contains the exact R (reader), B
(builder), and common opt-in diagnostic patches, build tools, original validated
driver, corrected timing driver and administrative provenance helpers. It needs
this research repository's `lab.harness`; pass that root explicitly as `--project`.
No artifact requires the original worker-directory name.

`U` means upstream behavior with common dormant diagnostics. Runtime allocator
settings and library hashes are pinned per reproduction. The published run used
`LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libtcmalloc.so.4`, all other listed allocator
variables unset. A different allocator/compiler is a new scoped reproduction;
do not claim its executable hashes equal the archived ones.

For replay using the existing immutable v1 evidence, pass its root to
`make_inputs_v2.py --prior-root PRIOR --output inputs.json`. The helper verifies
old source/executable hashes and payload identities, checks rw/cw dimensions and
records new sidecar hashes. Invoke `cache_timing_v2.py` with explicit
`--prior-root`, `--inputs-manifest`, `--inputs-sha256`, `--validation-sha256`,
`--output` (must not exist), and `--project`. Obtain the two requested hashes from
SHA256 of inputs.json and PRIOR/science_v1/validate_report_v1.json. Supply the
expected reviewed hashes to reproduce the exact archived run; do not disable
checks to make a mismatch pass.

For a fresh local reproduction:

1. Create a new RUNROOT and copy this package's files into it. Download and unpack
   the exact upstream tarball from
   https://codeload.github.com/cado-nfs/cado-nfs/tar.gz/70354d7a8d54e985e46ca0fb6fb64d10716ba8bc
   (archive SHA256 df86a1af2921ffbdf601066e8c824f0af40c2803bbcf43d319808fae338e32e0).
2. Use an isolated Python environment for CADO configuration, installing the
   exact `build_requirements_v1.txt` versions. Configure CMake with `MPI=0`,
   `-DCMAKE_BUILD_TYPE=Release -DENABLE_SHARED=OFF`,
   `-DBWC_GF2_ARITHMETIC_BACKENDS=b64 -DBWC_GF2_MATMUL_BACKENDS=bucket`,
   empty `BWC_GFP_ARITHMETIC_BACKENDS`/`BWC_GFP_MATMUL_BACKENDS`, and
   `-DPYTHON_EXECUTABLE` pointing into that environment. CADO supplies embedded
   fmt and gf2x. Save configuration logs and enforce the10-minute preparation
   window; no systemwide installation is needed.
3. Run `prepare_arms_v1.py --source SOURCE --output RUNROOT/arms`, followed by
   `build_arms_v1.py --source SOURCE --build BUILD --output RUNROOT`. This builds
   only build_matcache/bench_matcache, at2 jobs, saves separate static CADO arm
   executables, and hashes their resolved runtime libraries. Keep the allocator
   environment identical from here onward. Run
   `capture_reproduction_metadata.py --root RUNROOT --stage freeze`.
4. From the research root, run `cache_experiment_v1.py` with explicit
   `--source SOURCE --build BUILD --bin-root RUNROOT/bin --output RUNROOT/science_v1
   --project RESEARCH --phase validate`. Stop if any gate fails. This is the
   actual sparse-XOR/cache/reload/path/capacity/mutant validation, not a substitute
   kernel. Preserve its report and hash.
5. For a complete fresh recreation of the archived two-version history, run the
   same v1 driver with `--phase performance`. This creates precisely the three
   prescribed fixture files and baseline hashes, and retains the original20
   rows. Its wall values have the documented timeout-polling defect and are only
   coarse launcher intervals. Run `capture_reproduction_metadata.py --root RUNROOT
   --stage finalize` and then `make_inputs_v2.py` as above.
6. Run corrected `cache_timing_v2.py` into a fresh output directory using the
   newly frozen input/validation hashes. It reuses fixtures and passed validation;
   it does not rebuild, regenerate, or search for favorable parameters. It
   executes exactly20 fresh-child read/build/save measurements with4GiBAS,
   45sCPU,60swall child limits and a600s accounting ledger. The10s prior debit is
   administrative, not a certified bound on old interpreter setup.

Steps5–6 reproduce the preserved history deliberately. The current accepted
revision task reused existing v1 fixtures/reports and ran step6 only. Repository
research runs require their usual authorization, core gate, registered
predictions/run records and independent review. Stop on scientific failure;
never silently patch a frozen driver or overwrite an earlier output directory.

GNU time RSS is the CADO child high-water value; global page cache and Python
harness memory are additional. Corrected elapsed is launcher-inclusive blocking
wait, measured before watchdog cleanup; GNU-time elapsed/CPU are also retained
at0.01s print precision. These small synthetic matrices do not establish a whole
factorization memory or general speedup claim.
