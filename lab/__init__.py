"""Experiment engine for the PPS/Walsh research project.

Extracted from the patterns the first 35 experiment scripts kept rewriting.
Core instruments (walsh, perm_pps, toffoli_arith, ...) stay at the repo root;
this package holds the machinery *around* them:

    harness   -- the protocol made executable: predictions declared before
                 measurement, must-fail controls, PASS/FAIL reporting
    measure   -- pullback support/density/peak with a content-addressed cache
    gf2       -- rank/kernel, the linear+affine structure finder, coset splits
    modarith  -- order, beta*2^alpha split, Carmichael lambda, real bit tables
    variants  -- overridable modexp builders and the wrap registry
    nulls     -- null models with the known sampling traps designed out

Finished experiment scripts under experiments/ are lab-notebook records and
do NOT get retrofitted onto this engine; it exists for the next experiment.
"""
from lab.harness import Experiment
from lab.measure import (pullback, support, sparsity, density, peak_pps,
                         fn_spectrum, fn_support)
from lab.gf2 import (rank_kernel, structures, find_structures,
                     check_structure, coset_split, quadrant_counts, slice_fn)
from lab.modarith import (order, v2_split, ord2, carmichael, bit_table, tile)
from lab.variants import (register_wrap, build_modexp, verify_correctness,
                          VariantModExp, WRAPS)
from lab.nulls import random_boolean, random_table, planted_structure
