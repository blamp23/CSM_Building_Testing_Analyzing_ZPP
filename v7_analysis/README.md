# v7_analysis — active Stages 2, 4, 11

What's left here after the 2026-10-26 chapter reframe. Three active subdirs,
all still used:

```
preprocessing/
  fix_a_tightening.py       Stage 2 — LD_72 polytope conditioning (patched
                            2026-10-12 to fix a numerical degeneracy that
                            blocked polytope rounding on D_72 under default
                            PolyRound settings)

loop_removal/
  cff_gurobi.py             Stage 4 — CycleFreeFlux via Gurobi. Bit-equivalent
                            on L1 to the scipy-HiGHS reference implementation
                            of Desouki et al. 2015, ~2× faster. Warm-starts
                            between consecutive samples.

matched_floor/
  hopsy_sample_direct.py    Stage 3 (variant) — direct hopsy billiard-walk
                            sampler with explicit Chebyshev starting point
                            and feasibility shrink
  run_pipeline.sh           Stage 11 — matched-floor F2 resample for 9 cells
                            (local polytope prep + Grace sampling + CFF)
```

The W1 statistics (`statistics/`, `figures/`, `tier0/`, `positive_control/`,
`launchers/`, `plots/`) that used to live here are now in
`../legacy/v7_pipeline_deprecated/` — superseded by `v8_analysis/scripts/`.
