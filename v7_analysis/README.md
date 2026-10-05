# v7_analysis — downstream pipeline from §11 amendments

This directory holds the analysis code added after the `v7-frozen` tag.
Everything here corresponds to a dated §11 amendment in
`FAMILY_PREREG.md`; nothing here changes the frozen model-building
pipeline (`base_model/`, `extraction/`, `validation/`) or the pre-
registered §1-§4 statistical rules — those remain the record of what
was pre-registered.

```
v7_analysis/
├── preprocessing/
│   └── fix_a_tightening.py       — §9 LD_72 polytope conditioning (s150)
├── loop_removal/
│   └── cff_gurobi.py             — CycleFreeFlux via Gurobi (s51b;
│                                    bit-equivalent on L1 to scipy-HiGHS
│                                    reference, ~2× faster)
├── statistics/
│   ├── two_tier_fast.py          — §2 pre-reg as originally written
│   │                               (mean normaliser + EPS-shift zeros)
│   ├── two_tier_fast_median.py   — §11.12 sensitivity (median + drop zeros)
│   ├── two_tier_fast_mean_dropzero.py
│   │                             — §11.14 METHOD OF RECORD
│   │                               (mean + drop zeros, matches v5)
│   ├── findings_v7.py            — §1–§4 + §11.9 (Pool excluded) +
│   │                               §11.10 floor-suspect +
│   │                               §11.13 SNR + envelope tracking
│   ├── findings_v7_median.py     — §11.12 variant of findings_v7
│   └── findings_v7_mean_dropzero.py
│                                 — method-of-record variant
├── plots/
│   ├── density_mean.py           — KDE overlays, mean normaliser
│   ├── density_median.py         — KDE overlays, median normaliser
│   ├── cdf_median.py             — empirical CDFs (tail-divergence view)
│   └── trajectory.py             — per-reaction trajectory plot across
│                                    BL / D / LD (clearest per-reaction view)
├── positive_control/
│   └── build.py                  — §11.15 corrected: builds LD_72 KD (10%
│                                    cap) and KO (lb=ub=0) models from the
│                                    Fix A tightened LD_72, re-FVA tightens
├── launchers/
│   ├── run_contrasts_mean_eps.sh      — fires 42 contrasts, §2 pre-reg
│   ├── run_contrasts_median.sh        — fires 42 contrasts, §11.12
│   └── run_contrasts_mean_dropzero.sh — fires 42 contrasts, §11.14
```

Expected directory layout at runtime (not in this repo):

```
grace_v4/
├── results/
│   ├── pooled_v7/{cell}_pooled.csv
│   ├── subsets_v7/{cell}_subset_membership.csv
│   ├── stats_v7/                        — §2 pre-reg
│   ├── stats_v7_median/                 — §11.12
│   ├── stats_v7_mean_dropzero/          — §11.14 method of record
│   ├── findings_v7/                     — §2 pre-reg + §11.9/10/13
│   ├── findings_v7_median/              — §11.12
│   └── findings_v7_mean_dropzero/       — §11.14 method of record
```

See `FAMILY_PREREG.md` for the dated amendment trail and
`v7_methods_for_reviewer.md` for the methods walkthrough.
