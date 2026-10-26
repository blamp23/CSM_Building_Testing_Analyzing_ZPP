# pipeline/05_chapter archive

104 scripts preserved for provenance. **Not** needed to reproduce the chapter —
the linear pipeline is `pipeline/05_chapter/scripts/01..09`.

Keep for forensics only. If you need something here, it's probably because you
want to see how a decision got made, not to run anything.

## Categories

**Sampling + loop removal (used by Stages 0–5, superseded by `pipeline/02_sampling/`
and `pipeline/03_loop_removal/`):**
`s41d_hopsy_sample_direct.py`, `s41d_aggregate.py`, `s51*_cyclefreeflux*.py`,
`s52_post_sample_subsets.py`, `s60_split_chain*.py`, `s70_cff_verify.py`

**W1 findings pipeline (superseded by chapter reframe — scripts 04–07):**
`s80_findings_v7*.py` (3 normaliser variants), `s81*_density_plots*.py` (5
density-plot iterations before `09_density_plots.py`), `s82*_positive_control*.py`,
`s83_fva_capability.py`, `s85_export_tier_tables.py`, `s90*/s91*/s92*_fig*.py` (W1
heatmap + condition-axis + developmental-axis figures),
`two_tier_subsystem_test*.py` (4 variants), `run_contrasts_*.sh` (launchers).

**Tier 0 structural analysis (now at `legacy/v7_pipeline_deprecated/tier0/`;
these were earlier iterations):** `s150_tier0_2x2_v7.py`, `s151_tier0_fixups_bl96.py`,
`s16{0..9}_*.py`, `s17{0..9}_*.py`, `s179_tables.py`.

**Lipid-relay 2×2 (PARKED, see `PARKED_lipid_relay_state_20261016.md`):**
`s180_f1_relay.py`, `s180_figures_3456.py`, `s181_f2_late_relay.py`,
`s181_fig7_unknown_zone.py`, `s184_lipid_relay_swap_sample.py`,
`s185_lipid_relay_baseline.py`, `s186_lipid_relay_decompose.py`,
`s187*_build_polytopes*`, `s188_grace_tarball.py`,
`s189_table5b_developmental.py`, `s190_reaction_footprint.py`,
`s191_footprint_figures.py`, `run_lipid_relay.slurm`.

**Grace resample remediation (one-off fixes for the matched-floor F2 round):**
`s192_cff_swap_corners.py`, `s193_retry_round_LD120_LD96.py`,
`s194*.py` (Rhat/ESS checks), `s20b_rank_rhat.py`.

**Pre-reframe chapter figures (replaced by `pipeline/05_chapter/scripts/08`/`09`):**
`s200_figures_v7*.py`, `s201_condition_axis_curves.py`,
`s202_run_matched_f2_pipeline.sh`.

**ATP-sweep iterations (final = `01_substrate_yields.py`; these preceded it):**
`s251_atp_sweep_figure.py`, `s252_atp_sweep_transporter_check.py`,
`s253_atp_sweep_pfba_essentials.py`, `s254_atp_sweep_pfba_figure.py`,
`s255_diagnose_excess_atp.py`, `s256_loopless_excess_atp.py`,
`s257_atp_sweep_clean.py`.

**Medium-energy gate build-up (preceded `02`/`03`):**
`s258_greedy_addback.py`, `s259_universal_lost_addback.py`,
`s260_anaerobic_ldh_gap.py`, `s261_medium_energy_gate.py`,
`s262_inventory_LD72.py`, `s263_oxphos_reconcile.py`,
`s264*_set_gate*.py`, `s265_sampled_uptake_C5.py`,
`s266_compute_demand_bounds.py`, `s267_tier0_before_after.py`.

**Collision siblings of lifted scripts (same number, different purpose):**
- `s268_fva_before_after.py` — pre-reframe FVA sanity check
- `s269_carotenoid_route.py` — carotenoid pathway probe (unused in chapter)
- `s270_LD72_old_vs_new_contrast.py` — contrast diff between model versions
- `s270_per_reaction_deseq.py` — Mann-Whitney variant superseded by `04`/`07`
- `s270_resample_pipeline.sh` — one-off resample driver
