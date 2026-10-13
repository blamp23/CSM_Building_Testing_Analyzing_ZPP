"""
s85 -- export the final tier tables.
Method of record: mean + drop zeros (findings_v7_mean_dropzero).
Median variant used ONLY for the concordance column.
EPS variant NOT used (retracted per §11.14).

Outputs (overwritten):
  results/findings_v7_mean_dropzero/tierA_final.csv
  results/findings_v7_mean_dropzero/tierB_final.csv
  results/findings_v7_mean_dropzero/tierA_by_subsystem.csv
"""
import pandas as pd, numpy as np, re, h5py
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MOR  = ROOT/"results/findings_v7_mean_dropzero"
MED  = ROOT/"results/findings_v7_median"
MF96 = ROOT/"results/stats_v7_matched96"
MF2  = ROOT/"results/stats_v7_matched_f2"
MOD  = ROOT.parent/"reviewer_packet_v4_slim/v7_models"

# -------------------- load method-of-record findings
f1 = pd.read_csv(MOR/"findings_family1.csv"); f1['family']='F1'
f2 = pd.read_csv(MOR/"findings_family2.csv"); f2['family']='F2'
mor = pd.concat([f1,f2], ignore_index=True)

# -------------------- median variant for concordance
m1 = pd.read_csv(MED/"findings_family1.csv"); m1['family']='F1'
m2 = pd.read_csv(MED/"findings_family2.csv"); m2['family']='F2'
med = pd.concat([m1,m2], ignore_index=True)
median_tierA = med[(med.q_bh <= 0.10) & (med.tier == 'A')][['family','contrast','subsystem']].assign(median_concordant=True)
mor = mor.merge(median_tierA, on=['family','contrast','subsystem'], how='left')
mor['median_concordant'] = mor['median_concordant'].fillna(False)

# -------------------- matched-floor resolution for BL_96-involving (F1) contrasts
mat = {}
for pair in ['BL_96_vs_D_96', 'BL_96_vs_LD_96', 'D_96_vs_LD_96']:
    df = pd.read_csv(MF96/f"{pair}.csv")
    for _, r in df.iterrows(): mat[(pair, r.subsystem)] = (r.W1_obs, r.p_emp)
# -------------------- matched-floor resolution for F2 adjacent-hpf contrasts
mat_f2 = {}
if MF2.exists():
    for fp in sorted(MF2.glob("*_vs_*.csv")):
        pair = fp.stem
        df_m = pd.read_csv(fp)
        for _, r in df_m.iterrows():
            mat_f2[(pair, r.subsystem)] = (r.W1_obs, r.p_emp)
    print(f"[matched-f2] loaded {len(mat_f2)} (contrast × subsystem) matched-floor rows "
          f"from {len(list(MF2.glob('*_vs_*.csv')))} contrasts", flush=True)
def floor_resolution(r):
    # F1: use the §11.10 BL_96 matched-floor table (via floor_suspect)
    if r.family == 'F1':
        if not r.floor_suspect: return ""
        v = mat.get((r.contrast, r.subsystem))
        if v is None: return "not_tested"
        W1_m, p_m = v
        verdict = "survives" if p_m <= 0.05 else "artifact"
        return f"W1={W1_m:.3f}, p={p_m:.4f}, {verdict}"
    # F2: new 2026-10-13 matched-floor resample for adjacent-hpf contrasts
    if r.family == 'F2':
        v = mat_f2.get((r.contrast, r.subsystem))
        if v is None: return ""
        W1_m, p_m = v
        verdict = "survives" if p_m <= 0.05 else "artifact"
        return f"W1={W1_m:.3f}, p={p_m:.4f}, {verdict}"
    return ""
mor['floor_resolved'] = mor.apply(floor_resolution, axis=1)
# F2 matched-floor W1 and p as explicit columns on every F2 row (NaN if not tested).
def _mf2_W(r):
    if r.family != 'F2': return np.nan
    v = mat_f2.get((r.contrast, r.subsystem))
    return v[0] if v is not None else np.nan
def _mf2_p(r):
    if r.family != 'F2': return np.nan
    v = mat_f2.get((r.contrast, r.subsystem))
    return v[1] if v is not None else np.nan
mor['W1_matched_floor'] = mor.apply(_mf2_W, axis=1)
mor['p_matched_floor']  = mor.apply(_mf2_p, axis=1)

# -------------------- pre-correction retinoid flag (§11.17 pre-§11.19 bounds)
# Set from the subsystem list; survives any downstream in-place edits.
RETINOID_SUBS = {"Retinol metabolism", "Vitamin A metabolism",
                 "Carotenoid biosynthesis"}
mor['pre_correction_retinoid'] = mor.subsystem.isin(RETINOID_SUBS)

# -------------------- Tier 0 Frame A biomass_max per cell (child_biomax)
TIER0 = ROOT.parent/"reviewer_packet_v4_slim/s146_tier0_decomposition.csv"
bmax_t0 = pd.read_csv(TIER0).set_index('cell')['child_biomax'].to_dict()
print(f"[bmax Tier 0] loaded {len(bmax_t0)} cells. e.g. LD_48={bmax_t0.get('LD_48'):.2f}, LD_72={bmax_t0.get('LD_72'):.2f}", flush=True)

def floor_exposed(r):
    if r.family != 'F2': return ""
    bmA = bmax_t0.get(r.cell_A); bmB = bmax_t0.get(r.cell_B)
    if bmA is None or bmB is None: return ""
    ratio = max(bmA, bmB) / min(bmA, bmB)
    if ratio > 1.15:
        return f"yes (bmax ratio {ratio:.2f}: {r.cell_A}={bmA:.2f}, {r.cell_B}={bmB:.2f})"
    return "no"
mor['floor_exposed'] = mor.apply(floor_exposed, axis=1)

# -------------------- v5 reference
V5 = ROOT.parent/"reviewer_packet_v4_slim/s124_findings_annotated.csv"
if V5.exists():
    v5 = pd.read_csv(V5)
    def norm(c):
        m = re.match(r"([A-Z]+)(\d+)_vs_([A-Z]+)(\d+)", str(c))
        return f"{m.group(1)}_{m.group(2)}_vs_{m.group(3)}_{m.group(4)}" if m else str(c)
    v5['contrast'] = v5['contrast'].apply(norm)
    v5_set = set(zip(v5['contrast'], v5['subsystem']))
else:
    v5_set = set()
mor['v5_concordant'] = mor.apply(lambda r: (r.contrast, r.subsystem) in v5_set, axis=1)

# -------------------- parse contrast → cell_A, cell_B, hpf fields
def parse_contrast_full(r):
    m = re.match(r"([A-Z]+)_(\d+)_vs_([A-Z]+)_(\d+)", r.contrast)
    cA, hA, cB, hB = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    if r.family == 'F1':
        return pd.Series(dict(hpf=hA, hpf_A=hA, hpf_B=hB, condition=""))
    else:
        return pd.Series(dict(hpf="", hpf_A=hA, hpf_B=hB, condition=cA))
mor = mor.join(mor.apply(parse_contrast_full, axis=1))

# -------------------- p_type (emp vs gamma)
mor['p_type'] = np.where(mor.N_tier2.fillna(0) < 15, 'gamma', 'emp')

# -------------------- high_zero_frac flag (either cell > 0.25)
mor['high_zero_frac'] = (mor.zero_frac_A > 0.25) | (mor.zero_frac_B > 0.25)

# -------------------- null_median: pull from the raw stats CSVs
null_med = {}
for pair, cA, cB in [
    *[((f"{cA}_{h}_vs_{cB}_{h}"), f"{cA}_{h}", f"{cB}_{h}")
      for h in (24,48,72,96,120) for cA, cB in [("BL","D"),("BL","LD"),("D","LD")]],
    *[((f"{c}_{HPFS[i]}_vs_{c}_{HPFS[i+1]}"), f"{c}_{HPFS[i]}", f"{c}_{HPFS[i+1]}")
      for c in ("BL","D","LD") for i, HPFS in [(0,[24,48,72,96,120]),(1,[24,48,72,96,120]),
                                                (2,[24,48,72,96,120]),(3,[24,48,72,96,120])]]
]:
    fp = ROOT/f"results/stats_v7_mean_dropzero/{pair}.csv"
    if not fp.exists(): continue
    df = pd.read_csv(fp)
    for _, r in df.iterrows():
        null_med[(pair, r.subsystem)] = r.null_median
mor['null_median'] = mor.apply(lambda r: null_med.get((r.contrast, r.subsystem), np.nan), axis=1)

# -------------------- envelope_tracking (parse out triggering exchange + ratio)
def envelope_detail(r):
    if not r.envelope_tracking: return ""
    direction = r.direction
    thresh_hi, thresh_lo = 2.0, 0.5
    bits = []
    for item in str(r.ex_fc_summary).split(";"):
        if not item or ":" not in item: continue
        try:
            ex, rest = item.split(":", 1)
            parts = rest.split("=")
            ratio = float(parts[-1])
            if (direction > 0 and ratio >= thresh_hi) or (direction < 0 and ratio <= thresh_lo):
                bits.append(f"{ex}({rest})")
        except Exception:
            continue
    return "; ".join(bits) if bits else "flagged (details in ex_fc_summary)"
mor['envelope_trigger'] = mor.apply(envelope_detail, axis=1)

# -------------------- tier assignment (final headline)
# Tier A: passes every methodological filter + median concordance
# Tier B: passes q ≤ 0.10 + snr ≥ 2 but fails ≥ 1 Tier A filter (and is NOT envelope-tracked
#          or asymmetric-zeros, which would make it Tier C)
# Tier C: anything else that is still q ≤ 0.10 (flagged for interpretive caution)
sig = mor.q_bh <= 0.10
pass_q_snr2 = sig & (mor.snr >= 2)
not_flagged = (~mor.asymmetric_zeros) & (~mor.envelope_tracking)
mor['floor_safe'] = mor.apply(lambda r:
    (not r.floor_suspect) or ("survives" in str(r.floor_resolved)), axis=1)
tier_a_mask = (pass_q_snr2 & not_flagged & mor.floor_safe
               & (mor.n_sig_contrasts >= 2) & (mor.N_tier2 >= 15)
               & (mor.snr >= 3) & mor.median_concordant)
# Tier B explicitly: passes q + snr ≥ 2, NOT asymmetric_zeros, NOT envelope_tracking,
# but fails at least one Tier A requirement (snr<3, N<15, no replication, floor-unresolved,
# or median non-concordant).
tier_b_mask = pass_q_snr2 & not_flagged & ~tier_a_mask
mor['tier'] = np.where(tier_a_mask, 'A', np.where(tier_b_mask, 'B', 'C'))

# -------------------- replicated_in: other Tier A contrasts for same subsystem × family.
# Strict rule: a Tier A row must have ≥ 1 OTHER actual Tier A row for the same
# subsystem in the same family (not just n_sig_contrasts ≥ 2 across §1 eligibility
# — the other eligible contrasts may have been demoted by subsequent filters).
# Demote rows that fail this to Tier B.
tA = mor[mor.tier == 'A']
rep_lookup = tA.groupby(['family','subsystem'])[['contrast','direction']].apply(
    lambda g: list(zip(g.contrast, g.direction))).to_dict()
def replicated_in(r):
    if r.tier != 'A': return ""
    others = [f"{c}({'+' if d>0 else '-'})"
              for c, d in rep_lookup.get((r.family, r.subsystem), []) if c != r.contrast]
    return "; ".join(sorted(others))
mor['replicated_in'] = mor.apply(replicated_in, axis=1)
# Demote Tier A rows with no actual replicate in the same family
demote_mask = (mor.tier == 'A') & (mor.replicated_in == '')
if demote_mask.any():
    print(f"[demote] {int(demote_mask.sum())} Tier A rows had no actual same-family Tier A replicate; demoted to Tier B", flush=True)
mor.loc[demote_mask, 'tier'] = 'B'

# -------------------- §11.20 F2 matched-floor demotion
# Any F2 Tier A row whose floor_resolved verdict is "artifact" → Tier C with
# flag_reason = "floor_artifact". Survivors keep Tier A plus the matched-floor
# W1 and p columns.
mor['flag_reason'] = ""
f2_artifact_mask = ((mor.tier == 'A') & (mor.family == 'F2') &
                    mor.floor_resolved.str.contains("artifact", na=False))
n_demote = int(f2_artifact_mask.sum())
mor.loc[f2_artifact_mask, 'tier'] = 'C'
mor.loc[f2_artifact_mask, 'flag_reason'] = "floor_artifact"
if n_demote:
    print(f"[§11.20 demote] {n_demote} F2 Tier A rows failed matched-floor (p > 0.05); → Tier C 'floor_artifact'", flush=True)

# -------------------- final column order
cols = ['family', 'contrast', 'cell_A', 'cell_B', 'hpf', 'hpf_A', 'hpf_B', 'condition',
        'subsystem', 'N_rep', 'N_tier2', 'carry_A', 'carry_B', 'loop_flagged',
        'zero_frac_A', 'zero_frac_B', 'high_zero_frac', 'p_fisher_zeros',
        'W1_obs', 'null_median', 'z',
        'p_primary', 'p_type', 'q_bh',
        'noise_floor_W1', 'snr',
        'direction',
        'median_concordant',
        'envelope_tracking', 'envelope_trigger',
        'floor_suspect', 'floor_resolved', 'floor_exposed',
        'W1_matched_floor', 'p_matched_floor',
        'n_sig_contrasts', 'replicated_in',
        'v5_concordant', 'pre_correction_retinoid', 'tier', 'flag_reason']
# ensure every column exists
for c in cols:
    if c not in mor.columns: mor[c] = ""

# -------------------- flag_reason for Tier C rows (enumerate every failed gate)
def _tier_c_reason(r):
    if r.flag_reason == "floor_artifact": return "floor_artifact"
    if r.tier != 'C': return r.flag_reason
    rs = []
    if pd.isna(r.snr) or r.snr < 2:          rs.append("snr<2")
    if bool(r.asymmetric_zeros):             rs.append("asymmetric_zeros")
    if bool(r.envelope_tracking):            rs.append("envelope_tracking")
    if pd.notna(r.N_tier2) and r.N_tier2 < 15: rs.append("N_tier2<15")
    if not bool(r.floor_safe):               rs.append("floor_exposed")
    if pd.notna(r.n_sig_contrasts) and r.n_sig_contrasts < 2: rs.append("not_replicated")
    if not bool(r.median_concordant):        rs.append("not_median_concordant")
    return ", ".join(rs) if rs else "unknown"
# only compute for rows that are q ≤ 0.10 and tier == C
c_mask = (mor.q_bh <= 0.10) & (mor.tier == 'C')
mor.loc[c_mask, 'flag_reason'] = mor[c_mask].apply(_tier_c_reason, axis=1)

tier_a = mor[mor.tier == 'A'][cols].sort_values(['family','subsystem','contrast'])
tier_b = mor[mor.tier == 'B'][cols].sort_values(['family','subsystem','contrast'])
tier_c = mor[(mor.tier == 'C') & (mor.q_bh <= 0.10)][cols].sort_values(
    ['family','subsystem','contrast'])

tier_a.to_csv(MOR/"tierA_final.csv", index=False)
tier_b.to_csv(MOR/"tierB_final.csv", index=False)
tier_c.to_csv(MOR/"tierC_final.csv", index=False)
print(f"[write] tierA_final.csv  rows={len(tier_a)}", flush=True)
print(f"[write] tierB_final.csv  rows={len(tier_b)}", flush=True)
print(f"[write] tierC_final.csv  rows={len(tier_c)}", flush=True)

# -------------------- all_rows_in_family.csv (every row that entered the family)
all_cols = cols.copy()
# ensure tier is set to "" for rows that didn't make q ≤ 0.10
mor_all = mor.copy()
mor_all['tier'] = np.where(mor_all.q_bh <= 0.10, mor_all['tier'], "")
mor_all[all_cols].sort_values(['family','contrast','subsystem']).to_csv(
    MOR/"all_rows_in_family.csv", index=False)
print(f"[write] all_rows_in_family.csv  rows={len(mor_all)}", flush=True)

# -------------------- by-subsystem summary
def summarise(gr):
    hpfs = sorted(set(list(gr.hpf_A) + list(gr.hpf_B)))
    return pd.Series(dict(
        n_contrasts_F1=int((gr.family=='F1').sum()),
        n_contrasts_F2=int((gr.family=='F2').sum()),
        n_contrasts=len(gr),
        families="; ".join(sorted(set(gr.family))),
        hpf_list=", ".join(str(h) for h in hpfs),
        directions=", ".join(sorted(set(gr.direction.dropna().astype(int).astype(str)))),
        W1_min=float(gr.W1_obs.min()),
        W1_max=float(gr.W1_obs.max()),
        snr_min=float(gr.snr.min()),
        snr_max=float(gr.snr.max()),
    ))
by_sub = (tier_a.groupby('subsystem').apply(summarise)
          .sort_values('n_contrasts', ascending=False).reset_index())
by_sub.to_csv(MOR/"tierA_by_subsystem.csv", index=False)
print(f"[write] tierA_by_subsystem.csv  rows={len(by_sub)}  (sum of n_contrasts = {int(by_sub.n_contrasts.sum())})", flush=True)

# --- verification prints ---
print()
print(f"=== VERIFICATION ===")
print(f"Tier A rows:                                  {len(tier_a)}")
print(f"Tier B rows:                                  {len(tier_b)}")
print(f"Tier C (q ≤ 0.10, flagged):                  {int((mor.tier=='C').sum())}")
print(f"Tier A rows with empty replicated_in:        {(tier_a.replicated_in == '').sum()} "
      f"(should be 0 by Tier A rule n_sig_contrasts ≥ 2)")
print(f"Tier A rows with envelope_tracking=True:     {int(tier_a.envelope_tracking.sum())} "
      f"(should be 0)")
print(f"Tier A rows with null_median NaN:            {int(tier_a.null_median.isna().sum())} "
      f"(should be 0)")
print(f"Tier A rows with high_zero_frac=True:        {int(tier_a.high_zero_frac.sum())}")
print(f"Tier A rows with floor_exposed starting 'yes': {int(tier_a.floor_exposed.str.startswith('yes').fillna(False).sum())}")
print()
# Floor-exposed Tier A rows breakdown
fx = tier_a[tier_a.floor_exposed.str.startswith('yes').fillna(False)]
if len(fx):
    print(f"--- Floor-exposed F2 Tier A rows ({len(fx)}): pending matched-floor resample ---")
    print(fx[['family','contrast','subsystem','W1_obs','z','floor_exposed']].to_string(index=False))
