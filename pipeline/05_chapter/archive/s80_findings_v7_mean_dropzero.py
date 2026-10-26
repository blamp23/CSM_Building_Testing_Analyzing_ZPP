"""
s80_findings_v7 -- Implements FAMILY_PREREG §1-§4 and §11 verbatim on
the v7 two-tier contrast CSVs.

Inputs:
    results/stats_v7/{A}_vs_{B}.csv     primary contrasts (Family 1 + 2)
    results/stats_v7/neg_{cell}.csv     negative control split-chain
    results/pooled_v7/{cell}_pooled.csv per-cell CFF'd samples (for zero
                                        fractions, direction, gamma-null
                                        regeneration)
    results/subsets_v7/{cell}_subset_membership.csv  per-cell subsystem map
    reviewer_packet_v4_slim/s124_findings_annotated.csv   v5 reference

Outputs:
    results/findings_v7/findings_family1.csv
    results/findings_v7/findings_family2.csv
    results/findings_v7/negcontrol_summary.csv
    results/findings_v7/negcontrol_per_subsystem.csv
    results/findings_v7/report.txt

Pipeline per row:
  §1 filter:  carry_both >= 5 AND loop_flagged/N_tier2 < 0.30 AND
              subsystem not in {Exchange, Artificial, Transport}.
  §2 p_prim: p_gamma if N_tier2 < 15, else p_emp. Gamma MLE on strictly-
             positive null W1 draws (regenerated from pooled samples +
             membership + rng(0)).
  §3 BH:     global within each family (not per contrast).
  §4 report: W1_obs, z, direction, asymmetric_zeros flag, p_fisher_zeros,
             neg_control_suspect flag, concordance vs v5.
  §11 tiers: A (mid-N, zero-symmetric, not neg-flagged, >= 2 contrasts),
             B (passes family+q but single-contrast or small-N),
             C (flagged).
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, gamma as gamma_dist
import re, sys, h5py, scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]   # grace_v4/
REPO = ROOT.parent                            # Discrete_Models/
STATS = ROOT/"results/stats_v7_mean_dropzero"  # §11.14 method of record
POOL  = ROOT/"results/pooled_v7"
SUBS  = ROOT/"results/subsets_v7"
MODELS = REPO/"reviewer_packet_v4_slim/v7_models"
FC_FP  = ROOT/"data/phase2_exchange_bounds_long.csv"
V5    = REPO/"reviewer_packet_v4_slim/s124_findings_annotated.csv"
OUT   = ROOT/"results/findings_v7_mean_dropzero"
OUT.mkdir(parents=True, exist_ok=True)

SNR_TIER_A = 3.0
SNR_TIER_B = 2.0
FC_ENVELOPE_THRESH = 2.0   # |log2(fc_A/fc_B)| >= 1 (= ratio >= 2 or <= 0.5)

HPFS   = [24, 48, 72, 96, 120]
CONDS  = ["BL", "D", "LD"]
CELLS  = [f"{c}_{h}" for c in CONDS for h in HPFS]

# Family 1: between-condition at each hpf (15 contrasts)
FAM1 = []
for h in HPFS:
    for cA, cB in [("BL","D"), ("BL","LD"), ("D","LD")]:
        FAM1.append((f"{cA}_{h}", f"{cB}_{h}"))

# Family 2: within-condition adjacent hpf (12 contrasts)
FAM2 = []
for c in CONDS:
    for i in range(len(HPFS)-1):
        FAM2.append((f"{c}_{HPFS[i]}", f"{c}_{HPFS[i+1]}"))

NEG_CELLS = CELLS  # 15 split-chain contrasts

EXCLUDE_SUBS = {"Exchange reactions", "Artificial reactions",
                "Transport reactions", "Exchange", "Artificial", "Transport",
                "Isolated", "Pool reactions"}
# §1 excludes Exchange/Artificial/Transport; §11.9 (2026-10-03) adds "Pool
# reactions" on the same grounds (pseudo-reactions assembling lumped
# substrate pools, not biochemical pathway members).

CARRY_TOL = 1e-4
LOOP_FRAC_MAX = 0.30
ZERO_TOL = 1e-10
ZERO_ASYM_DELTA = 0.02
NEG_Z_SUSP = 3.0
Q_THRESH = 0.10

# ---------------------------------------------------------------- fc + exchange helpers

def load_fc_table():
    """Return dict {(cond, hpf, ex_rxn): fc}."""
    df = pd.read_csv(FC_FP)
    out = {}
    for _, r in df.iterrows():
        out[(r.condition, int(r.hpf), r.ex_rxn)] = float(r.fc)
    return out

_SUB_EX_CACHE = {}

def _resolve_h5_str(f, obj):
    """Decode a MATLAB-saved char array from an h5py object."""
    arr = np.asarray(obj).ravel()
    try:
        return "".join(chr(int(c)) for c in arr)
    except Exception:
        return ""

def build_sub_exchange_map(cell):
    """For `cell`, return dict {subsystem: set(ex_rxn_ids touching sub's mets)}.
    Cached per cell."""
    if cell in _SUB_EX_CACHE:
        return _SUB_EX_CACHE[cell]
    mat_fp = MODELS/f"trans_rfastcormics_{cell}.mat"
    with h5py.File(mat_fp, "r") as f:
        m = f["parent"]
        S_grp = m["S"]
        data = np.asarray(S_grp["data"])
        ir   = np.asarray(S_grp["ir"])
        jc   = np.asarray(S_grp["jc"])
        n_rxn = len(jc) - 1
        n_met = int(ir.max()) + 1
        S = sp.csc_matrix((data, ir, jc), shape=(n_met, n_rxn))
        rxn_refs = np.asarray(m["rxns"]).ravel()
        rxns = [_resolve_h5_str(f, f[r]) for r in rxn_refs]
        # subSystems can be nested; resolve to a per-rxn string
        sub_refs = np.asarray(m["subSystems"]).ravel()
        subs = []
        for ref in sub_refs:
            obj = f[ref]
            arr = np.asarray(obj).ravel()
            if arr.dtype == object and len(arr) > 0:
                sub = _resolve_h5_str(f, f[arr[0]])
            else:
                sub = _resolve_h5_str(f, obj)
            subs.append(sub)
    Snnz  = np.diff(S.indptr)
    is_ex = Snnz == 1
    # Metabolite of each exchange column (its one nonzero row)
    ex_met = {}
    for j in np.where(is_ex)[0]:
        met = int(S.indices[S.indptr[j]])
        ex_met[rxns[j]] = met
    # For each subsystem: union of metabolites touched by any NON-exchange rxn
    sub_mets = {}
    for j in range(n_rxn):
        if is_ex[j]: continue
        s = subs[j]
        if not s: continue
        m_set = sub_mets.setdefault(s, set())
        m_set.update(S.indices[S.indptr[j]:S.indptr[j+1]].tolist())
    # For each subsystem, find exchanges whose single metabolite is in sub's mets
    sub_exchanges = {}
    for sub, mets in sub_mets.items():
        touching = {rxn for rxn, met in ex_met.items() if met in mets}
        sub_exchanges[sub] = touching
    _SUB_EX_CACHE[cell] = sub_exchanges
    return sub_exchanges

def compute_envelope_flag(cell_a, cell_b, subsystem, direction, fc_lookup):
    """Return (envelope_tracking_bool, ex_fc_summary_string).
    direction = +1 (A higher), -1 (B higher)."""
    cond_a, hpf_a = cell_a.split("_"); hpf_a = int(hpf_a)
    cond_b, hpf_b = cell_b.split("_"); hpf_b = int(hpf_b)
    exs_a = build_sub_exchange_map(cell_a).get(subsystem, set())
    exs_b = build_sub_exchange_map(cell_b).get(subsystem, set())
    all_exs = exs_a | exs_b
    tracked = False
    summary_parts = []
    for ex in sorted(all_exs):
        fcA = fc_lookup.get((cond_a, hpf_a, ex))
        fcB = fc_lookup.get((cond_b, hpf_b, ex))
        if fcA is None or fcB is None or fcA <= 0 or fcB <= 0:
            continue
        ratio = fcA / fcB
        summary_parts.append(f"{ex}:{fcA:.2f}/{fcB:.2f}={ratio:.2f}")
        if direction > 0 and ratio >= FC_ENVELOPE_THRESH:
            tracked = True
        elif direction < 0 and ratio <= 1.0 / FC_ENVELOPE_THRESH:
            tracked = True
    return tracked, ";".join(summary_parts)

# ---------------------------------------------------------------- helpers

def load_contrast_csv(fp):
    df = pd.read_csv(fp)
    # some rows have NaN W1/z/p (subsystems with too few carry_both to run tier 2);
    # keep them for the §1 filter audit but drop before BH.
    return df

def load_samples(cell):
    df = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
    return df

def load_membership(cell):
    m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
    m = m.drop_duplicates("representative").set_index("representative")
    return m

def exclude_structural(sub):
    s = str(sub).strip()
    if s in EXCLUDE_SUBS: return True
    low = s.lower()
    return any(key in low for key in ("exchange", "artificial", "transport",
                                       "isolated"))

from collections import OrderedDict
_PAIR_CACHE = OrderedDict()
_PAIR_CACHE_MAX = 3   # ~320 MB per pair; cap at ~1 GB total

def _prepare_pair(cell_a, cell_b, pooled_cache, mem_cache):
    """Precompute AA/AB/log_AA/log_AB/subsys/bg/carry_both/loop/common once
    per (cell_a, cell_b) pair. LRU-cached — hundreds of subsystems within
    the same contrast reuse the same arrays, but we evict old pairs so
    peak RAM stays bounded as we move through contrasts."""
    key = (cell_a, cell_b)
    if key in _PAIR_CACHE:
        _PAIR_CACHE.move_to_end(key)
        return _PAIR_CACHE[key]
    XA = pooled_cache[cell_a]; XB = pooled_cache[cell_b]
    mem = mem_cache[cell_a]
    common = XA.index.intersection(XB.index)
    AA = np.abs(XA.loc[common].values)
    AB = np.abs(XB.loc[common].values)
    medA, medB = np.median(AA, 1), np.median(AB, 1)
    carryA = medA > CARRY_TOL; carryB = medB > CARRY_TOL
    carry_both = carryA & carryB
    loop = (AA.max(1) > 100.0) | (AB.max(1) > 100.0)
    gpr = mem.has_gpr.reindex(common).fillna(0).values == 1
    subsys = mem.subsystem.reindex(common).fillna("UNMAPPED")
    exch = subsys.str.contains("Exchange|Artificial|Transport", case=False).values
    # §11.14 method of record: mean normaliser on NON-ZERO values, zeros dropped
    AB_concat = np.concatenate([AA, AB], axis=1)
    mask_nz = AB_concat >= ZERO_TOL
    with np.errstate(all="ignore"):
        masked = np.where(mask_nz, AB_concat, np.nan)
        mu = np.nanmean(masked, axis=1)
    mu_safe = np.where(np.isfinite(mu) & (mu > 0), mu, 1.0)
    mu_col = mu_safe[:, None]
    log_AA = np.where(AA >= ZERO_TOL,
                       np.log2(np.maximum(AA, ZERO_TOL) / mu_col), np.nan)
    log_AB = np.where(AB >= ZERO_TOL,
                       np.log2(np.maximum(AB, ZERO_TOL) / mu_col), np.nan)
    bg = np.where(gpr & ~exch & ~loop & carry_both)[0]
    cache = dict(common=common, AA=AA, AB=AB, subsys=subsys.values,
                 carry_both=carry_both, loop=loop, log_AA=log_AA, log_AB=log_AB,
                 bg=bg)
    _PAIR_CACHE[key] = cache
    while len(_PAIR_CACHE) > _PAIR_CACHE_MAX:
        _PAIR_CACHE.popitem(last=False)
    return cache

Q_GRID = np.linspace(0, 1, 201)

def regenerate_null_for_small_n(cell_a, cell_b, subsystem, n_tier2,
                                 pooled_cache, mem_cache, ndraw=2000):
    """Fresh rng(0) null for one subsystem — statistically equivalent to the
    two_tier replay, which also uses the same statistic over 2000 draws."""
    pc = _prepare_pair(cell_a, cell_b, pooled_cache, mem_cache)
    carry_both = pc["carry_both"]; loop = pc["loop"]
    log_AA = pc["log_AA"]; log_AB = pc["log_AB"]; bg = pc["bg"]
    idx_sub = np.where(pc["subsys"] == subsystem)[0]
    idx2 = idx_sub[carry_both[idx_sub] & ~loop[idx_sub]]
    if len(idx2) < 1 or len(bg) <= len(idx2):
        return None
    rng = np.random.default_rng(0)
    draws = np.stack([rng.choice(bg, len(idx2), replace=False) for _ in range(ndraw)])
    la = log_AA[draws].reshape(ndraw, -1)
    lb = log_AB[draws].reshape(ndraw, -1)
    with np.errstate(all="ignore"):
        qa = np.nanquantile(la, Q_GRID, axis=1)
        qb = np.nanquantile(lb, Q_GRID, axis=1)
    return np.abs(qa - qb).mean(axis=0)

def gamma_p(w1_obs, null):
    pos = null[null > 0]
    if len(pos) < 10:
        return np.nan, np.nan, np.nan
    shape, loc, scale = gamma_dist.fit(pos, floc=0)
    p = 1.0 - gamma_dist.cdf(w1_obs, shape, loc=0, scale=scale)
    return shape, scale, max(p, 1.0/(len(null)+1))

def compute_sub_extras(cell_a, cell_b, subsystem,
                       pooled_cache, mem_cache):
    """zero_frac_A, zero_frac_B, p_fisher_zeros, direction for one row."""
    pc = _prepare_pair(cell_a, cell_b, pooled_cache, mem_cache)
    AA = pc["AA"]; AB = pc["AB"]
    idx_sub = np.where(pc["subsys"] == subsystem)[0]
    if len(idx_sub) == 0:
        return np.nan, np.nan, np.nan, np.nan
    sub_A = AA[idx_sub]; sub_B = AB[idx_sub]
    zero_A = int((sub_A < ZERO_TOL).sum())
    zero_B = int((sub_B < ZERO_TOL).sum())
    total_A = sub_A.size; total_B = sub_B.size
    nonzero_A = total_A - zero_A; nonzero_B = total_B - zero_B
    zfA = zero_A / total_A
    zfB = zero_B / total_B
    try:
        _, p_fish = fisher_exact([[zero_A, nonzero_A], [zero_B, nonzero_B]])
    except Exception:
        p_fish = np.nan
    EPS = 1e-30
    log_A = np.log2(sub_A + EPS).mean()
    log_B = np.log2(sub_B + EPS).mean()
    direction = float(np.sign(log_A - log_B))   # +1: A higher, -1: B higher
    return zfA, zfB, p_fish, direction

def bh_fdr(p):
    """Benjamini-Hochberg q-values. NaN inputs pass through as NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    ok = ~np.isnan(p)
    pp = p[ok]
    if len(pp) == 0: return q
    order = np.argsort(pp)
    ranked = pp[order]
    n = len(ranked)
    bh = ranked * n / np.arange(1, n+1)
    # enforce monotonicity from the right
    for i in range(n-2, -1, -1):
        bh[i] = min(bh[i], bh[i+1])
    bh = np.clip(bh, 0, 1)
    out = np.full(n, np.nan)
    out[order] = bh
    q[ok] = out
    return q

# ---------------------------------------------------------------- v5 ref

def load_v5():
    if not V5.exists():
        return pd.DataFrame(columns=["contrast", "subsystem", "z", "q_bh"])
    v5 = pd.read_csv(V5)
    # v5 contrast uses "BL72_vs_LD72" style (no underscore between cond and
    # hpf). Normalize to our "BL_72_vs_LD_72" style.
    def norm(c):
        m = re.match(r"([A-Z]+)(\d+)_vs_([A-Z]+)(\d+)", str(c))
        if m: return f"{m.group(1)}_{m.group(2)}_vs_{m.group(3)}_{m.group(4)}"
        return str(c)
    v5["contrast_norm"] = v5["contrast"].apply(norm)
    return v5[["contrast_norm","subsystem","z","q_bh","robustness_tier"]].rename(
        columns={"contrast_norm":"contrast","z":"v5_z","q_bh":"v5_q",
                 "robustness_tier":"v5_tier"})

v5 = load_v5()
print(f"[v5] loaded {len(v5)} reference rows", flush=True)
fc_lookup = load_fc_table()
print(f"[fc] loaded {len(fc_lookup)} (cond,hpf,ex_rxn) fc entries", flush=True)

# ---------------------------------------------------------------- load caches

pooled_cache = {}
mem_cache = {}
needed_cells = set()
for a, b in FAM1 + FAM2: needed_cells.update([a, b])
needed_cells.update(NEG_CELLS)
print(f"[cache] loading {len(needed_cells)} pooled + membership CSVs...", flush=True)
for c in sorted(needed_cells):
    pooled_cache[c] = load_samples(c)
    mem_cache[c] = load_membership(c)
print(f"[cache] done. pooled shapes example: {next(iter(pooled_cache.values())).shape}", flush=True)

# ---------------------------------------------------------------- neg control first (per-subsystem absz lookup)

print(f"\n=== Negative control (15 split-chain contrasts) ===", flush=True)
neg_rows = []
for cell in NEG_CELLS:
    fp = STATS/f"neg_{cell}.csv"
    if not fp.exists():
        print(f"  WARN missing {fp}"); continue
    df = pd.read_csv(fp)
    df["cell"] = cell
    neg_rows.append(df)
neg_df = pd.concat(neg_rows, ignore_index=True)
# §1 filter for negs (same rules)
neg_filt = neg_df.copy()
neg_filt["pass_struct"] = ~neg_filt["subsystem"].apply(exclude_structural)
neg_filt["pass_carry"]  = neg_filt["carry_both"] >= 5
neg_filt["pass_loop"]   = (neg_filt["loop_flagged"] / neg_filt["N_tier2"].replace(0, np.nan)) < LOOP_FRAC_MAX
neg_fam = neg_filt[neg_filt.pass_struct & neg_filt.pass_carry & neg_filt.pass_loop & neg_filt["N_tier2"].notna()].copy()

# p_primary for negs: use p_emp if N>=15, else fit gamma
small_neg = neg_fam[neg_fam.N_tier2 < 15]
print(f"[neg] {len(neg_fam)} rows pass §1 filter; {len(small_neg)} need gamma fit", flush=True)
# Neg contrasts use halfA vs halfB of same cell. For gamma regeneration we
# use the FULL cell's pooled samples -- a conservative approximation since
# background size is slightly smaller for halves but subsystem index invariant.
neg_fam = neg_fam.sort_values("cell").reset_index(drop=True)
g_shape = np.full(len(neg_fam), np.nan); g_scale = np.full(len(neg_fam), np.nan)
p_gamma_ = np.full(len(neg_fam), np.nan)
last_cell = None
for k, row in neg_fam.iterrows():
    if row.cell != last_cell:
        last_cell = row.cell
        print(f"  [neg] processing {row.cell} ...", flush=True)
    if row.N_tier2 < 15:
        null = regenerate_null_for_small_n(row.cell, row.cell, row.subsystem,
                                            row.N_tier2, pooled_cache, mem_cache)
        if null is not None:
            shape, scale, p = gamma_p(row.W1_obs, null)
            g_shape[k] = shape; g_scale[k] = scale; p_gamma_[k] = p
neg_fam["gamma_shape"] = g_shape
neg_fam["gamma_scale"] = g_scale
neg_fam["p_gamma"] = p_gamma_

neg_fam["p_primary"] = np.where(neg_fam.N_tier2 < 15, neg_fam.p_gamma, neg_fam.p_emp)
neg_fam["q_bh"] = bh_fdr(neg_fam.p_primary.values)

# §11.13: noise_floor_W1[sub] = 95th pct W1_obs across the 15 split-chain
# contrasts for that subsystem. snr = W1_obs / noise_floor_W1 replaces the
# older |z| > 3 neg-control-suspect flag in tiering.
noise_floor_W1 = (neg_fam.groupby("subsystem")["W1_obs"]
                  .quantile(0.95)
                  .rename("noise_floor_W1"))
# Legacy neg-control-suspect kept for provenance (we still report it as a
# column), computed as max |z| across the 15 split-chain contrasts > 3.
neg_suspect = (neg_fam.assign(absz=neg_fam.z.abs())
               .groupby("subsystem")["absz"].max()
               .rename("neg_control_absz_max"))
suspect_set = set(neg_suspect[neg_suspect > NEG_Z_SUSP].index)
print(f"[§11.13] noise_floor_W1 computed for {len(noise_floor_W1)} subsystems "
      f"(median noise floor = {float(noise_floor_W1.median()):.3f})", flush=True)

neg_fam.to_csv(OUT/"negcontrol_per_subsystem.csv", index=False)
neg_sum = pd.DataFrame({
    "median_absz":   [neg_fam.z.abs().median()],
    "p95_absz":      [neg_fam.z.abs().quantile(0.95)],
    "n_rows":        [len(neg_fam)],
    "n_q_le_0.10":   [int((neg_fam.q_bh <= Q_THRESH).sum())],
    "empirical_FDR": [float((neg_fam.q_bh <= Q_THRESH).sum()) / max(len(neg_fam), 1)],
    "n_suspect_subs": [len(suspect_set)],
})
neg_sum.to_csv(OUT/"negcontrol_summary.csv", index=False)
print(neg_sum.to_string(index=False))
print(f"[neg] {len(suspect_set)} subsystems flagged neg_control_suspect (|z|>{NEG_Z_SUSP})", flush=True)

# ---------------------------------------------------------------- process a family

def process_family(contrasts, label):
    print(f"\n=== {label} ({len(contrasts)} contrasts) ===", flush=True)
    rows = []
    for a, b in contrasts:
        fp = STATS/f"{a}_vs_{b}.csv"
        if not fp.exists():
            print(f"  WARN missing {fp}"); continue
        df = load_contrast_csv(fp)
        df["cell_A"] = a; df["cell_B"] = b
        df["contrast"] = f"{a}_vs_{b}"
        rows.append(df)
    full = pd.concat(rows, ignore_index=True)

    # §1 filter
    full["pass_struct"] = ~full["subsystem"].apply(exclude_structural)
    full["pass_carry"]  = full["carry_both"] >= 5
    full["pass_loop"]   = (full["loop_flagged"] / full["N_tier2"].replace(0, np.nan)) < LOOP_FRAC_MAX
    full["pass_tier2"]  = full["N_tier2"].notna()
    full["in_family"]   = full.pass_struct & full.pass_carry & full.pass_loop & full.pass_tier2
    fam = full[full.in_family].copy()
    print(f"[{label}] total rows: {len(full)}; in family: {len(fam)}", flush=True)

    # Single pass per row: compute gamma (if small-N) AND zero_frac / direction.
    # Reusing _PAIR_CACHE across both means each (cell_A, cell_B) pair is
    # prepared only once per contrast block.
    small = fam[fam.N_tier2 < 15]
    print(f"[{label}] gamma fits needed for {len(small)} small-N rows; "
          f"extras for all {len(fam)} rows", flush=True)
    fam = fam.sort_values(["contrast"]).reset_index(drop=True)   # cluster by pair
    g_shape = np.full(len(fam), np.nan); g_scale = np.full(len(fam), np.nan)
    p_gamma_ = np.full(len(fam), np.nan)
    zfA_ = np.full(len(fam), np.nan); zfB_ = np.full(len(fam), np.nan)
    pfz_ = np.full(len(fam), np.nan); dir_ = np.full(len(fam), np.nan)
    last_contrast = None
    for k, row in fam.iterrows():
        if row.contrast != last_contrast:
            last_contrast = row.contrast
            print(f"  [{label}] processing {row.contrast} ...", flush=True)
        zfA, zfB, pfz, dr = compute_sub_extras(row.cell_A, row.cell_B,
                                                row.subsystem,
                                                pooled_cache, mem_cache)
        zfA_[k] = zfA; zfB_[k] = zfB; pfz_[k] = pfz; dir_[k] = dr
        if row.N_tier2 < 15:
            null = regenerate_null_for_small_n(row.cell_A, row.cell_B,
                                                row.subsystem, row.N_tier2,
                                                pooled_cache, mem_cache)
            if null is not None:
                shape, scale, p = gamma_p(row.W1_obs, null)
                g_shape[k] = shape; g_scale[k] = scale; p_gamma_[k] = p
    fam["gamma_shape"] = g_shape; fam["gamma_scale"] = g_scale
    fam["p_gamma"] = p_gamma_
    fam["zero_frac_A"] = zfA_; fam["zero_frac_B"] = zfB_
    fam["p_fisher_zeros"] = pfz_; fam["direction"] = dir_
    fam["asymmetric_zeros"] = (fam.zero_frac_A - fam.zero_frac_B).abs() > ZERO_ASYM_DELTA
    fam["neg_control_suspect"] = fam.subsystem.isin(suspect_set)   # legacy, provenance only
    # §11.10 floor-suspect flag (BL_96's bmax 40 % below D_96/LD_96 under Frame A)
    fam["floor_suspect"] = fam.contrast.str.contains("BL_96")
    # §11.13 SNR (envelope-tracking deferred until after q_bh)
    fam = fam.merge(noise_floor_W1, on="subsystem", how="left")
    fam["snr"] = fam.W1_obs / fam.noise_floor_W1

    fam["p_primary"] = np.where(fam.N_tier2 < 15, fam.p_gamma, fam.p_emp)

    # §3 BH-FDR global within family
    fam["q_bh"] = bh_fdr(fam.p_primary.values)

    # §11.13 envelope-tracking (needs q_bh for sig-row early-skip)
    env_track = np.zeros(len(fam), dtype=bool)
    env_summary = [""] * len(fam)
    for k, row in enumerate(fam.reset_index(drop=True).itertuples(index=False)):
        if not (row.q_bh <= Q_THRESH): continue
        tracked, summary = compute_envelope_flag(
            row.cell_A, row.cell_B, row.subsystem, row.direction, fc_lookup)
        env_track[k] = tracked
        env_summary[k] = summary
    fam["envelope_tracking"] = env_track
    fam["ex_fc_summary"] = env_summary

    # concordance vs v5
    if len(v5):
        fam = fam.merge(v5, on=["contrast","subsystem"], how="left")
        fam["found_in"] = np.where(fam.v5_z.notna() & (fam.q_bh <= Q_THRESH), "both",
                           np.where(fam.v5_z.notna() & ~(fam.q_bh <= Q_THRESH), "v5_only",
                           np.where(fam.v5_z.isna() & (fam.q_bh <= Q_THRESH), "v7_only",
                                    "neither")))
    else:
        fam["v5_z"] = np.nan; fam["v5_q"] = np.nan; fam["v5_tier"] = np.nan
        fam["found_in"] = np.where(fam.q_bh <= Q_THRESH, "v7_only", "neither")

    # §11.13 tiering based on snr + envelope + asymmetric_zeros + N_tier2 +
    # replication (q <= 0.10 still a precondition).
    sig = fam.q_bh <= Q_THRESH
    tier = np.full(len(fam), "none", dtype=object)
    flagged = fam.asymmetric_zeros | fam.envelope_tracking
    tier_c_mask = sig & (fam.snr.fillna(0) < SNR_TIER_B)        # snr < 2 always C
    tier_c_mask |= sig & flagged                                 # any flag → C
    tier_c_mask |= sig & fam.snr.isna()                          # no neg baseline → C
    tier[tier_c_mask] = "C"
    tier_b_mask = sig & (fam.snr >= SNR_TIER_B) & ~flagged       # candidate B/A
    tier[tier_b_mask] = "B"
    fam["tier"] = tier

    # Tier A promotion: snr >= 3, mid-N, zero-symmetric, not envelope-tracked,
    # AND replicated across >= 2 contrasts in SAME family under the same rules.
    tier_a_eligible = (sig & (fam.snr >= SNR_TIER_A) & (fam.N_tier2 >= 15)
                       & ~fam.asymmetric_zeros & ~fam.envelope_tracking
                       & fam.snr.notna())
    reps = (fam[tier_a_eligible]
            .groupby("subsystem").size().rename("n_sig_contrasts"))
    fam = fam.merge(reps, on="subsystem", how="left")
    fam.loc[tier_a_eligible & (fam.n_sig_contrasts >= 2), "tier"] = "A"
    fam["n_sig_contrasts"] = fam["n_sig_contrasts"].fillna(0).astype(int)

    out_fp = OUT/f"findings_{label.lower()}.csv"
    cols_order = ["contrast","cell_A","cell_B","subsystem","N_tier2","N_rep",
                  "carry_A","carry_B","carry_both","loop_flagged",
                  "W1_obs","z","direction",
                  "noise_floor_W1","snr",
                  "zero_frac_A","zero_frac_B","asymmetric_zeros",
                  "p_fisher_zeros",
                  "p_emp","p_gamma","gamma_shape","gamma_scale",
                  "p_primary","q_bh","tier","n_sig_contrasts",
                  "envelope_tracking","ex_fc_summary",
                  "neg_control_suspect","floor_suspect",
                  "found_in","v5_z","v5_q","v5_tier"]
    fam = fam.reindex(columns=[c for c in cols_order if c in fam.columns])
    fam.to_csv(out_fp, index=False)

    n_sig = int((fam.q_bh <= Q_THRESH).sum())
    print(f"[{label}] {n_sig} rows with q<=0.10. Tiers: "
          f"A={int((fam.tier=='A').sum())}, "
          f"B={int((fam.tier=='B').sum())}, "
          f"C={int((fam.tier=='C').sum())}")
    return fam

fam1 = process_family(FAM1, "Family1")
fam2 = process_family(FAM2, "Family2")

# ---------------------------------------------------------------- report

tierA_1 = fam1[fam1.tier == "A"].sort_values("z", ascending=False)
tierA_2 = fam2[fam2.tier == "A"].sort_values("z", ascending=False)

with open(OUT/"report.txt", "w") as f:
    f.write("v7 two-tier findings report\n")
    f.write("="*60 + "\n\n")
    f.write(f"Family 1 (between-condition at hpf):  {len(FAM1)} contrasts\n")
    f.write(f"Family 2 (within-condition adj hpf):  {len(FAM2)} contrasts\n")
    f.write(f"Neg control (split-chain):            {len(NEG_CELLS)} contrasts\n\n")
    f.write("--- Family 1 ---\n")
    f.write(f"  rows passing §1 filter: {len(fam1)}\n")
    f.write(f"  q <= {Q_THRESH}:          {int((fam1.q_bh <= Q_THRESH).sum())}\n")
    f.write(f"  Tier A (headline):      {int((fam1.tier=='A').sum())}\n")
    f.write(f"  Tier B:                 {int((fam1.tier=='B').sum())}\n")
    f.write(f"  Tier C:                 {int((fam1.tier=='C').sum())}\n\n")
    f.write("--- Family 2 ---\n")
    f.write(f"  rows passing §1 filter: {len(fam2)}\n")
    f.write(f"  q <= {Q_THRESH}:          {int((fam2.q_bh <= Q_THRESH).sum())}\n")
    f.write(f"  Tier A (headline):      {int((fam2.tier=='A').sum())}\n")
    f.write(f"  Tier B:                 {int((fam2.tier=='B').sum())}\n")
    f.write(f"  Tier C:                 {int((fam2.tier=='C').sum())}\n\n")
    f.write("--- Neg control ---\n")
    f.write(neg_sum.to_string(index=False) + "\n\n")
    f.write("--- Tier A, Family 1 ---\n")
    if len(tierA_1):
        f.write(tierA_1[["contrast","subsystem","N_tier2","W1_obs","z","direction",
                         "q_bh","n_sig_contrasts","found_in","v5_z"]]
                .to_string(index=False) + "\n")
    else:
        f.write("(none)\n")
    f.write("\n--- Tier A, Family 2 ---\n")
    if len(tierA_2):
        f.write(tierA_2[["contrast","subsystem","N_tier2","W1_obs","z","direction",
                         "q_bh","n_sig_contrasts","found_in","v5_z"]]
                .to_string(index=False) + "\n")
    else:
        f.write("(none)\n")

print("\n" + "="*60)
print(open(OUT/"report.txt").read())
