"""
Greedy add-back for the ATP-sweep condition-specific failures.

For each single-failure probe (base-usable + 2/3 conditions usable + 1
infeasible):
  - Load the base at that hpf; apply clean probe bounds (class-1 closed,
    O2 at physiological OCR, substrate forced to -1).
  - SKO every reaction in base.rxns − cell.rxns to force the base into
    the cell's rxn-superset constraint.
  - Greedy add-back: at each step unlock the one reaction whose unlock
    produces the largest ATP gain. Stop when ATP > 0.5.
  - Record the ordered restoring set with GPR and per-gene zFPKM in the
    failing cell. Classify the probe:
      artefact        — every restoring rxn has no GPR OR all mapped
                        genes are above the exp threshold (model dropped
                        a reaction whose genes were available).
      transcriptional — at least one restoring rxn has mapped genes and
                        none are above threshold (model correctly
                        dropped an unexpressed reaction).

Writes atp_sweep_addback_probes.csv and atp_sweep_addback_summary.csv.
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from s250_atp_sweep import (
    load_model, class1_vitamin_exchanges, prep_bounds,
    V7_BASE_DIR, V7_EXTRACTED, ATP_RXN, O2_EX,
)
# Inline these to avoid re-running the clean sweep at import time.
OCR = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
def prep_bounds_clean(model, class1_set, probed_mar):
    lb, ub = prep_bounds(model, class1_set)
    for rid in class1_set:
        if rid == probed_mar: continue
        if rid in model["rxns"]:
            j = model["rxns"].index(rid)
            lb[j] = 0.0; ub[j] = 1000.0
    return lb, ub

OUT = REPO / "grace_v4/results/atp_sweep_v7"
CONDS = ("BL", "D", "LD")
MIN_ATP = 0.5
MAX_STEPS = 30

sweep = pd.read_csv(OUT / "atp_sweep_v7_clean.csv")
panel = pd.read_csv(OUT / "atp_sweep_substrates.csv")
diffs = pd.read_csv(OUT / "atp_sweep_v7_clean_condition_differences.csv")

zfp = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_scores_by_cell.csv")
zthr = pd.read_csv(REPO / "rfastcoreomics_trans/results/zfpkm_thresholds_by_cell.csv")
thr_map = dict(zip(zthr.iloc[:, 0], zthr["exp_threshold_zFPKM"]))
zfp = zfp.set_index(zfp.columns[0])

gmap = pd.read_csv(REPO / "results/qc/gene_map_ensdarg_to_model.csv")
SYM_TO_ENS = (gmap.groupby("model_gene")["ensembl_gene_id"]
              .apply(lambda s: sorted(set(s.dropna()))).to_dict())

class1 = class1_vitamin_exchanges()


def parse_gpr_genes(gpr):
    if not gpr: return []
    toks = re.findall(r"[A-Za-z0-9_.\-]+", gpr)
    return [t for t in toks if t.lower() not in ("and", "or", "not")]


def classify_rxn(gpr, cell, thr):
    """
    Returns (classification_tag, zfpkms_str, states_str) for one reaction.
    tag ∈ {no_gpr, all_expressed, none_expressed, mixed, no_mapping}
    """
    if gpr is None or (isinstance(gpr, float) and np.isnan(gpr)):
        return "no_gpr", "", ""
    g = str(gpr).replace("\x00", "").strip()
    if not g: return "no_gpr", "", ""
    gpr = g
    gpr = str(gpr).strip()
    symbols = parse_gpr_genes(gpr)
    vals, states = [], []
    for sym in symbols:
        ens = SYM_TO_ENS.get(sym, [])
        if not ens:
            vals.append(f"{sym}=NOMAP"); states.append("nomap"); continue
        zs = [float(zfp.loc[e, cell]) for e in ens if e in zfp.index and cell in zfp.columns]
        if not zs:
            vals.append(f"{sym}=NA"); states.append("na"); continue
        zmax = max(zs)
        vals.append(f"{sym}={zmax:.2f}")
        states.append("exp" if zmax >= thr else "noexp")
    exp_states = [s for s in states if s in ("exp", "noexp")]
    if not exp_states:
        return "no_mapping", ";".join(vals), ";".join(states)
    if all(s == "exp" for s in exp_states):
        return "all_expressed", ";".join(vals), ";".join(states)
    if all(s == "noexp" for s in exp_states):
        return "none_expressed", ";".join(vals), ";".join(states)
    return "mixed", ";".join(vals), ";".join(states)


def probe_lp_bounds(base, h, sub_mar, class1_set):
    lb, ub = prep_bounds_clean(base, class1_set, sub_mar)
    idx = {r: i for i, r in enumerate(base["rxns"])}
    lb[idx[sub_mar]] = -1.0; ub[idx[sub_mar]] = -1.0
    lb[idx[O2_EX]]   = -float(OCR[h]); ub[idx[O2_EX]] = 0.0
    return lb, ub, idx


base_cache = {}
def get_base(h):
    if h not in base_cache:
        base_cache[h] = load_model(V7_BASE_DIR / f"baked_ocr_{h}_v7.mat", "m")
    return base_cache[h]


cell_cache = {}
def get_cell(cell):
    if cell not in cell_cache:
        cell_cache[cell] = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
    return cell_cache[cell]


# ---- Build list of single-failure probes -----------------------------------
probes = []
for _, row in diffs.iterrows():
    flags = {c: bool(row[c]) for c in CONDS}
    if sum(flags.values()) != 2: continue
    h = int(row.hpf); sub_label = row.substrate
    br = sweep[(sweep.model == f"base_{h}") & (sweep.substrate == sub_label)]
    if br.empty or not bool(br.iloc[0].usable): continue
    failing = next(c for c, v in flags.items() if not v)
    probes.append((h, sub_label, f"{failing}_{h}"))

print(f"[info] single-failure probes: {len(probes)}", flush=True)


# ---- Greedy add-back engine ------------------------------------------------
all_detail_rows = []
summary_rows = []
for h, sub_label, cell in probes:
    print(f"\n[probe] hpf={h} sub={sub_label} failing={cell}", flush=True)
    sub = panel[panel.label == sub_label].iloc[0]
    base = get_base(h); cm = get_cell(cell)
    cell_rxns = set(cm["rxns"])
    lb0, ub0, idx = probe_lp_bounds(base, h, sub.mar, class1)
    # Pre-filter candidates to pFBA sparse route in the base: at each
    # step we only consider reactions that are on the parsimonious ATP
    # path. Reactions off this path cannot increase ATP under the probe
    # objective, so they are safely excluded from the greedy search.
    # This takes the per-probe search space from ~400 reactions down to
    # ~10-30 (the s253 "candidates" set).
    import scipy.sparse as sp
    S_mat = base["S"]; n_met, n_rxn = S_mat.shape
    env_p = gp.Env(empty=True); env_p.setParam("OutputFlag", 0); env_p.start()
    gm_p1 = gp.Model("fba", env=env_p)
    vp = gm_p1.addMVar(n_rxn, lb=lb0, ub=ub0)
    gm_p1.addMConstr(S_mat.tocsr(), vp, "=", np.zeros(n_met))
    cp = np.zeros(n_rxn); cp[idx[ATP_RXN]] = 1.0
    gm_p1.setMObjective(None, cp, 0.0, sense=GRB.MAXIMIZE)
    gm_p1.optimize()
    atp_star = gm_p1.ObjVal if gm_p1.Status == GRB.OPTIMAL else 0.0
    gm_p1.dispose()
    gm_p2 = gp.Model("pfba", env=env_p)
    vpp = gm_p2.addMVar(n_rxn, lb=0.0, ub=np.where(ub0 > 0, ub0, 0.0))
    vpn = gm_p2.addMVar(n_rxn, lb=0.0, ub=np.where(lb0 < 0, -lb0, 0.0))
    SS = sp.hstack([S_mat, -S_mat], format="csr")
    gm_p2.addMConstr(SS, gp.MVar.fromlist(list(vpp.tolist()) + list(vpn.tolist())),
                     "=", np.zeros(n_met))
    gm_p2.addConstr(vpp[idx[ATP_RXN]] - vpn[idx[ATP_RXN]] == atp_star)
    gm_p2.setMObjective(None, np.ones(2 * n_rxn), 0.0, sense=GRB.MINIMIZE)
    gm_p2.optimize()
    sparse_fluxes = np.array([vpp[j].X - vpn[j].X for j in range(n_rxn)])
    sparse_set = {base["rxns"][j] for j in range(n_rxn) if abs(sparse_fluxes[j]) > 1e-6}
    gm_p2.dispose(); env_p.dispose()
    print(f"  pFBA ATP*={atp_star:.2f}  sparse_route={len(sparse_set)}", flush=True)

    # Candidate pool = sparse route ∩ (base - cell)
    missing_all = [r for r in base["rxns"] if r not in cell_rxns]
    missing = [r for r in missing_all if r in sparse_set]
    print(f"  candidates (sparse ∩ missing-in-cell): {len(missing)} (of {len(missing_all)} missing)",
          flush=True)
    # SKO: zero EVERYTHING missing in cell (candidate pool + the rest), then
    # only the candidate pool is available to unlock in greedy.
    native_lb = {r: lb0[idx[r]] for r in missing}
    native_ub = {r: ub0[idx[r]] for r in missing}
    lb_cur = lb0.copy(); ub_cur = ub0.copy()
    for r in missing_all:
        j = idx[r]; lb_cur[j] = 0.0; ub_cur[j] = 0.0

    # Build persistent Gurobi model
    S = base["S"]; n_met, n_rxn = S.shape
    env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
    gm = gp.Model("addback", env=env)
    v = gm.addMVar(n_rxn, lb=lb_cur, ub=ub_cur)
    gm.addMConstr(S.tocsr(), v, "=", np.zeros(n_met))
    c = np.zeros(n_rxn); c[idx[ATP_RXN]] = 1.0
    gm.setMObjective(None, c, 0.0, sense=GRB.MAXIMIZE)
    gm.setParam("Method", 1)
    gm.optimize()
    atp0 = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    print(f"  cell-sim baseline (missing SKO'd): ATP = {atp0:.4f}", flush=True)
    if atp0 > MIN_ATP:
        print(f"  cell-sim already above {MIN_ATP}, probe not condition-specific infeasible under clean policy")
        summary_rows.append(dict(hpf=h, substrate=sub_label, failing_cell=cell,
                                 baseline_atp=atp0, final_atp=atp0, n_steps=0,
                                 restoring_set="", verdict="not_infeasible_under_clean"))
        gm.dispose(); env.dispose()
        continue

    # Fast probe: unlock ALL candidates at once. If ATP still < 0.5 the
    # pFBA pre-filter missed something; fall through to a looser search.
    for r in missing:
        j = idx[r]; v[j].lb = native_lb[r]; v[j].ub = native_ub[r]
    gm.optimize()
    atp_all_unlocked = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    print(f"  all pFBA-missing unlocked: ATP = {atp_all_unlocked:.4f}", flush=True)
    if atp_all_unlocked < MIN_ATP:
        # Fall through to full greedy over ALL missing (slow, used only for
        # probes the pFBA pre-filter couldn't close).
        for r in missing:
            j = idx[r]; v[j].lb = 0.0; v[j].ub = 0.0
        gm.optimize()
        missing = missing_all
        native_lb = {r: lb0[idx[r]] for r in missing}
        native_ub = {r: ub0[idx[r]] for r in missing}
        print(f"  [fallback] expanded candidate pool to {len(missing)}", flush=True)
    else:
        # Lock them back off; greedy will unlock in order of best gain.
        for r in missing:
            j = idx[r]; v[j].lb = 0.0; v[j].ub = 0.0
        gm.optimize()

    picks = []
    remaining = set(missing)
    step = 0
    cur = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
    while cur < MIN_ATP and remaining and step < MAX_STEPS:
        step += 1
        best_rxn = None; best_atp = cur
        for r in remaining:
            j = idx[r]
            ol, ou = v[j].lb, v[j].ub
            v[j].lb = native_lb[r]; v[j].ub = native_ub[r]
            gm.optimize()
            if gm.Status == GRB.OPTIMAL and gm.ObjVal > best_atp + 1e-9:
                best_atp = gm.ObjVal; best_rxn = r
            v[j].lb = ol; v[j].ub = ou
        if best_rxn is None:
            # No single reaction improves ATP individually, but all-at-once
            # unlock reaches ATP > 0.5 — the restoring set is a group of
            # reactions that must appear together (common for stoichiometric
            # couples: transporter + its reciprocal, or two halves of a
            # bidirectional carrier). Report the whole candidate pool as the
            # unordered restoring set with per-reaction classification.
            print(f"  step {step}: no single-rxn gain — reporting the pool as a group")
            for r in remaining:
                j = idx[r]; v[j].lb = native_lb[r]; v[j].ub = native_ub[r]
            gm.optimize()
            cur = gm.ObjVal if gm.Status == GRB.OPTIMAL else 0.0
            # Run pFBA on the fully-unlocked LP to find the actual active
            # restoring subset
            import scipy.sparse as sp
            lb_now = np.array([v[j].lb for j in range(n_rxn)])
            ub_now = np.array([v[j].ub for j in range(n_rxn)])
            env2 = gp.Env(empty=True); env2.setParam("OutputFlag", 0); env2.start()
            gm2 = gp.Model("pfba", env=env2)
            pp = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(ub_now > 0, ub_now, 0.0))
            pn = gm2.addMVar(n_rxn, lb=0.0, ub=np.where(lb_now < 0, -lb_now, 0.0))
            Sh = sp.hstack([S.tocsr(), -S.tocsr()], format="csr")
            gm2.addMConstr(Sh, gp.MVar.fromlist(list(pp.tolist()) + list(pn.tolist())),
                           "=", np.zeros(n_met))
            gm2.addConstr(pp[idx[ATP_RXN]] - pn[idx[ATP_RXN]] == cur)
            gm2.setMObjective(None, np.ones(2 * n_rxn), 0.0, sense=GRB.MINIMIZE)
            gm2.optimize()
            active_restoring = []
            for r in remaining:
                j = idx[r]
                f = pp[j].X - pn[j].X
                if abs(f) > 1e-6:
                    active_restoring.append((r, f))
            gm2.dispose(); env2.dispose()
            print(f"    group-restoring set (pFBA-active): {len(active_restoring)} rxns",
                  flush=True)
            for r, f in active_restoring:
                j = idx[r]
                gpr_raw = base["grRules"][j]
                gpr = "" if (gpr_raw is None or (isinstance(gpr_raw, float) and np.isnan(gpr_raw))) else str(gpr_raw).replace("\x00","").strip()
                subs_str = base["subs"][j] if j < len(base["subs"]) else ""
                tag, zfpkm_str, state_str = classify_rxn(gpr, cell, thr_map.get(cell, np.nan))
                picks.append(dict(step=step, rxn=r, atp_after=cur,
                                  gpr=gpr, subsystem=subs_str, flux=float(f),
                                  gene_zfpkms=zfpkm_str, gene_states=state_str,
                                  rxn_class=tag))
                print(f"       +{r:12s}  v={f:+7.3f}  cls={tag:14s}  gpr={gpr!r}",
                      flush=True)
            remaining.clear()
            break
        j = idx[best_rxn]
        v[j].lb = native_lb[best_rxn]; v[j].ub = native_ub[best_rxn]
        gm.optimize()
        cur = gm.ObjVal
        gpr_raw = base["grRules"][j]
        gpr = "" if (gpr_raw is None or (isinstance(gpr_raw, float) and np.isnan(gpr_raw))) else str(gpr_raw).strip()
        subs_str = base["subs"][j] if j < len(base["subs"]) else ""
        tag, zfpkm_str, state_str = classify_rxn(gpr, cell, thr_map.get(cell, np.nan))
        picks.append(dict(step=step, rxn=best_rxn, atp_after=cur,
                          gpr=gpr, subsystem=subs_str,
                          gene_zfpkms=zfpkm_str, gene_states=state_str,
                          rxn_class=tag))
        remaining.discard(best_rxn)
        print(f"  step {step}: +{best_rxn} ({subs_str[:35]}) → ATP={cur:.4f}  cls={tag}  gpr={gpr!r}",
              flush=True)
    gm.dispose(); env.dispose()

    # Probe-level verdict
    verdict = "artefact"
    reason = ""
    for p in picks:
        if p["rxn_class"] == "none_expressed":
            verdict = "transcriptional"
            reason = f"restoring rxn {p['rxn']} has mapped genes, all below threshold"
            break
    if not picks:
        verdict = "no_picks"; reason = "add-back could not reach 0.5"
    elif verdict == "artefact":
        # Determine sub-tag
        tags = [p["rxn_class"] for p in picks]
        if all(t == "no_gpr" for t in tags):
            reason = "all restoring reactions have no GPR (gene-less infrastructure)"
        elif all(t in ("no_gpr", "all_expressed") for t in tags):
            reason = "all restoring reactions either gene-less or had genes above threshold"
        else:
            reason = "mix; classified as artefact because no 'none_expressed' hit"
    summary_rows.append(dict(hpf=h, substrate=sub_label, failing_cell=cell,
                             baseline_atp=atp0, final_atp=cur, n_steps=len(picks),
                             restoring_set=";".join(p["rxn"] for p in picks),
                             verdict=verdict, reason=reason))
    for p in picks:
        p.update(dict(hpf=h, substrate=sub_label, failing_cell=cell))
        all_detail_rows.append(p)

detail = pd.DataFrame(all_detail_rows)
detail.to_csv(OUT / "atp_sweep_addback_detail.csv", index=False)
pd.DataFrame(summary_rows).to_csv(OUT / "atp_sweep_addback_summary.csv", index=False)
print(f"\n[write] atp_sweep_addback_detail.csv rows={len(detail)}", flush=True)
print(f"[write] atp_sweep_addback_summary.csv rows={len(summary_rows)}", flush=True)

# Summary verdict table
print("\n=== Probe verdicts ===")
print(pd.DataFrame(summary_rows)[["hpf","substrate","failing_cell","n_steps","baseline_atp","final_atp","verdict","reason"]].to_string())
