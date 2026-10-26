"""
B1: Under final apply_context_bounds for LD_72, list every exchange with
lb < 0. Columns: MAR, met_name, carbon, ledger class / action, lb, in_105.
Explicit check: FAD (MAR01939), PE-LD pool (MAR00656), and each MEDIUM
row (esp MAR13066 carbonate, MAR13072 Mg2+, MAR13073 Cu2+).
"""
import sys, re
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "functions"))
sys.path.insert(0, str(REPO / "grace_v4/scripts"))
from apply_context_bounds import apply_context_bounds, MEDIUM
from s250_atp_sweep import load_model, V7_EXTRACTED

OUT = REPO / "grace_v4/results/atp_sweep_v7"

cell = "LD_72"; cond = "LD"; hpf = 72
m = load_model(V7_EXTRACTED / f"trans_rfastcormics_{cell}.mat", "parent")
S = m["S"]
is_ex = (np.diff(S.indptr) == 1)
lb, ub, idx, _ = apply_context_bounds(m["lb"].copy(), m["ub"].copy(),
                                     m["rxns"], is_ex, cond, hpf)

led = pd.read_csv(REPO / "reviewer_packet_v4_slim/s130_v7_ledger.csv")
led_class = dict(zip(led.exchange_MAR_id.dropna(), led["class"]))
led_action = dict(zip(led.exchange_MAR_id.dropna(), led.action))

phase2 = pd.read_csv(REPO / "grace_v4/data/phase2_exchange_bounds_long.csv")
measured_105 = set(phase2.ex_rxn.unique())

def count_C(f):
    if not isinstance(f, str): return 0
    mm = re.search(r"C(\d*)(?![a-z])", f)
    if not mm: return 0
    return int(mm.group(1)) if mm.group(1) else 1

rows = []
for j in np.where(is_ex)[0]:
    if lb[j] >= 0: continue
    rid = m["rxns"][j]
    col = S.getcol(j).tocoo()
    met_i = int(col.row[0])
    rows.append(dict(
        mar=rid,
        met_name=m["metNames"][met_i],
        met_mam=m["mets"][met_i],
        formula=m["metFormulas"][met_i],
        carbon=count_C(m["metFormulas"][met_i]),
        ledger_class=led_class.get(rid, ""),
        ledger_action=led_action.get(rid, ""),
        is_MEDIUM=rid in MEDIUM,
        in_105_measured=rid in measured_105,
        lb=lb[j], ub=ub[j],
    ))

df = pd.DataFrame(rows).sort_values(["is_MEDIUM","in_105_measured","ledger_action","lb"],
                                    ascending=[False, True, True, True])
df.to_csv(OUT / "inventory_LD72_openUptake.csv", index=False)
print(f"[write] inventory_LD72_openUptake.csv  n_open_uptake={len(df)}", flush=True)
print(f"  of which: in-105 measured = {int(df.in_105_measured.sum())}")
print(f"  of which: MEDIUM           = {int(df.is_MEDIUM.sum())}")
print(f"  of which: ledger class 1   = {int((df.ledger_class==1).sum())}")

# Spot checks
for check in ["MAR13066","MAR13072","MAR13073","MAR01939","MAR00656",
              "MAR09276","MAR09361","MAR09152","MAR09154","MAR09147"]:
    row = df[df.mar == check]
    if row.empty:
        print(f"  SPOT: {check:10s} not open (lb>=0 or not in model)")
    else:
        r = row.iloc[0]
        print(f"  SPOT: {check:10s} {r.met_name:25s} C={r.carbon:3d} "
              f"lb={r.lb:+8.1f} action={r.ledger_action:20s} MEDIUM={r.is_MEDIUM} meas={r.in_105_measured}")

# MEDIUM breakdown explicitly
print("\n=== MEDIUM rows explicitly: ===")
med_rows = df[df.is_MEDIUM]
for _, r in med_rows.iterrows():
    print(f"  {r.mar}  {r.met_name:30s}  formula={r.formula:12s}  C={r.carbon:3d}  lb={r.lb:+8.1f}")
if len(med_rows) < len(MEDIUM):
    missing = MEDIUM - set(med_rows.mar)
    print(f"  [note] {len(missing)} MEDIUM rows have lb >= 0 in LD_72 (pruned or non-negative):")
    for mid in sorted(missing):
        if mid in idx:
            j = idx[mid]
            col = S.getcol(j).tocoo(); met_i = int(col.row[0])
            print(f"    {mid}  {m['metNames'][met_i]}  lb={lb[j]}  ub={ub[j]}  (pruned={mid not in m['rxns']})")
        else:
            print(f"    {mid}  (pruned from this model)")

# Ledger-closed rows that are mysteriously open?
print("\n=== Ledger-closed rows with lb < 0 (potential bug) ===")
closed_set = set(led.loc[led.action == "close_exchange", "exchange_MAR_id"].dropna())
viol = df[df.mar.isin(closed_set)]
if viol.empty:
    print("  (none — all ledger-closed rows are closed, as expected)")
else:
    print(viol.to_string())
