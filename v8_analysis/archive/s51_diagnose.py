"""
Diagnose CFF LP failure. Runs the exact CFF formulation from s51 on ONE sample
and prints the HiGHS status + message. Also measures the residual |S*v0|
that the LP is being asked to enforce exactly, and tries the fix (b_eq set
to S*v0 rather than zero, so the LP only needs to preserve what v0 already
satisfies).
"""
import argparse, time
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.optimize import linprog
from pathlib import Path

OCR_MAP = {24: 97, 48: 173, 72: 227, 96: 273, 120: 313}
BOUND_EX = 10
BIO_LB_FRAC = 0.9
O2_EX = "MAR09048"
BIOMASS_RXN = "MAR00021"

p = argparse.ArgumentParser()
p.add_argument("--cell", default="LD_72")
p.add_argument("--chain", type=int, default=1)
p.add_argument("--sample", type=int, default=0)
p.add_argument("--indir", default="hopsy_v4")
args = p.parse_args()

repo = Path.cwd()

def load_v73(mat_fp):
    with h5py.File(mat_fp, "r") as f:
        m = f["parent"]
        def deref_str(ref):
            out = []
            for r in np.asarray(ref).ravel():
                obj = f[r]; arr = np.asarray(obj).ravel()
                out.append("".join(chr(int(c)) for c in arr))
            return out
        S_grp = m["S"]
        data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
        return dict(
            S=sp.csc_matrix((data, ir, jc), shape=(int(ir.max())+1, len(jc)-1)),
            rxns=deref_str(m["rxns"]))

cell = args.cell
cond, hpf = cell.split("_"); hpf = int(hpf)
M = load_v73(repo/"models"/f"trans_rfastcormics_{cell}.mat")
S = M["S"]; rxn_ids = M["rxns"]; n_rxn = len(rxn_ids); n_met = S.shape[0]
idx = {r:i for i,r in enumerate(rxn_ids)}
Snnz = np.diff(S.tocsc().indptr); is_ex = Snnz == 1
j_bio = idx[BIOMASS_RXN]; ex_or_bio = is_ex.copy(); ex_or_bio[j_bio] = True
internal = ~ex_or_bio; n_int = int(internal.sum())
int_idx = np.where(internal)[0]

# load one sample
fp_in = repo/"results"/args.indir/f"flux_samples_{cell}_ch{args.chain}.csv"
df_in = pd.read_csv(fp_in, index_col=0).reindex(rxn_ids)
V_in = df_in.values.astype(float)
v0 = V_in[:, args.sample]
print(f"cell={cell} chain={args.chain} sample={args.sample}", flush=True)
print(f"  n_rxn={n_rxn}  n_met={n_met}  n_int={n_int}  n_ex_bio={ex_or_bio.sum()}", flush=True)

# 1. residual of v0
resid = S @ v0
print(f"  input residual |S*v0|_max        = {np.max(np.abs(resid)):.3e}", flush=True)
print(f"  input residual |S*v0|_median     = {np.median(np.abs(resid)):.3e}", flush=True)
print(f"  count |S*v0|_i > 1e-6            = {int((np.abs(resid) > 1e-6).sum())} / {n_met}", flush=True)
print(f"  count |S*v0|_i > 1e-4            = {int((np.abs(resid) > 1e-4).sum())} / {n_met}", flush=True)

# 2. distribution of v0 magnitudes
print(f"  max |v0|                         = {np.max(np.abs(v0)):.3g}", flush=True)
print(f"  count |v0| > 100                 = {int((np.abs(v0) > 100).sum())}", flush=True)
print(f"  count |v0| > 900                 = {int((np.abs(v0) > 900).sum())}", flush=True)

# ---- CURRENT LP (b_eq = 0) ----
def build_lp(v0, b_eq_full):
    n_var = n_rxn + 2*n_int
    c_obj = np.zeros(n_var); c_obj[n_rxn:n_rxn+n_int] = 1.0; c_obj[n_rxn+n_int:] = 1.0
    # S · v = b_eq_S  (first n_met rows); vp/vn identity rows follow
    S_coo = S.tocoo()
    rows = [S_coo.row]; cols = [S_coo.col]; data = [S_coo.data]
    for k, j in enumerate(int_idx):
        rows.append(np.array([n_met+k, n_met+k, n_met+k]))
        cols.append(np.array([j, n_rxn+k, n_rxn+n_int+k]))
        data.append(np.array([1.0, -1.0, 1.0]))
    A_eq = sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(n_met + n_int, n_var))
    b_eq = np.zeros(n_met + n_int)
    b_eq[:n_met] = b_eq_full
    lb_v = np.zeros(n_rxn); ub_v = np.zeros(n_rxn)
    for j in int_idx:
        if v0[j] >= 0: lb_v[j] = 0.0; ub_v[j] = max(v0[j], 0.0)
        else:          lb_v[j] = min(v0[j], 0.0); ub_v[j] = 0.0
    for j in np.where(ex_or_bio)[0]:
        lb_v[j] = v0[j]; ub_v[j] = v0[j]
    bounds = list(zip(lb_v, ub_v)) + [(0.0, None)]*n_int + [(0.0, None)]*n_int
    return c_obj, A_eq, b_eq, bounds

for tag, b_eq_S in [("CURRENT (b_eq=0)", np.zeros(n_met)),
                     ("FIX A (b_eq = S*v0 residual)", S @ v0)]:
    print(f"\n---- LP variant: {tag} ----", flush=True)
    c_obj, A_eq, b_eq, bounds = build_lp(v0, b_eq_S)
    t0 = time.time()
    res = linprog(c_obj, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                  method="highs", options={"disp": False,
                                            "presolve": True})
    dt = time.time() - t0
    print(f"  success       = {res.success}", flush=True)
    print(f"  status        = {res.status}", flush=True)
    print(f"  message       = {res.message}", flush=True)
    print(f"  fun (L1 min)  = {res.fun if res.success else 'N/A'}", flush=True)
    print(f"  wall          = {dt:.1f}s", flush=True)

    if res.success:
        v = res.x[:n_rxn]
        L1_before = float(np.abs(v0[internal]).sum())
        L1_after  = float(np.abs(v[internal]).sum())
        print(f"  L1 internal   : {L1_before:.4g} -> {L1_after:.4g} "
              f"({100*(1-L1_after/L1_before):.2f}% removed)", flush=True)
        exch_diff = float(np.max(np.abs(v[ex_or_bio] - v0[ex_or_bio])))
        print(f"  |v_ex - v0_ex|_max = {exch_diff:.2e}", flush=True)
        Sv_after = np.max(np.abs(S @ v))
        print(f"  |S*v|_max post = {Sv_after:.2e}", flush=True)

# ---- ALSO try relaxed equality (linprog no direct slack, so simulate with
# A_ub +/- tol). This is a sanity check on the tolerance angle.
print(f"\n---- LP variant: FIX B (b_eq = 0 but A_ub +/- 1e-4 slack on S*v) ----", flush=True)
n_var = n_rxn + 2*n_int
c_obj = np.zeros(n_var); c_obj[n_rxn:n_rxn+n_int] = 1.0; c_obj[n_rxn+n_int:] = 1.0
S_coo = S.tocoo()
# A_ub: S*v <= slack; -S*v <= slack  (i.e., -slack <= S*v <= slack)
Su = sp.csr_matrix((np.concatenate([S_coo.data, -S_coo.data]),
                    (np.concatenate([S_coo.row, S_coo.row + n_met]),
                     np.concatenate([S_coo.col, S_coo.col]))),
                    shape=(2*n_met, n_var))
b_ub = np.full(2*n_met, 1e-4)
# vp/vn identity as equality (needed for L1 formulation)
rows = []; cols = []; data = []
for k, j in enumerate(int_idx):
    rows.append(np.array([k, k, k]))
    cols.append(np.array([j, n_rxn+k, n_rxn+n_int+k]))
    data.append(np.array([1.0, -1.0, 1.0]))
A_eq_only = sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                           shape=(n_int, n_var))
b_eq_only = np.zeros(n_int)
lb_v = np.zeros(n_rxn); ub_v = np.zeros(n_rxn)
for j in int_idx:
    if v0[j] >= 0: lb_v[j] = 0.0; ub_v[j] = max(v0[j], 0.0)
    else:          lb_v[j] = min(v0[j], 0.0); ub_v[j] = 0.0
for j in np.where(ex_or_bio)[0]:
    lb_v[j] = v0[j]; ub_v[j] = v0[j]
bounds = list(zip(lb_v, ub_v)) + [(0.0, None)]*n_int + [(0.0, None)]*n_int
t0 = time.time()
res = linprog(c_obj, A_ub=Su, b_ub=b_ub, A_eq=A_eq_only, b_eq=b_eq_only,
              bounds=bounds, method="highs")
print(f"  success       = {res.success}", flush=True)
print(f"  status        = {res.status}", flush=True)
print(f"  message       = {res.message}", flush=True)
if res.success:
    v = res.x[:n_rxn]
    L1_before = float(np.abs(v0[internal]).sum())
    L1_after  = float(np.abs(v[internal]).sum())
    print(f"  L1 internal   : {L1_before:.4g} -> {L1_after:.4g} "
          f"({100*(1-L1_after/L1_before):.2f}% removed)", flush=True)
