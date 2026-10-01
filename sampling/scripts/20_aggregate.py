"""
s41d_aggregate -- After all chains for a cell are written (via --chain_only),
                  run the 4 assertions + compute Rhat/ESS + write diagnostics.
"""
import argparse, pickle
import numpy as np, pandas as pd
from pathlib import Path
from scipy.optimize import linprog

p = argparse.ArgumentParser()
p.add_argument("--cell", default="LD_72")
p.add_argument("--n_chains", type=int, default=4)
p.add_argument("--outdir_name", default="hopsy_v4")
args = p.parse_args()

repo = Path.cwd()
outdir = repo/"results"/args.outdir_name
meta_fp = outdir/f"{args.cell}_polytope_meta.pkl"
cache_fp = outdir/f"{args.cell}_polytope_cache.pkl"

with open(meta_fp,"rb") as fh: meta = pickle.load(fh)
rxn_ids = meta["rxn_ids"]; idx = meta["idx"]
n_rxn = meta["n_rxn"]; n_met = meta["n_met"]
lb = meta["lb"]; ub = meta["ub"]
bio_max_full = meta["bio_max_full"]; bio_lb = meta["bio_lb"]; frame = meta["frame"]

# Load all chain samples
chains = []
for ch in range(1, args.n_chains+1):
    fp = outdir/f"flux_samples_{args.cell}_ch{ch}.csv"
    if not fp.exists():
        raise FileNotFoundError(fp)
    df = pd.read_csv(fp, index_col=0)
    chains.append(df.reindex(rxn_ids).values.T)   # (n_samples, n_rxn)
full = np.stack(chains, axis=0)   # (n_chains, n_samples, n_rxn)
n_ch, n_smp, _ = full.shape
print(f"[{args.cell}] loaded {n_ch} chains x {n_smp} samples x {n_rxn} rxns")

# Rebuild S for assertions
import h5py, scipy.sparse as sp
mat_fp = repo/"models"/f"trans_rfastcormics_{args.cell}.mat"
with h5py.File(mat_fp,"r") as f:
    S_grp = f["parent"]["S"]
    data = np.asarray(S_grp["data"]); ir=np.asarray(S_grp["ir"]); jc=np.asarray(S_grp["jc"])
    S = sp.csc_matrix((data, ir, jc), shape=(n_met, n_rxn))

# ---- assertions ----
print("=== ASSERTIONS on 100 random samples ===")
rng = np.random.default_rng(0)
flat = full.reshape(-1, n_rxn)
check_idx = rng.choice(flat.shape[0], size=min(100, flat.shape[0]), replace=False)
X = flat[check_idx]

# 1. mass balance
massbal_max = float(np.abs(X @ S.toarray().T).max())
print(f"[1] |S·v|.max() = {massbal_max:.2e}  ({'PASS' if massbal_max < 1e-6 else 'FAIL'})")

# 2. bounds
lb_viol = float(np.max(lb[None,:] - X - 1e-6, initial=-np.inf))
ub_viol = float(np.max(X - ub[None,:] - 1e-6, initial=-np.inf))
print(f"[2] max lb violation = {lb_viol:.3e}, max ub violation = {ub_viol:.3e}  "
      f"({'PASS' if max(lb_viol,ub_viol) < 1e-4 else 'FAIL'})")

# 3. biomass
bio = X[:, idx["MAR00021"]]
bio_min = float(bio.min())
print(f"[3] biomass min={bio_min:.4f}, mean={bio.mean():.4f}, need >= {bio_lb:.4f}  "
      f"({'PASS' if bio_min >= bio_lb - 1e-4 else 'FAIL'})")

# 4. biomass_max reduced vs full
with open(cache_fp,"rb") as fh: polytope = pickle.load(fh)
T = polytope.transformation.values; shift = polytope.shift.values
c_eff = T[idx["MAR00021"]]; c_const = shift[idx["MAR00021"]]
lp = linprog(-c_eff, A_ub=polytope.A.values, b_ub=polytope.b.values,
             bounds=[(None,None)]*len(c_eff), method="highs")
bio_max_red = float(-lp.fun + c_const) if lp.success else float("nan")
diff = abs(bio_max_full - bio_max_red)
print(f"[4] biomass_max full={bio_max_full:.4f}, reduced={bio_max_red:.4f}, "
      f"diff={diff:.2e}  ({'PASS' if diff < 1e-2 else 'FAIL'})")

# ---- Rhat + ESS ----
half = n_smp // 2
X2 = np.concatenate([full[:, :half, :], full[:, half:, :]], axis=0)
M2 = X2.shape[0]; N2 = X2.shape[1]
chain_mean = X2.mean(axis=1)
grand_mean = chain_mean.mean(axis=0)
B = N2 * ((chain_mean - grand_mean) ** 2).sum(axis=0) / (M2 - 1)
W_ = ((X2 - chain_mean[:, None, :]) ** 2).sum(axis=(0, 1)) / (M2 * (N2 - 1))
rhat = np.full(n_rxn, np.nan)
ok = W_ > 1e-12
rhat[ok] = np.sqrt(((N2 - 1) / N2 * W_[ok] + B[ok] / N2) / W_[ok])
z = (full - full.mean(axis=1, keepdims=True)) / (full.std(axis=1, keepdims=True) + 1e-12)
lag1 = (z[:, :-1, :] * z[:, 1:, :]).mean(axis=1).mean(axis=0)
ess_bulk = n_ch * n_smp * np.clip(1 - lag1, 1e-6, None) / np.clip(1 + lag1, 1e-6, None)

finite = np.isfinite(rhat)
lines = [
    f"Cell: {args.cell}   frame: {frame}",
    f"Chains: {n_ch}   samples/chain: {n_smp}",
    f"Non-fixed reactions: {finite.sum()} / {n_rxn}",
    f"biomass_max_full = {bio_max_full:.4f}",
    f"biomass_max_reduced = {bio_max_red:.4f}",
    f"biomass_lb enforced = {bio_lb:.4f}",
    f"biomass min across all samples = {full[:,:,idx['MAR00021']].min():.4f}",
    f"biomass mean across all samples = {full[:,:,idx['MAR00021']].mean():.4f}",
    "",
    f"Rhat < 1.05: {(rhat[finite] < 1.05).mean()*100:.1f}%",
    f"Rhat < 1.10: {(rhat[finite] < 1.10).mean()*100:.1f}%",
    f"Rhat > 1.20: {(rhat[finite] > 1.20).mean()*100:.1f}%",
    f"ESS bulk median: {np.nanmedian(ess_bulk):.1f}",
    f"ESS bulk 10th %: {np.nanquantile(ess_bulk, .1):.1f}",
]
diag = "\n".join(lines)
print("\n" + diag)
(outdir/f"{args.cell}_diagnostics.txt").write_text(diag)
pd.DataFrame({"rxn": rxn_ids, "rhat": rhat, "ess_bulk": ess_bulk, "lag1": lag1}) \
    .to_csv(outdir/f"{args.cell}_rhat_ess.csv", index=False)
print(f"\nwrote {outdir/f'{args.cell}_diagnostics.txt'}")
