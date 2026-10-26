"""
s147 -- v6 vs v7 sizes table (n_rxns, n_mets, n_genes per cell).
"""
import h5py, numpy as np, pandas as pd
from pathlib import Path

REPO = Path("/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models")
V6_DIR = REPO / "ocr_anchored_extraction/models"
V7_DIR = REPO / "reviewer_packet_v4_slim/v7_models"
CELLS = [f"{c}_{h}" for c in ("BL","D","LD") for h in (24,48,72,96,120)]


def sizes(fp, root):
    try:
        with h5py.File(fp, "r") as f:
            m = f[root]
            n_rxn = np.asarray(m["rxns"]).shape[-1]
            n_met = np.asarray(m["mets"]).shape[-1]
            try:
                n_gene = np.asarray(m["genes"]).shape[-1]
            except Exception:
                n_gene = np.nan
    except Exception as e:
        return dict(n_rxn=np.nan, n_met=np.nan, n_gene=np.nan, err=str(e))
    return dict(n_rxn=n_rxn, n_met=n_met, n_gene=n_gene)


rows = []
for cell in CELLS:
    v6_fp = V6_DIR / f"trans_rfastcormics_{cell}.mat"
    v7_fp = V7_DIR / f"trans_rfastcormics_{cell}.mat"
    s6 = sizes(v6_fp, "parent") if v6_fp.exists() else {"n_rxn": np.nan, "n_met": np.nan, "n_gene": np.nan}
    s7 = sizes(v7_fp, "parent") if v7_fp.exists() else {"n_rxn": np.nan, "n_met": np.nan, "n_gene": np.nan}
    rows.append(dict(cell=cell,
                     v6_rxn=s6["n_rxn"], v7_rxn=s7["n_rxn"], drxn=s7["n_rxn"]-s6["n_rxn"] if not np.isnan(s6["n_rxn"]) else np.nan,
                     v6_met=s6["n_met"], v7_met=s7["n_met"], dmet=s7["n_met"]-s6["n_met"] if not np.isnan(s6["n_met"]) else np.nan,
                     v6_gene=s6["n_gene"], v7_gene=s7["n_gene"], dgene=s7["n_gene"]-s6["n_gene"] if not np.isnan(s6["n_gene"]) else np.nan))

df = pd.DataFrame(rows)
df.to_csv(REPO/"reviewer_packet_v4_slim/s147_sizes_v6_v7.csv", index=False)
print(df.to_string(index=False))
