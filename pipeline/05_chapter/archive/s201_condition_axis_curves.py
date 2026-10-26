"""
s201 — condition-axis CDF document, same-hpf comparisons on same page.

One page per (subsystem, hpf) that has ANY F1 Tier A or B row.
Page layout (landscape 14×6):
  Left margin : subsystem name + hpf
  Right       : 3 sub-panels side-by-side, one per F1 contrast at that hpf
                (BL vs D, BL vs LD, D vs LD). Each sub-panel has density
                strip on top + CDF below (shared x). Annotation box with
                W1 / snr / q / dir / tier / null_med / zf; shaded "W" area
                between the contrasted pair. If a contrast has no tiered
                row at that hpf, draw it anyway with "not significant".

Also emits two Tier A reader tables:
  tierA_table.pdf     visual heatmap-style table of Tier A hits
  tierA_table.md      compact, sorted, human-scannable markdown
"""
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import patches as mpatches

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

ROOT = Path(__file__).resolve().parents[2]
MOR  = ROOT / "grace_v4/results/findings_v7_mean_dropzero"
POOL = ROOT / "grace_v4/results/pooled_v7"
SUBS = ROOT / "grace_v4/results/subsets_v7"
OUT  = ROOT / "grace_v4/results/figures_v7"
OUT.mkdir(parents=True, exist_ok=True)

CCOL = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
CONDS = ["BL", "D", "LD"]
ZERO_TOL = 1e-10

F1_PAIRS_PER_HPF = [("BL","D"), ("BL","LD"), ("D","LD")]

allf = pd.read_csv(MOR/"all_rows_in_family.csv")
f1 = allf[allf.family=="F1"].copy()
# F1 contrasts are {cond}_{hpf}_vs_{cond}_{hpf} with same hpf on both sides
f1["hpf"] = f1.cell_A.str.split("_").str[1].astype(int)

tA = pd.read_csv(MOR/"tierA_final.csv")
tB = pd.read_csv(MOR/"tierB_final.csv")

f1_tiered = f1[f1.tier.isin(["A","B"])].copy()
# Unique (subsystem, hpf) targets, with "any_A" flag
targets = (f1_tiered.assign(any_A=(f1_tiered.tier=="A"))
           .groupby(["subsystem","hpf"], as_index=False)
           .agg(any_A=("any_A","max")))
targets["hpf"] = targets["hpf"].astype(int)
targets_tierA_first = targets.sort_values(["any_A","hpf","subsystem"],
                                          ascending=[False, True, True])
print(f"Target pages: {len(targets_tierA_first)} (any_A first)")

pool_cache = {}
def pool(cell):
    if cell not in pool_cache:
        pool_cache[cell] = pd.read_csv(POOL/f"{cell}_pooled.csv", index_col=0)
    return pool_cache[cell]
mem_cache = {}
def mem(cell):
    if cell not in mem_cache:
        m = pd.read_csv(SUBS/f"{cell}_subset_membership.csv")
        mem_cache[cell] = m.drop_duplicates("representative").set_index("representative")
    return mem_cache[cell]

SHORT = {
    "Phenylalanine, tyrosine and tryptophan biosynthesis": "Aromatic AA",
    "Valine, leucine, and isoleucine metabolism": "BCAA",
    "Arginine and proline metabolism": "Arg/Pro",
    "Glycolysis / Gluconeogenesis": "Glycolysis",
    "Fatty acid activation (cytosolic)": "FA activation (cyt)",
    "Fatty acid activation (endoplasmic reticular)": "FA activation (ER)",
}
def short(s):
    return SHORT.get(s, re.sub(r"\s+metabolism$","",str(s)))


def kde_ish(samples, xs, bw=0.35):
    s = samples[np.isfinite(samples)]
    if len(s) == 0:
        return np.zeros_like(xs)
    if len(s) > 20000:
        s = np.random.default_rng(0).choice(s, 20000, replace=False)
    diffs = (xs[:, None] - s[None, :]) / bw
    k = np.exp(-0.5 * diffs**2)
    return k.mean(axis=1) / (bw * np.sqrt(2*np.pi))


def build_curves(subsystem, hpf):
    cells = {c: f"{c}_{hpf}" for c in CONDS}
    rep_sets = []
    for c in CONDS:
        m = mem(cells[c])
        rep_sets.append(set(m.index[m.subsystem == subsystem]) & set(pool(cells[c]).index))
    common = sorted(set.intersection(*rep_sets))
    if len(common) < 3:
        return None
    Xs = {}
    for c in CONDS:
        Xs[c] = np.abs(pool(cells[c]).loc[common].values)
    all_concat = np.concatenate([Xs[c] for c in CONDS], axis=1)
    mu_r = all_concat.mean(axis=1)
    mu_r = np.where(mu_r <= 0, 1e-30, mu_r)
    log_curves = {}; zf = {}
    for c in CONDS:
        X = Xs[c]
        X_log = np.log2(X / mu_r[:, None] + 1e-30)
        flat = X_log.ravel()
        zf[c] = (X.ravel() < ZERO_TOL).mean()
        log_curves[c] = flat[X.ravel() >= ZERO_TOL]
    return dict(curves=log_curves, zf=zf, n_rep=len(common))


def contrast_meta(subsystem, hpf, cA, cB):
    """Return dict with W1/snr/q/dir/tier/null_med if row exists, else all NaN."""
    row = f1[(f1.subsystem == subsystem) & (f1.cell_A == cA) & (f1.cell_B == cB)]
    if len(row) == 0:
        return None
    r = row.iloc[0]
    # pull null_median from the tier-final tables if row is tiered
    nm = np.nan
    if r.tier == "A":
        hit = tA[(tA.subsystem==subsystem) & (tA.contrast==r.contrast)]
        if len(hit): nm = float(hit.iloc[0]["null_median"])
    elif r.tier == "B":
        hit = tB[(tB.subsystem==subsystem) & (tB.contrast==r.contrast)]
        if len(hit): nm = float(hit.iloc[0]["null_median"])
    return dict(W1=r.W1_obs, snr=r.snr, q=r.q_bh, direction=r.direction,
                tier=r.tier if r.tier in ("A","B","C") else "", null_med=nm)


def draw_sub(ax_dens, ax_cdf, cv, cA, cB, meta, show_legend=False):
    xs_pool = np.concatenate(list(cv["curves"].values()))
    xlo, xhi = np.percentile(xs_pool, [1, 99]) if len(xs_pool) else (-5, 5)
    xlo = max(xlo, -30); xhi = min(xhi, 30)
    xs = np.linspace(xlo, xhi, 400)

    # density strip
    for c in CONDS:
        if len(cv["curves"][c]) == 0: continue
        d = kde_ish(cv["curves"][c], xs)
        ax_dens.plot(xs, d, color=CCOL[c], lw=1.0)
    ax_dens.set_yticks([]); ax_dens.set_xticklabels([])
    ax_dens.set_xlim(xlo, xhi)
    for sp in ("top","right","left"): ax_dens.spines[sp].set_visible(False)
    ax_dens.spines["bottom"].set_lw(0.5)
    ax_dens.axvline(0, ls=":", color="#888", lw=0.5)
    ax_dens.tick_params(bottom=False)

    # CDF
    y_c = {}
    for c in CONDS:
        if len(cv["curves"][c]) == 0: continue
        srt = np.sort(cv["curves"][c])
        ys = np.searchsorted(srt, xs, side="right") / len(srt)
        ax_cdf.plot(xs, ys, color=CCOL[c], lw=1.3, label=c)
        y_c[c] = ys
    if cA in y_c and cB in y_c:
        ax_cdf.fill_between(xs, y_c[cA], y_c[cB], color="#E88A8A", alpha=0.75, lw=0)
    ax_cdf.set_xlim(xlo, xhi); ax_cdf.set_ylim(0, 1)
    ax_cdf.axvline(0, ls=":", color="#888", lw=0.5)
    ax_cdf.axhline(0.5, ls=":", color="#CCC", lw=0.5)
    ax_cdf.set_xlabel("log₂(|v|/μ_r)", fontsize=8)
    ax_cdf.set_ylabel("CDF", fontsize=8)
    for sp in ("top","right"): ax_cdf.spines[sp].set_visible(False)
    ax_cdf.tick_params(labelsize=7)

    # annotation box
    if meta is None:
        anno = "no data in family"
    else:
        arr = "↑" if meta["direction"]>0 else "↓"
        tier_lbl = f"Tier {meta['tier']}" if meta["tier"] else "not significant"
        nm_str = f"null med={meta['null_med']:.3f}" if np.isfinite(meta['null_med']) else ""
        anno = (f"W1={meta['W1']:.2f}  snr={meta['snr']:.1f}  q={meta['q']:.3f}\n"
                f"dir={arr}  {tier_lbl}  {nm_str}\n"
                f"zf " + "  ".join(f"{c} {cv['zf'][c]:.2f}" for c in CONDS))
    ax_cdf.text(0.02, 0.98, anno, ha="left", va="top", transform=ax_cdf.transAxes,
                fontsize=7,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
    ax_cdf.set_title(f"{cA} vs {cB}", fontsize=9, loc="left")
    if show_legend:
        ax_cdf.legend(fontsize=7, loc="lower right", frameon=False)


def render_pages(tgt_df, pdf_fp):
    with PdfPages(pdf_fp) as pdf:
        for _, t in tgt_df.iterrows():
            sh = short(t.subsystem)
            hpf = int(t.hpf)
            cv = build_curves(t.subsystem, hpf)
            fig = plt.figure(figsize=(14, 6))
            gs = fig.add_gridspec(2, 4, width_ratios=[0.8, 3, 3, 3],
                                  height_ratios=[1, 4],
                                  wspace=0.18, hspace=0.03,
                                  left=0.03, right=0.98, top=0.92, bottom=0.09)
            ax_lbl = fig.add_subplot(gs[:, 0]); ax_lbl.set_axis_off()
            ax_lbl.text(0.95, 0.55, f"{sh}\n{hpf} hpf",
                        ha="right", va="center", fontsize=11, color="#222",
                        transform=ax_lbl.transAxes)
            if cv is None:
                ax = fig.add_subplot(gs[:, 1:]); ax.set_axis_off()
                ax.text(0.5, 0.5, f"n_common_reps < 3 at hpf {hpf}", ha="center", va="center",
                        fontsize=10, transform=ax.transAxes)
                pdf.savefig(fig); plt.close(fig); continue
            for i, (cA_cond, cB_cond) in enumerate(F1_PAIRS_PER_HPF):
                cA = f"{cA_cond}_{hpf}"; cB = f"{cB_cond}_{hpf}"
                ax_dens = fig.add_subplot(gs[0, 1+i])
                ax_cdf  = fig.add_subplot(gs[1, 1+i], sharex=ax_dens)
                meta = contrast_meta(t.subsystem, hpf, cA, cB)
                draw_sub(ax_dens, ax_cdf, cv, cA, cB, meta, show_legend=(i==2))
            fig.suptitle(f"F1 — {sh} @ {hpf} hpf   (n_reps={cv['n_rep']})",
                         fontsize=10, y=0.98)
            pdf.savefig(fig); plt.close(fig)
    print(f"[write] {pdf_fp.name}  ({len(tgt_df)} pages)")


render_pages(targets_tierA_first, OUT/"condition_axis_curves.pdf")

# =========================================================================
# Tier A reader table — visual PDF + compact markdown
# =========================================================================
tA_short = tA.copy()
tA_short["short"] = tA_short.subsystem.apply(short)
# derive hpf robustly from contrast ("cond_h_vs_cond_h" for F1 or "cond_h1_vs_cond_h2" for F2)
tA_short["hpf"] = tA_short.cell_A.str.split("_").str[1].astype(int)
tA_short["arrow"] = tA_short.direction.apply(lambda d: "↑" if d>0 else "↓")

# Pivot visual: rows = subsystems (sorted by n_contrasts desc), cols = contrast
all_contrasts_f1 = sorted(set(tA_short[tA_short.family=="F1"].contrast),
                          key=lambda c: (int(c.split("_")[1]), c))
all_contrasts_f2 = sorted(set(tA_short[tA_short.family=="F2"].contrast),
                          key=lambda c: (c.split("_")[0], int(c.split("_")[1])))
all_contrasts = all_contrasts_f1 + all_contrasts_f2

sub_order = (tA_short.groupby("short").size()
             .sort_values(ascending=False).index.tolist())

# Build matrix of cell text (arrow + W1) and background color
n_r = len(sub_order); n_c = len(all_contrasts)
cell_txt = np.full((n_r, n_c), "", dtype=object)
cell_bg  = np.full((n_r, n_c), "#FFFFFF", dtype=object)
for _, r in tA_short.iterrows():
    i = sub_order.index(r["short"]); j = all_contrasts.index(r.contrast)
    cell_txt[i, j] = r.arrow
    cell_bg[i, j] = "#C8E0F4" if r.direction>0 else "#F4D0C8"

# Visual table PDF (landscape)
fig_h = max(5, 0.28 * n_r + 2)
fig_w = max(14, 0.5 * n_c + 3)
fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_axis_off()
# header row
col_hdr = [c.replace("_vs_", "\nvs\n") for c in all_contrasts]
row_hdr = sub_order
# draw
col_w = 1.0; row_h = 1.0
x_off = 3.0; y_off = n_r + 2
# column headers (grouped by family)
for j, c in enumerate(all_contrasts):
    fam = "F1" if c in all_contrasts_f1 else "F2"
    bg = "#DDE8F3" if fam=="F1" else "#F3EADD"
    ax.add_patch(mpatches.Rectangle((x_off+j*col_w, y_off), col_w, 1.5,
                                     facecolor=bg, edgecolor="black", lw=0.3))
    ax.text(x_off+j*col_w+col_w/2, y_off+0.75, c.replace("_vs_", "\nvs\n"),
            ha="center", va="center", fontsize=5.5, rotation=0)
# family labels
ax.text(x_off + len(all_contrasts_f1)*col_w/2, y_off+1.9, "Family 1 (between-condition)",
        ha="center", va="bottom", fontsize=9, fontweight="bold", color="#1F4E79")
ax.text(x_off + len(all_contrasts_f1)*col_w + len(all_contrasts_f2)*col_w/2, y_off+1.9,
        "Family 2 (within-condition, adjacent hpf)",
        ha="center", va="bottom", fontsize=9, fontweight="bold", color="#7F5416")
# rows
for i, sub in enumerate(row_hdr):
    y = y_off - 1 - i*row_h
    # row header
    ax.add_patch(mpatches.Rectangle((0, y), x_off, row_h, facecolor="#F8F8F8", edgecolor="black", lw=0.3))
    ax.text(x_off-0.1, y+row_h/2, sub, ha="right", va="center", fontsize=7)
    for j, c in enumerate(all_contrasts):
        bg = cell_bg[i,j]; txt = cell_txt[i,j]
        ax.add_patch(mpatches.Rectangle((x_off+j*col_w, y), col_w, row_h,
                                         facecolor=bg, edgecolor="black", lw=0.3))
        if txt:
            ax.text(x_off+j*col_w+col_w/2, y+row_h/2, txt, ha="center", va="center",
                    fontsize=9, color="#000")
ax.set_xlim(-0.1, x_off + n_c*col_w + 0.5)
ax.set_ylim(-0.5, y_off + 2.5)
ax.set_aspect("equal")

# legend
leg_y = -1.5
ax.add_patch(mpatches.Rectangle((x_off, leg_y), 0.6, 0.6, facecolor="#C8E0F4", edgecolor="black", lw=0.3))
ax.text(x_off+0.7, leg_y+0.3, "↑ earlier-letter / earlier-hpf cond HIGHER",
        va="center", fontsize=7)
ax.add_patch(mpatches.Rectangle((x_off+6, leg_y), 0.6, 0.6, facecolor="#F4D0C8", edgecolor="black", lw=0.3))
ax.text(x_off+6.7, leg_y+0.3, "↓ earlier-letter / earlier-hpf cond LOWER",
        va="center", fontsize=7)

fig.suptitle("Tier A findings — subsystem × contrast grid", fontsize=11, y=0.98)
fig.savefig(OUT/"tierA_table.pdf", bbox_inches="tight")
fig.savefig(OUT/"tierA_table.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"[write] tierA_table.pdf + png  ({n_r} subsystems × {n_c} contrasts)")

# Markdown: compact, grouped by subsystem
lines = ["# Tier A findings — compact reader",
         "",
         f"**62 rows across {n_r} unique subsystems and {n_c} unique contrasts.**",
         "",
         "Arrow ↑ = earlier-letter condition (F1) or earlier hpf (F2) is HIGHER; ↓ = LOWER.",
         ""]
for sub in sub_order:
    rows = tA_short[tA_short["short"]==sub].sort_values(["family","contrast"])
    lines.append(f"## {sub}  (n={len(rows)})")
    for _, r in rows.iterrows():
        lines.append(f"- {r.family}  {r.contrast}  {r.arrow}  "
                     f"W1={r.W1_obs:.2f}  snr={r.snr:.1f}  q={r.q_bh:.3f}")
    lines.append("")
(OUT/"tierA_table.md").write_text("\n".join(lines))
print(f"[write] tierA_table.md")

print("\nDONE")
