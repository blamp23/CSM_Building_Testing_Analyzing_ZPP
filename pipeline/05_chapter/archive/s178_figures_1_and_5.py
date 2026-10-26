"""
s178 -- Chapter Figures 1 and 5 (revision).

Fixes applied vs the first draft:
  - Fig 5B now uses post-§11.19 corrected numbers from
    tier0_C3_2x2_before_after.csv (not the uncorrected tier0_2x2_v7.csv).
    BL_96–D_96 shifts from topology-dominated (85%) to mixed (~55/45).
    BL_96–LD_96 stays topology-dominated.
  - Fig 1A tier B/C labels now match §2.12 definitions.
  - Fig 1B labels the three columns KO / KD(FVA-max — did not fire) /
    KD2(sampled-max — fired, W1=2.29).
  - Fig 1C caption: 976 within-cell split-chain tests, 5 at q ≤ 0.10,
    empirical FDR 0.5 %.
  - More horizontal gutter in both figures; shorter panel titles.
"""
import sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
R = REPO / "grace_v4/results"
OUT = R / "chapter_figures"
OUT.mkdir(parents=True, exist_ok=True)

COND_COLOR = {"BL": "#1F77B4", "D": "#111111", "LD": "#F0A800"}
HPFS = [24, 48, 72, 96, 120]

# ============================================================
# FIGURE 1
# ============================================================
fig1 = plt.figure(figsize=(17.5, 5.4))
gs = fig1.add_gridspec(1, 3, width_ratios=[1.0, 1.5, 1.0], wspace=0.85)

# ---------- Panel 1A: Funnel ----------
axA = fig1.add_subplot(gs[0, 0])
stages = [
    ("rows in family\n(F1 948, F2 739)",                                       1687, "#AAAAAA"),
    ("q ≤ 0.10\n(BH-FDR within family)",                                        146, "#4C72B0"),
    ("Tier A (50)\nq ≤ 0.10, snr ≥ 2,\nclean + replicated",                      50, "#2E7D32"),
    ("Tier B (63)\nq ≤ 0.10, snr ≥ 2,\nclean; fails replication / N",            63, "#B07000"),
    ("Tier C (33)\nq ≤ 0.10 but fails\n> 1 quality gate\n(12/33 floor-artefact)", 33, "#B00020"),
]
# format thousands with a comma for the n = labels
def _fmt(n): return f"{n:,}"
ys = np.arange(len(stages))[::-1]
max_n = max(s[1] for s in stages)
for i, (lab, n, c) in enumerate(stages):
    axA.barh(ys[i], n, color=c, edgecolor="black", linewidth=0.6, alpha=0.9)
    label = f"n = {_fmt(n)}"
    if n >= max_n * 0.4:
        axA.text(n - 25, ys[i], label, va="center", ha="right",
                 fontsize=9, fontweight="bold", color="white")
    else:
        axA.text(n + 25, ys[i], label, va="center", ha="left",
                 fontsize=9, fontweight="bold")
    axA.text(-50, ys[i], lab, va="center", ha="right", fontsize=8.5)
axA.set_xlim(-250, max_n * 1.08)
axA.set_ylim(-0.5, len(stages) - 0.5)
axA.set_xticks([]); axA.set_yticks([])
for spine in ("top", "right", "bottom", "left"):
    axA.spines[spine].set_visible(False)
axA.set_title("A — family funnel", fontsize=11, loc="left", pad=6)

# ---------- Panel 1B: Positive controls, grouped & labelled ----------
axB = fig1.add_subplot(gs[0, 1])
pos_specs = [
    ("KO\n(hard knockout)",          R / "stats_v7_positivecontrol/LD_72_vs_LD_72_KO_steroid.csv",  "#2E7D32"),
    ("KD\n(FVA-max cap —\ndid not fire)", R / "stats_v7_positivecontrol/LD_72_vs_LD_72_KD_steroid.csv",  "#4C72B0"),
    ("KD2\n(sampled-max cap —\nfired)", R / "stats_v7_positivecontrol/LD_72_vs_LD_72_KD2_steroid.csv", "#B07000"),
]
y_offset = 0
yticks, ylabels = [], []
group_bounds = []  # (y_start, y_end, label, color)
for name, fp, color in pos_specs:
    if not fp.exists(): continue
    d = pd.read_csv(fp)
    top = d.sort_values("W1_obs", ascending=False).head(5)
    y_start = y_offset
    for i, (_, r) in enumerate(top.iterrows()):
        y = y_offset + i
        axB.scatter(r.null_median, y, color="grey", s=36, alpha=0.6, marker="o")
        axB.scatter(r.W1_obs, y, color=color, s=70, edgecolor="black", linewidths=0.8, zorder=5)
        axB.plot([r.null_median, r.W1_obs], [y, y], color=color, alpha=0.5, lw=1.2)
        # Convert chapter's "Beta oxidation" to β-oxidation; wrap long labels.
        sub = r.subsystem.replace("Beta oxidation", "β-oxidation").replace("Beta-oxidation", "β-oxidation")
        if len(sub) > 32:
            mid = len(sub) // 2
            split = sub.rfind(" ", 0, mid + 8)
            if split > 0:
                sub = sub[:split] + "\n" + sub[split + 1:]
            else:
                sub = sub[:30] + "…"
        if r.get("p_emp", 1.0) <= 0.001:
            axB.text(r.W1_obs + 0.05, y, "***", fontsize=10, fontweight="bold",
                     va="center", color=color)
        yticks.append(y); ylabels.append(sub)
    group_bounds.append((y_start, y_offset + len(top) - 1, name, color))
    y_offset += len(top) + 1.5

axB.set_yticks(yticks); axB.set_yticklabels(ylabels, fontsize=8)
axB.set_xlabel("W1 (null median → observed)", fontsize=9)
axB.set_title("B — positive controls at LD_72 (steroid pathway)", fontsize=11, loc="left")
axB.grid(axis="x", alpha=0.3)
axB.invert_yaxis()

# Group labels on the right axis
ax_right = axB.twinx()
ax_right.set_ylim(axB.get_ylim())
ax_right.set_yticks([(gb[0] + gb[1]) / 2 for gb in group_bounds])
ax_right.set_yticklabels([gb[2] for gb in group_bounds], fontsize=8)
ax_right.tick_params(axis="y", length=0, pad=10)
for gb in group_bounds:
    ax_right.axhspan(gb[0] - 0.4, gb[1] + 0.4, xmin=0.0, xmax=0.015,
                     color=gb[3], alpha=0.7)

# Marker legend
handles = [plt.Line2D([0], [0], marker="o", color="grey", linestyle="",
                      markersize=7, label="null median"),
           plt.Line2D([0], [0], marker="o", color="black", linestyle="",
                      markersize=8, markeredgecolor="black", label="observed W1")]
axB.legend(handles=handles, loc="lower right", fontsize=8, frameon=False)

# ---------- Panel 1C: Negative-control noise floor ----------
axC = fig1.add_subplot(gs[0, 2])
nc = pd.read_csv(R / "v7_final/negcontrol_all_rows.csv")
axC.hist(nc.W1_obs, bins=35, color="grey", alpha=0.7,
         edgecolor="black", linewidth=0.4)
pct95 = nc.W1_obs.quantile(0.95)
axC.axvline(pct95, color="red", linestyle="--", linewidth=1.3,
            label=f"95th-pct noise floor = {pct95:.2f}")
axC.set_xlabel("W1 (split-chain null)", fontsize=9)
axC.set_ylabel("count", fontsize=9)
axC.set_title("C — split-chain negative control", fontsize=11, loc="left")
axC.legend(fontsize=7.5, loc="upper right")
axC.grid(axis="y", alpha=0.3)
# Headline annotation
n_sig_nc = int((nc.q_bh <= 0.10).sum())
axC.text(0.02, 0.72,
         f"n = {len(nc)} tests\n(15 cells × ~65 subsystems)\n"
         f"significant: {n_sig_nc}/{len(nc)} at q ≤ 0.10\n"
         f"empirical FDR ≈ {100*n_sig_nc/len(nc):.1f} %",
         transform=axC.transAxes, fontsize=8,
         bbox=dict(facecolor="white", edgecolor="grey", linewidth=0.5))

fig1.suptitle("Figure 1 — family funnel, positive controls, split-chain noise floor",
              y=1.03, fontsize=12)
fig1.text(0.5, -0.03,
          "Tier definitions (§2.12): A = clean + replicated; B = clean, fails replication/N; "
          "C = fails > 1 gate.  "
          "KO zeroes the target reactions so they drop out of the W1 comparison; "
          "only subsystems that still carry non-steroid flux (e.g. Protein assembly) remain testable, "
          "which is why the KO panel is sparse. "
          "KD2 fired at W1 = 2.29 (snr ≈ 38) once the knockdown was imposed on the sampled flux cap, "
          "not the FVA-max cap — see §3.5.",
          ha="center", fontsize=8, style="italic", wrap=True)
fig1.savefig(OUT / "fig1_funnel.pdf", bbox_inches="tight")
fig1.savefig(OUT / "fig1_funnel.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig1_funnel.pdf'} / .png")
plt.close(fig1)


# ============================================================
# FIGURE 5
# ============================================================
fig5 = plt.figure(figsize=(17, 5.2))
gs = fig5.add_gridspec(1, 3, width_ratios=[1.2, 1.2, 1.0], wspace=0.5)

# ---------- Panel 5A: Parent vs child dumbbell ----------
axA = fig5.add_subplot(gs[0, 0])
pc = pd.read_csv(R / "tier0_v7/tier0_parent_child_v7.csv")
pc = pc.sort_values(["hpf", "condition"]).reset_index(drop=True)
ys = np.arange(len(pc))
for i, (_, r) in enumerate(pc.iterrows()):
    color = COND_COLOR.get(r.condition, "#555555")
    axA.plot([r.child_biomax, r.parent_biomax], [i, i], color=color, lw=1.6, alpha=0.6)
    axA.scatter(r.child_biomax, i, color=color, s=70, zorder=5)
    axA.scatter(r.parent_biomax, i, color=color, s=70, marker="s", zorder=5,
                edgecolor="black", linewidths=0.8)
axA.set_yticks(ys)
axA.set_yticklabels([r.cell for _, r in pc.iterrows()], fontsize=8)
axA.set_xlabel("biomass_max (relative units)", fontsize=9)
axA.set_title("A — parent (base) vs child (extracted) biomass", fontsize=11, loc="left")
axA.grid(axis="x", alpha=0.3)
legend_items = [
    plt.Line2D([0], [0], marker="s", color="grey", linestyle="", markersize=9,
               markeredgecolor="black", label="parent"),
    plt.Line2D([0], [0], marker="o", color="grey", linestyle="", markersize=9, label="child"),
]
for cond, c in COND_COLOR.items():
    legend_items.append(plt.Line2D([0], [0], marker="o", color=c, linestyle="",
                                   markersize=9, label=cond))
axA.legend(handles=legend_items, loc="lower right", fontsize=8, frameon=False)

# ---------- Panel 5B: 2×2 attribution — horizontal stacked bar ----------
# Values are the Table 5 topology/envelope shares (ground truth from the
# chapter), paired with the signed total_gap from the C3 CSV.
# This avoids re-deriving the share formula from scratch — same source as
# the table, so the figure cannot drift from it.
axB = fig5.add_subplot(gs[0, 1])

# Shares now computed over |topology| + |envelope| + |interaction| (so bars
# are honest about how much of the raw gap the two additive terms describe).
# For the four small-gap pairs the interaction segment is 0; for the two
# BL_96 pairs the interaction is a majority of the gap.
TABLE5 = pd.DataFrame([
    # pair,              env_%, topo_%, inter_%, gap (raw AA−BB), interaction (signed, units)
    ("BL_96 – LD_96",     1.9,   19.5,   78.7,   -4.39,          -3.59),
    ("BL_96 – D_96",     17.8,   22.0,   60.3,   -4.04,          -2.43),
    ("BL_120 – LD_120", 92.9,    7.1,    0.0,   -0.84,           0.00),
    ("D_120 – LD_120",  95.1,    4.9,    0.0,   -0.73,           0.00),
    ("BL_72 – D_72",   100.0,    0.0,    0.0,   +1.57,           0.00),
    ("D_72 – LD_72",   100.0,    0.0,    0.0,   -1.35,           0.00),
], columns=["pair", "env", "topo", "inter", "gap", "interaction"])

# Three-segment stacked bar: envelope (green) + topology (red) + interaction (grey).
# Each segment is a share of |topo| + |env| + |interaction|, so the bar is
# honest about how much of the raw gap lives in the non-separable term.
ys = np.arange(len(TABLE5))
env_color   = "#2E7D32"
topo_color  = "#B00020"
inter_color = "#808080"

for i, (_, r) in enumerate(TABLE5.iterrows()):
    axB.barh(ys[i], r.env, color=env_color, edgecolor="black", linewidth=0.4, alpha=0.85)
    axB.barh(ys[i], r.topo, left=r.env, color=topo_color, edgecolor="black",
             linewidth=0.4, alpha=0.85)
    axB.barh(ys[i], r.inter, left=r.env + r.topo, color=inter_color,
             edgecolor="black", linewidth=0.4, alpha=0.85)
    # Mid-segment numeric labels, only when segment ≥ 10 % to stay legible.
    if r.env >= 10:
        axB.text(r.env / 2, ys[i], f"{r.env:.1f}", ha="center", va="center",
                 fontsize=7.5, color="white", fontweight="bold")
    if r.topo >= 10:
        axB.text(r.env + r.topo / 2, ys[i], f"{r.topo:.1f}", ha="center",
                 va="center", fontsize=7.5, color="white", fontweight="bold")
    if r.inter >= 10:
        axB.text(r.env + r.topo + r.inter / 2, ys[i], f"{r.inter:.1f}",
                 ha="center", va="center", fontsize=7.5, color="white",
                 fontweight="bold")
    # Signed raw gap + interaction magnitude on the right.
    if abs(r.interaction) > 0.1:
        gap_str = f"Δ = {r.gap:+.2f}\n(interaction {r.interaction:+.2f})"
    else:
        gap_str = f"Δ = {r.gap:+.2f}"
    axB.text(103, ys[i], gap_str, va="center", fontsize=7.5, color="black")

axB.set_yticks(ys)
axB.set_yticklabels(TABLE5.pair, fontsize=9)
axB.invert_yaxis()
axB.set_xlim(0, 130)
axB.set_xticks([0, 25, 50, 75, 100])
axB.set_xticklabels(["0", "25", "50", "75", "100"], fontsize=8)
axB.set_xlabel("attribution share (%)  [envelope → topology → interaction]",
               fontsize=9)
axB.set_title("B — Tier 0 2×2 attribution (shares of |env|+|topo|+|interaction|)",
              fontsize=10.5, loc="left")
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# Legend: three swatches, no counts (the verdict column lives in Table 5).
leg_handles = [
    plt.Rectangle((0, 0), 1, 1, color=env_color,  label="envelope"),
    plt.Rectangle((0, 0), 1, 1, color=topo_color, label="topology"),
    plt.Rectangle((0, 0), 1, 1, color=inter_color, label="interaction (non-separable)"),
]
axB.legend(handles=leg_handles, loc="lower right", fontsize=8, frameon=False,
           bbox_to_anchor=(1.0, -0.26), ncol=3)

# ---------- Panel 5C: BL_96 greedy add-back ----------
axC = fig5.add_subplot(gs[0, 2])
ab = pd.read_csv(R / "tier0_v7/bl96_greedy_addback.csv")
ab_plot = ab.head(15).copy()
colors_step = ["#B00020" if r.rxn == "MAR00466" else "#1F77B4"
               for _, r in ab_plot.iterrows()]
axC.bar(range(len(ab_plot)), ab_plot.gain_this_step, color=colors_step,
        edgecolor="black", linewidth=0.5)
mar_idx = ab_plot[ab_plot.rxn == "MAR00466"].index.tolist()
if mar_idx:
    i = mar_idx[0]
    # The step-1 "86 %" is of the greedy ceiling (≈3.72 units), NOT the raw
    # parent-vs-child gap (≈18.1 units) and NOT Table 5's attribution gaps.
    # Make both denominators visible on the panel.
    frac = ab.frac_of_full_recovery.iloc[0]
    greedy_ceiling = ab.gain_this_step.iloc[0] / frac   # ≈ 3.72
    axC.annotate(f"MAR00466\n(geneless Transport)\n"
                 f"{frac:.0%} of greedy\nceiling ({greedy_ceiling:.1f} u)",
                 xy=(i, ab_plot.loc[i, "gain_this_step"]),
                 xytext=(i + 1.5, ab_plot.loc[i, "gain_this_step"] * 1.02),
                 fontsize=8, ha="left", va="top",
                 arrowprops=dict(arrowstyle="->", color="#B00020", lw=0.9))
axC.set_xlabel("add-back step", fontsize=9)
axC.set_ylabel("biomass_max gain", fontsize=9)
axC.set_title("C — BL_96 add-back", fontsize=11, loc="left")
axC.grid(axis="y", alpha=0.3)

fig5.suptitle("Figure 5 — Tier 0: capability loss, corrected attribution, and the BL_96 lesson",
              y=1.04, fontsize=12)
fig5.text(0.5, -0.03,
          "A: every extracted cell has biomass_max below its parent (base); BL_96 the extreme. "
          "B: shares of |env|+|topo|+|interaction| (interaction = raw gap − topology − envelope, "
          "the non-separable part of a 2×2 swap). The four small-gap pairs are explained by envelope "
          "alone. For the two BL_96 pairs the interaction segment carries 60–79 % of the deficit — "
          "topology and envelope are complementary, not independent, so the gap is not separable "
          "into additive effects. Δ is raw AA−BB (Table 5 total gap column); the parenthetical "
          "is the signed interaction. "
          "C: step-1 reaction MAR00466 (geneless Transport) adds 3.20 units; this is 86 % of the "
          "greedy-add-back ceiling (3.72 u) and about 18 % of the raw parent-vs-child gap (18.1 u).",
          ha="center", fontsize=8, style="italic", wrap=True)
fig5.savefig(OUT / "fig5_tier0.pdf", bbox_inches="tight")
fig5.savefig(OUT / "fig5_tier0.png", dpi=200, bbox_inches="tight")
print(f"[write] {OUT/'fig5_tier0.pdf'} / .png")
plt.close(fig5)

# Print the Table 5 ground-truth values for the record
print("\n--- Fig 5B — Table 5 shares (ground truth, hard-coded) ---")
print(TABLE5.to_string(index=False))
