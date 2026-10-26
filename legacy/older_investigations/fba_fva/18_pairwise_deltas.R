# s06_pairwise_deltas.R -- DESCRIPTIVE pairwise comparison of ll-FVA
# capacity between conditions at matched hpf.
#
# Two products:
#   1. Jaccard trajectories per condition-pair over hpf
#      (three set-versions: required, blocked, and combined-class agreement)
#   2. Subsystem-level delta heatmaps: median |Delta_range| per
#      subsystem x hpf, for LD-vs-D and LD-vs-BL.
#
# Reads: loopless_fva/results/analysis_v2/reps_all.csv  (from s05)
# Writes: loopless_fva/results/analysis_v2/pairwise/*.csv, *.png
#
# Descriptive only. No p-values. Deltas combine topology + envelope effects.

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

here    <- "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
res_dir <- file.path(here, "loopless_fva", "results")
in_dir  <- file.path(res_dir, "analysis_v2")
out_dir <- file.path(in_dir, "pairwise")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DROP_SUBSYS <- c("Transport", "Transport reactions", "Exchange/demand reactions",
                 "Artificial reactions", "Isolated", "")

message("Loading reps_all.csv ...")
reps <- read_csv(file.path(in_dir, "reps_all.csv"), show_col_types = FALSE) |>
  mutate(hpf = factor(hpf, levels = c(24,48,72,96,120)),
         condition = factor(condition, levels = c("BL","D","LD")),
         class = case_when(
           blocked  ~ "blocked",
           required ~ "required",
           TRUE     ~ "optional"))

reps_bio <- reps |>
  filter(!subsystem %in% DROP_SUBSYS) |>
  # collapse the small number of reps that appear under multiple subsystems
  # (e.g. MAR05282 in both Pool reactions and Protein degradation) to a single
  # canonical row per (cell, representative) with the largest-subsystem tag.
  group_by(cell, representative) |>
  slice_head(n = 1) |>
  ungroup()
message(sprintf("  reps total: %d ; biology-facing: %d",
                nrow(reps), nrow(reps_bio)))

# ============================================================================
# 1. Jaccard trajectories
# ============================================================================
message("\nJaccard trajectories per condition-pair...")

# pair a cell's reps by representative name; intersect
pair_at_hpf <- function(h, cA, cB) {
  a <- reps_bio |> filter(condition == cA, hpf == h) |>
       select(representative, subsystem, class_A = class,
              range_A = range, max_A = max_flux, min_A = min_flux)
  b <- reps_bio |> filter(condition == cB, hpf == h) |>
       select(representative, class_B = class,
              range_B = range, max_B = max_flux, min_B = min_flux)
  inner_join(a, b, by = "representative")
}

pairs_def <- tribble(
  ~pair,     ~cA,  ~cB,
  "LD_vs_D", "LD", "D",
  "LD_vs_BL","LD", "BL",
  "D_vs_BL", "D",  "BL"
)

HPF_LEVELS <- c("24","48","72","96","120")
jaccard_traj <- pairs_def |>
  rowwise() |>
  mutate(
    dat = list(map(HPF_LEVELS, function(h) {
      pd <- pair_at_hpf(h, cA, cB)
      req_A <- pd$representative[pd$class_A == "required"]
      req_B <- pd$representative[pd$class_B == "required"]
      blk_A <- pd$representative[pd$class_A == "blocked"]
      blk_B <- pd$representative[pd$class_B == "blocked"]
      tibble(
        hpf = h,
        n_shared_reps = nrow(pd),
        jaccard_required = length(intersect(req_A, req_B)) /
                            max(1, length(union(req_A, req_B))),
        jaccard_blocked  = length(intersect(blk_A, blk_B)) /
                            max(1, length(union(blk_A, blk_B))),
        class_agreement  = mean(pd$class_A == pd$class_B))
    }) |> bind_rows())
  ) |>
  unnest(dat) |>
  ungroup() |>
  mutate(hpf = factor(hpf, levels = HPF_LEVELS))

write_csv(jaccard_traj, file.path(out_dir, "01_jaccard_trajectories.csv"))

# ---- Figure: Jaccard trajectories ----
pj <- jaccard_traj |>
  pivot_longer(c(jaccard_required, jaccard_blocked, class_agreement),
               names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                         jaccard_required = "Jaccard(required)",
                         jaccard_blocked  = "Jaccard(blocked)",
                         class_agreement  = "class agreement (rep-level)")) |>
  ggplot(aes(x = hpf, y = value, color = pair, group = pair)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = c(LD_vs_D  = "#111111",
                                LD_vs_BL = "#F0A800",
                                D_vs_BL  = "#1F77B4"),
                     labels = c(LD_vs_D  = "Dark vs Light-Dark",
                                LD_vs_BL = "Blue vs Light-Dark",
                                D_vs_BL  = "Blue vs Dark")) +
  labs(x = "hpf", y = "similarity",
       title = "Pairwise similarity of ll-FVA capacity classification",
       subtitle = "1 = identical between conditions; descriptive, single-model") +
  theme_bw(11) + theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig01_jaccard_trajectories.png"),
       pj, width = 10, height = 4.5, dpi = 150)

# ============================================================================
# 2. Per-representative delta table (for the heatmap + top-N tables)
# ============================================================================
message("Building per-rep deltas...")

delta_rep <- pairs_def |>
  rowwise() |>
  mutate(dat = list(
    map(HPF_LEVELS, function(h) {
      pd <- pair_at_hpf(h, cA, cB)
      pd |>
        mutate(hpf = h,
               d_range = range_A - range_B,
               d_max   = max_A - max_B,
               d_min   = min_A - min_B,
               abs_d_range = abs(d_range),
               abs_d_max   = abs(d_max),
               # Relative shift: |Δrange| / mean(range across both cells).
               # Dimensionless -> comparable across subsystems with very
               # different absolute flux magnitudes. Uses a small epsilon
               # to avoid division-by-zero when both cells have range=0.
               mean_range  = 0.5 * (range_A + range_B),
               rel_d_range = abs_d_range / (mean_range + 1e-6),
               class_change = ifelse(class_A == class_B, NA_character_,
                                     paste0(class_B, "->", class_A)))
    }) |> bind_rows()
  )) |>
  unnest(dat) |>
  ungroup() |>
  mutate(hpf = factor(hpf, levels = HPF_LEVELS))

write_csv(delta_rep, file.path(out_dir, "02_rep_deltas_by_pair_hpf.csv.gz"))

# ============================================================================
# 3. Subsystem-level medians of |Delta_range| and |Delta_max|
# ============================================================================
sub_delta <- delta_rep |>
  group_by(pair, hpf, subsystem) |>
  summarise(n_rep_shared  = n(),
            n_class_diff  = sum(!is.na(class_change)),
            # UNSIGNED (magnitude): how big is the change, regardless of direction
            median_abs_dr = median(abs_d_range),
            median_abs_dm = median(abs_d_max),
            max_abs_dr    = max(abs_d_range),
            median_rel_dr = median(rel_d_range),
            # SIGNED (direction): positive => A (=cA) > B; negative => B > A.
            # For LD_vs_D: positive = LD has more range = "D lost capability."
            # For LD_vs_BL: positive = LD has more range = "BL lost capability."
            # For D_vs_BL: positive = D has more range = "BL lost capability."
            median_signed_dr     = median(d_range),
            median_signed_rel_dr = median(d_range / (0.5*(range_A+range_B) + 1e-6)),
            .groups = "drop")

write_csv(sub_delta, file.path(out_dir, "03_subsystem_deltas.csv"))

# ---- Figure: subsystem x hpf heatmap of chosen delta metric, per pair ----
# metric_col: which column of sub_delta to plot.
# signed: TRUE for diverging palette (direction: A > B is red, B > A is blue);
#         FALSE for sequential palette (magnitude only).
make_heatmap <- function(pair_name, out_png, metric_col, title_metric,
                         legend_name, signed = FALSE) {
  dat <- sub_delta |>
    filter(pair == pair_name, n_rep_shared >= 3)  # drop tiny subsystems
  if (nrow(dat) == 0) return(invisible(NULL))
  mat <- dat |>
    select(subsystem, hpf, value = all_of(metric_col)) |>
    pivot_wider(names_from = hpf, values_from = value, values_fill = NA) |>
    as.data.frame()
  rownames(mat) <- mat$subsystem; mat$subsystem <- NULL
  mat <- as.matrix(mat[, c("24","48","72","96","120")])
  # drop subsystems whose reps are missing in >1 hpf (else hclust NAs)
  keep_row <- rowSums(is.na(mat)) <= 1
  mat <- mat[keep_row, , drop = FALSE]
  # cap extreme outliers for visualisation
  if (signed) {
    cap <- quantile(abs(mat), 0.98, na.rm = TRUE)
    mat_c <- pmin(pmax(mat, -cap), cap)
    pal <- colorRamp2(c(-cap, 0, cap), c("#3B4CC0","white","#B40426"))
  } else {
    cap <- quantile(mat, 0.98, na.rm = TRUE)
    mat_c <- pmin(mat, cap)
    pal <- colorRamp2(c(0, cap/2, cap), c("#f7fbff","#fdbb84","#7f0000"))
  }
  # fill remaining single-cell NAs with row median for clustering
  for (i in seq_len(nrow(mat_c))) {
    rm <- median(mat_c[i, ], na.rm = TRUE)
    mat_c[i, is.na(mat_c[i, ])] <- rm
  }

  # size column: total shared reps at that hpf across all subsystems (context)
  col_size <- sub_delta |> filter(pair == pair_name) |>
    group_by(hpf) |> summarise(n = sum(n_rep_shared), .groups = "drop") |>
    arrange(hpf) |> pull(n)

  col_ha <- HeatmapAnnotation(
    hpf = as.character(c(24,48,72,96,120)),
    col = list(hpf = c("24"="#f7fbff","48"="#c6dbef","72"="#6baed6",
                       "96"="#2171b5","120"="#08306b")),
    n_shared = anno_barplot(col_size, height = unit(1, "cm"))
  )

  png(out_png, width = 1400, height = 1600, res = 130)
  print(Heatmap(mat_c,
      name = legend_name,
      col  = pal,
      cluster_columns = FALSE, cluster_rows = TRUE,
      top_annotation = col_ha,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 9),
      show_row_dend = FALSE,
      column_title = sprintf("Subsystem-level %s (%s), capped at 98%%tile",
                             title_metric, pair_name),
      cell_fun = function(j, i, x, y, width, height, fill) {
        v <- dat$n_class_diff[match(paste(rownames(mat)[i], colnames(mat)[j]),
                                    paste(dat$subsystem, dat$hpf))]
        if (!is.na(v) && v > 0)
          grid.text(v, x, y, gp = gpar(fontsize = 6, col = "black"))
      }))
  dev.off()
}

# ABSOLUTE-shift heatmaps (biased toward high-flux subsystems)
make_heatmap("LD_vs_D",  file.path(out_dir, "fig02a_subsystem_delta_abs_LD_vs_D.png"),
             "median_abs_dr", "|Δrange| (mmol/gDW/hr)", "median |Δrange|")
make_heatmap("LD_vs_BL", file.path(out_dir, "fig02b_subsystem_delta_abs_LD_vs_BL.png"),
             "median_abs_dr", "|Δrange| (mmol/gDW/hr)", "median |Δrange|")
make_heatmap("D_vs_BL",  file.path(out_dir, "fig02c_subsystem_delta_abs_D_vs_BL.png"),
             "median_abs_dr", "|Δrange| (mmol/gDW/hr)", "median |Δrange|")

# RELATIVE-shift heatmaps (dimensionless, comparable across subsystems)
make_heatmap("LD_vs_D",  file.path(out_dir, "fig02d_subsystem_delta_rel_LD_vs_D.png"),
             "median_rel_dr", "|Δrange| / mean(range) (dimensionless)",
             "rel |Δrange|")
make_heatmap("LD_vs_BL", file.path(out_dir, "fig02e_subsystem_delta_rel_LD_vs_BL.png"),
             "median_rel_dr", "|Δrange| / mean(range) (dimensionless)",
             "rel |Δrange|")
make_heatmap("D_vs_BL",  file.path(out_dir, "fig02f_subsystem_delta_rel_D_vs_BL.png"),
             "median_rel_dr", "|Δrange| / mean(range) (dimensionless)",
             "rel |Δrange|")

# SIGNED RELATIVE-shift heatmaps (diverging palette; direction matters)
# For LD_vs_D: red = LD has more range (D lost capability),
#              blue = D has more range (D gained capability).
# For LD_vs_BL: red = LD > BL; blue = BL > LD.
# For D_vs_BL:  red = D > BL;  blue = BL > D.
make_heatmap("LD_vs_D",  file.path(out_dir, "fig02g_subsystem_delta_signed_LD_vs_D.png"),
             "median_signed_rel_dr",
             "signed relative Δrange -- red: LD>D (D lost), blue: D>LD (D gained)",
             "signed rel Δrange", signed = TRUE)
make_heatmap("LD_vs_BL", file.path(out_dir, "fig02h_subsystem_delta_signed_LD_vs_BL.png"),
             "median_signed_rel_dr",
             "signed relative Δrange -- red: LD>BL, blue: BL>LD",
             "signed rel Δrange", signed = TRUE)
make_heatmap("D_vs_BL",  file.path(out_dir, "fig02i_subsystem_delta_signed_D_vs_BL.png"),
             "median_signed_rel_dr",
             "signed relative Δrange -- red: D>BL, blue: BL>D",
             "signed rel Δrange", signed = TRUE)

# ============================================================================
# 6. Publication-quality faceted heatmap (all three pairs on one page)
#    - filter to subsystems where max|signed_rel_dr| >= threshold in ANY pair
#      (union) so all three panels share the same subsystem axis and can be
#      compared column-by-column
#    - landscape orientation: hpf rows, subsystems columns
#    - subsystem labels angled 45 degrees, truncated at 42 chars
#    - signed diverging palette (blue = pair's second condition has more range,
#      red = pair's first condition has more range)
#    - cell text overlay = n representatives whose CLASS (blocked / required at
#      90% biomass / optional) also changed between the two conditions
# ============================================================================

THRESHOLD <- 0.2
MIN_N_REP <- 5
TRUNC_LEN <- 42

`%||%` <- function(a, b) if (!is.null(a)) a else b

# labels for each pair
PAIR_TITLE <- c(
  LD_vs_D  = "Dark  vs  Light-Dark",
  LD_vs_BL = "Blue  vs  Light-Dark",
  D_vs_BL  = "Blue  vs  Dark"
)

# build a signed_rel_dr matrix per pair: rows = subsystem, cols = hpf
build_pair_mat <- function(pair_name) {
  d <- sub_delta |> filter(pair == pair_name, n_rep_shared >= MIN_N_REP)
  m <- d |>
    select(subsystem, hpf, value = median_signed_rel_dr) |>
    pivot_wider(names_from = hpf, values_from = value, values_fill = NA) |>
    as.data.frame()
  rownames(m) <- m$subsystem; m$subsystem <- NULL
  as.matrix(m[, c("24","48","72","96","120")])
}
build_overlay_mat <- function(pair_name) {
  d <- sub_delta |> filter(pair == pair_name, n_rep_shared >= MIN_N_REP)
  m <- d |>
    select(subsystem, hpf, value = n_class_diff) |>
    pivot_wider(names_from = hpf, values_from = value, values_fill = 0) |>
    as.data.frame()
  rownames(m) <- m$subsystem; m$subsystem <- NULL
  as.matrix(m[, c("24","48","72","96","120")])
}

mat_LD_D  <- build_pair_mat("LD_vs_D")
mat_LD_BL <- build_pair_mat("LD_vs_BL")
mat_D_BL  <- build_pair_mat("D_vs_BL")

# union of subsystems passing threshold in ANY pair
pass_row <- function(mat) {
  s <- apply(abs(mat), 1, max, na.rm = TRUE)
  names(s)[!is.na(s) & s >= THRESHOLD]
}
sub_union <- unique(c(pass_row(mat_LD_D), pass_row(mat_LD_BL), pass_row(mat_D_BL)))
message(sprintf("Union of subsystems passing |max| >= %.2f in any pair: %d",
                THRESHOLD, length(sub_union)))

# align every matrix to the union
align_mat <- function(mat, subs) {
  out <- matrix(NA_real_, nrow = length(subs), ncol = 5,
                dimnames = list(subs, c("24","48","72","96","120")))
  common <- intersect(subs, rownames(mat))
  out[common, ] <- mat[common, ]
  out
}
mLD_D  <- align_mat(mat_LD_D,  sub_union)
mLD_BL <- align_mat(mat_LD_BL, sub_union)
mD_BL  <- align_mat(mat_D_BL,  sub_union)

# cluster columns (subsystems) once, using the sum of |signed| across all pairs,
# so the same order is used in all three panels.
cluster_mat <- pmax(abs(mLD_D), abs(mLD_BL), abs(mD_BL), na.rm = TRUE)
cluster_mat[is.na(cluster_mat)] <- 0
sub_order <- hclust(dist(cluster_mat))$order
sub_ordered <- rownames(cluster_mat)[sub_order]

reorder_mat <- function(m) m[sub_ordered, , drop = FALSE]
mLD_D  <- reorder_mat(mLD_D)
mLD_BL <- reorder_mat(mLD_BL)
mD_BL  <- reorder_mat(mD_BL)
oLD_D  <- align_mat(build_overlay_mat("LD_vs_D"),  sub_union)[sub_ordered, ]
oLD_BL <- align_mat(build_overlay_mat("LD_vs_BL"), sub_union)[sub_ordered, ]
oD_BL  <- align_mat(build_overlay_mat("D_vs_BL"),  sub_union)[sub_ordered, ]

# shared color scale across all 3 panels
cap <- max(quantile(abs(c(mLD_D, mLD_BL, mD_BL)), 0.98, na.rm = TRUE), THRESHOLD)
cap_v <- function(m) { m[is.na(m)] <- 0; pmin(pmax(m, -cap), cap) }
tLD_D  <- t(cap_v(mLD_D))
tLD_BL <- t(cap_v(mLD_BL))
tD_BL  <- t(cap_v(mD_BL))
tolo_LD_D  <- t(oLD_D)
tolo_LD_BL <- t(oLD_BL)
tolo_D_BL  <- t(oD_BL)

# condition color anchors
COL_LD <- "#F0A800"   # Light-Dark = gold/yellow
COL_D  <- "#111111"   # Dark = black
COL_BL <- "#1F77B4"   # Blue

# Per-panel diverging palettes: positive => first condition (cA) dominates
# (has more range); negative => second (cB) dominates. The color in each
# cell IS the condition with more capacity at that subsystem x hpf.
pal_LD_D  <- colorRamp2(c(-cap, 0, cap), c(COL_D,  "white", COL_LD))
pal_LD_BL <- colorRamp2(c(-cap, 0, cap), c(COL_BL, "white", COL_LD))
pal_D_BL  <- colorRamp2(c(-cap, 0, cap), c(COL_BL, "white", COL_D))

# truncate long subsystem labels
trunc_name <- function(s) ifelse(nchar(s) > TRUNC_LEN,
                                 paste0(substr(s, 1, TRUNC_LEN-1), "…"), s)
colnames(tLD_D)  <- trunc_name(colnames(tLD_D))
colnames(tLD_BL) <- trunc_name(colnames(tLD_BL))
colnames(tD_BL)  <- trunc_name(colnames(tD_BL))

# figure size scales with number of subsystem columns
n_cols <- ncol(tLD_D)
w_in <- max(11, 0.32 * n_cols + 4)
h_in <- 8.5   # room for 3 panels vertically + caption

# cell_fun factories
cf_factory <- function(overlay_mat) {
  function(j, i, x, y, width, height, fill) {
    v <- overlay_mat[i, j]
    if (!is.na(v) && v > 0)
      grid.text(v, x, y, gp = gpar(fontsize = 7, col = "black"))
  }
}

# build the three Heatmap objects with row_title = the pair-explanation string
ht_LD_D <- Heatmap(tLD_D,
    name = "Dark ↔ Light-Dark",
    col = pal_LD_D,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_dend = FALSE, show_row_dend = FALSE,
    row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    row_title = PAIR_TITLE["LD_vs_D"], row_title_gp = gpar(fontsize = 10),
    row_title_side = "right", row_title_rot = 0,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                labels_gp = gpar(fontsize = 8),
                                at = c(-cap, 0, cap),
                                labels = c("Dark", "0", "Light-Dark")),
    show_column_names = FALSE,
    cell_fun = cf_factory(tolo_LD_D))

ht_LD_BL <- Heatmap(tLD_BL,
    name = "Blue ↔ Light-Dark",
    col = pal_LD_BL,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_dend = FALSE, show_row_dend = FALSE,
    row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    row_title = PAIR_TITLE["LD_vs_BL"], row_title_gp = gpar(fontsize = 10),
    row_title_side = "right", row_title_rot = 0,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                labels_gp = gpar(fontsize = 8),
                                at = c(-cap, 0, cap),
                                labels = c("Blue", "0", "Light-Dark")),
    show_column_names = FALSE,
    cell_fun = cf_factory(tolo_LD_BL))

ht_D_BL <- Heatmap(tD_BL,
    name = "Blue ↔ Dark",
    col = pal_D_BL,
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_column_dend = FALSE, show_row_dend = FALSE,
    row_names_side = "left", row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    row_title = PAIR_TITLE["D_vs_BL"], row_title_gp = gpar(fontsize = 10),
    row_title_side = "right", row_title_rot = 0,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                                labels_gp = gpar(fontsize = 8),
                                at = c(-cap, 0, cap),
                                labels = c("Blue", "0", "Dark")),
    cell_fun = cf_factory(tolo_D_BL))

png(file.path(out_dir, "fig03_pub_faceted_signed_deltas.png"),
    width = w_in * 130, height = h_in * 130, res = 130)
draw(ht_LD_D %v% ht_LD_BL %v% ht_D_BL,
     column_title = sprintf(
       "Signed relative Δrange between conditions  ·  |max| ≥ %.2f, n ≥ %d reps  ·  cell number = n reps changing capacity class (blocked / required / optional)",
       THRESHOLD, MIN_N_REP),
     column_title_gp = gpar(fontsize = 10, fontface = "italic"),
     padding = unit(c(6, 12, 8, 12), "mm"))
dev.off()

# ============================================================================
# 4. Top-N "widest capacity" per cell (biology-facing) -- what the paper
#    would call "which reactions have the most room to vary?"
# ============================================================================
top_capacity <- reps_bio |>
  filter(!blocked) |>
  group_by(cell) |>
  slice_max(range, n = 25, with_ties = FALSE) |>
  arrange(cell, desc(range)) |>
  select(cell, condition, hpf, representative, subsystem, min_flux, max_flux, range)
write_csv(top_capacity, file.path(out_dir, "04_top25_capacity_per_cell.csv"))

# Which representatives appear in top-25 in LD but NOT in D (or vice versa)
# at each hpf -- "highest-capacity reps that differ between conditions"
top_diff <- lapply(HPF_LEVELS, function(h) {
  ld <- top_capacity |> filter(hpf == h, condition == "LD") |> pull(representative)
  d  <- top_capacity |> filter(hpf == h, condition == "D")  |> pull(representative)
  bl <- top_capacity |> filter(hpf == h, condition == "BL") |> pull(representative)
  tibble(hpf = h,
         LD_only_vs_D  = paste(setdiff(ld, d),  collapse = ";"),
         D_only_vs_LD  = paste(setdiff(d,  ld), collapse = ";"),
         LD_only_vs_BL = paste(setdiff(ld, bl), collapse = ";"),
         BL_only_vs_LD = paste(setdiff(bl, ld), collapse = ";"))
}) |> bind_rows()
write_csv(top_diff, file.path(out_dir, "04_top25_capacity_diff_names.csv"))

# ============================================================================
# 7. Compensation analysis: do Blue and Dark diverge from Light-Dark in the
#    SAME direction across subsystems, or in OPPOSITE directions?
#
#    Per hpf, correlate signed_rel_dr of LD_vs_D vs LD_vs_BL across
#    subsystems. Report:
#      * Pearson r  (magnitude of linear agreement)
#      * sign concordance (fraction of subsystems where both deltas share
#        the same sign) -- a "same-side" statistic
#      * median |LD_vs_D| and median |LD_vs_BL| for context
#    Also make a per-hpf scatter plot.
# ============================================================================
message("\nCompensation analysis: are Blue and Dark on the same side of LD?")

comp_dat <- sub_delta |>
  filter(pair %in% c("LD_vs_D","LD_vs_BL"), n_rep_shared >= MIN_N_REP) |>
  select(pair, subsystem, hpf, median_signed_rel_dr) |>
  pivot_wider(names_from = pair, values_from = median_signed_rel_dr) |>
  drop_na()

comp_summary <- comp_dat |>
  group_by(hpf) |>
  summarise(
    n_subsystems      = n(),
    pearson_r         = cor(LD_vs_D, LD_vs_BL, method = "pearson"),
    spearman_r        = cor(LD_vs_D, LD_vs_BL, method = "spearman"),
    sign_concordance  = mean(sign(LD_vs_D) == sign(LD_vs_BL)),
    n_both_positive   = sum(LD_vs_D > 0 & LD_vs_BL > 0),
    n_both_negative   = sum(LD_vs_D < 0 & LD_vs_BL < 0),
    n_opposite        = sum(sign(LD_vs_D) * sign(LD_vs_BL) < 0),
    median_abs_LD_D   = median(abs(LD_vs_D)),
    median_abs_LD_BL  = median(abs(LD_vs_BL)),
    .groups = "drop")
write_csv(comp_summary, file.path(out_dir, "06_compensation_summary.csv"))
message("\nCompensation summary (LD_vs_D  vs  LD_vs_BL, per hpf):")
print(comp_summary)

# Per-hpf scatter -- each dot = one subsystem
p_comp <- comp_dat |>
  mutate(quadrant = case_when(
    LD_vs_D > 0 & LD_vs_BL > 0 ~ "both lost vs LD",
    LD_vs_D < 0 & LD_vs_BL < 0 ~ "both gained vs LD",
    TRUE                       ~ "opposite / one flat")) |>
  ggplot(aes(x = LD_vs_D, y = LD_vs_BL)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(intercept = 0, slope = 1, color = "grey50",
              linetype = "dashed") +
  geom_point(aes(color = quadrant), alpha = 0.75, size = 2) +
  facet_wrap(~ hpf, nrow = 1) +
  scale_color_manual(values = c("both lost vs LD"     = "#B40426",
                                "both gained vs LD"   = "#3B4CC0",
                                "opposite / one flat" = "grey60")) +
  labs(x = "signed rel Δrange  (Dark - Light-Dark)",
       y = "signed rel Δrange  (Blue - Light-Dark)",
       color = NULL,
       title = "Do Blue and Dark shift the same way against Light-Dark?",
       subtitle = "each dot = subsystem;  dashed line = perfect agreement (Dark and Blue identical)") +
  theme_bw(10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig04_compensation_scatter.png"),
       p_comp, width = 14, height = 4.5, dpi = 150)

# ============================================================================
# 8. Index
# ============================================================================
outputs <- list.files(out_dir, pattern = "\\.csv$|\\.csv\\.gz$|\\.png$", full.names = FALSE)
writeLines(sort(outputs), file.path(out_dir, "INDEX.txt"))
message("\nDONE. Outputs in ", out_dir)
message("Files: ", length(outputs))
