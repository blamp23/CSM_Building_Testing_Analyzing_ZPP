# s05_analysis_v2.R -- DESCRIPTIVE loopless-FVA summary across 15 cells.
#
# What this script IS:
#   * Descriptive summaries of per-cell capability classification
#   * Descriptive summaries of per-cell flexibility (rel_range)
#   * Enzyme-subset-representative aggregation (not raw reactions)
#   * Descriptive comparisons: differences in counts / medians between
#     conditions within an hpf, and across hpf within a condition
#   * Figures that visualise those descriptive comparisons
#
# What this script IS NOT:
#   * A statistical test between conditions. Every ll-FVA interval here
#     is a deterministic function of one model, so there is no replication
#     axis for hypothesis testing. To get p-values that mean something,
#     an ensemble over RNA-seq replicates x metabolomics replicates is
#     required; that ensemble does not yet exist and this script does
#     not pretend it does.
#   * Any inference about gene essentiality. The renamed column
#     `required_at_90pct_biomass` is a network-level constraint at a
#     specific operating point, not a gene-KO prediction.
#
# Reads:
#   loopless_fva/results/augmented/llfva_aug_{cell}.csv  (from s04)
#
# Writes:
#   loopless_fva/results/analysis_v2/*.csv
#   loopless_fva/results/analysis_v2/*.png

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

here    <- "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
res_dir <- file.path(here, "loopless_fva", "results")
aug_dir <- file.path(res_dir, "augmented")
out_dir <- file.path(res_dir, "analysis_v2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

CELLS <- expand.grid(condition = c("BL","D","LD"),
                     hpf       = c(24,48,72,96,120)) |>
         mutate(cell = paste0(condition, "_", hpf)) |>
         arrange(condition, hpf)

DROP_SUBSYS <- c("Transport", "Transport reactions", "Exchange/demand reactions",
                 "Artificial reactions", "Isolated", "")

# ============================================================================
# 1. Load augmented data
# ============================================================================
message("Loading augmented llfva CSVs...")

# per-cell biomass_max (from llfva_summary.csv, filled in the FVA run)
bio_max_tbl <- read_csv(file.path(res_dir, "llfva_summary.csv"),
                        show_col_types = FALSE) |>
  select(cell, biomass_max) |>
  # BL_120 has NA in summary; fill with the value logged in the yesterday run
  mutate(biomass_max = ifelse(cell == "BL_120" & is.na(biomass_max),
                              8.5619, biomass_max))

llfva <- CELLS |>
  rowwise() |>
  mutate(dat = list(read_csv(file.path(aug_dir, sprintf("llfva_aug_%s.csv", cell)),
                             show_col_types = FALSE))) |>
  unnest(dat) |>
  ungroup() |>
  left_join(bio_max_tbl, by = "cell") |>
  mutate(hpf = factor(hpf, levels = c(24,48,72,96,120)),
         condition = factor(condition, levels = c("BL","D","LD")),
         blocked = as.logical(blocked),
         required_at_90pct_biomass = as.logical(required_at_90pct_biomass),
         range_gt_1 = as.logical(range_gt_1),
         binding_lb = as.logical(binding_lb),
         binding_ub = as.logical(binding_ub))

message(sprintf("  %d rxn-cell rows loaded", nrow(llfva)))

# ============================================================================
# 2. Reduce to enzyme-subset REPRESENTATIVES
#    Rationale: reactions in the same subset share GPR + flux always; counting
#    them all inflates counts for parts of the model (isozymes / compartmental
#    duplicates) at the expense of others. One row per (cell, representative)
#    is the biology-faithful unit.
# ============================================================================
message("Collapsing to enzyme-subset representatives...")

rep_data <- llfva |>
  group_by(cell, condition, hpf, representative, subsystem) |>
  summarise(
    n_members     = n(),
    subset_size   = first(subset_size),
    has_gpr       = first(has_gpr),
    # Use median across subset members for min/max. For a well-formed subset,
    # all members share flux exactly, so median == any member's value. This
    # is robust to the v4_slim-subset / v6-model reaction-list drift (some
    # named representatives aren't in the v6 model, and to the retracted-
    # collapse case where members can disagree.
    min_flux      = median(min_flux),
    max_flux      = median(max_flux),
    range         = max_flux - min_flux,
    midpoint      = (min_flux + max_flux)/2,
    biomass_max   = first(biomass_max),
    # rel_range (bounded [0,1]) is degenerate for irreversible forward
    # reactions (always 2/3), so we keep it but rely on rel_range_bio
    # -- range as fraction of the cell's biomass_max -- as the primary
    # biology-comparable flexibility measure.
    rel_range     = ifelse(abs(midpoint)+range > 0, range/(abs(midpoint)+range), 0),
    rel_range_bio = range / biomass_max,
    blocked       = all(blocked),
    required      = any(required_at_90pct_biomass),
    binding_ub    = any(binding_ub),
    binding_lb    = any(binding_lb),
    .groups = "drop"
  )

# biology-facing (drop structural subsystems)
rep_bio <- rep_data |> filter(!subsystem %in% DROP_SUBSYS)

message(sprintf("  reps total: %d ; biology-facing: %d",
                nrow(rep_data), nrow(rep_bio)))

write_csv(rep_data, file.path(out_dir, "reps_all.csv"))
write_csv(rep_bio,  file.path(out_dir, "reps_biology.csv"))

# ============================================================================
# 3. Descriptive counts per cell (reps, biology-facing)
# ============================================================================
counts <- rep_bio |>
  group_by(condition, hpf, cell) |>
  summarise(n_reps       = n(),
            n_blocked    = sum(blocked),
            n_required   = sum(required),
            n_binding_ub = sum(binding_ub),
            n_binding_lb = sum(binding_lb),
            median_rel_range_bio = median(rel_range_bio[!blocked]),
            .groups = "drop") |>
  mutate(frac_blocked  = n_blocked / n_reps,
         frac_required = n_required / n_reps)
write_csv(counts, file.path(out_dir, "01_counts_per_cell.csv"))

# ---- Figure 1: fraction required and fraction blocked over hpf, by condition ----
p1a <- counts |>
  ggplot(aes(x = hpf, y = n_required, color = condition, group = condition)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_color_manual(values = c(BL="#1F77B4", D="#111111", LD="#F0A800"),
                     labels = c(BL="Blue", D="Dark", LD="Light-Dark")) +
  labs(x = "hpf", y = "n representatives required at 90% biomass",
       title = "Required-representative count over development",
       subtitle = "descriptive; no replication -> no p-values") +
  theme_bw(11) + theme(legend.position = "bottom")

p1b <- counts |>
  ggplot(aes(x = hpf, y = median_rel_range_bio, color = condition, group = condition)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  scale_color_manual(values = c(BL="#1F77B4", D="#111111", LD="#F0A800"),
                     labels = c(BL="Blue", D="Dark", LD="Light-Dark")) +
  labs(x = "hpf", y = "median (range / biomass_max), unblocked reps",
       title = "Network flexibility over development",
       subtitle = "flexibility per rep, normalised by biomass flux") +
  theme_bw(11) + theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig01_capability_trajectories.png"),
       p1a / p1b, width = 9, height = 7, dpi = 150)

# ============================================================================
# 4. Class change between conditions at each hpf (descriptive counts, no test)
# ============================================================================
# For each pair (A,B) x hpf: how many reps change class from A to B?
class_of <- function(blocked, required) {
  case_when(blocked  ~ "blocked",
            required ~ "required",
            TRUE     ~ "optional")
}

rep_class <- rep_bio |>
  mutate(class = class_of(blocked, required)) |>
  select(representative, subsystem, condition, hpf, class)

pairs_mat <- function(a_cond, b_cond) {
  d <- rep_class |>
    filter(condition %in% c(a_cond, b_cond)) |>
    pivot_wider(id_cols = c(representative, subsystem, hpf),
                names_from = condition, values_from = class) |>
    drop_na()
  d <- d |> rename(A = !!a_cond, B = !!b_cond)
  d |> mutate(change = paste0(A, "->", B),
              pair   = paste0(a_cond, " vs ", b_cond))
}

class_changes <- bind_rows(
  pairs_mat("D",  "LD"),
  pairs_mat("BL", "LD"),
  pairs_mat("BL", "D")
)
write_csv(class_changes |> count(pair, hpf, change) |>
            pivot_wider(names_from = change, values_from = n, values_fill = 0),
          file.path(out_dir, "02_class_transitions_by_pair_hpf.csv"))

# List of reps whose class *actually differs* between conditions -- these are
# the biology-interesting rows (small set, name-level, no p-values).
class_diffs <- class_changes |>
  filter(A != B) |>
  arrange(pair, hpf, subsystem, representative)
write_csv(class_diffs, file.path(out_dir, "02_class_diff_reps_by_pair_hpf.csv"))

# ---- Figure 2: heatmap of D-vs-LD class changes per subsystem x hpf ----
d_vs_ld <- class_changes |>
  filter(pair == "D vs LD", A != B) |>
  count(subsystem, hpf, change)

if (nrow(d_vs_ld) > 0) {
  # subsystems with any D-vs-LD change
  sub_lst <- d_vs_ld |> group_by(subsystem) |>
    summarise(total = sum(n)) |> arrange(-total) |> pull(subsystem)
  mat <- d_vs_ld |>
    filter(change %in% c("required->optional","optional->required",
                          "blocked->optional","optional->blocked",
                          "required->blocked","blocked->required")) |>
    mutate(direction = case_when(
      change %in% c("required->optional","required->blocked") ~ "D lost capability",
      change %in% c("optional->required","blocked->required") ~ "D gained capability",
      TRUE ~ "other")) |>
    group_by(subsystem, hpf, direction) |>
    summarise(n = sum(n), .groups = "drop")
  p2 <- ggplot(mat, aes(x = hpf, y = subsystem, fill = direction, size = n)) +
    geom_point(shape = 21, color = "grey40") +
    scale_size_area(max_size = 8) +
    scale_fill_manual(values = c("D lost capability"="#B40426",
                                  "D gained capability"="#3B4CC0",
                                  "other"="grey")) +
    labs(x = "hpf", y = NULL,
         title = "D vs LD class transitions (biology-facing subsystems)",
         subtitle = "size = n representatives changing class; descriptive only") +
    theme_bw(9)
  ggsave(file.path(out_dir, "fig02_D_vs_LD_class_transitions.png"),
         p2, width = 9, height = min(14, 2 + 0.25*length(sub_lst)), dpi = 150)
}

# ============================================================================
# 5. Subsystem-level required fractions and flexibility, descriptive heatmap
# ============================================================================
sub_summary <- rep_bio |>
  group_by(subsystem, condition, hpf, cell) |>
  summarise(n_rep         = n(),
            n_required    = sum(required),
            frac_required = n_required / n_rep,
            median_rr_bio = median(rel_range_bio[!blocked]),
            .groups = "drop")
write_csv(sub_summary, file.path(out_dir, "03_subsystem_summary.csv"))

# ---- Figure 3: heatmap of frac_required, subsystem x cell ----
mat_req <- sub_summary |>
  select(subsystem, cell, frac_required) |>
  pivot_wider(names_from = cell, values_from = frac_required, values_fill = 0) |>
  as.data.frame()
rownames(mat_req) <- mat_req$subsystem; mat_req$subsystem <- NULL
mat_req <- as.matrix(mat_req[, CELLS$cell])
# keep subsystems with size >= 3 reps in any cell (else noisy)
sizes <- rep_bio |> group_by(subsystem, cell) |> summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = cell, values_from = n, values_fill = 0) |>
  as.data.frame()
rownames(sizes) <- sizes$subsystem; sizes$subsystem <- NULL
keep <- apply(as.matrix(sizes[, CELLS$cell]), 1, max) >= 3
mat_req <- mat_req[keep[rownames(mat_req)], , drop = FALSE]

col_ha <- HeatmapAnnotation(
  condition = CELLS$condition, hpf = as.character(CELLS$hpf),
  col = list(condition = c(BL="#1F77B4", D="#111111", LD="#F0A800"),
             hpf       = c("24"="#f7fbff","48"="#c6dbef","72"="#6baed6",
                           "96"="#2171b5","120"="#08306b")))

png(file.path(out_dir, "fig03_frac_required_heatmap.png"),
    width = 1400, height = 1600, res = 130)
Heatmap(mat_req,
        name = "frac_required",
        col = colorRamp2(c(0, 0.5, 1), c("#f7fbff","#6baed6","#08306b")),
        cluster_columns = FALSE, cluster_rows = TRUE,
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 9),
        show_row_dend = FALSE,
        column_title = "Fraction of representatives required at 90% biomass")
dev.off()

# ============================================================================
# 6. Key pathway trajectories (retinoid, GAG, xenobiotic)
# ============================================================================
PATHWAYS <- c(
  "Vitamin A metabolism", "Retinol metabolism",
  "Xenobiotics metabolism", "ROS detoxification",
  "Chondroitin / heparan sulfate biosynthesis",
  "Chondroitin sulfate degradation",
  "Heparan sulfate degradation",
  "Keratan sulfate biosynthesis",
  "Keratan sulfate degradation",
  "N-glycan metabolism", "O-glycan metabolism")

pathway_summary <- rep_data |>
  filter(subsystem %in% PATHWAYS) |>
  group_by(subsystem, condition, hpf) |>
  summarise(n_rep         = n(),
            n_required    = sum(required),
            n_blocked     = sum(blocked),
            median_rr_bio = median(rel_range_bio[!blocked]),
            .groups = "drop")
write_csv(pathway_summary, file.path(out_dir, "04_key_pathway_summary.csv"))

# ---- Figure 4: pathway-specific required + rel_range trajectories ----
p4a <- pathway_summary |>
  ggplot(aes(x = hpf, y = n_required, color = condition, group = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~ subsystem, scales = "free_y") +
  scale_color_manual(values = c(BL="#1F77B4", D="#111111", LD="#F0A800"),
                     labels = c(BL="Blue", D="Dark", LD="Light-Dark")) +
  labs(x = "hpf", y = "n representatives required",
       title = "Required-representative counts in key pathways") +
  theme_bw(10) + theme(legend.position = "bottom")
p4b <- pathway_summary |>
  ggplot(aes(x = hpf, y = median_rr_bio, color = condition, group = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~ subsystem, scales = "free_y") +
  scale_color_manual(values = c(BL="#1F77B4", D="#111111", LD="#F0A800"),
                     labels = c(BL="Blue", D="Dark", LD="Light-Dark")) +
  labs(x = "hpf", y = "median rel_range",
       title = "Flexibility (rel_range) in key pathways") +
  theme_bw(10) + theme(legend.position = "bottom")
ggsave(file.path(out_dir, "fig04_key_pathway_required.png"),
       p4a, width = 12, height = 8, dpi = 150)
ggsave(file.path(out_dir, "fig04_key_pathway_relrange.png"),
       p4b, width = 12, height = 8, dpi = 150)

# ============================================================================
# 7. Binding-constraint audit -- how many reps hit input lb/ub per cell
#    Useful for the methods: how much of the reported range is set by the box
#    vs by the network?
# ============================================================================
binding_audit <- rep_bio |>
  group_by(condition, hpf, cell) |>
  summarise(n_reps            = n(),
            n_binding_ub      = sum(binding_ub),
            n_binding_lb      = sum(binding_lb),
            n_binding_any     = sum(binding_lb | binding_ub),
            frac_binding_any  = mean(binding_lb | binding_ub),
            .groups = "drop")
write_csv(binding_audit, file.path(out_dir, "05_binding_constraint_audit.csv"))

# ============================================================================
# 8. INDEX
# ============================================================================
outputs <- list.files(out_dir, pattern = "\\.csv$|\\.png$", full.names = FALSE)
writeLines(sort(outputs), file.path(out_dir, "INDEX.txt"))

message("\nDONE. Descriptive outputs in ", out_dir)
message("Files: ", length(outputs))
message("\nRemember: no p-values from this script are legitimate. ",
        "To turn descriptive counts into significance claims, ",
        "an ensemble of models (RNA reps x metab reps) is required.")
