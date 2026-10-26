# s04_augment_llfva.R -- add lb/ub, rel_range, and clearer column names
# to the per-cell llfva CSVs. Does NOT require re-running MATLAB.
#
# Reads:
#   loopless_fva/results/llfva_{cell}.csv
#   loopless_fva/results/pfba_{cell}.csv          (source of lb/ub)
#   reviewer_packet_v4_slim/subsets/{cell}_subset_membership.csv
#
# Writes:
#   loopless_fva/results/augmented/llfva_aug_{cell}.csv
#
# Column changes vs original llfva CSVs:
#   +  lb, ub                     : model bounds after sampling override
#   +  binding_lb, binding_ub     : min/max sit on the input box (diagnostic)
#   +  rel_range                  : range / (|midpoint| + range)  in [0,1]
#   +  representative, subset_id  : enzyme-subset assignment
#   +  subset_size, has_gpr       : subset metadata
#   ~  essential  -> required_at_90pct_biomass   (renamed to avoid confusion
#                                                  with gene essentiality)
#   ~  flexible   -> range_gt_1     (kept; note this is the *absolute* metric,
#                                    use rel_range for biology-comparable)
#
# Note on lb/ub source: pfba CSVs store the bounds actually in force at
# pFBA solve time, which is the same sampling-style regime as the FVA.
# So joining lb/ub from pfba is equivalent to reading them from the FVA
# model without re-running MATLAB.

suppressPackageStartupMessages({
  library(tidyverse)
})

here <- "/Users/lamp_b/Library/CloudStorage/OneDrive-TexasA&MUniversity/Hala, David's files - Benji_COBRA/Tanguay_Data/Discrete_Models"
res_dir <- file.path(here, "loopless_fva", "results")
subset_dir <- file.path(here, "reviewer_packet_v4_slim", "subsets")
aug_dir <- file.path(res_dir, "augmented")
dir.create(aug_dir, showWarnings = FALSE, recursive = TRUE)

CELLS <- expand.grid(condition = c("BL","D","LD"),
                     hpf       = c(24,48,72,96,120)) |>
         mutate(cell = paste0(condition, "_", hpf))

TOL_ZERO  <- 1e-6
TOL_BOUND <- 1e-3   # loose enough to catch numerical near-touches

drift_log <- list()

for (i in seq_len(nrow(CELLS))) {
  cell_i <- CELLS$cell[i]
  message(sprintf("[%d/%d] %s ...", i, nrow(CELLS), cell_i))

  fva <- read_csv(file.path(res_dir, sprintf("llfva_%s.csv", cell_i)),
                  show_col_types = FALSE) |>
    mutate(blocked   = as.logical(blocked),
           essential = as.logical(essential),
           flexible  = as.logical(flexible))

  pfba <- read_csv(file.path(res_dir, sprintf("pfba_%s.csv", cell_i)),
                   show_col_types = FALSE) |>
    select(rxn, lb, ub, v_pfba)

  subs <- read_csv(file.path(subset_dir, sprintf("%s_subset_membership.csv", cell_i)),
                   show_col_types = FALSE) |>
    select(rxn, representative, subset_id, subset_size, has_gpr)

  n_fva <- nrow(fva); n_pfba <- nrow(pfba); n_subs <- nrow(subs)

  aug <- fva |>
    left_join(pfba, by = "rxn") |>
    left_join(subs, by = "rxn") |>
    rename(required_at_90pct_biomass = essential,
           range_gt_1                = flexible) |>
    mutate(
      # rel_range: bounded [0,1]. 0 = pinned reaction, 1 = maximally flexible.
      midpoint  = (min_flux + max_flux) / 2,
      rel_range = ifelse((abs(midpoint) + range) > 0,
                          range / (abs(midpoint) + range),
                          0),
      # binding diagnostics
      binding_lb = !is.na(lb) & (min_flux <= lb + TOL_BOUND),
      binding_ub = !is.na(ub) & (max_flux >= ub - TOL_BOUND),
      # subset representative fallback (some rxns not in subset file)
      representative = coalesce(representative, rxn),
      subset_size    = coalesce(subset_size, 1L),
      has_gpr        = coalesce(has_gpr, NA)
    ) |>
    select(rxn, name, subsystem, representative, subset_id, subset_size, has_gpr,
           lb, ub, min_flux, max_flux, midpoint, range, rel_range,
           binding_lb, binding_ub,
           blocked, required_at_90pct_biomass, range_gt_1)

  n_missing_lbub  <- sum(is.na(aug$lb) | is.na(aug$ub))
  n_missing_subs  <- sum(is.na(aug$subset_id))

  drift_log[[cell_i]] <- tibble(cell = cell_i,
                                n_fva = n_fva, n_pfba = n_pfba, n_subs = n_subs,
                                n_missing_lbub = n_missing_lbub,
                                n_missing_subset = n_missing_subs)

  write_csv(aug, file.path(aug_dir, sprintf("llfva_aug_%s.csv", cell_i)))
}

drift_tbl <- bind_rows(drift_log)
write_csv(drift_tbl, file.path(aug_dir, "join_drift_report.csv"))

message("\nJoin drift report:")
print(drift_tbl, n = 15)
message(sprintf("\nDONE. Augmented CSVs in: %s", aug_dir))
