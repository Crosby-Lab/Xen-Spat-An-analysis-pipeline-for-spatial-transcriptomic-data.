# Check gap_profile_from_dice against the published Table S8 values, and check
# dice_gap_profile end to end on a synthetic cell table.
source("dice_gap_profile.R")

tableS8 <- list(
  "046"    = c(0.02, 0.12, 0.03, 0.05, 0.09, 0.83),
  "047"    = c(0.04, 0.14, 0.22, 0.60, 0.15, 0.08),
  "049"    = c(0.05, 0.70, 0.19, 0.03, 0.06, 0.06),
  "089"    = c(0.02, 0.74, 0.08, 0.04, 0.02, 0.01, 0.01, 0.04, 0.01, 0.08),
  "097"    = c(0.07, 0.06, 0.04, 0.10, 0.02, 0.08, 0.04, 0.01, 0.05, 0.67),
  "100"    = c(0.13, 0.64, 0.13, 0.08, 0.05, 0.03, 0.09, 0.01, 0.07, 0.02),
  "S1_TOP" = c(0.01, 0.00, 0.01, 0.02, 0.45, 0.00, 0.00, 0.75, 0.01, 0.09, 0.00),
  "S3_TOP" = c(0.01, 0.00, 0.77, 0.00, 0.09, 0.00, 0.34),
  "S4_TOP" = c(0.25, 0.06, 0.00, 0.39, 0.03, 0.02, 0.00, 0.00, 0.63)
)
# published number of niches in the tumor compartment
expected_n_tumor <- c("046"=1,"047"=1,"049"=1,"089"=1,"097"=1,"100"=1,
                      "S1_TOP"=2,"S3_TOP"=1,"S4_TOP"=1)

cat(sprintf("%-8s %6s %6s %7s %9s %9s\n",
            "sample","gap1","gap2","ratio","n_tumor","min_dice"))
ok <- TRUE
for (s in names(tableS8)) {
  g <- gap_profile_from_dice(tableS8[[s]])
  cat(sprintf("%-8s %6.2f %6.2f %7.2f %9d %9.2f\n",
              s, g$largest_gap, g$second_largest_gap, g$gap_ratio,
              g$split_rank, g$min_dice_above_split))
  if (g$split_rank != expected_n_tumor[[s]]) {
    cat("   MISMATCH: expected", expected_n_tumor[[s]], "tumor niches\n"); ok <- FALSE
  }
}
cat(if (ok) "\nAll assignments match Table S8.\n" else "\nMISMATCH FOUND.\n")

# ---- end to end on a synthetic cell table ---------------------------------
set.seed(1)
mk <- function(sample, niche, n_cells, n_tumor)
  data.frame(sample = sample, niche = niche,
             cell_type = c(rep("Tumor Cells", n_tumor),
                           rep("Fibroblasts", n_cells - n_tumor)),
             stringsAsFactors = FALSE)
cells <- rbind(
  mk("X1","N1", 1000, 800), mk("X1","N2", 1200,  40),
  mk("X1","N3",  900,  20), mk("X1","N4", 1100,  30),
  # a section with no tumor cells at all
  mk("X2","N1", 500, 0), mk("X2","N2", 600, 0), mk("X2","N3", 700, 0)
)
res <- dice_gap_profile(cells, tumor_types = "Tumor Cells")
print(res$samples[, c("sample","n_niches","gap_ratio","n_tumor_niches",
                      "min_dice_above_split","flag")])
print(res$niches[, c("sample","niche","n_cells","n_tumor_cells","dice","rank",
                     "gap_below","compartment","compartment_at_second_gap")])
