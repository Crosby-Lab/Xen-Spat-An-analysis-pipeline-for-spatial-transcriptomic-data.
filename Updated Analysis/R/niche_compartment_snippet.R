# ---------------------------------------------------------------------------
# Drop-in replacement for the block that assigned niches to the tumor or
# stromal compartment. Put this where the old largest-gap code sat, inside
# the niche analysis script, after BuildNicheAssay has been run and the
# metadata table (one row per cell, with the niche and cell type columns)
# exists.
#
# The assignment is unchanged. What is new is that the gap profile is written
# out, which is what Section 2.7 of the manuscript says the pipeline reports.
# ---------------------------------------------------------------------------

source("R/dice_gap_profile.R")        # adjust the path if the scripts move

# `meta` is the per-cell table: sample, niche, cell_type.
# `TUMOR_LABELS` is the cell type label or labels counted as tumor cells in
# this dataset. The labels differ between datasets, so set them per run:
#   discovery cohort        c("Tumor Cells")
#   validation S1_TOP/S4_TOP c("Tumor_Cells", "Proliferative_Tumor_Cells")
#   validation S3_TOP        c("Epithelial_&_Tumor_Cells")
#   validation S2_TOP        character(0)   # normal breast, no tumor label

prof <- dice_gap_profile(meta,
                         sample_col  = "sample",
                         niche_col   = "niche",
                         type_col    = "cell_type",
                         tumor_types = TUMOR_LABELS)

# the compartment assignment the rest of the pipeline consumes
meta$compartment <- prof$niches$compartment[
  match(paste(meta$sample, meta$niche), paste(prof$niches$sample, prof$niches$niche))]

# the gap profile, written next to the other per-sample outputs
write.csv(prof$niches,  file.path(outdir, "dice_gap_profile_niches.csv"),  row.names = FALSE)
write.csv(prof$samples, file.path(outdir, "dice_gap_profile_samples.csv"), row.names = FALSE)

# print the two diagnostics so a weak split is visible in the run log
for (i in seq_len(nrow(prof$samples))) with(prof$samples[i, ],
  message(sprintf("%-8s  gap ratio %-6s  tumor niches %-3s  min Dice above split %-6s  %s",
                  sample,
                  ifelse(is.na(gap_ratio), "NA", sprintf("%.2f", gap_ratio)),
                  ifelse(is.na(n_tumor_niches), "NA", n_tumor_niches),
                  ifelse(is.na(min_dice_above_split), "NA",
                         sprintf("%.2f", min_dice_above_split)),
                  ifelse(is.na(flag) | flag == "", "", paste("FLAG:", flag)))))

# A gap ratio near 1, or a low min Dice above the split, means the binary
# assignment should not be taken at face value. Section 2.7 and the Discussion
# describe what to do in that case: use the Dice coefficient as a continuous
# score, or repeat the downstream comparison using
# prof$niches$compartment_at_second_gap as a sensitivity check.
