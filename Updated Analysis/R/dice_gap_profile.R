# ---------------------------------------------------------------------------
# Xen-Spat: Dice coefficient between tumor cells and each niche, with the
# gap profile used to place the tumor/stromal compartment boundary.
#
# This replaces the step that selected the largest gap and returned only the
# resulting assignment. It returns the same assignment plus the quantities a
# user needs to judge whether that assignment is well supported:
#   - every niche's Dice coefficient, in rank order
#   - the gap below each rank
#   - the largest and second-largest gap and their ratio (the "gap ratio")
#   - the partition that the second-largest gap would have produced
#   - the lowest Dice coefficient above the split
#
# Dice is computed from the cell-level labels as
#     2 * (tumor cells in the niche) / (tumor cells in the sample + cells in the niche)
#
# Two diagnostics are reported because they fail independently:
#   gap_ratio            near 1  -> the coefficients are distributed smoothly,
#                                  the split is set by small differences and the
#                                  partition is not stable.
#   min_dice_above_split near 0  -> the split sits among uniformly low
#                                  coefficients, which happens in sections with
#                                  little tumor; most niches are then called
#                                  tumor although none overlaps it substantially.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Core arithmetic, exposed so that it can be checked against a published table
# of Dice coefficients without needing the cell-level data.
#   dice: numeric vector, one value per niche, any order
# Returns the sorted values and the two candidate splits.
# ---------------------------------------------------------------------------
gap_profile_from_dice <- function(dice) {
  d <- sort(dice[!is.na(dice)], decreasing = TRUE)
  k <- length(d)
  if (k < 3)
    return(list(dice_sorted = d, gaps = if (k == 2) -diff(d) else numeric(0),
                split_rank = if (k >= 1) 1L else NA_integer_,
                split_rank_second = NA_integer_,
                largest_gap = if (k == 2) -diff(d) else NA_real_,
                second_largest_gap = NA_real_, gap_ratio = NA_real_,
                min_dice_above_split = if (k >= 1) d[1] else NA_real_))
  gaps <- -diff(d)
  ord  <- order(gaps, decreasing = TRUE)
  g1 <- gaps[ord[1]]; g2 <- gaps[ord[2]]
  list(dice_sorted = d, gaps = gaps,
       split_rank = ord[1], split_rank_second = ord[2],
       largest_gap = g1, second_largest_gap = g2,
       gap_ratio = if (g2 > 0) g1 / g2 else Inf,
       min_dice_above_split = d[ord[1]])
}

#' @param cells      data.frame, one row per cell, after quality filtering
#' @param sample_col column holding the sample identifier
#' @param niche_col  column holding the niche label (BuildNicheAssay output)
#' @param type_col   column holding the uniform cell type label
#' @param tumor_types cell type labels counted as tumor cells
#' @param min_gap_ratio  gap ratios below this are flagged as weakly separated
#' @param min_dice       values of min_dice_above_split below this are flagged
#' @return list(niches = per-niche data.frame, samples = per-sample data.frame)
dice_gap_profile <- function(cells,
                             sample_col = "sample",
                             niche_col  = "niche",
                             type_col   = "cell_type",
                             tumor_types = c("Tumor Cells"),
                             min_gap_ratio = 2,
                             min_dice = 0.2) {

  stopifnot(all(c(sample_col, niche_col, type_col) %in% names(cells)))
  cells <- cells[!is.na(cells[[niche_col]]), , drop = FALSE]
  if (!any(cells[[type_col]] %in% tumor_types))
    warning("No cell matches tumor_types: ", paste(tumor_types, collapse = ", "))

  niche_rows  <- list()
  sample_rows <- list()

  for (s in sort(unique(as.character(cells[[sample_col]])))) {
    d <- cells[as.character(cells[[sample_col]]) == s, , drop = FALSE]
    n_tumor_total <- sum(d[[type_col]] %in% tumor_types)

    nl    <- as.character(d[[niche_col]])
    n_cls <- tapply(nl, nl, length)
    n_tum <- tapply(d[[type_col]] %in% tumor_types, nl, sum)
    niches <- names(n_cls)

    dice <- if (n_tumor_total == 0) rep(NA_real_, length(niches)) else
      as.numeric(2 * n_tum[niches] / (n_tumor_total + n_cls[niches]))

    tab <- data.frame(sample = s, niche = niches,
                      n_cells = as.integer(n_cls[niches]),
                      n_tumor_cells = as.integer(n_tum[niches]),
                      n_tumor_cells_in_sample = n_tumor_total,
                      dice = dice, stringsAsFactors = FALSE)

    # ---- no tumor cells: no partition exists -------------------------------
    if (n_tumor_total == 0) {
      tab$rank <- NA_integer_; tab$gap_below <- NA_real_
      tab$compartment <- NA_character_
      tab$compartment_at_second_gap <- NA_character_
      niche_rows[[s]] <- tab
      sample_rows[[s]] <- data.frame(
        sample = s, n_niches = nrow(tab), n_tumor_cells_in_sample = 0L,
        largest_gap = NA_real_, second_largest_gap = NA_real_,
        gap_ratio = NA_real_, split_rank = NA_integer_,
        n_tumor_niches = NA_integer_, n_tumor_niches_at_second_gap = NA_integer_,
        min_dice_above_split = NA_real_, max_dice_below_split = NA_real_,
        flag = "no tumor cells; no compartment assigned",
        stringsAsFactors = FALSE)
      next
    }

    tab <- tab[order(-tab$dice), , drop = FALSE]
    tab$rank <- seq_len(nrow(tab))
    k <- nrow(tab)
    tab$gap_below <- c(-diff(tab$dice), NA_real_)   # gap between rank i and i+1

    # ---- fewer than three niches: a gap ratio is not defined ---------------
    if (k < 3) {
      split1 <- 1L
      tab$compartment <- ifelse(tab$rank <= split1, "tumor", "stromal")
      tab$compartment_at_second_gap <- NA_character_
      niche_rows[[s]] <- tab
      sample_rows[[s]] <- data.frame(
        sample = s, n_niches = k, n_tumor_cells_in_sample = n_tumor_total,
        largest_gap = if (k == 2) tab$gap_below[1] else NA_real_,
        second_largest_gap = NA_real_, gap_ratio = NA_real_,
        split_rank = split1, n_tumor_niches = split1,
        n_tumor_niches_at_second_gap = NA_integer_,
        min_dice_above_split = tab$dice[split1],
        max_dice_below_split = if (k == 2) tab$dice[2] else NA_real_,
        flag = "fewer than three niches; gap ratio not defined",
        stringsAsFactors = FALSE)
      next
    }

    gp     <- gap_profile_from_dice(tab$dice)
    split1 <- gp$split_rank; split2 <- gp$split_rank_second
    g1 <- gp$largest_gap;    g2 <- gp$second_largest_gap
    gap_ratio <- gp$gap_ratio

    tab$compartment <- ifelse(tab$rank <= split1, "tumor", "stromal")
    tab$compartment_at_second_gap <- ifelse(tab$rank <= split2, "tumor", "stromal")

    flags <- character(0)
    if (is.finite(gap_ratio) && gap_ratio < min_gap_ratio)
      flags <- c(flags, sprintf("gap ratio %.2f below %.2f: coefficients distributed smoothly, split not stable",
                                gap_ratio, min_gap_ratio))
    if (tab$dice[split1] < min_dice)
      flags <- c(flags, sprintf("lowest Dice above the split is %.2f: split sits among low coefficients",
                                tab$dice[split1]))
    if (split1 >= k - 1)
      flags <- c(flags, "split at the lowest rank: all but one niche assigned to the tumor compartment")

    niche_rows[[s]] <- tab
    sample_rows[[s]] <- data.frame(
      sample = s, n_niches = k, n_tumor_cells_in_sample = n_tumor_total,
      largest_gap = g1, second_largest_gap = g2, gap_ratio = gap_ratio,
      split_rank = split1, n_tumor_niches = split1,
      n_tumor_niches_at_second_gap = split2,
      min_dice_above_split = tab$dice[split1],
      max_dice_below_split = tab$dice[split1 + 1],
      flag = if (length(flags)) paste(flags, collapse = "; ") else "",
      stringsAsFactors = FALSE)
  }

  list(niches  = do.call(rbind, c(niche_rows,  list(make.row.names = FALSE))),
       samples = do.call(rbind, c(sample_rows, list(make.row.names = FALSE))))
}


#' Convenience wrapper for a Seurat object carrying a niche assay.
#' @param obj      Seurat object for one sample
#' @param sample_id identifier written into the output
#' @param niche_col  meta.data column with the niche labels (BuildNicheAssay
#'                   writes "niches" by default)
#' @param type_col   meta.data column with the uniform cell type labels
dice_gap_profile_seurat <- function(obj, sample_id,
                                    niche_col = "niches",
                                    type_col  = "cell_type",
                                    tumor_types = c("Tumor Cells"), ...) {
  md <- obj@meta.data
  stopifnot(all(c(niche_col, type_col) %in% names(md)))
  cells <- data.frame(sample = sample_id,
                      niche = as.character(md[[niche_col]]),
                      cell_type = as.character(md[[type_col]]),
                      stringsAsFactors = FALSE)
  dice_gap_profile(cells, tumor_types = tumor_types, ...)
}
