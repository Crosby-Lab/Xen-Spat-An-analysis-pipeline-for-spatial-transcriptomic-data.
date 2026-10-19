# Dice gap profile for the tumor/stromal compartment boundary

`dice_gap_profile.R` replaces the step that chose the largest gap in the Dice
coefficients and returned only the resulting compartment assignment. It returns
the same assignment plus the quantities needed to judge whether that assignment
is well supported.

## Why

The compartment boundary is placed at the largest gap between consecutive sorted
Dice coefficients, so no threshold has to be chosen for a new dataset. That rule
assumes the coefficients separate into two groups. Two things can go wrong, and
they fail independently:

| Diagnostic | Failure signalled |
|---|---|
| `gap_ratio` near 1 | The coefficients are distributed smoothly rather than bimodally. The split is set by small differences between adjacent niches and is not stable. |
| `min_dice_above_split` near 0 | The split sits among uniformly low coefficients, which happens in a section with little tumor. Most niches get called tumor although none overlaps it substantially. |

Both are reported for every sample, together with the partition the
second-largest gap would have produced, so a sensitivity check needs no rerun.

## Use

```r
source("dice_gap_profile.R")

# cells: one row per cell after quality filtering, with sample, niche and
# cell type columns. The niche column is the BuildNicheAssay output.
res <- dice_gap_profile(cells,
                        sample_col  = "sample",
                        niche_col   = "niche",
                        type_col    = "cell_type",
                        tumor_types = c("Tumor Cells"))

res$samples   # one row per sample: gap_ratio, split_rank, flags
res$niches    # one row per niche: dice, rank, gap_below, compartment

write.csv(res$niches,  "dice_gap_profile_niches.csv",  row.names = FALSE)
write.csv(res$samples, "dice_gap_profile_samples.csv", row.names = FALSE)
```

For a Seurat object carrying a niche assay:

```r
res <- dice_gap_profile_seurat(obj, sample_id = "097",
                               niche_col = "niches", type_col = "cell_type")
```

If a sample has more than one tumor cell label, pass them all, for example
`tumor_types = c("Tumor Cells", "Tumor Cells_MKI67")`. Labels not listed are
treated as non-tumor.

## Output columns

`res$niches`, sorted by sample then by descending Dice:

| Column | Meaning |
|---|---|
| `n_cells`, `n_tumor_cells` | cells in the niche, and how many are tumor cells |
| `n_tumor_cells_in_sample` | tumor cells in the whole sample, the other set in the Dice coefficient |
| `dice` | 2 × `n_tumor_cells` / (`n_tumor_cells_in_sample` + `n_cells`) |
| `rank` | 1 is the highest Dice coefficient in that sample |
| `gap_below` | this niche's Dice minus the next rank's |
| `compartment` | tumor or stromal, at the largest gap |
| `compartment_at_second_gap` | the same at the second-largest gap, for a sensitivity check |

`res$samples`, one row per sample:

| Column | Meaning |
|---|---|
| `largest_gap`, `second_largest_gap` | the two largest gaps |
| `gap_ratio` | largest / second largest. `Inf` when the second is zero |
| `split_rank`, `n_tumor_niches` | number of niches above the split |
| `n_tumor_niches_at_second_gap` | number above the second-largest gap |
| `min_dice_above_split` | lowest Dice coefficient assigned to the tumor compartment |
| `max_dice_below_split` | highest Dice coefficient assigned to the stromal compartment |
| `flag` | empty when nothing was flagged, otherwise what was flagged and why |

## Edge cases handled

- **No tumor cells in the section.** Dice is undefined; no compartment is
  assigned and the flag says so. This is the normal breast section S2_TOP.
- **Fewer than three niches.** Only one gap exists, so no ratio is defined. The
  top niche is assigned to the tumor compartment and the flag says the ratio is
  unavailable.
- **Ties in the Dice coefficients.** Tied values sit at adjacent ranks with a gap
  of zero between them, so a split is never placed inside a tie.

## Thresholds

`min_gap_ratio` (default 2) and `min_dice` (default 0.2) only set when a flag is
raised; they never change the assignment. They are deliberately loose. In the
data reported in the manuscript the gap ratio was 3.9 to 28.5 in the six
discovery biopsies and 1.2 to 1.7 in the validation sections, and the lowest Dice
coefficient above the split was 0.45.

## Checking it

`test_dice_gap_profile.R` runs the gap arithmetic on the published Table S8
coefficients and confirms it reproduces the number of niches assigned to the
tumor compartment in all nine sections, then runs the full function on a small
synthetic cell table including a section with no tumor cells.

```
Rscript test_dice_gap_profile.R
```
