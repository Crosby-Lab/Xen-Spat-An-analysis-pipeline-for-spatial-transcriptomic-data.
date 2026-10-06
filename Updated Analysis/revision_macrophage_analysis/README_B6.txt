Xen-Spat revision, task B6: macrophage near vs margin reanalysis (run 2026-10-04)

Input: Integrated_Data_CelltypeAnnotated_9225.h5ad (Jesu, Box). 348,246 cells x 329 genes.
Sample map (batch -> sample): BC1=046, BC2=047, BC3=049, BC4=089, BC5=097, BC6=100.
Counts: integer counts reconstructed exactly from adata.raw (log1p of library-size-normalized data);
reconstructed per-cell totals equal obs.total_counts for every cell.

Cell labels: integrated (BBKNN) annotation in obs.leiden, because per-sample uniform labels are not
stored in this file. Tumor cells = all "Tumor Cells_*" clusters.

Distance: for each sample, Euclidean distance (um) from each non-tumor cell centroid to the nearest
tumor cell centroid. Near < 25 um, margin 25 to 75 um, far > 75 um.

1. Epithelial transcripts vs distance (B6_epithelial_vs_distance.csv, B6_Fig_epithelial_vs_distance):
   genes EPCAM, KRT7, KRT8, KRT19, CDH1, GRHL2, DSP, AR, ESR1, ERBB2, VTCN1 (KRT18 not in panel);
   mean counts per cell in 5 um bins (bins with < 10 cells omitted), macrophages and fibroblasts.

2. Tumor-enriched genes (B6_tumor_enriched_genes.csv): log2[(mean in tumor cells + 0.01) /
   (mean in macrophages > 50 um from tumor + 0.01)] > 1 in at least 4 of 6 samples, using counts
   per 1,000 transcripts. 87 genes flagged. (> 50 um rather than > 75 um because some samples have
   very few far macrophages, e.g. 10 in 047.)

3. Pseudobulk DE (B6_pseudobulk_DE_all_genes.csv, B6_pseudobulk_DE_excluded.csv): counts summed per
   patient for near and margin macrophages (12 profiles; cell numbers in B6_pseudobulk_metadata.csv),
   PyDESeq2 0.5.4, design ~ patient + zone, contrast near vs margin, Benjamini-Hochberg.
   Run with all 329 genes and with the 87 tumor-enriched genes removed.
