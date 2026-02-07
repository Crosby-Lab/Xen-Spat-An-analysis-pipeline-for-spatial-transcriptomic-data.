###############################################################################
### Generating Celltype annotataion from Xenium Raw Data
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date :
### Environment : Python 3.9.0
### Packages  : numpy 1.23.4, panda 2.3.3, matplotlib 3.7.2, seaborn 0.12.2, scanpy 1.9.3, squidpy 1.6.1
#######################################################################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


######## Use  custom folder

adata = sc.read_10x_h5(
    filename="tutorial_data/xenium_data/output_XETG00054_0006530_046_20230623_154919/cell_feature_matrix.h5"
)


df = pd.read_csv(
    "tutorial_data/xenium_data/output_XETG00054_0006530_046_20230623_154919/cells.csv"
)

print(df)


df.set_index(adata.obs_names, inplace=True)
adata.obs = df.copy()

adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

print(adata.obsm["spatial"])

print(" Size of DataFrame:", adata.var_names )


sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)

print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")


fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)


axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)

plt.show()


sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=5)


# #IF YOU ARE DOING MOUSE YOU MIGHT NEED TO CHANGE MT- to Mt. Always double check you actually labeld MT

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'

print(adata.var)

print(adata.var[adata.var.mt == True])


sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, .02)
print(f'{lower_lim} to {upper_lim}')


adata = adata[(adata.obs.n_genes_by_counts < upper_lim) & (adata.obs.n_genes_by_counts > lower_lim)]

sc.pp.normalize_total(adata, inplace=True)

sc.pp.log1p(adata) #change to log counts


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5) #these are default values

sc.pl.highest_expr_genes(adata)

sc.pp.regress_out(adata, ['total_counts']) #Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed

sc.pp.scale(adata, max_value=10) #scale each gene to unit variance

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

sc.pl.pca_loadings(adata,include_lowest = True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)

sc.tl.umap(adata)
sc.pl.umap(adata)

sc.tl.leiden(adata, resolution = 0.25)
sc.pl.umap(adata, color=['leiden'])

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, fontsize =14)


############ Labelling Leiden Clusters  ###########################################################################

#046
new_cluster_names = [
    'Myeloid Cells', 'Breast Epithelium_1','Breast Epithelium_2',
    'Endothelial Cells','Breast Epithelium_3','Lymphoid Cells_1',
    'Lymphoid Cells_2', 'Fibroblasts','Breast Epithelium_4']


adata.rename_categories('leiden', new_cluster_names)

sc.pl.umap(adata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')


sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],

    wspace=0.4,
    )

# ############################ SUBSet ANALYSIS  046 ###########################################################

#046

print(adata.obs.leiden.unique().tolist())

M = ['Myeloid Cells','Lymphoid Cells_1','Lymphoid Cells_2']

bdata = adata[adata.obs.leiden.isin(M)]

sc.pl.umap(bdata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')

sc.pp.neighbors(bdata)
sc.tl.umap(bdata)
sc.tl.leiden(bdata, resolution = 0.3, key_added="leiden")


sc.pl.umap(bdata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')


sc.tl.rank_genes_groups(bdata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(bdata, n_genes=10, sharey=False,fontsize =14)


new_cluster_names1 = [
    'Macrophages','CD8 T Cells', 'Macrophages(MMP12)','Macrophages(CXCL10',
    'Tregs','Macrophages(IL2RA)','Breast Epithelial Cells','Epithelial Cells(IFNg)',
    'Epithelial Cells(MZB1)','Macrophages(CXCL5)']


bdata.rename_categories('leiden', new_cluster_names1)


sc.pl.umap(bdata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')

sc.pl.rank_genes_groups_dotplot(bdata, groupby="leiden", standard_scale="var", n_genes=5)

sq.pl.spatial_scatter(
    bdata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],

    wspace=0.4,
    )



M = [ 'Breast Epithelium_1','Breast Epithelium_2',
      'Endothelial Cells','Breast Epithelium_3', 'Fibroblasts','Breast Epithelium_4']

cdata = adata[adata.obs.leiden.isin(M)]

alldata = bdata.concatenate([cdata])  ## Combining all clusters

sc.pl.umap(alldata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')

sc.tl.rank_genes_groups(alldata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(alldata, n_genes=10, sharey=False,fontsize =14)

sc.tl.dendrogram(alldata, 'leiden')
sc.pl.dendrogram(alldata, 'leiden')

sc.pl.rank_genes_groups_dotplot(alldata, n_genes=5, groupby = 'leiden')


######################## Uniform Clusteriing     ######################################################################
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Breast Epithelial Cells', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Breast Epithelium_1', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Breast Epithelium_2', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Breast Epithelium_3', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Breast Epithelium_4', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Epithelial Cells(IFNg)', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Epithelial Cells(MZB1)', 'Tumor cells')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Macrophages(MMP12)', 'Macrophages')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Macrophages(CXCL10', 'Macrophages')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Macrophages(CXCL5)', 'Macrophages')
alldata.obs['leiden'] = alldata.obs['leiden'].replace('Macrophages(IL2RA)', 'Macrophages')
############################################################################################################


sc.pl.umap(alldata, color='leiden',  legend_fontsize= 'large',legend_loc='on data')

sc.tl.rank_genes_groups(alldata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(alldata, n_genes=10, sharey=False,fontsize =14)

sc.tl.dendrogram(alldata, 'leiden')
sc.pl.dendrogram(alldata, 'leiden')


sc.pl.rank_genes_groups_heatmap(alldata, n_genes=5, groupby = 'leiden', show_gene_labels = True)
sc.pl.rank_genes_groups_dotplot(alldata, n_genes=5, groupby = 'leiden')

alldata = alldata[alldata.obs_names.sort_values(), :]

alldata.obs.to_csv('tutorial_data/046metadataUniform.csv')  # Exporting METADATA file

alldata.write('tutorial_data/XETG00054_0006530_046_Uniform.h5ad')  # String procesed AnnData

###########################   End of Analysis ######################################################


