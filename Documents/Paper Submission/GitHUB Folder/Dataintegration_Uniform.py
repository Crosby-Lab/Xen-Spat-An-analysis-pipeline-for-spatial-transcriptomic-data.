###############################################################################
### Multiple Data Integration
### Author : Jesuchristopher Joseph (jesu.joseph@duke.edu)
### Date :
### Environment : Python 3.9.0
### Packages  : scanpy 1.9.3, matplotlib 3.7.2, bbknn 1.6.0
#######################################################################################################################

import scanpy as sc
import matplotlib.pyplot as pl
import bbknn


sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor= 'white')

adata1 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006530_046_Uniform_72925.h5ad')

adata1.var_names_make_unique()

adata2 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006674_047_Uniform_72925.h5ad')

adata2.var_names_make_unique()

adata3 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006530_049_Uniform_72925.h5ad')

adata3.var_names_make_unique()

adata4 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006674_089_Uniform_72925.h5ad')

adata4.var_names_make_unique()

adata5 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006530_097_Uniform_72925.h5ad')

adata5.var_names_make_unique()

adata6 = sc.read_h5ad('tutorial_data/xenium_data/BC_Uniform/XETG00054_0006674_100_Uniform_72925.h5ad')

adata6.var_names_make_unique()


adata = adata1.concatenate([adata2, adata3, adata4, adata5, adata6], batch_categories = ['BC1','BC2','BC3','BC4','BC5','BC6'])

sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca_variance_ratio(adata, log=True)

sc.external.pp.bbknn(adata,batch_key = 'batch')  # Data Integration

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution = 1.0)

sc.pl.umap(adata, color=['leiden','batch'])

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, fontsize =14)


new_cluster_names = [
    'Tumor Cells_TACSTD2_HER2_3','Macrophages','T Cells','Fibroblasts','Tumor Cells_GATA3','Endothelial Cells','Tumor Cells_SMAD3',
    'Basal/Myoepithelial Cells','Tumor Cells_TOP2A','Tregs','Plasma Cells', 'Adipocytes','B Cells','CD8 T Cells','Mast Cells','Plasmacytoid DCs','Dendritic Cells']


adata.rename_categories('leiden', new_cluster_names)

sc.pl.rank_genes_groups_dotplot(adata, n_genes=8, groupby = 'leiden')


adata.obs.to_csv('tutorial_data/DataIntegration_PostmetadataUniform_103025.csv')
adata.write('tutorial_data/Integrated_Data_CelltypeAnnotated_9225.h5ad')