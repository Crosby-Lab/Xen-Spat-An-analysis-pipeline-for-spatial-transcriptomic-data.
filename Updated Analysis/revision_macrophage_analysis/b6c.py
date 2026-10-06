import numpy as np, pandas as pd, scipy.sparse as sp, warnings, logging
warnings.filterwarnings('ignore'); logging.disable(logging.WARNING)
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
C=sp.load_npz('counts_recon.npz').tocsr(); obs=pd.read_csv('obs_with_dist.csv.gz'); obs['sample']=obs['sample'].astype(str).str.zfill(3)
genes=pd.read_csv('/mnt/user-data/uploads/Downloads/xenspat_B6_export/genes.csv',header=None)[0].tolist()
samples=['046','047','049','089','097','100']
# tumor-enriched flag: mean normalized expression in tumor cells vs macrophages >50 um
lib=np.asarray(C.sum(1)).ravel(); N=sp.diags(1e3/lib)@C
flags=[]
for s in samples:
    t=(obs['sample']==s)&obs.leiden.str.startswith('Tumor Cells')
    m=(obs['sample']==s)&(obs.leiden=='Macrophages')&(obs.dist>50)
    mt=np.asarray(N[np.where(t)[0]].mean(0)).ravel(); mm=np.asarray(N[np.where(m)[0]].mean(0)).ravel()
    flags.append(np.log2((mt+0.01)/(mm+0.01))>1)
nflag=np.sum(flags,0); tumor_enriched=[g for g,n in zip(genes,nflag) if n>=4]
pd.DataFrame({'gene':genes,'n_samples_tumor_enriched':nflag}).sort_values('n_samples_tumor_enriched',ascending=False).to_csv('B6_tumor_enriched_genes.csv',index=False)
print('tumor-enriched genes (>=4/6):',len(tumor_enriched))
# pseudobulk
rows=[];meta=[]
for s in samples:
    for z in ['near','margin']:
        m=(obs['sample']==s)&(obs.leiden=='Macrophages')&(obs.zone==z)
        rows.append(np.asarray(C[np.where(m)[0]].sum(0)).ravel()); meta.append((f'{s}_{z}',s,z,m.sum()))
cnt=pd.DataFrame(rows,columns=genes,index=[m[0] for m in meta]).astype(int)
md=pd.DataFrame(meta,columns=['id','patient','zone','n_cells']).set_index('id')
md.to_csv('B6_pseudobulk_metadata.csv')
def run(cols,tag):
    dds=DeseqDataSet(counts=cnt[cols],metadata=md[['patient','zone']],design='~patient + zone',refit_cooks=True,inference=DefaultInference(n_cpus=4),quiet=True)
    dds.deseq2()
    st=DeseqStats(dds,contrast=['zone','near','margin'],inference=DefaultInference(n_cpus=4),quiet=True); st.summary()
    r=st.results_df.sort_values('padj'); r.index.name='gene'; r.to_csv(f'B6_pseudobulk_DE_{tag}.csv'); return r
allr=run(genes,'all_genes'); exc=run([g for g in genes if g not in tumor_enriched],'excluded')
for name,r in [('ALL',allr),('EXCLUDED',exc)]:
    sig=r[(r.padj<0.05)]
    print(name,'n sig',len(sig),'up near',(sig.log2FoldChange>0).sum(),'down',(sig.log2FoldChange<0).sum())
    print(sig.head(25)[['baseMean','log2FoldChange','padj']].round(3).to_string())
for g in ['GRHL2','AR','DSP','VTCN1','MRC1','CD163','CD14','CXCL10','CCND1','PRLR']:
    print(g,'tumor_enriched' if g in tumor_enriched else '', allr.loc[g,['log2FoldChange','padj']].round(4).to_dict() if g in allr.index else '', exc.loc[g,['log2FoldChange','padj']].round(4).to_dict() if g in exc.index else 'excluded')
print('tumor enriched list:',tumor_enriched)
