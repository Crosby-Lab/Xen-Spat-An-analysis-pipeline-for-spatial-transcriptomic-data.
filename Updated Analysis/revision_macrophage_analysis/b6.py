import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.spatial import cKDTree
D='/mnt/user-data/uploads/Downloads/xenspat_B6_export/'
C=sp.load_npz('counts_recon.npz').tocsr(); obs=pd.read_csv(D+'obs.csv.gz'); genes=pd.read_csv(D+'genes.csv',header=None)[0].tolist()
smap={'BC1':'046','BC2':'047','BC3':'049','BC4':'089','BC5':'097','BC6':'100'}
obs['sample']=obs.batch.map(smap)
obs['is_tumor']=obs.leiden.str.startswith('Tumor Cells')
obs['dist']=np.nan
for s,df in obs.groupby('sample'):
    t=df[df.is_tumor]; tree=cKDTree(t[['x_centroid','y_centroid']].values)
    nt=df[~df.is_tumor]; d,_=tree.query(nt[['x_centroid','y_centroid']].values,k=1)
    obs.loc[nt.index,'dist']=d
obs['zone']=pd.cut(obs.dist,[-0.01,25,75,np.inf],labels=['near','margin','far'])
obs.to_csv('obs_with_dist.csv.gz',index=False)
# counts of macrophages by zone
tab=obs[obs.leiden=='Macrophages'].groupby(['sample','zone'],observed=False).size().unstack()
print(tab)
# enrichment check like Fig 6e for macrophages
print((tab.T/tab.sum(1)).T.round(2))
