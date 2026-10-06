import numpy as np, pandas as pd, scipy.sparse as sp, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})
C=sp.load_npz('counts_recon.npz').tocsr(); obs=pd.read_csv('obs_with_dist.csv.gz'); obs['sample']=obs['sample'].astype(str).str.zfill(3)
genes=pd.read_csv('/mnt/user-data/uploads/Downloads/xenspat_B6_export/genes.csv',header=None)[0].tolist(); gi={g:i for i,g in enumerate(genes)}
EPI=[g for g in ['EPCAM','KRT7','KRT8','KRT19','CDH1','GRHL2','DSP','AR','ESR1','ERBB2','VTCN1'] if g in gi]
samples=['046','047','049','089','097','100']
bins=np.arange(0,105,5); rows=[]
for ct in ['Macrophages','Fibroblasts']:
    for s in samples:
        m=(obs.leiden==ct)&(obs['sample']==s)
        idx=np.where(m)[0]; d=obs.dist.values[idx]
        sub=C[idx][:,[gi[g] for g in EPI]].toarray()
        b=np.digitize(d,bins)-1
        for k in range(len(bins)-1):
            sel=b==k
            if sel.sum()<10: continue
            mm=sub[sel].mean(0)
            for g,v in zip(EPI,mm): rows.append((s,ct,bins[k],sel.sum(),g,v))
            rows.append((s,ct,bins[k],sel.sum(),'epithelial_sum',sub[sel].sum(1).mean()))
ep=pd.DataFrame(rows,columns=['sample','cell_type','distance_bin_um','n_cells','gene','mean_count']); ep['mean_count']=ep.mean_count.astype(float); ep['distance_bin_um']=ep.distance_bin_um.astype(float)
ep.to_csv('B6_epithelial_vs_distance.csv',index=False)
fig,axs=plt.subplots(2,2,figsize=(6.5,4.6),dpi=300,sharex=True)
cols=dict(zip(samples,['#1B7F79','#6A4C93','#C8553D','#3D7A2E','#1F3864','#BF9000']))
for j,ct in enumerate(['Macrophages','Fibroblasts']):
    for i,g in enumerate(['epithelial_sum','VTCN1']):
        ax=axs[i,j]
        for s in samples:
            q=ep[(ep.cell_type==ct)&(ep.gene==g)&(ep['sample']==s)]
            ax.plot(q.distance_bin_um+2.5,q.mean_count,color=cols[s],lw=1,label=s)
        ax.set_title(f'{ct}: '+('sum of epithelial genes' if g=='epithelial_sum' else 'VTCN1'),fontsize=7.5)
        ax.tick_params(labelsize=6); ax.axvline(25,color='grey',ls=':',lw=0.7)
        for sp_ in ['top','right']: ax.spines[sp_].set_visible(False)
        if j==0: ax.set_ylabel('Mean transcripts per cell',fontsize=6.5)
        if i==1: ax.set_xlabel('Distance to nearest tumor cell (µm)',fontsize=6.5)
axs[0,0].legend(fontsize=5.5,frameon=False,ncol=2)
fig.tight_layout(); fig.savefig('B6_Fig_epithelial_vs_distance.png'); fig.savefig('B6_Fig_epithelial_vs_distance.svg')
# summary: ratio near(0-10) vs >50
for ct in ['Macrophages','Fibroblasts']:
    q=ep[(ep.cell_type==ct)&(ep.gene=='epithelial_sum')]
    a=q[q.distance_bin_um<10].groupby('sample').mean_count.mean(); b=q[q.distance_bin_um>=50].groupby('sample').mean_count.mean()
    print(ct,'0-10um',a.round(2).to_dict(),' >=50um',b.round(2).to_dict())
print('epi genes',EPI)
