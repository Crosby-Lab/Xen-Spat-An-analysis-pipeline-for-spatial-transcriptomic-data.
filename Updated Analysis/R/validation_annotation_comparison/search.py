import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

TGT = np.array(sorted([0.01,0.00,0.01,0.02,0.45,0.00,0.00,0.75,0.01,0.09,0.00], reverse=True))
df = pd.read_pickle('s1top.pkl')
xy = df[['Cellordinate.x','Cellordinate.y']].to_numpy(float)
lab = df['leiden'].astype(str).to_numpy()
cats = pd.unique(lab)
TUM = ['Tumor_Cells','Proliferative_Tumor_Cells']
is_t = np.isin(lab, TUM); T = is_t.sum()

_, idx = cKDTree(xy).query(xy, k=25)
onehot = np.zeros((len(df), len(cats)), np.float32)
code = {c:i for i,c in enumerate(cats)}
onehot[np.arange(len(df)), [code[c] for c in lab]] = 1
summ = onehot[idx].sum(1)
z = (summ - summ.mean(0))/summ.std(0, ddof=1)
z = np.clip(z, -10, 10)                     # Seurat ScaleData scale.max

def dvec(cl):
    out=[]
    for c in np.unique(cl):
        m = cl==c
        out.append(2*(m&is_t).sum()/(T+m.sum()))
    return np.array(sorted(out, reverse=True))

best=None
for seed in range(40):
    km = KMeans(n_clusters=11, n_init=10, random_state=seed).fit(z)
    v = dvec(km.labels_)
    err = np.abs(np.round(v,2)-TGT).sum()
    if best is None or err < best[0]:
        best=(err, v, km.labels_.copy(), seed, km.inertia_)
        print('seed %2d  L1 err %.3f  inertia %.0f  %s' % (seed, err, km.inertia_, np.round(v,2)))
print('\nbest L1 error %.3f at seed %d' % (best[0], best[3]))
print('recovered', np.round(best[1],2))
print('Table S8 ', TGT)
np.save('clusters_best.npy', best[2])
df['niche_recovered']=best[2]+1
df.to_pickle('s1top_with_niches.pkl')
