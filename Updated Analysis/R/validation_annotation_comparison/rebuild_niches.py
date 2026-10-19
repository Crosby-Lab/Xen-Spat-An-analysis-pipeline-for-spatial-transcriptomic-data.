"""Reimplementation of Seurat v5 BuildNicheAssay, to recover the S1_TOP niche
assignment from the saved metadata, and check it against Table S8."""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

TARGET = {  # Table S8, S1_TOP, niches 1..11
 1:0.01, 2:0.00, 3:0.01, 4:0.02, 5:0.45, 6:0.00,
 7:0.00, 8:0.75, 9:0.01, 10:0.09, 11:0.00}

df = pd.read_pickle('s1top.pkl')
xy = df[['Cellordinate.x','Cellordinate.y']].to_numpy(float)
lab = df['leiden'].astype(str).to_numpy()
cats = pd.unique(lab)                       # Seurat uses unique() order
print('%d cells, %d label values' % (len(df), len(cats)))

K = 25
_, idx = cKDTree(xy).query(xy, k=K)         # includes self, as Seurat's does
onehot = np.zeros((len(df), len(cats)), dtype=np.float32)
code = {c:i for i,c in enumerate(cats)}
onehot[np.arange(len(df)), [code[c] for c in lab]] = 1
summ = onehot[idx].sum(1)                   # neighbour cell-type counts
print('niche feature matrix', summ.shape)

# Seurat ScaleData: centre and scale each feature across cells
z = (summ - summ.mean(0)) / summ.std(0, ddof=1)

def dice_vector(clusters, tumor_labels):
    is_t = np.isin(lab, tumor_labels)
    T = is_t.sum()
    out = {}
    for c in np.unique(clusters):
        m = clusters == c
        out[c] = 2*(m & is_t).sum() / (T + m.sum())
    return out

for tum in (['Tumor_Cells'], ['Tumor_Cells','Proliferative_Tumor_Cells']):
    best = None
    for seed in range(6):
        km = KMeans(n_clusters=11, n_init=30, random_state=seed).fit(z)
        d = dice_vector(km.labels_, tum)
        v = sorted(d.values(), reverse=True)
        if best is None or km.inertia_ < best[0]:
            best = (km.inertia_, v, km.labels_.copy(), d)
    v = [round(x,2) for x in best[1]]
    print('\ntumor labels = %s' % tum)
    print('  recovered Dice (sorted):', v)
    print('  Table S8     (sorted):', sorted([round(x,2) for x in TARGET.values()], reverse=True))
    np.save('clusters_%d.npy' % len(tum), best[2])
