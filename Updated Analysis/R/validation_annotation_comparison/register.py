import numpy as np, pandas as pd, cv2, json
from scipy import ndimage as ndi
from matplotlib.path import Path

df = pd.read_pickle('s1top_with_niches.pkl')
xy = df[['Cellordinate.x','Cellordinate.y']].to_numpy(float)
UM = 8.0                                  # microns per raster pixel
x0,y0 = xy.min(0) - 40
W = int((xy[:,0].max()+40 - x0)/UM)+1; H = int((xy[:,1].max()+40 - y0)/UM)+1
cell = np.zeros((H,W), np.uint8)
ix = ((xy[:,0]-x0)/UM).astype(int); iy = ((xy[:,1]-y0)/UM).astype(int)
cell[iy,ix] = 255
cell = ndi.binary_closing(cell>0, np.ones((7,7)))
cell = ndi.binary_fill_holes(cell)
lab,n = ndi.label(cell); sizes = ndi.sum(cell,lab,range(1,n+1))
cell = (lab == np.argmax(sizes)+1).astype(np.uint8)*255
print('cell raster', cell.shape, 'area', (cell>0).sum())

he = (np.load('../ma.npy')*255).astype(np.uint8)
print('H&E mask', he.shape, 'area', (he>0).sum())

HH = max(cell.shape[0], he.shape[0])+60; WW = max(cell.shape[1], he.shape[1])+60
def pad(m):
    o = np.zeros((HH,WW), np.uint8); o[:m.shape[0],:m.shape[1]] = m; return o
C, Hm = pad(cell), pad(he)
def bb(m):
    ys,xs = np.nonzero(m); return xs.min(),xs.max(),ys.min(),ys.max()
cx0,cx1,cy0,cy1 = bb(C); hx0,hx1,hy0,hy1 = bb(Hm)
sx = (cx1-cx0)/(hx1-hx0); sy = (cy1-cy0)/(hy1-hy0)
M = np.array([[1/sx,0,hx0-cx0/sx],[0,1/sy,hy0-cy0/sy]], np.float32)   # cell-frame -> H&E frame
F = cv2.WARP_INVERSE_MAP|cv2.INTER_NEAREST
iou = lambda M: (((cv2.warpAffine(Hm,M,(WW,HH),flags=F)>0)&(C>0)).sum() /
                 ((cv2.warpAffine(Hm,M,(WW,HH),flags=F)>0)|(C>0)).sum())
print('IoU at bbox init %.4f  (aniso sx %.3f sy %.3f)'%(iou(M),sx,sy))
for sig,bl in ((9,31),(4,21),(2,11)):
    a = cv2.GaussianBlur(C.astype(np.float32)/255,(bl,bl),sig)
    b = cv2.GaussianBlur(Hm.astype(np.float32)/255,(bl,bl),sig)
    cc,M = cv2.findTransformECC(a,b,M.copy(),cv2.MOTION_AFFINE,
          (cv2.TERM_CRITERIA_EPS|cv2.TERM_CRITERIA_COUNT,3000,1e-8),None,3)
    print('  sigma %2d ecc %.4f IoU %.4f'%(sig,cc,iou(M)))
# M maps cell-raster px -> H&E px; invert to go H&E px -> cell px -> microns
Mi = cv2.invertAffineTransform(M)
T = np.zeros((2,3))
T[0,0]=UM*Mi[0,0]; T[0,1]=UM*Mi[0,1]; T[0,2]=UM*Mi[0,2]+x0
T[1,0]=UM*Mi[1,0]; T[1,1]=UM*Mi[1,1]; T[1,2]=UM*Mi[1,2]+y0
print('H&E px -> microns:\n', np.round(T,4))
np.save('T_he_to_um.npy', T); np.save('final_iou.npy', np.array([iou(M)]))

j = json.load(open('../overlay/annotation_polygons_S1_TOP.json'))
for r in j['regions']:
    v = np.array(r['vertices_panelA_px'], float)
    r['vertices_xenium_um'] = np.round((T[:,:2]@v.T).T + T[:,2], 1).tolist()
j['transform_panelA_px_to_xenium_um'] = T.tolist()
j['registration'] = {'method':'affine ECC, H&E tissue outline to the rasterised cell coordinates',
                     'tissue_outline_IoU': round(float(iou(M)),4),
                     'raster_um_per_px': UM}
json.dump(j, open('annotation_polygons_S1_TOP_um.json','w'), indent=1)

region = np.full(len(df), 'Unannotated tissue', dtype=object)
order = ['Dilated ducts','Normal breast epithelium','Ductal hyperplasia',
         'Solid DCIS','DCIS-comedo','Tumor with immune infiltrates']
for cls in order:
    for r in j['regions']:
        if r['class'] != cls: continue
        p = Path(np.array(r['vertices_xenium_um']))
        inside = p.contains_points(xy)
        region[inside] = cls
df['region'] = region
print(df['region'].value_counts())
df.to_pickle('s1top_final.pkl')
