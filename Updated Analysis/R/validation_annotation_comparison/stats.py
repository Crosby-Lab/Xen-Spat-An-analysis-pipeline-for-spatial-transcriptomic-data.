import numpy as np, pandas as pd
df = pd.read_pickle('s1top_final.pkl')
TUM = ['Tumor_Cells','Proliferative_Tumor_Cells']
is_t = df['leiden'].isin(TUM).to_numpy(); T = is_t.sum()
nic = df['niche_recovered'].to_numpy()

rows=[]
for c in np.unique(nic):
    m = nic==c
    rows.append((c, m.sum(), (m&is_t).sum(), 2*(m&is_t).sum()/(T+m.sum())))
d = pd.DataFrame(rows, columns=['niche','n_cells','n_tumor','dice']).sort_values('dice',ascending=False)
d['rank']=range(1,len(d)+1)
gaps = -np.diff(d['dice'].to_numpy())
o = np.argsort(-gaps); s1,s2 = o[0],o[1]
print(d.to_string(index=False, float_format=lambda v:'%.3f'%v))
print('\nlargest gap %.3f after rank %d | second %.3f after rank %d | gap ratio %.2f'
      %(gaps[s1],s1+1,gaps[s2],s2+1,gaps[s1]/gaps[s2]))
tum_n = set(d['niche'].to_numpy()[:s1+1]); alt_n = set(d['niche'].to_numpy()[:s2+1])
print('tumor compartment niches:', sorted(tum_n), ' at second gap:', sorted(alt_n))
df['compartment'] = np.where(np.isin(nic,list(tum_n)),'tumor','stromal')
df['compartment_alt'] = np.where(np.isin(nic,list(alt_n)),'tumor','stromal')

order=['Solid DCIS','DCIS-comedo','Tumor with immune infiltrates','Ductal hyperplasia',
       'Dilated ducts','Normal breast epithelium','Unannotated tissue']
res=[]
for r in order:
    g = df[df.region==r]
    res.append(dict(region=r, n_cells=len(g),
        pct_tumor_compartment=100*(g.compartment=='tumor').mean(),
        pct_at_second_gap=100*(g.compartment_alt=='tumor').mean(),
        pct_tumor_labelled=100*g['leiden'].isin(TUM).mean()))
res=pd.DataFrame(res)
print('\n', res.to_string(index=False, float_format=lambda v:'%.1f'%v))
print('\ntumor-labelled cells captured by the tumor compartment: %.1f%%'
      %(100*df.loc[is_t,'compartment'].eq('tumor').mean()))
print('tumor compartment is %.1f%% of all cells'%(100*df.compartment.eq('tumor').mean()))
comp = df.groupby(['compartment','leiden']).size().unstack(fill_value=0)
comp = (100*comp.T/comp.sum(1)).T
print('\ncomposition (%):\n', comp.round(1).to_string())
res.to_csv('S1_TOP_region_vs_compartment.csv',index=False)
d.to_csv('S1_TOP_dice_gap_profile.csv',index=False)
comp.round(2).to_csv('S1_TOP_compartment_composition.csv')
df[['Cell_ID','Cellordinate.x','Cellordinate.y','leiden','niche_recovered','region','compartment','compartment_alt']].to_csv('S1_TOP_cells_annotated.csv',index=False)
