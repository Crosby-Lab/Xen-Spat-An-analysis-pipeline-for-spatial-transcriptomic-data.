import pandas as pd, numpy as np, json
U='/mnt/user-data/uploads/Local backup/Documents/RCode/Brestcancer/'
SEC={
 'S1_TOP':dict(src='recon', tum=['Tumor_Cells','Proliferative_Tumor_Cells']),
 'S2_TOP':dict(src=U+'XeniumdataNiche_S2_TOP_82726_25.csv', tum=None),
 'S3_TOP':dict(src=U+'XeniumdataNiche_S3_TOP_82726_25.csv', tum=['Epithelial_&_Tumor_Cells']),
 'S4_TOP':dict(src=U+'XeniumdataNiche_S4_TOP_82726_25.csv', tum=['Tumor_Cells','Proliferative_Tumor Cells']),
}
def profile(dice):
    s=sorted(dice.items(), key=lambda kv:-kv[1])
    v=[x[1] for x in s]; g=[v[i]-v[i+1] for i in range(len(v)-1)]
    o=sorted(range(len(g)), key=lambda i:-g[i])
    return s,g,o[0],o[1],(g[o[0]]/g[o[1]] if g[o[1]]>0 else np.inf)

rows=[]
for sec,cfg in SEC.items():
    if cfg['src']=='recon':
        d=pd.read_pickle('s1top_final.pkl').rename(columns={'niche_recovered':'niches'})
    else:
        d=pd.read_csv(cfg['src'])
    labs=d.leiden.value_counts()
    print('===',sec,' cells %d  niches %d  labels %d'%(len(d),d.niches.nunique(),len(labs)))
    print('   labels:', list(labs.index))
    if cfg['tum'] is None:
        print('   NO TUMOR LABEL -> no partition\n'); 
        rows.append(dict(section=sec,n_cells=len(d),n_niches=d.niches.nunique(),tumor_label='none',
                         gap_ratio=np.nan,n_tumor_niches=np.nan,capture_sel=np.nan,purity_sel=np.nan,
                         capture_alt=np.nan,purity_alt=np.nan)); continue
    it=d.leiden.isin(cfg['tum']); T=it.sum()
    dice={n: 2*(it&(d.niches==n)).sum()/(T+(d.niches==n).sum()) for n in sorted(d.niches.unique())}
    s,g,i1,i2,ratio=profile(dice)
    sel=[k for k,_ in s[:i1+1]]; alt=[k for k,_ in s[:i2+1]]
    cs=d.niches.isin(sel); ca=d.niches.isin(alt)
    r=dict(section=sec,n_cells=len(d),n_niches=d.niches.nunique(),
           tumor_label='+'.join(cfg['tum']), n_tumor_cells=int(T),
           dice_sorted=[round(x[1],3) for x in s],
           gap1=round(g[i1],3), gap2=round(g[i2],3), gap_ratio=round(ratio,2),
           n_tumor_niches=len(sel), n_tumor_niches_alt=len(alt),
           capture_sel=round(100*it[cs].sum()/T,1), purity_sel=round(100*it[cs].mean(),1),
           capture_alt=round(100*it[ca].sum()/T,1), purity_alt=round(100*it[ca].mean(),1),
           pct_cells_in_tumor=round(100*cs.mean(),1))
    rows.append(r)
    print('   tumor label(s):', r['tumor_label'], '| %d tumor cells (%.1f%%)'%(T,100*T/len(d)))
    print('   dice:', r['dice_sorted'])
    print('   gap1 %.3f gap2 %.3f ratio %.2f | %d niche(s) selected, %d at second gap'
          %(g[i1],g[i2],ratio,len(sel),len(alt)))
    print('   capture %.1f%% purity %.1f%%   |  alt capture %.1f%% purity %.1f%%\n'
          %(r['capture_sel'],r['purity_sel'],r['capture_alt'],r['purity_alt']))
pd.DataFrame(rows).to_csv('validation_cohort_gap_summary.csv',index=False)
json.dump(rows,open('validation_cohort_gap_summary.json','w'),indent=1,default=str)
