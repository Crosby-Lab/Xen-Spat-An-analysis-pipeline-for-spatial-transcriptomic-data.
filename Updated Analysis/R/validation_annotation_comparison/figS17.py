import numpy as np, pandas as pd, json, sys, matplotlib
matplotlib.use('Agg')
sys.path.insert(0,'/home/claude/r2/s1')
from PAPER_STYLE import use_paper_style, INK, TUMOR, STROMA, CELLTYPE
use_paper_style()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

df=pd.read_pickle('s1top_final.pkl')
TUM=['Tumor_Cells','Proliferative_Tumor_Cells']
df['compartment']=np.where(df.niche_recovered.isin([7,2]),'tumor','stromal')
df['compartment_alt']=np.where(df.niche_recovered.isin([7]),'tumor','stromal')
j=json.load(open('annotation_polygons_S1_TOP_um.json'))
V=json.load(open('validation_cohort_gap_summary.json'))

def panel_letter(ax,s,pad=6):
    ax.set_title(s, loc='left', fontsize=13, fontweight='bold', color=INK, pad=pad)

fig=plt.figure(figsize=(7.2,8.1),dpi=300)
gs=fig.add_gridspec(3,2,height_ratios=[1.05,1.0,0.95],width_ratios=[1.12,1.0],
                    hspace=0.78,wspace=0.42)

# ---- a: spatial ----
ax=fig.add_subplot(gs[0,:])
st=df[df.compartment=='stromal']; tu=df[df.compartment=='tumor']
ax.scatter(st['Cellordinate.x'],st['Cellordinate.y'],s=0.12,c=STROMA,lw=0,alpha=.45,rasterized=True)
ax.scatter(tu['Cellordinate.x'],tu['Cellordinate.y'],s=0.12,c=TUMOR,lw=0,alpha=.55,rasterized=True)
for r in j['regions']:
    v=np.array(r['vertices_xenium_um']); v=np.vstack([v,v[:1]])
    ax.plot(v[:,0],v[:,1],color=INK,lw=0.6)
ax.set_aspect('equal'); ax.invert_yaxis()
ax.set_xlabel('x (µm)',fontsize=8.5); ax.set_ylabel('y (µm)',fontsize=8.5)
ax.tick_params(labelsize=8)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=TUMOR,ms=4,label='Tumor compartment (niches 7 and 2)'),
                   Line2D([],[],marker='o',ls='',color=STROMA,ms=4,label='Stromal compartment'),
                   Line2D([],[],color=INK,lw=0.9,label='Annotated region (Figure 7a)')],
          fontsize=7.5,frameon=False,loc='lower left',bbox_to_anchor=(0.085,1.0),ncol=3,
          handletextpad=0.4,columnspacing=1.6,labelcolor=INK)
panel_letter(ax,'a')

# ---- b: region coverage ----
order=['DCIS-comedo','Solid DCIS','Tumor with immune infiltrates','Ductal hyperplasia',
       'Dilated ducts','Normal breast epithelium']
sel=[100*df[df.region==r].compartment.eq('tumor').mean() for r in order]
alt=[100*df[df.region==r].compartment_alt.eq('tumor').mean() for r in order]
ax2=fig.add_subplot(gs[1,0]); y=np.arange(len(order))[::-1]
ax2.barh(y+0.19,sel,height=0.36,color=TUMOR,label='Selected partition (largest gap)')
ax2.barh(y-0.19,alt,height=0.36,color='none',edgecolor=TUMOR,lw=0.9,
         label='Alternative (second-largest gap)')
for yy,a,b in zip(y,sel,alt):
    ax2.text(a+1.5,yy+0.19,'%.0f'%a,va='center',fontsize=7,color=INK)
    ax2.text(b+1.5,yy-0.19,'%.0f'%b,va='center',fontsize=7,color=INK)
ax2.set_yticks(y); ax2.set_yticklabels(order,fontsize=7.5)
ax2.set_xlim(0,100); ax2.set_xlabel('Cells in the tumor compartment (%)',fontsize=8)
ax2.tick_params(axis='x',labelsize=8)
for s in ('top','right'): ax2.spines[s].set_visible(False)
ax2.legend(fontsize=7,frameon=False,loc='lower right',bbox_to_anchor=(1.02,-0.52),labelcolor=INK)
panel_letter(ax2,'b')

# ---- c: composition by published cell type label ----
ax3=fig.add_subplot(gs[1,1])
keep=['Tumor_Cells','Proliferative_Tumor_Cells','T_Lymphocytes','Macrophages',
      'Stromal_Cells','Perivascular_Cells','Endothelial_Cells','Normal_Epithelial_Cells']
groups=[('Niche 7\nDice 0.75',df[df.niche_recovered==7]),
        ('Niche 2\nDice 0.46',df[df.niche_recovered==2]),
        ('Stromal\ncompartment',df[df.compartment=='stromal'])]
yy=np.arange(len(groups))[::-1]
for gi,(nm,g) in zip(yy,groups):
    pct=(100*g['leiden'].value_counts()/len(g)); left=0
    for k in keep+['Other']:
        v=pct.get(k,0) if k!='Other' else 100-sum(pct.get(z,0) for z in keep)
        if v<=0: continue
        ax3.barh(gi,v,left=left,height=0.52,color=CELLTYPE.get(k,'#DCDCDC'),
                 edgecolor='white',lw=0.4)
        if v>=12:
            dark = k in ('T_Lymphocytes','Normal_Epithelial_Cells','Plasma_Cells','Endothelial_Cells')
            ax3.text(left+v/2,gi,'%.0f'%v,ha='center',va='center',fontsize=7,
                     color='white' if dark else INK)
        left+=v
ax3.set_yticks(yy); ax3.set_yticklabels([n for n,_ in groups],fontsize=7.5)
ax3.set_xlim(0,100); ax3.set_xlabel('Cells (%)',fontsize=8); ax3.tick_params(axis='x',labelsize=8)
for s in ('top','right'): ax3.spines[s].set_visible(False)
lbl={'Tumor_Cells':'Tumor cells','Proliferative_Tumor_Cells':'Proliferative tumor cells',
     'T_Lymphocytes':'T lymphocytes','Macrophages':'Macrophages','Stromal_Cells':'Stromal cells',
     'Other':'Other'}
ax3.legend(handles=[plt.Rectangle((0,0),1,1,color=CELLTYPE.get(k,'#DCDCDC')) for k in lbl],
           labels=list(lbl.values()),fontsize=7,frameon=False,loc='lower left',
           bbox_to_anchor=(-0.05,-0.74),ncol=2,handlelength=1.1,handletextpad=0.4,
           columnspacing=1.0,labelcolor=INK)
panel_letter(ax3,'c')

# ---- d: capture vs purity ----
ax4=fig.add_subplot(gs[2,:])
for r in V:
    if r.get('tumor_label')=='none' or r.get('gap_ratio') is None: continue
    cs,ps=float(r['capture_sel']),float(r['purity_sel'])
    ca,pa=float(r['capture_alt']),float(r['purity_alt'])
    ax4.annotate('',xy=(cs,ps),xytext=(ca,pa),
                 arrowprops=dict(arrowstyle='-|>',color=INK,lw=0.9,alpha=.55,shrinkA=4,shrinkB=4))
    ax4.plot(ca,pa,'o',ms=6,mfc='white',mec=TUMOR,mew=1.3,zorder=3)
    ax4.plot(cs,ps,'o',ms=6.5,color=TUMOR,mec='white',mew=0.6,zorder=4)
    ax4.annotate('%s  (gap ratio %.2f)'%(r['section'],float(r['gap_ratio'])),(cs,ps),
                 textcoords='offset points',xytext=(8,-3),fontsize=7.5,color=INK)
ax4.set_xlim(35,100); ax4.set_ylim(60,100)
ax4.set_xlabel('Tumor-labeled cells captured by the tumor compartment (%)',fontsize=8)
ax4.set_ylabel('Tumor-labeled cells in the\ntumor compartment (%)',fontsize=8)
ax4.tick_params(labelsize=8)
for s in ('top','right'): ax4.spines[s].set_visible(False)
ax4.legend(handles=[Line2D([],[],marker='o',ls='',color=TUMOR,ms=6,label='Selected partition (largest gap)'),
                    Line2D([],[],marker='o',ls='',mfc='white',mec=TUMOR,mew=1.3,ms=6,label='Alternative (second-largest gap)')],
           fontsize=7.5,frameon=False,loc='lower left',handletextpad=0.5,labelcolor=INK)
panel_letter(ax4,'d')

fig.savefig('Figure_S17.png',bbox_inches='tight',facecolor='white')
fig.savefig('Figure_S17.pdf',bbox_inches='tight',facecolor='white')
print('saved')
