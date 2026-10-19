"""Visual conventions taken from the published Xen-Spat figures."""
import matplotlib
from matplotlib import font_manager as fm

INK   = '#2A3F5F'           # axis text / labels in Figures 3, 4, 7, 8
TUMOR = '#EC2227'           # red, as used for tumor/overlap marks in Figures 4h and 7e-g
STROMA= '#A0A0A0'

# cell type colours sampled from the Figure 7 legend
CELLTYPE = {
 'Tumor_Cells'                        : '#F3A460',
 'Proliferative_Tumor_Cells'          : '#F2634B',
 'Tumor_Adjacent_Myoepithelial_Cells' : '#F9C0CB',
 'Endothelial_Cells'                  : '#EC2227',
 'Perivascular_Cells'                 : '#7F8133',
 'Normal_Epithelial_Cells'            : '#3954A4',
 'Normal_Myoepithelial_Cells'         : '#FAA41C',
 'Macrophages'                        : '#9CD089',
 'Plasma_Cells'                       : '#7F479C',
 'T_Lymphocytes'                      : '#7F1518',
 'Stromal_Cells'                      : '#BB8D8E',
 'not_assigned'                       : '#DCDCDC',
}
# niche colour vector, identical across samples in the published figures
NICHE = {1:'#F4A460',2:'#808080',3:'#FF0000',4:'#808000',5:'#FFFF00',6:'#90EE90',
         7:'#9400D3',8:'#800000',9:'#0000FF',10:'#FF6347',11:'#FFC0CB'}

def use_paper_style():
    for f in ('Carlito-Regular.ttf','Carlito-Bold.ttf'):
        try: fm.fontManager.addfont('/usr/share/fonts/truetype/crosextra/'+f)
        except Exception: pass
    matplotlib.rcParams.update({
        'font.family':'Carlito',
        'text.color':INK, 'axes.labelcolor':INK,
        'xtick.color':INK, 'ytick.color':INK,
        'axes.edgecolor':INK, 'axes.linewidth':0.8,
        'xtick.major.width':0.8,'ytick.major.width':0.8,
        'figure.facecolor':'white','axes.facecolor':'white',
    })
