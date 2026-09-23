"""Render genomic panels from regenerated tables, preserving mutation identity.

Notebook selections are retained; source-verified annotation corrections are
applied to labels. A short label is never an
aggregation key: variants at distinct coordinates remain separate rows.
"""
import ast
import hashlib
import json
import os
import re
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path.home() / '.cache' / 'matplotlib'))
os.environ['MPLBACKEND'] = 'Agg'
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/genomics/generated_tables'
OUT = ROOT/'data/genomics/plot_tables'
FIG = ROOT/'working/analysis/stats-rework/out/genomic-panels'
NB = ROOT/'working/figures'
ENDPOINTS = {
 'Fig2f': ('MG_AMI','Figure 2/F-H heatmaps/05-PA.ipynb',13,'Blues'),
 'Fig2g': ('MG_LEV','Figure 2/F-H heatmaps/07-PL.ipynb',14,'Reds'),
 'Fig2h': ('MG_CEF','Figure 2/F-H heatmaps/06-PC.ipynb',18,'Greens'),
 'S5a': ('MG_AMI','Figure 2/F-H heatmaps/05-PA.ipynb',12,'Blues'),
 'S5b': ('MG_LEV','Figure 2/F-H heatmaps/07-PL.ipynb',13,'Reds'),
 'S5c': ('MG_CEF','Figure 2/F-H heatmaps/06-PC.ipynb',16,'Greens'),
 'S3e': ('MG_CEF_R','Supplemental Figure 3 - CefR/E - genetics/08-PCr.ipynb',11,'Purples'),
 'S4f': ('PbEc_CEF','Supplemental Figure 4 - Pb/F - Mutations/10-ATEC12.ipynb',16,'Greys'),
 'Fig4c': ('PLAC','Figure 4/02-PLAC.ipynb',29,'Greens'),
}
CROSSWALK=[]
CHANGES=[]

def source(path,cell):
 return ''.join(json.loads((NB/path).read_text())['cells'][cell]['source'])

def literals(code):
 """Read only literal list/dictionary assignments from the original plotting cell."""
 result={}
 for node in ast.parse(code).body:
  if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name):
   value=node.value
   if isinstance(value,ast.Call) and isinstance(value.func,ast.Attribute) and value.func.attr=='array': value=value.args[0]
   try: result[node.targets[0].id]=ast.literal_eval(value)
   except (ValueError,TypeError): pass
 return result

def table(name):
 d=pd.read_json(DATA/f'{name}.json',orient='table',precise_float=True)
 if 'position' in d:
  numeric=pd.to_numeric(d.position)
  assert np.equal(numeric,numeric.astype(int)).all()
  d['coordinate']=numeric.astype(int).astype(str)
 else: d['coordinate']=d.abs_position.astype(str)
 # Keep source coordinates and alleles, including contig identity for PbEc.
 cols=['coordinate']+[c for c in ['ref_seq','new_seq','aa_new_seq'] if c in d]
 d['mutation_id']=d[cols].fillna('').astype(str).agg('|'.join,axis=1)
 return d

def plain(label):
 label=re.sub(r'\$\\(?:it|mathrm)\{([^}]+)\}\$',r'\1',str(label))
 return label.replace('$\\Delta$','Δ').replace('small_indel','indel').replace('snp_intergenic','intergenic SNP').replace('mobile_element_insertion','IS insertion')

def matrix(d,panel,order,aliases,pop='Pop',value='last_freq'):
 d=d.copy()
 d['source_label']=d.label
 # U00096.3:4185113 is codon 1290 (ATG→AGG) in the original breseq call.
 # Keep the notebook label in the crosswalk, but use the correct residue in outputs.
 corrected=d.coordinate.eq('4185113') & d.label.eq('rpoB M129R')
 if corrected.any():
  assert d.loc[corrected,'aa_pos'].astype(float).eq(1290).all()
  d.loc[corrected,'label']='rpoB M1290R'
  order=['rpoB M1290R' if label=='rpoB M129R' else label for label in order]
  aliases={**aliases,'rpoB M1290R':'RpoB M1290R'}
 if d.duplicated([pop,'mutation_id']).any():
  raise ValueError(f'{panel}: repeated population/mutation identity')
 # Resolve all site differences, even if they occur in different cultures.
 counts=d.groupby('label').mutation_id.nunique()
 ambiguous=set(counts[counts>1].index)
 d['matrix_label']=[f'{label} [{coord}]' if label in ambiguous else label for label,coord in zip(d.label,d.coordinate)]
 numbers={label:{key:i+1 for i,key in enumerate(sorted(g.mutation_id.unique()))} for label,g in d.groupby('label') if label in ambiguous}
 d['display_label']=[plain(aliases.get(label,label))+(f' ({numbers[label][key]})' if label in ambiguous else '') for label,key in zip(d.label,d.mutation_id)]
 if d.duplicated([pop,'matrix_label']).any():
  # Same position, distinct alleles require the allele identity as well.
  dupl=set(d.loc[d.duplicated([pop,'matrix_label'],False),'matrix_label'])
  mask=d.matrix_label.isin(dupl)
  d.loc[mask,'matrix_label']=d.loc[mask,'matrix_label']+' '+d.loc[mask,'mutation_id']
  d.loc[mask,'display_label']=d.loc[mask,'display_label']+' '+d.loc[mask,'mutation_id']
 assert not d.duplicated([pop,'matrix_label']).any()
 m=d.pivot(index='matrix_label',columns=pop,values=value).fillna(0)
 rank={s:i for i,s in enumerate(order)}
 labels=d.drop_duplicates('matrix_label').set_index('matrix_label')
 row_order=sorted(m.index,key=lambda s:(rank.get(labels.loc[s,'label'],len(rank)),labels.loc[s,'label'],labels.loc[s,'mutation_id']))
 m=m.reindex(row_order)
 # Preserve the original mean-based values for comparison, without using them to plot.
 old=d.pivot_table(index='label',columns=pop,values=value,fill_value=0)
 collisions=d[d.duplicated([pop,'label'],False)]
 for (population,label),g in collisions.groupby([pop,'label']):
  CHANGES.append(dict(panel=panel,population=str(population),source_label=label,
   coordinates=g.coordinate.tolist(),frequencies=g[value].tolist(),historical_mean=float(old.loc[label,population])))
 for _,r in d.iterrows():
  CROSSWALK.append(dict(panel=panel,population=str(r[pop]),mutation_id=r.mutation_id,coordinate=r.coordinate,
   source_label=r.source_label,display_label=r.display_label,matrix_label=r.matrix_label,frequency=float(r[value])))
 m.to_csv(OUT/f'{panel}_matrix.csv',float_format='%.15g')
 labels[['label','display_label','coordinate','mutation_id']].to_csv(OUT/f'{panel}_labels.csv')
 return m,labels.loc[m.index,'display_label'].tolist()

def endpoint(panel):
 name,path,cell,cmap=ENDPOINTS[panel]; cfg=literals(source(path,cell));d=table(name)
 aliases=cfg.get('alt_labels',cfg.get('alternate_labels',cfg.get('alternate_labels_dict',{})))
 order=cfg.get('manual_label_order',cfg.get('custom_order',[]))
 if panel=='S4f':
  d=d[d.gene_product.isin(cfg['gene_product_names'])].copy()
  # Preserve gene-product selection while keeping each variant's identity.
  original_labels=d.label.copy()
  d['label']=d.gene_product.map(aliases).fillna(d.label)
  # The gene-product alias applied S342R to a separate atpD indel.
  atpd=d.gene_product.eq('ATP synthase beta chain (EC 3.6.3.14)') & original_labels.str.contains('small_indel')
  d.loc[atpd,'label']='atpD indel'
  aliases={};order=sorted(d.label.unique())
 elif panel!='S5a': d=d[d.label.isin(order)].copy()
 if not len(d): raise ValueError((panel,cfg.keys()))
 m,labels=matrix(d,panel,order,aliases)
 count=7 if panel=='S4f' else 6 if panel=='S3e' else 10
 m=m.reindex(columns=range(1,count+1),fill_value=0)
 if panel=='S4f':
  idx=np.argsort(-m[7].to_numpy(),kind='stable'); m=m.iloc[idx]; labels=[labels[i] for i in idx]
 m.to_csv(OUT/f'{panel}_matrix.csv',float_format='%.15g')
 return m,labels,cmap

def s6():
 name='S6d';d=table('MG_untreated');cfg=literals(source('Supplemental Figure 6/B - heatmap/05-Punt.ipynb',17));aliases=cfg['alt_labels']
 meta=pd.read_csv(ROOT/'data/genomics/reference/metadata_complete.csv')
 records=[]
 for _,r in d.iterrows():
  days=sorted(meta[(meta.Strain=='P')&(meta.Culture==r.Pop)].Day.unique())
  if len(days)!=len(r.freq):raise ValueError((r.Pop,days,r.freq))
  for day,f in zip(days,r.freq):records.append({**r.to_dict(),'sample':f'{r.Pop}_Day{day}','frequency':f})
 d=pd.DataFrame(records);d=d[d.label.isin(aliases)]
 m,labels=matrix(d,name,list(aliases),aliases,pop='sample',value='frequency')
 m=m.reindex(columns=sorted(m.columns,key=lambda s:(int(s.split('Day')[1]),int(s.split('_')[0]))))
 m.to_csv(OUT/f'{name}_matrix.csv',float_format='%.15g')
 return m,labels,'Oranges'

def s9():
 path='Supplemental Figure 9 - big heatmap/13-bigheatmap.ipynb';code=source(path,7)
 marker='    # Replace the index with numeric portions and transpose the heatmap'
 assert code.count(marker)==1
 code=code.split(marker)[0]+'\n    return combined_df, heatmap_data\n'
 ns={'np':np,'pd':pd,'plt':plt};exec(compile(code,path+':cell7-selection','exec'),ns)
 d,historical=ns['plot_combined_heatmap'](table('combined_lineages'),None)
 historical.to_csv(OUT/'S9_historical_matrix.csv',float_format='%.15g')
 order=historical.columns.tolist()
 # Notebook PLAC cell 14 distinguishes the insertion from the clade deletion.
 original=d.label.copy()
 d.loc[(d.position.astype(int)==3326602)&(d.gene_name=='ftsH'),'label']='ftsH small_ins'
 aliases={label:plain(label) for label in d.label.unique()}
 m,labels=matrix(d,'S9',order,aliases,pop='population')
 m=m.reindex(columns=historical.index)
 m.to_csv(OUT/'S9_matrix.csv',float_format='%.15g')
 return m,labels,'viridis'

def heatmap(m,labels,cmap,width,height,transpose=False,annot=False,bottom_override=None,top_override=None):
 """Dimensions and text sizes are physical points, independent of saved DPI."""
 plt.rcParams.update({'font.family':'Times New Roman','font.size':8.5,'pdf.fonttype':42,'svg.fonttype':'none'})
 arr=m.T if transpose else m
 y=list(map(str,m.columns)) if transpose else labels
 x=labels if transpose else list(map(str,m.columns))
 # Measure text in points rather than rely on raster-derived margins.
 temp=plt.figure(figsize=(1,1));renderer=temp.canvas.get_renderer()
 text_width=[]
 for s in y:
  t=temp.text(0,0,s,fontsize=8.5);text_width.append(t.get_window_extent(renderer).width*72/temp.dpi);t.remove()
 left=max(text_width,default=15)+(21 if transpose else 6)
 if transpose:
  xwidth=[]
  for item in x:
   t=temp.text(0,0,item,fontsize=8.5);xwidth.append(t.get_window_extent(renderer).width*72/temp.dpi);t.remove()
  left=max(left,max(w*.67-j*(width-left-5)/len(x) for j,w in enumerate(xwidth))+3)
 plt.close(temp)
 bottom=(max(len(str(s)) for s in x)*3.8*.74+18) if transpose else 30
 if bottom_override is not None:bottom=bottom_override
 top=25 if top_override is None else top_override;right=5
 if not transpose and any('Day' in str(s) for s in x):bottom=48
 if width-left-right<60 or height-bottom-top<50:raise ValueError(f'Panel too small: {width,height,left,bottom}')
 fig=plt.figure(figsize=(width/72,height/72))
 ax=fig.add_axes([left/width,bottom/height,(width-left-right)/width,(height-bottom-top)/height])
 colors=plt.get_cmap(cmap)(np.linspace(0,1,256));colors[0]=[1,1,1,1]
 mesh=ax.pcolormesh(arr.to_numpy(),cmap=ListedColormap(colors),vmin=0,vmax=1,edgecolors='#dddddd',linewidth=.25,rasterized=False)
 bar=fig.add_axes([(width-61)/width,(height-7)/height,56/width,3/height])
 cb=fig.colorbar(mesh,cax=bar,orientation='horizontal',ticks=[0,.5,1]);cb.ax.tick_params(labelsize=7.5,length=1.5,pad=1,width=.5);cb.outline.set_linewidth(.4);cb.solids.set_rasterized(False)
 fig.text(3/width,(height-3)/height,'Frequency',fontsize=8.5,va='top')
 ax.set_ylim(len(y),0);ax.set_xlim(0,len(x))
 ax.set_xticks(np.arange(len(x))+.5,x,rotation=48 if transpose else 90 if 'Day' in str(x[0]) else 0,ha='right' if transpose else 'center',rotation_mode='anchor')
 ax.set_yticks(np.arange(len(y))+.5,y)
 ax.tick_params(length=0,pad=3,labelsize=8.5)
 for s in ax.spines.values():s.set_linewidth(.65)
 if not transpose and 'Day' in str(x[0]):
  ax.set_xticklabels([v.split('_')[0] for v in x],rotation=0)
  for j in range(0,len(x),4):
   ax.text(j+2,-.065,'Day '+x[j].split('Day')[1],transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8.5)
   ax.axvline(j,color='black',lw=.7)
  ax.set_xlabel('Culture',fontsize=9.5,labelpad=20)
 elif not transpose:ax.set_xlabel('Culture',fontsize=9.5,labelpad=3)
 else:ax.set_ylabel('Culture',fontsize=9.5,labelpad=3)
 if annot and (width-left-right)/len(x)>=11:
  for i,row in enumerate(arr.to_numpy()):
   for j,v in enumerate(row):
    if v>0:ax.text(j+.5,i+.5,f'{v:.2f}',ha='center',va='center',fontsize=8,color='white' if v>.7 else 'black')
 return fig

def validate_bounds(fig,name):
 from matplotlib.text import Text
 fig.canvas.draw();renderer=fig.canvas.get_renderer();page=fig.bbox
 bad=[]
 for text in fig.findobj(Text):
  if not text.get_visible() or not text.get_text().strip():continue
  box=text.get_window_extent(renderer)
  if box.x0<-.5 or box.y0<-.5 or box.x1>page.width+.5 or box.y1>page.height+.5:bad.append(text.get_text())
 if bad:raise ValueError((name,'text outside page',bad))

def save(fig,name):
 FIG.mkdir(parents=True,exist_ok=True)
 validate_bounds(fig,name)
 fig.savefig(FIG/f'{name}.pdf',metadata={'CreationDate':None,'ModDate':None})
 fig.savefig(FIG/f'{name}.png',dpi=300)
 plt.close(fig)

def trajectories():
 d=table('PLAC');days=[0,3,16,19,31,34,47,55,63,82]
 code=source('Figure 4/02-PLAC.ipynb',19)
 groups={}
 for group in ['gyrA','rpoB','fusA','trkH','clade']:
  ns={'np':np,'sns':sns};exec(code.replace("group = 'fusA'",f'group = {group!r}'),ns)
  groups[group]=list(ns['select_mutations'])
 aliases=literals(source('Figure 4/02-PLAC.ipynb',29)).get('alt_labels',{})
 pal={group:dict(zip(values,sns.color_palette('muted',len(values)))) for group,values in groups.items()}
 def plot(ax,pop,group):
  sub=d[d.Pop==pop]; selected=groups[group]
  for _,r in sub.iterrows():
   assert len(r.freq)==len(days)
   is_selected=r.label in selected
   ax.plot(range(10),r.freq,color=pal[group].get(r.label,'.55'),lw=1.0 if is_selected else .45,alpha=1 if is_selected else .15,
           marker='o' if is_selected else None,ms=2,label=plain(aliases.get(r.label,r.label)) if is_selected else None)
  ax.set_xlim(-.15,9.15);ax.set_ylim(-.04,1.04);ax.set_yticks([0,.5,1]);ax.set_xticks(range(10),days)
  ax.tick_params(labelsize=8.5,width=.65,length=2)
  for s in ax.spines.values():s.set_linewidth(.65)
  ax.set_xlabel('Day',fontsize=9.5);ax.set_ylabel('Frequency',fontsize=9.5)
 def legend(ax,group,pop=None):
  from matplotlib.lines import Line2D
  return [Line2D([0],[0],color=pal[group][l],lw=1,label=plain(aliases.get(l,l))) for l in groups[group] if pop is None or l in set(d.loc[d.Pop==pop,'label'])]
 # Main panel: culture 4, four selected mutation groups, consistent with caption.
 fig,axes=plt.subplots(2,2,figsize=(240/72,220/72))
 fig.subplots_adjust(left=.135,right=.985,bottom=.37,top=.85,hspace=1.5,wspace=.48)
 fig.text(.135,.985,r'MG$^{\mathrm{LEV,AMI,CEF}}$-4 (persistent)',fontsize=9,va='top')
 for i,(ax,group) in enumerate(zip(axes.flat,['gyrA','fusA','trkH','clade'])):
  plot(ax,4,group);ax.set_title(['i  gyrA','ii  fusA','iii  trkH','iv  Convergent alleles'][i],loc='left',fontsize=9)
  ax.set_xticks([0,2,4,6,9],[0,16,31,47,82])
  ax.set_xlabel('Day' if i>=2 else '',fontsize=9,labelpad=2)
  ax.set_ylabel('Frequency' if i%2==0 else '',fontsize=9,labelpad=2)
  if i%2:ax.tick_params(labelleft=False)
  if i<2:ax.legend(handles=legend(ax,group,4),loc='upper left',bbox_to_anchor=(-.06,-.27),ncol=1,fontsize=8.5,frameon=False,handlelength=1.0,columnspacing=.5,labelspacing=.05,borderpad=0)
 fig.legend(handles=legend(None,'clade',4),loc='lower left',bbox_to_anchor=(.07,.005),ncol=2,fontsize=8.5,frameon=False,handlelength=1.2,columnspacing=.65,labelspacing=.1)
 save(fig,'Fig4a')
 # Supplementary trajectories split by culture into two readable, matching grids.
 with PdfPages(FIG/'Supp_Figure_8.pdf',metadata={'CreationDate':None}) as pdf:
  for page,pops in enumerate([range(1,6),range(6,11)],1):
   fig,axes=plt.subplots(5,5,figsize=(10.6,7.6));fig.subplots_adjust(left=.055,right=.99,bottom=.29,top=.92,hspace=.58,wspace=.24)
   for row,pop in enumerate(pops):
    for col,group in enumerate(groups):
     ax=axes[row,col];plot(ax,pop,group)
     if row==0:ax.set_title(f'{chr(97+col)}  '+('Convergent alleles' if group=='clade' else group),loc='left',fontsize=10)
     if col>0:ax.set_ylabel('');ax.set_yticklabels([])
     if row<4:ax.set_xlabel('');ax.set_xticklabels([])
     if col==0:ax.text(-.30,.5,f'Culture {pop}',rotation=90,transform=ax.transAxes,va='center',ha='center',fontsize=9)
     ax.tick_params(axis='x',labelsize=7.5)
   for col,group in enumerate(groups):
    axes[-1,col].legend(handles=legend(axes[-1,col],group),loc='upper left',bbox_to_anchor=(-.08,-.65),fontsize=7.5,frameon=False,handlelength=1.2,labelspacing=.1)
   fig.suptitle(f'Supplementary Figure 8 · cultures {pops.start}–{pops.stop-1}',fontsize=11,x=.055,ha='left')
   validate_bounds(fig,f'S8 page {page}')
   pdf.savefig(fig);fig.savefig(FIG/f'Supp_Figure_8_page{page}.png',dpi=300);plt.close(fig)
 # Presence from first to last observed nonzero sampling day, as in the notebook.
 clade=groups['clade'];fig,axes=plt.subplots(2,5,figsize=(322/72,220/72),sharey=True)
 fig.subplots_adjust(left=.30,right=.99,bottom=.13,top=.93,wspace=.18,hspace=.32)
 emergence=[]
 for pop,ax in enumerate(axes.flat,1):
  for j,label in enumerate(clade):
   rows=d[(d.Pop==pop)&(d.label==label)]
   assert len(rows)<=1,(pop,label)
   if len(rows):
    nonzero=np.where(np.array(rows.iloc[0].freq)>0)[0]
    if len(nonzero):
     start,end=int(nonzero[0]),int(nonzero[-1]);ax.plot([start,end],[j,j],lw=4,color=pal['clade'][label],solid_capstyle='butt')
     ax.scatter([start,end],[j,j],s=6,color=pal['clade'][label]);emergence.append(dict(culture=pop,label=label,first_observed_day=days[start],last_observed_day=days[end]))
  ax.set_title(f'Culture {pop}',fontsize=8.5);ax.set_xlim(-.5,9.5);ax.set_ylim(8.5,-.5);ax.set_xticks([0,4,9],[0,31,82]);ax.tick_params(labelsize=8.5,length=2,width=.65)
  ax.set_yticks(range(9),[plain(aliases.get(l,l)) for l in clade]);ax.set_xlabel('Day',fontsize=8.5)
  if pop<=5:ax.set_xlabel('');ax.set_xticklabels([])
  for s in ax.spines.values():s.set_linewidth(.65)
 save(fig,'Fig4b');pd.DataFrame(emergence).to_csv(OUT/'Fig4b_observed_spans.csv',index=False)
 rows=[]
 for _,r in d.iterrows():
  for day,f in zip(days,r.freq):rows.append(dict(culture=r.Pop,mutation_id=r.mutation_id,coordinate=r.coordinate,source_label=r.label,day=day,frequency=f))
 pd.DataFrame(rows).to_csv(OUT/'PLAC_trajectories.csv',index=False,float_format='%.15g')

def main():
 OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.family':'Times New Roman','font.size':8.5,'pdf.fonttype':42,'svg.fonttype':'none'})
 for panel in ENDPOINTS:
  m,labels,cmap=endpoint(panel)
  fig=heatmap(m,labels,cmap,579 if panel=='Fig4c' else 464.4,125 if panel=='Fig4c' else max(180,len(m)*10+45),annot=len(m)<16,bottom_override=18 if panel=='Fig4c' else None,top_override=20 if panel=='Fig4c' else None)
  if panel=='Fig4c':
   from matplotlib.lines import Line2D
   for culture,label in enumerate(fig.axes[0].get_xticklabels(),1):
    label.set_color('teal' if culture in [2,5] else 'maroon' if culture in [3,8] else '#bb40bb');label.set_weight('bold')
   fig.axes[0].set_xlabel('')
   fig.text(.105,.015,'Culture',fontsize=9,ha='right',va='bottom')
   fig.legend(handles=[Line2D([],[],marker='s',ls='',ms=4,color=color,label=label) for color,label in [('teal','Tolerant'),('#bb40bb','Persistent'),('maroon','Intermediate')]],loc='upper center',bbox_to_anchor=(.55,1.02),ncol=3,fontsize=8.5,frameon=False,handlelength=.6,columnspacing=1)
  save(fig,panel)
 m,labels,cmap=s6();save(heatmap(m,labels,cmap,464.4,len(m)*10+70),'S6d')
 m,labels,cmap=s9()
 with PdfPages(FIG/'Supp_Figure_9.pdf',metadata={'CreationDate':None}) as pdf:
  for page,start in enumerate(range(0,len(m),43),1):
   sub=m.iloc[start:start+43];fig=heatmap(sub,labels[start:start+43],cmap,720,len(sub)*10+80)
   ax=fig.axes[0];ax.set_xticklabels([s.split('_')[-1].lstrip('0') for s in sub.columns],rotation=0)
   groups=[]
   for j,c in enumerate(sub.columns):
    prefix=c.split('_')[0]
    if not groups or groups[-1][0]!=prefix:groups.append([prefix,j,j])
    else:groups[-1][2]=j
   for group,a,b in groups:
    ax.text((a+b+1)/2,1.01,group,transform=ax.get_xaxis_transform(),ha='center',va='bottom',fontsize=9)
    ax.axvline(a,color='black',lw=.7)
   ax.set_xlabel('Culture; allele frequency 0 (white) to 1 (yellow)',fontsize=9.5)
   fig.text(.01,.018,f'Supplementary Figure 9 · {page}',va='bottom',fontsize=9)
   validate_bounds(fig,f'S9 page {page}')
   pdf.savefig(fig);fig.savefig(FIG/f'Supp_Figure_9_page{page}.png',dpi=300);plt.close(fig)
 trajectories()
 pd.DataFrame(CROSSWALK).to_csv(OUT/'panel_source_crosswalk.csv',index=False,float_format='%.15g')
 (OUT/'aggregation_corrections.json').write_text(json.dumps(CHANGES,indent=2)+'\n')
 (OUT/'manifest.json').write_text(json.dumps({'source':'Regenerated genomic tables; original notebook selections and aliases',
   'mutation_identity':'coordinate plus source alleles; distinct variants are not averaged',
   'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.glob('*.json'))},
   'plot_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
 print(f'Rendered genomic plots; {len(CROSSWALK)} matrix entries; {len(CHANGES)} historical duplicate groups.')

if __name__=='__main__':main()
