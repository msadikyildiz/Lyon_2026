"""Assemble Figure 1 and Supplementary Figure 1 with supplied schematic artwork."""
from pathlib import Path
import os,sys,json,hashlib
os.environ.setdefault('MPLCONFIGDIR',str(Path.home()/'.cache'/'matplotlib'));os.environ['MPLBACKEND']='Agg'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fitz
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'working/figures-assembled';DATA=ROOT/'data/figure-assets'
def main():
 manifest=json.loads((DATA/'schematic_sources.json').read_text())
 for name,item in manifest.items():assert hashlib.sha256((DATA/name).read_bytes()).hexdigest()==item['sha256']
 plt.rcParams.update({'font.family':'Times New Roman','font.size':9.5,'mathtext.fontset':'stix','pdf.fonttype':42})
 fig=plt.figure(figsize=(6.45,4.25));fig.text(.02,.965,'Figure 1',weight='bold',fontsize=12)
 ax=fig.add_axes([.035,.52,.94,.37]);ax.imshow(plt.imread(DATA/'Figure1_schematics.png'));ax.axis('off')
 for i,label in enumerate('abc'):fig.text(.02+i*.333,.905,label,fontsize=12,weight='bold')
 points=[]
 for i,(sheet,color,drug) in enumerate([('PA','royalblue','Amikacin'),('PL','firebrick','Levofloxacin'),('PC','forestgreen','Cefepime')]):
  raw=pd.read_excel(ROOT/'working/figures/Figure 1/D-F/SurvivalData.xlsx',sheet_name=sheet)
  d=raw[(raw.Strain=='MG1655')&raw.Culture.isin(range(1,11))].copy();d['day']=np.floor(d.Day).astype(int);d['phase']=np.where(d.Day==d.day,'before','after')
  assert not d.duplicated(['Culture','day','phase']).any()
  d=d.pivot(index=['Culture','day'],columns='phase',values='CFU').reset_index();assert d[['before','after']].notna().all().all()
  d['survival_percent']=100*d.after/d.before;d['antibiotic']=drug;points.append(d)
  ax=fig.add_axes([.085+i*.331,.125,.245,.33]);fig.text(.02+i*.333,.48,'def'[i],fontsize=12,weight='bold')
  for _,g in d.groupby('Culture'):ax.plot(g.day,g.survival_percent,color='.7',lw=.7)
  mean=d.groupby('day').survival_percent.mean();ax.plot(mean.index,mean.values,color=color,lw=1.1)
  ax.set_yscale('log');ax.set_ylim(.0005,3000);last=int(d.day.max());ax.set_xlim(1,last);ax.set_xticks(sorted(set([1,last]+[t for t in [5,10,15] if t<=last-2])));ax.set_yticks([.001,.1,10,1000]);ax.tick_params(labelsize=8.5,length=2.5,width=.65,pad=2)
  ax.set_title(drug,fontsize=9.5,pad=5);ax.set_xlabel('Day',fontsize=9.5);ax.set_ylabel('Survival (%)',fontsize=9.5,labelpad=2)
  ax.spines[['top','right']].set_visible(False)
  for s in ax.spines.values():s.set_linewidth(.65)
 pd.concat(points).to_csv(ROOT/'working/analysis/stats-rework/out/fig1_all_survival_records.csv',index=False)
 for ext,dpi in [('pdf',300),('png',300)]:fig.savefig(OUT/f'Figure_1.{ext}',dpi=dpi)
 fig.savefig(OUT/'Figure_1_preview.png',dpi=100);plt.close(fig)
 with fitz.open(DATA/'Supplementary_Figure_1_schematic.pdf') as src:
  doc=fitz.open();page=doc.new_page(width=464.4,height=src[0].rect.height+30);page.insert_text((8,16),'Supplementary Figure 1',fontsize=12,fontname='hebo');page.show_pdf_page(fitz.Rect(8,30,456.4,page.rect.height),src)
  doc.save(OUT/'Supplementary_Figure_1.pdf',garbage=4,deflate=True);page.get_pixmap(dpi=300).save(OUT/'Supplementary_Figure_1.png');page.get_pixmap(dpi=100).save(OUT/'Supplementary_Figure_1_preview.png')
 print('Rendered Figure 1 from numerical trajectories and supplied schematics; retained S1 illustration.')
if __name__=='__main__':main()
