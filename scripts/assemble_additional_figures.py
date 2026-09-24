"""Compose numerical genomic figures and retain the supplied schematic explicitly."""
from pathlib import Path
import json,hashlib
import fitz
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'working/figures-assembled';PANELS=ROOT/'working/analysis/stats-rework/out/genomic-panels'
RECORD=[]
def place(page,path,rect,label=None):
 with fitz.open(path) as src:
  page.show_pdf_page(fitz.Rect(rect),src,keep_proportion=True)
  scale=min((rect[2]-rect[0])/src[0].rect.width,(rect[3]-rect[1])/src[0].rect.height)
  RECORD.append({'source':path.relative_to(ROOT).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'scale':scale,'label':label})
 if label:page.insert_text((rect[0],rect[1]-3),label,fontsize=12,fontname='hebo')
def finish(doc,name):
 doc.set_metadata({'title':name.replace('_',' '),'subject':'Source data and code: https://github.com/msadikyildiz/Lyon_2026'})
 doc.save(OUT/f'{name}.pdf',garbage=4,deflate=True)
 for i,page in enumerate(doc):
  suffix='' if len(doc)==1 else f'_page{i+1}'
  page.get_pixmap(dpi=300).save(OUT/f'{name}{suffix}.png');page.get_pixmap(dpi=100).save(OUT/f'{name}{suffix}_preview.png')
 doc.close()
def main():
 # Numerical panels use their final physical sizes; retain the vector schematic.
 doc=fitz.open();p=doc.new_page(width=595.28,height=782.36)
 p.insert_text((8,16),'Figure 4',fontsize=12,fontname='hebo')
 place(p,PANELS/'Fig4a.pdf',[8,36,248,256],'a')
 place(p,PANELS/'Fig4b.pdf',[265,36,587,256],'b')
 place(p,PANELS/'Fig4c.pdf',[8,274,587,399],'c')
 asset=json.loads((ROOT/'data/figure-assets/Figure4_schematic.json').read_text());source=ROOT/asset['source']
 assert hashlib.sha256(source.read_bytes()).hexdigest()==asset['sha256']
 clip=fitz.Rect(asset['crop'])
 width=464.4;h=width*clip.height/clip.width;x=(p.rect.width-width)/2;y=415
 p.insert_text((8,y-3),'d',fontsize=12,fontname='hebo')
 with fitz.open(source) as schematic:
  p.show_pdf_page(fitz.Rect(x,y,x+width,y+h),schematic,pno=asset['page'],clip=clip,keep_proportion=True)
 p.insert_text((8,777),asset['credit'],fontsize=8,fontname='helv')
 RECORD.append({**asset,'figure':'Figure 4','page':1,'rect':[x,y,x+width,y+h],'retained_schematic':True});finish(doc,'Figure_4')
 # Complete endpoint heatmaps retain a fixed physical label size across three pages.
 doc=fitz.open()
 for letter in 'abc':
  with fitz.open(PANELS/f'S5{letter}.pdf') as src:
   page=doc.new_page(width=480.4,height=src[0].rect.height+35)
   page.insert_text((8,16),'Supplementary Figure 5'+(' (continued)' if letter!='a' else ''),fontsize=12,fontname='hebo')
   place(page,PANELS/f'S5{letter}.pdf',[8,35,472.4,src[0].rect.height+35],letter)
 finish(doc,'Supplementary_Figure_5')
 for number in [8,9]:
  with fitz.open(PANELS/f'Supp_Figure_{number}.pdf') as src:
   doc=fitz.open();doc.insert_pdf(src)
  RECORD.append({'source':(PANELS/f'Supp_Figure_{number}.pdf').relative_to(ROOT).as_posix(),'figure':f'Supplementary Figure {number}','pages':len(doc)})
  finish(doc,f'Supplementary_Figure_{number}')
 (OUT/'additional_assembly_manifest.json').write_text(json.dumps(RECORD,indent=2)+'\n')
 files=sorted([p for p in OUT.glob('*.pdf') if p.stem.startswith(('Figure_','Supplementary_Figure_'))],key=lambda p:(p.stem.startswith('Supplementary'),int(p.stem.split('_')[-1])))
 with fitz.open() as combined:
  toc=[]
  for file in files:
   toc.append([1,file.stem.replace('_',' '),len(combined)+1])
   with fitz.open(file) as src:combined.insert_pdf(src)
  combined.set_toc(toc);combined.save(OUT/'All_figures.pdf',garbage=4,deflate=True)
 print(f'Combined {len(files)} figures; multi-page panels retain readable physical sizes.')
if __name__=='__main__':main()
