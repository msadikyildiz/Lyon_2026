"""Render assembly-specific panels with typography specified at manuscript print size.

Numerical helpers and cached observations are shared with plot_final.py. MDK panels retain confirmed zero-count observations as censored records.
No dose-response statistics are refitted.
"""
import hashlib
import fitz
import json
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_final as original
from fig5a_doubling import XLSX as GROWTH_SOURCE
from plot_final import HERE, LAB, COL, PANELS, DRUGS, spread, geo_ci, nice_limits, brackets_for

ROOT = HERE.parents[2]
OUT = HERE / 'out/print-panels'
TICK = 8.5
AXIS = 9.5
ANNOTATION = 8.0
SPINE = .75
LINE = 1.05
ERROR = .65


def canvas(width, height, bottom=36, top=6, right=5, left=47):
    fig = plt.figure(figsize=(width/72, height/72), dpi=144)
    ax = fig.add_axes([left/width, bottom/height, (width-left-right)/width,
                       (height-bottom-top)/height])
    return fig, ax


def style(ax, log=True):
    if log:
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10, numticks=20))
        ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext(base=10))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=range(2, 10), numticks=100))
        ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    for spine in ax.spines.values():
        spine.set_linewidth(SPINE)
    ax.tick_params(axis='both', which='major', length=3, width=SPINE, labelsize=TICK, pad=2)
    ax.tick_params(axis='both', which='minor', length=1.7, width=.5)
    ax.xaxis.label.set_size(AXIS)
    ax.yaxis.label.set_size(AXIS)
    ax.yaxis.labelpad = 3
    ax.xaxis.labelpad = 3


def fit_panel(asset, width, height, final):
    stem, drug, value, _ = asset.split('_')
    spec = next(s for s in PANELS if s[3] == stem and s[4] == value)
    panel, ds, order, _, _, interval, use_brackets = spec
    gf = original.comparison_cohort(ds, original.split_strain(pd.read_pickle(HERE/'cache'/f'{ds}.pkl')))
    sub = gf[gf.Antibiotic == drug]
    present = [g for g in order if g in set(sub['group'])]
    br = brackets_for(final, panel, drug, present) if use_brackets else []
    rotation = 65 if len(present) > 4 else (0 if len(present) == 2 else 50)
    bottom = 48 if len(present) > 4 else (22 if len(present) == 2 else 56 if 'PLAC' in present else 42)
    fig, ax = canvas(width, height, bottom=bottom)
    summaries, points, tops = [], [], []
    for i, group in enumerate(present):
        values = sub[sub['group'] == group][value].dropna().values
        gm, lo, hi = geo_ci(values)
        summaries.append(dict(group=group,mean=gm,low=lo if interval else None,high=hi if interval else None))
        points.extend(dict(group=group,value=float(v)) for v in values)
        if interval and len(values) >= 10:
            ax.boxplot([values], positions=[i], widths=.48, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=COL[group], alpha=.12, linewidth=.5),
                       medianprops=dict(color='grey', linewidth=.5), whiskerprops=dict(color='grey',linewidth=.5),
                       capprops=dict(color='grey',linewidth=.5), zorder=3, manage_ticks=False)
        ax.scatter(spread(len(values),i),values,s=10,facecolor=COL[group],edgecolor='black',linewidth=.35,zorder=10)
        ax.hlines(gm,i-.42,i+.42,color='black',linewidth=LINE,zorder=7)
        if interval and np.isfinite(lo) and lo > 0:
            ax.errorbar(i,gm,yerr=[[gm-lo],[hi-gm]],fmt='none',ecolor='black',elinewidth=ERROR,capsize=2,capthick=ERROR,zorder=7)
            tops.append(hi)
        tops.append(values.max())
    lows = [r['low'] for r in summaries if r['low'] is not None and np.isfinite(r['low'])]
    highs = [r['high'] for r in summaries if r['high'] is not None and np.isfinite(r['high'])]
    ylim = nice_limits([p['value'] for p in points],len(br),max(highs) if highs else None,min(lows) if lows else None)
    style(ax)
    # Reserve actual point-sized bracket rows instead of scaling tiny annotations with the panel.
    axis_height = ax.get_position().height * height
    first_gap, step = 4/axis_height, 11/axis_height
    loglo, loghi = np.log10(ylim)
    if br:
        reserve = first_gap + len(br)*step + 2/axis_height
        if reserve >= .82:
            raise ValueError(f'{asset}: increase panel height for readable brackets')
        loghi = max(loghi, loglo+(np.log10(max(tops))-loglo)/(1-reserve))
    ax.set_ylim(10**loglo,10**loghi)
    ax.set_xlim(-.7,len(present)-.3)
    ax.set_xticks(range(len(present)),[LAB.get(g,g) for g in present],rotation=rotation,
                  ha='right' if rotation else 'center',rotation_mode='anchor')
    ax.set_ylabel(f"{drug}\n{'IC$_{50}$' if value=='IC50' else 'MIC'} (μg/mL)")
    for tick in ax.get_xticklabels():tick.set_fontsize(TICK)
    start=(np.log10(max(tops))-loglo)/(loghi-loglo)+first_gap
    labels=[]
    for j,r in enumerate(br):
        x1,x2=sorted((present.index(r['Group1']),present.index(r['Group2'])))
        y=start+j*step
        ax.plot([x1,x1,x2,x2],[y,y+1.5/axis_height,y+1.5/axis_height,y],
                color='black',lw=.6,transform=ax.get_xaxis_transform(),clip_on=False)
        text=original.p_text(r['p_holm']);labels.append(text)
        ax.text((x1+x2)/2,y+2/axis_height,text,fontsize=ANNOTATION,ha='center',va='bottom',transform=ax.get_xaxis_transform())
    return fig,dict(kind='susceptibility',asset=asset,points=points,summaries=summaries,brackets=labels,
                    ylim=list(ax.get_ylim()),tick_pt=TICK,axis_pt=AXIS,spine_pt=SPINE)


def growth_panel(width,height):
    fig,ax=canvas(width,height,bottom=45,left=37)
    g=original.load_doubling_times();order=['MG','gatA','glvC','hipA','selB','rpoZ','ftsH','fimE']
    keys=['MG','gata','glvc','hipa','selb','rpoz','ftsh','fime'];payload=[]
    for i,name in enumerate(order):
        v=g[g.Culture==name].doubling_time_min.dropna().values
        ax.scatter(spread(len(v),i),v,s=10,c='grey',edgecolor='black',linewidth=.35,zorder=10)
        ax.hlines(v.mean(),i-.42,i+.42,color='black',lw=LINE,zorder=7)
        ax.errorbar(i,v.mean(),yerr=v.std(ddof=1),fmt='none',ecolor='black',elinewidth=ERROR,capsize=2,capthick=ERROR)
        payload.append(dict(group=name,values=v.tolist(),mean=v.mean(),sd=v.std(ddof=1)))
    style(ax,False);ax.set_ylim(0,25);ax.set_xlim(-.7,7.7)
    ax.set_xticks(range(8),[LAB[k] for k in keys],rotation=55,ha='right',rotation_mode='anchor')
    ax.set_ylabel('Doubling time (min)')
    return fig,dict(kind='growth',groups=payload,tick_pt=TICK,axis_pt=AXIS,spine_pt=SPINE)


def displaced_end_labels(ax, entries, height):
    """Minimize label displacement with measured text gaps and monotone leaders."""
    fig=ax.figure
    ordered=sorted(entries,key=lambda e:e['y'])
    loglo,loghi=np.log10(ax.get_ylim())
    axis_height=ax.get_position().height*height
    axis_width=ax.get_position().width*fig.get_figwidth()*72
    desired=np.array([(np.log10(e['y'])-loglo)/(loghi-loglo)*axis_height for e in ordered])
    texts=[ax.text(1+14/axis_width,y/axis_height,e['label'],transform=ax.transAxes,
                   fontsize=8,ha='left',va='center',clip_on=False,
                   color=tuple(.7*c for c in mcolors.to_rgb(e['color'])))
           for e,y in zip(ordered,desired)]
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    sizes=np.array([t.get_window_extent(renderer).height*72/fig.dpi for t in texts])
    gap=1.6
    offsets=np.r_[0,np.cumsum((sizes[:-1]+sizes[1:])/2+gap)]
    lower=1+sizes[0]/2
    upper=axis_height-1-sizes[-1]/2-offsets[-1]
    if upper<lower:
        raise ValueError(f'Endpoint labels need {axis_height+lower-upper:.2f} pt; axes have {axis_height:.2f} pt')
    # Pool-adjacent-violators solves the least-squares ordered-position problem.
    blocks=[]
    for i,target in enumerate(desired-offsets):
        blocks.append([i,i+1,float(target),1])
        while len(blocks)>1 and blocks[-2][2]/blocks[-2][3]>blocks[-1][2]/blocks[-1][3]:
            right=blocks.pop();left=blocks.pop()
            blocks.append([left[0],right[1],left[2]+right[2],left[3]+right[3]])
    fitted=np.empty(len(texts))
    for start,stop,total,count in blocks:fitted[start:stop]=np.clip(total/count,lower,upper)
    positions=fitted+offsets
    endpoint=(7-ax.get_xlim()[0])/(ax.get_xlim()[1]-ax.get_xlim()[0])
    result=[]
    for e,text,target,pos in zip(ordered,texts,desired,positions):
        text.set_position((1+14/axis_width,pos/axis_height))
        # All diagonal segments stay in the right-hand gutter, away from data.
        ax.plot([endpoint+2/axis_width,1+2/axis_width,1+11/axis_width,1+12/axis_width],
                [target/axis_height,target/axis_height,pos/axis_height,pos/axis_height],
                transform=ax.transAxes,clip_on=False,color='.5',lw=.45,zorder=2)
        result.append(dict(label=e['label'],endpoint_y=e['y'],label_y_axes=pos/axis_height,
                           displacement_pt=float(pos-target)))
    fig.canvas.draw()
    boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in texts]
    actual_gaps=[(b.y0-a.y1)*72/fig.dpi for a,b in zip(boxes,boxes[1:])]
    assert all(g>=gap-.05 for g in actual_gaps),actual_gaps
    assert boxes[0].y0>=ax.bbox.y0 and boxes[-1].y1<=ax.bbox.y1
    return dict(labels=result,minimum_text_gap_pt=min(actual_gaps),leader_diagonals_outside_axes=True)


def legacy_mdk(asset,width,height):
    mapping={
        'original_image8':('Fig2e_MDK','PC','MG','CEF',10,{2,4,5,6,10},{7,8}),
        'original_image17':('Fig3e_MDK','PLAC','MG','LEV,AMI,CEF',10,{2,5},{3,8}),
        'original_image41':('SuppFig4e_MDK','ATECc','Pb','CEF',6,{1,5,6},set()),
    }
    stem,prefix,base,exponent,n,teal,maroon=mapping[asset]
    parent='ATEC1' if prefix=='ATECc' else 'P1'
    summ=pd.read_csv(HERE/'out'/f'{stem}_summary.csv')
    # Summaries retain confirmed zero plates and validate their censoring intervals.
    order=[parent]+[f'{prefix}{i}' for i in range(1,n+1)]
    colors={parent:'black',**{f'{prefix}{i}':('teal' if i in teal else 'maroon' if i in maroon else 'violet') for i in range(1,n+1)}}
    right=89 if prefix=='PLAC' else 57
    fig,ax=canvas(width,height,bottom=29,top=15.5 if prefix=='PC' else 6,right=right,left=35)
    payload=[];entries=[]
    for name in order:
        s=summ[summ.Name==name].sort_values('Time')
        med=s['median'].to_numpy();mad=s.mad_scaled.to_numpy()
        ax.errorbar(s.Time,med,yerr=mad,fmt='-o',color=colors[name],lw=LINE,markersize=2.8,
                    elinewidth=ERROR,capsize=1.5,capthick=ERROR,alpha=.8,zorder=3)
        label=base+'-1' if name==parent else base+rf'$^{{\mathrm{{{exponent}}}}}$-'+name[len(prefix):]
        entries.append(dict(label=label,y=med[-1],color=colors[name]))
        payload.extend(dict(Name=name,Time=float(t),median=float(m),mad=float(a),basis='confirmed records; censored median and normal-scaled MAD') for t,m,a in zip(s.Time,med,mad))
    style(ax);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    ax.set_ylim(1e-8,2);ax.set_xlim(-.12,7.15);ax.set_xticks([0,1,2,3,5,7]);ax.set_xlabel('Time (h)');ax.set_ylabel('Survivor fraction')
    if prefix=='PLAC':
        records=pd.read_csv(HERE/'out'/f'{stem}_records.csv')
        censored=records[records.below_detection & records.included_in_original_panel]
        ax.scatter(censored.Time,censored.detection_limit,marker='v',s=24,facecolor='white',edgecolor='black',linewidth=.75,zorder=10)
    end_labels=displaced_end_labels(ax,entries,height)
    if prefix=='PC':
        handles=[mlines.Line2D([],[],color=c,lw=LINE,label=t) for c,t in [('teal','Tolerant'),('violet','Persistent'),('maroon','Intermediate')]]
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,1),frameon=False,fontsize=7.5,ncol=3,
                   handlelength=1.4,handletextpad=.4,columnspacing=1,labelspacing=.1,borderaxespad=0)
    return fig,dict(kind='legacy_mdk',source_summary=stem+'_summary.csv',records=payload,end_labels=end_labels,
                    tick_pt=TICK,axis_pt=AXIS,spine_pt=SPINE,line_pt=LINE)


def modern_mdk(asset,width,height):
    if asset=='SuppFig11_MDK_double_mutants':
        stem=asset;rep=pd.read_csv(HERE/'out/supp11_recomputed.csv');summ=pd.read_csv(HERE/'out/supp11_summary.csv')
        order=['hipA','hipA_selB','hipA_rpoZ','hipA_ftsH','wt']
        col=dict(zip(order,['violet','lightsalmon','blueviolet','steelblue','black']))
        lab={'wt':'MG','hipA':LAB['hipa'],**{g:'MG'+rf'$^{{\mathrm{{{g.replace("_","+")}}}}}$' for g in order[1:-1]}}
    elif asset=='Fig5e_MDK':
        stem=asset;order=['gata','glvc','hipa','selb','rpoz','ftsh','fimE','wt']
        col=dict(zip(order,['teal','blue','violet','maroon','orange','purple','red','black']))
        lab={**{g:LAB[g] for g in order[:6]},'fimE':LAB['fime'],'wt':'MG'}
    else:
        stem=asset;order=['hipa','plac','prs','pc','wt'];col=dict(zip(order,['violet','purple','teal','blue','black']))
        lab={'hipa':LAB['hipa'],'plac':LAB['PLAC']+'-7','prs':r'MG$^{\mathrm{prs}}$','pc':LAB['PC']+'-5','wt':'MG'}
    if asset!='SuppFig11_MDK_double_mutants':
        rep=pd.read_csv(HERE/'out'/f'{stem}_records.csv');summ=pd.read_csv(HERE/'out'/f'{stem}_summary.csv')
    fig,ax=canvas(width,height,bottom=29,top=62 if len(order)>5 else 52,left=36,right=5)
    offsets=dict(zip(order,np.linspace(-.25,.25,len(order))));payload=[]
    for name in order:
        s=summ[summ.Name==name].sort_values('Time');med=s['median'].values;mad=s.mad_scaled.values;t=s.Time.values+offsets[name]
        ax.errorbar(t,med,yerr=[np.where(med-mad>0,mad,0),mad],color=col[name],lw=LINE,marker='o',markersize=2.7,
                    capsize=1.5,capthick=ERROR,elinewidth=ERROR,alpha=.8,label=lab[name],zorder=4)
        r=rep[rep.Name==name];x=np.empty(len(r))
        for time in r.Time.unique():
            mask=(r.Time==time).values;x[mask]=spread(mask.sum(),time+offsets[name],s=.035,gap=.012)
        censored=r.below_detection.to_numpy() if 'below_detection' in r else np.zeros(len(r),bool)
        values=r.fraction_plot.to_numpy() if 'fraction_plot' in r else r.fraction_observed.to_numpy()
        ax.scatter(x[~censored],values[~censored],s=8,facecolor=col[name],edgecolor='black',linewidth=.3,zorder=8)
        if censored.any():ax.scatter(x[censored],values[censored],s=24,facecolor='white',edgecolor=col[name],linewidth=.75,marker='v',zorder=10)
        if 'mad_scaled_max' in s:
            for x0,mid,lo,hi in zip(t,med,s.mad_scaled_min,s.mad_scaled_max):
                if hi>lo:
                    ax.plot([x0,x0],[mid+lo,mid+hi],ls='--',color=col[name],lw=ERROR)
                    if mid-hi>0:ax.plot([x0,x0],[mid-hi,mid-lo],ls='--',color=col[name],lw=ERROR)
        payload.extend(dict(Name=name,Time=float(time),median=float(m),mad=float(a)) for time,m,a in zip(s.Time,med,mad))
    style(ax);ax.set_ylim(1e-9,1.5);ax.set_xlim(-.5,7.5);ax.set_xticks([0,1,2,3,5,7]);ax.set_xlabel('Time (h)');ax.set_ylabel('Survivor fraction')
    axis_height = ax.get_position().height * height
    ax.legend(fontsize=8,frameon=False,loc='lower center',bbox_to_anchor=(.5,1+4/axis_height),ncol=2,
              handlelength=1.4,columnspacing=1,handletextpad=.4,labelspacing=.4,borderaxespad=0)
    return fig,dict(kind='mdk',asset=asset,records=payload,replicate_points=len(rep),tick_pt=TICK,axis_pt=AXIS,spine_pt=SPINE,line_pt=LINE)


def other_numeric(asset,width,height):
    if asset=='original_image5':
        raw=pd.read_pickle(HERE/'cache/fig2_paplpc__df_analysis.pkl')
        data=raw[(raw.Strain=='PA5')&(raw.Antibiotic=='Levofloxacin')]
        fit=pd.read_pickle(HERE/'cache/fig2_paplpc.pkl')
        fit=fit[(fit.Strain=='PA5')&(fit.Antibiotic=='Levofloxacin')].iloc[0]
        fig,ax=canvas(width,height,left=36,bottom=50)
        # Stored observations and fitted curve are exported in Source Data.
        x=json.loads(fit.x_fit) if isinstance(fit.x_fit,str) else fit.x_fit
        y=json.loads(fit.y_fit) if isinstance(fit.y_fit,str) else fit.y_fit
        positive=sorted(data.loc[data.Dose>0,'Dose'].unique()); pseudo=min(positive)/2
        ax.plot(x,y,color='black',lw=LINE,ls='--',alpha=.7)
        ax.scatter(data.Dose.where(data.Dose>0,pseudo),data.OD_final,s=10,color='grey',marker='x',linewidth=.7)
        ax.set_xscale('log');ax.set_xlim(pseudo/1.05,max(positive)*1.05);ax.set_ylim(-.01,1.01)
        ax.set_xticks(positive[::2],[f'{v:.2g}' for v in positive[::2]],rotation=45,ha='right',rotation_mode='anchor')
        ax.set_yticks([0,.25,.5,.75,1]);ax.set_xlabel('Levofloxacin (μg/mL)');ax.set_ylabel('OD$_{600}$');style(ax,False)
        for j,(key,color) in enumerate([('IC50','magenta'),('MIC','blue')]):
            value=float(fit[key]);threshold=float(fit[key.lower()+'_threshold'])
            ax.plot([pseudo,value,value],[threshold,threshold,0],color=color,lw=.7,ls=':')
            ax.text(.98,.93-j*.14,f'{key}: {value:.3g} μg/mL',transform=ax.transAxes,ha='right',fontsize=7.5,color=color)
    elif asset=='original_image32':
        data=pd.read_excel(ROOT/'working/figures/Supplemental Figure 3 - CefR/A - ResistanceEvo/pcr_evolution.xlsx')
        g=data.groupby('Day').Concentration.agg(['mean','std'])
        fig,ax=canvas(width,height,left=52,bottom=29)
        ax.errorbar(g.index,g['mean'],yerr=g['std'],fmt='o-',color='purple',lw=LINE,ms=2.5,elinewidth=ERROR,capsize=1.5)
        ax.set_xlabel('Day');ax.set_ylabel('Maximum concentration\n(μg/mL)');style(ax)
        ax.set_xticks(g.index[g.index%2==1])
    else:
        rel,sheet,strain=('Figure 3/A - survival','PLAC','MG1655') if asset=='original_image15' else ('Supplemental Figure 4 - Pb/A - Survival','ATEC','ATEC')
        raw=pd.read_excel(ROOT/'working/figures'/rel/'SurvivalData.xlsx',sheet_name=sheet)
        raw=raw[(raw.Strain==strain)&~raw.Culture.isin(['A','B','C','D'])].copy()
        raw['treated']=raw.Day%1!=0;raw['day']=raw.Day.where(~raw.treated,raw.Day-.5)
        grouping=['day','Strain','Culture', 'Drug' if sheet=='PLAC' else 'Evolution']
        before=raw[~raw.treated].set_index(grouping).CFU
        after=raw[raw.treated].set_index(grouping).CFU
        assert before.index.is_unique and after.index.is_unique
        data=(after/before*100).dropna().rename('survival_percent').reset_index()
        fig,ax=canvas(width,height,left=43,bottom=29)
        for _,g in data.groupby('Culture'):
            g=g.sort_values('day');ax.plot(g.day,g.survival_percent,color='.6',lw=.65,alpha=.5)
        mean=data.groupby('day').survival_percent.mean()
        ax.plot(mean.index,mean,color='forestgreen',lw=LINE)
        if sheet=='PLAC':
            for cutoff,color in [(32,'blue'),(17,'red')]:
                m=mean[mean.index<cutoff];ax.plot(m.index,m,color=color,lw=LINE)
        ax.set_xlabel('Day');ax.set_ylabel('Survival (%)');style(ax);ax.set_ylim(.0005,3000)
        ax.set_xticks(range(0,int(data.day.max())+1,15 if sheet=='PLAC' else 2))
    return fig,dict(kind='numerical_reconstruction',asset=asset,source_records=len(data),tick_pt=TICK,axis_pt=AXIS,spine_pt=SPINE)


def main():
    import sys
    sys.path.insert(0,str(ROOT/'scripts'))
    import reproduce_genomic_plots as gp
    gp.OUT.mkdir(parents=True,exist_ok=True)
    genomic={'original_image12':'Fig2f','original_image13':'Fig2g','original_image11':'Fig2h',
             'original_image34':'S3e','original_image44':'S4f','original_image49':'S6d'}
    OUT.mkdir(exist_ok=True)
    layout=json.loads((HERE/'assembly_layout.json').read_text());final=pd.read_csv(HERE/'out/final_statistics.csv')
    replacements={};audit=[]
    for figure in layout['figures']:
        scale=min(6.45*72/figure['width'],7.5*72/figure['height'])
        for panel in figure['panels']:
            asset=panel['asset'];rect=panel['rect'];width=(rect[2]-rect[0])*scale;height=(rect[3]-rect[1])*scale
            if asset.startswith('original_image') and asset not in ['original_image8','original_image17','original_image41','original_image5','original_image15','original_image32','original_image38',*genomic]:continue
            if asset in genomic:
                key=genomic[asset]; m,labels,cmap=gp.s6() if key=='S6d' else gp.endpoint(key)
                fig=gp.heatmap(m,labels,cmap,width,height,transpose=key not in ['S4f','S6d'],annot=key.startswith('Fig2'))
                details=dict(kind='genomic',matrix=key,mutation_identity_preserved=True,tick_pt=8.5)
            elif asset in ['original_image5','original_image15','original_image32','original_image38']:fig,details=other_numeric(asset,width,height)
            elif asset.endswith('_labels') and 'doubling' not in asset:fig,details=fit_panel(asset,width,height,final)
            elif 'doubling' in asset:fig,details=growth_panel(width,height)
            elif asset.startswith('original_image'):fig,details=legacy_mdk(asset,width,height)
            else:fig,details=modern_mdk(asset,width,height)
            stem=figure['name']+'__'+asset
            fig.canvas.draw()
            renderer=fig.canvas.get_renderer();bounds=[]
            legends = list(fig.legends)
            for ax in fig.axes:
                if ax.get_legend() is not None:
                    legends.append(ax.get_legend())
                for text in [*ax.texts,*ax.get_xticklabels(),*[t for t in ax.get_yticklabels() if ax.get_ylim()[0]<=t.get_position()[1]<=ax.get_ylim()[1]],ax.xaxis.label,ax.yaxis.label]:
                    if text.get_visible() and text.get_text():
                        box=text.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
                        bounds.append([text.get_text(),*[v*72 for v in box.extents]])
            for artist in [*legends, *fig.texts]:
                if artist.get_visible():
                    box=artist.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
                    bounds.append(['legend' if artist in legends else artist.get_text(),*[v*72 for v in box.extents]])
            bad=[x for x in bounds if x[1]<-.4 or x[2]<-.4 or x[3]>width+.4 or x[4]>height+.4]
            if bad:raise ValueError(f'{stem}: text outside panel: {bad}')
            for ext in ['pdf','svg','png']:
                fig.savefig(OUT/f'{stem}.{ext}',dpi=600 if ext=='png' else None)
            plt.close(fig)
            path=OUT/f'{stem}.pdf'
            with fitz.open(path) as saved:
                pdf_width,pdf_height=saved[0].rect.width,saved[0].rect.height
            replacements[figure['name']+'::'+asset]=dict(source=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                native_width=pdf_width,native_height=pdf_height,content_box=[0,0,pdf_width,pdf_height],content_aspect=pdf_width/pdf_height)
            audit.append(dict(figure=figure['name'],panel=panel['panel'],print_width_pt=width,print_height_pt=height,
                              text_bounds=bounds,boundary_issues=bad,**details))
    inputs = [Path(__file__), ROOT/'scripts/reproduce_genomic_plots.py', *sorted((ROOT/'data/genomics/generated_tables').glob('*.json')), HERE/'plot_final.py', HERE/'assembly_layout.json',
              *sorted((ROOT/'working/figures').rglob('SurvivalData.xlsx')),
              *sorted((ROOT/'working/figures').rglob('*.ipynb')), HERE/'fig5a_doubling.py', GROWTH_SOURCE, HERE/'stats.py', HERE/'report_results.py',
              *sorted((HERE/'cache').glob('*.pkl')), *sorted((HERE/'out').glob('*.csv'))]
    (OUT/'manifest.json').write_text(json.dumps({'panels':replacements,
        'layout_sha256':hashlib.sha256((HERE/'assembly_layout.json').read_bytes()).hexdigest(),
        'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}},indent=2)+'\n')
    (OUT/'render_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(f'Rendered {len(audit)} panels at manuscript size.')


if __name__=='__main__':main()
