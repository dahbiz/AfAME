#!/usr/bin/env python3
# AfAME - parallel-tempering search and algebraic certification of
# absolutely maximally entangled (AME) states.
# Copyright (C) 2026 Zakaria Dahbi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
"""Regenerate the four external figures loaded by mmes.tex.
Run: python3 data/plot_paper_figures.py
Requires Python 3, NumPy and Matplotlib; no TeX installation for plotting.
Input files live beside this script. The CRT diagram is drawn inside mmes.tex.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/afame-pra-matplotlib')
from pathlib import Path
from collections import Counter
import csv, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
BLUE='#0072B2'; ORANGE='#D55E00'; GREEN='#009E73'; PURPLE='#7A5195'; GREY='#6B7280'
plt.rcParams.update({'font.family':'serif','font.serif':['STIXGeneral'], 'mathtext.fontset':'stix',
 'font.size':9,'axes.labelsize':9,'axes.titlesize':9.5,'xtick.labelsize':8,'ytick.labelsize':8,
 'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(ROOT/(name+'.pdf'),bbox_inches='tight',pad_inches=.04)
 fig.savefig(ROOT/(name+'.png'),dpi=220,bbox_inches='tight',pad_inches=.04)
 plt.close(fig)
def rows(name):
 with (HERE/name).open(newline='') as f:return list(csv.DictReader(f))

def tradeoffs(data):
 fig,ax=plt.subplots(1,2,figsize=(7.1,3.15),layout='constrained')
 x=np.array([d['tails'][1]/70 for d in data]);y=np.array([d['tails'][0]/70 for d in data])
 assert len(data)==817
 a=ax[0];a.scatter(x,y,s=12,c='#B9BEC5',alpha=.8,edgecolors='none')
 xx=np.linspace(66/70,1,50);a.plot(xx,132/35-3*xx,'--',color=GREY,lw=1,zorder=2)
 for label,u,v,m,c,offset in [('A',1,54/70,'o',BLUE,(-14,-2)),('B',68/70,60/70,'s',ORANGE,(8,0)),('C',66/70,66/70,'^',PURPLE,(9,0))]:
  a.scatter(u,v,s=45,marker=m,color=c,edgecolors='white',linewidth=.5,zorder=4)
  a.annotate(label,(u,v),xytext=offset,textcoords='offset points',weight='bold',va='center')
 a.set(xlim=(.925,1.012),ylim=(.74,.966),xticks=[.94,.96,.98,1],yticks=[.75,.8,.85,.9,.95],
       xlabel=r'$p_3=\Pr(r\geq 3)$',ylabel=r'$p_4=\Pr(r=4)$',title='(a) Eight-qutrit entanglement yield')
 a.grid(alpha=.17);a.set_axisbelow(True)
 inset=a.inset_axes([.11,.14,.39,.34]);inset.scatter(x,y,s=3,color='#ABB1B9')
 inset.plot([1,66/70],[54/70,66/70],color=BLUE,lw=.8)
 inset.set(xlim=(-.04,1.05),ylim=(-.04,1.05),xticks=[0,1],yticks=[0,1])
 inset.tick_params(labelsize=6,pad=1,length=2);inset.set_title('All 817 classes',fontsize=7,pad=3)
 inset.set_xlabel(r'$p_3$',fontsize=7,labelpad=-6);inset.set_ylabel(r'$p_4$',fontsize=7,labelpad=-3)
 a=ax[1]
 # Upper envelope only; the segment is not a continuum of deterministic states.
 a.plot([44,48],[48,48],':',color=GREY,lw=1.1)
 a.plot([48,52],[56,56],'--',color=GREY,lw=1.1)
 a.scatter([48],[48],facecolors='white',edgecolors=GREY,s=24,zorder=3)
 a.scatter([44],[48],s=45,color=BLUE,marker='o',zorder=4)
 a.scatter([48],[56],s=45,color=PURPLE,marker='s',zorder=4)
 a.scatter([44],[34],s=35,color=GREY,marker='x',zorder=4)
 a.annotate('(44, 48)',(44,48),xytext=(7,8),textcoords='offset points')
 a.annotate('(48, 56)',(48,56),xytext=(-4,-15),textcoords='offset points',ha='right')
 a.annotate('Earlier (44, 34)',(44,34),xytext=(8,1),textcoords='offset points',fontsize=8)
 a.annotate(r'$F_6\leq56$',(51.7,56),xytext=(0,6),textcoords='offset points',ha='right',fontsize=8)
 a.set(xlim=(42.8,52),ylim=(30,60),xticks=[44,46,48,50,52],yticks=[32,40,48,56],
       xlabel=r'Sector cost $\mathcal{C}_6$',ylabel=r'Maximally mixed four-party subsets $F_6$',title='(b) Binary–ternary product trade-off')
 a.grid(alpha=.17);a.set_axisbelow(True)
 return fig

def main():
 save(tradeoffs(json.loads((HERE/'ternary_class_profiles.json').read_text())),'fig_tradeoffs')
 records=rows('recombination_benchmark.csv')
 assert len(records)==150
 rng=np.random.default_rng(20260923);summary=[]
 fig,ax=plt.subplots(figsize=(7.1,2.95),layout='constrained')
 for j,n in enumerate((10,11,12)):
  z=[r for r in records if int(r['N'])==n];g=np.array([int(r['gain']) for r in z])
  assert len(g)==50 and all(int(r['joint_cost'])-int(r['recombined_cost'])==int(r['gain']) and int(r['proposals'])==2000 for r in z)
  ci=np.quantile(rng.choice(g,(50000,50)).mean(axis=1),[.025,.975]);counts=Counter(g)
  # Deterministic horizontal stacking prevents coincident observations disappearing.
  for gain,count in sorted(counts.items()):
   offsets=(np.arange(count)-(count-1)/2)*.016
   ax.scatter(j+offsets,[gain]*count,s=9,color=BLUE,alpha=.85,edgecolors='none',zorder=3)
  ax.errorbar(j+.38,g.mean(),yerr=[[g.mean()-ci[0]],[ci[1]-g.mean()]],fmt='D',ms=4,color=ORANGE,capsize=3,lw=1,zorder=5)
  ax.text(j,37.2,f'{sum(g>0)}/50 improved; mean {g.mean():.2f}',ha='center',fontsize=8)
  summary.append({'N':n,'joint_mean':np.mean([int(r['joint_cost']) for r in z]),'recombined_mean':np.mean([int(r['recombined_cost']) for r in z]),'mean_gain':g.mean(),'ci95':ci.tolist(),'strict_gains':int(sum(g>0))})
 ax.axhline(0,lw=.6,color=GREY);ax.set(xlim=(-.47,2.62),ylim=(-1.7,40),xticks=[0,1,2],xticklabels=[r'$N=10$',r'$N=11$',r'$N=12$'],yticks=[0,10,20,30],ylabel=r'Paired cost reduction $\Delta\mathcal{C}_6$',xlabel='Party number (50 matched seeds per group)')
 ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
 ax.legend(handles=[Line2D([],[],marker='o',ls='',color=BLUE,ms=3,label='One seed'),Line2D([],[],marker='D',ls='-',color=ORANGE,ms=4,label='Mean and 95% bootstrap interval')],loc='upper left',bbox_to_anchor=(0,1.18),ncol=2,frameon=False)
 save(fig,'fig_recombination')
 (HERE/'figure_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
 near=rows('near_ame_entropies.csv')
 pairs=sorted((r for r in near if int(r['N'])==4 and int(r['k'])==2),key=lambda r:r['subset'])
 assert len(pairs)==6 and sum(r['passing']=='False' for r in pairs)==2
 a2=[];a3=[]
 for r in pairs:
  u,v=map(int,r['sector_ranks'].split('+'));a2.append(u);a3.append(v*math.log2(3))
 fig,ax=plt.subplots(figsize=(3.45,2.55),layout='constrained');x=np.arange(6)
 ax.bar(x,a2,width=.65,color=BLUE,label=r'$\mathbb{F}_2$')
 ax.bar(x,a3,bottom=a2,width=.65,color='#E69F00',hatch='///',lw=.3,edgecolor='white',label=r'$\mathbb{F}_3$')
 ax.axhline(2*math.log2(6),color=GREY,ls='--',lw=.8,label='Maximum')
 for i,r in enumerate(pairs):
  if r['passing']=='False':ax.text(i,float(r['entropy_bits'])+.08,'−1 bit',ha='center',fontsize=7,color=ORANGE)
 ax.set(xticks=x,xticklabels=['{'+','.join(r['subset'])+'}' for r in pairs],ylim=(0,5.7),xlabel=r'Two-party subset $S$',ylabel=r'$S(\rho_S)$ (bits)')
 ax.legend(ncol=3,frameon=False,loc='lower center',bbox_to_anchor=(.5,1.02),handlelength=1.1,columnspacing=.8)
 save(fig,'fig_four_six')
 seven=[r for r in near if int(r['N'])==7];assert len(seven)==63
 fig,axs=plt.subplots(1,2,figsize=(7.1,2.7),layout='constrained');ax=axs[0]
 for k in (1,2,3):
  z=[r for r in seven if int(r['k'])==k];counts=Counter(round(float(r['entropy_bits']),7) for r in z)
  for entropy,count in sorted(counts.items()):
   ax.scatter(k,entropy,s=20+2*count,color=BLUE if abs(entropy-k)<1e-8 else ORANGE,marker='o' if abs(entropy-k)<1e-8 else 's',zorder=3)
   ax.annotate(f'{count}/{len(z)}',(k,entropy),xytext=(8,0),textcoords='offset points',va='center',fontsize=8)
 ax.plot([1,2,3],[1,2,3],ls='--',color=GREY,lw=.8,label='AME limit')
 ax.set(xlim=(.7,3.55),ylim=(.75,3.3),xticks=[1,2,3],yticks=[1,2,3],xlabel=r'Subset size $|S|$',ylabel=r'$S(\rho_S)$ (bits)',title='(a) Seven-qubit entropy profile');ax.legend(frameon=False,loc='upper left')
 ax=axs[1];counts=Counter(round(float(r['deficit_bits']),7) for r in seven if int(r['k'])==3);assert counts=={0.:32,1.:3}
 for x,num,color,hatch in [(0,32,BLUE,''),(1,3,ORANGE,'///')]:
  ax.bar(x,num,width=.55,color=color,hatch=hatch,edgecolor='white');ax.text(x,num+.8,str(num),ha='center',fontsize=9)
 ax.set(xlim=(-.6,1.6),ylim=(0,37),xticks=[0,1],yticks=[0,8,16,24,32],xlabel=r'Entropy deficit $3-S(\rho_S)$ (bits)',ylabel='Three-party subsets',title='(b) Balanced-cut deficits')
 save(fig,'fig_seven_qubit')
 print('Four manuscript figures regenerated; 150 paired benchmark rows and profile counts validated.')
if __name__=='__main__':main()
