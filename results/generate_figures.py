#!/usr/bin/env python3
"""Generate all manuscript figures and a seeded finite-shot numerical experiment."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/finite-controller-mpl')
from pathlib import Path
import json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
ROOT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.bbox':'tight'})
NAVY='#17324d';TEAL='#00817a';RED='#b84a48';GREY='#74828c'
def save(fig,name):
    fig.savefig(ROOT/'figures'/f'{name}.pdf');fig.savefig(ROOT/'figures'/f'{name}.png',dpi=180);plt.close(fig)
# Exact four-frame example: plot derived probabilities, no measured data.
a=(2-math.sqrt(2))/4;p=np.array([a,1-a]);q0=np.array([0,.5]);q1=np.array([.5,1.])
fig,ax=plt.subplots(figsize=(5.2,4.0))
ax.plot([0,.5],[.5,1],color=GREY,lw=1.5,label='Mixtures of the two tables')
for q,label in [(q0,'Identity table'),(q1,'Shift table')]:
    ax.scatter(*q,s=55,color=NAVY,zorder=4);ax.annotate(label,q,xytext=(8,10) if label=='Identity table' else (8,-16),textcoords='offset points',color=NAVY)
ax.add_patch(Rectangle(p-.25,.5,.5,facecolor=TEAL,alpha=.10,edgecolor=TEAL,lw=1.5))
ax.scatter(*p,s=65,marker='D',color=TEAL,zorder=5,label='Ideal quantum pair')
ax.scatter(.25,.75,s=32,color=RED,zorder=5,label='Best table mixture')
ax.set(xlim=(-.16,.65),ylim=(.4,1.17),xlabel=r'$p(G)$',ylabel=r'$p(GA)$')
ax.text(-.12,.425,'Shaded square: tolerance 1/4 about the ideal pair',fontsize=8,color=TEAL)
ax.legend(loc='upper left',frameon=False);ax.set_aspect('equal',adjustable='box');save(fig,'joint_test_geometry')
# Repeated local benchmark.
w=json.loads((ROOT/'data/repeated_witness.json').read_text());from fractions import Fraction
seq=np.array([float(Fraction(s)) for s in w['sequence']]);n=np.arange(194)
theta=.5*math.acos(-(2+math.sqrt(2))/4);ideal=(6-math.sqrt(2))/17*(1-np.cos(2*n*theta));finite=seq[n%17];err=np.abs(finite-ideal)
fig,axs=plt.subplots(2,1,figsize=(6.7,4.9),gridspec_kw={'height_ratios':[1,1.05]})
axs[0].plot(n[:52],ideal[:52],'-',color=TEAL,lw=1.2,label='Ideal $(T,H)^n$');axs[0].plot(n[:52],finite[:52],'o--',color=NAVY,ms=2.7,lw=.8,label='Certified period-17 table')
axs[0].set(ylabel='Outcome-1 probability',xlabel='Repeated blocks (first 51 shown)',ylim=(-.03,.82));axs[0].legend(frameon=False,loc='upper right',ncol=2)
axs[1].plot(n,err,color=NAVY,lw=.8);axs[1].axhline(1/3,color=RED,ls='--',lw=1,label='Threshold 1/3');axs[1].scatter([193],[err[193]],color=RED,s=27,zorder=4);axs[1].annotate('First crossing: 193',(193,err[193]),xytext=(110,.385),arrowprops={'arrowstyle':'->','color':RED},fontsize=8,color=RED)
axs[1].set(xlabel='Repeated blocks',ylabel='Absolute probability error',ylim=(-.01,.43),xlim=(0,197));axs[1].legend(frameon=False,loc='upper left');fig.tight_layout(h_pad=1.4);save(fig,'repeated_horizon')
# Finite shot simulation under the four-frame target and both null models.
# Decision uses distance to the two fixed-table model vectors, with Hoeffding radius.
rng=np.random.default_rng(20261003);alpha=.05;reps=20000
Ns=[25,50,100,200,400,800];rows=[]
for N in Ns:
    r=math.sqrt(math.log(4/alpha)/(2*N))
    row={'shots_per_circuit':N,'confidence_radius':r}
    for name,truth in [('quantum',p),('identity',q0),('shift',q1)]:
        freq=rng.binomial(N,truth,size=(reps,2))/N
        distance=np.minimum(np.max(np.abs(freq-q0),axis=1),np.max(np.abs(freq-q1),axis=1))
        rate=float(np.mean(distance>r));row[name+'_rejection_rate']=rate
    rows.append(row)
(ROOT/'data/sampling_experiment.json').write_text(json.dumps({'seed':20261003,'replications_per_setting':reps,'alpha':alpha,'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__,'rows':rows},indent=2))
fig,ax=plt.subplots(figsize=(5.8,3.2))
for name,col,lab,mark in [('quantum',TEAL,'Ideal quantum probabilities','o'),('identity',NAVY,'Identity-table model','s'),('shift',RED,'Shift-table model','^')]:
    ax.plot(Ns,[r[name+'_rejection_rate'] for r in rows],marker=mark,color=col,label=lab,ms=4,lw=1.2)
ax.axhline(.05,color=GREY,lw=.8,ls=':');ax.set_xscale('log',base=2);ax.set_xticks(Ns,labels=[str(x) for x in Ns]);ax.set(xlabel='Shots per circuit',ylabel='Fraction rejecting the model class',ylim=(-.04,1.04));ax.legend(frameon=False,loc='center right');fig.tight_layout();save(fig,'sampling_power')
print('Generated three vector figures and finite-shot experiment data.')
