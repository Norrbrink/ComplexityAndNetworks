# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 11:23:11 2022

@author: Alexander Norrbrink
"""

from Oslo_fast import Ricepile
from Oslo_algorithm import Sandpile
import matplotlib.pyplot as plt
import numpy as np

parameters = {
    # Use LaTeX to write all text
    "text.usetex": False,
    "font.family": "Times New Roman",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 16,
    "font.size": 16,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.figsize": (8,6)
}

plt.rcParams.update(parameters)

#%%
systems = [4, 8, 16, 32, 64, 128, 256, 512, 1024]
slowpiles = [Sandpile(i, 0.5) for i in systems]
piles = [Ricepile(i, 0.5) for i in systems]
#%% Drive systems 1100000 grains
for i in range(len(piles)):
    print(systems[i])
    piles[i].Oslo(1100000) 
#%%Test if Oslo Algorithm has the correct averaged height at i = 1
#piles[2].heighttest()
#piles[3].heighttest()
#plt.show()
#%% Test if Oslo Algorithm has the correct form and functionality of the steady state solutions for extrema values of p
sandp0 = Sandpile(16, 0)
sandp1 = Sandpile(16, 1)
sandp0.Oslo(1000)
sandp1.Oslo(1000)

sandp0.onnomtest()
sandp1.onnomtest()

#%% test tc for p = 0 and 1, tc= zth*L(L+1)/2
p0piles = [Sandpile(i, 0) for i in systems]
p1piles = [Sandpile(i, 1) for i in systems]
for i in range(len(systems)):
    p0piles[i].Oslo(4*systems[i]**2+10)
    p1piles[i].Oslo(systems[i]**2+10)
    print('P0: tc={}'.format(p0piles[i].steadstatetime))
    print('P1: tc={}'.format(p1piles[i].steadstatetime))
#%% test set of recurrent systems
configsystems = [1,2]#,3,4,5] 
colours = ['black', 'red', 'navy', 'green', 'purple']
configpiles = [Ricepile(i, 0.5) for i in configsystems]
for j in range(len(configpiles)):
    print(j)
    configpiles[j].Oslo(1100000, configs=True) 
    time = np.linspace(0, len(configpiles[j].numbofconfig), len(configpiles[j].numbofconfig))
    plt.plot(time, configpiles[j].numbofconfig, '.', c=colours[j])

plt.xlabel(r'time, (grains)')
plt.ylabel('Number of Configurations')
plt.savefig('testforconfigurations.png', dpi=1000)
plt.show()
        

#%% test Average avalanche size 
for i in range(len(systems)):
    print('Average Avalanche Size = {:.0f}'.format(np.mean(piles[i].avalanches[piles[i].steadstatetime + 1:])))
#%% Height at i=1 over time

ht = [j.height1time for j in piles]
times = [range(len(k.hoft), len(k.hoft)+len(k.height1time)) for k in piles]
for i in range(0, len(ht)):
    plt.plot(times[i], ht[i], '.', c='red', Zorder=10)#, label='steady state {}'.format(systems[i]))
plt.plot(range(0, len(piles[-1].hoft)), piles[-1].hoft, '.',  c='black', Zorder=0, label='Transient')
plt.xlabel('time, (grains)')
plt.ylabel('Height at i=1, (grains)')
plt.xlim(0, 1100000)
plt.ylim(0, 2000)
plt.plot([0], [0], '.', c='red', label='Steady states')
plt.legend()
plt.savefig('Heightat1.png', dpi=1000)
plt.show()

#%% Height Data Collapse

heightpiles = [[Ricepile(i, 0.5) for i in systems] for j in range(1000)]
heighttimes = np.arange(0, 1000, 1)

averagedheights = [np.array([0 for i in range(1000)]) for j in range(len(systems))]

for i in range(1000):
    print(i)
    for j in range(len(systems)):
        heightpiles[i][j].skiptosteadystate()
        heightpiles[i][j].Oslo(1000)
        averagedheights[j] = averagedheights[j] + np.array(heightpiles[i][j].height1time)
#%% Task 2e-g

logmeanht = [np.log(np.mean(averagedheights[j]/1000)) for j in range(len(systems))]
logsys = [np.log(systems[i]) for i in range(len(logmeanht))]
meanht = [(np.mean(averagedheights[j]/1000)) for j in range(len(systems))]
stdht = [np.log(np.sqrt(np.mean((averagedheights[j]/1000)**2) - (np.mean(averagedheights[j]/1000))**2)) for j in range(len(systems))]
errht = [np.std(averagedheights[j]/1000, ddof=1)/(np.sqrt(len(averagedheights[j]/1000))*systems[j]) for j in range(len(systems))]
#plt.errorbar(logsys, logmeanht, yerr=errht, fmt='.', c='black', capsize=4, elinewidth=1,  ecolor='black', label='<$h_1$>')
plt.xlim(0, 1050)
plt.ylim(0, 1600)
plt.errorbar(systems, meanht, yerr=errht, fmt='.', c='black', capsize=4, elinewidth=1,  ecolor='black', label='<$h_1$>')
fit = np.polyfit(systems[-3:], meanht[-3:], 1)
fit_fn = np.poly1d(fit)
a0 = fit[0] 
plt.plot(np.linspace(0, 1050, 1000), fit_fn(np.linspace(0, 1050, 1000)), c='red', label='Linear fit: y={:.2f}x + {:.2f}'.format(fit_fn[1], fit_fn[0]))
#plt.plot(range(9), np.arange(0.413, 9.413, 1), c='black', label='y = x + 0.45')
plt.xlabel('L, (grains)')
plt.ylabel(r'Average Height after $t_c$')
plt.plot([], [], ' ', label=r'$a_0$ = {:.2f}'.format(a0))
plt.legend()
plt.savefig('a0determination.png', dpi=1000)
plt.show()

plt.plot(logsys, stdht, '.', c='black')
fit = np.polyfit(logsys[-3:], stdht[-3:], 1)
fit_fn = np.poly1d(fit)
plt.plot(np.linspace(0, 8, 1000), fit_fn(np.linspace(0, 8, 1000)), c='red', label='Linear fit: y={:.2f}x + {:.2f}'.format(fit_fn[1], fit_fn[0]))
plt.xlabel('ln(L)')
plt.ylabel(r'ln(Standard Deviation)')
plt.legend()
plt.xlim(0, 8)
plt.ylim(0, 4)
plt.savefig('STDscaling.png', dpi=1000)
plt.show()
 
correctionht = [np.log(np.mean(averagedheights[j]/1000) - a0*systems[j]) for j in range(len(systems))]
plt.errorbar(logsys, correctionht, yerr=errht, fmt='.', c='black', capsize=4, elinewidth=1,  ecolor='black', label='<$h_1$>')
plt.xlim(0, 8)
plt.ylim(0, 4)
#plt.plot(systems, np.exp(0.413)*np.array(systems), c='black')
fit = np.polyfit(logsys[-3:], correctionht[-3:], 1)
fit_fn = np.poly1d(fit)
w1 = 1 - fit[0]
a1 = np.exp(fit[1])/a0
plt.plot(np.linspace(0, 8, 1000), fit_fn(np.linspace(0, 8, 1000)), c='red', label='Linear fit: y={:.2f}x + {:.2f}'.format(fit_fn[1], fit_fn[0]))
plt.plot([], [], ' ', label=r"$\omega_1$ = {:.2f}".format(w1))
plt.xlabel('ln(L)')
plt.ylabel(r'ln(Average Height after $t_c$ minus $a_0$L)')
plt.legend()
plt.savefig('w1determination.png', dpi=1000)
plt.show()

#%% Task 2d
plt.show()

for j in range(len(systems)):
    plt.plot(heighttimes/systems[j]**2, np.log(averagedheights[j]/(1000*systems[j])), '.', c='r')
    

plt.xlim(0, 70)
#plt.ylim(0, np.log(2000))
plt.xlabel(r'time since $t_c$, (grains)')
plt.ylabel('Height at i=1, (grains)')
plt.plot([0], [0], '.', c='red', label='Steady states')
plt.legend()
plt.savefig('heightcollapse.png', dpi=1000)
plt.show()      
         
#%% Task 2g
normalisations = [1099989, 1099939, 1099790, 1099019, 1096539, 1085740, 1043585, 874629, 201479]
for j in range(5, len(systems)):
    d = piles[j].probablities()
    plt.plot(np.array(list(d.keys()))/(a0*systems[j]), np.array(list(d.values()))/(normalisations[j]), '.', color='red')
              

plt.xlim(1.1, 1.25)
plt.ylim(0, 0.3)
plt.xlabel('h, (grains)')
plt.ylabel('P(h; L)')
plt.plot([0], [0], '.', c='red', label='Steady states')
plt.legend()
plt.savefig('heightcollapse.png', dpi=1000)
plt.show()      


#%% <tc> determination
tcs = []
for i in systems:
    tccurr = []
    for j in range(1000):
        current = Sandpile(i, 0.5)
        tccurr.append(current.tc())
    tcs.append(tccurr)
#%%
meantcs = [np.log(np.mean(i)) for i in tcs]
errtcs = [np.std(tcs[i], ddof=1)/(np.sqrt(len(tcs[i]))*systems[i]) for i in range(len(systems))]
logsys = [np.log(systems[i]) for i in range(len(meantcs))]
plt.errorbar(logsys, meantcs, yerr=errtcs, fmt='.', c='red', capsize=4, elinewidth=1,  ecolor='black', label='<$t_c$>')
plt.xlim(0, 8)
plt.ylim(0, 16)
plt.plot(range(9), np.arange(0, 17, 2), c='black', label='Power law - 2nd order')
plt.xlabel('ln(L)')
plt.ylabel('ln(<$t_c$>)')

plt.legend()
plt.savefig('tcagainstL.png', dpi=1000)
plt.show()




#%% Avalanche size probabilities
for i in range(len(piles)):
    x, y = piles[i].avalancheProb()
    logx = [np.log(j) for j in x]
    logy = [np.log(j) for j in y]
    plt.plot(logx, logy, '.', color='black')
fit = np.polyfit(logx[10:21], logy[10:21], 1)
fit_fn = np.poly1d(fit)
plt.plot(np.linspace(0, 16, 1000), fit_fn(np.linspace(0, 16, 1000)), c='red', label='Linear fit: y={:.2f}x {:.2f}'.format(fit_fn[1], fit_fn[0]))
    
plt.xlim(0, 16)
plt.xlabel('s')
plt.ylabel('P(s; L)')
plt.legend()
plt.savefig('probplot.png', dpi=1000)
plt.show()
#%% Finite Size scaling
linestyles = [(0, (1,10)), 'dotted', 'dashed', 'solid']
for i in range(5, len(piles)):
    x, y = piles[i].avalancheProb()
    logx = [(j/(systems[i]**2.25)) for j in x]
    logy = [(y[j]*x[j]**abs(fit_fn[1])) for j in range(len(y))]
    plt.plot(logx, logy, linestyle = linestyles[i-5], color='black', label='L = {}'.format(systems[i]))
#fit = np.polyfit(logx, logy, 1)
#fit_fn = np.poly1d(fit)
#plt.plot(np.linspace(0, 16, 1000), fit_fn(np.linspace(0, 16, 1000)), c='red', label='Linear fit: y={:.2f}x {:.2f}'.format(fit_fn[1], fit_fn[0]))
    
plt.ylim(10**-3, 10**0)
plt.xlim(10**-7, 10**0)
plt.xlabel(r's/$L^D$')
plt.ylabel(r'$s^{\tau_s}$ P(s; L)')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('finitesizesacling.png', dpi=1000)
plt.show()
#%% Measuring the kth moment
k = [1,2,3,4]
colours = ['navy', 'red', 'green', 'black']
moments = []
for i in k:
    kmom = []
    for j in piles:
        kmom.append(np.log(j.kthmoment(i)))
    moments.append(kmom)

exponents = []

for i in range(4):
    plt.plot(logsys, moments[i], '.',  c= colours[i], label='Moment k={}'.format(k[i]))
    grad, inter = np.polyfit(logsys, moments[i], 1)
    exponents.append(grad)
    if inter < 0:
        plt.plot(range(9), [grad*i + inter for i in range(9)], c=colours[i], label='y = {:.2f} x - {:.2f}'.format(grad, abs(inter)))
    else:
        plt.plot(range(9), [grad*i + inter for i in range(9)], c=colours[i], label='y = {:.2f} x + {:.2f}'.format(grad, inter))
plt.xlim(0, 8)
plt.ylim(0, 60)
plt.xlabel('ln(L)')
plt.ylabel(r'ln(<$s^k$>)')
plt.legend()
plt.savefig('kthmoments.png', dpi=1000)
plt.show()

plt.plot(k, exponents, '.', c='black')
plt.xlim(0, 5)
plt.ylim(0, 10)
plt.xlabel('k')
plt.ylabel(r'D(1+k-$\tau_s$)')
fit = np.polyfit(k, exponents, 1)
fit_fn = np.poly1d(fit)
plt.plot(np.linspace(0, 5, 1000), fit_fn(np.linspace(0, 5, 1000)), c='red', label='Linear fit: y={:.2f}x {:.2f}'.format(fit_fn[1], fit_fn[0]))
plt.plot([], [], ' ', label=r"D = {:.2f}, $\tau_s$ = {:.2f}".format(fit_fn[1], 1-fit_fn[0]/fit_fn[1]))
plt.legend()
plt.savefig('Ddetermination.png', dpi=1000)
plt.show()
