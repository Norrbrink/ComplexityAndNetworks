# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 09:32:11 2022

@author: Alexander
"""

from Network import Network
import matplotlib.pyplot as plt
import numpy as np
import random
from sklearn.metrics import r2_score

mvalues = [8, 16, 32, 64, 128]
distributions = []
Ns=[]
kbins = []
for i in range(len(mvalues)):
    print(i)
    evalue = (mvalues[i]*(mvalues[i]-1))/2
    a = Network(rand=True, N=mvalues[i], E=mvalues[i]*(mvalues[i]-1)//2)
    #a.status()
    a.drive(50000, mvalues[i], 'P1')
    k, counts, N = a.degreeprobabilities()
    distributions.append(np.array(counts))
    kbins.append(k)
    Ns.append(N)
    #plt.plot(k, np.log(np.array(counts)))
#plt.xlim(0, 0.7)
#%%


colours = ['red', 'navy', 'orange', 'green', 'purple']

for i in range(len(mvalues)):
    x = np.log(kbins[i]/(mvalues[i]+1))
    y = -3*x
    plt.plot(np.log(kbins[i]/(mvalues[i]+1)), np.log(distributions[i]*(mvalues[i])), '.', c=colours[i], label='m={}'.format(mvalues[i]))
    print(r2_score(y, np.log(distributions[i]*(mvalues[i]))))
plt.plot(x, y, c='black', label=r'Power Law: $\gamma$ =3 ')
plt.xlim(0)
plt.legend()
plt.xlabel('ln(k/m)')
plt.ylabel(r'ln(mP$\infty$(k))')
plt.savefig('P1distribution.png')
plt.show()
for i in range(len(mvalues)):
    x = kbins[i]
    y = (2*mvalues[i]*(mvalues[i]+1))/(x*(x+1)*(x+2))
    plt.plot(x, y-distributions[i], '.', c=colours[i], label='m={}'.format(mvalues[i]))   
#plt.plot(x, y, c='black', label=r'Power Law: $\gamma$ =3 ')
plt.xlim(0)
plt.legend()
plt.xlabel('k')
plt.ylabel(r'theoretical minus numerical P$\infty$(k)')
plt.savefig('P1comp.png')
#%%
attachednode = []
bins = [-0.5, 0.5, 0.5, 1.5, 1.5, 2.5, 2.5, 3.5, 3.5, 4.5]
for j in range(10000):
    b = Network(rand=False, N=5, connections=[[2, 1], [2, 3], [3, 0], [1, 3], [0, 2], [0, 1], [3, 4]])
    b.drive(6, 1, 'P1')
    attachednode.append(b.adjacencylist[-1][0])
plt.hist(attachednode, bins=bins,  density=True, facecolor='r', edgecolor='black')
plt.xlabel('Attached Node')
plt.ylabel('Fraction of total attachments')

#%%
Nvalues = [500, 1000, 5000, 10000, 25000, 50000]
mvalues = 32
k1distributions = []
k1Ns=[]
k1bins = []
k1s = []
for i in range(len(Nvalues)):
    print(i)
    evalue = (mvalues*(mvalues-1))/2
    d = Network(rand=True, N=mvalues, E=mvalues*(mvalues-1)//2)
    #a.status()
    d.drive(Nvalues[i], mvalues, 'P1')
    k1s.append(max([len(i) for i in d.adjacencylist]))
#%%

plt.plot(Nvalues, k1s, '.',  c='r', label='Measured')
x = np.linspace(0, 50000, 1000)
y = -0.5 + np.sqrt(1+4*x*32*(33))/2
#plt.xlim(0)
#plt.ylim(0)
plt.plot(x,y, c='black', label='Theoretical')
plt.legend()
plt.savefig('P1K1.png')
plt.show()
#plt.xlim(0, 0.7)

#%%
mvalues = [8, 16, 32, 64, 128]
P2distributions = []
P2Ns=[]
P2kbins = []
for i in range(len(mvalues)):
    print(i)
    evalue = (mvalues[i]*(mvalues[i]-1))/2
    e = Network(rand=True, N=mvalues[i], E=mvalues[i]*(mvalues[i]-1)//2)
    #a.status()
    e.drive(50000, mvalues[i], 'P2')
    k, counts, N = e.degreeprobabilities()
    P2distributions.append(np.array(counts))
    P2kbins.append(k)
    P2Ns.append(N)
    #plt.plot(k, np.log(np.array(counts)))
#plt.xlim(0, 0.7)
#%%
colours = ['red', 'navy', 'orange', 'green', 'purple']

for i in range(len(mvalues)):
    x = P2kbins[i]
    y = mvalues[i]**(x-mvalues[i])/(1+mvalues[i])**(1+x-mvalues[i])
    plt.plot(P2kbins[i], P2distributions[i], '.', c=colours[i], label='m={}'.format(mvalues[i]))
    #print(r2_score(y, np.log(distributions[i]*(mvalues[i]))))
    plt.plot(x, y, c='black', label=r'Theoretical Solution')
plt.xlim(0)
plt.legend()
plt.xlabel('k)')
plt.ylabel(r'P$\infty$(k))')
plt.savefig('P2distribution.png')
plt.show()

#%%
Nvalues = [500, 1000, 5000, 10000, 50000]
mvalues = 32
k1distributions = []
k1Ns=[]
k1bins = []
k1s = []
for i in range(len(Nvalues)):
    print(i)
    evalue = (mvalues*(mvalues-1))/2
    d = Network(rand=True, N=mvalues, E=mvalues*(mvalues-1)//2)
    #a.status()
    d.drive(Nvalues[i], mvalues, 'P2')
    k1s.append(max([len(i) for i in d.adjacencylist]))
#%%

plt.plot(Nvalues, k1s, '.',  c='r', label='Measured')
x = np.linspace(0, 50000, 1000)
y = 32 - np.log(x)/(np.log(32)-np.log(33))
#plt.xlim(0)
#plt.ylim(0)
plt.plot(x,y, c='black', label='Theoretical')
plt.legend()
plt.savefig('P2K1.png')
plt.show()
#plt.xlim(0, 0.7)
