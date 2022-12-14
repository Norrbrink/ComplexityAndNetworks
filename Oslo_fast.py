# -*- coding: utf-8 -*-
"""
Created on Sat Feb 12 16:46:36 2022

@author: Alexander
"""

import random
import numpy as np
from collections import Counter
from logbin_2020 import logbin
from numba import jit

@jit(nopython=True)
def numbarelax(L, z, thresholds, steadystate, options):
    '''Relaxes all sites with slopes larger than their threshold slopes'''
    relaxations = 0
    for i in range(0, L):
        if z[i] > thresholds[i]: #Relaxes Site i=1
            if i == 0:
                z[i] -= 2
                z[i+1] += 1
                relaxations += 1
            elif i < L - 1: #Relaxes Site 1<i<L
                    z[i] -= 2
                    z[i+1] += 1
                    z[i-1] += 1
                    relaxations += 1
            else: #Relaxes Site i=L
                steadystate = True
                z[i] -= 1
                z[i-1] += 1
                relaxations += 1
            thresholds[i] = np.random.choice(options) #Changes the new threshold slopes
    return relaxations, z, thresholds, steadystate

@jit(nopython=True)
def numbaconfig(config, configurations):
    if config.tolist() not in configurations.tolist():
        return True
    else:
        return False
    
class Ricepile:
    '''An Oslo Model for ricegrains, this is a quick version to run a longer simulation'''
    def __init__(self, L, p):
        self.L = L
        self.z = np.array([0 for i in range(L)])
        self.height = sum(self.z)
        self.options = np.array([1, 2])
        self.weighting= [p, 1 - p]
        self.thresholds = np.array([random.choices(self.options, weights=self.weighting, k=1)[0] for i in range(L)])
        self.avalanches = []
        self.height1time = [] #Array of Heights for Steady State
        self.hoft = [0] #Array of Heights for Transient State
        self.steadystate = False
        self.t = 0
        self.configurations = []
        self.numbofconfig = []
        
    def __drive(self):
        '''Adds one rice-grain to Site i=1'''
        self.z[0] += 1
        self.t += 1
        
    def __relax(self):
        relaxations, self.z, self.thresholds, self.steadystate = numbarelax(self.L, self.z, self.thresholds, self.steadystate, self.options)
        return relaxations
    
    def Oslo(self, N, configs=False):
        '''Drives and Relaxes N times'''
        for i in range(N):
            self.__drive()
            s = 0
            relaxed = self.thresholds - self.z
            while len([*filter(lambda x: x < 0, relaxed)]) > 0: 
                s += self.__relax()
                relaxed = self.thresholds - self.z
            self.height = sum(self.z)    
            if self.steadystate:
                self.height1time.append(self.height)
                if configs:
                    if self.z.tolist() not in self.configurations:
                        self.configurations.append(self.z.tolist())
                    self.numbofconfig.append(len(self.configurations))
            else:
                self.hoft.append(self.height)
                self.steadstatetime = self.t + 1
            self.avalanches.append(s)
        
    def heighttest(self):
        '''Function that calculates the time averaged height at site 1 for all time after tc'''
        mean = sum(self.height1time)/len(self.height1time)
        std = np.std(self.height1time, ddof=1)
        print(mean, u'\u00B1', std/np.sqrt(len(self.height1time)))
    
    def onnomtest(self, skipped = False):
        '''Function that calculates the average avalanche size for all time after tc'''
        if skipped:
            mean= np.mean(self.avalanches)
            std = np.std(self.avalanches, ddof=1)
            print(mean, u'\u00B1', std/np.sqrt(len(self.avalanches)))
        else:
            mean = np.mean(self.avalanches[self.steadstatetime + 1:])
            std = np.std(self.avalanches[self.steadstatetime + 1:], ddof=1)
            print(mean, u'\u00B1', std/np.sqrt(len(self.avalanches)+ 1 - self.steadstatetime))
        
    def tc(self):
        '''Function that calculates the number of grains before reaching Steady State'''
        SUM = 0
        for i in range(self.L):
            SUM += self.thresholds[i]*(i + 1)
        return SUM
    
    def timeAveragedHeight(self):
        '''Function that calculates the time averaged height at site 1 for all time after tc'''
        return np.mean(self.height1time)
    
    def timeHeightSTD(self):
        '''Function that calculates the STD of the height at site 1 for all time after tc'''
        return np.std(self.height1time, ddof=1)
    
    def probablities(self):
        '''Function that calculates the frequency of heights'''
        print('Total number of configurations considered: {}'.format(len(self.height1time)))
        return Counter(self.height1time)

    def avalancheProb(self, scale = 1.2, zeroes = False, skipped = False):
        '''Function that applies the logbin function on the Ricepile'''
        if skipped:
            x, y = logbin(self.avalanches, scale, zeroes)
        else:
            x, y = logbin(self.avalanches[self.steadstatetime:], scale, zeroes)
        return x, y
    
    def kthmoment(self, k, skipped = False):
        '''Function that calculates the kth moment'''
        SUM = 0
        if skipped:
            T = len(self.avalanches)
            tstart = 0
        else:
            T = len(self.avalanches[self.steadstatetime:])
            tstart = self.steadstatetime
        for i in self.avalanches[tstart:]:
            SUM += (i**k)/T
        return SUM
    
    def skiptosteadystate(self):
        '''Function that calculates the tc and then moves the simulation there'''
        if self.steadystate:
            print('Already Steady State')
        else:
            self.t = self.tc()
            self.z = self.thresholds.copy()
            self.height = sum(self.z)
            self.steadystate = True
            self.steadstatetime = self.t.copy()