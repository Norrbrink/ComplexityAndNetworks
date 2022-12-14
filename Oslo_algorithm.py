import random
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from logbin_2020 import logbin

class Sandpile:
    '''An Oslo Model for ricegrains, this is a slower version but is more robust'''
    def __init__(self, L, p):
        self.L = L 
        self.pile = np.array([0 for i in range(L)])
        self.z = np.array([0 for i in range(L)])
        self.options = np.array([1, 2])
        self.weighting= [p, 1 - p]
        self.thresholds = np.array([random.choices(self.options, weights=self.weighting, k=1)[0] for i in range(L)])
        self.avalanches = []
        self.height1time = []
        self.hoft = []
        self.steadystate = False
        
    def __drive(self):
        '''Adds one rice-grain to Site i=1'''
        self.pile[0] += 1
        self.z[0] += 1
        
    def __relax(self):
        '''Relaxes all sites with slopes larger than their threshold slopes'''
        relaxations = 0
        for i in range(0, self.L):
            if self.z[i] > self.thresholds[i]:
                if i == 0:
                    self.z[i] -= 2
                    self.z[i+1] += 1
                    self.pile[i] -= 1
                    self.pile[i+1] += 1
                    relaxations += 1
                elif i < self.L - 1:
                    self.z[i] -= 2
                    self.z[i+1] += 1
                    self.z[i-1] += 1
                    self.pile[i] -= 1
                    self.pile[i+1] += 1
                    relaxations += 1
                else:
                    if self.steadystate == False:
                        self.steadstatetime = len(self.avalanches)
                        self.steadystate = True
                    self.z[i] -= 1
                    self.z[i-1] += 1
                    self.pile[i] -= 1
                    relaxations += 1
                self.thresholds[i] = random.choices(self.options, weights=self.weighting, k=1)[0]
        return relaxations
    
    def Oslo(self, N):
        '''Drives and Relaxes N times'''
        for i in range(N):
            self.__drive()
            s = 0
            relaxed = self.thresholds - self.z
            while len([*filter(lambda x: x < 0, relaxed)]) > 0:
                s += self.__relax()
                relaxed = self.thresholds - self.z
            if self.steadystate:
                self.height1time.append(self.pile[0])
            else:
                self.hoft.append(self.pile[0])
            self.avalanches.append(s)
        self.visualise()
          
    def visualise(self, Print=False):
        '''Visualises the current configuration of the Pile'''
        x = [i+1 for i in range(self.L)]
        xcon = [i for i in range(self.L+2)]
        container = [1 for i in range(self.L+1)]
        #bottom = [1 for i in range(self.L)]
        container[0] = 2*self.L
        container.append(0)
        if Print:
            print('Heights:', self.pile)
            print('Slopes:', self.z)
            print('Threshold Slopes', self.thresholds)
            print('Avalanches:', self.avalanches)
        fig, ax = plt.subplots()
        ax.bar(xcon, container, width = 1, align ='edge', color = 'black')
        ax.bar(x, self.pile, width=1, align='edge', edgecolor= 'black', color='gray', bottom = container[1::-2])
        ax.set_xlim(0, self.L+2)
        ax.set_ylim(0, 2*self.L) 
        ax.set_ylabel('Height of at Site, (Grains)')
        ax.set_xlabel('Site')
        
    def heighttest(self):
        '''Function that calculates the time averaged height at site 1 for all time after tc'''
        mean = sum(self.height1time)/len(self.height1time)
        std = np.std(self.height1time, ddof=1)
        print(mean, u'\u00B1', std/np.sqrt(len(self.height1time)))
    
    def onnomtest(self):
        '''Function that calculates the average avalanche size for all time after tc'''
        mean = np.mean(self.avalanches[self.steadstatetime:])
        std = np.std(self.avalanches[self.steadstatetime:], ddof=1)
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

    def avalancheProb(self, scale = 1., zeroes = False):
        '''Function that applies the logbin function on the Ricepile'''
        x, y = logbin(self.avalanches[self.steadstatetime:], scale, zeroes)
        return x, y
    
    def kthmoment(self, k):
        '''Function that calculates the kth moment'''
        SUM = 0
        T = len(self.avalanches[self.steadstatetime:])
        for i in self.avalanches[self.steadstatetime:]:
            SUM += (i**k)/T
        return SUM