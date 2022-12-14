# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 09:21:41 2022

@author: Alexander
"""
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from logbin_2020 import logbin

class Network:
    def __init__(self, rand=True, N=random.choice(range(6, 10)), connections=[], E=random.choice(range(1, 7))):
        self.t = 0
        #self.edgemat = [[0 for i in range(N)] for j in range(N)]
        self.adjacencylist = [[] for i in range(N)]
        if rand:
            while len(connections) < E:
                connection = random.sample(range(N), 2)
                if connection not in connections and connection[::-1] not in connections:
                    connections.append(connection)
        else: 
            E = len(connections)
        print(connections)
        for i in range(E):
            #self.edgemat[connections[i][0]][connections[i][1]] += 1
            #self.edgemat[connections[i][1]][connections[i][0]] += 1
            self.adjacencylist[connections[i][0]].append(connections[i][1])
            self.adjacencylist[connections[i][1]].append(connections[i][0])
        for i in self.adjacencylist:
            i.sort()
        
    def drive(self, N, m, phase = 'P1'):
        while len(self.adjacencylist) < N:
            self.t += 1
            #for i in self.edgemat:
             #   i.append(0)
            #self.edgemat.append([0 for i in range(len(self.edgemat) + 1)])
            self.adjacencylist.append([])
            if phase == 'P1':
                self.__add_edges(self.__probabilities('PA'), m)
            elif phase == 'P2': 
                self.__add_edges(self.__probabilities('Random'), m)
            elif phase == 'P3':
                r = m//2
                self.__add_edges(self.__probabilities('Random'), r)
                self.__add_edges(self.__probabilities('PA'), m-r, existing=True)
        #self.status()
        
    def __add_edges(self, probabilities, m, existing=False):
        options = range(len(self.adjacencylist))
        if existing:
            i = sum([len(j) for j in self.adjacencylist])
            M = i + m
            while i < M:
                connection = random.sample(range(len(self.adjacencylist)), 2) 
                if connection[0] not in self.adjacencylist[connection[1]]:
                    self.adjacencylist[connection[0]].append(connection[1])
                    self.adjacencylist[connection[1]].append(connection[0])  
                i = sum([len(j) for j in self.adjacencylist])   
        else:
            choices = np.random.choice(options, m, replace=False, p=probabilities)
            choices2 = [len(self.adjacencylist)-1 for i in range(m)]
            for i in range(m):
            #self.edgemat[choices[i]][-1] += 1
            #self.edgemat[-1][choices[i]] += 1
                self.adjacencylist[choices[i]].append(choices2[i])
                self.adjacencylist[choices2[i]].append(choices[i])    
            
    def __probabilities(self, type='PA'):
        probabilitieslist = []
        #probabilities = []
        if type == 'PA':
            for i in range(len(self.adjacencylist)):
                probabilitieslist.append(len(self.adjacencylist[i]))
        elif type == 'Random':
            probabilitieslist = [1/len(self.adjacencylist) for i in range(len(self.adjacencylist))]
        else:
            print('Invalid probability type')
            #probabilities.append(sum([self.edgemat[j][i] for j in range(len(self.edgemat))]))
        #if probabilitieslist == probabilities:
        return np.array(probabilitieslist)/sum(probabilitieslist)
        #else:
         #   print('Probabilities unequal, debug required')
    
    def degreeprobabilities(self):
        N = len(self.adjacencylist)
        dist = [len(i) for i in self.adjacencylist]
        ks, distribution = logbin(dist)
        return ks, distribution, N
            
    
    def status(self):
        print('Number of Nodes:', len(self.adjacencylist))
        print('Adjacency list:', self.adjacencylist)
        #print('Edge matrix:',  self.edgemat)
        self.visualise()
    
    def visualise(self):
        plt.axis('off')
        R = len(self.adjacencylist)
        angles = np.linspace(0, 2*np.pi, R, endpoint=False)
        locations = [R*np.cos(angles), R*np.sin(angles)]
        plt.plot(locations[0], locations[1], '.', ms=8, c='r' )
        for i in range(len(self.adjacencylist)):
            plt.text(locations[0][i], locations[1][i], str(i), color='Black', fontsize=12)
            for j in self.adjacencylist[i]:
                linex = [locations[0][i], locations[0][j]]
                liney = [locations[1][i], locations[1][j]]
                plt.plot(linex, liney, c='black')
        
        plt.show()
    
    