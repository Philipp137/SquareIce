#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 14:27:19 2019

@author: phil
"""
import itertools
from numpy import *
from numpy.linalg import svd,norm
import numpy as np
from numpy import exp, mod,meshgrid
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# %%
def gauss_law_is_satisfied(config,boundary,Lx,Ly):
    sigma =[0,0,0,0]
    for ix in range(Lx):
        for iy in range(Ly):
            
            index = ix*Ly + iy
            link_up = index*2+1
            link_right = index*2
            sigma[0] =  config[link_up]
            sigma[1] =  config[link_right]
            
            # link  down
            if iy == 0: 
                index_down = ix *Ly+ (Ly - 1)
                link_down = index_down*2+1
            else:
                index_down = ix *Ly + (iy - 1)
                link_down = index_down*2+1
            
            
            sigma[2] = config[link_down]
            
            if ix == 0:
                sigma[3] = boundary[iy]
                
            else:
                index_left = (ix - 1)*Ly + iy
                link_left = index_left*2
                sigma[3] = config[link_left]
            
            
            charge = -sigma[0]-sigma[1]+sigma[2]+sigma[3];
            #print("ix",ix,"iy",iy,"links ^>v<",sigma,"charge",charge)
            
            if charge !=0 :
                return False

    return True

# %%
Lx = 5
Ly = 2
boundary=[-1,1]
Nlinks = 2*Lx*Ly
# %%

states=[-1,1]
list_of_lattice_config = list(itertools.product([-1,1],repeat=Nlinks))
count = 0;
list_of_gauss_invariant_states =[]
for n,config in enumerate(list_of_lattice_config):
    #print(n," Checking config:" ,config)
 
    if gauss_law_is_satisfied(config,boundary,Lx,Ly):
        count = count + 1
        list_of_gauss_invariant_states.append(config)
        
print("Number of states fulfilling the gauss law: ",count)             