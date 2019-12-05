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
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
    return p
# %% 
def draw_lattice(config,Lx,Ly,boundary_left):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    
    ax.set_xlim([-1,Lx])
    ax.set_ylim([-1,Ly])
    ax.axis('off')
    ax.set_title('$L_x ='+str(Lx)+ ', L_y='+str(Ly)+'$')
    
    # draw incoming arrows
    for iy in range(Ly):
        if boundary_left[iy]==1:
                arrow_bound = plt.arrow(0-1,iy,0.9,0, head_width=0.08, 
                                        length_includes_head=True, 
                                        alpha =0.3,linestyle='dashed') 
        else:
                arrow_bound = plt.arrow(0,iy,-0.9,0, linestyle='dashed',
                                        head_width=0.08, length_includes_head=True, 
                                        alpha = 0.3)     
        ax.add_artist(arrow_bound)
    
    
    for ix in range(Lx):
        for iy in range(Ly):
            
            node_id = ix*Ly + iy
            circle = plt.Circle((ix,iy), 0.1, color='r')
            
            link_up_id = node_id*2+1
            if config[link_up_id]==1:
                arrow_up = plt.arrow(ix,iy,0,0.9,
                                     head_width=0.08, length_includes_head=True) 
            else:
                arrow_up = plt.arrow(ix,iy+1,0,-0.9,
                                     head_width=0.08, length_includes_head=True) 
            
            link_right_id = node_id*2
            
            if ix < Lx - 1:
                if config[link_right_id]==1:
                    arrow_right = plt.arrow(ix+0.1,iy,0.8,0,
                                        head_width=0.08, length_includes_head=True) 
                else:
                        arrow_right = plt.arrow(ix+1,iy,-0.9,0, 
                                        head_width=0.08, length_includes_head=True)     
            else:
                if config[link_right_id]==1:
                    arrow_bound = plt.arrow(ix+0.1,iy,0.9,0, alpha = 0.3, 
                                            linestyle='dashed',
                                            head_width=0.08, length_includes_head=True) 
                else:
                    arrow_bound = plt.arrow(ix+1,iy,-0.9,0, alpha =0.3,
                                                linestyle='dashed',
                                                head_width=0.08, length_includes_head=True)     
                ax.add_artist(arrow_bound)
        # copy periodic links at iy = 0 from iy = Ly - 1
            if iy == 0:
                # copy link of the upper bound
                node_id = ix*Ly + Ly - 1   
                link_up_id = node_id*2+1
                if config[link_up_id]==1:
                     arrow_down = plt.arrow(ix,iy-1,0,0.9, alpha = 0.3, linestyle='dashed',
                                     head_width=0.08, length_includes_head=True) 
                else:
                    arrow_down = plt.arrow(ix,iy-0.1,0,-0.9,alpha = 0.3, linestyle='dashed',
                                     head_width=0.08, length_includes_head=True) 
                ax.add_artist(arrow_down)
             
            if ix == Lx - 1:
                  arrow_right.set_alpha=0.3
                
            ax.add_artist(arrow_up)
            ax.add_artist(arrow_right)
            ax.add_artist(circle)
    
    plt.legend([arrow_bound,arrow_right], ['boundary link','active link'],
               bbox_to_anchor=(1.15, 0.2),
               handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),})
    
    
# %%
def gauss_law_is_satisfied(config,boundary_left,boundary_right,Lx,Ly):
    sigma =[0,0,0,0]
    
    for iy in range(Ly):
        ix = Lx - 1
        node_id = ix*Ly + iy
        linkindex_right = node_id*2
        if config[linkindex_right] != boundary_right[iy]:
            return False
        
    for ix in range(Lx):
        for iy in range(Ly):
            
            node_id = ix*Ly + iy
            linkindex_up = node_id*2+1
            linkindex_right = node_id*2
            sigma[0] =  config[linkindex_up]
            sigma[1] =  config[linkindex_right]
            
            # link  down
            if iy == 0: 
                node_id_down = ix *Ly+ (Ly - 1)
                linkindex_down = node_id_down*2+1
            else:
                node_id_down = ix *Ly + (iy - 1)
                linkindex_down = node_id_down*2+1
            
            
            sigma[2] = config[linkindex_down]
            
            if ix == 0:
                sigma[3] = boundary_left[iy]
                
            else:
                node_id_left = (ix - 1)*Ly + iy
                linkindex_left = node_id_left*2
                sigma[3] = config[linkindex_left]
            
            
            charge = -sigma[0]-sigma[1]+sigma[2]+sigma[3];
            #print("ix",ix,"iy",iy,"links ^>v<",sigma,"charge",charge)
            
            if charge !=0 :
                return False
            

    return True

# %%
Lx =2
Ly = 2
boundary=[-1, 1,-1] # spins in vertical direction of boundary at ix = 0 and ix = Lx
Nlinks = 2*Lx*Ly
# %%

states=[-1,1]
list_of_lattice_config = list(itertools.product(states,repeat=Nlinks))


# %%
count = 0;
list_of_gauss_invariant_states =[]
for n,config in enumerate(list_of_lattice_config):
    #print(n," Checking config:" ,config)
 
    if gauss_law_is_satisfied(config,boundary,boundary,Lx,Ly):
        count = count + 1
        list_of_gauss_invariant_states.append(config)
        draw_lattice(config,Lx,Ly,boundary)

        
print(list_of_lattice_config)

print("Number of states fulfilling the gauss law: ",count)             