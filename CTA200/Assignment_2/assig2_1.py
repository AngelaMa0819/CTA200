#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 16:03:51 2020

@author: xma
"""


import numpy as np
import matplotlib.pyplot as plt

def draw (x_i,x_f,y_i,y_f):
    """
    
    Draw the graph for the c according to their divergence index
    
    Parameters: 
        x_i: start value for x range
        x_f: end value for x range
        y_i: start value for y range
        y_f: end value for y range
    
    Output:
       graph for color density according to divergence index
       where converge point c is denote as white color
        
        
     """
    
    #set up the array for x,y
    x = np.linspace(x_i,x_f,100)
    y = np.linspace(y_i,y_f,100)
    
    
    #setup array for the diverge index
    where_diverge = np.zeros((100,100))
    bound_x = []
    bound_y = []
    
    #assign values to c array
    for i in range(len(x)):
        for j in range(len(y)):
            c = x[i] + y[j] * 1j
            
            
            #setup the array for z 
            z = np.zeros(100,dtype=complex)
            
            #assign value to z and find diverge index
            for k in range (len(z)-1):
                
                if np.isinf(np.abs(z[k])**2) == True:
                    where_diverge[j,i] = k
                    break
                else:
                    z[k+1] = np.abs(z[k]) ** 2 + c
           
            
            #assign value diverge index for converge ones as 100
            if where_diverge[j,i] == 0:
                where_diverge[j,i] = len(z)
                bound_x.append(x[i])
                bound_y.append(y[j])
    
    
                
    #plot the result with a color bar        
    X,Y = np.meshgrid(x,y)
    plt.pcolormesh(X,Y,where_diverge)
    plt.colorbar()
    plt.scatter(bound_x,bound_y,color='w')
    plt.xlabel("Real Axis")
    plt.ylabel("Imagery Axis")
    plt.title("Level of Divergence")
    plt.show()
    
    
#call the function for different range
draw(-2,2,-2,2)
draw(-1,1,-1,1)
draw(-0.5,0.5,-0.5,0.5)
draw(0.225,0.3,-0.1,0.1)
                
                