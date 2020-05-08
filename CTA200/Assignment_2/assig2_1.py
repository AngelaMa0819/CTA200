# -*- coding: utf-8 -*-
"""
Created on Wes May 6 15:34:00 2020

@author: Xiaoyi Ma
"""

import numpy as np
import matplotlib.pyplot as plt




def draw(x_i,x_f,y_i,y_f):
    
    """
    
    Draw the color code graph for the ability of convergence
    of sequence that made of c
    
    Parameters: 
        the start and end value for x and y
    
    Output: 
        the graph for convergence
        
    """
    
    #setup the array for x,y and c
    x = np.linspace(x_i,x_f,100)
    y = np.linspace(y_i,y_f,100)
   
    #setup the array for bound and diverge c
    
    
    where_diverge = np.zeros((100,100))
    
    #asign values to the c array
    for i in range(len(x)):
        for j in range(len(y)):
            c = x[i] + y[j] * 1j
            
            
            #setup the array for z and its magnitude
            z = np.zeros(100,dtype=complex)
            
            
            
            #asign value to the z array
            for k in range(len(z)-1):
                
                if np.isinf(np.abs(z[k])**2) == True:
                            
                        where_diverge[j,i]=k
                        
                        break
                else:
        
                    z[k+1] = np.abs(z[k]) ** 2 + c
        
            
            #exam where the array diverge
            if where_diverge[j,i] == 0:
                
                where_diverge[j,i]=len(z)
            
    #plot the result with a color bar for all points
    
    X,Y = np.meshgrid(x,y) 
    plt.pcolormesh(X,Y,where_diverge) 
    plt.colorbar()
    plt.xlabel("Real Axis")
    plt.ylabel("Imagery Axis")
    plt.title("Level of Divergence of the Magnitude")
    plt.show()
    
     
   
#call the function for different ranges
draw(-2,2,-2,2)
draw(-1,1,-1,1)
draw(-0.5,0.5,-0.5,0.5)
draw(0.225,0.3,-0.1,0.1)
    

        
        