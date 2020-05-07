# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import math



def draw(x_i,x_f,y_i,y_f):
    
    #setup the array for x,y and c
    x = np.linspace(x_i,x_f,100)
    y = np.linspace(y_i,y_f,100)
   
    #setup the array for bound and diverge c
    
    c_bound = []
    c_diverge = []
    where_diverge = np.zeros((100,100))
    
    #asign values to the c array
    for i in range(len(x)):
        for j in range(len(y)):
            c = x[i] + y[j] * 1j
            
           
            #setup the array for z and its magnitude
            z = np.zeros(100,dtype=complex)
            z_mag = np.zeros(100,dtype=complex)
    
        
            #asign value to the z array
            for k in range(len(z)-1):
        
                z[k+1] = (np.real(z[k]) ** 2 + np.imag(z[k]) ** 2) + c
        
            #asign value to the z_mag array
            for k in range(len(z)):
            
                z_mag[k] = np.real(z[k]) ** 2 + np.imag(z[k]) ** 2
          
            
            #exam where the array diverge
            if math.isfinite(z_mag[len(z)-1])== True:
                
                c_bound.append(c)
                where_diverge[j,i]=len(z_mag)
                
            
            else:
                c_diverge.append(c)
                for k in range(len(z_mag)):
                    if np.isinf(z[k]) == True:
                        
                        where_diverge[j,i]=k
                        
                        break
    
    
    
    #plot the result with a color bar for all points
    
    X,Y = np.meshgrid(x,y) 
    plt.pcolormesh(X,Y,where_diverge) 
    plt.colorbar()
    plt.xlabel("Real Axis")
    plt.ylabel("Imagery Axis")
    plt.title("Level of Divergence of the Magnitude")
    plt.show()
    
    
            
        
        
        
        
   
#call the function for different range
draw(-2,2,-2,2)
draw(-1,1,-1,1)
draw(-0.5,0.5,-0.5,0.5)
draw(0.225,0.3,-0.1,0.1)
    

        
        