#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 18:47:00 2020

@author: Xiaoyi Ma
"""

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

#Setup constants for time
t_0 = 0
t_max = 200
tpoints = np.linspace(t_0, t_max,1000 )


#Setup rate constants
beta = 0.2
gamma = 0.02

#Initial Conditions
I_0 = 1
S_0 = 999
R_0 = 0
N = 1000
y_0 = [S_0,I_0,R_0]


#Define differential equations for SIR models
def SIR(t,y):
    
    """
    
   Define the differential equation for S,I,R parameter
    
    Parameters: 
        t:time for the function,float
        y:value of S,I,R parameter,np.array
    
    Output:
        fS,fI,fR: result of differential functions in np.array
        
        
    """
    
    #Define the inital value for S,I,R,N
    S = y[0]
    I = y[1]
    R = y[2]
    
    
    #function for S,I,R
    fS = - beta * S * I / N 
    fI = beta * S * I / N - gamma * I
    fR = gamma * I 
    return ([fS, fI, fR]) 

#solve the ODE with Scipy
sol= solve_ivp(lambda t,y: SIR(t,y),[t_0,t_max],y_0,t_eval=tpoints)
S,I,R = sol.y


#plot the graph for S,I,R
plt.plot(tpoints, np.array(S, float)/N,'b-', label = 'Susceptible fraction')
plt.plot(tpoints, np.array(I, float)/N,'y-' ,label = 'Infected fraction')
plt.plot(tpoints, np.array(R, float)/N, 'r-',label = 'Recovered fraction')
plt.xlabel('Time (day)')
plt.ylabel('Population Fractions(people/day)')
plt.title('Populaton Fractions vs Time')
plt.legend()
plt.grid()
plt.show()



#Bonus: new model invole the death parameter

#define a new variable, inital death parameter as zero, the
D_0 = 0
#death rate as one percent
delta = 0.006

y_0_d = [S_0,I_0,R_0,D_0]


#Define differential equations for SIR models
def SIRD(t,y):
    
    """
    
   Define the differential equation for S,I,R,D parameter
    
    Parameters: 
        t:time for the function,float
        y:value of S,I,R parameter,np.array
    
    Output:
        fS,fI,fR,fD: result of differential functions in np.array
        
        
    """
    
    #Define the inital value for S,I,R,N
    S = y[0]
    I = y[1]
    R = y[2]
    D = y[3]
    
    
    #function for S,I,R
    fS = - beta * S * I / N 
    fI = beta * S * I / N - gamma * I - delta * I
    fR = gamma * I 
    fD = delta * I
    return ([fS, fI, fR, fD]) 

#solve the ODE with Scipy
sol= solve_ivp(lambda t,y: SIRD(t,y),[t_0,t_max],y_0_d,t_eval=tpoints)
S,I,R,D = sol.y


#plot the graph for S,I,R
plt.plot(tpoints, np.array(S, float)/N,'b-', label = 'Susceptible fraction')
plt.plot(tpoints, np.array(I, float)/N,'y-' ,label = 'Infected fraction')
plt.plot(tpoints, np.array(R, float)/N, 'r-',label = 'Recovered fraction')
plt.plot(tpoints, np.array(D, float)/N, 'k-',label = 'Death fraction')
plt.xlabel('Time (day)')
plt.ylabel('Population Fractions(people/day)')
plt.title('Populaton Fractions vs Time')
plt.legend()
plt.grid()
plt.show()
