#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 22:26:43 2020

@author: angelama
"""

import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np
from sympy import *


#PART 1 

#Question 1


def constant(q,s):
    
    """
    
   Find varaibles in the critical curve equation with s,q
    
    Parameters: 
        s:seperation between two lens
        q:mass ratio between two lens
    
    Output:
        z_1c,z_2c:postion of two lens
        esp_1,esp_2:mass ratio of two lens and the total mass
        
    """
    
    #define variable in the polynomial with given parameter
    esp_1 = 1 / (q+1)
    esp_2 = q / (q+1)
    z_1c = s * q / (1+q)
    z_2c = -s / (q+1)
    
    return  esp_1,esp_2,z_1c,z_2c
    
   

def coeff(esp_1,esp_2,z_1c,z_2c,phi):
    
    """
    
      Find coefficients of the critical curve equation 
    
    Parameters: 
        z_1c,z_2c:postion of two lens
        esp_1,esp_2:mass ratio of two lens and the total mass
    
    Output:
        coeffs:coefficients of the critical curve equation 
        
        
    """
    
    #define the conjugate of z
    z_c = symbols('z_c')
    
    #define the numerator and denominator on LHS
    A = esp_2 * (z_c - z_1c) ** 2 + esp_1 * (z_c - z_2c) ** 2
    B = ((z_c - z_2c) ** 2) * ((z_c - z_1c) ** 2)
    
    func = Poly( A-np.exp(1j * phi) * B, z_c)
    
    coeffs = np.array(func.coeffs())
    
    return coeffs.astype(complex)

#find the coefficient and roots of the equation for values of phi 
phi = np.linspace(0,2*np.pi,100)
z = []

#find constants first since it is independent of phi
esp_1,esp_2,z_1c,z_2c = constant(0.2,1.2)

for i in range(len(phi)):
  
    coeffs = coeff(esp_1,esp_2,z_1c,z_2c,phi[i])
    
    roots = np.roots(coeffs)
    z.append(np.conjugate(roots))
    
z_lin = np.reshape(z,len(z)*4)


#find the position on source plane
zeta = z_lin - esp_1 / (np.conjugate(z_lin)-z_1c) - esp_2 / (np.conjugate(z_lin)-z_2c)


#find the casutic wave with Mulens Model

model = mm.Caustics(0.2,1.2)#input as (q,s)
X,Y = model.get_caustics(n_points=400)

#plot the caustic wave
plt.scatter(-np.real(zeta),np.imag(zeta),label='Python')
#inverse of the x value since the mulensmodel take the heavy mass on the positive axis
plt.scatter(X,Y,label='Mulens')
plt.xlabel('X coordinates (X/$r_E$)')
plt.ylabel('Y coordinates (Y/$r_E$)')
plt.title('Caustic Wave for Gravitational Microlensing')
plt.grid()
plt.legend()
plt.savefig('2_1.pdf')
plt.show()

#From the graph, we can see that shape of the caustics curve
#from the code the MulensModel is the same




   


    
    


    
    
