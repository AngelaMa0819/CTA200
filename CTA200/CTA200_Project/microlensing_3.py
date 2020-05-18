#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 23:58:43 2020

@author: angelama
"""


import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np
import datetime
import julian
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u

#PART 2

#Question 2
CB_c = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
                  
#assign variable of the model
t_0 = 2459031.5 #JD
t_E = 100 #days
u_0 = 0.1
alpha = 0
rho = 1/1000
s = 1.2
q = 0.01
coord = SkyCoord('18:00:00 -30:00:00',unit=(u.hourangle, u.deg))

#set up the time array for t_0-2t_E to t_0+2t_E
t_gc = julian.from_jd(t_0,'jd')
t_start = julian.to_jd(t_gc - datetime.timedelta(days=2*t_E),fmt='jd')
t_end = julian.to_jd(t_gc + datetime.timedelta(days=2*t_E),fmt='jd')
times = np.linspace(t_start,t_end,100)


#change the unit for the times array to t_E
times_plot = (times-t_0)/t_E

#define parameters
params = mm.modelparameters.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho':rho ,'alpha':alpha,'s':s,'q':q})


def magnif (q,s,rho):
    
    """
    
  Find the magnification, trajectory  and caustics curve of the binary lens system
    
    Parameters: 
        s:seperation between two lens
        q:mass ratio between two lens
        rho:source  scaled  size
    
    Output:
        light: magnification amplitude of of the binary lens system
        traj.x,traj,y: trajectory of the binary lens system in x and y direction
        X,Y: x and y coordination of points on caustics curve
        
    """
    #find mass with q
    m_1 = 1 / (q+1)
    m_2 = q / (q+1)
    
    #reset parameters
    params.q = q
    params.s = s
    params.rho = rho
    
    #find trajectory
    traj = mm.Trajectory(times,params,coords=coord)
    
    #find magnification
    light = []
    model_1 = mm.BinaryLens(m_1,m_2,s)
    
    for i in range(len(traj.x)):    
         light.append(model_1.vbbl_magnification(traj.x[i],traj.y[i],rho))
    
    #find caustics curves
    model_2 = mm.Caustics(q,s)
    X,Y = model_2.get_caustics()
    return light,traj.x,traj.y,X,Y


#Part a)
    
#call the function for different q
result_1 = magnif(0.01,1.2,0.001)
result_2 = magnif(0.1,1.2,0.001)
result_3 = magnif(1,1.2,0.001)

#plot the result
plt.plot(times_plot,2.5 * np.log10(result_1[0]),c=CB_c[1],label='q=0.01')
plt.plot(times_plot,2.5 * np.log10(result_2[0]),'b',label='q=0.10')
plt.plot(times_plot,2.5 * np.log10(result_3[0]),c='g',label='q=1.00')
plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Binary Lense')
plt.grid()
plt.legend()
plt.savefig('2_2a_l.pdf')
plt.show()



#plot the source Trajectory
plt.plot(result_1[1],result_1[2],c=CB_c[1],label='q=0.01')
plt.plot(result_2[1],result_2[2],'b-',label='q=0.10')
plt.plot(result_3[1],result_3[2],'g-',label='q=1.00')
plt.legend()
#plot the caustic wave
plt.scatter(result_1[3],result_2[4],s=0.05,color=CB_c[1])
plt.scatter(result_2[3],result_2[4],s=0.05,color='b')
plt.scatter(result_3[3],result_3[4],s=0.05,color='g')

plt.xlabel('X coordinates (X/$r_E$)')
plt.ylabel('Y coordinates (Y/$r_E$)')
plt.title('Caustic Wave for Gravitational Microlensing')
plt.grid()
plt.savefig('2_2a_c.pdf')
plt.show()




#Part b)
#call the function for different s
result_4 = magnif(0.1,0.8,0.001)
result_5 = magnif(0.1,1.0,0.001)
result_6 = magnif(0.1,1.2,0.001)

#plot the result
plt.plot(times_plot,2.5 * np.log10(result_4[0]),c=CB_c[1],label='s=0.8')
plt.plot(times_plot,2.5 * np.log10(result_5[0]),'b-',label='s=1.0')
plt.plot(times_plot,2.5 * np.log10(result_6[0]),'g-',label='s=1.2')
plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Binary Lense')
plt.grid()
plt.legend()
plt.savefig('2_2b_l.pdf')
plt.show()


#add the legend to the graph
plt.plot(result_4[3][0],result_4[4][0],c=CB_c[1],label='s=0.8')
plt.plot(result_5[3][0],result_5[4][0],'b-',label='s=1.0')
plt.plot(result_6[3][0],result_6[4][0],'g-',label='s=1.2')
plt.legend()

#plot the caustic wave
plt.scatter(result_4[3],result_4[4],s=0.05,color=CB_c[1])
plt.scatter(result_5[3],result_5[4],s=0.05,color='b')
plt.scatter(result_6[3],result_6[4],s=0.05,color='g')
plt.xlabel('X coordinates (X/$r_E$)')
plt.ylabel('Y coordinates (Y/$r_E$)')
plt.title('Caustic Wave for Gravitational Microlensing')
plt.grid()
plt.savefig('2_2b_c.pdf')
plt.show()




#Part c)
#call the function for different rho
result_7 = magnif(0.1,1.2,0.001)
result_8 = magnif(0.1,1.2,0.01)
result_9 = magnif(0.1,1.2,0.1)



#plot the result
plt.plot(times_plot,2.5 * np.log10(result_7[0]),c=CB_c[1],label='$\\rho$=0.001')
plt.plot(times_plot,2.5 * np.log10(result_8[0]),'b-',label='$\\rho$=0.001')
plt.plot(times_plot,2.5 * np.log10(result_9[0]),'g-',label='$\\rho$=0.100')
plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Binary Lense')
plt.grid()
plt.legend()
plt.savefig('2_2c_l.pdf')
plt.show()

#add the legend to the graph
plt.plot(result_7[3][0],result_7[4][0],c=CB_c[1],label='$\\rho$=0.001')
plt.plot(result_8[3][0],result_8[4][0],'b-',label='$\\rho$=0.010')
plt.plot(result_9[3][0],result_9[4][0],'g-',label='$\\rho$=0.100')
plt.legend()

#plot the caustic wave
plt.scatter(result_7[3],result_7[4],s=0.05,color=CB_c[1])
plt.scatter(result_8[3],result_8[4],s=0.05,color='b')
plt.scatter(result_9[3],result_9[4],s=0.05,color='g')
plt.xlabel('X coordinates (X/$r_E$)')
plt.ylabel('Y coordinates (Y/$r_E$)')
plt.title('Caustic Wave for Gravitational Microlensing')
plt.grid()
plt.savefig('2_2c_c.pdf')
plt.show()



    