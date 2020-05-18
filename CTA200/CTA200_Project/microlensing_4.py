#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 02:27:12 2020

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

#Question 3

#assign variable of the model
t_0 = 2459031.5 #JD
t_E = 100 #days
u_0 = 0.1
alpha = 0
rho = 1/1000
s = 1.2
q = 1/10
m_1  = 10
m_2  = 1
r_E = 10 #in the unit of AU
coord = SkyCoord('18:00:00 -30:00:00',unit=(u.hourangle, u.deg))
paral = {'earth_orbital': True, 'satellite': False, 'topocentric': False}

#define gravtational constant
G = 39.478  #AU^-3 * yr^-2 * M_solar^-1

#set up the time array for t_0-2t_E to t_0+2t_E
t_gc = julian.from_jd(t_0,'jd')
t_start = julian.to_jd(t_gc - datetime.timedelta(days=2*t_E),fmt='jd')
t_end = julian.to_jd(t_gc + datetime.timedelta(days=2*t_E),fmt='jd')
times = np.linspace(t_start,t_end,100)


#change the unit for the times array to t_E
times_plot = (times-t_0)/t_E



#Part a
#find period of the orbital motion of the binary axis
T = 2 * np.pi * np.sqrt((s * r_E) ** 3/(G * (m_1 + m_2)))
print("The period of the orbital motion is",T,"yr")




#Part b
#find rate of change of alpha in deg/year
dalpha_dt = 360/T
print(dalpha_dt)

#since he binary system is exactly face-on,ds_dt=0
ds_dt=0

#define parameters
#mm.modelparameters.which_parameters('lens orbital motion')

params_1 = mm.modelparameters.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho':rho,'alpha':alpha,'s':s,'q':q,'ds_dt':ds_dt,'dalpha_dt':dalpha_dt})
params_2 = mm.modelparameters.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'rho':rho,'alpha':alpha,'s':s,'q':q})
params_3 = mm.modelparameters.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E_N':0.2 ,'pi_E_E':0,'rho':rho,'alpha':alpha,'s':s,'q':q})
#get the trajectory of the source
traj_1 = mm.Trajectory(times,params_1,coords=coord)
traj_2 = mm.Trajectory(times,params_2,coords=coord)
traj_3 = mm.Trajectory(times,params_3,parallax=paral,coords=coord)



#get the magnification of the model
light_1 = []
light_2 = []
light_3 = []
model_1 = mm.BinaryLens(m_1/11,m_2/11,s)


    
for i in range(len(traj_1.x)):    
    light_1.append(model_1.vbbl_magnification(traj_1.x[i],traj_1.y[i],rho))
    light_2.append(model_1.vbbl_magnification(traj_2.x[i],traj_2.y[i],rho))
    light_3.append(model_1.vbbl_magnification(traj_3.x[i],traj_3.y[i],rho))
    
#plot the light curve
    
plt.plot(times_plot,2.5 * np.log10(light_2),'r--',label='Static')
plt.plot(times_plot,2.5 * np.log10(light_1),'y-',lw=2,label='Orbital Motion')
plt.plot(times_plot,2.5 * np.log10(light_3),'b--',label='Parallax')
plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Binary Lense')
plt.grid()
plt.legend()
plt.savefig('2_3b.pdf')
plt.show()




#Part c
#find caustics waves
model_2 = mm.Caustics(q,s)
X,Y = model_2.get_caustics()

#get the trajectory of the source
plt.plot(traj_1.x,traj_1.y,'y-',label='Orbital Motion')
plt.plot(traj_2.x,traj_2.y,'r-',label='Static')
plt.plot(traj_3.x,traj_3.y,'b-',label='Parallax')
plt.plot(X[0],X[0],'g-',label='Static Caustics')
plt.scatter(X,Y,s=0.05,color='g')
plt.xlabel('X coordinates (X/$r_E$)')
plt.ylabel('Y coordinates (Y/$r_E$)')
plt.legend()
plt.title('Trajectory for Gravitational Microlensing')
plt.grid()
plt.savefig('2_3c.pdf')
plt.show()
   
