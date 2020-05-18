#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:17:55 2020

@author: angelama

Part 1 of CTA-200H Computing Course Project
"""

import MulensModel as mm
import matplotlib.pyplot as plt
import numpy as np
import datetime
import julian
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy import units as u

#PART 1 

#Question 2



#assign variable of the model
t_0 = 2459031.5 #JD
t_E = 100 #days
u_0 = 0.1
params = mm.modelparameters.ModelParameters({'t_0': t_0, 'u_0': u_0, 't_E': t_E, 'pi_E_N':0 ,'pi_E_E':0})



#set up the time array for t_0-2t_E to t_0+2t_E
t_gc = julian.from_jd(t_0,'jd')
t_start = julian.to_jd(t_gc - datetime.timedelta(days=2*t_E),fmt='jd')
t_end = julian.to_jd(t_gc + datetime.timedelta(days=2*t_E),fmt='jd')
times = np.linspace(t_start,t_end,100)


#change the unit for the times array to t_E
times_plot = (times-t_0)/t_E
    


#set up PSPL model 
#model_PSPL= mm.MagnificationCurve(times,params,coords='18:00:00 -30:00:00')

model_PSPL_1= mm.MagnificationCurve(times,params,coords='18:00:00 -30:00:00')
params.u_0 = 0.3
model_PSPL_2= mm.MagnificationCurve(times,params,coords='18:00:00 -30:00:00')
params.u_0 = 1.0
model_PSPL_3= mm.MagnificationCurve(times,params,coords='18:00:00 -30:00:00')



#get light curve data for different u_0
curve_1 = model_PSPL_1.get_point_lens_magnification()
curve_2 = model_PSPL_2.get_point_lens_magnification()
curve_3 = model_PSPL_3.get_point_lens_magnification()
"""
model_PSPL.parameters.u_0 = 0.3
print(model_PSPL.parameters)
curve_2 = model_PSPL.get_point_lens_magnification()


model_PSPL.parameters.u_0 = 1.0
print(model_PSPL.parameters)
curve_3 = model_PSPL.get_point_lens_magnification()
"""



#plot the light curve
plt.plot(times_plot,2.5 * np.log10(curve_1),'b-',label='$u_0$=0.1')
plt.plot(times_plot,2.5 * np.log10(curve_2),'g-',label='$u_0$=0.3')
plt.plot(times_plot,2.5 * np.log10(curve_3),'r-',label='$u_0$=1.0')

plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Point-source Lense Without Parallax')
plt.grid()
plt.legend()
plt.savefig('1.2.pdf')
plt.show()




#Question 3


#assign new values to the parameters
params.u_0 = 0.3
params.pi_E_E= 0.5

#set up parallax and coordinates of the event
paral = {'earth_orbital': True, 'satellite': False, 'topocentric': False}
coord = SkyCoord('18:00:00 -30:00:00',unit=(u.hourangle, u.deg))


#set up PSPL model and get light curve data for different pi_E
model_PSPL_4= mm.MagnificationCurve(times,params,parallax=paral,coords=coord)
curve_4 = model_PSPL_4.get_point_lens_magnification() 

#change the value for pi_E
params.pi_E_N= 0.5
params.pi_E_E= 0

#get new model and data for new parameters
model_PSPL_5= mm.MagnificationCurve(times,params,parallax=paral,coords=coord)
curve_5 = model_PSPL_5.get_point_lens_magnification() 

#plot the curve

plt.plot(times_plot,2.5 * np.log10(curve_2),'g--',label='Non-parallax')


plt.plot(times_plot,2.5 * np.log10(curve_4),'b-',label='$\pi_E$=[0,0.5]')
plt.plot(times_plot,2.5 * np.log10(curve_5),'r-',label='$\pi_E$=[0.5,0]')

plt.xlabel('Time ($(t-t_0)/t_E$)')
plt.ylabel('Magnification in $2.5log_{10}$ scale')
plt.title('Light Curve of Point-source Lense')
plt.grid()
plt.legend()
plt.savefig('1.3.1.pdf')
plt.show()




#find the trajectory
traj_1 = mm.Trajectory(times,model_PSPL_2.parameters,coords=coord)
traj_2 = mm.Trajectory(times,model_PSPL_2.parameters,coords=coord)
traj_3 = mm.Trajectory(times,model_PSPL_3.parameters,coords=coord)
traj_4 = mm.Trajectory(times,model_PSPL_4.parameters,parallax=paral,coords=coord)
traj_5 = mm.Trajectory(times,model_PSPL_5.parameters,parallax=paral,coords=coord)


#plot the trajectory
plt.plot(traj_2.x,traj_2.y,'g--',label='Non-parallax')


plt.plot(traj_4.x,traj_4.y,'b-',label='$\pi_E$=[0,0.5]')
plt.plot(traj_5.x,traj_5.y,'r-',label='$\pi_E$=[0.5,0]')

plt.xlabel('X coordinates of the Trajectory')
plt.ylabel('Y coordinates of the Trajectory')
plt.title('Source Trajectory on the Source Plane')
plt.grid()
plt.legend()
plt.savefig('1.3.2.pdf')
plt.show()

