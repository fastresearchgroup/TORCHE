"""
Created on Mon Feb 13 20:25:03 2017

@author: Jonah Haefner

Title: master_example.py

The purpose of this script is to ensure the four functions included in this
package are functioning properly
dP_1 = 21.94 kPa
dP 2 = 25.67 kPa
Nu1 = 142.52
Nu2 = 147.31
"""

import TORCHE as TE

# Geometric parameters
d = 0.0254		# Outside diameter of tube or cylinder (m)
a = 1.25 		# Transverse pitch to diameter ratio
b = 1.25		# Longitudinal pitch to diameter ratio
geom = 'inline' # Tube layout 
N_rows = 10 	# Number of tube rows

# Fluid thermo-physical properties
rho = 1940 		# Density of the working fluid - FLiBe salt (kg/m^3)
mu = 0.0056 	# Dynamic visocity of the working fluid - FLiBe salt (Pa-s)
Pr = 1			# Prandtl number of the working fluid
Pr_w = 1 		# Prandtl number of the working fluid based on the wall film temperature

# Flow behavior
vel	= 0.5		# Free-stream velocity before interacting with the tube bank (m/s)
Re = 22000		# Reynolds number of the flow based on the maximium velocity in the minimum area between tubes

dP_1 = TE.dP_Zu(rho,a,b,geom,N_rows,vel,Re)
print('The Pressure Drop calculated by Zukauskas is',dP_1/1000,'kPa')
dP_2 = TE.dP_GG(rho,a,b,geom,N_rows,vel,Re,Return="")
print('The Pressure Drop calculated by Gaddis-Gnielinski is',dP_2/1000,'kPa')
Nu_1 = TE.HT_Zu(rho,Pr,Pr_w,a,b,d,geom,N_rows,vel,Re)
print('The Nusselt Number calculated by Zukauskas is', Nu_1)
Nu_2 = TE.HT_GG(rho,Pr,a,b,d,geom,N_rows,vel,Re)
print('The Nusselt Number calculated by Gnielinski is', Nu_2)