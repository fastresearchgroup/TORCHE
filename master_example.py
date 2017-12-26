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

import TORCHE

# Geometric parameters
d = 0.0254		# Outside diameter of tube or cylinder (m)
a = 1.25 		# Transverse pitch to diameter ratio
b = 1.25		# Longitudinal pitch to diameter ratio

# Fluid thermo-physical properties
rho = 1940 		# Density of the working fluid - FLiBe (kg/m^3)

# Flow behavior
vel	= 0.5		# Free-stream velocity before interacting with the tube bank (m/s)


dP_1 = TORCHE.dP_Zu(rho,a,b,vel,.0056,10,22000,'inline')
print('The Pressure Drop calculated by Zukauskas is',dP_1/1000,'kPa')
dP_2 = TORCHE.dP_GG(rho,a,b,vel,22000,10,'inline')
print('The Pressure Drop calculated by Gaddis-Gnielinski is',dP_2/1000,'kPa')
Nu_1 = TORCHE.HT_Zu(rho,a,b,d,vel,10,1,1,1.787e-6,'inline',22000)
print('The Nusselt Number calculated by Zukauskas is', Nu_1)
Nu_2 = TORCHE.HT_GG(rho,a,b,d,vel,22000,10,'inline',1)
print('The Nusselt Number calculated by Gnielinski is', Nu_2)