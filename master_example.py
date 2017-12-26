"""
Created on Mon Feb 13 20:25:03 2017

@author: Jonah Haefner

Title: master_example.py

The purpose of this script is to ensure the four functions included in this
package are functioning properly
dP_1 = 27.12 kPa   - This needs to be investigated. I think it is due to the angle of attack.
dP 2 = 25.67 kPa
Nu1 = 142.52
Nu2 = 147.31
"""

import TORCHE

d = 0.0254		# Outside diameter of tube or cylinder (m)
a = 1.25 	# transverse pitch to diameter ratio
b = 1.25	# longitudinal pitch to diameter ratio


rho = 1940 	# Density of the working fluid - FLiBe (kg/m^3)

dP_1 = TORCHE.dP_Zu(a,b,.5,1940,.0056,10,.0254,22000,'inline')
print('The Pressure Drop calculated by Zukauskas is',dP_1/1000,'kPa')
dP_2 = TORCHE.dP_GG(a,b,.5,1940,22000,10,'inline')
print('The Pressure Drop calculated by Gaddis-Gnielinski is',dP_2/1000,'kPa')
Nu_1 = TORCHE.HT_Zu(a,b,d,.5,1940,10,1,1,1.787e-6,'inline',22000)
print('The Nusselt Number calculated by Zukauskas is', Nu_1)
Nu_2 = TORCHE.HT_GG(a,b,d,.5,1940,22000,10,'inline',1)
print('The Nusselt Number calculated by Gnielinski is', Nu_2)