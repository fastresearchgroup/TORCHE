# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 20:25:03 2017

@author: Jonah Haefner

Title: master_example.py

The purpose of this script is to ensure the four functions included in this
package are functioning properly
dP_1 = 27.12 kPa  
dP 2 = 25.67 kPa
Nu1 = 142.52
Nu2 = 147.31
"""

import TORCHE
dP_1 = TORCHE.dP_Zu(1.25,1.25,.5,1940,.0056,10,.0254,22000,'inline')
print('The Pressure Drop calculated by Zukauskas is',dP_1/1000,'kPa')
dP_2 = TORCHE.dP_GG(1.25,1.25,.5,1940,22000,10,'inline')
print('The Pressure Drop calculated by Gaddis-Gnielinski is',dP_2/1000,'kPa')
Nu_1 = TORCHE.HT_Zu(1.25,1.25,.5,1940,10,1,1,.0254,1.787e-6,'inline',22000)
print('The Nusselt Number calculated by Zukauskas is', Nu_1)
Nu_2 = TORCHE.HT_GG(1.25,1.25,.0254,.5,1940,22000,10,'inline',1)
print('The Nusselt Number calculated by Gnielinski is', Nu_2)