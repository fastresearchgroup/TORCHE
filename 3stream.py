'''
Author: Lane Carasik
Purpose: Three-stream heat exchanger temperature distribution (1-D)
'''

#----------------------------------------------------------------------------------#
# import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#----------------------------------------------------------------------------------#
# Define functions needed

def Prandtl(c_p,mu,k):
	# Prandtl number calculation
	Pr = c_p*mu/k
	
	return Pr
    
def Reynolds(rho,vel,D_h,mu):
	# Reynolds number calculation
	Re = rho*vel*D_h/mu
	
	return Re

#----------------------------------------------------------------------------------#
# Options/Inputs

# Solver Settings
nodes = 100  # Number of nodes for 1-D calculation

# Input Conditions
# Temperatures
T_1in = 30 # Inlet Temperature (Stream 1)
T_2in = 50 # Inlet Temperature (Stream 2)
T_3in = 80 # Inlet Temperature (Stream 3)
# Flow rates
V_flow_1 = 1 # Flow rate (Stream 1) 
V_flow_2 = 5 # Flow rate (Stream 2) 
V_flow_3 = 3 # Flow rate (Stream 3) 

# Call Fluid properties
'''
Currently assuming liquids and not gases.
'''


#----------------------------------------------------------------------------------#
# Basic Global Calulations

# Call mass flow rate calculation 

# Call area calculation
A_1 = .010 # Flow area of stream 1
A_2 = .015 # Flow area of stream 2
A_3 = .025 # Flow area of stream 3

# Call Hydraulic diameter calculation

# Call Reynolds number calculation

# Call Prandtl number calcation
Pr_1 = Prandtl(c_p,mu,k)
Pr_2 = Prandtl(c_p,mu,k)
Pr_2 = Prandtl(c_p,mu,k)
# Call Nusselt number calculation
# Call HTC Calculation

#----------------------------------------------------------------------------------#
# Boundary Value Problem Calculation

#----------------------------------------------------------------------------------#
# Output Formatted Text

#----------------------------------------------------------------------------------#
# Output relevant plots