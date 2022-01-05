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
# Define or import material properties
    
def thermophys_Na(T):
    # Where T is in Kelvin
    Cp = -3.001*(10**6)*T**-2 + 1658 - 0.8479*T + 4.454*(10**-4)*T**2
    k = 104+0.047*T
    ln_mu = 556.835/T - 0.3958*np.log(T) - 6.4406
    mu = np.exp(ln_mu)
    rho = 1014-0.235*T
    return rho,Cp,mu,k
 
def thermophys_FLiBe(T):
    # Where T is in Kelvin
    rho = 2413-0.488*(T); #kg/m^3, error +- 2%
    Cp = 2386
    k = 1.1
    mu = (0.116*np.exp(3755/T))*0.001
    return rho,Cp,mu,k
 
def thermophys_H20(T):
    # Where T is in Kelvin
    rho = 713.8
    Cp = 5750
    k = 0.548
    mu = 0.086e-3
    return rho,Cp,mu,k

#----------------------------------------------------------------------------------#
# Options/Inputs

# Solver Settings
nodes = 100  # Number of nodes for 1-D calculation

# Input Conditions
# Temperatures
T_1in = 30 # Inlet Temperature (Stream 1) - C
T_2in = 50 # Inlet Temperature (Stream 2) - C
T_3in = 80 # Inlet Temperature (Stream 3) - C
# Flow rates
V_flow_1 = 1 # Flow rate (Stream 1) - m3/sec
V_flow_2 = 5 # Flow rate (Stream 2) - m3/sec
V_flow_3 = 3 # Flow rate (Stream 3) - m3/sec

# Call Fluid properties
'''
Currently assuming liquids and not gases.
'''
[rho_1,cp_1,mu_1,k_1] = thermophys_H20(T_1in+273.15) # Thermophysical (Stream 1)
[rho_2,cp_2,mu_2,k_2] = thermophys_H20(T_2in+273.15) # Thermophysical (Stream 2)
[rho_3,cp_3,mu_3,k_3] = thermophys_H20(T_3in+273.15) # Thermophysical (Stream 3)

#----------------------------------------------------------------------------------#
# Basic Global Calulations

# Call mass flow rate calculation (mdot = rho*area*vel = rho*vol_flow)
mdot_1 = rho_1*V_flow_1 # Mass Flow Rate of Stream 1 - kg/sec
mdot_2 = rho_2*V_flow_2 # Mass Flow Rate of Stream 2 - kg/sec
mdot_3 = rho_3*V_flow_3 # Mass Flow Rate of Stream 3 - kg/sec

# Call heat capacity rates (C = mdot*cp)
C_1 = mdot_1*cp_1 # Heat Capacity of Stream 1
C_2 = mdot_2*cp_2 # Heat Capacity of Stream 1
C_3 = mdot_3*cp_3 # Heat Capacity of Stream 1

# Call area calculation
A_1 = .010 # Flow area of stream 1
A_2 = .015 # Flow area of stream 2
A_3 = .025 # Flow area of stream 3

# Call Hydraulic diameter calculation

# Call Reynolds number calculation

# Call Prandtl number calcation
Pr_1 = Prandtl(cp_1,mu_1,k_1) # Prandtl (Stream 1)
Pr_2 = Prandtl(cp_2,mu_2,k_2) # Prandtl (Stream 2)
Pr_2 = Prandtl(cp_3,mu_3,k_3) # Prandtl (Stream 3)
# Call Nusselt number calculation
# Call HTC Calculation

#----------------------------------------------------------------------------------#
# Boundary Value Problem Calculation

#----------------------------------------------------------------------------------#
# Output Formatted Text

#----------------------------------------------------------------------------------#
# Output relevant plots