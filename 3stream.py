'''
Author: Lane Carasik
Purpose: Three-stream heat exchanger temperature distribution (1-D)
'''

#----------------------------------------------------------------------------------#
# import modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from iapws import IAPWS95

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

def f_Blasius(Re):
    # Friction factor for calculating smooth pipe 
    # friction factor
    f = 0.3164*Re**(-0.25)
    
    return f

def Nu_DB(Re,Pr):
    # Nusselt number correlation (Dittus Boelter)
    # for smooth pipe flow
    # cooling so n = 0.3 or n = 0.4 for heating
    Nu = 0.023*(Re**0.8)*(Pr**0.4)

    return Nu   

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
    # Need to add the ability to change temperature and pressure
    # Likely need to use the IAPWS to make this easy.
    rho = 713.8
    Cp = 5750
    k = 0.548
    mu = 0.086e-3
    return rho,Cp,mu,k

def NaNO2_60wt_KNO2_40wt(T):
	'''
	Temp range: 300-600 C (573.16-873.16 K)
	'''
	rho = 1000*(2.0772-0.0006*T) # Temp in C - results in kg/m3
	mu = 0.001*5.103*np.exp(-2.575*(10**3)*(T+273.15)**(-1)+1.305*(10**6)*((T+273.15)**(-2))) # Temp in K - results in Pa-s
	k = 0.459 + T*0.0003 # Temp in C - results in W/m-K
	Cp = 1570 # result in J/kg-K

	return rho,mu,k,Cp

def NaNO3_NaNO2_KNO3(T):
	'''
	HITEC salt

	Temp range: 174.84-499.84 C (448-773 K)

	Ref: Serrano-Lopez 2013 - Molten salts database for energy applications
	Density SAM-NREL 2012 (0.07–0.49–0.44)
	'''
	rho = 2279.799-0.7324*(T+273.15) # Temp in K - results in kg/m3
	mu = np.exp(-4.343 - 2.0143*(np.log(T + 273.15-273) - 5.011)) # Temp in K - results in Pa-s
	k = 0.45 # results in W/m-K
	Cp = 1560 # result in J/kg-K

	return rho,mu,k,Cp

def thermophys_H20_iapws95(T,P):

    T_inK = T+273.15 # Temperature in Kelvin
    P_inMPa = P/1e6  # Pressure in MPa

    H20 = IAPWS95(T=T_inK,P=P_inMPa) # Thermophysical (Stream 1)
    rho = H20.rho
    cp = H20.cp
    mu = H20.mu
    k = H20.k

    return rho, cp, mu, k    

#----------------------------------------------------------------------------------#
# Options/Inputs

# Solver Settings
nodes = 100  # Number of nodes for 1-D calculation

# Input Conditions
# Temperatures
T_1in = 30 # Inlet Temperature (Stream 1) - C
T_2in = 50 # Inlet Temperature (Stream 2) - C
T_3in = 80 # Inlet Temperature (Stream 3) - C
# Pressures
P_1in = 101325 # Inlet Pressure (Stream 1) - Pa
P_2in = 101325 # Inlet Pressure (Stream 2) - Pa
P_3in = 101325 # Inlet Pressure (Stream 3) - Pa
# Flow rates
V_flow_1 = 1 # Flow rate (Stream 1) - m3/sec
V_flow_2 = 5 # Flow rate (Stream 2) - m3/sec
V_flow_3 = 3 # Flow rate (Stream 3) - m3/sec

# Call Fluid properties
'''
Currently assuming liquids and not gases.
It currently calls the IAPWS95 database using the IAPWS python module. I created a
wrapper function to reduce overhead for the user. It should be portable enough for
future users if other fluid properties are needed.
'''
[rho_1,cp_1,mu_1,k_1] = thermophys_H20_iapws95(T_1in,P_1in) # Thermophysical (Stream 1)
[rho_2,cp_2,mu_2,k_2] = thermophys_H20_iapws95(T_2in,P_2in) # Thermophysical (Stream 2)
[rho_3,cp_3,mu_3,k_3] = thermophys_H20_iapws95(T_3in,P_3in) # Thermophysical (Stream 3)

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
Pr_3 = Prandtl(cp_3,mu_3,k_3) # Prandtl (Stream 3)

# Call Nusselt number calculation

# Call HTC Calculation

#----------------------------------------------------------------------------------#
# Boundary Value Problem Calculation

#----------------------------------------------------------------------------------#
# Output Formatted Text

print(Pr_1)
print(Pr_2)
print(Pr_3)

#----------------------------------------------------------------------------------#
# Output relevant plots