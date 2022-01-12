'''
Author: Lane Carasik
Purpose: Three-stream heat exchanger temperature distribution (1-D)
'''

#----------------------------------------------------------------------------------#
# import modules

import numpy as np
from scipy.integrate import solve_bvp
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

def Nu_DB(Re,Pr,horc):
    # Nusselt number correlation (Dittus Boelter)
    # for smooth pipe flow
    # cooling so n = 0.3 or n = 0.4 for heating
    if horc == 'h':
    	n = 0.4
    elif horc == 'c':
    	n = 0.3

    Nu = 0.023*(Re**0.8)*(Pr**n)

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
    cp = H20.cp*1000
    mu = H20.mu
    k = H20.k

    return rho, cp, mu, k    

#----------------------------------------------------------------------------------#
# Options/Inputs
'''
This assumes a three stream concentenric heat exchanger 
Stream 1 - Inner most tube (D_h = D_1in)
Stream 2 - Intermediate Annulus (D_h = D_2in - D_1out)
Stream 3 - Outer Annulus (D_h = D_3in - D_2out)
'''

# Solver Settings
nodes = 100  # Number of nodes for 1-D calculation

# Input Conditions
# Geometry of Heat Exchanger
D_1in = 10/100 # m - Inner Diameter (Tube 1)
D_1out = 11/100 # m - Outer Diamter (Tube 1)
D_2in = 13/100 # m - Inner Diamter (Tube 2)
D_2out = 14/100 # m - Outer Diamter (Tube 2)
D_3in = 16/100 # m - Inner Diameter (Tube 3)
D_3out = 19/100 # m - Outer Diameter (Tube 3)
L_HX = 200/100 # m - Length (HX)

i_2 = -1 # Add stream definition
i_3 = -1

# Temperatures
T_1in = 30 # Inlet Temperature (Stream 1) - C
T_2in = 50 # Inlet Temperature (Stream 2) - C
T_3in = 80 # Inlet Temperature (Stream 3) - C
# Pressures
P_1in = 101325 # Inlet Pressure (Stream 1) - Pa
P_2in = 101325 # Inlet Pressure (Stream 2) - Pa
P_3in = 101325 # Inlet Pressure (Stream 3) - Pa
# Flow rates
V_flow_1 = 0.00454609 # Flow rate (Stream 1) - m3/sec
V_flow_2 = 0.0227304  # Flow rate (Stream 2) - m3/sec
V_flow_3 = 0.0136383  # Flow rate (Stream 3) - m3/sec

# Call Fluid properties
'''
Currently assuming liquids and not gases.
It currently calls the IAPWS95 database using the IAPWS python module. 
I created a wrapper function to reduce overhead for the user. 
It should be portable enough for future users if other fluid properties 
are needed.
'''
[rho_1,cp_1,mu_1,k_1] = thermophys_H20_iapws95(T_1in,P_1in) # Thermophysical (Stream 1)
[rho_2,cp_2,mu_2,k_2] = thermophys_H20_iapws95(T_2in,P_2in) # Thermophysical (Stream 2)
[rho_3,cp_3,mu_3,k_3] = thermophys_H20_iapws95(T_3in,P_3in) # Thermophysical (Stream 3)

# Call heat exchanger tube wall properties
'''
This just assumes Stainless Steel 316L annealed sheets
http://www.matweb.com/search/datasheet_print.aspx?matguid=1336be6d0c594b55afb5ca8bf1f3e042.
Future work would be to make this more general.
'''
kw_12 = 15.9 # Thermal conductivity of wall (Stream 1-2) - W/m-K
kw_23 = 15.9 # Thermal conductivity of wall (Stream 1-2) - W/m-K

#----------------------------------------------------------------------------------#
# Basic Global Calulations

# Call area calculation
A_1 = np.pi*(D_1in/2)**2 					   # Flow area of stream 1 - m2
A_2 = np.pi*(D_2in/2)**2 - np.pi*(D_1out/2)**2 # Flow area of stream 2 - m2
A_3 = np.pi*(D_3in/2)**2 - np.pi*(D_2out/2)**2 # Flow area of stream 3 - m2

# Call Hydraulic diameter calculation
Dh_1 = D_1in 			# Hydraulic Diameter (Stream 1) - m
Dh_2 = D_2in - D_1out	# Hydraulic Diameter (Stream 2) - m
Dh_3 = D_3in - D_2out	# Hydraulic Diameter (Stream 3) - m

# Call velocity calculation (velocity = flow rate/flow area)
vel_1 = V_flow_1/A_1 # Velocity (Stream 1) - m/sec
vel_2 = V_flow_2/A_2 # Velocity (Stream 2) - m/sec
vel_3 = V_flow_3/A_3 # Velocity (Stream 3) - m/sec

# Call mass flow rate calculation (mdot = rho*area*vel = rho*vol_flow)
mdot_1 = rho_1*V_flow_1 # Mass Flow Rate of Stream 1 - kg/sec
mdot_2 = rho_2*V_flow_2 # Mass Flow Rate of Stream 2 - kg/sec
mdot_3 = rho_3*V_flow_3 # Mass Flow Rate of Stream 3 - kg/sec

# Call heat capacity rates (C = mdot*cp)
C_1 = mdot_1*cp_1 # Heat Capacity of Stream 1
C_2 = mdot_2*cp_2 # Heat Capacity of Stream 2
C_3 = mdot_3*cp_3 # Heat Capacity of Stream 3

# Call Reynolds number calculation
Re_1 = Reynolds(rho_1,vel_1,Dh_1,mu_1) # Re Number (Stream 1)
Re_2 = Reynolds(rho_2,vel_2,Dh_2,mu_2) # Re Number (Stream 2)
Re_3 = Reynolds(rho_3,vel_3,Dh_3,mu_3) # Re Number (Stream 3)

# Call Prandtl number calcation
Pr_1 = Prandtl(cp_1,mu_1,k_1) # Prandtl (Stream 1)
Pr_2 = Prandtl(cp_2,mu_2,k_2) # Prandtl (Stream 2)
Pr_3 = Prandtl(cp_3,mu_3,k_3) # Prandtl (Stream 3)

# Call Nusselt number calculation
Nu_1 = Nu_DB(Re_1,Pr_1,horc='h') # Nu Number (Stream 1)
Nu_2 = Nu_DB(Re_2,Pr_2,horc='h') # Nu Number (Stream 2)
Nu_3 = Nu_DB(Re_3,Pr_3,horc='c') # Nu Number (Stream 3)

# Call Heat Transfer Coefficient Calculation
HTC_1 = Nu_1*k_1/Dh_1	# HTC (Stream 1) - W/m2-K
HTC_2 = Nu_2*k_2/Dh_2	# HTC (Stream 2) - W/m2-K
HTC_3 = Nu_3*k_3/Dh_3	# HTC (Stream 3) - W/m2-K

# Call Overall Heat Transfer Coefficient Calculation
# Stream 1 to 2
UA_12 = 1/(1/(HTC_2*A_2)+1/(2*np.pi*kw_12*L_HX)*np.log(D_1out/D_1in)+1/(HTC_1*A_1))
# Stream 2 to 3
UA_32 = 1/(1/(HTC_3*A_3)+1/(2*np.pi*kw_23*L_HX)*np.log(D_2out/D_2in)+1/(HTC_2*A_2))

# Call NTU calculation
NTU_1 = UA_12/C_1

# Call Heat Capacity Stream Ratio Calculation
Cs_12 = C_1/C_2 # Heat capacity ratio (1->2)
Cs_32 = C_3/C_2 # Heat capacity ratio (3->2)

# Call Conductance Ratio calculation
R_star = UA_32/UA_12

# Call inlet temperature ratio (Dimensionless)
theta_3in = (T_3in-T_1in)/(T_2in-T_1in) 

#----------------------------------------------------------------------------------#
# Boundary Value Problem Calculation

# See notes from Selic three stream heat exchanger design

# Call Function to Solve
def fun(x,y):
	return np.vstack((NTU_1*(y[1]-y[0]),
		i_2*NTU_1*Cs_12*(y[0]-y[1])+i_2*NTU_1*R_star*Cs_12*(y[2]-y[1]),
		i_3*(Cs_12/Cs_32)*R_star*NTU_1*(y[1]-y[2])
	))

#Call Boundary Conditions to use
def bc(ya,yb):
	if i_2 == 1:
		res_2 = ya[1]-theta_3in # theta_3in might be the wrong one here.
	elif i_2 == -1:
		res_2 = yb[1]-theta_3in
	if i_3 == 1:
		res_3 = ya[2]
	elif i_3 == -1:
		res_3 = yb[2]
	return np.array([ya[0]-1,res_2,res_3])

# Solve and make sure it can handle different directions of the 
# parallel/counter flow heat exchanger
x = np.linspace(0,1,nodes)
y = np.zeros((3,x.size))
 
#Call bvp solver
sol = solve_bvp(fun,bc,x,y,max_nodes = 1e6,verbose=2)

Theta_1 = sol.sol(x)[0] # Stream 1
Theta_2 = sol.sol(x)[1] # Stream 2
Theta_3 = sol.sol(x)[2] # Stream 3

print(Theta_1)
print(Theta_2)
print(Theta_3)

#----------------------------------------------------------------------------------#
# Output Formatted Text

print(Pr_1)
print(Pr_2)
print(Pr_3)
print(Re_1)
print(Re_2)
print(Re_3)
print(Nu_1)
print(Nu_2)
print(Nu_3)
print(HTC_1)
print(HTC_2)
print(HTC_3)

#----------------------------------------------------------------------------------#
# Output relevant plots

h = 10
w = 8

lw = 2.5
fs = 14

k = 1
plt.figure(k, figsize=(h,w))
plt.plot(x,Theta_1,'r--',linewidth = lw,label=r'$\theta_{1}$')
plt.plot(x,Theta_2,'k--',linewidth = lw,label=r'$\theta_{2}$')
plt.plot(x,Theta_3,'b--',linewidth = lw,label=r'$\theta_{3}$')
#plt.axhline(y=Tboil, xmin = 0, xmax = 1, color = 'r',linewidth = lw, label='Coolant Boiling Temp')
#plt.axhline(y=TmeltCoolant, xmin = 0, xmax = 1, color = 'b',linewidth = lw, label='Coolant Melting Temp')
plt.legend(loc='center left')
plt.xlabel('Length - m',fontsize = fs)
plt.ylabel('Theta - ',fontsize = fs)
plt.grid()
# fig = TCinput.Reactor_Title + '_Axial_Temperatures'
# if print_logic == 0:
#     plt.savefig(fig + '.png', dpi = 300, format = "png",bbox_inches="tight")
#     plt.savefig(fig + '.eps', dpi = 300, format = "eps",bbox_inches="tight")
#     plt.savefig(fig + '.svg', dpi = 300, format = "svg",bbox_inches="tight")
k = k + 1

plt.show()