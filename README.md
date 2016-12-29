# PressureDrop
This is a collection of Python scripts for pressure drop calculations for crossflow tube bundles for both inline and staggered configurations

Author: Jonah Haefner
Contributing and previous work done by: Julien Clayton and Lane Carasik

These models are simple to use pressure drop calculators based off the work of Zukauskas et al. and Gnielinski and Gaddis.

###Zukauskas Model

Inputs:
- a = dimensionless transverse pitch
- b = dimensionless longitudinal pitch
- v = the free stream fluid velocity in m/s 
- rho = density in g/cm^3
- u = dynamic viscosity in pa*s
- N = number of rows in the tube bundle 
- Dtube = The diamter of the tubes in the bundle in m

Outputs:
- Pressure Drop in kPa

Limitations:
- Will work for values of b = 1.25, 1.5, or 2
- Reynolds number needs to be below 150000

This work is discussed in detail in the following work:

1. J. B. Haefner, L. B. Carasik, Y. A. Hassan, Development and Verification of MATLAB Functions for Pressure Drop Estimation in Crossflow Tube Bundles. In Transactions of the American Nuclear Society (Vol. 115(1), pp. 1535-1537). Las Vegas, Nevada

References:
1. A. ZUKAUSKAS, R. ULINSKAS, Heat Transfer in Banks of Tubes in Crossflow, Hemisphere Publishing, United States (1988).

###Gnielinski & Gaddis Model
Inputs
- Rho [kg/m^3]
- u = free stream velocity without tubes [m/s]
- a = transverse pitch ratio
- b = longitudinal pitch ratio
- - a and b [m]
- - Reynolds number calculated based on velcoity at narrowest point
- N = # of tubes

Calculated Variables
- D_lam = Drag coefficient due to laminar flow
- D_turb = Drag coefficient due to turbulent flow
- f_nt = Coefficient for influence of inlet and outlet pressure losses
- f_ti, f_ts = geometric arragement factor  

Notes:
- Not limited to specific pitch-diameter ratios
- Calculates the velcoity through the narrowest area
Validity:
* The available experimental measurements used to develop the above equations covered the following ranges: Reynolds number: 1 ≤ Re ≤ 3 /* 10E5 and number of tube rows: nR ≥  5.
* The tube bundles used in the experimental investigations had transverse and longitudinal pitch ratios as follows: 
 * In the range Re < 10E3: For in-line tube arrangement: a x b = 1:25 x 1:25; 1:5 x 1:5; 2:0 x 2:0. For staggered tube arrangement: a x b = 1:25 x 1:0825; 1:5 x 1:299; 1:768 x 0:884.
 * In the range Re ≥ 10E3: For in-line tube arrangement: 1:25 ≤ a ≤ 3:0; 1:2 ≤  b ≤  3:0. For staggered tube arrangement: 1:25 ≤  a ≤  3:0; 0:6 ≤  b ≤  3:0; c ≥ 1:25.
 Citation: VDI. VDI Heat Atlas. Berlin, Heidelberg: VDI-Buch-Springer, 2010. Chapter L1
Output:
- Pressure Drop in kPa
