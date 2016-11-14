# PressureDrop
Python script for pressure drop calculations for crossflow tube bundles for both inline and staggered configurations

Author: Jonah Haefner
Contributing and previous work done by: Julien Clayton and Lane Carasik

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
