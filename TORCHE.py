# -*- coding: utf-8 -*-
"""
Date of last edit: 2 APR 2017
Author: Jonah Haefner
Texas A&M University
Thermal-hydraulics Group

Title: Heat Transfer Models for Crossflow Heat Exchangers

Functions included:
    Zukauskas
    Gnielinski
    PressureDrop
"""

import TORCHE
import numpy
import sys

def dP_Zu(a,b,v,rho,u,N,Dtube, Re,geom,eta_wall=1, eta=1,alpha=0,beta=90):
''' #####Units####
# #a = dimensionless transverse pitch
# #b = dimensionless longitudinal pitch
# #v = the free stream fluid velocity in m/s 
# #rho = density in kg/m^3
# #u = dynamic viscosity in pa*s
# #N = number of rows in the tube bundle 
# #Dtube = The diameter of the tubes in the bundle in m
# #Will work for pitch diamters of 1.25, 1.5, or 2
# #Reynolds number needs to be below 150000

# Output Pressure Drop in kPa
# Rotated Crossflow (k4) and Inclined Crossflow(k5)
    #  angle of attack, α, may not be 0°, i.e., the bank may 
    # be rotated at some arbitrary angle to the flow 
    # The angle of incidence β is not 90°, i.e. the flow is 
    # not exactly perpendicular to the tube bank. 
    
    Return dP''''
    x = (a-1)/(b-1);
    v_mean = v*(a/(a-1));
    if geom == 'inline':
        # These are valid for Reynolds numbers between 2E3 and 2E6 %
        if Re > 2E3:
            c_0 = numpy.array([.267, .235, .247]);      #b = 1.25, 1.5, and 2
            c_1 = numpy.array([.249E4, .197E4, -.595]) #b = 1.25, 1.5, and 2
            c_2 = numpy.array([-.927E7, -.124E8, .15]); #b = 1.25, 1.5, and 2
            c_3 = numpy.array([.1E11, .312E11, -.137]); #b = 1.25, 1.5, and 2
            c_4 = numpy.array([0, -.274E14, .396]);     #b = 1.25, 1.5, and 2
        elif Re <= 800 and Re > 3:
            c_0 = numpy.array([.272, .263, .188]);      #b = 1.25, 1.5, and 2
            c_1 = numpy.array([.207E3, .867E2, 56.6]);  #b = 1.25, 1.5, and 2
            c_2 = numpy.array([.102E3, -.202, -646]);   #b = 1.25, 1.5, and 2
            c_3 =numpy.array([-.286E3, 0, 6010]);      #b = 1.25, 1.5, and 2
            c_4 = numpy.array([0, 0, -18300]);          #b = 1.25, 1.5, and 2
        elif Re <= 2E3 and Re > 800:
            c_0 =numpy.array([.272, .263, .247]);      #b = 1.25, 1.5, and 2
            c_1 = numpy.array([.207E3, .867E2, 56.6-(4.766E-2*Re)]); #b = 1.25, 1.5, and 2
            c_2 = numpy.array([.102E3, -.202, .15]); #b = 1.25, 1.5, and 2
            c_3 = numpy.array([-.286E3, 0, -.137]); #b = 1.25, 1.5, and 2
            c_4 = numpy.array([0, 0, .396]);     #b = 1.25, 1.5, and 2
        #There may be errors when b=2 and Re < 800 
        
        ##########Correction Factors #################
        
        # Non-rectangular bundle # (valid up to Reynolds numbers of 150000)
        
        k_1 = 1;
        
        if x != 1.00:
            
            if (Re >= 1000 and Re < 10000):
            
                k_1 = .9849*x**(-.8129);
            
            elif (Re >= 10000 and Re < 70000):  
            
                k_1 = .9802*x**(-.7492);
            
            elif (Re >= 70000 and Re < 150000):
                            
                k_1 = .988*x**(-.6388);
                        
         # Entry Losses #
                        
        # Entry loss coefficients (Re > 1E4 but < 1E6) 
        el_1 = numpy.array([1.9, 1.1, 1, 1, 1, 1, 1]);
        
        #Entry loss coefficients (Re > 1E6)
        el_2 = numpy.array([2.7, 1.8, 1.5, 1.4, 1.3, 1.2, 1.2]);
         
        
        if (b == 1.25):
            
            i = 0;
            
        elif (b == 1.5):
                
            i = 1;
                
        elif (b == 2):
                    
            i = 2;
                    
        k_3 = 1;
        
        if (N < 7 and N > 0 and Re < 1E6 and Re > 1E2):
               
            k_3 = el_1[N-1];
            
        elif (N < 9 and N > 0 and Re > 1E6):
            
            k_3 = el_2[N-1];
               
        elif (Re <= 1E2):
            print('Reynolds number needs to be greater than 100')
    
    if geom == 'staggered':
            #Mean Velocity Calculations. Found on same website as above
            if a <= (2*b**2 -.5):
                v_mean = v*(a/(a-1));
            if a > (2*b**2 -.5):
                v_mean = v*(a/(numpy.sqrt(4*b^2+a^2)-2))                     
            if Re < 1E3:
                c_0 = [.795,.683,.343];      #b = 1.25, 1.5, and 2
                c_1 = [.247E3,.111E3,.303E3]; #b = 1.25, 1.5, and 2
                c_2 = [.335E3,-.973E2,-.717E5]; #b = 1.25, 1.5, and 2
                c_3 = [-.155E4,-.426E3,.88E7]; #b = 1.25, 1.5, and 2
                c_4 = [.241E4,.574E3,-.38E9];     #b = 1.25, 1.5, and 2
            elif Re >= 1E3 and Re < 1E4:
                c_0 = [.245,.203,.343];      #b = 1.25, 1.5, and 2
                c_1 = [.339E4, .248E4, .303E3 ]; #b = 1.25, 1.5, and 2
                c_2 = [-.984E7, -.758E7, -.717E5]; #b = 1.25, 1.5, and 2
                c_3 = [.132E11, .104E11, .88E7]; #b = 1.25, 1.5, and 2
                c_4 = [-.599E13, -.482E13, -.38E9];     #b = 1.25, 1.5, and 2
            elif Re >= 1E4:
                c_0 = [.245, .203, .162  ];   #b = 1.25, 1.5, and 2
                c_1 = [.339E4, .248E4, .181E4 ]; #b = 1.25, 1.5, and 2
                c_2 = [-.984E7, -.758E7, .792E8]; #b = 1.25, 1.5, and 2
                c_3 = [.132E11, .104E11, -.165E13]; #b = 1.25, 1.5, and 2
                c_4 = [-.599E13, -.482E13, .872E16];  #b = 1.25, 1.5, and 2
            
            if (b == 1.25):                
                i = 0;
                
            elif (b == 1.5):
                i = 1;
                    
            elif (b == 2):
                i = 2;
            
            #Entry loss coefficients (Re > 1E2 but < 1E4)
            el_1 = [1.4, 1.3, 1.2, 1.1, 1, 1, 1];
            #Entry loss coefficients (Re > 1E4 but < 1E6) 
            el_2 = [1.1, 1.05, 1, 1, 1, 1, 1];
            #Entry loss coefficients (Re > 1E6)
            el_3 = [.25, .45, .6, .65, .7, .75, .8];

            k_3 = 1; #Setting for tubes > 7
            if (N < 7 and N > 0 and Re < 1E4 and Re > 1E2):
                k_3 = el_1(N);
            elif (N < 7 and N > 0 and Re < 1E6 and Re >= 1E4):
                k_3 = el_2(N);
            elif (N < 7 and N > 0 and Re > 1E6):
                k_3 = el_3(N);
   #Has not been added in yet: k_1 for unequal pitch-diameter ratios
   # Can be found in graphical form in Zukauskas' High Performance Single Phase Heat Exchangers
   # Electronically found in cubic splines based of Re at: 
   # http://trace.tennessee.edu/cgi/viewcontent.cgi?article=1949&context=utk_gradthes
            k_1 = 1;
    '''Influence of temperature on fluid properties
    From: Heat Exchanger Design Handbook 1983'''
    p=1;
    if eta_wall>eta:
        p=0.776*numpy.exp(-0.545*Re**0.256)
    elif eta_wall < eta and Re < 10E3:
        p=0.968*numpy.exp(-1.076*Re**0.196)
    k_2=(eta_wall/eta)**p
    
    '''
    Deviation form normal incidence
    '''
    k_4 = numpy.cos(numpy.deg2rad(alpha))
    k_5 = 1
    if beta != 90:
        k_5 = -1E-06*beta**3 + 1E-04*beta**2 + 0.0132*beta - 0.033
    #R² = 0.9994. Note this is a curve fitted to experimental data and not exact

  #Power series for Euler number per row. From same website as above.
    Eu_p = (c_0[i])+(c_1[i]/Re**1)+(c_2[i]/Re**2)+(c_3[i]/Re**3)+(c_4[i]/Re**4);
    
    Eu = Eu_p*k_1*k_2*k_3*k_4*k_5;    
             
    #Using the relation Eu = dP/((1/2)*rho*v^2)
    dP = Eu*((rho*v_mean**2)/2); #Pressure drop per row 
    dP_total = dP*N/1000; #pressure drop across N rows in kPa   
    return dP_total


def inline(a,b,v,rho,u,N,Dtube, Re):
    geom = 'inline'
    PD= TORCHE.dP_Zu(a,b,v,rho,u,N,Dtube, Re,geom)
    return PD
    
def staggered(a,b,v,rho,u,N,Dtube, Re):
    geom = 'staggered'
    PD= TORCHE.dP_Zu(a,b,v,rho,u,N,Dtube, Re,geom)
    return PD 
    
    
    
def HT_Zu(a,b,u,rho,N,Pr,PrS,d,v,geom,Re):
    '''
        Zukauskas Model
        Cite: Zhukauskas, A., R. Ulinskas. Heat Transfer in Tube Banks in Crossflow.
        Hemisphere Publishing Corporation. New York, NY. 1988.
        
        Inputs:
            # Rho = density [kg/m^3]
            # u = free stream velocity without tubes [m/s]
            # a = transverse pitch ratio
            # b = longitudinal pitch ratio
            # a and b [m]
            # Reynolds number calculated based on velcoity at narrowest point
            # N = # of tubes
            # Pr = Prandtl Number
            # d = Tube Diameter
            # PrS = Wall Prandtl Number
            # v = kinematic viscosity (m^2/s)
            
            NOTE: Currently only applicable for staggered geometry
            Validity:10 < Re < 2e6
                     0.7 < Pr < 500
        Return Nusselt Number
        '''
    c = numpy.sqrt(b**2+(a/2)**2);
    if c > (a+d)/2:
        Vmax= a/(a-d)*u;
    else:
        Vmax= a/(2(b-d))*u;
#    Re= d*Vmax/v
    

    if geom == 'inline':  
        n= 0.36;
        ###Correction for less than 20 tube bundles
        if N <=14:
            c2v = [.7,.8,.85,.9,.92,.94,.95,.96,.97,.97,.98,.98,.99,.99] ;
            c2 = c2v[N-1];
    ### Coefficients based on range of Reynolds Number
        if Re > 1.6 and Re <= 100:
            c1 = 0.90;
            m =  0.40;
        elif Re > 100 and Re <= 1e3:
            c1=.52;
            m= .5;
        elif Re > 1e3 and Re <= 2e5:
            c1 = .27;
            m =  .63;
        elif Re > 2e5:
            c1 = .033;
            m = .83;
        else:
            print( 'Re out of Range')
    
    if geom == 'staggered':
            ###Correction for less than 20 tube bundles
        if Re <= 1000 and N < 20:
            c2= 1-numpy.exp(-numpy.sqrt(3*N**(1/numpy.sqrt(2))));    
        if Re > 1000 and N < 20:
            c2 = 1-numpy.exp(-N**(1/numpy.sqrt(3)));
         ### Coefficients based on range of Reynolds Number
        if Re <= 40 and Re > 1.6:
            c1 = 1.04;
            m= .4;
            n=.36;
        elif Re > 40 and Re <= 1000:
            c1 = 0.52;
            m = 0.5;
            n = 0.36;
        elif Re > 1000 and Re <= 2e5:
            if a/b < 2:
                c1 = 0.35*(a/b)**0.2;
            else:
                c1 = 0.40;
            m = 0.6;
            n=0.36;
        elif Re > 2e5 and Re <= 2e6:
            c1 = 0.022;
            m = 0.84;
            n = 0.36;
        else:
            print( 'Re out of Range')
    
    Nu = c1*c2*(Re**m)*(Pr**n)*(Pr/PrS)**.25;
    return Nu
    
def dP_GG(a,b,u,rho,Re,N,geom,Return=""):
	'''
	Description:
		Calculating the velcoity through the narrowest area
	Inputs:
		rho = density of fluid [kg/m^3]
        u = free stream velocity [m/s]
        a = transverse pitch ratio
        b = longitudinal pitch ratio
        Re = Reynolds number calculated based on velcoity at narrowest point
        N = Number of tubes
        D_lam = Drag coefficient due to laminar flow
        D_turb = Drag coefficient due to turbulent flow
        f_nt = Coefficient for influence of inlet and outlet pressure losses
        f_ti, f_ts = geometric arragement factor  
	Outputs:
		dP = Pressure drop across tube bundle [Pa]
	Warnings:
		Ranges of validity have yet to be included in this script. Information soon to come.
        '''
   if geom == 'inline':
        u0 = u*(a/(a-1))
        D_lam= 280*numpy.pi*((b**(.5)-0.6)**2+0.75)/(a**(1.6)*(4*a*b-numpy.pi)*Re);
        if (N > 5 and N <= 10):
            f_nt= (1/a**2)*(1/N-.1);
        elif N > 10:
            f_nt = 0
        elif N <=5:
            sys.exit('Model only valid for greater than 5 tubes')
        f_ti = (0.22+1.2*(1-(0.94/b))**(0.6)/(a-0.85)**1.3)*10**(0.47*(b/a-1.5))+0.03*(a-1)*(b-1);
        D_turb = f_ti/(Re**(.1*b/a));
        D_tot = D_lam+(D_turb+f_nt)*(1-numpy.exp(-(Re+1000)/2000));
   if geom == "staggered":
        N = N-1;
        if  b < .5*(2*a+1)**(.5):
            c=((a/2)**2+b**2)**(1/2);
            u0 = u*a/(2*(c-1));
            D_lam= 280*numpy.pi*((b**.5-0.6)**2+0.75)/((c**1.6)*(4*a*b-numpy.pi)*Re);
            if (N > 5 and N <= 10):
                f_nt = (2*(c-1)/a*(a-1)**2)*(1/N-.1);
            elif N > 10:
                f_nt = 0;
        else:
            u0 = u*(a/(a-1))
            D_lam= 280*numpy.pi*((b**(.5)-0.6)**2+0.75)/(a**(1.6)*(4*a*b-numpy.pi)*Re);
            if (N > 5 and N <= 10):
                f_nt= (1/a**2)*(1/N-.1);
            elif N > 10:
                f_nt = 0
        f_ts = 2.5+1.2/(a-0.85)**1.08+0.4*(b/a-1)**3-0.01*(a/b-1)**3;
        D_turb = f_ts/Re**0.25;
        D_tot = D_lam+(D_turb+f_nt)*(1-numpy.exp(-(Re+200)/1000));
   if Return == "D_tot":
        return D_tot    
   PressDrop = .5*D_tot*(N*rho*u0**2); # in kPa
   return PressDrop
def HT_GG(a,b,d,u,rho,Re,N,geom,Pr):
    '''
    Gnielinski and Martin Model
    Leveque Analogy
    Cite: Martin, H., 2002, “The Generalized Lévêque Equation and its practical use for the prediction of heat and mass transfer rates from pressure drop,”
    Chem. Eng. Sci., vol. 57, pp. 3217-3223.
        # Rho [kg/m^3]
        # u = free stream velocity without tubes [m/s]
        # a = transverse pitch ratio
        # b = longitudinal pitch ratio
        # a and b [m]
        # Reynolds number calculated based on velcoity at narrowest point
        # N = # of tubes
        #geom = tube geometry (inline or staggered)
        #d= tube diameter
        # dh = hydraulic diameter
        # Pr = Prandtl Number
        
    Return Nusselt Number
    '''
    xi = dP_GG(a,b,u,rho,Re,N,geom,"D_tot") # Pressure Drop Coefficient from Gaddis-Gnielinski Pressure Drop Model
    if Re > 2.5e5: #Correction made my Holger Martin for use in Heat Transfer calculation 
        xi = xi*(1+(Re-2.5e5)/3.25e5); # Does the xi need to be in kPa?
    xi_f = .5*xi; #Original total drag coefficient used to calculate drag coefficient due to friction ~.5
    if b > 1:
       dh=((4*a/numpy.pi)-1)*d;
    else:
       dh=((4*a*b/numpy.pi)-1)*d;
    if geom == 'staggered':
        c=((a/2)**2+b**2)**.5;
        L=c*d;
    if geom == 'inline':
        L=b*d;
    NuPr=0.404*(xi_f*Re**2*dh/L)**(1/3); #Nu/Pr**(1/3)
    Nu = NuPr*(Pr)**(1/3);
    return Nu
    
