# -*- coding: utf-8 -*-
# Author: Jonah Haefner
#Contributing and previous work done by: Julien Clayton and Lane Carasik
import numpy
#==============================================================================
# #####Units####
# #a = dimensionless transverse pitch
# #b = dimensionless longitudinal pitch
# #v = the free stream fluid velocity in m/s 
# #rho = density in kg/m^3
# #u = dynamic viscosity in pa*s
# #N = number of rows in the tube bundle 
# #Dtube = The diamter of the tubes in the bundle in m
# #Will work for pitch diamters of 1.25, 1.5, or 2
# #Reynolds number needs to be below 150000
# Output Pressure Drop in kPa
#==============================================================================
    
import PressureDrop



def Zukauskas(a,b,v,rho,u,N,Dtube, Re,geom):
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
                i = 1;
                
            elif (b == 1.5):
                i = 2;
                    
            elif (b == 2):
                i = 3;
            
            #Entry loss coefficients (Re < 1E2 but < 1E4)
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

            k_1 = 1
  #Power series for Euler number per row. From same website as above.
    Eu_p = (c_0[i])+(c_1[i]/Re**1)+(c_2[i]/Re**2)+(c_3[i]/Re**3)+(c_4[i]/Re**4);
    
    Eu = Eu_p*k_3*k_1;                
    #Using the relation Eu = dP/((1/2)*rho*v^2)
    dP = Eu*((rho*v_mean**2)/2); #Pressure drop per row 
    dP_total = dP*N/1000; #pressure drop across N rows in kPa   
    return dP_total


def inline(a,b,v,rho,u,N,Dtube, Re):
    geom = 'inline'
    PD= PressureDrop.Zukauskas(a,b,v,rho,u,N,Dtube, Re,geom)
    return PD
    
def staggered(a,b,v,rho,u,N,Dtube, Re):
    geom = 'staggered'
    PD= PressureDrop.Zukauskas(a,b,v,rho,u,N,Dtube, Re,geom)
    return PD 
        

    


