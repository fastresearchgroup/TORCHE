# -*- coding: utf-8 -*-
# Author: Jonah Haefner
#Contributing and previous work done by: Julien Clayton and Lane Carasik
import numpy
#####Units####
#a = dimensionless transverse pitch
#b = dimensionless longitudinal pitch
#v = the free stream fluid velocity in m/s 
#rho = density in g/cm^3
#u = dynamic viscosity in pa*s
#N = number of rows in the tube bundle 
#Dtube = The diamter of the tubes in the bundle in m
#Will work for values of b = 1.25, 1.5, or 2
#Reynolds number needs to be below 150000

# Output Pressure Drop in kPa


def integer(a,b,v,rho,u,N,Dtube, Re):
    rho = rho *1000;
    x = (a-1)/(b-1);
    v_mean = v*(a/(a-1));
    
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
    #There may be errors when b=2 and Re < 800 but didn't input for code
    #simplicity.
    
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

  #Power series for Euler number per row. From same website as above.
    Eu_p = (c_0[i]/Re**0)+(c_1[i]/Re**1)+(c_2[i]/Re**2)+(c_3[i]/Re**3)+(c_4[i]/Re**4);
    Eu = Eu_p*k_3*k_1;
                    
    #Using the relation Eu = dP/((1/2)*rho*v^2)
    dP = Eu*((rho*v_mean**2)/2); #Pressure drop per row 
    dP_total = dP*N/1000; #pressure drop across N rows    
    return dP_total
    
def array(a,b,v,rho,u,N,Dtube, Re_A):
    rho = rho * 1000;
    dPList = numpy.zeros(len(Re_A));
    for t in range (0,len(Re_A)):
        Re = Re_A[t]
        x = (a-1)/(b-1);
        v_mean = v*(a/(a-1));
        
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
        #There may be errors when b=2 and Re < 800 but didn't input for code
        #simplicity.
        
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
       
      #Power series for Euler number per row. From same website as above.
        Eu_p = (c_0[i]/Re**0)+(c_1[i]/Re**1)+(c_2[i]/Re**2)+(c_3[i]/Re**3)+(c_4[i]/Re**4);
        Eu = Eu_p*k_3*k_1;
                        
        #Using the relation Eu = dP/((1/2)*rho*v^2)
        dP = Eu*((rho*v_mean**2)/2); #Pressure drop per row 
        dP_total = dP*N/1000; #pressure drop across N rows
        dPList[t]= dP_total
    return dPList

def inline(a,b,v,rho,u,N,Dtube, Re):
    if type(Re) == list:
        Re = numpy.array(Re,dtype=float)
    if type(Re) == int:
        PD= testtype.integer(a,b,v,rho,u,N,Dtube, Re)
    elif type(Re) == numpy.ndarray:
        PD=testtype.array(a,b,v,rho,u,N,Dtube, Re)
    return PD 
        
