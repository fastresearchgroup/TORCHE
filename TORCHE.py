"""
Date of last edit: 2 APR 2017
Author(s): Jonah Haefner* and Lane Carasik^
*Texas A&M University
*Thermal-hydraulics Research Group

^Kairos Power LLC.
^Modeling and Simulation Group

TOolbox for Reactor Cross-Flow Heat Exchangers: 
Python Scripts for calculation of Pressure drop and Heat Transfer in crossflow tube bundles based on models found across the literature.

Functionality:
	Pressure Drop Models
		- Zhukauskas (Inline and Staggered)
		- Gaddis-Gnielinski (Inline and Staggered)
	Heat Transfer Models
		- Zhukauskas (Staggered)
		- Gaddis-Gnielinski (Inline and Staggered)
"""

import numpy as np
import sys

def dP_Zu(rho,a,b,geom,N,u,Re,eta_wall=1, eta=1,alpha=0,beta=90):
	'''
	Description:
		Calculating the presssure drop using the Zhukauskas correlation
	Inputs:
		rho = density of the fluid [kg/m^3]
		a = transverse pitch to diameter ratio
		b = longitudinal pitch to diameter ratio
		geom = tube geometry (inline or staggered) [string]
		N = Number of tubes
		u = free stream velocity [m/s]
        Re = Reynolds number calculated based on velcoity at narrowest point and tube outer diameter
	Output:
		dP_total = Pressure Drop across tubes [Pa]
	Warnings:
	
	Validity:	
		Reynolds number:
			10 ≤ Re ≤ 10e6
		Prandtl number:
			0.7 < Pr < 500
		Number of rows of tubes:
			Nr >= 1
		For in-line tube arrangement: 	a x b = 1.25 x 1.25;   1.5 x 1.5;   2.0 x 2.0. 
		For staggered tube arrangement: a x b = 1.25 x ----;   1.5 x ----;  2.0 x ----.

	Citation: Zhukauskas, A., R. Ulinskas. Heat Transfer in Tube Banks in Crossflow.
        Hemisphere Publishing Corporation. New York, NY. 1988.
	'''
	
	if Re < 10:
		print('The provided Reynolds number is out of the lower bounds (10) of validity at',Re)
	elif Re > 10E6:
		print('The provided Reynolds number is out of the higher bounds (10e6) of validity at',Re)
	
	Pr = 10 # Hardcoded for now.
	if Pr < 0.7:
		print('The provided Prandtl number is out of the lower bounds (0.7) of validity at')
	elif Pr > 500:
		print('The provided Prandtl number is out of the higher bounds (500) of validity at')
	
	if a != b:
		print('Has not been added in yet: k_1 for unequal pitch-diameter ratios Can be found in graphical form in Zukauskas High Performance Single Phase Heat Exchangers Electronically found in cubic splines based of Re at: http://trace.tennessee.edu/cgi/viewcontent.cgi?article=1949&context=utk_gradthes')
	
	x = (a-1)/(b-1)
	v_max = u*(a/(a-1))
	if geom in ['inline','INLINE','Inline','square','SQUARE','Square']:
        # These are valid for Reynolds numbers between 2E3 and 2E6
		if Re > 2E3:
			c_0 = np.array([.267, .235, .247])      #b = 1.25, 1.5, and 2
			c_1 = np.array([.249E4, .197E4, -.595]) #b = 1.25, 1.5, and 2
			c_2 = np.array([-.927E7, -.124E8, .15]) #b = 1.25, 1.5, and 2
			c_3 = np.array([.1E11, .312E11, -.137]) #b = 1.25, 1.5, and 2
			c_4 = np.array([0, -.274E14, .396])     #b = 1.25, 1.5, and 2
		elif Re <= 800 and Re > 3:
			c_0 = np.array([.272, .263, .188])      #b = 1.25, 1.5, and 2
			c_1 = np.array([.207E3, .867E2, 56.6])  #b = 1.25, 1.5, and 2
			c_2 = np.array([.102E3, -.202, -646])   #b = 1.25, 1.5, and 2
			c_3 =np.array([-.286E3, 0, 6010])		#b = 1.25, 1.5, and 2
			c_4 = np.array([0, 0, -18300])          #b = 1.25, 1.5, and 2
			if b == 2:
				print('There may be errors when b=2 and Re < 800')
		elif Re <= 2E3 and Re > 800:
			c_0 = np.array([.272, .263, .247])      #b = 1.25, 1.5, and 2
			c_1 = np.array([.207E3, .867E2, 56.6-(4.766E-2*Re)]) #b = 1.25, 1.5, and 2
			c_2 = np.array([.102E3, -.202, .15]) #b = 1.25, 1.5, and 2
			c_3 = np.array([-.286E3, 0, -.137]) #b = 1.25, 1.5, and 2
			c_4 = np.array([0, 0, .396])     #b = 1.25, 1.5, and 2
        
        ##########Correction Factors #################
        
        # Non-rectangular bundle # (valid up to Reynolds numbers of 150000)
        
		k_1 = 1;
		if x != 1.00:            
			if (Re >= 1000 and Re < 10000):
				k_1 = .9849*x**(-.8129)
			elif (Re >= 10000 and Re < 70000):  
				k_1 = .9802*x**(-.7492)
			elif (Re >= 70000 and Re < 150000):                
				k_1 = .988*x**(-.6388)
        # Entry Losses #                   
        # Entry loss coefficients (Re > 1E4 but < 1E6) 
		el_1 = np.array([1.9, 1.1, 1, 1, 1, 1, 1])
        #Entry loss coefficients (Re > 1E6)
		el_2 = np.array([2.7, 1.8, 1.5, 1.4, 1.3, 1.2, 1.2])
		if (b == 1.25):
			i = 0
		elif (b == 1.5):
			i = 1
		elif (b == 2):  
			i = 2
		k_3 = 1
        
		if (N < 7 and N > 0 and Re < 1E6 and Re > 1E2):
			k_3 = el_1[N-1]
		elif (N < 9 and N > 0 and Re > 1E6):
			k_3 = el_2[N-1]
		elif (Re <= 1E2):
			print('Reynolds number needs to be greater than 100')
			
	if geom in ['staggered','STAGGERED','Staggered','triangular','TRIANGULAR','Triangular']:
		if a <= (2*b**2 -.5):
			v_max = u*(a/(a-1))
		if a > (2*b**2 -.5):
			v_max = u*(a/(np.sqrt(4*b^2+a^2)-2))                     
		if Re < 1E3:
			c_0 = [.795,.683,.343]				#b = 1.25, 1.5, and 2
			c_1 = [.247E3,.111E3,.303E3]		#b = 1.25, 1.5, and 2
			c_2 = [.335E3,-.973E2,-.717E5]		#b = 1.25, 1.5, and 2
			c_3 = [-.155E4,-.426E3,.88E7]		#b = 1.25, 1.5, and 2
			c_4 = [.241E4,.574E3,-.38E9]		#b = 1.25, 1.5, and 2
		elif Re >= 1E3 and Re < 1E4:
			c_0 = [.245,.203,.343]			    #b = 1.25, 1.5, and 2
			c_1 = [.339E4, .248E4, .303E3]		#b = 1.25, 1.5, and 2
			c_2 = [-.984E7, -.758E7, -.717E5]	#b = 1.25, 1.5, and 2
			c_3 = [.132E11, .104E11, .88E7] 	#b = 1.25, 1.5, and 2
			c_4 = [-.599E13, -.482E13, -.38E9]	#b = 1.25, 1.5, and 2
		elif Re >= 1E4:
			c_0 = [.245, .203, .162  ]			#b = 1.25, 1.5, and 2
			c_1 = [.339E4, .248E4, .181E4 ]		#b = 1.25, 1.5, and 2
			c_2 = [-.984E7, -.758E7, .792E8]	#b = 1.25, 1.5, and 2
			c_3 = [.132E11, .104E11, -.165E13]	#b = 1.25, 1.5, and 2
			c_4 = [-.599E13, -.482E13, .872E16]	#b = 1.25, 1.5, and 2
		
		if (b == 1.25):                
			i = 0
		elif (b == 1.5):
			i = 1
		elif (b == 2):
			i = 2
	   
		#Entry loss coefficients (Re > 1E2 but < 1E4)
		el_1 = [1.4, 1.3, 1.2, 1.1, 1, 1, 1]
		#Entry loss coefficients (Re > 1E4 but < 1E6) 
		el_2 = [1.1, 1.05, 1, 1, 1, 1, 1]
		#Entry loss coefficients (Re > 1E6)
		el_3 = [.25, .45, .6, .65, .7, .75, .8]

		k_3 = 1; 								# Setting for tubes > 7
		if (N < 7 and N > 0 and Re < 1E4 and Re > 1E2):
			k_3 = el_1(N)
		elif (N < 7 and N > 0 and Re < 1E6 and Re >= 1E4):
			k_3 = el_2(N)
		elif (N < 7 and N > 0 and Re > 1E6):
			k_3 = el_3(N)

		k_1 = 1;
	'''
	Influence of temperature on fluid properties
	From: Heat Exchanger Design Handbook 1983
	'''
	p = 1
	if eta_wall>eta:
		p=0.776*np.exp(-0.545*Re**0.256)
	elif eta_wall < eta and Re < 10E3:
		p=0.968*np.exp(-1.076*Re**0.196)
	k_2=(eta_wall/eta)**p
    
	'''
	Deviation form normal incidence
	This needs to be better explained. It almost suggests that alpha CAN not be 0 and that beta CAN not be 90. 
	Rotated Crossflow (k4) and Inclined Crossflow(k5) angle of attack, α, may not be 0°, i.e., the bank maybe rotated at some arbitrary angle to the flow.
	The angle of incidence β is not 90°, i.e. the flow is not exactly perpendicular to the tube bank. 
	'''
	
	k_4 = np.cos(np.deg2rad(alpha))
	k_5 = 1
	if beta != 90:
		k_5 = -1E-06*beta**3 + 1E-04*beta**2 + 0.0132*beta - 0.033
    #R² = 0.9994. Note this is a curve fitted to experimental data and not exact

	#Power series for Euler number per row. From same website as above.
	Eu_p = (c_0[i])+(c_1[i]/Re**1)+(c_2[i]/Re**2)+(c_3[i]/Re**3)+(c_4[i]/Re**4)
	Eu = Eu_p*k_1*k_2*k_3*k_4*k_5
	
    #Using the relation Eu = dP/((1/2)*rho*v^2)
	dP = Eu*((rho*v_max**2)/2)		# Pressure drop per row [Pa]
	dP_total = dP*N					# Pressure drop across tube bundle [Pa]
	
	return dP_total

def HT_Zu(rho,Pr,Pr_w,a,b,d,geom,N,u,Re):
	'''
	Description:
		Calculate Nusselt number using the Zukauskas correlation
	Inputs:
		rho = density of the fluid [kg/m^3]
		Pr = Prandtl Number of the fluid
		Pr_w = Prandtl Number of the fluid at wall film temperature
		a = transverse pitch to diameter ratio
		b = longitudinal pitch to diameter ratio
		d = Tube Diameter [m]
		geom = inline or staggered
		N = Number of tubes
		u = free stream velocity without tubes [m/s]
		Re = Reynolds number calculated based on velcoity at narrowest point
	Outputs:
		Nu = Nusselt number
	Warnings:
		NOTE: Currently only applicable for staggered geometry 
	Validity:	
		Number of rows of tubes:
			Nr >= 2
		Reynolds number:
			100 ≤ Re ≤ 3x10e5
		For in-line tube arrangement: 	1.25 ≤ a ≤ 3.0; 1.2 ≤ b ≤ 3.0. 
		For staggered tube arrangement: 1.008 ≤ a ≤ 2.0; 1.008 ≤ b ≤ 2.0
			10 < Re < 2e6 ???
			0.7 < Pr < 500 ???
	Citation: Zhukauskas, A., R. Ulinskas., 1988, "Heat Transfer in Tube Banks in Crossflow,"
		Hemisphere Publishing Corporation. New York, NY. 1988.
		Zukauskas, A., 1972, "Heat Transfer from Tubes in Cross Flow," Adv. in Heat Trans, vol. 8, Academic Press, New York.
	'''
	
	if Re < 100:
		print('')
	c = np.sqrt(b**2+(a/2)**2)
	if c > (a+d)/2:
		Vmax= a/(a-d)*u
	else:
		Vmax= a/(2(b-d))*u
		
	if geom in ['inline','INLINE','Inline','square','SQUARE','Square']:  
		n = 0.36
		###Correction for less than 20 tube rows
		if N <=14:
			c2v = [.7,.8,.85,.9,.92,.94,.95,.96,.97,.97,.98,.98,.99,.99]
			c2 = c2v[N-1]
		### Coefficients based on range of Reynolds Number
		if Re > 1.6 and Re <= 100:
			c1 = 0.90
			m =  0.40
		elif Re > 100 and Re <= 1e3:
			c1=.52
			m= .5
		elif Re > 1e3 and Re <= 2e5:
			c1 = .27
			m =  .63
		elif Re > 2e5:
			c1 = .033
			m = .83
		else:
			print( 'Re out of Range')
			
	if geom in ['staggered','STAGGERED','Staggered','triangular','TRIANGULAR','Triangular']:
		###Correction for less than 20 tube rows
		if Re <= 1000 and N < 20:
			c2= 1-np.exp(-np.sqrt(3*N**(1/np.sqrt(2))))    
		if Re > 1000 and N < 20:
			c2 = 1-np.exp(-N**(1/np.sqrt(3)))
		### Coefficients based on range of Reynolds Number
		if Re <= 40 and Re > 1.6:
			c1 = 1.04
			m = .4
			n =.36
		elif Re > 40 and Re <= 1000:
			c1 = 0.52
			m = 0.5
			n = 0.36
		elif Re > 1000 and Re <= 2e5:
			if a/b < 2:
				c1 = 0.35*(a/b)**0.2
			else:
				c1 = 0.40
			m = 0.6
			n = 0.36
		elif Re > 2e5 and Re <= 2e6:
			c1 = 0.022
			m = 0.84
			n = 0.36
		else:
			print( 'Re out of Range')
    
	Nu = c1*c2*(Re**m)*(Pr**n)*(Pr/Pr_w)**.25
	
	return Nu
    
def dP_GG(rho,a,b,geom,N,u,Re,Return=""):
	'''
	Description:
		Calculate the pressure drop or total drag coefficient through a tube bundle using the Gaddis-Gnielinski correlation
	Inputs:
		rho = density of fluid [kg/m^3]
        a = transverse pitch ratio
        b = longitudinal pitch ratio
		geom = tube geometry (inline or staggered) [string]
		N = Number of tubes
		u = free stream velocity [m/s]
        Re = Reynolds number calculated based on velcoity at narrowest point and tube outer diameter
		Return = Blank returns pressure drop and "D_tot" returns total drag coefficient
	Outputs:
		dP_total = Pressure drop across tube bundle [Pa]
	Warnings:
		Not limited to specific pitch-diameter ratios within limits after a high enough Reynolds numbers. 
		The Reynolds number is based on the maximum velocity in the narrowest area between tubes.
	Validity:
		Reynolds number:
			1 =< a =< 3x10e5
		Number of rows of tubes:
			Nr >= 5
		In the range Re < 10E3: 
			For in-line tube arrangement: 	a x b = 1.25 x 1.25;   1.5 x 1.5;   2.0 x 2.0. 
			For staggered tube arrangement: a x b = 1.25 x 1.0825; 1.5 x 1.299; 1.768 x 0.884.
		In the range Re ≥ 10E3: 
			For in-line tube arrangement: 	1.25 ≤ a ≤ 3.0; 1.2 ≤ b ≤ 3.0. 
			For staggered tube arrangement: 1.25 ≤ a ≤ 3.0; 0.6 ≤ b ≤ 3.0; c ≥ 1.25.
	Citation: VDI. VDI Heat Atlas. Berlin, Heidelberg: VDI-Buch-Springer, 2010. Chapter L1
		Gaddis, E., Gnielinski, V., 1985, "Pressure drop in cross flow across tube bundles," Int. Chem. Eng., vol. 25(1), pp. 1-15
    '''
	if geom in ['inline','INLINE','Inline','square','SQUARE','Square']:
		u0 = u*(a/(a-1))
		D_lam= 280*np.pi*((b**(.5)-0.6)**2+0.75)/(a**(1.6)*(4*a*b-np.pi)*Re)	# Drag coefficient due to laminar flow
		if (N > 5 and N <= 10):
			f_nt= (1/a**2)*(1/N-.1)						# Coefficient for influence of inlet and outlet pressure losses
		elif N > 10:
			f_nt = 0									# Coefficient for influence of inlet and outlet pressure losses
		elif N <=5:
			sys.exit('Model only valid for greater than 5 tubes')
		f_ti = (0.22+1.2*(1-(0.94/b))**(0.6)/(a-0.85)**1.3)*10**(0.47*(b/a-1.5))+0.03*(a-1)*(b-1) # Coefficient for geometric arrangement factor  
		D_turb = f_ti/(Re**(.1*b/a))					# Drag coefficient due to turbulent flow
		D_tot = D_lam+(D_turb+f_nt)*(1-np.exp(-(Re+1000)/2000))
	if geom in ['staggered','STAGGERED','Staggered','triangular','TRIANGULAR','Triangular']:
		N = N-1
		if  b < .5*(2*a+1)**(.5):
			c=((a/2)**2+b**2)**(1/2)
			u0 = u*a/(2*(c-1))
			D_lam= 280*np.pi*((b**.5-0.6)**2+0.75)/((c**1.6)*(4*a*b-np.pi)*Re)
			if (N > 5 and N <= 10):
				f_nt = (2*(c-1)/a*(a-1)**2)*(1/N-.1)	# Coefficient for influence of inlet and outlet pressure losses
			elif N > 10:
				f_nt = 0								# Coefficient for influence of inlet and outlet pressure losses
		else:
			u0 = u*(a/(a-1))
			D_lam= 280*np.pi*((b**(.5)-0.6)**2+0.75)/(a**(1.6)*(4*a*b-np.pi)*Re)
			if (N > 5 and N <= 10):
				f_nt= (1/a**2)*(1/N-.1)					# Coefficient for influence of inlet and outlet pressure losses
			elif N > 10:
				f_nt = 0								# Coefficient for influence of inlet and outlet pressure losses
				
		f_ts = 2.5+1.2/(a-0.85)**1.08+0.4*(b/a-1)**3-0.01*(a/b-1)**3	# Coefficient for geometric arrangement factor  
		D_turb = f_ts/Re**0.25											# Drag coefficient due to turbulent flow
		D_tot = D_lam+(D_turb+f_nt)*(1-np.exp(-(Re+200)/1000))			# Total drag coefficient based on laminar and turbulent contributions

	if Return == "D_tot":
		return D_tot
	
	dP_total = .5*D_tot*(N*rho*u0**2) 					# Pressure drop across tube bundle [Pa]
   
	return dP_total
   
def HT_GG(rho,Pr,a,b,d,geom,N,u,Re):
	'''
	Description:
		Calculating the nusselt number using the Gnielinski and Martin Model based on the Leveque Analogy
	Inputs:
		rho = density of fluid [kg/m^3]
		Pr = Prandtl Number of the fluid
        a = transverse pitch ratio
        b = longitudinal pitch ratio
		d = tube diameter [m]
		geom = tube geometry (inline or staggered) [string]
		N = Number of tubes
		u = free stream velocity [m/s]
        Re = Reynolds number calculated based on velcoity at narrowest point and outside tube diameter
	Outputs:
		Nu = Nusselt Number
			Validity:	10 < Re < 2e6
				0.7 < Pr < 500
	Warnings:
		Ranges of validity have yet to be included in this script. Information soon to come.
	Citation: Martin, H., 2002, “The Generalized Lévêque Equation and its practical use for the prediction of heat and mass transfer rates from pressure drop,”
		Chem. Eng. Sci., vol. 57, pp. 3217-3223.
    '''
	
	if b > 1:
	   dh=((4*a/np.pi)-1)*d						# Hydraulic diameter [m]
	else:
	   dh=((4*a*b/np.pi)-1)*d					# Hydraulic diameter [m]
		
	if geom in ['inline','INLINE','Inline','square','SQUARE','Square']:
		L=b*d
		
	if geom in ['staggered','STAGGERED','Staggered','triangular','TRIANGULAR','Triangular']:
		c=((a/2)**2+b**2)**.5
		L=c*d

	xi = dP_GG(rho,a,b,geom,N,u,Re,"D_tot") 	# Drag Coefficient from Gaddis-Gnielinski Pressure Drop Model
	
	if Re > 2.5e5: 								 
		xi = xi*(1+(Re-2.5e5)/3.25e5)			# Correction made Holger Martin for use in Heat Transfer calculation
		
	xi_f = .5*xi 								# Total drag coefficient used to calculate drag coefficient due to friction ~.5
		
	NuPr = 0.404*(xi_f*Re**2*dh/L)**(1/3) #Nu/Pr**(1/3)
	Nu = NuPr*(Pr)**(1/3)
	
	return Nu
    
