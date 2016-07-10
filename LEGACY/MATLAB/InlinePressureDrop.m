%Original Author: Julien Clayton
%Last Modified: 7/9/2016
%Most Recent Author: Jonah Haefner
%All Authors: Jonah Haefner, Julien Clayton, Lane Carasik, ...

%%%Units%%%
%a = dimensionless transverse pitch
%b = dimensionless longitudinal pitch
%v = the free stream fluid velocity in m/s 
%rho = density in g/cm^3
%u = dynamic viscosity in pa*s
%N = number of rows in the tube bundle 
%D_tube = The diamter of the tubes in the bundle in m
%Will work for values of b = 1.25, 1.5, or 2
%Reynolds number needs to be below 150000

function [dP_total, Re, v_mean,Eu] = InlinePressureDrop(a,b,v,rho,u,N,D_tube, Re)
x = double((a-1)./(b-1));
rho = rho*1000;               %g/cm^3 ---> kg/m^3(density)
v_mean = v.*(a./(a-1));   %mean flow velocity in min XSection of tube bank (Zukauska's book assumes avg flow velocity is approx. max velocity)
%Re = double((rho.*D_tube.*v_mean)./u);  %MinXS Reynolds number (if not already known)

%coefficients for Euler number (calculated using the power series)
%coefficients come from here: http://www.thermopedia.com/content/1211/#TUBE_BANKS_CROSSFLOW_OVER_FIG2

% These are valid for Reynolds numbers between 2E3 and 2E6 %
if Re > 2E3
c_0 = [.267, .235, .247];      %b = 1.25, 1.5, and 2
c_1 = [.249E4, .197E4, -.595]; %b = 1.25, 1.5, and 2
c_2 = [-.927E7, -.124E8, .15]; %b = 1.25, 1.5, and 2
c_3 = [.1E11, .312E11, -.137]; %b = 1.25, 1.5, and 2
c_4 = [0, -.274E14, .396];     %b = 1.25, 1.5, and 2
elseif Re <= 800 && Re > 3
c_0 = [.272, .263, .188];      %b = 1.25, 1.5, and 2
c_1 = [.207E3, .867E2, 56.6]; %b = 1.25, 1.5, and 2
c_2 = [.102E3, -.202, -646]; %b = 1.25, 1.5, and 2
c_3 = [-.286E3, 0, 6010]; %b = 1.25, 1.5, and 2
c_4 = [0, 0, -18300];     %b = 1.25, 1.5, and 2
elseif Re <= 2E3 && Re > 800
c_0 = [.272, .263, .247];      %b = 1.25, 1.5, and 2
c_1 = [.207E3, .867E2, 56.6-(4.766E-2*Re)]; %b = 1.25, 1.5, and 2
c_2 = [.102E3, -.202, .15]; %b = 1.25, 1.5, and 2
c_3 = [-.286E3, 0, -.137]; %b = 1.25, 1.5, and 2
c_4 = [0, 0, .396];     %b = 1.25, 1.5, and 2
%There may be errors when b=2 and Re < 800 but didn't input for code
%simplicity.
end
%%%%%%%%% Correction Factors %%%%%%%%%% 

% Non-rectangular bundle % (valid up to Reynolds numbers of 150000)

k_1 = 1;

if x ~= 1.00
    
    if (Re >= 1000 && Re < 10000)
    
        k_1 = .9849.*x.^(-.8129);
    
        else if (Re >= 10000 && Re < 70000)  
    
            k_1 = .9802.*x.^(-.7492);
    
            else if (Re >= 70000 && Re < 150000)
                    
                k_1 = .988.*x.^(-.6388);
             
                end 
                
              end
            
    end
    
end 

% Entry losses % 

%Entry loss coefficients (Re > 1E4 but < 1E6) 
el_1 = [1.9, 1.1, 1, 1, 1, 1, 1];

%Entry loss coefficients (Re > 1E6) 
el_2 = [2.7, 1.8, 1.5, 1.4, 1.3, 1.2, 1.2];
 

if (b == 1.25)
    
    i = 1;
    
    else if (b == 1.5)
        
        i = 2;
        
        else if (b == 2)
            
            i = 3;
            
            end 
        
          end 
    
end 


k_3 = 1;

if (N < 7 && N > 0 && Re < 1E6 && Re > 1E2)
       
    k_3 = el_1(N);
    
    else if (N < 7 && N > 0 && Re >= 1E6)
    
        k_3 = el_2(N);
       
        else if (Re <= 1E2)
                %disp('Reynolds number needs to be greater than 100')
                %break 
            end 
            
        end 
                
    
end

%Power series for Euler number per row. From same website as above.
Eu_p = (c_0(i)./Re.^0)+(c_1(i)./Re.^1)+(c_2(i)./Re.^2)+(c_3(i)./Re.^3)+(c_4(i)./Re.^4);


    
Eu = Eu_p.*k_3.*k_1;
%Using the relation Eu = dP/((1/2)*rho*v^2)
dP = Eu.*((rho.*v_mean.^2)./2); %Pressure drop per row 
dP_total = dP.*N/1000; %pressure drop across 10 rows expressed in kPa


end

