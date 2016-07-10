%Original Author: Jonah Haefner
%Last Modified: 10/27/2015
%Most Reecent Author: Jonah Haefner
%References: Julien Clayton, Lane Carasik, ...

%%%Units%%%
%a = dimensionless transverse pitch
%b = dimensionless longitudinal pitch
%v = the free stream fluid velocity in m/s 
%rho = density in kg/m^3
%u = dynamic viscosity in pa*s
%N = number of rows in the tube bundle 
%D_tube = The diamter of the tubes in the bundle in m
%Will work for values of b = 1.25, 1.5, or 2
%Reynolds number needs to be below 150000 and greater than E2


%coefficients for Euler number (calculated using the power series)
%coefficients come from here: http://www.thermopedia.com/content/1211/#TUBE_BANKS_CROSSFLOW_OVER_FIG2

% These are valid for Reynolds numbers between 2E3 and 2E6 %
function [dP_total, Re, v_mean,Eu] = StaggeredPressureDrop(a,b,v,rho,u,N,D_tube, Re)
x = double((a)./(b));      %Used for correction factor k_1 which I don't have data for.
%Mean Velocity Calculations. Found on same website as above
if a <= (2*b^2 -.5)
v_mean = v.*(a./(a-1));
elseif a > (2*b^2 -.5)
v_mean = v.*(a./(sqrt(4*b^2+a^2)-2))
end

% Coefficients, c_i, to generate pressure drop coefficients for equilateral triangle banks.
if Re < 1E3
c_0 = [.795,.683,.343];      %b = 1.25, 1.5, and 2
c_1 = [.247E3,.111E3,.303E3]; %b = 1.25, 1.5, and 2
c_2 = [.335E3,-.973E2,-.717E5]; %b = 1.25, 1.5, and 2
c_3 = [-.155E4,-.426E3,.88E7]; %b = 1.25, 1.5, and 2
c_4 = [.241E4,.574E3,-.38E9];     %b = 1.25, 1.5, and 2
elseif Re >= 1E3 && Re < 1E4
c_0 = [.245,.203,.343];      %b = 1.25, 1.5, and 2
c_1 = [.339E4, .248E4, .303E3 ]; %b = 1.25, 1.5, and 2
c_2 = [-.984E7, -.758E7, -.717E5]; %b = 1.25, 1.5, and 2
c_3 = [.132E11, .104E11, .88E7]; %b = 1.25, 1.5, and 2
c_4 = [-.599E13, -.482E13, -.38E9];     %b = 1.25, 1.5, and 2
elseif Re >= 1E4
c_0 = [.245, .203, .162  ];   %b = 1.25, 1.5, and 2
c_1 = [.339E4, .248E4, .181E4 ]; %b = 1.25, 1.5, and 2
c_2 = [-.984E7, -.758E7, .792E8]; %b = 1.25, 1.5, and 2
c_3 = [.132E11, .104E11, -.165E13]; %b = 1.25, 1.5, and 2
c_4 = [-.599E13, -.482E13, .872E16];     %b = 1.25, 1.5, and 2
end 
% Assigning correct values to each tube spacing
if (b == 1.25)
    
    i = 1;
    
    else if (b == 1.5)
        
        i = 2;
        
        else if (b == 2)
            
            i = 3;
            
            end 
        
          end 
    
end 

%%%%%%%%% Correction Factors %%%%%%%%%% 

%k_1 is the influence of pitch ratio

        %Will add later if data can be found

%k_2 =  Influence of Temperature on Fluid properties. Neglected 

% k_3 is Entry length effects 

%Entry loss coefficients (Re < 1E2 but < 1E4)
el_1 = [1.4, 1.3, 1.2, 1.1, 1, 1, 1];
%Entry loss coefficients (Re > 1E4 but < 1E6) 
el_2 = [1.1, 1.05, 1, 1, 1, 1, 1];
%Entry loss coefficients (Re > 1E6)
el_3 = [.25, .45, .6, .65, .7, .75, .8];

k_3 = 1; %Setting for tubes > 7
if (N < 7 && N > 0 && Re < 1E4 && Re > 1E2)
        k_3 = el_1(N);
elseif (N < 7 && N > 0 && Re < 1E6 && Re >= 1E4)
      k_3 = el_2(N);
elseif (N < 7 && N > 0 && Re > 1E6)
      k_3 = el_3(N);
end

%Power series for Euler number per row. From same website as above.
Eu_p = (c_0(i)./Re.^0)+(c_1(i)./Re.^1)+(c_2(i)./Re.^2)+(c_3(i)./Re.^3)+(c_4(i)./Re.^4);


    
Eu = Eu_p.*k_3; %.*k_1;
%Using the relation Eu = dP/((1/2)*rho*v^2)
dP = Eu.*((rho.*v_mean.^2)./2); %Pressure drop per row 
dP_total = dP.*N; %pressure drop across 10 rows


end
