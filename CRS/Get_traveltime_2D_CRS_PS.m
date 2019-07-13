function T = Get_traveltime_2D_CRS_PS(MM, HH, model)

% Calculate DSR-PS traveltime approximation
% Abakumov Ivan
% 1st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

%T = zeros(size(MM));  

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
Vp   = model(4);                 % km/s 
Vs   = model(5);                 % km/s  
gamma = Vp/Vs; 
sigma = (gamma-1)/(gamma+1);
upsilon = 2/(gamma+1);

% define Vps (see eq. 4.3)
Vps = 2*Vp*Vs/(Vp+Vs); 
t0 = 2*Rnip/Vps; 
w = 2*sin(alpha)/Vps; 
M = (cos(alpha))^2/(Vps*Rnip); 
N = (cos(alpha))^2/(Vps*Rn); 

% for i=1:size(MM,1); 
%     for j=1:size(MM,2); 
%         
%         % define gamma-CMP coordinates (see eq. D.31)
%         m = MM(i,j) + sigma*HH(i,j); 
%         h = upsilon*HH(i,j);
% 
%         % CRS-PS stacking operator (see eq. 4.12)
%         T(i,j) = sqrt((t0+w*m)^2 + 2*t0*(N*m^2 + gamma*M*h^2));
%     end
% end


% define gamma-CMP coordinates (see eq. D.31)
m = MM + sigma*HH; 
h = upsilon*HH;

% CRS-PS stacking operator (see eq. 4.12)
T = sqrt((t0+w*m).^2 + 2*t0*(N*m.^2 + gamma*M*h.^2));


