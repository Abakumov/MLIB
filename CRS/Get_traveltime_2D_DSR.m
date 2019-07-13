function T = Get_traveltime_2D_DSR_PS(MM, HH, model)

% Calculate DSR-PS traveltime approximation
% Abakumov Ivan
% 1st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

%T = zeros(size(MM)); 

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
V0   = model(4);                 % km/s 

t0 = 2*Rnip/V0; 
w = 2*sin(alpha)/V0; 
M = (cos(alpha))^2/(V0*Rnip); 
N = (cos(alpha))^2/(V0*Rn); 

% for i=1:size(MM,1); 
%     for j=1:size(MM,2); 
%         
%         % define CMP coordinates 
%         m = MM(i,j); 
%         h = HH(i,j);
%         
%         % DSR-PS stacking operator (see eq. 4.10)
%         Ss = sqrt((t0+w*(m-h))^2 + 2*t0*(N*m^2 - 2*N*m*h + M*h^2));
%         Sg = sqrt((t0+w*(m+h))^2 + 2*t0*(N*m^2 + 2*N*m*h + M*h^2));
%         
%         T(i,j) = (Ss + Sg)/2; 
%     end
% end

Ss = sqrt((t0+w*(MM-HH)).^2 + 2*t0*(N*MM.^2 - 2*N*MM.*HH + M*HH.^2));
Sg = sqrt((t0+w*(MM+HH)).^2 + 2*t0*(N*MM.^2 + 2*N*MM.*HH + M*HH.^2));
        
T = (Ss + Sg)/2; 



