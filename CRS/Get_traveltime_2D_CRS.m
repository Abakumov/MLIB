function T = Get_traveltime_2D_CRS(MM, HH, model)

% Calculate CRS traveltime approximation
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
%         m = MM(i,j); 
%         h = HH(i,j); 
% 
%         % CRS stacking operator (see eq. 1.9)
%         T(i,j) = sqrt((t0+w*m).^2 + 2*t0*(N*m.^2 + M*h.^2));
%     end
% end


% CRS stacking operator (see eq. 1.9)
T = sqrt((t0+w*MM).^2 + 2*t0*(N*MM.^2 + M*HH.^2));
