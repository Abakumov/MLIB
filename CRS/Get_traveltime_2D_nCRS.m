function T = Get_traveltime_2D_nCRS(M, H, model)

% Calculate n-CRS traveltime approximation
% See paper Fomel, Kazinnin, 2013, non-hyperbolic Common Reflection surface
% Geophysical prospecting, 61, pp 21-27
% Abakumov Ivan
% 1st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

%T = zeros(size(M));  

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
V0   = model(4);                 % km/s 
        
t0 = 2*Rnip/V0; 
a1 = 2*sin(alpha)/V0; 
a2 = 2*t0*(cos(alpha))^2/(V0*Rn); 
b2 = 2*t0*(cos(alpha))^2/(V0*Rnip); 
c = 2*b2 + a1^2 - a2; 

% for i=1:size(M,1); 
%     for j=1:size(M,2); 
%         
%         m = M(i,j); 
%         h = H(i,j); 
%         
%         F=(t0+a1*m)^2+a2*m^2;
%         F1=(t0+a1*(m-h))^2+a2*(m-h)^2;
%         F2=(t0+a1*(m+h))^2+a2*(m+h)^2;
%         T(i,j)=sqrt((F+c*h^2+sqrt(F1*F2))/2.);
%     end
% end

F=(t0+a1*M).^2+a2*M.^2;
F1=(t0+a1*(M-H)).^2+a2*(M-H).^2;
F2=(t0+a1*(M+H)).^2+a2*(M+H).^2;
T=sqrt((F+c*H.^2+sqrt(F1.*F2))/2.);