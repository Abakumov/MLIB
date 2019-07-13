function T = Get_traveltime_2D_nCRS_PS(MM, HH, model)

% Calculate n-CRS type traveltime approximation for converted waves
% Abakumov Ivan
% 10th April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

%T = zeros(size(M));  

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
Vp   = model(4);                 % km/s 
Vs   = model(5);                 % km/s  
gamma = Vp/Vs; 
upsilon = 2/(gamma+1);

% define Vps (see eq. 4.3)
Vps = 2*Vp*Vs/(Vp+Vs); 
t0 = 2*Rnip/Vps; 
w = 2*sin(alpha)/Vps; 
M = (cos(alpha))^2/(Vps*Rnip); 
N = (cos(alpha))^2/(Vps*Rn); 

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

Fs = sqrt((t0+w*(MM-HH)).^2 + 2*t0*N*(MM-HH).^2);
Fg = sqrt((t0+w*(MM+HH)).^2 + 2*t0*N*(MM+HH).^2);
Fsg = 2*t0*gamma*(M-N)*upsilon^2*HH.^2;
T = sqrt(((Fs + gamma*Fg)/(1+gamma)).^2 + Fsg); 