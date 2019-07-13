function T = Get_traveltime_2D_iCRS(MM, HH, model)

% Calculate i-CRS traveltime approximation
% Abakumov Ivan
% 7st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg
% Original code was provided by Benjamin Schwarz

T = zeros(size(MM));  
niso=1;

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
v0   = model(4);                 % km/s 

t0 = 2*Rnip/v0; 
vnmo=sqrt(2.*v0*Rnip/(t0*(cos(alpha))^2));
lambda=1./sqrt(1.+vnmo^2/v0^2*(sin(alpha))^2);
V=vnmo*lambda;
xc=-Rn*sin(alpha)/(cos(alpha))^2*lambda^2;
H=v0*Rn/(vnmo*(cos(alpha))^2)*lambda^2;
R=(v0*Rn/(vnmo*(cos(alpha))^2)-vnmo*t0/2)*lambda;

for j=1:size(MM,2); 
    tantheta0=(MM(1,j)-xc)/H;
    tantheta=tantheta0;
    costheta=1./sqrt(1.+tantheta^2);
    sintheta=tantheta*costheta;
    for i=1:size(MM,1); 
        for iso=1:niso
            a=MM(i,j)-xc-R*sintheta;
            b=(H-R*costheta)^2;
            ts=(a-HH(i,j))^2+b;
            tg=(a+HH(i,j))^2+b;
            tantheta=tantheta0+HH(i,j)*(ts-tg)/(H*(ts+tg+2.*sqrt(ts*tg)));
            costheta=1./sqrt(1.+tantheta^2);
            sintheta=tantheta*costheta;
        end
        a=MM(i,j)-xc-R*sintheta;
        b=(H-R*costheta)^2;
        ts=sqrt((a-HH(i,j))^2+b);
        tg=sqrt((a+HH(i,j))^2+b);
        T(i,j)=1./V*(ts+tg);
    end
end