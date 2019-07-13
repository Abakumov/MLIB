function [T, Xr] = Get_traveltime_2D_sixth_order(M, H, model)

% Calculate exact traveltime for PS or PP wave
% reflected from the circle by solution of sixth-order equation
% Abakumov Ivan
% 1st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg

T = zeros(size(M));  
Xr = zeros(size(M));  

alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
Vp   = model(4);                 % km/s 
Vs   = model(5);                 % km/s  
R = Rn-Rnip; 
xc = -Rn*sin(alpha);              % center of circle 
zc = Rn*cos(alpha);
gamma = Vp/Vs; 
D = (1-gamma^2)/(1+gamma^2);    % see eq. E5
Vp = Vp/R; 
Vs = Vs/R; 

for i=1:size(M,1); 
    for j=1:size(M,2); 
        m = M(i,j); 
        h = H(i,j); 

        % Make model scaling (unit circle, xs = -hnew, xg = + hnew)
        m1 = (xc-m)/R; 
        m2 = zc/R; 
        h  = h/R; 

        xs = -h; 
        xg =  h; 

        % Coefficients of 6th order equation
        A1 = m1^2 + m2^2 + 1 + h^2 -2*m1;
        A2 = m1^2 + m2^2 + 1 + h^2 +2*m1;
        B1 = 4*D*h*m2 + 4*m1*m2;
        B2 = 4*h*m2 + 4*D*m1*m2;
        C1 = 8*D*h*m1 + 4*h^2 - 2*m2^2 + 4*m1^2;
        C2 = 8*h*m1 + D*(4*h^2 - 2*m2^2 + 4*m1^2);

        b6 = (2*h*(m1-1) - D*A1)*m2^2;
        b5 = 2*h*(m1-1)*B1 - A1*B2 - 4*D*m2^3;
        b4 = 2*h*(m1+1)*m2^2 + 2*h*(m1-1)*C1 - A2*D*m2^2 - 4*m2*B2 - A1*C2;
        b3 = -8*m2*(4*h*m1 + D*(4*m1^2 - m2^2));
        b2 = 2*h*(m1-1)*m2^2 + 2*h*(m1+1)*C1 - A1*D*m2^2 + 4*m2*B2 - A2*C2;
        b1 = -2*h*(m1+1)*B1 + A2*B2 - 4*D*m2^3;
        b0 = (2*h*(m1+1)-D*A2)*m2^2;

        u = roots([b6,b5,b4,b3,b2,b1,b0]); 
        thetar = 2*atan(u); 
        % find only real roots
        eps = 1e-8; 
        ind = (imag(thetar)<eps); 
        thetar = real(thetar(ind)); 

        xref = m1 + cos(thetar); 
        zref = m2 + sin(thetar); 

        [t, ind] = min(sqrt((xs-xref).^2+zref.^2)/Vp + sqrt((xg-xref).^2+zref.^2)/Vs);
        T(i,j) = t; 
        Xr(i,j) = thetar(ind); 
    end
end
