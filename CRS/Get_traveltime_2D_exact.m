function [T, Xr] = Get_traveltime_2D_exact(M, H, model)

% Calculate exact traveltime for PS or PP wave
% reflected from the circle by dihotomia
% Abakumov Ivan
% 1st April 2016
% abakumov_ivan@mail.ru
% University of Hamburg


T = zeros(size(M));  
Xr = zeros(size(M));  

% alpha = model(1);                % deg 
% Rnip = model(2);                 % km 
% Rn   = model(3);                 % km
% Vp   = model(4);                 % km/s 
% Vs   = model(5);                 % km/s   

for i=1:size(M,1); 
    for j=1:size(M,2); 
        xg = M(i,j) + H(i,j); 
        xs = M(i,j) - H(i,j); 
        thetar = linspace(-pi, pi, 21); 
        for iter = 1:10
            [xref, zref] = Get_2D_circle(model, thetar); 
            [t, ind] = min(Get_2D_traveltime(xs, xg, xref, zref, model)); 
            minind = max(ind-1, 1); 
            maxind = min(ind+1, 21); 
            thetar = linspace(thetar(minind), thetar(maxind), 21); 
        end 
        T(i,j) = t; 
        Xr(i,j) = thetar(11); 
    end
end

end
 

function t = Get_2D_traveltime(xs, xg, xref, zref, model)
% see eq. D.9
Vp   = model(4);                 % km/s 
Vs   = model(5);                 % km/s 

t = sqrt((xs-xref).^2+zref.^2)/Vp + sqrt((xg-xref).^2+zref.^2)/Vs;
end


function [xref, zref] = Get_2D_circle(model, thetar)
% see eq. D.11
alpha = model(1);                % deg 
Rnip = model(2);                 % km 
Rn   = model(3);                 % km
R = Rn-Rnip; 

xref = -Rn*sin(alpha) + R*sin(alpha - thetar); 
zref = Rn*cos(alpha) - R*cos(alpha - thetar);
end
