function Uz = get_subsidence_vertical_Dyskin(D,R,r)
%
% see formula 32 in Dyskin, A., Pasternak, E., and Shapiro, S. A. (2020).  
% Subsidence, uplift and shift due to fluid extraction and production in a 
% finite reservoir. International Journal for Numerical and Analytical Methods in
% Geomechanics (in review)

Uz = zeros(size(r)); 

Uz(r<R) = 1  - D^3/(R^2 + D^2)^1.5 ...
    - (15*r(r<R).^2*D^3*R^2)/(4*(R^2 + D^2)^(7/2));

Uz(r>=R) = (3*D^3*R^2)./(2*(r(r>=R).^2 + D^2).^(5/2)) ...
    - (15*D^3*R^4*(2*D^2 - 5*r(r>=R).^2))./(16*(r(r>=R).^2 + D^2).^(9/2)); 

