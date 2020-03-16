function Ur = get_subsidence_radial_Dyskin(D,R,r)
%
% see formula 33 in Dyskin, A., Pasternak, E., and Shapiro, S. A. (2020).  
% Subsidence, uplift and shift due to fluid extraction and production in a 
% finite reservoir. International Journal for Numerical and Analytical Methods in
% Geomechanics (in review)
 
Ur = zeros(size(r)); 

Ur(r<R) = -3/4 * ((r(r<R)*D^2*R^2)/((R^2 + D^2)^(5/2)) ...
- (5*r(r<R).^3*D^2*R^2*(4*D^2 - 3*R^2))/(8*(R^2+D^2)^(9/2)));

Ur(r>=R) = -3/4 * ((r(r>=R)*D^2*R^2)./((r(r>=R).^2 + D^2).^(5/2)) ...
    - (5*r(r>=R)*D^2*R^4.*(4*D^2-3*r(r>=R).^2))./(8*(r(r>=R).^2+D^2).^(9/2))); 


