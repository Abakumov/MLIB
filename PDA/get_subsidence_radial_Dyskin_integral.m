function Ur = get_subsidence_radial_Dyskin_integral(D,R,r)
%
% see formula 20 in Dyskin, A., Pasternak, E., and Shapiro, S. A. (2020).  
% Subsidence, uplift and shift due to fluid extraction and production in a 
% finite reservoir. International Journal for Numerical and Analytical Methods in
% Geomechanics (in review)
 
nu = D/R; 
rho = r/R; 

Ur = zeros(size(rho)); 

for i=1:length(rho)

    I11 = @(a,nu,rho) a.*exp(-nu*a).*besselj(1,a).*besselj(1,a*rho);
    Ur(i) = nu*integral(@(a)I11(a,nu,rho(i)),0,Inf);

end

