function Uz = get_subsidence_vertical_Dyskin_integral(D,R,r)
%
% see formula 21 in Dyskin, A., Pasternak, E., and Shapiro, S. A. (2020).  
% Subsidence, uplift and shift due to fluid extraction and production in a 
% finite reservoir. International Journal for Numerical and Analytical Methods in
% Geomechanics (in review)
 
nu = D/R; 
rho = r/R; 

Uz = zeros(size(rho)); 

for i=1:length(rho)
    
    I00 = @(a,nu,rho) exp(-nu*a).*besselj(1,a).*besselj(0,a*rho);
    I01 = @(a,nu,rho) a.*exp(-nu*a).*besselj(1,a).*besselj(0,a*rho);

    Uz(i) = (integral(@(a)I00(a,nu,rho(i)),0,Inf) + nu*integral(@(a)I01(a,nu,rho(i)),0,Inf));


end

