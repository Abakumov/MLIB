function Ur = get_subsidence_radial_Geertsma(D,R,r)
%
% see formula 7 in Geertsma, J. (1973). Land subsidence above compacting oil 
% and gas reservoirs. Journal of Petroleum Technology, 25(06), 734-744.
 
nu = D/R; 
rho = r/R; 

Ur = zeros(size(rho)); 

for i=1:length(rho)

    f_B = @(a,nu,rho) exp(-nu*a).*besselj(1,a).*besselj(1,a*rho);
    Ur(i) = integral(@(a)f_B(a,nu,rho(i)),0,Inf);

end

