function Uz = get_subsidence_vertical_Geertsma(D,R,r)
%
% see formula 6 in Geertsma, J. (1973). Land subsidence above compacting oil 
% and gas reservoirs. Journal of Petroleum Technology, 25(06), 734-744.
 
nu = D/R; 
rho = r/R; 

Uz = zeros(size(rho)); 

for i=1:length(rho)

    f_A = @(a,nu,rho) exp(-nu*a).*besselj(1,a).*besselj(0,a*rho);
    Uz(i) = integral(@(a)f_A(a,nu,rho(i)),0,Inf);

end

