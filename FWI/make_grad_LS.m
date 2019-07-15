function [ Grad ] = make_grad_LS( G, U, Q )
%MAKE_GRAD_LS compute gradient for Least Squares functional
% See R.-E. Plessix, (2006) A review of the adjoint-state method for computing
% the gradient of a functional with geophysical applications
% Geophys. J. Int. V. 167, pp. 495–503, formula (36)
% NOTE: formula gives dJ/dm where m = 1/v^2
% => dJ/dv = (dJ/dm)*(dm/dv) 
%   Abakumov Ivan
%   University of Hamburg
%   25.06.2015
    [nx, nz, nt, nshot] = size(U); 

    Grad = zeros(nx, nz);
    for s=1:nshot
        grad = zeros(nx, nz); 
        for i=1:nx
            for j=1:nz
                q = squeeze(Q(i,j,end:-1:1,s)); 
                u = squeeze(U(i,j,:,s)); 
                u(2:end-1) = diff(u,2);                                             % second derivative
                grad(i,j) = - sum(q.*u)/G.dt;                                       % devide by (dt)^2 in diff and multiply by dt in integral 
            end
        end
        Grad = Grad + grad; 
    end
end

