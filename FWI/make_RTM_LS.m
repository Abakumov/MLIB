function [ RTM ] = make_RTM_LS( G, U, Q )
%MAKE_RTM compute RTM image
% NOTE: formula gives dJ/dm where m = 1/v^2
% => dJ/dv = (dJ/dm)*(dm/dv) 
%   Abakumov Ivan
%   University of Hamburg
%   25.06.2015

    [nx, nz, nt, nshot] = size(U); 

    RTM = zeros(nx, nz);
    for s=1:nshot
        rtm = zeros(nx, nz); 
        for i=1:nx
            for j=1:nz
                q = squeeze(Q(i,j,end:-1:1,s)); 
                u = squeeze(U(i,j,:,s)); 
                rtm(i,j) = - sum(q.*u)*G.dt;
            end
        end
        RTM = RTM + rtm; 
    end
end