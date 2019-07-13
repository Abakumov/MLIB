function [ w, M, N ] = my_P2A(v0, alpha, beta, KNIP, KN)
%my_P2A get stacking parameters from wavefield attributes 
% see equation 1.20

    R = Get_Rmatrix_2x2(alpha, beta);

    w = - 2/v0*sin(alpha)*[cos(beta); sin(beta)];
    
    M = 1/v0*R*KNIP*R';
    
    N = 1/v0*R*KN*R';
end

