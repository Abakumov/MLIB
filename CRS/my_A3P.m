function [ alpha, beta, KNIP, KN ] = my_A3P(v0, w, M, N)
%my_A3P get wavefield attributes from the stacking parameters
% see equation 1.20

    beta = atan(w(2)/w(1)); 
    
    alpha = asin(norm(w(1:2))*v0/2*sign(-w(1)));
   
    R = Get_Rmatrix_3x3(alpha, beta);
    
    KNIP = v0*R'*M*R;
    
    KN   = v0*R'*N*R;
    
end