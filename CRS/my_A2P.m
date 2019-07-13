function [ alpha, beta, KNIP, KN ] = my_A2P(v0, w, M, N)
    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru
    %my_A2P get stacking parameters from wavefield attributes 
    % see equation 1.20
    % Abakumov, I. (2017). Systematic analysis of double-square-root-based stacking operators

    beta = atan(w(2)/w(1)); 
    
    alpha = asin(norm(w)*v0/2*sign(-w(1)));
   
    R = Get_Rmatrix_2x2(alpha, beta);
    
    Q = inv(R); 
    
    KNIP = v0*Q*M*Q';
    
    KN   = v0*Q*N*Q';
    
end