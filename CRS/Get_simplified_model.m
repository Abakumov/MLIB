function smodel = Get_simplified_model(CRS_param)
% This function make simplified model from CRS stacking parameters
%
% Abakumov Ivan
% 28th August 2016

    x0 = CRS_param.x0;
    t0 = CRS_param.t0; 
    v0 = CRS_param.v0; 
    w  = CRS_param.w; 
    M  = CRS_param.M; 
    N  = CRS_param.N; 

    if length(M(:)) == 9 
        [ alpha, beta, KNIP3, KN3 ] = my_A3P(v0, w, M, N);
        KNIP = KNIP3(1:2,1:2); 
        KN   = KN3(1:2,1:2); 

    else
        [ alpha, beta, KNIP, KN ] = my_A2P(v0, w, M, N);
    end
    
    %% Equalize non-diagonal elements
    KNIP(1,2) = (KNIP(1,2) + KNIP(2,1))/2; 
    KNIP(2,1) = KNIP(1,2); 
    KN(1,2) = (KN(1,2) + KN(2,1))/2; 
    KN(2,1) = KN(1,2); 

    %% Find eigen values and eigen vectors of matrix KNIP

    % Rotation matrix from the general Cartesian system L 
    % to the standard ray-centered system
    R1 = Get_Rmatrix_3x3(alpha, beta);    
 
    % Additional rotation for angle delta (eq. 2.6)
    [R_delta, KNIP_diag] = eig(KNIP);                                  
    R2 =  zeros(3,3); 
    R2(3,3) = 1; 
    R2(1:2,1:2) = R_delta; 
 
    % Rotation matrix from the general Cartesian system L
    % to the special ray-centered system L' (eq. 2.6)
    % x_src = R'*(x-x0);                   % from L to L'
    % x = R*x_src + x0;                    % from L' to L
    R = R1*R2;           
 
    % Curvatures of NIP and normal waves in the system L'
    KNIP_src = R_delta'*KNIP*R_delta; 
    KN_src   = R_delta'*KN*R_delta; 
 
    % Diagonal elements of matrix KNIP (eq. 2.8)
    knip11 = KNIP_diag(1,1); 
    knip22 = KNIP_diag(2,2); 
  
    % Auxiliary anisotropic medium (eq. 2.13)
    A11 = 2*v0/(t0*knip11); 
    A22 = 2*v0/(t0*knip22); 
    A33 = v0^2; 
 
    % Depth of the NIP point (eq. 2.11)
    RNIP = t0*v0/2; 

    % Curvature of the reflector
    KR_src = inv(inv(KN_src) - inv(KNIP_src)); 

    %% Make simplified model
    smodel.x0  = x0; 
    smodel.R   = R; 
    smodel.A11 = A11;
    smodel.A22 = A22; 
    smodel.A33 = A33; 
    smodel.RNIP = RNIP; 
    smodel.KR_src  = KR_src; 
end
