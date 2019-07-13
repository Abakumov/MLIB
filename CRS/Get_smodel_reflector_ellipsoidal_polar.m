function ref = Get_smodel_reflector_ellipsoidal_polar(smodel, THETA, PHI)
% Get ellipsoidal reflector (polar angles) in the simplified model
%
% Abakumov Ivan
% 28th August 2016

    %% Make ellipsoidal reflector in the simplified model
    
    RNIP   = smodel.RNIP; 
    KR_src = smodel.KR_src; 
    R      = smodel.R; 
    X0     = smodel.x0;
    
    %
    Rel = 1/sqrt(det(KR_src)); 
    RN = RNIP + Rel; 
    
    % see equations (B.7)
    dphi = asin(KR_src(1,2)/sqrt(KR_src(1,1)*KR_src(2,2)));
    A_src = sqrt(Rel*KR_src(2,2)/det(KR_src));
    B_src = sqrt(Rel*KR_src(1,1)/det(KR_src));
    
    ref.theta = THETA; 
    ref.phi   = PHI; 
    
    ref.xx_src = zeros(size(ref.theta)); 
    ref.yy_src = zeros(size(ref.theta)); 
    ref.zz_src = zeros(size(ref.theta)); 
    
    % see equations (B.6)
    for i=1:size(ref.theta,1); 
        for j=1:size(ref.theta,2)
            theta = ref.theta(i,j);
            phi   = ref.phi(i,j);
            ref.xx_src(i,j) = A_src*sin(theta)*cos(phi); 
            ref.yy_src(i,j) = B_src*sin(theta)*sin(phi-dphi); 
            ref.zz_src(i,j) = RN - Rel*cos(theta);
        end
    end
    
    ref.xx = zeros(size(ref.theta)); 
    ref.yy = zeros(size(ref.theta)); 
    ref.zz = zeros(size(ref.theta)); 
    
    for i=1:size(ref.xx_src,1); 
        for j=1:size(ref.xx_src,2)
            x_src = [ ref.xx_src(i,j);...
                      ref.yy_src(i,j);...
                      ref.zz_src(i,j)];
            x = R*x_src + X0;
            ref.xx(i,j) = x(1); 
            ref.yy(i,j) = x(2); 
            ref.zz(i,j) = x(3); 
        end
    end
end