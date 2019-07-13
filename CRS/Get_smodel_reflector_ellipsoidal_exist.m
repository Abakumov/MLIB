function ind = Get_smodel_reflector_ellipsoidal_exist(smodel, XX, YY, ZZ)
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
   
    ind = (XX == 0); 
                             
                                     
    for i=1:size(XX,1); 
        for j=1:size(XX,2)
            for k=1:size(XX,3)
                x = [ XX(i,j,k);...
                      YY(i,j,k);...
                      ZZ(i,j,k)];
                      x_src = R'*(x - X0);
            ind(i,j,k) =   (RN - x_src(3))^2 + Rel*(KR_src(1,1)*x_src(1)^2 ...
                                                +   KR_src(2,2)*x_src(2)^2 ...
                                                + 2*KR_src(1,2)*x_src(1)*x_src(2)) < Rel^2 ;
            end
        end
    end
end
    
    
    
    
    
    
    
    
    
    
    