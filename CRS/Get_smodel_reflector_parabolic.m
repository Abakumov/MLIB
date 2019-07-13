function ref = Get_smodel_reflector_parabolic(smodel, XX, YY)
% Get parabolic reflector in the simplified model
%
% Abakumov Ivan
% 28th August 2016

    %% Make parabolic reflector in the simplified model
    
    RNIP   = smodel.RNIP; 
    KR_src = smodel.KR_src; 
    R      = smodel.R; 
    X0     = smodel.x0;
    
    ref.xx_src = XX; 
    ref.yy_src = YY; 
 
    % see eq. 2.14 
    ref.zz_src = RNIP + 0.5*(KR_src(1,1)*ref.xx_src.^2 ...
                         +   KR_src(2,2)*ref.yy_src.^2 ...
                         + 2*KR_src(1,2)*ref.xx_src.*ref.yy_src); 
    
    ref.xx = zeros(size(ref.xx_src)); 
    ref.yy = zeros(size(ref.yy_src)); 
    ref.zz = zeros(size(ref.zz_src));              

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
