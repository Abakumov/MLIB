function cray = Get_smodel_central_ray(smodel)
% Get central ray in simplified model
%
% Abakumov Ivan
% 28th August 2016

    %% Make central ray in the simplified model

    X0 = smodel.x0; 
    R = smodel.R; 
    RNIP = smodel.RNIP; 

    cray.zz_src = 0:1:RNIP; 
    cray.xx_src = zeros(size(cray.zz_src)); 
    cray.yy_src = zeros(size(cray.zz_src));
    cray.xx = zeros(size(cray.xx_src)); 
    cray.yy = zeros(size(cray.yy_src)); 
    cray.zz = zeros(size(cray.zz_src)); 

    for i=1:size(cray.xx_src,1); 
        for j=1:size(cray.xx_src,2)
            x_src = [ cray.xx_src(i,j);... 
                      cray.yy_src(i,j);... 
                      cray.zz_src(i,j)];
            x = R*x_src + X0;
            cray.xx(i,j) = x(1); 
            cray.yy(i,j) = x(2); 
            cray.zz(i,j) = x(3); 
        end
    end
end
