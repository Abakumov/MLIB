function [T, Xref] = Get_model_exact_traveltime(Xs, Xg, model)
    % Calculate exact traveltime in analytic models
    % Abakumov Ivan
    % 27th August 2016
    % abakumov_ivan@mail.ru
    % Hamburg University

    Get_model_parameters;
    
    T = zeros(1,size(Xs,2)); 
    Xref = zeros(size(Xs)); 
    
    for i = 1:size(Xs,2); 
        xs = Xs(:,i);
        xg = Xg(:,i);
    
        %% User parameters
        eps = 10^-12;                       % Accuracy of traveltime (sec)
        nxsteps = 50;                       % Number of search steps in x direction 
        nysteps = 50;                       % Number of search steps in y direction

        %% No changes below
        xmin = G.x0; 
        xmax = G.mx; 
        ymin = G.y0;
        ymax = G.my; 
        xstep = (xmax-xmin)/nxsteps; 
        ystep = (ymax-ymin)/nysteps; 
        told  = 100; 
        dt    = 100; 
        nxsteps = 50; 
        nysteps = 50;

        while abs(dt) > eps;  % 10+2, 10-4, 10-6, 10-8, 10-11, 10-13, 10-14

            [XX, YY] = meshgrid(xmin:xstep:xmax, ymin:ystep:ymax);

            % find zref
            [zref, ind] = Get_model_surface(XX,YY,model); 
            XX = XX(ind); 
            YY = YY(ind); 
            ZZ = zref(ind); 

            % find traveltime
            TT = Get_model_traveltime(xs,xg,XX,YY,ZZ,model);

            % find min traveltime
            [tnew,  I] = min(TT(:));

            % update search conditions
            xmin = max(G.x0, XX(I)-(2+0.1*rand(1))*xstep);  
            xmax = min(G.mx, XX(I)+(2+0.1*rand(1))*xstep); 
            ymin = max(G.y0, YY(I)-(2+0.1*rand(1))*ystep);
            ymax = min(G.my, YY(I)+(2+0.1*rand(1))*ystep); 
            xstep = (xmax-xmin)/nxsteps; 
            ystep = (ymax-ymin)/nysteps; 
            dt = told - tnew; 
            told  = tnew; 
        end
        T(i) = tnew; 
        Xref(:,i) = [XX(I); YY(I); ZZ(I)]; 
    end
end