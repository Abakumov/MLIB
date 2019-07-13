function t = Get_traveltime_3D_CRS(MM, HH, t0, w, M, N)

    % Hyperbolic 3D CRS formula
    % see equation 1.19
    % size(w) = 2 x 1
    % size(M) = 2 x 2
    % size(N) = 2 x 2
    % size(MM) = 2 x fold
    % size(HH) = 2 x fold
    % size(t)  = 1 x fold

    fold = size(MM,2); 
    t = zeros(1, fold); 
    for i=1:fold
        m = MM(:,i); 
        h = HH(:,i); 
        t2 = (t0 + w'*m)^2 + 2*t0*(m'*N*m + h'*M*h);
        t(i) = sqrt(t2);
    end

  