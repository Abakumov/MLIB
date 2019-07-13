function t = Get_traveltime_3D_nCRS(MM, HH, t0, w, M, N)
    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru
    
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
        Fs = (t0 + w'*(m-h))^2 + 2*t0*(m'-h')*N*(m-h);
        Fg = (t0 + w'*(m+h))^2 + 2*t0*(m'+h')*N*(m+h);
        t(i) = sqrt( 0.25*(sqrt(Fs) + sqrt(Fg))^2 + 2*t0*h'*(M-N)*h);
    end
end