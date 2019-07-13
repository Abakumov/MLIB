function t = Get_traveltime_3D_iCRS_el_LIA(MM, HH, t0, v0, w, M, N, niter)
    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru
    % 3D iCRS formula
    % ellipsoidal reflector
    % Linearized iterative approach
    % size(w) = 2 x 1
    % size(M) = 2 x 2
    % size(N) = 2 x 2
    % size(MM) = 2 x fold
    % size(HH) = 2 x fold
    % size(t)  = 1 x fold

    fold = size(MM,2); 
    t = zeros(niter, fold); 
    K = M/(M-N)*N;
    tau = 1/(v0^2*sqrt(det(K)));
    
    for i=1:fold
        dxs = MM(:,i) - HH(:,i);
        dxg = MM(:,i) + HH(:,i);
        
        dxr = [0; 0];
        khi = sqrt(1 - dxr'*K*dxr/tau); 
        zeta = t0+2*tau*(1-khi);
        
        ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
        tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
        
        for iter=1:niter
             
            dxsg  = (tg0*dxs + ts0*dxg)/(ts0 + tg0);

            b = M*dxsg; 
            A = (zeta + w'*dxsg)/(t0*khi)*K + M;

            dxr = A\b; 
            khi = sqrt(1 - dxr'*K*dxr/tau); 
            zeta = t0+2*tau*(1-khi);
            
            ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
            tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
            
            t(iter,i) = ts0 + tg0;
        end
    end
end