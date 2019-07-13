function t = Get_traveltime_3D_iCRS_par_LIA(MM, HH, t0, w, M, N, niter)

    % 3D iCRS formula
    % parabolic reflector
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
    
    for i=1:fold
        dxs = MM(:,i) - HH(:,i);
        dxg = MM(:,i) + HH(:,i);
        
        dxr = [0; 0];
        zeta = t0+dxr'*K*dxr;
        
        ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
        tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
        
        for iter=1:niter
             
            dxsg  = (tg0*dxs + ts0*dxg)/(ts0 + tg0);
          
            b = M*dxsg; 
            A = (zeta + w'*dxsg)/t0*K + M;

            dxr = A\b; 
           
            zeta = t0+dxr'*K*dxr;
        
            ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
            tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
            
            t(iter,i) = ts0 + tg0;
        end
    end
end