function t = Get_traveltime_3D_iCRS_el_TIA(MM, HH, t0, v0, w, M, N, niter)

    % Abakumov Ivan
    % University of Hamburg
    % e-mail: abakumov_ivan@mail.ru
    % 3D iCRS formula
    % ellipsoidal reflector
    % Trigonometric iterative approach
    % size(w) = 2 x 1
    % size(M) = 2 x 2
    % size(N) = 2 x 2
    % size(MM) = 2 x fold
    % size(HH) = 2 x fold
    % size(t)  = 1 x fold

    fold = size(MM,2); 
    t = zeros(niter, fold); 
    K = M/(M-N)*N;
    tau = 1/(v0*sqrt(det(K*v0)));
    
    Ae = sqrt(tau*K(2,2)/det(K));
    Be = sqrt(tau*K(1,1)/det(K));
    sdphi = K(1,2)/sqrt(K(1,1)*K(2,2));
    dphi  = asin(sdphi);
    cdphi = cos(dphi);  
    
    for i=1:fold
        dxs = MM(:,i) - HH(:,i);
        dxg = MM(:,i) + HH(:,i);
        
        theta = 0; 
        phi   = 0;
        
        dxr = [Ae*sin(theta)*cos(phi) ; Be*sin(theta)*sin(phi-dphi)];
        zeta = t0+2*tau*(1-cos(theta));
        
        ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
        tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
        
        for iter=1:niter
             
            dxsg  = (tg0*dxs + ts0*dxg)/(ts0 + tg0);

            b = M*dxsg; 
            delta = M*(dxr-dxsg); 
           
            kappa = (zeta + w'*dxsg)*tau/t0;
            omega = M(1,1)*(Ae*cos(phi))^2 +2*M(1,2)*(Ae*cos(phi))*(Be*sin(phi+dphi)) +M(2,2)*(Be*sin(phi+dphi))^2;

            if iter == 1
                phi   = atan(    delta(2)*Be*cdphi/( delta(1)*Ae - delta(2)*Be*sdphi) ); 
            else
                dxrdphi  = [-Ae*sin(theta)*sin(phi) ;  Be*sin(theta)*cos(phi-dphi)];
                dxrddphi = [-Ae*sin(theta)*cos(phi) ; -Be*sin(theta)*sin(phi-dphi)];
                
                 phi = phi - delta'*dxrdphi/(delta'*dxrddphi + dxrdphi'*M*dxrdphi);
            end
            
            theta = atan( ( b(1)*Ae*cos(phi) + b(2)*Be*sin(phi-dphi) )/(kappa+omega*cos(theta)) );
 
            dxr = [Ae*sin(theta)*cos(phi) ; Be*sin(theta)*sin(phi-dphi)];
            zeta = t0+2*tau*(1-cos(theta));
 
            ts0 = 0.5*sqrt( (zeta+w'*dxs)^2 + 2*t0*(dxs-dxr)'*M*(dxs-dxr) );
            tg0 = 0.5*sqrt( (zeta+w'*dxg)^2 + 2*t0*(dxg-dxr)'*M*(dxg-dxr) );
            
            t(iter,i) = ts0 + tg0;
        end
    end
end