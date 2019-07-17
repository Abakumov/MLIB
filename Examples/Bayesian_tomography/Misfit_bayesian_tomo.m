function F = Misfit_bayesian_tomo(X,G,XX,ZZ,acq,texact,sigma)
% Computes L2 misfit 
    nshot = length(acq.sx); 
    nrece = length(acq.rx); 
    
    Gold=G.oldGrid;
    velmod = voronoi2velmod(X,G,XX,ZZ);
    tt = zeros(nshot,nrece); 
    
    for s=1:nshot
        S = [acq.sx(s) acq.sz(s)];
        tti  = FSM2D(Gold,S,velmod); 
        for r=1:nrece
            tt(s,r) = tti(acq.grx(r), acq.grz(r));  
        end
    end
    F = exp(-rms(tt(:)-texact(:))/sigma);
end

