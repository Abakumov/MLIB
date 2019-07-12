function [p, ux, uz, uxx, uzz, uxz, uzx] = exel2d(rho,vp,vs,xs,zs,xl,zl,wf,fr,stype)
    nfr = length(fr); 
    ifr = find(fr>0);
    nsh = length(xs);    
    om  = fr(ifr).*(2*pi);
    
    if isvector(xl)
        nxl = length(xl);
        nzl = length(zl);
        p  = single(zeros(nfr,nxl,nzl,nsh));
        ux = p; uz = p; uxx = p; uzz = p; uxz = p; uzx = p;
        for ish=1:nsh
            for izl=1:nzl
                for ixl=1:nxl                               
                    x  = xl(ixl)-xs(ish);            
                    z  = -(zl(izl)-zs(ish));
                    [pf,uxf,uzf,uxxf,uzzf,uxzf,uzxf] = wavefield(x,z,rho,vp,vs,om,wf(ifr),stype);
                    p  (ifr,ixl,izl,ish) = pf;
                    ux (ifr,ixl,izl,ish) = uxf;
                    uz (ifr,ixl,izl,ish) = uzf;
                    uxx(ifr,ixl,izl,ish) = uxxf;
                    uzz(ifr,ixl,izl,ish) = uzzf;
                    uxz(ifr,ixl,izl,ish) = uxzf;
                    uzx(ifr,ixl,izl,ish) = uzxf;
                end
            end
        end
    elseif size(xl,1)==2
        nrec = size(xl,2);
        p  = single(zeros(nfr,nrec,nsh));
        ux = p; uz = p; uxx = p; uzz = p; uxz = p; uzx = p;
        for ish=1:nsh
            for irec=1:nrec                    
                x  = xl(1,irec)-xs(ish);
                z  = -(xl(2,irec)-zs(ish));
                [pf,uxf,uzf,uxxf,uzzf,uxzf,uzxf] = wavefield(x,z,rho,vp,vs,om,wf(ifr),stype);
                p  (ifr,irec,ish) = pf;
                ux (ifr,irec,ish) = uxf;
                uz (ifr,irec,ish) = uzf;
                uxx(ifr,irec,ish) = uxxf;
                uzz(ifr,irec,ish) = uzzf;
                uxz(ifr,irec,ish) = uxzf;
                uzx(ifr,irec,ish) = uzxf;
            end
        end
    end
end

function [p,ux,uz,uxx,uzz,uxz,uzx] = wavefield(x,z,rho,vp,vs,om,wf,stype)
    x2 = x^2;
    z2 = z^2;
    r2 = x2+z2;
    r  = sqrt(r2);
    vp2 = vp^2;
    vs2 = vs^2;
    
    if r<eps(r)
        p = 0.*om;
        ux = p; uz = p; uxx = p; uzz = p; uxz = p; uzx = p;
        return;
    end
    
    rkp  = om.*r./vp;
    irkp = rkp.^-1;
    
    rks  = om.*r./vs;
    irks = rks.^-1;
    
    h0p = besselh(0,1,rkp);
    h1p = besselh(1,1,rkp);
    h2p = besselh(2,1,rkp);
    
    h0s = besselh(0,1,rks);
    h1s = besselh(1,1,rks);
    h2s = besselh(2,1,rks);
    
    gp1 = h1p.*irkp./vp2;
    gp2 = h0p./vp2;
    
    gs1 = h1s.*irks./vs2;
    gs2 = h0s./vs2;
    
    gp1x = -h2p.*(x/r2/vp2);
    gp1z = -h2p.*(z/r2/vp2);
    gp2x = -h1p.*rkp.*(x/r2/vp2);
    gp2z = -h1p.*rkp.*(z/r2/vp2);
    
    gs1x = -h2s.*(x/r2/vs2);
    gs1z = -h2s.*(z/r2/vs2);
    gs2x = -h1s.*rks.*(x/r2/vs2);
    gs2z = -h1s.*rks.*(z/r2/vs2);
    
    C = wf.*(1i/4/rho);
    
    switch stype        
        case 1
            ux  = C.*(gp1 + (x2/r2).*gp2 + ((x2-z2)/r2).*gs1 + (z2/r2).*gs2);
            uz  = (C.*(x*z/r2)).*(2.*gp1 - gp2 - 2.*gs1 + gs2);
            uzz = -C.*((x*z/r2).*(2.*gp1z - gp2z - 2.*gs1z + gs2z) + (x*(x2-z2)/r2^2).*(2.*gp1 - gp2 - 2.*gs1 + gs2));
            uzx = C.*((x*z/r2).*(2.*gp1x - gp2x - 2.*gs1x + gs2x) + (z*(z2-x2)/r2^2).*(2.*gp1 - gp2 - 2.*gs1 + gs2));
            uxz = -C.*(gp1z - (2*z*x2/r2^2).*gp2 + (x2/r2).*gp2z - (4*z*x2/r2^2).*gs1 + ((x2-z2)/r2).*gs1z + (2*z*x2/r2^2).*gs2 + (z2/r2).*gs2z);
            uxx = C.*(gp1x - (2*x*z2/r2^2).*gp2 + (x2/r2).*gp2x + (4*x*z2/r2^2).*gs1 + ((x2-z2)/r2).*gs1x - (2*x*z2/r2^2).*gs2 + (z2/r2).*gs2x);
            p   = -(rho*vp2).*(uxx + uzz);
        case 2
            ux  = (C.*(x*z/r2)).*(2.*gp1 - gp2 - 2.*gs1 + gs2);
            uz  = C.*(gp1 + (z2/r2).*gp2 + ((z2-x2)/r2).*gs1 + (x2/r2).*gs2);
            uxx = C.*((x*z/r2).*(2.*gp1x - gp2x - 2.*gs1x + gs2x) + (z*(z2-x2)/r2^2).*(2.*gp1 - gp2 - 2.*gs1 + gs2));
            uxz = -C.*((x*z/r2).*(2.*gp1z - gp2z - 2.*gs1z + gs2z) + (x*(x2-z2)/r2^2).*(2.*gp1 - gp2 - 2.*gs1 + gs2));
            uzx = C.*(gp1x - (2*x*z2/r2^2).*gp2 + (z2/r2).*gp2x - (4*x*z2/r2^2).*gs1 + ((z2-x2)/r2).*gs1x + (2*x*z2/r2^2).*gs2 + (x2/r2).*gs2x);
            uzz = -C.*(gp1z - (2*z*x2/r2^2).*gp2 + (z2/r2).*gp2z + (4*z*x2/r2^2).*gs1 + ((z2-x2)/r2).*gs1z - (2*z*x2/r2^2).*gs2 + (x2/r2).*gs2z);
            p   = -(rho*vp2).*(uxx + uzz);
        otherwise
            error('Unknown/non-implemented source type!');
    end
    
end