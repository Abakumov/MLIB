function [maxi, maxj, maxk, maxx, maxy, maxz] = get_max_of_pdf(F, G) 

    % Apply spline interpolation near maximum of the pdf function and update
    % maximum location
    %
    % input: F - 3D array, G - GridClass
    % output - (maxi maxj maxk) - index of nearest to maximum gridpoint
    %          (maxx maxy maxz) - interpolated coordinates of maximum
    %
    % *Author*: Ivan Abakumov
    %
    % *Publication date*: 25th October 2017
    %    
    % *E-mail*: ivan.abakumov@fu-berlin.de

    [~, maxind] = max(F(:));
    [maxi,maxj,maxk] = ind2sub(size(F),maxind);

    indx = max((maxi-2),1):min((maxi+2),G.nx);
    indy = max((maxj-2),1):min((maxj+2),G.ny);
    indz = max((maxk-2),1):min((maxk+2),G.nz);
    
    V = F(indx,indy,indz); 

    [x, y, z] = ndgrid(G.xx(indx),G.yy(indy),G.zz(indz)); 

    indxx = min(G.xx(indx)):(G.dx/10):max(G.xx(indx));
    indyy = min(G.yy(indy)):(G.dy/10):max(G.yy(indy));
    indzz = min(G.zz(indz)):(G.dz/10):max(G.zz(indz));

    [xq, yq, zq] = ndgrid(indxx, indyy, indzz);  

    % spline interpolation requires at least 4 points in each dimension
    % cubic interpolation requires at least 2 points in each dimension
    if min(size(V))>3
        Vq = interpn(x,y,z,V,xq,yq,zq, 'spline');
    else
        Vq = interpn(x,y,z,V,xq,yq,zq, 'cubic');
    end

    [~, maxind] = max(Vq(:));
    [vqx,vqy,vqz] = ind2sub(size(Vq),maxind);

    maxx = indxx(vqx); 
    maxy = indyy(vqy); 
    maxz = indzz(vqz);
    
end
