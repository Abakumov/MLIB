function tt = get_tti(tti,G,acq)
%GET_TTI find traveltimes at the position of receivers
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019

inisize = size(tti);
ndims = length(inisize);


% 2D case
  
if ndims == 2
   
    tt = zeros(size(acq.rx));
    acq.grx = x2grid_floor(acq.rx,G.x0,G.dx,G.nx);
    acq.grz = x2grid_floor(acq.rz,G.z0,G.dz,G.nz);
    

    for recenum = 1:length(acq.rx)
    
        xm = acq.grx(recenum);  
        xp = min(acq.grx(recenum)+1,G.nx);   
        xx = acq.rx(recenum) - G.xx(xm);        
        zm = acq.grz(recenum);  
        zp = min(acq.grz(recenum)+1,G.nz);   
        zz = acq.rz(recenum) - G.zz(zm);                
        dx = G.dx;
        dz = G.dz;

        tt(recenum) =  tti(xm, zm)*(dx-xx)*(dz-zz)/(dx*dz)+...
                                      tti(xm, zp)*(dx-xx)*zz/(dx*dz)+...
                                      tti(xp, zm)*xx*(dz-zz)/(dx*dz)+...
                                      tti(xp, zp)*xx*zz/(dx*dz);
    end  
 
% more dimensions
elseif ndims>3
    disp('Error, algorithm is designed for arrays with numer of dimentions less then 3')
end



end

