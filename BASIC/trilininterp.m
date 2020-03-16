function interp = trilininterp(initial, fragm)

% Ivan Abakumov

% for description of realized here trilinear interpolation
% see for details http://en.wikipedia.org/wiki/Trilinear_interpolation
% (formula in Russian version)

inisize = size(initial);
ndims = length(inisize);
newsize = (inisize-1)*fragm + 1;
interp = zeros(newsize);

  % 1D & 2D case
  
if ndims == 2
  fragm2 = fragm^2;
  for i=1:newsize(1)
    for j=1:newsize(2);
      im = 1 + floor((i-1)/fragm);  
      jm = 1 + floor((j-1)/fragm);  
      ip = 1 + ceil((i-1)/fragm);   
      jp = 1 + ceil((j-1)/fragm);   
      ii = i - 1 - (im-1)*fragm;        
      jj = j - 1 - (jm-1)*fragm;        
   
      interp(i, j) =  initial(im, jm)/fragm2*(fragm-ii)*(fragm-jj)+...
                      initial(im, jp)/fragm2*(fragm-ii)*jj+...
                      initial(ip, jm)/fragm2*ii*(fragm-jj)+...
                      initial(ip, jp)/fragm2*ii*jj;
    end
  end  
  
  % 3D case  
elseif ndims == 3
  fragm3 = fragm^3;
  for i=1:newsize(1)
    for j=1:newsize(2)
      for k=1:newsize(3)   
                                      % k: 1 2 3 4 5 6  
        im = 1 + floor((i-1)/fragm);  
        jm = 1 + floor((j-1)/fragm);  
        km = 1 + floor((k-1)/fragm);  % 1 1 1 1 1 2
        ip = 1 + ceil((i-1)/fragm);   
        jp = 1 + ceil((j-1)/fragm);   
        kp = 1 + ceil((k-1)/fragm);   % 1 2 2 2 2 2
        ii = i - 1 - (im-1)*fragm;        
        jj = j - 1 - (jm-1)*fragm;        
        kk = k - 1 - (km-1)*fragm;   % 0 1 2 3 4 0 

        interp(i, j, k) =  initial(im, jm, km)/fragm3*(fragm-ii)*(fragm-jj)*(fragm-kk)+...
                           initial(im, jm, kp)/fragm3*(fragm-ii)*(fragm-jj)*kk+...
                           initial(im, jp, km)/fragm3*(fragm-ii)*jj*(fragm-kk)+...
                           initial(im, jp, kp)/fragm3*(fragm-ii)*jj*kk+...
                           initial(ip, jm, km)/fragm3*ii*(fragm-jj)*(fragm-kk)+...
                           initial(ip, jm, kp)/fragm3*ii*(fragm-jj)*kk+...
                           initial(ip, jp, km)/fragm3*ii*jj*(fragm-kk)+...
                           initial(ip, jp, kp)/fragm3*ii*jj*kk;
      end
    end
  end
  
  % more dimensions
elseif ndims>3
    disp('Error, algorithm is designed for arrays with numer of dimentions less then 3')
end