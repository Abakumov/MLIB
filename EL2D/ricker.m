function a = ricker(fpeak,t)

%RICKER 
% Ricker wavelet in time
% fpeak is peak frequency in Hz
% t is time array
    tcut = 1.5/fpeak;    
    s = (t-tcut).*fpeak;
    a = (1-2*pi^2*s.^2).*exp(-pi^2*s.^2);
    a(abs(s)>4.) = 0;    
end

