function a = ricker8(fpeak,t)
%RICKER
% 8th derivative of Ricker wavelet in time
% fpeak is peak frequency in Hz
% t is time array
    tcut = 1.5/fpeak;    
    s = (t-tcut).*fpeak;
    %a = (1-2*pi^2*s.^2).*exp(-pi^2*s.^2);
    %f = sym('(1-2*w^2)*exp(-w^2)');
    %simplify(diff(f,8))
    w = pi*s; 
    a = -16*exp(-w.^2).*(32*w.^10 - 720*w.^8 + 5040*w.^6 - 12600*w.^4 + 9450*w.^2 - 945);
    a(abs(s)>4.) = 0;    
end

