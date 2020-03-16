 function y = get_ricker_2_diff(t,tau,fc)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% Ricker wavelet in time
% fc is peak frequency in Hz
% t is time array

    tcut = 1.5/fc;    
    s = (t-tau-tcut).*fc;
    a = - pi^2*s.^2; 
    a1 =  2*(pi*fc)^2*(t-tau-tcut);
    a2 =  - 2*(pi*fc)^2;
    x0 = (1+2*a).*exp(a);
    x1 = (3+2*a).*exp(a);
    x2 = (5+2*a).*exp(a);
    
    y = x2.*(a1).^2 + x1*a2; 
    y = y/(-6*pi^2*fc^2); % normalize
    
    %x(abs(s)>4.) = 0;    
end
