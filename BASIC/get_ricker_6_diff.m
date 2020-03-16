 function y = get_ricker_6_diff(t,tau,fc)

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
    x3 = (7+2*a).*exp(a);
    x4 = (9+2*a).*exp(a);
    x5 = (11+2*a).*exp(a);
    x6 = (13+2*a).*exp(a);

    y = x6.*(a1).^6 + 15*x5.*a1.^4*a2 + 45*x4.*a1.^2*a2^2 + 15*x3*(a2)^3; 
    y = y/(15*7*a2^3); % normalize
    
    %x(abs(s)>4.) = 0;    
end

