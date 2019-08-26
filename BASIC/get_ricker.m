function x = get_ricker(t,tau,fc)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% Ricker wavelet in time
% fc is peak frequency in Hz
% t is time array

    tcut = 1.5/fc;    
    s = (t-tau-tcut).*fc;
    x = (1-2*pi^2*s.^2).*exp(-pi^2*s.^2);
    %x(abs(s)>4.) = 0;    
end

