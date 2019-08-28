function beta = get_spectral_variance(signal,dt)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% Calculate spectral variance of the signal (angular bandwidth of the pulse) 
% in time domain

beta = sqrt(sum((diff(signal,1)/dt).^2)/sum(signal.^2));


end

