function [f, Psd] = get_doublesided_PSD(trace,Fs)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 29th August 2019

% Calculate double sided power spectral density

nfft = 2^nextpow2(length(trace));
Psd = abs(fft(trace,nfft)).^2/length(trace)/Fs;
Psd = Psd(1:length(Psd)/2);         % Doublesided
%PSDdb = 10*log10(Psd);
f = Fs/2*linspace(0,1,nfft/2);






end

