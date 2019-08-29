function vartau = get_SZZB(SNR,Td,Ta)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% Shapiro bound type II


vartau =    Ta^2/12*erfc(sqrt(SNR/4)) ...
    + (Td^2/(4*pi^2)./SNR).*gammainc(SNR/4,3/2) ...
    - (Td^2/(4*pi^2)./SNR).^1.5.*(32/(3*Ta*sqrt(2*pi))).*gammainc(SNR/4,2);





