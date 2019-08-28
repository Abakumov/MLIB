function vartau = get_CRB(SNR,beta)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

% Cramer Rao Bound



vartau = 1./beta.^2./SNR;


end

