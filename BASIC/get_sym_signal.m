function x = get_sym_signal(t,tau,fc,damp)
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

    x = cos(2*pi*fc*(t-tau)).*exp(-damp*abs(t-tau));

end

