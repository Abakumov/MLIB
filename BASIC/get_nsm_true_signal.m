function x = get_nsm_true_signal(t,t0,fc,damp)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019


    x = myheaviside(t-t0).*sin(2*pi*fc*(t-t0) ).*exp(-damp*(t-t0));


end

