function x = get_AWG_noise(t,sigma)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 26th August 2019

    x = sigma*randn(size(t));
end

