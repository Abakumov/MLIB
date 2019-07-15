function y = heaviside(x)
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019


    y = zeros(size(x));
    y(x>0) = x;

end

