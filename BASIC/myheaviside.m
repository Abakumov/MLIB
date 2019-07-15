function y = myheaviside(x)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019

    y = zeros(size(x));
    ind = (x>=0);
    y(ind) = x(ind);

end

