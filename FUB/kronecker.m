function delta = kronecker(i,j)

% Kronecker delta function
% A function of two variables, usually just positive integers. 
% The function is 1 if the variables are equal, and 0 otherwise
%
% *Author*: Ivan Abakumov
%
% *Publication date*: 22nd September 2017
%    
% *E-mail*: abakumov_ivan@mail.ru

if i==j
    delta=1; 
else
    delta=0;
end
