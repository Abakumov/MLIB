function xmax = fithyp(imax, y1, y2, y3)

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019

% fit hyperbola for given three point

% Y    y1      y2      y3
% X  imax-1   imax   imax+1
% x0 - point where function value is maximum
% 
%
%     fit y=y0 + a*(x-x0)**2
%     let x0=imax+d
%     then
%     y1=y0 + a*(1+2d+d*d)
%     y2=y0 + a*(d*d)
%     y3=y0 + a*(1-2d+d*d)
%     y2-y1 = -a*(1+2d)
%     y2-y3 = -a*(1-2d)
%     r=(y2-y1)/(y2-y3) = (1+2d)/(1-2d)
%     d=(r-1)/(r+1)/2

if y1==y3
    xmax=imax;
elseif y1==y2
    xmax=imax-0.5;
elseif y3==y2
    xmax=imax+0.5;
else
	ratio=(y2-y1)/(y2-y3);
	d=(ratio-1)/(ratio+1)/2;
	xmax=imax+d;
end
