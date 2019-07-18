function outstr = mynum2str(varargin)

% MYNUM2STR - returns string with additional zeros:
% e.g. if innum = 23 result is '000023'
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 18th July 2019

if nargin == 1
    innum = varargin{1};
    strlength = 5; 
end
if nargin == 2
   innum = varargin{1};  
   strlength = varargin{2}; 
end




innum = round(abs(innum));

if strlength == 5

    if innum < 10
        outstr = ['0000' num2str(innum)];
    else
        if innum < 100
            outstr = ['000' num2str(innum)];
        else 
            if innum < 1000
                outstr = ['00' num2str(innum)];
            else
                if innum < 10000
                    outstr = ['0' num2str(innum)];
                else
                    disp('Program input is integer positive laues less that 10000');
                end
            end
        end
    end
elseif strlength == 4
    if innum < 10
        outstr = ['000' num2str(innum)];
    else
        if innum < 100
            outstr = ['00' num2str(innum)];
        else 
            if innum < 1000
                outstr = ['0' num2str(innum)];
            else
                if innum < 10000
                    outstr = num2str(innum);
                else
                    disp('Program input is integer positive laues less that 1000');
                end
            end
        end
    end
else
    disp('Variable strlength not properly set (possible values strlength=4 or = 5)');
    
end