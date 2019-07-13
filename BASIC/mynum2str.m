function outstr = mynum2str(innum)

% MYNUM2STR - returns string with additional zeros:
% e.g. if innum = 23 result is '000023'
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019


innum = round(abs(innum));

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
% 
% if innum < 10
%     outstr = ['000' num2str(innum)];
% else
%     if innum < 100
%         outstr = ['00' num2str(innum)];
%     else 
%         if innum < 1000
%             outstr = ['0' num2str(innum)];
%         else
%             if innum < 10000
%                 outstr = num2str(innum);
%             else
%                 disp('Program input is integer positive laues less that 1000');
%             end
%         end
%     end
% end