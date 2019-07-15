function c = myxcorr(a,b)
%MYXCORR cross-correlation function estimates.
%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019
    
    corrLength=length(a)+length(b)-1;
    c=fftshift(ifft(fft(a,corrLength).*conj(fft(b,corrLength))));
    

end

