function [f, amp]=myfft(trace, dt)
% Frendly Fourier transform 
% Abakumov Ivan
% 9 Sept 2014 



Fs = 1/dt;
L = length(trace);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(trace,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
amp = 2*abs(Y(1:NFFT/2+1));

% Plot single-sided amplitude spectrum.
%plot(f,amp) 
%title('Single-Sided Amplitude Spectrum of y(t)')
%xlabel('Frequency (Hz)')
%ylabel('|Y(f)|')