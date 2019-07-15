function output=add_time_shift_fft(trace,tshift,Fs)
% Adds time shift "tshift" to trace

%% Author: Ivan Abakumov
% Freie Universität Berlin
% E-mail: abakumov_ivan@mail.ru
% Publication date: 12th July 2019


L = length(trace);
[Sx, Sy] = size(trace);
% check that trace is a vector
if min(Sx, Sy)~=1
    disp('Error, trace shold be vector!');
end
% chech if Sy = L and transport
if Sy==L
  trace=trace';
end

i = sqrt(-1);
NFFT = 2^nextpow2(L);
f = Fs/2*linspace(0,1,NFFT/2+1);
Ftrace = fft(trace,NFFT);
Ftrace(1:NFFT/2+1) = Ftrace(1:NFFT/2+1).*exp(-i*2*pi*f'*tshift);
Ftrace((NFFT/2+2):NFFT)=conj(Ftrace(NFFT/2:-1:2));
n1=ceil(abs(tshift)*Fs);
traceShifted=real( ifft(Ftrace) );
traceShifted=traceShifted(1:L);

output = zeros(size(trace));
if tshift>0
    output(n1:L)=traceShifted(n1:L);
else
    output(1:L-n1)=traceShifted(1:L-n1);
end

output=reshape(output, Sx, Sy);

