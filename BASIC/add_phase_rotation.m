function output = add_phase_rotation(trace, phaselag)
%% ADD_PHASE_ROTATION.m
% adds phase shift "phaselag" to the "trace" using DFT
% 
% Note: size(trace) = [L, 1]

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

% add phase shift using DFT
i = sqrt(-1);
NFFT = 2^nextpow2(L);
Ftrace = fft(trace,NFFT);
pFtrace(1:NFFT/2+1) = Ftrace(1:NFFT/2+1)*exp(i*phaselag);
pFtrace((NFFT/2+2):NFFT)=conj(pFtrace(NFFT/2:-1:2));
Strace=real(ifft(pFtrace));
output=Strace(1:L);
output=reshape(output, Sx, Sy);
