function f = ComputeObjectiveSample(x)

% Test the "lbfgs" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.

alpha = x(1); 
delta = x(2); 
epsilon = x(3); 
dtheta = x(4); 

Sample = MLD('/home/ivan/Desktop/MLIB/UMM/TempSample.mat'); 

ym = get_Vqp_VTI_weak(Sample,alpha,delta,epsilon,dtheta);
yo = Sample.Vqp;
dy = ym - yo;    
m = length(dy);
f = sum(dy.^2)/m/2;