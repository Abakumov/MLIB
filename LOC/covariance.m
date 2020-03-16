function M=covariance(x,y,z)

% calculate mean values to remove possible DC shifts

x = x - mean(x);
y = y - mean(y);
z = z - mean(z);

X = [x; y; z]; 

M = X*X'; 

M=M./length(x);