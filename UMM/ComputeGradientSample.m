function grad = ComputeGradientSample(x)

grad = zeros(size(x));

dx = 1e-6; 

for i = 1:length(grad) 
    xp = x; 
    xp(i) = xp(i) + dx;  
    xm = x; 
    xm(i) = xm(i) - dx;  
    fp = ComputeObjectiveSample(xp); 
    fm = ComputeObjectiveSample(xm);
    grad(i) = (fp - fm)/2/dx; 
end

