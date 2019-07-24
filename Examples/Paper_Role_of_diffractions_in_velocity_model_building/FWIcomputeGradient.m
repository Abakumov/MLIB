function Grad = FWIcomputeGradient (dv)

while exist('Grad.mat', 'file') ~= 2
    pause(1); 
    disp('Attantion with gradient!')
end

Grad = MLD('Grad.mat'); 
delete('Grad.mat');
  
  
  
