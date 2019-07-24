function Grad = TOMOcomputeGradient (dv)

while exist('TOMOGrad.mat', 'file') ~= 2
    pause(1); 
    disp('Attantion with gradient!')
end

Grad = MLD('TOMOGrad.mat'); 
    
! rm -f TOMOGrad.mat
  
  
  
