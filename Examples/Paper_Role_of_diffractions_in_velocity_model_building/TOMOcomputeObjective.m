function Funct = TOMOcomputeObjective (dv)

%function [Funct, Grad] = TOMOcomputeObjective (dv)
 
INV_settings; 
 

dv = reshape(dv, length(target.xx), length(target.zz));
bvelmod(target.xx, target.zz) =  bvelmod(target.xx, target.zz) + dv; 
  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
callTOMOcomputeObjective=MLD('callTOMOcomputeObjective.mat'); 

Grad = zeros(length(target.xx), length(target.zz)); 
Funct = 0; 

bfb = zeros(length(acq.sx), length(acq.rx)); 
mfb = zeros(length(acq.sx), length(acq.rx)); 

for shotnum = 1:1:length(acq.sx);
    
    S = [acq.sx(shotnum), acq.sz(shotnum)];

    % baseline
    btti = FSM2D(Gold,S,bvelmod);
    for recenum=1:length(acq.rx)
        bfb(shotnum, recenum) = btti(acq.grx(recenum),acq.grz(recenum));
    end
     
    % monitor 
    mtti = FSM2D(Gold,S,mvelmod);
    for recenum=1:length(acq.rx)
        mfb(shotnum, recenum) = mtti(acq.grx(recenum),acq.grz(recenum));
    end
         
    % residual timeshifts        
    dt = mfb(shotnum,:) - bfb(shotnum,:);
    dcr = zeros(size(mvelmod)); 
    for recenum=1:length(acq.rx)
        dcr(acq.grx(recenum),acq.grz(recenum)) = dt(recenum);
    end
    
    % adjoint state
    grad = ASM2D(Gold,dcr,btti);
          
    % functional and gradient
    funct = 0.5*sum(dt.^2); 
    Grad = Grad + grad(target.xx, target.zz);
    Funct = Funct + funct; 
end

Grad = Grad.*weight;
     
%% Calculate gradient and functional
% dJ/dm - is result of correlation
% now we multiply it by dm/dv   (m=1/v^2)
vel = bvelmod(target.xx, target.zz); 
Grad = Grad./vel.^3; 
Const = (1e+10)/29;     


Grad = Const*Grad;
Funct = Const*Funct; 

% apply regularization
% funct 10^10 - 10^5
% regularization - 10^7

Funct = Funct + lambda*0.5*sum(dv(:).^2);
Grad = Grad + lambda*dv;
     
     
%% save 
fname_G = [outfolder tomobasename '_G_' num2str(callTOMOcomputeObjective) '.mat'];
fname_F = [outfolder tomobasename '_F_' num2str(callTOMOcomputeObjective) '.mat'];
fname_dV = [outfolder tomobasename '_dv_' num2str(callTOMOcomputeObjective) '.mat'];
save(fname_G, 'Grad');
save(fname_F, 'Funct');
save(fname_dV, 'dv');
     
callTOMOcomputeObjective = callTOMOcomputeObjective + 1; 
save('callTOMOcomputeObjective.mat', 'callTOMOcomputeObjective'); 
   
Grad = reshape(Grad,1,length(target.xx)*length(target.zz));
save('TOMOGrad.mat', 'Grad'); 
     
end
