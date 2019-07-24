function Funct = FWIcomputeObjective (dv)

INV_settings; 

%dv  = zeros(1,length(target.xx)*length(target.zz));
dv = reshape(dv, length(target.xx), length(target.zz));
bvelmod(target.xx, target.zz) =  bvelmod(target.xx, target.zz) + dv; 

Grad = zeros(length(target.xx), length(target.zz)); 
Funct = 0; 

%% Calculate traveltimes of direct waves

tdir = zeros(length(acq.sx), length(acq.rx)); 

for s=1:length(acq.sx)
    S = [acq.sx(s), acq.sz(s)];
    tti = FSM2D(Gold,S,bvelmod);
    
    for r=1:length(acq.rx)
        tdir(s,r) = tti(acq.grx(r),acq.grz(r));
    end
end

tdir = round((tdir)/G.dt + 1) + 25; 
%% Make sin filter

filter = ones(1,101);
filter(1:26) = sin(linspace(0, pi/2, 26)).^2;
filter(76:101) = sin(linspace(pi/2, 0, 26)).^2;

%%


for ishot=1:nthreads:length(acq.sx);

    shotnum = ishot:min(length(acq.sx), ishot+nthreads-1); 
 
    % baseline
    [f09,~,~] = Get_FD_forward(G, acq, shotnum, fpeak, bvelmod);
    s09 = zeros(length(acq.rx), G.nt, length(shotnum));
    for s=1:length(shotnum); 
        for i=1:length(acq.rx)
            minind = max(1,tdir(shotnum(s),i));
            maxind = min(G.nt, tdir(shotnum(s),i) + 100);
            ind = minind:maxind;
            filt = zeros(G.nt,1);
            filt(ind) = filter(1:length(ind));
            s09(i, :, s) = squeeze(f09(acq.grx(i),acq.grz(i),:,s)).*filt;
        end
    end
    U = f09(target.xx,target.zz,:,:);  
    clear f09; 
 
    % monitor 
    [f10,~,~] = Get_FD_forward(G, acq, shotnum, fpeak, mvelmod);
    s10 = zeros(length(acq.rx), G.nt,length(shotnum));
    for s=1:length(shotnum); 
        for i=1:length(acq.rx)
            minind = max(1,tdir(shotnum(s),i));
            maxind = min(G.nt, tdir(shotnum(s),i) + 100);
            ind = minind:maxind;
            filt = zeros(G.nt,1);
            filt(ind) = filter(1:length(ind));
            s10(i,:,s) = squeeze(f10(acq.grx(i),acq.grz(i),:,s)).*filt;
        end
    end
    clear f10; 
         
    % residual seismograms        
    ds = s10-s09; 
    ds = ds(:, end:-1:1,:); 
    
    
          
    % adjoint state
    [q,~,~] = Get_FD_backward(G, acq, shotnum, fpeak, ds, bvelmod);
    Q = q(target.xx,target.zz,:,:);
    clear q
         
    % functional and gradient
    funct = 0.5*sum(sum(sum(ds.^2))).*G.dt; 
    grad=make_grad_LS(G,U,Q);
    clear U Q

    Grad = Grad + grad;
    Funct = Funct + funct; 
end

Grad = Grad.*weight;
     
%% Calculate gradient and functional
% dJ/dm - is result of correlation
% now we multiply it by dm/dv   (m=1/v^2)
vel = bvelmod(target.xx, target.zz); 
Grad = -2*Grad./vel.^3; 
Const = (1e+18)/29;      
Grad = Const*Grad;
Funct = Const*Funct; 
     
     
%% save 
callFWIcomputeObjective=MLD('callFWIcomputeObjective.mat'); 
fname_G = [outfolder fwibasename '_G_' num2str(callFWIcomputeObjective) '.mat'];
fname_F = [outfolder fwibasename '_F_' num2str(callFWIcomputeObjective) '.mat'];
fname_dV = [outfolder fwibasename '_dv_' num2str(callFWIcomputeObjective) '.mat'];
fname_fr = [outfolder fwibasename '_fr_' num2str(callFWIcomputeObjective) '.mat'];
save(fname_G, 'Grad');
save(fname_F, 'Funct');
save(fname_dV, 'dv');
save(fname_fr, 'fpeak');
     
callFWIcomputeObjective = callFWIcomputeObjective + 1; 
save('callFWIcomputeObjective.mat', 'callFWIcomputeObjective'); 
   
Grad = reshape(Grad,1,length(target.xx)*length(target.zz));
save('Grad.mat', 'Grad'); 
     
end
  
  
  
  
  
  
  
  
  