
thisfolder = pwd;

%% Load files

acq = MLD([thisfolder '/data/OBN_DIFF_acquisition.mat']);

% 12.5 x 12.5 m grid
G = MLD([thisfolder '/data/OBN_DIFF_G_file.mat']);
bvelmod = MLD([thisfolder '/data/OBN_DIFF_bvelmod.mat']);
mvelmod = MLD([thisfolder '/data/OBN_DIFF_mvelmod.mat']);

tomobasename = 'TOMO_First';
avaxbasename = 'AVAX_First';
fwibasename = 'FWI_First';

outfolder = [thisfolder '/result/TOMO2/'];

% make old G
Gold = oldGrid(G);

%% Frequency (for FWI)
% 1 m 
%fpeak=100; 
%fpeak=250;     
%fpeak=500;

% 12.5 m 
%fpeak = 8; 
fpeak = 20; 
%fpeak = 40; 

% 25 m 
%fpeak = 4; 
%fpeak = 10; 
%fpeak = 20; 

%% Regularization factor
lambda = 0.0000001; 

%% Number of threads
nthreads = 13; 

%% Define target zone

target.xx = 151:351;
target.zz = 121:181;

%% Make weight function
%weight = MLD([thisfolder '/data/OBN_DIFF_weight.mat']);
% 
% weight = ones(G.nx, G.nz); 
% wind = 80; % in meter 
% pwind = 28;  
% for s=1:length(acq.sx)
%     for i=x2grid(acq.sx(s)-wind,G.x0,G.dx,G.nx):x2grid(acq.sx(s)+wind,G.x0,G.dx,G.nx)
%         for j=x2grid(acq.sz(s)-wind,G.z0,G.dz,G.nz):x2grid(acq.sz(s)+wind,G.z0,G.dz,G.nz)
%             f = 1-exp(-sqrt((acq.gsx(s)-i)^2 + (acq.gsz(s)-j)^2)/pwind); 
%             weight(i,j) = min(weight(i,j),f); 
%         end
%     end
% end
% for r=1:length(acq.rx)
%     for i=x2grid(acq.rx(r)-wind,G.x0,G.dx,G.nx):x2grid(acq.rx(r)+wind,G.x0,G.dx,G.nx)
%         for j=x2grid(acq.rz(r)-wind,G.z0,G.dz,G.nz):x2grid(acq.rz(r)+wind,G.z0,G.dz,G.nz)
%             f = 1-exp(-sqrt((acq.grx(r)-i)^2 + (acq.grz(r)-j)^2)/pwind); 
%             weight(i,j) = min(weight(i,j),f); 
%         end
%     end
% end
% weight = weight(target.xx,target.zz);  
% weight = SME(weight,1); 
% imagesc(weight')