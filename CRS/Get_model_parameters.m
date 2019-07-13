%% This script defines parameters of the model
% Abakumov Ivan
% 27th August 2016

%% Parameters of the reflector

% 1. Flat reflector
if     model == 11 || model == 12 || model == 13 || model == 14
    zc = 1000; 
    %xc = 2000;  yc = 2000;  zc = 6401000;
    %Rc = 6400000;  

% 2. Plane dipping reflector
elseif model == 21 || model == 22 || model == 23 || model == 24
     xc = 2000;  yc = 2000;  zc = 1000;
     dip = -pi/12; az = pi/6;
%    xc = 2500;  yc = 1500;  zc = 10500;
%    Rc = 10000;  

% 3. Point diffractor: R = 10 m
elseif model == 31 || model == 32 || model == 33 || model == 34
    xc = 1520;  yc = 2720;  zc = 1010;
    Rc = 10;  

% 4. Sphere: R = 1 km
elseif model == 41 || model == 42 || model == 43 || model == 44
    xc = 1750;  yc = 2400;  zc = 2000;
    Rc = 1000;

% 5. Ellipsoid
elseif model == 51 || model == 52 || model == 53 || model == 54 
    xc = 1750;  yc = 2400;  zc = 2000;
    ae = 1000;   be = 1500;  ce = 2000;  de = 1000;

% 6. Complex analytical surface
elseif model == 61 || model == 62 || model == 63 || model == 64
    f = MLD('/home/ivan/Desktop/MLIB/CRS/models/model_6x_reflector.mat');

% 0. Model to check 3D i-CRS formulas
elseif model == 0   
    xc = 2000;  yc = 2000;  zc = 2000;
    Ae = 900; Be = 1600; de = 1200;
    phi = 0; 
    ae = 1/sqrt((cos(phi)/Ae)^2 + (sin(phi)/Be)^2); 
    be = 1/sqrt((cos(phi)/Be)^2 + (sin(phi)/Ae)^2); 
    ce = 1/sqrt((1/Ae^2 - 1/Be^2)*sin(phi)*cos(phi)); 
    
    Rc = 1000; 
end
    
%% Parameters of the velocity model 

% 1. V = const
if     model == 11 || model == 21 || model == 31 || model == 41 || model == 51 || model == 61 || model == 0
    v0 = 1500;

% 2. V = v0 + k*z
elseif model == 12 || model == 22 || model == 32 || model == 42 || model == 52 || model == 62
    v0 = 1500;   
    k = 0.5;
% 3. With water layer
% V = v0              z < z0
% V = v0 + k(z-z0)    z > z0
elseif model == 13 || model == 23 || model == 33 || model == 43 || model == 53 || model == 63
    v0 = 1500;   
    k = 0.5;
    z0 = 250; 
    
% 4. Elliptical anisotropy  
elseif model == 14 || model == 24 || model == 34 || model == 44 || model == 54 || model == 64
    A11 = 1500*1500; 
    A22 = 1700*1700; 
    A33 = 1800*1800; 
    v0 = sqrt(A33);
    
elseif model > 100
    target_model = model - floor(model/100)*100; 
    CRS_param = MLD(['/home/ivan/Desktop/MLIB/CRS/models/model_' num2str(target_model) '_CRS_param.mat']); 
    smodel = Get_simplified_model(CRS_param);
    A11 = smodel.A11; 
    A22 = smodel.A22; 
    A33 = smodel.A33; 
    v0 = sqrt(A33); 
    clear target_model smodel R CRS_param
end

%% G-file

if model < 100
     G = MLD('/home/ivan/Desktop/MLIB/CRS/models/model_G_file.mat');
    
elseif model > 100
     G = MLD('/home/ivan/Desktop/MLIB/CRS/models/model_G_file.mat');
     G.x0 = G.x0 - 2000; 
     G.y0 = G.y0 - 2000; 
     G.setGrid;
end 

