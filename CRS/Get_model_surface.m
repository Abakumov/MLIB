function [zref, ind] = Get_model_surface(XX,YY,model)
% This function returns depth of the reflector surface
% zref = f(XX, YY)
% and indicator that surface exists in the point
% Abakumov Ivan
% 27th August 2016

Get_model_parameters;

eps = 1e-20; 

% 1. Flat reflector
if model == 11 || model == 12 || model == 13 || model == 14
    zref = zc*ones(size(XX)); 
    ind = (zref>0); 
%    ze = real(sqrt(Rc^2 - (XX-xc).^2 - (YY-yc).^2));
%    zref = zc - ze;
%    ind = 1 - ze./(ze+eps);  
%    ind = (ind==0); 

% 2. Plane dipping reflector
elseif model == 21 || model == 22 || model == 23 || model == 24
    zref = zc*ones(size(XX)) + (XX-xc)*sin(dip)*cos(az) + (YY-yc)*sin(dip)*sin(az);
    ind = (zref>0); 
%    ze = real(sqrt(Rc^2 - (XX-xc).^2 - (YY-yc).^2));
%    zref = zc - ze;
%    ind = 1 - ze./(ze+eps);  
%    ind = (ind==0); 

% 3. Point diffractor: R = 10 m  and 4. Sphere: R = 1 km
elseif model == 31 || model == 32 || model == 33 || model == 34 || model == 41 || model == 42 || model == 43 || model == 44
    ze = real(sqrt(Rc^2 - (XX-xc).^2 - (YY-yc).^2));
    zref = zc - ze;
    ind = 1 - ze./(ze+eps);  
    ind = (ind==0); 

% 5. Ellipsoid
elseif model == 51 || model == 52 || model == 53 || model == 54 || model == 0 
    ze = de*real(sqrt(1 - (XX-xc).^2/ae^2 - (YY-yc).^2/be^2 - 2*(XX-xc).*(YY-yc)/ce^2));
    zref = zc - ze;
    ind = 1 - ze./(ze+eps); 
    ind = (ind==0); 
    
% 6. Complex analytical surface
elseif model == 61 || model == 62 || model == 63 || model == 64
    zref = f(XX, YY); 
    ind = (zref>0); 
    
% 7. Paraboloidal reflector in the simplified model
elseif model > 100 && model < 200
    target_model = model - floor(model/100)*100; 
    CRS_param = MLD(['/home/zmaw/u250128/Desktop/MLIB/CRS/models/model_' num2str(target_model) '_CRS_param.mat']); 
    smodel = Get_simplified_model(CRS_param);
    ref_par = Get_smodel_reflector_parabolic(smodel, XX, YY); 
    zref = ref_par.zz_src; 
    ind = (zref>0); 

% 8. ellipsoidal reflector in the simplified model
elseif model > 200
    target_model = model - floor(model/100)*100; 
    CRS_param = MLD(['/home/zmaw/u250128/Desktop/MLIB/CRS/models/model_' num2str(target_model) '_CRS_param.mat']); 
    smodel = Get_simplified_model(CRS_param);
    ref_elc = Get_smodel_reflector_ellipsoidal_cart(smodel, XX, YY);
    zref = ref_elc.zz_src; 
    ind = ref_elc.ind; 
end