%% Table 200 - table with traveltime errors
% 
% *Author*: Abakumov Ivan
% *Publication date*: 31st August 2016

%% Define working folder, add links to Library and SeisLab 

clear; close all; clc;
mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
path(path, mlibfolder);
addmypath;

%% Introduction

% Note:
% Accuracy of traveltime 10e-12
% Accuracy of attributes 10e-10

for model = [11 12 13 21 22 23 31 32 33 41 42 43 51 52 53 61 62 63 ]
%for model = 21:21
    acquisition = 3;

    % Get model parameters
    Get_model_parameters; 
    
    % Get acquisition geometry
    Get_model_acquisition_geometry; 
    
    % Load stacking parameters

    CRS_param = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_CRS_param.mat']); 
      
    X0 = CRS_param.x0; 
    t0 = CRS_param.t0; 
    v0 = CRS_param.v0; 
    w3 = CRS_param.w; 
    M3 = CRS_param.M; 
    N3 = CRS_param.N; 
   
    % take (x,y) components 
     % (reality)
     w = w3(1:2, 1); 
     N = N3(1:2, 1:2);
     M = M3(1:2, 1:2);  
    
    % Download exact traveltimes
    tti_ex = MLD([mlibfolder '/CRS/models/model_' num2str(model) '_traveltimes_for_acq_' num2str(acquisition) '.mat']);
    
    %% Traveltime approximations
    
     HH = (Xg(1:2, :) - Xs(1:2,:))/2;
     MM = (Xg(1:2, :) + Xs(1:2,:))/2;
     MM(1,:) = MM(1,:) - X0(1); 
     MM(2,:) = MM(2,:) - X0(2);
     
     fold = size(HH,2); 
     
     %
     tti_crs  = real(Get_traveltime_3D_CRS (MM, HH, t0, w, M, N));
     tti_dsr  = real(Get_traveltime_3D_DSR (MM, HH, t0, w, M, N));
     tti_ncrs = real(Get_traveltime_3D_nCRS(MM, HH, t0, w, M, N));
     tti_icrs_par_LIA = real(Get_traveltime_3D_iCRS_par_LIA(MM, HH, t0, w, M, N, 4));
     tti_icrs_el_LIA = real(Get_traveltime_3D_iCRS_el_LIA(MM, HH, t0, v0, w, M, N, 4));
     tti_icrs_el_TIA = real(Get_traveltime_3D_iCRS_el_TIA(MM, HH, t0, v0, w, M, N, 4));

     err_rms_crs  = sqrt(sum(((tti_crs  - tti_ex)./tti_ex).^2)/fold)*100;
     err_rms_dsr  = sqrt(sum(((tti_dsr  - tti_ex)./tti_ex).^2)/fold)*100;
     err_rms_ncrs = sqrt(sum(((tti_ncrs - tti_ex)./tti_ex).^2)/fold)*100;
     err_rms_icrs_par_LIA = sqrt(sum(((tti_icrs_par_LIA(4,:) - tti_ex)./tti_ex).^2)/fold)*100;
     err_rms_icrs_el_LIA  = sqrt(sum(((tti_icrs_el_LIA(4,:)  - tti_ex)./tti_ex).^2)/fold)*100;
     err_rms_icrs_el_TIA  = sqrt(sum(((tti_icrs_el_TIA(4,:)  - tti_ex)./tti_ex).^2)/fold)*100;
    
     disp(['Model #' num2str(model)]);
     disp('Err: CRS, DSR, nCRS, iCRS_par_LIA, i-CRS_el_LIA, i-CRS_el_TIA'); 
     err = [err_rms_crs, err_rms_dsr, err_rms_ncrs, err_rms_icrs_par_LIA, err_rms_icrs_el_LIA, err_rms_icrs_el_TIA   ]

end

