% This script returns acquisition geometry
% Abakumov Ivan
% 27th August 2016

% 2D line
if acquisition == 1 
    
    X0 = [2000; 2000; 0]; 
    azimuth = pi/4; 
    
    xs = (-2000:25:2000)*cos(azimuth) + X0(1);
    ys = (-2000:25:2000)*sin(azimuth) + X0(2);
    zs = zeros(size(xs));

    xg = (-2000:25:2000)*cos(azimuth) + X0(1);
    yg = (-2000:25:2000)*sin(azimuth) + X0(2);
    zg =  zeros(size(xg));
    
    Xs = zeros(3,length(xs)*length(xg)); 
    Xg = zeros(3,length(xs)*length(xg));
    
    for i=1:length(xs)
        for j=1:length(xg)
            k = j + (i-1)*length(xg);
            Xs(1,k) = xs(i);
            Xs(2,k) = ys(i);
            Xs(3,k) = zs(i);
            
            Xg(1,k) = xg(j);
            Xg(2,k) = yg(j);
            Xg(3,k) = zg(j);
         end
    end
    clear xs xg ys yg zs zg azimuth
    clear i j k
    
elseif acquisition == 2 
    
    X0 = [2000; 2000; 0]; 
    azimuth = pi/4; 
    
    hx = (0:25:2000)*cos(azimuth);
    hy = (0:25:2000)*sin(azimuth);
    
    mx = (-500:25:500)*cos(azimuth);
    my = (-500:25:500)*sin(azimuth);
    
    Xs = zeros(3,length(mx)*length(hx)); 
    Xg = zeros(3,length(mx)*length(hx));
    
    for i=1:length(mx)
        for j=1:length(hx)
            k = j + (i-1)*length(hx);
            Xs(1,k) = mx(i) - hx(j) + X0(1);
            Xs(2,k) = my(i) - hy(j) + X0(2);
            Xs(3,k) = 0;
            
            Xg(1,k) = mx(i) + hx(j) + X0(1);
            Xg(2,k) = my(i) + hy(j) + X0(2);
            Xg(3,k) = 0;
         end
    end
    clear xs xg ys yg zs zg azimuth
    clear i j k 
    
    
    
    
elseif acquisition == 3
    X0 = [2000; 2000; 0]; 
%     fold = 100;
%     hmax = 1000; 
%     mmax = 500; 
%     hmod = rand(1,fold)*hmax;
%     mmod = rand(1,fold)*mmax; 
%     haz  = rand(1,fold)*2*pi; 
%     maz  = rand(1,fold)*2*pi; 
%     hx = hmod.*cos(haz);
%     hy = hmod.*sin(haz); 
%     mx = mmod.*cos(maz); 
%     my = mmod.*sin(maz); 
%     Xs = zeros(3,fold);
%     Xg = zeros(3,fold);
% 
%     Xs(1,:) = mx - hx + X0(1);
%     Xs(2,:) = my - hy + X0(2);
%     Xg(1,:) = mx + hx + X0(1);
%     Xg(2,:) = my + hy + X0(2);
    
    acq = MLD('/home/zmaw/u250128/Desktop/MLIB/CRS/models/acquisition_3.mat');
    Xs = acq.Xs; 
    Xg = acq.Xg;
    clear hx hy mx my mmax mmod hmod hmax haz maz fold
    
    
elseif acquisition == 4
    % CMP geometry
    azimuth = 0:pi/18:2*pi;
    offset = 0:100:2000; 
    
    X0 = [2000; 2000; 0]; 
    Xs = zeros(3,length(azimuth)*length(offset)); 
    Xg = zeros(3,length(azimuth)*length(offset));
    
    for i = 1:length(offset);
        for j = 1:length(azimuth); 
            k = j + (i-1)*length(azimuth);
            Xs(1,k) = offset(i)*cos(azimuth(j)) + X0(1);
            Xs(2,k) = offset(i)*sin(azimuth(j)) + X0(2);
            Xs(3,k) = 0;
            
            Xg(1,k) = - offset(i)*cos(azimuth(j)) + X0(1);
            Xg(2,k) = - offset(i)*sin(azimuth(j)) + X0(2);
            Xg(3,k) = 0;
        end
    end

   
    clear i j k
    
    elseif acquisition == 5
    % ZO geometry
    azimuth = 0:pi/18:2*pi;
    offset = 0:50:1000; 
    
    X0 = [2000; 2000; 0]; 
    Xs = zeros(3,length(azimuth)*length(offset)); 
    Xg = zeros(3,length(azimuth)*length(offset));
    
    for i = 1:length(offset);     % size(MM) = 2 x fold
        for j = 1:length(azimuth); 
            k = j + (i-1)*length(azimuth);
            Xs(1,k) = offset(i)*cos(azimuth(j)) + X0(1);
            Xs(2,k) = offset(i)*sin(azimuth(j)) + X0(2);
            Xs(3,k) = 0;
            
            Xg(1,k) = offset(i)*cos(azimuth(j)) + X0(1);
            Xg(2,k) = offset(i)*sin(azimuth(j)) + X0(2);
            Xg(3,k) = 0;
        end
    end
    
    
    clear i j k
    
end

% additional rotation to L' if model > 100

if model > 100
    target_model = model - floor(model/100)*100; 
    CRS_param = MLD(['/home/zmaw/u250128/Desktop/MLIB/CRS/models/model_' num2str(target_model) '_CRS_param.mat']); 
    smodel = Get_simplified_model(CRS_param);
    R = smodel.R; 
    X0 = smodel.x0; 
    for i = 1:size(Xs,2); 
        Xs(:,i) = R'*(Xs(:,i) - X0);
        Xg(:,i) = R'*(Xg(:,i) - X0); 
    end
    X0 = [0;0;0]; 
end
clear target_model smodel i R CRS_param






