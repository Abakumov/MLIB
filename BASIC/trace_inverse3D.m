function ray = trace_inverse3D(Gold, S, alpha, beta, t0, velmod)    

    % Ivan Abakumov

    % Ray tracer in 2D isotropic medium
    % Input:    G, S, velmod, alpha    
    % Output:   ray = matrix[5xGnt]
    %           arrays: [x, z, px, pz, t]
    % Method Runge-Kutt, 4th order

    %% tmp
    
%     clear all; close all; clc;
%     mlibfolder = '/home/zmaw/u250128/Desktop/MLIB';
%     path(path, mlibfolder);
%     addmypath;
%     current_folder = pwd;
%     
%     G=GridClass;
%     
%     %   [m]             [m]          [m]             [s]
%     G.x0 = -2475;   G.y0 = 0;    G.z0 = -5000;   G.t0 = 0.00;            % initial point
%     G.nx =  1000;   G.ny = 1;    G.nz =  800;    G.nt = 1000;             % grid size
%     G.dx =  25;     G.dy = 25;   G.dz =  25;     G.dt = 0.007;         % grid step (meter)
%     
%     
%     G.gridInfo;
%     G.setGrid;
%     Gold = oldGrid(G);
%     
%     velmod = zeros(G.nx, G.nz);
%     
%     velmod(:, 1:240)=1500;
%     for i=241:390
%         velmod(:, i)=velmod(:,i-1)+8.0;
%     end
%     for i=391:G.nz
%         velmod(:,i)=velmod(:,390);
%     end
%     
%     S = [5000 0];
%     alpha = pi/6;

    %% Read input

    Gox=Gold(1); Goy=Gold(2);  Goz=Gold(3);  Got=Gold(10);
    Gnx=Gold(4); Gny=Gold(5);  Gnz=Gold(6);  Gnt=Gold(11);
    Gdx=Gold(7); Gdy=Gold(8);  Gdz=Gold(9);  Gdt=Gold(12);

    Gmx = Gox + (Gnx-1)*Gdx;  
    Gmy = Goy + (Gny-1)*Gdy; 
    Gmz = Goz + (Gnz-1)*Gdz; 
    Gmt = Got + (Gnt-1)*Gdt; 

 

    %% Find derivatives of velocity

    dvdx = zeros(size(velmod));
    dvdy = zeros(size(velmod));
    dvdz = zeros(size(velmod)); 

    dvdx(1:end-1,:,:) = diff(velmod,1,1)/Gdx; 
    dvdx(end,:,:) = dvdx(end-1,:,:); 
    dvdy(:,1:end-1,:) = diff(velmod,1,2)/Gdy; 
    dvdy(:,end,:) = dvdy(:,end-1,:); 
    dvdz(:,:,1:end-1) = diff(velmod,1,3)/Gdz;
    dvdz(:,:,end) = dvdz(:,:,end-1);

    %% Find Hamiltonian

    y=zeros(7, Gnt+1);

    % find velocity in the origin of the ray

    v0 = velmod(x2grid(S(1),Gox,Gdx,Gnx), x2grid(S(2),Goy,Gdy,Gny),x2grid(S(3),Goz,Gdz,Gnz)); 

    % set initial values
    y(1, 1) = S(1);
    y(2, 1) = S(2);
    y(3, 1) = S(3); 
    y(4, 1) = sin(alpha)*cos(beta)/v0;            % y(4:6) is slowness components
    y(5, 1) = sin(alpha)*sin(beta)/v0;            % only direction is important
    y(6, 1) = cos(alpha)/v0; 
    y(7 ,1) = t0;

    i=1;   t=t0;
    while  t>0 && y(1, i)>=Gox && y(1, i)<=Gmx && y(2, i)>=Goy && y(2, i)<=Gmy && y(3, i)>=Goz && y(3, i)<=Gmz 
            k1 = funct(y(1:6, i),           velmod, dvdx, dvdy, dvdz, Gold);
            k2 = funct(y(1:6, i)+Gdt*k1/2,  velmod, dvdx, dvdy, dvdz, Gold);
            k3 = funct(y(1:6, i)+Gdt*k2/2,  velmod, dvdx, dvdy, dvdz, Gold);
            k4 = funct(y(1:6, i)+Gdt*k3,    velmod, dvdx, dvdy, dvdz, Gold);
            y(1:6, i+1) = y(1:6, i) + Gdt*(k1+2*k2+2*k3+k4)/6;
            y(7, i+1) = y(7, i) - Gdt;
            t = t - Gdt;
            i = i+1;
    end;
    if abs(y(7,i-1)) < Gdt
        dy = (y(:,i) - y(:,i-1))/Gdt;
        y(:, i) = y(:,i-1) + y(7,i-1)*dy; 
        ray = y(1:7, 1:i);
    else
        ray = y(1:7, 1:i-1);
    end

end
   
function f=funct(y, velmod, dvdx, dvdy, dvdz, Gold)
     
    f = zeros(6,1);
    Gox=Gold(1); Goy=Gold(2);  Goz=Gold(3);
    Gnx=Gold(4); Gny=Gold(5);  Gnz=Gold(6);
    Gdx=Gold(7); Gdy=Gold(8);  Gdz=Gold(9);


    indx = x2grid(y(1),Gox,Gdx,Gnx);
    indy = x2grid(y(2),Goy,Gdy,Gny); 
    indz = x2grid(y(3),Goz,Gdz,Gnz);
    v0  = velmod(indx,indy,indz); 
    dvdx0 = dvdx(indx,indy,indz);
    dvdy0 = dvdy(indx,indy,indz);
    dvdz0 = dvdz(indx,indy,indz); 

    % see Cherveny, Seismic ray theory, Eq. 3.1.14
    % dx
    f(1) = v0^2*y(4);         % dx/dt = v^2*p_x 
    f(2) = v0^2*y(5);         % dy/dt = v^2*p_y
    f(3) = v0^2*y(6);         % dz/dt = v^2*p_z

    % dp
    f(4) = -dvdx0/v0;         % dp/dx = -1/v*dv/dx
    f(5) = -dvdy0/v0;         % dp/dy = -1/v*dv/dy
    f(6) = -dvdz0/v0;         % dp/dz = -1/v*dv/dz
end