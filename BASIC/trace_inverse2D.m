function ray = trace_inverse2D(Gold, S, alpha, t0, velmod)    

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

    xx = Gox:Gdx:Gmx; 
    yy = Goy:Gdy:Gmy; 
    zz = Goz:Gdz:Gmz; 
    tt = Got:Gdt:Gmt; 

    %% Find derivatives of velocity

    dvdx = zeros(size(velmod));
    dvdz = zeros(size(velmod)); 

    dvdx(1:end-1,:) = diff(velmod,1,1)/Gdx; 
    dvdx(end,:) = dvdx(end-1,:); 
    dvdz(:,1:end-1) = diff(velmod,1,2)/Gdz;
    dvdz(:,end) = dvdz(:,end-1);

    %% Find Hamiltonian

    y=zeros(5, Gnt+1);

    % find velocity in the origin of the ray

    v0 = velmod(x2grid(S(1),Gox,Gdx,Gnx), x2grid(S(2),Goz,Gdz,Gnz)); 

    % set initial values
    y(1, 1) = S(1);
    y(2, 1) = S(2);
    y(3, 1) = sin(alpha)/v0;            % y(3:4) is slowness components
    y(4, 1) = cos(alpha)/v0;            % only direction is important
    y(5 ,1) = t0;

    i=1;   t=t0;
    while  t>0 && y(1, i)>=Gox && y(1, i)<=Gmx && y(2, i)>=Goz && y(2, i)<=Gmz 
            k1 = funct(y(1:4, i),           velmod, dvdx, dvdz, Gold);
            k2 = funct(y(1:4, i)+Gdt*k1/2,  velmod, dvdx, dvdz, Gold);
            k3 = funct(y(1:4, i)+Gdt*k2/2,  velmod, dvdx, dvdz, Gold);
            k4 = funct(y(1:4, i)+Gdt*k3,    velmod, dvdx, dvdz, Gold);
            y(1:4, i+1) = y(1:4, i) + Gdt*(k1+2*k2+2*k3+k4)/6;
            y(5, i+1) = y(5, i) - Gdt;
            t = t - Gdt;
            i = i+1;
    end;
    if abs(y(5,i-1)) < Gdt
        dy = (y(:,i) - y(:,i-1))/Gdt;
        y(:, i) = y(:,i-1) + y(5,i-1)*dy; 
        ray = y(1:5, 1:i);
    else
        ray = y(1:5, 1:i-1);
    end

end
   
function f=funct(y, velmod, dvdx, dvdz, Gold)
     
    f = zeros(4,1);
    Gox=Gold(1); Goy=Gold(2);  Goz=Gold(3);  Got=Gold(10);
    Gnx=Gold(4); Gny=Gold(5);  Gnz=Gold(6);  Gnt=Gold(11);
    Gdx=Gold(7); Gdy=Gold(8);  Gdz=Gold(9);  Gdt=Gold(12);


    indx = x2grid(y(1),Gox,Gdx,Gnx);
    indz = x2grid(y(2),Goz,Gdz,Gnz);
    v0 = velmod(indx,indz); 
    dvdx0 = dvdx(indx,indz); 
    dvdz0 = dvdz(indx,indz); 

    % see Cherveny, Seismic ray theory, Eq. 3.1.14
    % dx
    f(1) = v0^2*y(3);         % dx/dt = v^2*p_x 
    f(2) = v0^2*y(4);         % dz/dt = v^2*p_z

    % dp
    f(3) = -dvdx0/v0;         % dp/dx = -1/v*dv/dx
    f(4) = -dvdz0/v0;         % dp/dz = -1/v*dv/dz
 
end