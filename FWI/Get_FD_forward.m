function [ p, ux, uz ] = Get_FD_forward( G, acq, shotnum, fpeak, vp, rho, vs)
%GET_FD_OMAN Returns wavefield calculated for particular velocity model and
%acquisition geometry
%   This function is based on EL2D code
%   Abakumov Ivan
%   University of Hamburg
%   25.06.2015

    %% Check input

    switch nargin
        case 5
            fprintf(1,'Constant dencity acoustic modeling: rho=const, vs=0.\n');
            rho=1000*ones(G.nx,G.nz); 
            vs = zeros(G.nx,G.nz); 
        case 6
            fprintf(1,'Acoustic modeling: vs=0.\n');
            vs = zeros(G.nx,G.nz);
        case 7
            fprintf(1,'Elastic modeling: vs=0.\n');
        otherwise
            error('Incorrect number of input parameters!\n');
    end
    
    %% Setting number of threads
    numthreads = [];
    % set number of threads here, otherwise number will be chosen automatically):
    %numthreads = 4; 
    if isempty(numthreads),
        import java.lang.*;
        numthreads = Runtime.getRuntime.availableProcessors;
    end

    setenv('OMP_NUM_THREADS',num2str(numthreads));
    fprintf(1,'el2d will be executed with %s threads.\n',getenv('OMP_NUM_THREADS'));

    %% Set up model and geometry
    x    = (G.x0:G.dx:G.mx);                       % x dimensions for model (in meters)
    z    = (G.z0:G.dz:G.mz);                       % z dimensions for model (in meters)
    xs   = acq.sx(shotnum);                        % x coordinates of shots (in meters)
    zs   = acq.sz(shotnum);                        % z coordinates of shots (in meters)
    xr   = acq.rx;                                 % x coordinates of receivers (in meters)
    zr   = acq.rz;                                 % z coordinates of receivers (in meters)
    xl   = x;                                   % x dimensions of local zone for snapshot output (in meters)
    zl   = z;                                   % z dimensions of local zone for snapshot output (in meters)
    nsht = length(xs);                          % number of shots
    nrec = length(xr);                          % number of receivers
    nx   = length(x);                           % x grid size
    nz   = length(z);                           % z grid size
    dx   = abs(x(2)-x(1));                      % x sampling
    dz   = abs(z(2)-z(1));                      % z sampling

    %% Check size of velocity model and dencity
    if (size(rho,1)~=nx || size(rho,2)~=nz),
        error('Incorrect size of rho!');
    end
    
    if (size(vp,1)~=nx || size(vp,2)~=nz),
        error('Incorrect size of vp!');
    end
    
    if (size(vs,1)~=nx || size(vs,2)~=nz),
        error('Incorrect size of vs!');
    end
    
    %% Absorbing parameters
    damp = [30 0.075]; %size of absrobing layer in grid points and coefficient
    % damp = 30; %just size of absorbing layer in grid points (default coefficient 0.075 will be used)
    % damp = []; %no absorbing boundaries

    %% Source
    stype = 0;               % source type (0 - explosion; 1 - shear source; 2 - horizontal force; 3 - vertical force)% stype = int32(0);
    tmin  = G.t0;            
    tmax  = G.mt;
    dt    = G.dt;             % time sampling
    t     = tmin:dt:tmax;    % time array
    nt    = length(t);       % number of time samples
    wt    = ricker(fpeak,t); % construct ricker wavelet

    %% Absorbing layer parameters
    % Indices of local axes inside full axes
    ilx = arrayfun(@(a)find(x==a,1),xl);
    ilz = arrayfun(@(a)find(z==a,1),zl);

    % Indices of shot positions inside full axes
    isx = arrayfun(@(a)find(x==a,1),xs);
    isz = arrayfun(@(a)find(z==a,1),zs);

    % Indices of receiver positions inside full axes
    irx = arrayfun(@(a)find(x==a,1),xr);
    irz = arrayfun(@(a)find(z==a,1),zr);

    % Check shot and receiver positions
    if (length(isx)~=nsht || length(isz)~=nsht),
        error('Incorrect shot positions!');
    end
    if (length(irx)~=nrec || length(irz)~=nrec),
        error('Incorrect receiver positions!');
    end

    %% Prepare input data for MEX
    irx   = int32(irx - 1); %C-style x indices of receivers
    irz   = int32(irz - 1); %C-style z indices of receivers
    isx   = int32(isx - 1); %C-style x indices of shots
    isz   = int32(isz - 1); %C-style z indices of shots
    ilx   = int32(ilx - 1); %C-style x indices of local zone
    ilz   = int32(ilz - 1); %C-style z indices of local zone
    stype = int32(stype);
    nt    = int32(nt);
    wt    = single(wt);
    fpeak = single(fpeak);
    dt    = single(dt);
    dx    = single(dx);
    dz    = single(dz);
    damp  = single(damp);
    rho   = single(rho);
    vp    = single(vp);
    vs    = single(vs);

    %% Check parameters
    if min(min(abs(vs)))==0
        vmin = min(min(vp));
    else
        vmin = min([min(min(vs)) min(min(vp))]);
    end
    vmax = max([max(max(vs)) max(max(vp))]);
    
    gamma = sqrt(12);
    crit1 = (vmin/(2.*fpeak))/max(dx,dz);                                   % wavelength for double peak frequency, minimum velocity
    crit2 = gamma/(dt*pi*vmax*sqrt(1./(dx^2)+1./(dz^2)));
    err = 0;
    if crit1 >= 2,
        fprintf(1,'Minimum number of grid points per shortest wavelength: %g >= 2.\n',crit1);
    else
        fprintf(1,'Minimum number of grid points per shortest wavelength: %g < 2!\n',crit1);
        err = 1;
    end
    if crit2 > 1,
        fprintf(1,'Stability parameter: %g > 1.\n',crit2);
    else
        fprintf(1,'Stability parameter: %g <= 1!\n',crit2);
        err = 1;
    end
    if err==1,
        error('Problems with stability!');
    end

    %% Run
    tic;
    fprintf(1,'Computing wavefields...');
    % The code can output spatial derivatives of displacement as well:
    % [p, ux, uz, uxx, uzz, uxz, uzx] = ...
    % where uxz stands for z derivative of x component of displacement etc
    [p, ux, uz] = el2d(0,rho,vp,vs,dx,dz,damp,ilx,ilz,stype,isx,isz,fpeak,dt,nt,wt);
    fprintf(1,'Done. ');
    toc;
end

