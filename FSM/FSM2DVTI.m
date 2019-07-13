%% Main function
function out = FSM2DVTI(G, S, alpha, betta, gamma, delta, epsilon)

% Find minimum values
alpmin = min(min(alpha));
betmin = min(min(betta));
delmin = min(min(delta));
gammin = min(min(gamma));
epsmin = min(min(epsilon));

% Definition of elastic parameters
[Gnx, Gnz] = size(alpha);

A33 = ones(Gnx+2, Gnz+2)*alpmin^2;
A44 = ones(Gnx+2, Gnz+2)*betmin^2;
A11 = ones(Gnx+2, Gnz+2)*alpmin^2*(1+2*epsmin);
A66 = ones(Gnx+2, Gnz+2)*betmin^2*(1+2*gammin);
A13 = ones(Gnx+2, Gnz+2)*(alpmin^2 - betmin^2)*sqrt(1 + 2*delmin/(1-(betmin^2/alpmin^2))) - betmin^2;

for i=1:Gnx
    for j=1:Gnz
        alp2 = alpha(i,j)^2; 
        bet2 = betta(i,j)^2;
        del = delta(i,j);
        gam = gamma(i,j);
        eps = epsilon(i,j);
    
        A33(i+1, j+1) = alp2;
        A44(i+1, j+1) = bet2;
        A11(i+1, j+1) = alp2*(1+2*eps);      
        A66(i+1, j+1) = bet2*(1+2*gam);      
        A13(i+1, j+1) = (alp2 - bet2)*sqrt(1 + 2*del/(1-(bet2/alp2))) - bet2;
    end
end

% Initialize times in source zone manually
tti = initialize_VTI_qP(G,S,A11,A33,A13,A44);

% Update times with FSM algorithm
out = fsmalgorithm_VTI_qP(G,A11,A33,A13,A44,tti);


%% Initialization

function tti = initialize_VTI_qP(G,S,A11,A33,A13,A44)

tmax = 1.e20;
tti = ones(size(A11))*tmax;

Gox=G(1); Goz=G(3);
Gnx=G(4); Gnz=G(6); 
Gdx=G(7); Gdz=G(9); 

Gmx = Gox + (Gnx-1)*Gdx; 
Gmz = Goz + (Gnz-1)*Gdz; 

	
if(S(1)<Gox || S(1)>Gmx || S(2)<Goz || S(2)>Gmz)
    disp('Shot coordinates out of range');
end
Gsx = x2grid(S(1),Gox,Gdx,Gnx);
Gsz = x2grid(S(2),Goz,Gdz,Gnz);
tti(Gsx+1, Gsz+1)=0;        % +1 due to general shift
  
for i=-5:1:5
    for j=-5:1:5
        indx = max(1, Gsx+1+i); 
        indx = min(Gnx+2, indx); 
        indz = max(1, Gsz+1+j); 
        indz = min(Gnz+2, indz); 
        a11 = A11(Gsx+1, Gsz+1); 
        a33 = A33(Gsx+1, Gsz+1); 
        a13 = A13(Gsx+1, Gsz+1); 
        a44 = A44(Gsx+1, Gsz+1); 
        x = max(Gox, Gox + (Gsx+i-1)*Gdx);
        x = min(Gmx, x);
        z = max(Goz, Goz + (Gsz+j-1)*Gdz); 
        z = min(Gmz, z); 
        dx = x - S(1); 
        dz = z - S(2); 
        dist = sqrt(dx.^2 + dz.^2); 
        if dist == 0 
            tti(indx, indz)=0;
        else
            ST = dx/dist;
            Vg = GetVg(a11,a33,a13,a44,ST);
            tti(indx, indz)=dist/Vg;    
        end
   end
end

%% General FSM solver

function tvltm = fsmalgorithm_VTI_qP(G,A11,A33,A13,A44,tti)

step = G(7);
tmax = 1.e20;
[sx, sz] = size(tti);
xmin=2; xmax=sx-1;
zmin=2; zmax=sz-1;

sq1 = sqrt(2)/2;
sq2 = sqrt(2);
tr(1,:) = [-1.0, -1.0,  0.0, -1.0,  sq2,  1.0,  sq1,  sq1,  0.0,  1.0];
tr(2,:) = [ 0.0, -1.0,  1.0, -1.0,  1.0,  sq2,  0.0,  1.0, -sq1,  sq1];
tr(3,:) = [ 1.0, -1.0,  1.0,  0.0,  sq2,  1.0, -sq1,  sq1, -1.0,  0.0];
tr(4,:) = [ 1.0,  0.0,  1.0,  1.0,  1.0,  sq2, -1.0,  0.0, -sq1, -sq1];
tr(5,:) = [ 1.0,  1.0,  0.0,  1.0,  sq2,  1.0, -sq1, -sq1,  0.0, -1.0];
tr(6,:) = [ 0.0,  1.0, -1.0,  1.0,  1.0,  sq2,  0.0, -1.0,  sq1, -sq1];
tr(7,:) = [-1.0,  1.0, -1.0,  0.0,  sq2,  1.0,  sq1, -sq1,  1.0,  0.0];
tr(8,:) = [-1.0,  0.0, -1.0, -1.0,  1.0,  sq2,  1.0,  0.0,  sq1,  sq1];

for iter=1:2
    for i=xmin:xmax
        for j=zmin:zmax
            a11 = A11(i,j); 
            a33 = A33(i,j); 
            a13 = A13(i,j); 
            a44 = A44(i,j);
            K = GetK(a11,a33,a13,a44);
            for triang = 1:8
                Ta = tti(i+tr(triang, 1), j+tr(triang, 2));
                Tb = tti(i+tr(triang, 3), j+tr(triang, 4));
                b = tr(triang, 5)*step;
                a = tr(triang, 6)*step;
                N = tr(triang, 7:10);
                TC = GetTC(Ta, Tb, a, b, N, K);
                for r = 1:4
                    TCr = TC(r);  
                    if (TCr<tmax)
                        Vg = CheckVg(Ta, Tb, a, b, N, K, TCr); 
                        if Vg==0
                            TC(r) = tmax;
                        end
                    end
                end
                TCn = min(real(TC));
                Vga = GetVg(a11,a33,a13,a44,N(1));          % N(1) = n11 = sin(theta)
                TCa = Ta + b/Vga;
                Vgb = GetVg(a11,a33,a13,a44,N(3));          % N(3) = n21 = sin(theta)
                TCb = Tb + a/Vgb;
                TCo = tti(i,j);
                tti(i, j) = min([TCn, TCa, TCb, TCo]);        
            end
        end
    end
    %2
    for i=xmin:xmax
        for j=zmax:-1:zmin
        a11 = A11(i,j); 
            a33 = A33(i,j); 
            a13 = A13(i,j); 
            a44 = A44(i,j);
            K = GetK(a11,a33,a13,a44);
            for triang = 1:8
                Ta = tti(i+tr(triang, 1), j+tr(triang, 2));
                Tb = tti(i+tr(triang, 3), j+tr(triang, 4));
                b = tr(triang, 5)*step;
                a = tr(triang, 6)*step;
                N = tr(triang, 7:10);
                TC = GetTC(Ta, Tb, a, b, N, K);
                for r = 1:4
                    TCr = TC(r);  
                    if (TCr<tmax)
                        Vg = CheckVg(Ta, Tb, a, b, N, K, TCr); 
                        if Vg==0
                            TC(r) = tmax;
                        end
                    end
                end
                TCn = min(real(TC));
                Vga = GetVg(a11,a33,a13,a44,N(1));          % N(1) = n11 = sin(theta)
                TCa = Ta + b/Vga;
                Vgb = GetVg(a11,a33,a13,a44,N(3));          % N(3) = n21 = sin(theta)
                TCb = Tb + a/Vgb;
                TCo = tti(i,j);
                tti(i, j) = min([TCn, TCa, TCb, TCo]);        
            end
        end
    end
    % 3
    for i=xmax:-1:xmin
        for j=zmin:zmax
        a11 = A11(i,j); 
            a33 = A33(i,j); 
            a13 = A13(i,j); 
            a44 = A44(i,j);
            K = GetK(a11,a33,a13,a44);
            for triang = 1:8
                Ta = tti(i+tr(triang, 1), j+tr(triang, 2));
                Tb = tti(i+tr(triang, 3), j+tr(triang, 4));
                b = tr(triang, 5)*step;
                a = tr(triang, 6)*step;
                N = tr(triang, 7:10);
                TC = GetTC(Ta, Tb, a, b, N, K);
                for r = 1:4
                    TCr = TC(r);  
                    if (TCr<tmax)
                        Vg = CheckVg(Ta, Tb, a, b, N, K, TCr); 
                        if Vg==0
                            TC(r) = tmax;
                        end
                    end
                end
                TCn = min(real(TC));
                Vga = GetVg(a11,a33,a13,a44,N(1));          % N(1) = n11 = sin(theta)
                TCa = Ta + b/Vga;
                Vgb = GetVg(a11,a33,a13,a44,N(3));          % N(3) = n21 = sin(theta)
                TCb = Tb + a/Vgb;
                TCo = tti(i,j);
                tti(i, j) = min([TCn, TCa, TCb, TCo]);        
            end
        end
    end
    %4
    for i=xmax:-1:xmin
        for j=zmax:-1:zmin
        a11 = A11(i,j); 
            a33 = A33(i,j); 
            a13 = A13(i,j); 
            a44 = A44(i,j);
            K = GetK(a11,a33,a13,a44);
            for triang = 1:8
                Ta = tti(i+tr(triang, 1), j+tr(triang, 2));
                Tb = tti(i+tr(triang, 3), j+tr(triang, 4));
                b = tr(triang, 5)*step;
                a = tr(triang, 6)*step;
                N = tr(triang, 7:10);
                TC = GetTC(Ta, Tb, a, b, N, K);
                for r = 1:4
                    TCr = TC(r);  
                    if (TCr<tmax)
                        Vg = CheckVg(Ta, Tb, a, b, N, K, TCr); 
                        if Vg==0
                            TC(r) = tmax;
                        end
                    end
                end
                TCn = min(real(TC));
                Vga = GetVg(a11,a33,a13,a44,N(1));          % N(1) = n11 = sin(theta)
                TCa = Ta + b/Vga;
                Vgb = GetVg(a11,a33,a13,a44,N(3));          % N(3) = n21 = sin(theta)
                TCb = Tb + a/Vgb;
                TCo = tti(i,j);
                tti(i, j) = min([TCn, TCa, TCb, TCo]);        
            end
        end
    end
end

tvltm = tti(2:end-1, 2:end-1);
   
%% Define coefficients

function K = GetK(A11,A33,A13,A44)

K = zeros(5, 1);

K(1) = -A11*A44;
K(2) = -A33*A44;
K(3) = -A11*A33 - A44^2 + (A13+A44)^2;
K(4) =  A11 + A44;
K(5) =  A33 + A44;
                
%% Find roots of Hamiltonian

function TC = GetTC(Ta, Tb, a, b, N, K)
        
cosg = sqrt(2)/2;  
tmax = 1.e20;
W = zeros(5, 1);
    
n11 = N(1); n12 = N(2); n21 = N(3); n22 = N(4);
k1 = K(1); k2 = K(2); k3 = K(3); k4 = K(4); k5 = K(5);
        
p11 = 2*(n11 - n21*cosg);       % 2 == 1/sin^2(gamma)
p12 = 2*(n21 - n11*cosg);
p21 = 2*(n12 - n22*cosg);
p22 = 2*(n22 - n12*cosg);
    
g1 = p11/b + p12/a;
g2 = -(p11/b*Ta + p12/a*Tb);
g3 = p21/b + p22/a;
g4 = -(p21/b*Ta + p22/a*Tb);

W(1) =   g1*g1*g1*g1*k1 +   g3*g3*g3*g3*k2 +    g1*g1*g3*g3*k3;                                             
W(2) = 4*g1*g1*g1*g2*k1 + 4*g3*g3*g3*g4*k2 + (2*g1*g1*g3*g4 + 2*g1*g2*g3*g3)*k3;
W(3) = 6*g1*g1*g2*g2*k1 + 6*g3*g3*g4*g4*k2 + (  g1*g1*g4*g4 + 4*g1*g2*g3*g4 + g2*g2*g3*g3)*k3 + g1*g1*k4 + g3*g3*k5;
W(4) = 4*g1*g2*g2*g2*k1 + 4*g3*g4*g4*g4*k2 + (2*g1*g2*g4*g4 + 2*g2*g2*g3*g4)*k3 + 2*g1*g2*k4 + 2*g3*g4*k5;
W(5) =   g2*g2*g2*g2*k1 +   g4*g4*g4*g4*k2 +    g2*g2*g4*g4*k3 + g2*g2*k4 + g4*g4*k5 - 1;
    
%disp(['W(1) =' num2str(W(1)), ' W(2) =' num2str(W(2)) ' W(3) =' num2str(W(3)), ' W(4) =' num2str(W(4)), ' W(5) =' num2str(W(5))])
TC = roots(W);
Tm = min(Ta, Tb);
for i=1:4
    iTC = imag(TC(i));
    rTC = real(TC(i));
    if (iTC~=0 || rTC<Tm)
        TC(i) = tmax;
    end
end
          
%% Check group velocity
    
function Vg = CheckVg(Ta, Tb, a, b, N, K, Tc)

Vg=0;
cosg = sqrt(2)/2;      
        
n11 = N(1); n12 = N(2); n21 = N(3); n22 = N(4);
k1 = K(1); k2 = K(2); k3 = K(3); k4 = K(4); k5 = K(5);
        
p11 = 2*(n11 - n21*cosg);       % 2 == 1/sin^2(gamma)
p12 = 2*(n21 - n11*cosg);
p21 = 2*(n12 - n22*cosg);
p22 = 2*(n22 - n12*cosg);

g1 = p11/b + p12/a;
g2 = -(p11/b*Ta + p12/a*Tb);
g3 = p21/b + p22/a;
g4 = -(p21/b*Ta + p22/a*Tb);
    
p1 = g1*Tc + g2;
p3 = g3*Tc + g4;
  
V(1) = 4*k1*p1^3 + 2*k3*p1*p3^2 + 2*k4*p1;
V(2) = 4*k2*p3^3 + 2*k3*p1^2*p3 + 2*k5*p3;
    
V =  V/norm(V); 
if (min(n11, n21) < V(1) && max(n11, n21)> V(1) && min(n12, n22) < V(2) && max(n12, n22) > V(2))
    Vg=1;
end
    
%% Find group velocity with group angle

function Vg = GetVg(a11,a33,a13,a44,ST)

if abs(ST) == 0 || abs(ST) == 1;
    % phase and group velocities coinside in horizontal and vertical
    % directions, we can use formula for phase velocity
    
    D = sqrt( (a33-a44)^2 + 2*(2*(a13+a44)^2 - (a33-a44)*(a11+a33-2*a44))*ST.^2 ...
         + ((a11+a33-2*a44)^2 - 4*(a13+a44)^2)*ST.^4);
    Vg = sqrt(0.5*(a33 + a44  + (a11 - a33)*ST.^2 + D));                    % in this case group == phase
else
    % See script "Analize Relationship between phase and group velocity"
    % See paper of Vladimir Grechka (2013), Ray-direction velocities in VTI
    % media, Geophysics, V. 78, N. 1, pp. F1-F5
    STini = ST; 
    psi = asin(STini);
    PSI = tan(psi); 
    
    % Thomsen parameters
    eps = (a11-a33)/(2*a33); 
    del = ((a13+a44)^2 - (a33-a44)^2)/(2*a33*(a33-a44)); 

    % F matrix, formulas A1-A6 in paper
    F11 = 1;
    F12 = -a33 -a44; 
    F13 = a33*a44; 
    F21 = -a11 - a44; 
    F22 = 2*a33*((eps-del)*a33 + (1+del)*a44);
    F31 = a11*a44; 

    % M matrix, formulas B1-B10 in paper
    A = F12^2 - 4*F11*F13; 
    B = F21^2 - 4*F11*F31; 

    M = zeros(7,3);
    M(1,3) = - A*F13; 
    M(2,2) = A*F22; 
    M(3,1) = F11*F22^2 + F13*F21^2 - F12*F21*F22;
    M(3,3) = 4*F11*F13*F22 - 2*F12*F13*F21; 
    M(4,2) = 2*(F12^2*F31 + F13*(F21^2-4*F11*F31)-F11*F22^2);
    M(5,1) = 4*F11*F22*F31 - 2*F12*F21*F31; 
    M(5,3) = F11*F22^2 + F12^2*F31 - F12*F21*F22; 
    M(6,2) = B*F22; 
    M(7,1) = -B*F31; 
    
    Mpsi =  M*[1; PSI; PSI^2]; 
    THETA = roots(Mpsi(end:-1:1));      % different order in roots function
    theta = atan(THETA); 
    theta = theta(imag(theta)==0);      % take only real roots
   
    VpQP  = zeros(size(theta)); 
    VgQP  = zeros(size(theta)); 
    psiQP = zeros(size(theta)); 
    
    for i=1:length(theta)   % for each real root find
        
        % phase velocity
        dtheta = 0.001; 
        ST = sin(theta(i)-dtheta); 
        D = sqrt( (a33-a44)^2 + 2*(2*(a13+a44)^2 - (a33-a44)*(a11+a33-2*a44))*ST.^2 ...
            + ((a11+a33-2*a44)^2 - 4*(a13+a44)^2)*ST.^4);
        VpQPm = sqrt(0.5*(a33 + a44  + (a11 - a33)*ST.^2 + D));
    
        ST = sin(theta(i)+dtheta); 
        D = sqrt( (a33-a44)^2 + 2*(2*(a13+a44)^2 - (a33-a44)*(a11+a33-2*a44))*ST.^2 ...
            + ((a11+a33-2*a44)^2 - 4*(a13+a44)^2)*ST.^4);
        VpQPp = sqrt(0.5*(a33 + a44  + (a11 - a33)*ST.^2 + D));
    
        VpQP(i) = (VpQPp+VpQPm)/2; 
    
        dVpQP = (VpQPp-VpQPm)/(2*dtheta); 

        % group angle
        psiQP(i) = (atan(dVpQP./VpQP(i))+theta(i));
    
        % group velocity
        VgQP(i) = sqrt(VpQP(i).^2 + dVpQP.^2); 
    end
    
    Vg = VgQP(abs(sin(psiQP)-STini)<0.00001); 
    if length(Vg)>1 
        Vg = Vg(1);
    end
end







