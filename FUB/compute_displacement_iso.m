function u_p = compute_displacement_iso(M,varargin)
%% Computes displacement for isotropic propagation of P-Wave from Moment-tensor

% © Nepomuk Boitz, FU Berlin (boitz@geophysik.fu-berlin.de) March 2016

if nargin == 1
    a = 6; 
    b = 6; 
end
if nargin == 2
   a = varargin{1}(1);  
   b = varargin{1}(2); 
end
u_p = zeros(3,360/a,180/b); 
for theta = 1:180/b+1
    for phi= 1:360/a+1               
        x(1) = sin(b*theta/180*pi)*cos(a*phi/180*pi- pi/2);
        x(2) = sin(b*theta/180*pi)*sin(a*phi/180*pi- pi/2);
        x(3) = -cos(b*theta/180*pi);       
        pp = x; 
        temp = 0; 
        for j=1:3
            for k=1:3
                temp = temp + (x(j)*pp(k) + x(k)*pp(j)) * M(j,k);           
             end
        end
        u_p(:,phi,theta) = x * temp; 
    end
end
