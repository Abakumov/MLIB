function p = compute_potency(normal,slip,varargin)
%% Computes potency tensor from 
% nargin == 2: normal- and slipvector
% nargin == 3: dip, azimuth, slip direction
% nargin == 4: dip, azimuth, slip direction, angle of deviation from 90° 
% between n and s

% © Nepomuk Boitz, Oktober 2016
% boitz@geophysik.fu-berlin.de


Ry = @(alpha) [cosd(alpha) 0 sind(alpha); 0 1 0 ; -sind(alpha) 0 cosd(alpha)];
Rx = @(alpha) [1 0 0; 0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];
Rz = @(alpha) [cosd(alpha) -sind(alpha) 0; sind(alpha) cosd(alpha) 0; 0 0 1];

if nargin == 2
    n = normal; 
    s = slip; 
end
if nargin > 2
    theta = normal; 
    phi = slip; 
    slip = varargin{1}; 
    n = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)];
    s = (Rz(phi)*Ry(theta)*[cosd(slip) sind(slip) 0]');    
end
if nargin > 3
    non_dc = varargin{2}; 
    M_rot = makehgtform('axisrotate',cross(n,s),non_dc*pi/180);
    s = M_rot(1:3,1:3)*s; 
end
p = zeros(3,3); 
for i=1:3
	for j=1:3
		p(i,j) = 0.5*(n(i)*s(j)+n(j)*s(i)); 
    end
end

