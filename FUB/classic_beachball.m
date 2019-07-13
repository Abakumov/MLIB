function classic_beachball(u,varargin)
%% Function plots a beachball of the amplitude-field up. 
% The amplitude-field should have the form: three components of displacement, 
% phi =  azimuth = [0:b:360]° and dip = theta = [0:a:180]° .  
% The parameters a and b should be positive integers
% Optional parameters are the position of the beachball and the beachball-radius

% © Nepomuk Boitz, FU Berlin (boitz@geophysik.fu-berlin.de) July 2017

if nargin == 1
    pos = [0,0,0]; 
    scale = 1; 
end
if nargin == 2
    pos = varargin{1} ;
    scale = 1; 
end
if nargin == 3
    pos = varargin{1} ; 
    scale = varargin{2}; 
end
x_z_scale = 1; 
if nargin >=4
    pos = varargin{1} ; 
    scale = varargin{2}; 
    x_z_scale = varargin{3}; 
end
%thet =  varargin{3} 
%phi_t = varargin{4} 
%%
a = 180/size(u,3);
b = 360/size(u,2);
hold on; 
for theta=1:1:size(u,3)
    for phi=1:1:size(u,2)
        x = zeros(3,1); 
        x(1) = sin(a*theta/180*pi)*cos(b*phi/180*pi -pi/2);
        x(2) = sin(a*theta/180*pi)*sin(b*phi/180*pi -pi/2);
        x(3) = -cos(a*theta/180*pi);   
        %sqrt((u(1,phi,theta)./norm(u(:,phi,theta))+x(1)).^2 + (u(2,phi,theta)./norm(u(:,phi,theta))+x(2)).^2 + (u(3,phi,theta)./norm(u(:,phi,theta))+x(3)).^2)
        if sqrt((u(1,phi,theta)./norm(u(:,phi,theta))+x(1)).^2 + (u(2,phi,theta)./norm(u(:,phi,theta))+x(2)).^2 + (u(3,phi,theta)./norm(u(:,phi,theta))+x(3)).^2) < 1
            xq(theta,phi,:) = -1; 
        else
            xq(theta,phi,:) = 1; 
        end
        %plot3([0 x(1)],[0 x(2)],[0 x(3)],'r-')
    end
end
[x,y,z] = sphere(size(u,3)-1);
% x(37:361,:) = []; 
% y(37:361,:) = []; 
% z(37:361,:) = []; 
%xq(50:360,:) = []; 


sh = surf(pos(1) + scale * -x * x_z_scale ,pos(2) + scale * y * x_z_scale,pos(3) + scale * z,xq,'linestyle','none'); 
map = [1 1 0;0 0 0];
colormap(map)
axis equal
xlabel('Easting [m]')
ylabel('Northing [m]')
%view([0 0 1])