function []=ray_visualization3D(wavetype, param, stations_names, sta)

if wavetype==1
    fp22 = fopen('tmpTTp.dat', 'rb');
    aa = fread(fp22, inf, 'double');
    fclose(fp22);
    fignum=6;
end

if wavetype==2
    fp33 = fopen('tmpTTs.dat', 'rb');
    aa = fread(fp33, inf, 'double');
    fclose(fp33);
    fignum=7;
end


G3D = param(1:12);

Gox = G3D(1);
Goy = G3D(2);
Goz = G3D(3);
Gnx = G3D(4);
Gny = G3D(5);
Gnz = G3D(6);
step =  G3D(7);
tti = reshape(aa,Gnx,Gny,Gnz);

x = zeros(3,Gnx);
xax = (0:(Gnx-1))*step;
yax = (0:(Gny-1))*step;
zax = (0:(Gnz-1))*step;


%imagesc(xax, zax, tti');
% colorbar;
% if wavetype==1
%     title('P-waves');
% end
% if wavetype==2
%     title('S-waves');
% end
% xlabel('Offset [m]');
% ylabel('Depth [m]');
% axis('equal');
figure;
hold on
set(gca, 'ZDir', 'reverse');
for stan=1:length(sta);
xx = param(stan+700);
yy = param(stan+800);
zz = param(stan+900);
xi = round((xx-Gox)/step)+1;
yi = round((yy-Goy)/step)+1;
zi = round((zz-Goz)/step)+1;
T = tti(xi,yi,zi);
dt = T/Gnx/1.002;
plot3(xx,yy,zz,'r^','MarkerSize',10)
text(xx,zz-max(zax)*0.02,stations_names(stan),'color','r');
for i=1:Gnx
    x(:,i) = [xx;yy;zz];
    xi = x2grid(xx, Gox, step, Gnx);
    yi = x2grid(yy, Goy, step, Gny);
    zi = x2grid(zz, Goz, step, Gnz);
    px = (tti(xi+1,yi,zi) - tti(xi,yi,zi))/step;
    py = (tti(xi,yi+1,zi) - tti(xi,yi,zi))/step;
    pz = (tti(xi,yi,zi+1) - tti(xi,yi,zi))/step;
    v2=1/(px^2 + py^2 + pz^2);
    xx = xx - v2*px*dt;
    yy = yy - v2*py*dt;
    zz = zz - v2*pz*dt;
end
plot3(x(1,:),x(2,:),x(3,:),'black');
end
axis equal;


