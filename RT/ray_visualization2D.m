function []=ray_visualization2D(tti, velmod, G, sta, color)

% Abakumov Ivan 
% FU Berlin

x = zeros(2,G.nx);


imagesc(G.xx, G.zz, velmod');
colorbar;
hold on
set(gca, 'ZDir', 'reverse');

for stan=1:size(sta,1)
    xx = sta(stan, 1);
    zz = sta(stan, 2);
    xi = x2grid(xx, G.x0, G.dx, G.nx);
    zi = x2grid(zz, G.z0, G.dz, G.nz);
    T = tti(xi,zi);
    dt = T/G.nx/1.002;
    for i=1:G.nx
        x(:,i) = [xx;zz];
        xi = x2grid(xx, G.x0, G.dx, G.nx);
        zi = x2grid(zz, G.z0, G.dz, G.nz);
        px = (tti(xi+1,zi) - tti(xi,zi))/G.dx;
        pz = (tti(xi,zi+1) - tti(xi,zi))/G.dy;
        v2=1/(px^2 + pz^2);
        xx = xx - v2*px*dt;
        zz = zz - v2*pz*dt;
  end
  plot(x(1,:),x(2,:),color);
end



