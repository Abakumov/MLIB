function text2line(h,ksi,z,T)
% Inserts text T in/near line with handle h
%  ksi - relative distance from the beginning of curve,
%  z - shift along normal to curve
%
set(gcf, 'CurrentObject', h)
x=h.XData;
y=h.YData;
i = round(ksi*numel(x));
% Get the local slope
dy=y(i+1)-y(i-1);
dx=x(i+1)-x(i-1);
d = dy/dx;
X = diff(get(gca, 'xlim'));
Y = diff(get(gca, 'ylim'));
p = pbaspect;
a = atan(d*p(2)*X/p(1)/Y)*180/pi;
% Display the text
switch z==0
    case 1
        text(x(i), y(i), T,'HorizontalAlignment','center', 'BackgroundColor', 'w', 'rotation', a);
    case 0
        ez=[dy,-dx]/norm([dy,-dx]); % unit normal vector
        text(x(i)+z*ez(1), y(i)+z*ez(2), T, 'HorizontalAlignment','center', 'rotation', a);
end