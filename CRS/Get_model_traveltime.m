function [TT] = Get_model_traveltime(xs,xg,XX,YY,ZZ,model)
% This function returns traveltime
% Abakumov Ivan
% 27th August 2016

    Get_model_parameters;

    % constant velocity
    if  model == 11 || model == 21 || model == 31 || model == 41 || model == 51 || model == 61 || model == 0
        ts = sqrt((xs(1)-XX).^2+(xs(2)-YY).^2+(xs(3)-ZZ).^2)/v0;
        tg = sqrt((xg(1)-XX).^2+(xg(2)-YY).^2+(xg(3)-ZZ).^2)/v0;
        TT = ts + tg; 

    % linear gradient
    elseif model == 12 || model == 22 || model == 32 || model == 42 || model == 52 || model == 62
       Vs = v0 + k*xs(3);                                   % velocity at the source point 
       Vg = v0 + k*xg(3);                                   % velocity at the receiver point
       Vr = v0 + k*ZZ;                                      % velocity at the point on reflector
       Rs = (xs(1)-XX).^2+(xs(2)-YY).^2+(xs(3)-ZZ).^2;      % distance from the source to the reflector point
       Rg = (xg(1)-XX).^2+(xg(2)-YY).^2+(xg(3)-ZZ).^2;      % distance from the receiver to the reflector point
       ts = 1/k*acosh(1 + k^2*Rs./(2*Vs*Vr));               % traveltime in the medium with linear gradient of velocity 
       tg = 1/k*acosh(1 + k^2*Rg./(2*Vg*Vr)); 
       TT  = ts + tg;
       
    % constant layer + linear gradient   
    elseif model == 13 || model == 23 || model == 33 || model == 43 || model == 53 || model == 63
        TT = zeros(size(XX)); 
        Rs = sqrt((xs(1)-XX).^2+(xs(2)-YY).^2);
        Rg = sqrt((xg(1)-XX).^2+(xg(2)-YY).^2);
        for i=1:size(XX,1)
            for j=1:size(XX,2); 
                zz = ZZ(i,j); 
                Vr = v0 + k*(zz-z0);                              % velocity at the point on reflector

                ts = @(r)(sqrt( r.^2  + (xs(3)-z0).^2))/v0 + ...
                    1/k*acosh(1 + k^2*((Rs(i,j)-r).^2 + (z0-zz).^2)./(2*v0*Vr));            % traveltime in the medium
                
                tg = @(r)(sqrt( r.^2  + (xg(3)-z0).^2))/v0 + ...
                    1/k*acosh(1 + k^2*((Rg(i,j)-r).^2 + (z0-zz).^2)./(2*v0*Vr));            % traveltime in the medium
                
                TT(i,j) = ts(fminbnd(ts, 0, Rs(i,j))) + tg(fminbnd(tg, 0, Rg(i,j)));
            end
        end
    elseif model == 14 || model == 24 || model == 34 || model == 44 || model == 54 || model == 64 || model > 100
        ts = sqrt((xs(1)-XX).^2/A11 + (xs(2)-YY).^2/A22  + (xs(3)-ZZ).^2/A33);
        tg = sqrt((xg(1)-XX).^2/A11 + (xg(2)-YY).^2/A22  + (xg(3)-ZZ).^2/A33);
        TT = ts + tg; 
    end
end