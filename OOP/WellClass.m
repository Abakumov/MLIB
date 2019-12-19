classdef WellClass < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        wellhead_X
        wellhead_Y
        wellhead_Z
        wellhead_lon
        wellhead_lat
        MD              % Measured Depth
        INC             % Inclination
        AZM             % Azimuth
        TVD             % True vertical depth
        NS              % +N / -S
        EW              % +E / -W
        VS              % Ver. Sect.
        DLS             % dogleg severity 
        RHOB            % Density
        DT              % P wave velocity
        SDT             % S wave velocity
        xx
        yy
        zz
        dd
        lat
        lon
     end
    
    methods
        %function valid = checkEvent(obj)
        %    if obj.ID ~=0
        %        valid = 1; 
        %    else 
        %        valid = 0; 
        %    end
        %end
              
        function setParameters(obj)
             obj.dd = 1:max(round(obj.MD));
             md = obj.MD;
             x =  obj.wellhead_X + obj.EW;
             y =  obj.wellhead_Y + obj.NS;
             z = -obj.wellhead_Z + obj.TVD;
             obj.xx = interp1(md,x,obj.dd);
             obj.yy = interp1(md,y,obj.dd);
             obj.zz = interp1(md,z,obj.dd);
             
        end
    end
end

