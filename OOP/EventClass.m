classdef EventClass < handle
    %EVENTCLASS defines class for storage of microseismic events
    % Author: Abakumov Ivan
    % Freie Universität Berlin
    % E-mail: abakumov_ivan@mail.ru
    % Publication date: 23rd of August, 2019
    properties
        ID
        nsens
        nt
        tracesx   % [nt x nsens]
        tracesy   % [nt x nsens]
        tracesz   % [nt x nsens]
        Htraces   % [nt x nsens]
        STALTA
        seg2filename
        DateTime
        TimeZero
        PTime 
        STime 
        Azimuth 
        Dip 
        SensorID
        SensorNorthing
        SensorEasting
        SensorDepth
        Stage
        StageTime 
        Year    
    	Month    
    	Day    
    	Hours    
    	Minutes    
        Seconds    
    	Northing    
    	Easting    
    	Depth     
        % ivan's localization
        TimeZero_calc
        Ptime_calc
        Stime_calc
        MisfitP
        MisfitS
        RMSE
        MEAN
        SIGMA
    	MomMag    
    	SeiMoment
    end
    
    methods
        function valid = checkEvent(obj)
            if obj.ID ~=0
                valid = 1; 
            else 
                valid = 0; 
            end
        end
              
        function setParameters(obj)
             % check some mandatory fields
             obj.tracesx = zeros(obj.nt,obj.nsens);
             obj.tracesy = zeros(obj.nt,obj.nsens); 
             obj.tracesz = zeros(obj.nt,obj.nsens); 
          %   obj.tracesabs = zeros(obj.nt,obj.nsens); 
          %   obj.tracesLTASTA = zeros(obj.nt,obj.nsens); 
        end
    end
end

