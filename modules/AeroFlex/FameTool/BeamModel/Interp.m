classdef Interp < hgsetget
    %INTERP Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties that reference objects
    properties
        Sets   = [];
        Spline  = [];
    end
    
    methods
        
        function obj = Interp(~)
            obj.Sets   = Sets;
            obj.Spline = Spline;
        end
    end
end
