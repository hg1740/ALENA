%% OPTIONS Input options
%
% This object contains sub-objects that can set the Geometry, Structure,
% Aerodynamic or Output options.

%   Copyright 2016 University of Bristol
%   Private function.
classdef Options
    
    properties
        Geom
        Struct
        Aero
        Outputs
    end
    
    methods
        function obj = Options(~)
            obj.Geom    = Optsgeom;
            obj.Struct  = Optstruct;
            obj.Aero    = Optsaero;
            obj.Outputs = Optsout;
        end
    end
    
end

