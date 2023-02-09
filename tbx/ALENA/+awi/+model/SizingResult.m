classdef (ConstructOnLoad) SizingResult < awi.model.AeroelasticResult
    %SizingResult Defines a set of sizing results over a beam.
    %
    % 'awi.model.SizingResult' derives from 'awi.model.AeroelasticResult',
    % therefore, it contains all results quantities in that superclass (and 
    % any other superclasses) as well as the following additional
    % quantities:
    %
    %   - USknT   : Box upper skin thickness
    %   - LSknT   : Box lower skin thickness
    %   - UShellT : Box upper shell thickness (skin thickness + smeared
    %               area of the stringers)
    %   - LShellT : Box lower shell thickness (skin thickness + smeared
    %               area of the stringers)
    %   - FSparT  : Front spar thickness
    %   - RSparT  : Rear spar thickness
        
    methods % construction
        function obj = SizingResult(varargin)            
           
            %Pass it on
            obj@awi.model.AeroelasticResult(varargin{:});
            
            %Add results quantities
            addBeamProperty(obj, 'USknT'  , 'Type', 'Box Dimensions');
            addBeamProperty(obj, 'LSknT'  , 'Type', 'Box Dimensions');
            addBeamProperty(obj, 'UShellT', 'Type', 'Box Dimensions');
            addBeamProperty(obj, 'LShellT', 'Type', 'Box Dimensions');
            addBeamProperty(obj, 'FSparT' , 'Type', 'Box Dimensions');
            addBeamProperty(obj, 'RSparT' , 'Type', 'Box Dimensions');
 
        end
    end
    
end

