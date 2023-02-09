classdef (ConstructOnLoad) Slat < awi.model.ControlSurface
    %Slat Defines an instance of a leading edge control surface known as a
    %"Slat".
    %
    % Very thin wrapper around the 'awi.model.ControlSurface' class. The
    % only reason this class exists is to allow control surfaces of
    % different types to be grouped together in the 'Tree' view.
    %
    % The parameterisation of the 'Slat' object is no different to the 
    % generic 'ControlSurface' class. This class simply defines a different
    % hinge line location and different colours for visualisation purposes.

    methods % construction
        function obj = Slat(varargin)
            
            %Pass it on to the superclass
            obj@awi.model.ControlSurface(varargin{:});
            
            %Set up some new defaults for plotting & hinge location
            obj.HingeLine = 'TE';
            obj.FaceColor = 'b';
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht)
           
            tag = 'Slats';
            
            %Pass it on
            hg = drawElement@awi.model.ControlSurface(obj, ht, tag);
            
        end
    end
    
end

