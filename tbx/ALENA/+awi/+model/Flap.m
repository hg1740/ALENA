classdef (ConstructOnLoad) Flap < awi.model.ControlSurface
    %Flap Defines an instance of a trailing edge control surface known as a
    %"Flap".
    %
    % Very thin wrapper around the 'awi.model.ControlSurface' class. The
    % only reason this class exists is to allow control surfaces of
    % different types to be grouped together in the 'Tree' view.
    %
    % The parameterisation of the 'Flap' object is no different to the 
    % generic 'ControlSurface' class. This class simply defines a different
    % hinge line location and different colours for visualisation purposes.

    methods % construction
        function obj = Flap(varargin)
            
            %Pass it on to the superclass
            obj@awi.model.ControlSurface(varargin{:});
            
            %Set up some new defaults for plotting & hinge location
            obj.HingeLine = 'LE';
            obj.FaceColor = 'g';
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht)
           
            tag = 'Flaps';
            
            %Pass it on
            hg = drawElement@awi.model.ControlSurface(obj, ht, tag);
            
        end
    end
    
end

