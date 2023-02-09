classdef BeamPropertyObject < matlab.mixin.SetGet
    %BeamPropertyObject Describes an object that defines properties of a
    %beam that must be calculated/derived, instead of simply prescribed as
    %a distribution.
    %
    % The actual beam properties will be defined at subclass level. This
    % object just handles the relationship between the 'BeamPropertyObject'
    % and the parent 'Beam'.    
    
    properties
        %Reference to the beam handle which this object belongs to.
        BeamHandle
        %Reference to the non-dimensional positon of this object along the
        %beam.
        BeamEta
    end
    
    methods % set / get
        function set.BeamHandle(obj, val) %set.BeamHandle
            %set.BeamHandle Set method for the property 'BeamHandle'.
            
            validateattributes(val, {'awi.model.Beam'}, {'scalar'}, ...
                class(obj), 'BeamHandle');
            obj.BeamHandle = val;
        end
        function set.BeamEta(obj, val)    %set.BeamEta
            %set.BeamEta Set method for the property 'BeamEta'.
            
            validateattributes(val, {'numeric'}, {'scalar', '<=', 1, ...
                'nonnegative'}, class(obj), 'BeamEta');
            obj.BeamEta = val;
        end
    end
    
end

