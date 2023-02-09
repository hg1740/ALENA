classdef ScalarProp < awi.fe.FEBaseClass 
    %ScalarProp Defines a property of a scalar mass, damper or spring.
    
   %Primary properties
    properties
        %Property identification number
        PID
    end
    
    %Shadown properties
    properties (Hidden)
        %Actual property value is hidden away
        PropVal
        %Name of the property
        
        %Nastran property token
    end
    
    methods % set / get        
        function set.PID(obj, val)  %set.PID
            %set.PID Set method for the property 'PID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.PropVal(obj, val) %set.PropVal
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'real', 'nonnan'}, class(obj), 'K');
            obj.PropVal = val;
        end
    end
    
end

