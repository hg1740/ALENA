classdef (ConstructOnLoad) TransientResultSet < awi.model.ResultSet
    %TransientResult Describes a results quantity that is defined in the
    %time or frequency domain.
    
    properties %DynamicUnits
        %Vector of time or frequency values
        TimeFreq
        %Units of meausre for the dynamic quantity
        DynamicUnits = 's';
    end
    
    properties (Dependent) %Label
        %Label for the x-axis
        Label
    end
    
    properties (Constant) %ValidDynamicUnits
        %Valid inputs for the property 'DynamicUnits'.
        ValidDynamicUnits = {'s', 'hz', 'rads'};
    end
    
    methods % set / get
        function set.DynamicUnits(obj, val) %set.DynamicUnits
            %set.DynamicUnits Set method for the property 'DynamicUnits'.
            %
            %   - 'DynamicUnits' must be one of the valus defined by
            %     'ValidDynamicUnits'.
            
            validatestring(val, obj.ValidDynamicUnits, class(obj), 'DynamicUnits');
            obj.DynamicUnits = val;
        end        
        function val = get.Label(obj)
            %get.Label Get method for the dependent property 'Label'.
            
            val = [' [', obj.DynamicUnits, ']'];
            
        end
    end
    
    methods % constructor
        function obj = TransientResultSet(varargin)
            %TransientResultSet Constructor for the class
            %'TransientResultSet'.
            
            %Pass it on
            obj@awi.model.ResultSet(varargin{:});
            
        end
    end
    
end

