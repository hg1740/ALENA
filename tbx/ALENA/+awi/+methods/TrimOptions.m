classdef TrimOptions < awi.methods.Options
    %TrimOptions Provides the basic simulation parameters for a aircraft
    %trim analysis using any of the aeroelastic methods in the AWI package.
    
    %Solver parameters
    properties
        %Incremement in trim variable for calculating numerical jacobian of
        %the trim system
        TrimVariableStepSize = 1e-3;        
    end
    
    %Trim parameters
    properties
        %Vector defining the global loads which will be balanced in the
        %trim solution
        ForceBalanceIndex = [];
        %ControlSurfDOF
        %RigidBodyDOF   
    end
    
    properties (Dependent)
        %A [N, 6] matrix of zeros and ones which extracts the trim loads
        %from the [6, 1] global loads matrix - N is number of trim ls.
        TrimLoadsIndexMatrix
    end
    
    methods % set / get
        function set.TrimVariableStepSize(obj, val)     %set.TrimVariableStepSize
            %set.TrimVariableStepSize Set method for the property
            %'TrimVariableStepSize'.
        
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'nonnan', 'real'}, class(obj), 'TrimVariableStepSize');
            obj.TrimVariableStepSize = val;            
        end
        function set.ForceBalanceIndex(obj, val)        %set.ForceBalanceIndex
            %set.ForceBalanceIndex Set method for the property
            %'ForceBalanceIndex'.
            
            validateattributes(val, {'numeric'}, {'vector', 'integer', ...
                'positive', '<' 7}, class(obj), 'ForceBalanceIndex');
            obj.ForceBalanceIndex = val;
        end
        function val = get.TrimLoadsIndexMatrix(obj)    %get.TrimLoadsIndexMatrix
            %get.TrimLoadsIndexMatrix Get method for the dependent property
            %'TrimLoadsIndexMatrix'.
            
            val = [];
            
            %Handle the case where the user has not defined the trim loads            
            trimForceIndex = obj.ForceBalanceIndex;
            if isempty(trimForceIndex) %Escape route
                return
            end
            
            %Remove unwanted rows
            val   = eye(6);
            index = setdiff(1 : 6, trimForceIndex);
            val(index, :) = []; 
            
        end
    end
    
end

