classdef Options < matlab.mixin.SetGet
    %Options Provides the basic simulation parameters for any generic
    %analysis in the AWI package.
    
    %Top level properties
    properties
        %Type of analysis being performed (TODO - Get rid of this! Should
        %be defined in the actual analysis object)
        AnalysisType = 'static';
        %Type of structural response being modelled 
        %{'linear', 'quasi-nonlinear', 'nonlinear'}
        StructuralResponse = 'linear';
        
        % Option for controlling load iteration in SOL 400
        %{'', 'QLINEAR', 'MILDLY', 'SEVERELY'}
        CTRLDEF = '';
    end
    
    %Gravity options
    properties
        %Flag to turn gravity effects on or off
        bModelGravity = true;
        %Acceleration due to gravity in global frame
        GravityVector = [0 ; 0 ; -9.80665 ; 0 ; 0 ; 0];
    end
    
    %Newton Raphson Parameters
    properties
        %Maximum number of iterations for the Newton-Raphson method
        MaxIterations = 100;
        %Tolerance on the residual of the static equations of motion when
        %solving using the Newton-Raphson method.
        ResidualTolerance = 1e-6;
        %Number of iterations between updating the tangent matrix for the
        %Newton-Raphson method
        TangentMatrixUpdateFrequency = 5;
        %Numerical damping for controlling convergence of Newton-Raphson
        NumericalDamping = 1;        
    end
    
    %Helper properties
    properties (Constant, Hidden)
        %Valid analysis types
        ValidAnalysisTypes = {'static'};
    end
    
    methods % set / get
        function set.StructuralResponse(obj, val)           %set.StructuralResponse
            %set.StructuralResponse Set method for the property
            %'StructuralResponse'.
            
            validatestring(val, {'linear', 'quasi-nonlinear', 'nonlinear'}, ...
                class(obj), 'StructuralResponse');
            obj.StructuralResponse = val;
        end
        function set.bModelGravity(obj, val)                %set.bModelGravity
            %set.bModelGravity Set method for the property 'bModelGravity'.
            validateattributes(val, {'logical'}, {'scalar'}, ...
                class(obj), 'bModelGravity');
            obj.bModelGravity = val;
        end
        function set.GravityVector(obj, val)                %set.GravityVector
            %set.GravityVector Set method for the property 'GravityVector'.
            validateattributes(val, {'numeric'}, {'column', 'nrows', 6, ...
                'nonnan', 'finite', 'real', 'nonempty'}, class(obj), 'GravityVector');
            obj.GravityVector = val;
        end
        function set.MaxIterations(obj, val)                %set.MaxIterations
            %set.MaxIterations Set method for the property
            %'MaxIterations'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'positive'}, class(obj), 'MaxIterations');
            obj.MaxIterations = val;
        end
        function set.ResidualTolerance(obj, val)            %set.ResidualTolerance
            %set.ResidualTolerance Set method for the property
            %'ResidualTolerance'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'positive'}, ...
                class(obj), 'ResidualTolerance');
            obj.ResidualTolerance = val;
        end
        function set.TangentMatrixUpdateFrequency(obj, val) %set.TangentMatrixUpdateFrequency
            %set.TangentMatrixUpdateFrequency Set method for the property
            %'TangentMatrixUpdateFrequency'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'positive'}, class(obj), 'TangentMatrixUpdateFrequency');
            obj.TangentMatrixUpdateFrequency = val;
        end
        function set.NumericalDamping(obj, val)             %set.NumericalDamping
            %set.NumericalDamping Set method for the property
            %'NumericalDamping'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'positive'}, ...
                class(obj), 'NumericalDamping');
            obj.NumericalDamping = val;
        end
        function set.AnalysisType(obj, val)                 %set.AnalysisType
            %set.AnalysisType Set method for the property 'AnalysisType'.
            
            str = lower(val);
            validatestring(str, obj.ValidAnalysisTypes);
            obj.AnalysisType = str;
        end
        function set.CTRLDEF(obj, val)  
            % Option for controlling load iteration in SOL 400
            %{'', 'QLINEAR', 'MILDLY', 'SEVERELY'}
            validatestring(val, {'', 'QLINEAR', 'MILDLY', 'SEVERELY'});
            obj.CTRLDEF = val;
        end
    end
    
    %Each analysis should have a top-level method that is invoked by the
    %corresponding 'Analysis' view from the GUI.
    methods
        function trim(obj)
            
        end
        function discreteGust(obj)
            
        end
        function continuousGust(obj)
            
        end
        function flutter(obj)
            
        end
        function sizing(obj)
            
        end
        function optimisation(obj)
            
        end
    end
    
end

