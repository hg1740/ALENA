classdef System < matlab.mixin.SetGet
    %System Describes a collection of awi.methods.isfe.SlenderBody objects in
    %terms of its states and orientation.
    %
    % Detailed Description:
    %   - A system can comprise of a number of "bodies" which themselves
    %     are described by a set of beam elements - an example of a body is 
    %     a fuselage, strut or lifting surface element.
    %   - The 'System' object handles the 'global' properties of the
    %     system. These are properties that affect the collection of
    %     SlenderBody objects that describe the system. An example of this
    %     would be the global-to-aircraft rotation matrix.
    
    %System states
    properties (Dependent)        
        %Column vector of forces and moments at each beam node
        BeamLoads
    end
    
    %System states (internal)
    properties (Hidden = true)
        %Column vector of forces and moments
        x_f
    end
    
    %System model
    properties (SetAccess = private)
        %Collection of 'awi.methods.isfe.Model' describing the system
        Model
    end
    
    %Transformation matrices
    properties (SetAccess = private)
        %Global-to-aircraft transformation matrix
        CGa = eye(3);
    end
        
    %Indexing the system matrices
    properties (SetAccess = private)
        %Upper and lower bounds for indexing the different bodies in the
        %system
        BodyLocalLoadIndex
    end
    
    methods % set / get
        function val = get.BeamLoads(obj) %get.BeamLoads 
            val = obj.x_f;
        end
    end
    
    methods % construction
        function obj = System(ModelData, Options)
            %System Constructor for the 'awi.methods.isfe.System' class.
            %
            % Actions performed:
            %   1. If an instance of ''awi.fe.isfe.Model'' has been passed 
            %   as the first argument then it will be passed to the
            %   ''initialise'' method.
            
            if nargin < 2
                ModelData = [];
                Options   = [];
            end
            if isa(ModelData, 'awi.methods.isfe.SlenderBody') && isa(Options, 'awi.methods.Options')
                initialise(obj, ModelData, Options);
            end
            
            %Stash the model
            obj.Model = ModelData;
            
        end
    end
    
    methods % setting up the default system states
        function initialise(obj, ModelData, AnalysisOptions)
            %initialise Sets up the initial system states.
            
            %Get the number of beam elements in each part
            nElem = [ModelData.NumElem];
            assert(numel(nElem) == numel(ModelData), ['Missing data in ', ...
                'the model description. Check that the ''awi.methods.isfe.Model'' ', ...
                'object(s) have been correctly initialised.']);
            
            %Define bounds for indexing
            ub = cumsum(nElem * 6);
            lb = [1, ub(1 : end - 1) + 1];
            obj.BodyLocalLoadIndex = [lb ; ub];            
            
            %Initialise states
            initialiseStates(ModelData, AnalysisOptions);   
            
            %Set up global state vector
            obj.x_f = vertcat(ModelData.x_f);

        end
    end
    
end

