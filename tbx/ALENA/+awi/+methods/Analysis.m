classdef Analysis < matlab.mixin.SetGet %(Abstract) --> If the object just holds methods then we probably want it to be abstract...
    %Analysis Provides the basic method blocks for any generic analysis
    %in the AWI package.
        
    %Top level properties
    properties
        %Handle to the 'awi.fe.FEModel' that contains the data for analysis
        AnalysisModel
    end
    
    properties (Dependent, Hidden = true)
        %Flatlist of the analysis model hierarchy
        FlatModel
    end
    
    methods % set / get
        function set.AnalysisModel(obj, val) %set.AnalysisModel 
            %set.AnalysisModel Set method for the property 'AnalysisModel'.
            %
            % 'AnalysisModel' must be a valid, scalar instance of the
            % 'awi.fe.FEModel' class.
            validateattributes(val, {'awi.fe.FEModel'}, {'scalar', ...
                'nonempty'}, class(obj), 'AnalysisModel');
            obj.AnalysisModel = val;            
        end
        function val = get.FlatModel(obj)    %get.FlatModel 
           %get.FlatModel Get method for the dependent property 'FlatModel'.
           
           val = [];
           if isempty(obj.AnalysisModel)
               return
           end
           val = flatlist(obj.AnalysisModel);
        end
    end
    
    methods % static analysis
        function StaticResults = static(obj, StaticOptions)
            %static Executes a static analysis using the analysis model
            %assigned to this method object and the options defined in
            %'StaticOptions'.
            
            StaticResults = [];
            
            %Check we can perform the analysis
            parse(obj, 'Analysis', 'static');
            
        end
    end
    
    methods % aeroelastic analysis
        function TrimResults = trim(obj, Aircraft, LoadCases)
            %trim Checks that the inputs to the 'trim' method are valid and
            %returns a scalar logical indicating whether the analysis can
            %continue. The actual trim analysis is conducted at the
            %subclass level.
            
            TrimResults = [];
            
            assert(isa(Aircraft, 'awi.model.Aircraft'), ['Unable to ', ...
                'perform the trim analysis without a valid '         , ...
                '''awi.model.Aircraft'' object.']);
            assert(isa(LoadCases, 'awi.model.LoadCase'), ['Unable to ', ...
                'perform the trim analysis without a valid set of'    , ...
                '''awi.model.LoadCase'' objects.']);
            
            %Check that the objects have the correct data. Invoke 'canTrim'
            
        end
        function discreteGust(obj)
            
        end
        function continuousGust(obj)
            
        end
        function flutter(obj)
            
        end
    end
    
    methods % multi-disciplinary analysis
        function sizing(obj)
            
        end
        function optimisation(obj)
            
        end        
    end
    
    methods (Access = protected)
        function parse(obj, varargin)
            %parse Checks the object has all the necessary information to
            %proceeed with the analysis.
            %
            % Generic Requirements:
            %   1. The ''AnalysisModel'' property must be populated.
            %   2. The analysis model must have the following properties
            %   defined:
            %       - Nodes
            %       - Beams
            %       - BeamProps
            %       - Materials
            %
            % Static Requirements
            %   1. The FEM must have a load defined.
            
            %Special options
            p = inputParser;
            addParameter(p, 'Analysis', @(x)any(validatestring(x, {'static'})));
            parse(p, varargin{:});
            
            %Get the flattened FEM
            AllFEM = obj.FlatModel;
            
            %Generic 
            assert(~isempty(obj.AnalysisModel), ['The ''AnalysisModel'' ', ...
                'property has not been defined. Assign a valid '         , ...
                '''awi.fe.FEModel'' object and continue with the analysis.']);
            i_checkFEObj(AllFEM, 'Nodes'    , 'awi.fe.Node');
            i_checkFEObj(AllFEM, 'Beams'    , 'awi.fe.Beam');
            i_checkFEObj(AllFEM, 'BeamProps', 'awi.fe.BeamProp');
            i_checkFEObj(AllFEM, 'Materials', 'awi.fe.Material');

            %Solution specific            
            switch lower(p.Results.Analysis)
                
                case 'static'
                    i_checkFEObj(AllFEM, 'PointLoads', 'awi.fe.PointLoad');
                    
            end
            
            %Local functions
            function i_checkFEObj(FEM, prp, nam) 
                assert(~isempty([FEM.(prp)]), sprintf(['No ''%s'' ', ...
                    'objects found in the analysis model. Update the ', ...
                    'model and run the analysis again.', nam]));
            end
            
        end
    end
    
end

