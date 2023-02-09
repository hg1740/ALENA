classdef TrimAnalysis < awi.view.Analysis
    
    properties
        
        MaxIterations = 10; %passed on to getTrim
        BeamModel;
        BeamModels;
        
    end
    
    properties (Transient, SetAccess = protected )
        
        hMaxIterations
        hBeamModel;
        
    end
    
    methods % get/set
                        
        function val = get.BeamModel(obj)
            
            %Get from uicontrol
            val = get(obj.hBeamModel, 'Value');
            
            %Valid ?
            val(val < 0 | val > numel(obj.BeamModels)) = [];
            
            %Return the actual selection, rather than the index of the selection
            val = obj.BeamModels(val);
                    
        end
        
        function set.BeamModel(obj, val)
        
            %Make a note
            obj.BeamModel = val;
            
            %No selection is valid from user's point of view
            if isempty(val)
                
                %But we need something to put into popup
                idx = 1;
                
            elseif isnumeric(val)
                
                %TODO: ensure within range ?
                idx = val;
                
            elseif isa(val, 'awi.model.BeamModel')

                %Allow selection to be set with an actual member of the list
                [b, idx] = ismember(val, obj.BeamModels); %#ok<MCSUP>

                %Ensure valid
                assert(b, 'invalid beam model');
                
            end
            
            %Assign in uicontrol
            set(obj.hBeamModel, 'Value', idx); %#ok<MCSUP>
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.BeamModels(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no LoadCases
                val = [];
                
            else
                
                %What have we got ?
                val = findall(obj.Model, 'Type', 'BeamModel');
                
            end
            
        end
        
    end
    
    methods % constructor
        
        function obj = TrimAnalysis( model, varargin )
            
            %Call superclass constructor - do not pass the model in, better to wait
            obj@awi.view.Analysis( model, ...                'Title', 'Trim Analysis', ...
                'OfferLoadCaseSelection', true, ...
                'OfferMethodSelection'  , true, ...
                'PushButtons'           , {'Trim...', @analyse}, ...
                varargin{:} );
            
            %Add a label, and popup (to control selection)
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Beam Model:')
            obj.hBeamModel = uicontrol('Parent', obj.hControls, 'Style', 'popup', ...
                'Callback', @obj.onBeamModelChange);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;            
            
            % Create a edit box for MaxIterations
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'max iteration:')
            obj.hMaxIterations = uicontrol('Parent', obj.hControls,...
                'Style', 'edit',...
                'String', num2str(obj.MaxIterations), ...
                'Callback', @obj.cbMaxIterations);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;
            
            %Now store the model - which will in turn trigger an update
            obj.Model = model;
            
        end
        
    end
    
    methods %Analysis
        
        function varargout = analyse(obj, varargin)
            %analyse Runs a trim analysis using the chosen aeroelastic
            %solver and analysis options as dictated by the user.
            
            %Invoke superclass method
            %   - Clears up the view by removing any extant progress
            %     windows etc.
            analyse@awi.view.Analysis(obj); %Should return a generic method object that we can call 
                        
            %Initialise the panel for the view
            hb = uiextras.HBox('Parent', obj.hTabPanel, 'Tag', 'View');
            obj.hTabPanel.TabTitles{end} = 'View';
            obj.hTabPanel.Selection = numel(obj.hTabPanel.TabTitles);
            
            %Check for FEM - If we don't have one then generate it
            
            %Add payload & fuel masses to the FEM
            
            %Using what method/formulation?
            idxAM = contains(varargin, 'AnalysisMethod');
            if any(idxAM)
                %Go with the user value
                ind = find(idxAM) + 1;
                analysisMethod = varargin{ind};
                validatestring(analysisMethod, obj.AnalysisMethods, 'analyse', 'AnalysisMethod');
                varargin([ind - 1, ind]) = [];
            else
                %Grab value from the view ui-objects
                analysisMethod = obj.AnalysisMethod;
            end
            
            %Run the analysis
            switch analysisMethod
                case 'MSC.Nastran'
                    
                    %Make the method object
                    Nas = awi.methods.Nastran;                    
                    
                    fprintf('Generating FEM...');
                    %Generate the FEM & assign it to the method object
                    FEM = convertToFE(obj.Aircraft);
                    fprintf('Complete!\n');
                    Nas.AnalysisModel = FEM;

                    %Do the analysis
                    TrimResult = Nas.trim(obj.Aircraft, obj.LoadCase);
                                        
                otherwise
                    %Use Robbie Cook's Nonlinear Intrinsic Strain-Based FE
                    ISBNFE = convertAWI2ISBNFE(obj);
                    
                    %Using what beam model?
                    res = obj.BeamModel;
                    assert(~isempty(res), 'no beam model');
                    
                    %Draw something in the view
                    ha = axes('Parent', uicontainer('Parent', hb), 'NextPlot', 'add', 'Box', 'on');
                    scatter3(ha, res.BM.Node.Coord(:,1), res.BM.Node.Coord(:,2),res.BM.Node.Coord(:,3),'k.')
                    axis(ha,'equal')
                    hold(ha,'on')
                    %PlotAircraft(x_dyn1,Matrices,MassProps)
                    %PlotAircraftv2([], ISBNFE.Matrices,0,0,[])
                    
                    %Do the trim analysis
                    [TrimResult, x_trim] = trim_ISBNFE(obj, ISBNFE);
                    
            end           
           
            %Add to session
            obj.BeamModel.add(TrimResult); 
            
            if nargout == 1
                varargout{1} = x_trim;
            end
            
        end
        
    end
    
    methods ( Access = protected )
        
        function val = getLoadCases(obj)
            
            %Start here
            val = getLoadCases@awi.view.Analysis(obj);
            
            %Anything ?
            if ~isempty(val)
                
                %Down-select
                val(~ismember({val.LoadCaseType},{'Cruise Design Case','Manoeuvre'})) = [];
                
            end
            
        end
        
        function update(obj)
            
            %Start with base class
            update@awi.view.Analysis(obj);
            
            %Any beam models available ?
            if isempty(obj.BeamModels)
                str = {};
            else
                str = {obj.BeamModels.Name};
            end
            
            %Get current selection
            sel = get(obj.hBeamModel, 'Value');
            
            %Ensure selection popup remains valid
            sel = min(sel, numel(str));
            
            %In case of no loadcases
            if sel == 0
                set(obj.hBeamModel, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
            else
                set(obj.hBeamModel, 'String', str, 'Value', sel, 'Enable', 'on');
            end
            
            
            %TODO: Display the result of the trim analysis (if it exists)
            
        end
        
        function cbMaxIterations(obj,~,~)
            
            obj.MaxIterations = str2double(get(obj.hMaxIterations, 'String'));
            
        end
         
        function onBeamModelChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
        end
       
    end
    
end
