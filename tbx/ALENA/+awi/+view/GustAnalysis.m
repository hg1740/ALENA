classdef GustAnalysis < awi.view.Analysis
    
    properties
        
        TrimResult
        TrimResults
        RBflag        = 0;
        EndTime       = 1;
        TimeStep      = 0.01;
        MaxIterations = 10; %passed on to getTrim
        BeamModel;
        BeamModels;
        
    end
    
    properties (Transient, SetAccess = protected )
        
        hTrimResult
        hFFflag
        hCflag
        hEndTime
        hTimeStep
        hMaxIterations
        hBeamModel;
        
    end
    
    properties (Dependent)
        %Names of the 'awi.model.ResultSet' objects
        GustResultNames
    end
    
    methods % set / get
        
        function val = get.TrimResult(obj)
            
            %Get from uicontrol
            val = get(obj.hTrimResult, 'Value');
            
            %Valid ?
            val(val < 0 | val > numel(obj.TrimResults)) = [];
            
            %Return the actual selection, rather than the index of the selection
            val = obj.TrimResults(val);
            
        end
        
        function set.TrimResult(obj, val)
            
            %Make a note
            obj.TrimResult = val;
            
            %No selection is valid from user's point of view
            if isempty(val)
                
                %But we need something to put into popup
                idx = 1;
                
            elseif isnumeric(val)
                
                %TODO: ensure within range ?
                idx = val;
                
            elseif isa(val, 'awi.model.TrimResult')
                
                %Allow selection to be set with an actual member of the list
                [b, idx] = ismember(val, obj.TrimResults); %#ok<MCSUP>
                
                %Ensure valid
                assert(b, 'invalid trim result');
                
            end
            
            %Assign in uicontrol
            set(obj.hTrimResult, 'Value', idx); %#ok<MCSUP>
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.TrimResults(obj)
            
            %No model ?
            if isempty(obj.BeamModel)
                
                %Then no LoadCases
                val = [];
                
            else
                
                %What have we got ?
                val = findall(obj.BeamModel, 'Type', 'TrimResult');
                
            end
            
        end
        
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
        
        function val = get.GustResultNames(obj)
            %get.GustResultNames Get method for the dependent property
            %'GustResultNames'.
            %
            % The results objects are named based on the value of gust
            % gradient.
            
            val = [];
            
            if isempty(obj.LoadCase)
                return
            end
            
            val = arrayfun(@(h) sprintf('One-Minus-Cosine (H = %im)', ...
                ceil(h)), obj.LoadCase.GustLength, 'Unif', false);
            
        end
        
    end
    
    methods % constructor
        function obj = GustAnalysis( model, varargin )
            
            %Call superclass constructor - do not pass the model in, better to wait
            obj@awi.view.Analysis( model, ...                'Title', 'Trim Analysis', ...
                'OfferLoadCaseSelection', true, ...
                'PushButtons', {'Analyse...', @analyse}, ...
                varargin{:} );
            
            %Update the 'Configure Analysis' page
            %   - New prompt
            obj.hPrompt.String = [obj.hPrompt.String ; {} ; {['If there ' , ...
                'is no ''One-Minus-Cosine'' load case present in the '    , ...
                'session then you can create one in one of three ways:']} ; ...
                {['1. In the tree view, right-click on ''Load Cases''>'   , ...
                'Edit>Add>Load Case, then change the type to ''One-Minus-Cosine''.']} ; ...
                {['2. In the tree view, right-click on an existing load ', ...
                'case and select Edit>Duplicate>This. This will create ' , ...
                'a new load case object with the same properties (e.g. ' , ...
                'altitude, mach, etc.). Next, change the type to ''One-Minus-Cosine''.']}
                {'3. In the tree view, select an existing load case and change the type to ''One-Minus-Cosine''.'}];
            obj.hPrompt.Parent.Heights = [200 -1];
            
            %Add a label, and popup (to control selection)
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Beam Model:')
            obj.hBeamModel = uicontrol('Parent', obj.hControls, 'Style', 'popup', ...
                'Callback', @obj.onBeamModelChange);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;
            
            %Add a label, and popup (to control selection)
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Trim Result:')
            obj.hTrimResult = uicontrol('Parent', obj.hControls, 'Style', 'popup', ...
                'Callback', @obj.onTrimResultChange);
            
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
            
            % Create a edit box for EndTime
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'end time:')
            obj.hEndTime = uicontrol('Parent', obj.hControls,...
                'Style', 'edit',...
                'String', num2str(obj.EndTime), ...
                'Callback', @obj.cbEndTime);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;
            
            % Create a edit box for TimeStep
            uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'time step:')
            obj.hTimeStep = uicontrol('Parent', obj.hControls,...
                'Style', 'edit',...
                'String', num2str(obj.TimeStep), ...
                'Callback', @obj.cbTimeStep);
            
            %Good height
            obj.hControls.Heights(end-1:end) = 25;
            
            % Create a check box for RBflag
            %             uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'time step:')
            uicontrol('Parent', obj.hControls, ...
                'Style', 'text', ...
                'String', ['Boundary Conditions', char(10)]);
            obj.hFFflag = uicontrol('Parent', obj.hControls, ...
                'Style'        , 'checkbox' , ... 
                'String'       , 'Free-Free', ...
                'Value'        , 1          , ...
                'ToolTipString', 'Run analysis using free-free boundary conditions', ...
                'Callback'     , {@obj.cbBCflag, 'free-free'});
            obj.hCflag = uicontrol('Parent', obj.hControls, ...
                'Style'        , 'checkbox', ...
                'String'       , 'Fixed'   , ...
                'ToolTipString', 'Run analysis using clamped root conditions', ...
                'Callback'     , {@obj.cbBCflag, 'clamped'});
                
            
            %Good height
            obj.hControls.Heights(end-2:end) = 25;
            
            %Now store the model - which will in turn trigger an update
            obj.Model = model;
            
        end        
    end
    
    methods
        
        function analyse(obj, varargin)
            %analyse Runs a gust analysis using the chosen aeroelastic
            %solver and analysis options as dictated by the user.
            
            %Invoke superclass method
            %   - Clears up the view by removing any extant progress
            %     windows etc.
            %AO = anlyse@awi.view.Analysis --> SHOULD BE THIS!
            analyse@awi.view.Analysis(obj);
            
            %Using what method/formulation? - Move this to superclass
            % - Something like this...
            %switch obj.Solver
            %   case Nonlinear Strain-Based
            %       AO = ... (AO = Analysis Object)
            %   case Nonlinear Beam Shapes
            %   case Nonlinear Displacement-Based 
            %
            % - And then something like this...
            %   Results = fn(AO)
            %
            ISBNFE = convertAWI2ISBNFE(obj);
            
            % Make sure we have a gust length
            assert(~isempty(obj.LoadCase.GustLength),'No Gust Length Specified!');
            
            %Using what beam model ?
            res = obj.BeamModel;
            assert(~isempty(res), 'no beam model');         
            
            % Initialise the panel for the view
            hb = uiextras.HBox('Parent', obj.hTabPanel, 'Tag', 'View');
            obj.hTabPanel.TabTitles{end} = 'View';
            obj.hTabPanel.Selection = numel(obj.hTabPanel.TabTitles);
            
            %Draw something in the view
            ha = axes('Parent', uicontainer('Parent', hb), 'NextPlot', 'add', 'Box', 'on');            
            scatter3(ha, res.BM.Node.Coord(:,1), res.BM.Node.Coord(:,2),res.BM.Node.Coord(:,3),'k.')
            axis(ha,'equal')
            hold(ha,'on')
            %PlotAircraft(x_dyn1,Matrices,MassProps)
            %PlotAircraftv2([], ISBNFE.Matrices,0,0,[]) 
            
            %Run the gust analysis
            [GustResult, ~] = gust_ISBNFE(obj, ISBNFE);
            
            %'GustResult' is of class 'awi.model.TransientResultSet'
            %   - We want to return both transient and 'static' results
            %   - For every 'GustResult' create an equivalent static result
            %     that contains the max/min of the time-domain response
            StaticGustResult = cell(1, numel(GustResult));               
            for iG = 1 : numel(GustResult) %For every gust...
                
                %Make the 'awi.model.ResultSet' object
                StaticGustResult{iG} = awi.model.ResultSet('Name', [GustResult(iG).Name, ' (Envelope)']);
                
                %Grab the data from the underlying 'BeamResults'
                br  = GustResult(iG).BeamResults;
                br_ = cell(1, numel(br));
                
                %For each beam result, grab the beam handle and find the
                %max/min loads across all time steps                
                for iB = 1 : numel(br) %For every 'awi.model.Beam'...
                
                    %Make the 'awi.model.BeamResult' (or subclass) object
                    fn = str2func(class(br(iB)));
                    br_{iB} = fn();
                    
                    %Assign the handle to the associated beam
                    br_{iB}.Beam = br(iB).Beam;
                    
                    %Grab the hidden 'BeamProperties' from the gust beam
                    %results obejcts
                    BeamProp = br(iB).BeamProperties;
                    
                    %Make sure we can select the 'Internal Loads'
                    types = {BeamProp.Type};
                    utype = unique(types);
                    if ~any(ismember(utype, 'Internal Loads'))
                        continue
                    end
                    
                    %For now, we just care about the internal loads
                    BP = selectBeamProp(BeamProp, 'Internal Loads');
                    
                    %Calculate the max/min
                    maxmin = arrayfun(@(bp) [max(bp.Value, [], 1) ; min(bp.Value, [], 1)], BP, 'Unif', false);
                    
                    %Loop through beam properties and assign data to the
                    %new instance of 'awi.model.BeamResult'
                    for iBP = 1 : numel(BP) %For every results quantity...                     
                        %What quantity?
                        q = BP(iBP).Quantity;                        
                        %Assign the value
                        br_{iB}.(q) = maxmin{iBP};                        
                        %Assign the distribution
                        br_{iB}.([q, '_eta']) = BP(iBP).Distribution;                        
                    end
                    
                end
                
                %Assign the 'awi.model.BeamResult' objects to the
                %'awi.model.ResultSet' object
                StaticGustResult{iG}.BeamResults = horzcat(br_{:});
                
            end
            
            %Collect the data
            StaticGustResult = horzcat(StaticGustResult{:});
            
            %Add to session
            obj.BeamModel.add(GustResult);
            obj.BeamModel.add(StaticGustResult);
            
        end
        
    end
    
    methods ( Access = protected )
        
        function val = getLoadCases(obj)
            
            val = getLoadCases@awi.view.Analysis(obj);
            val(~ismember({val.LoadCaseType},{'One-Minus-Cosine','Continuous-Turbulence'})) = [];
            
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
            
            %Any trim results available ?
            if isempty(obj.TrimResults)
                str = {};
            else
                str = {obj.TrimResults.Name};
            end
            str{end+1} = 'Calculate';
            
            %Get current selection
            sel = get(obj.hTrimResult, 'Value');
            
            %Ensure selection popup remains valid
            sel = min(sel, numel(str));
            
            %In case of no loadcases
            if sel == 0
                set(obj.hTrimResult, 'String', str, 'Value', 1, 'Enable', 'off');
            else
                set(obj.hTrimResult, 'String', str, 'Value', sel, 'Enable', 'on');
            end
            
            %TODO: Display the result of the trim analysis (if it exists)
            
        end
        
        function cbMaxIterations(obj,~,~)
            
            obj.MaxIterations = str2double(get(obj.hMaxIterations, 'String'));
            
        end
        
        function cbEndTime(obj,~,~)
            
            obj.EndTime = str2double(get(obj.hEndTime, 'String'));
            
        end
        
        function cbTimeStep(obj,~,~)
            
            obj.TimeStep = str2double(get(obj.hTimeStep, 'String'));
            
        end
        
        function cbBCflag(obj,src,~, varargin)
            
            %If obj.RBflag is true then the model is free-free
            %   - TODO : Update the property name to something more
            %   descriptive.
            
            %Toggle the status of the other checkbox
            switch varargin{1}
                case 'clamped'
                    if src.Value
                        set(obj.hFFflag, 'Value', 0);
                    else
                        set(obj.hFFflag, 'Value', 1);
                    end
                case 'free-free'
                    if src.Value
                        set(obj.hCflag, 'Value', 0);
                    else
                        set(obj.hCflag, 'Value', 1);
                    end
            end
            
            %Pass flag to underlying property
            if obj.hFFflag.Value
                %Free-Free conditions
                obj.RBflag = true;
            elseif obj.hCflag.Value
                %Clamped conditions
                obj.RBflag = false;
            end
            
%             obj.RBflag = get(obj.hRBflag, 'Value');
            
        end
        
        function onBeamModelChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
        end
        
        function onTrimResultChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
        end
        
    end
    
end
