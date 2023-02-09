classdef Analysis < mvc.view.Container
    %
    % Generic awi analysis, provides a starting point from which specific analyses can be sub-classed.
    
    properties
        %The model to which we are attached
        Model;
        %Object describing the analysis method we are using 
        AnalysisObject
        %The load case with which we are associated (may or may not be relevant, depending on analyis)
        LoadCase;
        %Function to call on selection change (not used, but required for compatibility with ViewManager)
        SelectionChangeFcn;
        %Track currently selected item(s) (not used, but required for compatibility with ViewManager)
        Selection;
        %Helpful message displayed on creation of the view associated with this Analysis
        Description = 'Configure Analysis using controls on the right hand side of this view';
    end
    
    properties (Dependent) %Useful shortcuts
        %Handle to the 'awi.model.Aircraft' object
        Aircraft;
        %Handle to the 'awi.model.LoadCases' objects in the session
        LoadCases;
        %Currently selected analysis method
        AnalysisMethod
    end
    
    properties %(Access = protected)
        %Does this analysus require a choice between methods?
        OfferMethodSelection = '';
        %Does this analysis require a choice between load cases ?  One of 'none', 'single', or 'multiple'
        OfferLoadCaseSelection = 'single';
        %List of pushbutton labels and functions, to be displayed in view
        Pushbuttons = {};
        %Load case selection (may or may not be relevant, depending on analyis)
        hLoadCase;
        %Method selection (may or may not be relevant, depending on
        %analysis)
        hMethod
        %The panel in which we display the analysis
        hPanel;
        %The tab panel which offers different views of the results.
        hTabPanel
        %Container for the options/configurables
        hOptions
        %A helpful prompt
        hPrompt;
        %Control the analysis (may or may not be relevant, depending on analyis)
        hControls;
        %Placeholder for listeners
        Listeners;
    end
    
    properties (Constant)
       %List of available analysis methods
       AnalysisMethods = {'MSC.Nastran', 'Instrinsic Strain-Based FE', 'Finite Volume Beam FE'};
    end
    
    methods % set / get
        
        function val = get.LoadCase(obj)
            
            %Get from uicontrol
            val = get(obj.hLoadCase, 'Value');
            
            %Valid ?
            val(val < 0 | val > numel(obj.LoadCases)) = [];
            
            %Return the actual selection, rather than the index of the selection
            val = obj.LoadCases(val);
            
        end
        
        function set.LoadCase(obj, val)
            
            %Make a note
            obj.LoadCase = val;
            
            %No selection is valid from user's point of view
            if isempty(val)
                
                %But we need something to put into popup
                idx = 1;
                
            elseif isnumeric(val)
                
                %TODO: ensure within range ?
                idx = val;
                
            elseif isa(val, 'awi.model.LoadCase')
                
                %Allow selection to be set with an actual member of the list
                [b, idx] = ismember(val, obj.LoadCases); %#ok<MCSUP>
                
                %Ensure valid
                assert(b, 'invalid load case');
                
            end
            
            %Assign in uicontrol
            set(obj.hLoadCase, 'Value', idx); %#ok<MCSUP>
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.LoadCases(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no LoadCases
                val = [];
                
            else
                
                %What have we got ?
                val = getLoadCases(obj);
                
            end
            
        end
        
        function val = get.Aircraft(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no Aircraft
                val = [];
                
            else
                
                %Pass it on
                val = obj.Model.Aircraft;
                
            end
            
        end
        
        function val = get.AnalysisMethod(obj)
            %get.AnalysisMethod Retrieves the name of the currently
            %selected analysis method.
           
            val = '';
            
            if isempty(obj.hMethod) %Escape route
                return
            end
            
            val = obj.hMethod.String{obj.hMethod.Value};
            
        end
        
        function set.Model(obj, val)
            
            %Make a note
            obj.Model = val;
            
            %If setting to nothing
            if isempty(val)
                
                %That's all
                return;
                
            end
            
            %Update content
            update(obj);
            
            %Listen for any subsequent change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged); %#ok<MCSUP>
            
        end
        
    end
    
    methods % construction
        
        function obj = Analysis(model, varargin)
            %Analysis Class constructor for the 'awi.view.Analysis' class.            
            
            %% Parse inputs
            
            %Caller may be supplying properties we need to grab at this level
            prp = {'OfferLoadCaseSelection', 'OfferMethodSelection', 'PushButtons', 'AnalysisOptions'};
            
            %Look for match in varargin
            idx = zeros(size(varargin));
            b = cellfun(@ischar, varargin);
            [~, idx(b)] = cellfun(@(x)ismember(x, prp), varargin(b));
            
            %Anything ?
            if any(idx)
                
                %Pull out corresponding values
                val = varargin(find(idx > 0) + 1);
                
                %Eliminate from varargin
                varargin(sort([find(idx > 0), find(idx > 0) + 1])) = [];
                
            end
            
            %So the properties we found in varargin are
            prp = prp(idx(idx > 0));
            
            %Call superclass constructor
            obj@mvc.view.Container(varargin{:});
            
            %Assign values at this level, if any were passed in
            if ~isempty(prp)
                pav = [prp; val];
                set(obj, pav{:});
            end            
            
            %% Build the GUI elements
            
            %Start with an HBox
            hh = uiextras.HBoxFlex('Parent', obj.UIContainer);
            
            %Containing a panel in which the analysis is displayed
            obj.hPanel = uiextras.Panel('Parent', hh, 'Padding', 6);
            
            %Populate the panel with a tab panel ready for some options
            obj.hTabPanel = uiextras.TabPanel('Parent', obj.hPanel, 'TabWidth', 100);
            hb = uiextras.VBox('Parent', obj.hTabPanel, 'Tag', 'Configure Analysis');
            obj.hTabPanel.TabTitles{end} = 'Configure Analysis';
            
            %Create a thin box at the top of the page to provide some help
            %to the user
            obj.hPrompt = uicontrol('Parent', hb, ... %Display a helpful message
                'Style', 'text', ...
                'String', {' ', ' ', ' ', obj.Description});
            
            %Create a container for the analysis options
            obj.hOptions = uix.Grid('Parent', hb);
            hb.Heights = [100, -1];
            
            %And a container for additional controls
            obj.hControls = uiextras.VBox('Parent', hh);
            hh.Widths(end) = 120;
            
            %Does this analysis require the user choose a Load Case ?
            switch obj.OfferLoadCaseSelection
                
                case {false, 'off', 'none'}
                    
                    %Nothing to do
                    
                case {true, 'one', 'single'}
                    
                    %Add a label, and popup (to control selection)
                    uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Load Case:')
                    obj.hLoadCase = uicontrol('Parent', obj.hControls, ...
                        'String'  , '-none-', ...
                        'Style'   , 'popup' , ...                        
                        'Callback', @obj.onLoadCaseChange);
                    
                    %Good height
                    obj.hControls.Heights(end-1:end) = 25;
                    
                case 'multiple'
                    
                    %Add a label, and listbox (to control selection)
                    uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Load Cases:')
                    obj.hLoadCase = uicontrol('Parent', obj.hControls, ...
                        'String', '-none-' , ...
                        'Style' , 'listbox', ...
                        'Min'   , 0        , ...
                        'Max'   , 2        , ...
                        'Callback', @obj.onLoadCaseChange);
                    obj.hControls.Heights(end-1:end) = [25, -1];
                    
                otherwise
                    error('bad option for ''OfferLoadCaseSelection''');
            end
            
            %Does this analysis require the user choose an analysis method?
            switch obj.OfferMethodSelection
                
                case {false, 'off', 'none'}
                    
                    %Nothing to do
                    
                case {true, 'on'}
                    
                    %Add a label, and popup (to control selection)
                    uicontrol('Parent', obj.hControls, 'Style', 'Text', 'String', 'Analysis Method:')
                    obj.hMethod = uicontrol('Parent', obj.hControls, ...
                        'String', '-none-', ...
                        'Style' , 'popup', ...
                        'Callback', @obj.onMethodChange);
                    
                    %Good height
                    obj.hControls.Heights(end-1:end) = 25;
                    
                otherwise
                    
            end
            
            %Add the buttons
            for i = 1:2:numel(obj.Pushbuttons)
                
                %Add the button
                uicontrol('Parent', obj.hControls, ...
                    'Style', 'Pushbutton', ...
                    'String', obj.Pushbuttons{i}, ...
                    'Callback', {@obj.onAnalyse, obj.Pushbuttons{i + 1}});
                
                %Good height
                obj.hControls.Heights(end) = 25;
                
            end
            
            %% Store the model
            obj.Model = model;
            
        end
        
    end
    
    methods % analysis
        
        function analyse(obj)
            %analyse Clears any extant views from the Tab Panel object
            %(obj.hTabPanel)
            
            %Keep children of the tab panel whose tag matches these tokens
            keepviews = {'Configure Analysis'};
            
            %What children does the tab panel have?
            ch = obj.hTabPanel.Children;
            
            %Escape route
            if isempty(ch)
                return
            end
            
            %Get the tags of the children of the tab panel
            tags = {ch.Tag};
            
            %Filter & deleete extant views
            idx  = ismember(tags, keepviews);
            delete(ch(~idx));
            
        end
        
    end
    
    methods ( Access = protected )
        
        function val = getLoadCases(obj)
            
            %Simple search
            val = findall(obj.Model, 'isa', 'awi.model.LoadCase');
            
        end
        
        function update(obj)
            
            %Anything to do ?
            if isempty(obj.Model)
                
                %Bail
                return;
                
            end
            
            %Need to show Load Case selection in view ?
            if obj.OfferLoadCaseSelection
                
                %Yes - list of loadcases (by name)
                if isempty(obj.LoadCases)
                    str = {};
                else
                    str = {obj.LoadCases.Name};
                end
                
                %Get current load case selection
                sel = get(obj.hLoadCase, 'Value');
                
                %Ensure selection popup remains valid
                sel = min(sel, numel(str));
                
                %In case of no loadcases
                if sel == 0
                    set(obj.hLoadCase, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
                else
                    set(obj.hLoadCase, 'String', str, 'Value', sel, 'Enable', 'on');
                end
                
            end
            
            %Need to show method selection in view?
            if obj.OfferMethodSelection
               
                %Yes - list of methods (by name)   
                str = obj.AnalysisMethods;
                
                %Get current selection
                sel = get(obj.hMethod, 'Value');
                
                %Ensure selection popup remains valid
                sel = min(sel, numel(str));
                
                %In case of no loadcases
                if sel == 0
                    set(obj.hMethod, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
                else
                    set(obj.hMethod, 'String', str, 'Value', sel, 'Enable', 'on');
                end
                
            end
            
        end
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onLoadCaseChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
        end
        
        function onMethodChange(obj, ~, ~)
           
            %Update content
            update(obj);
            
        end
        
        function onAnalyse(obj, hc, ~, fcn)
            
            %Might take a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Do it
                fcn(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
            end
            
            %Update content
            update(obj);
            
            %If we get this far and there's nothing being shown
            if isempty(obj.hPanel.Children)
                
                %Reinstate helpful message, if it has gone away
                % (TODO: not an ideal approach, must be a better way)
                if isempty(obj.hPrompt) || ~ishghandle(obj.hPrompt)
                    obj.hPrompt = uicontrol('Parent', obj.hPanel, ...
                        'Style', 'text', ...
                        'String', {'', '', '', obj.Description});
                end
                
            end
            
        end
        
    end
    
    methods %Method block for the intrinsic strain-based nonlinear FE formulation
        
        function ISBNFE = convertAWI2ISBNFE(obj)
            %convertAWI2ISBNFE Converts an AWI model into an analysis model
            %that is compatible with the Intrinsic, Strain-Based,
            %Nonlinear, Finite-Element (ISBNFE) formulation.
            %
            %   - TODO : Eventually this will be moved into a dedicated
            %            object that will handle all of this!            
            %   - TODO : The structure 'ISBNFE' should contain enough
            %            information such that it can be plotted on its own.
            %   - TODO : Should be able to initialise the ISBNFE
            %            formulation from an AEROFLEX object/structure
            
            %Start with a blank structure
            ISBNFE = struct();
            
            %Grab the AeroFlex beam model
            res = obj.BeamModel;
            assert(~isempty(res), 'no beam model');
            
            %Attach the load case fuel masses to the beam model! 
            %   - Create a new partID and append the data to the Conm2s
            BM = FuelObj2BeamModel(obj.LoadCase.FuelDistribution,res.BM);
  
            %Convert the AeroFlex beam model to the method specific
            %analysis model
            [ISBNFE.Matrices, ISBNFE.Aero,~] = InitialiseRobbie(BM);
            
            %Grab the names of all parts in the beam model
            ISBNFE.AllPartIDs = unique({res.BM.PartId(:).Part}, 'stable');
            
        end
        
        function [TrimResult, x_trim] = trim_ISBNFE(obj, ISBNFE, varargin)
            %trim_ISBNFE Performs a trim analysis using the Intrinsic,
            %Strain-Based, Nonlinear, Finite-Element (ISBNFE) formulation
            %and returns the results data (awi.model.ResultSet) and the
            %formulation specific data
            
            if nargin < 2
                ISBNFE = convertAWI2ISBNFE(obj);
            end
            
            p = inputParser;
            addParameter(p, 'SkipResultsGeneration', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %% Run the trim analysis
            x_trim = i_trim_ISBNFE(obj, ISBNFE);
            
            %% Plot results
            % Initialise the panel for the results view
            hb = uiextras.Grid('Parent', obj.hTabPanel);
            obj.hTabPanel.TabTitles{end} = 'Results';
            obj.hTabPanel.Selection = numel(obj.hTabPanel.Children);
            
            for i = [1 4 2 5 3 6]
                ha = axes('Parent', uicontainer('Parent',hb));
                plot(ha,x_trim.x_f{1}(i:6:end))
                xlabel(ha,'Span (m)')
                ylabel(ha,getLabels(i,'f'))
                grid(ha,'on')
                grid(ha,'minor')
                set(ha,'FontSize',14)
            end
            hb.Widths = [-1 -1 -1];
            
            %% Generate the 'awi.model.ResultSet' object 
            
            %Do we need to bother?
            if p.Results.SkipResultsGeneration
                TrimResult = [];
                return
            end
            
            %What do we want to call them?
            setNames = cellfun(@(x) sprintf('Nonlinear Trim (%s)', x), {obj.LoadCase.Name}, 'Unif', false);
            
            %Pass to generic method...
            TrimResult = generateResultsObjects(obj, {x_trim}, ISBNFE, setNames);
            
            
        end
        
        function x_trim = i_trim_ISBNFE(obj, ISBNFE)
            %i_trim_ISBNFE Performs a trim analysis using the Intrinsic,
            %Strain-Based, Nonlinear, Finite-Element (ISBNFE) formulation.
            
            if nargin < 2
                ISBNFE = convertAWI2ISBNFE(obj);
            end
            
            Sim.grav_fact               = obj.LoadCase.LoadFactor;
            [ISBNFE.Aero.rho,~,~,SoS,~] = ISA_h(obj.LoadCase.Altitude);
            ISBNFE.Aero.stall_angle     = 180;
            x_trim0.x_vg                = [-obj.LoadCase.Mach*SoS;0.0;0.0;0.0;0.0;0.0];
            
            %Assign control surface deflections as Trim Variables
            TrimVars = i_assignCSDeflections(obj.LoadCase, ISBNFE);
            
            % Specify the rigid body degree of freedom
            TrimVars.rb_vars = {2,1*pi/180,'Pitch Angle'};
            TrimVars.Thrust  = [];
            
            % Specify what it is we are trying to trim for (Lift = 3,
            % Pitching Moment = 5)
            if isempty(TrimVars.Flaps)
                TrimVars.ForceBalance = 3;
            else
                TrimVars.ForceBalance = [3 5];
            end
            
            % Numerical Damping for the static solver
            TrimVars.Damping      = 0.75;
            Sim.damping           = 1;
            
            % Initialise the panel for the trim view
            hb = uiextras.VBox('Parent', obj.hTabPanel);
            obj.hTabPanel.TabTitles{end} = 'Progress';
            obj.hTabPanel.Selection = numel(obj.hTabPanel.Children);
            
            %Do the TRIM analysis
            x_trim = getTrim( ...
                x_trim0, ...
                ISBNFE.Matrices, ...
                ISBNFE.Aero    , ...
                Sim            , ...
                TrimVars       , ...
                hb             , ...
                obj.MaxIterations);
            
            function TrimVars = i_assignCSDeflections(LoadCase, ISBNFE)
                %i_assginCSDeflections Assigns the control surface
                %deflections from the LoadCase to the analysis model.
                
                TrimVars.Flaps = cell(0,4);
                
                % Search for free control surfaces
                cs_free_idx = strcmpi(LoadCase.CsDeflections.Type,'free');
                ParentNames = LoadCase.CsDeflections.Parent;
                lc_cs_eta   = vertcat(LoadCase.ControlSurfaces.Eta);
                lc_cs_name  = {LoadCase.ControlSurfaces.Name};
                lc_cs_def   = -[LoadCase.CsDeflections.Value];
                
                fixedcount = 0;
                for i = 1:size(LoadCase.CsDeflections,1)
                    
                    member_id  = find(ismember(ISBNFE.AllPartIDs, ParentNames{i}));
                    member_eta = ISBNFE.Matrices.s{member_id}/ISBNFE.Matrices.s{member_id}(end);
                    cs_eta     = lc_cs_eta(i,:);
                    [~,node_id_in]  = min(abs(cs_eta(1) - member_eta));
                    [~,node_id_out] = min(abs(cs_eta(2) - member_eta));
                    beam_ids   = (node_id_in:node_id_out-1);
                    
                    if isempty(regexpi(lc_cs_name{i},'^p'))
                        lc_cs_def(i) = -lc_cs_def(i);
                    end
                    
                    if cs_free_idx(i)
                        fixedcount = fixedcount + 1;
                        TrimVars.Flaps{fixedcount,1} = member_id;
                        TrimVars.Flaps{fixedcount,2} = {beam_ids};
                        TrimVars.Flaps{fixedcount,3} = lc_cs_def(i);
                        TrimVars.Flaps{fixedcount,4} = lc_cs_name{i};
                    else
                        ISBNFE.Aero.delta_flap{member_id}(beam_ids) = pi/180*lc_cs_def(i);
                    end
                end
                
            end
            
        end
        
        function [GustResult, x_out] = gust_ISBNFE(obj, ISBNFE, varargin)
            %gust_ISBNFE Performs a gust analysis using the Intrinsic,
            %Strain-Based, Nonlinear, Finite-Element (ISBNFE) formulation.
            
            if nargin < 2
                ISBNFE = convertAWI2ISBNFE(obj);
            end
            
            p = inputParser;
            addParameter(p, 'SkipResultsGeneration', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %% Do a trim analysis
            
            %Prepare flight parameters
            [ISBNFE.Aero.rho,~,~,SoS,~] = ISA_h(obj.LoadCase.Altitude);
            ISBNFE.Aero.stall_angle = 180;
            x_in.x_vg     = [-obj.LoadCase.Mach*SoS;0.0;0.0;0.0;0.0;0.0];
            Sim.aero_flag = 1;
            
            if isempty(obj.TrimResult) %Do we need to run a trim analysis first?
                
                %Trim the aircraft
                [~, x_trim] = trim_ISBNFE(obj, ISBNFE, 'SkipResultsGeneration', true);

                %Update the CoG location and the control surface deflection
                ISBNFE.Matrices.CGa    = x_trim.CGa;
                ISBNFE.Aero.delta_flap = x_trim.delta_flap;
                
            else
                
                %Grab the data from the stored result
                %   - TODO : This ability is not currently catered for by
                %            the current results set.
                ISBNFE.Matrices.CGa    = obj.TrimResult.RotMat;
                ISBNFE.Aero.delta_flap = obj.TrimResult.ControlDeflection;
                
            end
            
            %% Prep the plots
            
            % Initialise the panel for the results view
            hb = uiextras.Grid('Parent',obj.hTabPanel, 'Tag', 'Results');
            obj.hTabPanel.TabTitles{end} = 'Results';
            obj.hTabPanel.Selection      = numel(obj.hTabPanel.Children);
            
            %Populate the axes
            ha = gobjects(1, 6);
            for i = [1 4 2 5 3 6]
                ha(i) = axes('Parent',uicontainer('Parent',hb));
                hold(ha(i),'on')
                xlabel(ha(i),'Span (m)')
                ylabel(ha(i), getLabels(i,'f'))
                grid(ha(i),'on')
                grid(ha(i),'minor')
                set(ha(i),'FontSize',14)
            end
                     
            %% Run the gust analysis
            
            switch lower(obj.LoadCase.LoadCaseType)
                
                case 'one-minus-cosine'
                    
                    %Set up simulation environment
                    %   - TODO : All this should be moved into the GUI
                    Sim.residual_flag  = 1;
                    Sim.sys_flag       = 1;
                    
                    Sim.Soln           = 2;
                    Sim.time.limits    = [0.00 obj.EndTime];
                    Sim.time.dt        = obj.TimeStep;
                    Sim.rb_flag        = obj.RBflag;
                    Sim.speed_up       = 1;
                    
                    Sim.eps            = 1e-3;
                    Sim.rbloads_flag   = 0;
                    Sim.rho_inf        = 1;
                    
                    ISBNFE.Aero.Leishman.a1   = 0.165;
                    ISBNFE.Aero.Leishman.a2   = 0.335;
                    ISBNFE.Aero.Leishman.b1   = 0.041;
                    ISBNFE.Aero.Leishman.b2   = 0.320;
                    
                    wgmax_vec = GustFamilyGen([],...
                        obj.LoadCase.GustLength,...
                        obj.LoadCase.Altitude,...
                        obj.LoadCase.Mach*SoS,...
                        Sim.time.limits(1):Sim.time.dt:Sim.time.limits(2),0);
                    
                    ISBNFE.Aero.gust.H        = obj.LoadCase.GustLength(1); % Gust gradient (gust length/2)
                    ISBNFE.Aero.gust.xg       = 0;
                    ISBNFE.Aero.gust.yg       = 0;
                    ISBNFE.Aero.gust.zg       = 0;
                    ISBNFE.Aero.gust.wgmax    = wgmax_vec.Gust.Amp(1);
                    ISBNFE.Aero.gust.theta    = 0*pi/180; % Gust orientation
                    ISBNFE.Aero.gust_switch   = 1;
                    ISBNFE.Aero.inflow_switch = 0;
                    ISBNFE.Aero.inflow_vel    = [0;0;0];
                    
                    Sim.grav_fact      = 1;
                    
                    x_dyn     = initDynSim(x_in,ISBNFE.Matrices,Sim);
                    
                    %% Run Dynamic Code
                    
                    %Tell the user how far through the simulation we are
                    hwb = waitbar(0, 'Running gust simulations...please wait, or hit ''Cancel'' to interrupt', ...
                        'CreateCancelBtn', 'delete(gcbf)', ...%@i_cancel
                        'Name', 'AWI Framework');
                    clu = onCleanup(@()delete(hwb(ishghandle(hwb))));
                    
                    %Preallocate
                    %   - One result for every gust...
                    x_out = cell(1, length(obj.LoadCase.GustLength));                    
                    
                    %Loop through gust gradient(s) and run simulation
                    for ii = 1:length(obj.LoadCase.GustLength)
                        
                        if ~ishandle(hwb) %Has the user cancelled?
                            break
                        end
                        
                        %Update the waitbar
                        waitbar(ii/length(obj.LoadCase.GustLength), hwb);
                        
                        %Grab gust gradient and velocity vector field
                        ISBNFE.Aero.gust.H     = obj.LoadCase.GustLength(ii);
                        ISBNFE.Aero.gust.wgmax = wgmax_vec.Gust.Amp(ii);
                        
                        %Run analysis
                        x_out{ii} = getDynamic( ...
                            x_dyn          , ...
                            ISBNFE.Matrices, ...
                            ISBNFE.Aero    , ...
                            Sim);
                        
                        %Plot the envelope of loads
                        for i = [1 4 2 5 3 6]
                            hp = plot(ha(i),max(x_out{ii}.x_f{1}(i:6:end,:),[],2));
                            plot(ha(i),min(x_out{ii}.x_f{1}(i:6:end,:),[],2),'Color',hp.Color)
                        end
                        hb.Widths = [-1 -1 -1];
                        drawnow
                    end
                    
                case 'continuous-turbluence'
                    
                otherwise
                    
                    error('Invalid Load Case Type!')
                    
            end
            
            %% Generate the 'awi.model.ResultSet' objects
            
            %Do we need to bother?
            if p.Results.SkipResultsGeneration
                GustResult = [];
                return
            end
            
            %Gust names - Combination of gust gradient and parent load case   
            gustNames = strcat(obj.GustResultNames, [' (', obj.LoadCase.Name, ')']);
            
            %Pass it to the generic method
            %   - 1 x 'awi.model.ResultSet' object per gust gradient
            GustResult = generateResultsObjects(obj, x_out, ISBNFE, gustNames);
            
        end
        
        function AWIBeamResult = generateResultsObjects(obj, x_data, ISBNFE, setNames)
            %generateResultsObjects Converts the 'x_data' structure that is
            %returned by the Intrinsic Strain-Based, Nonlinear Finite
            %Element (ISBNFE) solver into a series of AWI
            %'awi.model.ResultSet' objects.
            %
            %   - 1 x 'awi.model.ResultSet' object is created per cell in 
            %     the cell array 'x_data'.
            %   - 1 x 'awi.model.AeroelasticResult' is generated per
            %     element in the 'x_data.x_f' cell array.
            %   - 'x_data' is anologous to to the 'x_trim' and 'x_out'
            %     variables in the original 'FALCON' code.

            if nargin < 3 %ISBNFE method data
                ISBNFE = convertAWI2ISBNFE(obj);
            end
            if nargin < 4 %Names of the results sets
                setNames = [];
            end
            
            %Need access to the 'awi.model.Aircraft' object so that we can
            %then search the aircraft for the different beam objects.
            ac = findall(obj.Model, 'isa', 'awi.model.Aircraft');
                        
            %How many results sets?
            nSet = numel(x_data);
            
            %Check/assign results set names
            if isempty(setNames)
                setNames = arrayfun(@(i) sprintf('Set %i', i), 1 : nSet, 'Unif', false);
            else
                assert(numel(setNames) == nSet, ['The number of '     , ...
                    'results sets and the number of provided results ', ...
                    'set names must be the same.'])
            end
            
            %Take a sneak-peak to see if the results are static or dynamic 
            fname = {'x_f', 'x_disp', 'Cl', 'Cll'}; %fields we can interrogate
            idx   = cellfun(@(x) isfield(x_data{1}, x), fname);
            ind   = find(idx == true, 1);   %what fields are present?
            
            %Error out
            if isempty(ind)
                error('Unable to extract any results from the analysis. Check the output.');
            end
            
            %Results are transient if the output data has more than one
            %column
            if size(x_data{1}.(fname{ind}){1}, 2) > 1
                fn = @awi.model.TransientResultSet;
            else
                fn = @awi.model.ResultSet;
            end

            %Preallocate
            AWIBeamResult = cell(1, nSet);
            
            for iS = 1 : nSet %For each results set...
                
                %Create a single instance of 'awi.model.ResultSet'
                AWIBeamResult{iS} = fn('Name', setNames{iS});
                
                %Do not assign results to for the 'VTP' as it is currently only
                %modelled with one element and is therefore imcompatible with
                %the 'awi.mixin.Beamable' class.
                %   - TODO : Finesse this!
                partidx = find(~strcmp(ISBNFE.AllPartIDs,'VTP'));
                
                for i = partidx %For each part in the model...
                    
                    %Find the AWI LiftingSurface object based on the name in
                    %the AeroFlex beam model.
                    %   - TODO : Should store the handle instead!
                    LS_Results = awi.model.AeroelasticResult;
                    LS = findall(ac, 'Name', ISBNFE.AllPartIDs{i});
                    set(LS_Results, 'Beam', LS);
                    
                    %Calculate the eta position for the element
                    eta_mid    = ISBNFE.Matrices.s_out{i}(ISBNFE.Matrices.s_mp_ind{i})/ISBNFE.Matrices.s_out{i}(end);
                    eta_nodes  = ISBNFE.Matrices.s_out{i}/ISBNFE.Matrices.s_out{i}(end);
                                        
                    if isfield(x_data{iS}, 'x_f') %Internal Loads
                        
                        % Assign INTERNAL LOADS in the analysis local frame
                        LS_Results.Fx     = x_data{iS}.x_f{i}(1:6:end, :)';
                        LS_Results.Fx_eta = eta_mid;
                        LS_Results.Fy     = x_data{iS}.x_f{i}(2:6:end, :)';
                        LS_Results.Fy_eta = eta_mid;
                        LS_Results.Fz     = x_data{iS}.x_f{i}(3:6:end, :)';
                        LS_Results.Fz_eta = eta_mid;
                        LS_Results.Mx     = x_data{iS}.x_f{i}(4:6:end, :)';
                        LS_Results.Mx_eta = eta_mid;
                        LS_Results.My     = x_data{iS}.x_f{i}(5:6:end, :)';
                        LS_Results.My_eta = eta_mid;
                        LS_Results.Mz     = x_data{iS}.x_f{i}(6:6:end, :)';
                        LS_Results.Mz_eta = eta_mid;
                        %
                        %Populate the FAME loads
                        % RMatrix = [0, -1, 0 ; 1, 0, 0 ; 0, 0, 1];
                        %   (Beam local frame to FAME local frame)
                        
                        % Assign the INTERNAL LOADS in the FAME Coordinate frame
                        LS_Results.FxFAME     = x_data{iS}.x_f{i}(3:6:end, :)'; %Swap is due to FAME wing axis convention
                        LS_Results.FxFAME_eta = eta_mid;
                        LS_Results.FyFAME     = x_data{iS}.x_f{i}(1:6:end, :)'; %Swap is due to FAME wing axis convention
                        LS_Results.FyFAME_eta = eta_mid;
                        LS_Results.FzFAME     = x_data{iS}.x_f{i}(2:6:end, :)';
                        LS_Results.FzFAME_eta = eta_mid;
                        LS_Results.MxFAME     = x_data{iS}.x_f{i}(6:6:end, :)'; %Swap is due to FAME wing axis convention
                        LS_Results.MxFAME_eta = eta_mid;
                        LS_Results.MyFAME     = x_data{iS}.x_f{i}(4:6:end, :)'; %Swap is due to FAME wing axis convention
                        LS_Results.MyFAME_eta = eta_mid;
                        LS_Results.MzFAME     = x_data{iS}.x_f{i}(5:6:end, :)';
                        LS_Results.MzFAME_eta = eta_mid;
                        
                    end
                                        
                    if isfield(x_data{iS}, 'Cl') && isfield(x_data{iS}, 'Cll') %Aero Loads
                        
                        % Assign AERODYNAMIC LOADS
                        LS_Results.Cl         = x_data{iS}.Cl{i}(:,3);
                        LS_Results.Cl_eta     = eta_mid;
                        LS_Results.Cl_l       = x_data{iS}.Cll{i}(:,3);
                        LS_Results.Cl_l_eta   = eta_mid;
                        
                    end
                    
                    if isfield(x_data{iS}, 'x_disp') %Displacements
                        
                        % Assign DISPLACEMENTS
                        LS_Results.Tx         = squeeze(x_data{iS}.x_disp{i}(1,1,:));
                        LS_Results.Tx_eta     = eta_nodes;
                        LS_Results.Ty         = squeeze(x_data{iS}.x_disp{i}(2,1,:));
                        LS_Results.Ty_eta     = eta_nodes;
                        LS_Results.Tz         = squeeze(x_data{iS}.x_disp{i}(3,1,:));
                        LS_Results.Tz_eta     = eta_nodes;
                        
                    end
                    
                    %Assign data to the results set
                    if isempty(AWIBeamResult{iS}.BeamResults)
                        AWIBeamResult{iS}.BeamResults = LS_Results;
                    else
                        AWIBeamResult{iS}.BeamResults(end + 1) = LS_Results;
                    end
                    
                end
                
            end
            
            %Return a object array
            AWIBeamResult = horzcat(AWIBeamResult{:});
            
        end
        
    end
    
    methods %Method block for the displacement-based nonlinear FE formulation
        function convertAWI2DBNFE(obj)
            %convertAWI2DBNFE Converts an AWI model into an analysis model
            %that is compatible with the Displacement-Based, Nonlinear,
            %Finite-Element (DBNFE) formulation.
            
        end
    end
    
    methods %Method block for the nonlinear beam shapes formulation
        
        function convertAWI2NBS(obj)
            %convertAWI2NBS Converts an AWI model into an analysis model
            %that is compatible with the Nonlinear Beam Shapes (NBS)
            %formulation.
            
        end
        
    end
    
    
end
