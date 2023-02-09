classdef BeamResultsViewer < mvc.view.Container
    %BeamResultsViewer Generic results viewer for results quantities
    %defined over a beam.
    %
    % TODO - Allow the exporting to be extendable...
    
    properties
        %The model to which we are attached
        Model;
        %Function to call on selection change (not used, but required for compatibility with ViewManager)
        SelectionChangeFcn;
        %Track currently selected item(s) (not used, but required for compatibility with ViewManager)
        Selection;
        %Number of axes to display
        NumAxes = 1;
        %Numbers of columns into which results are arranged
        Cols = 1;
        %Lists of the extra variables being plotted on the left y-axis
        VariablesYLeft  = {};
        %Lists of the extra variables being plotted on the right y-axis
        VariablesYRight = {};
        %List of the results sets whose data is being plotted on the right
        %y-axis
        ResultsLeft  = {};
        %List of the results sets whose data is being plotted on the left
        %y-axis
        ResultsRight = {};
        %What arrangements of results are supported by this viewer?
        SupportedArrangements;
        %Which arrangement is currently selected?
        CurrentArrangement;
        %What styles are supported ?
        SupportedStyles = {'line', 'scatter', 'table'};
        %What style is currently selected?
        CurrentStyle;
        %What layouts are supported ?
        SupportedLayouts = {'grid', 'tabs'};
        %What layout is currently selected?
        CurrentLayout;
        %What tab is currently selected (if any) ?
        CurrentTab = 1;
    end
    
    properties (Dependent)
        %Numbers of rows into which results are arranged
        Rows
    end
    
    %Controlling export process
    properties
        %Path to folder for exporting
        ExportPath
        %Array of 'awi.model.ResultsSet' objects
        Res2Export
        %Cell array of beam names that are used to filter the exported
        %results
        BeamNames2Export
        %Cell array of results types that are used to filter the exported
        %results
        ResTypes2Export
        %Logical flag indicating whether the results should be exported in
        %the AWI file format
        ExportAWI
        %Logical flag indicating whether the results should be exported in
        %the CSV file format
        ExportCSV
        %Logical flag indicating whether the results should be exported in
        %the Excel file format
        ExportExcel
    end
    
    %Shortcuts
    properties (Dependent)
        %Handle to all 'awi.model.ResultSet' objects
        Results
        %Handle to current results
        CurrentResults
        %Current list of beam results
        CurrentBeamResults
        %Name of all available result sets
        AvailableResultsNames
        %Name of the current results set
        CurrentResultsNames
        %Handle to all 'awi.model.Beam' objects within the current
        %selection
        AvailableBeams
        %Handle to the currently selected 'awi.model.Beam' object
        CurrentBeam
        %Name of the 'awi.model.Beam' objects within the current selection
        AvailableBeamNames
        %Name of the current 'awi.model.Beam' object
        CurrentBeamName
        %Name of all beam property 'Types' in the current selection
        AvailableResultsTypes
        %Name of the current results type
        CurrentResultsType
        %Name of all the results quantities in the current selection
        AvailableResultsQuantities
        %Name of the current results quantities being plotted
        CurrentResultsQuantities
        
        Arrangement;
        CurrentLayoutIndex;
        
    end
    
    %Helper property to point towards the underlying BeamProperty objects
    properties (Dependent, Hidden = true)
        %2d array of all hidden BeamProperty objects from the resulting
        %selection
        CurrentBeamProperties
    end
    
    %To help persist view arrangements
    properties (Access = protected)
        ArrangementFileMask = {'*.raf', 'Results viewer arrangement files (*.raf)'; ...
            '*.mat', 'MATLAB data files (*.mat)'; ...
            '*.*', 'All files (*.*)'};
        ArrangementFile;
    end
    
    %Handles to GUI elements
    properties (Access = protected)
        %List of pushbutton labels and functions, to be displayed in view
        Pushbuttons = {};
        %Results selection
        hResults;
        %Beam results selection
        hBeamResults
        %List of results types that can be plotted
        hResultsType
        %List of results quantities that can be plotted
        hResultsList
        %Push button for calculating the loads envelope
        hLoadsEnv
        %The panel in which we display the results
        hPanel;
        %Axes embedded within panel
        hAxes;
        %Pushbutton(s) to control the analysis (may or may not be relevant, depending on analyis)
        hButtons;
        %Placeholder for listeners
        Listeners;
        %Controls displaying numbers of rows / columns, layout
        hNumAxes;
        hCols;
        hLayout;
        %Somewhere to attach a single legend - avoiding duplication and
        %freeing up real-estate)
        hLegend;
    end
    
    %Selection values for the various listboxes/popups
    properties (Access = protected)
        %Index number of the current beam results selection
        BeamSelection
        %Index number of the current results type selection
        PlotTypeSelection
        %Index number of the current results quantity selection
        PlotQuantitySelection
    end
    
    methods % set / get
        
        function set.NumAxes(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnan', 'finite', 'real', 'positive'}, class(obj), 'NumAxes');
            obj.NumAxes = val;
        end
        
        function set.Cols(obj, val)
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnan', 'finite', 'real', 'positive'}, class(obj), 'Cols');
            obj.Cols = val;
        end
        
        function val = get.Rows(obj)            
            val = ceil(obj.NumAxes / obj.Cols);
        end
        
        %TODO - CONSIDER WHETHER WE NEED TO VALIDATE THE SELECTION VALUES
        
        function set.Selection(obj, val)
            
            %             validateattributes(val, {'numeric'}, {'nonnegative'}, ...
            %                 class(obj), 'Selection');
            %
            %             %Selection stores as index into list of available results
            %             if isnumeric(val)
            %
            %                 %Ensure valid
            %                 assert(isempty(val) || all(ismember(val, 1:numel(obj.Results))), 'invalid selection'); %#ok<MCSUP>
            %
            %             elseif isa(val, 'awi.model.ResultSet')
            %
            %                 %Ensure valid
            %                 [b,idx] = ismember(val, obj.Results); %#ok<MCSUP>
            %                 assert(all(b), 'invalid selection');
            %                 val = idx(b);
            %                 val = val(:).';
            %
            %             else
            %
            %                 %This type of view can't present whatever is being selected.
            %                 % COULD issue a warning and continue
            %                 % warning(['selection set to class ''', class(val), '''']);
            %                 %
            %                 % but better I think to just continue, with no change, silently
            %                 return;
            %
            %             end
            
            %Make a note
            obj.Selection = val;
            
        end
        
        function set.BeamSelection(obj, val)
            
            validateattributes(val, {'numeric'}, {'nonnegative'}, ...
                class(obj), 'Selection');
            
            %             %Selection stores as index into list of available results
            %             if isnumeric(val)
            %
            %                 %Ensure valid
            %                 assert(isempty(val) || all(ismember(val, 1:numel(obj.AvailableBeamNames))), 'invalid selection'); %#ok<MCSUP>
            %
            %             elseif isa(val, 'awi.model.ResultSet')
            %
            %                 %Ensure valid
            %                 [b,idx] = ismember(val, obj.Results); %#ok<MCSUP>
            %                 assert(all(b), 'invalid selection');
            %                 val = idx(b);
            %                 val = val(:).';
            %
            %             else
            %
            %                 %This type of view can't present whatever is being selected.
            %                 % COULD issue a warning and continue
            %                 % warning(['selection set to class ''', class(val), '''']);
            %                 %
            %                 % but better I think to just continue, with no change, silently
            %                 return;
            %
            %             end
            
            %Make a note
            obj.BeamSelection = val;
            
        end
        
        function set.PlotTypeSelection(obj, val)
            
            validateattributes(val, {'numeric'}, {'nonnegative'}, ...
                class(obj), 'Selection');
            
            %             %Selection stores as index into list of available results
            %             if isnumeric(val)
            %
            %                 %Ensure valid
            %                 assert(isempty(val) || all(ismember(val, 1:numel(obj.AvailableResultsTypes))), 'invalid selection'); %#ok<MCSUP>
            %
            %             elseif isa(val, 'awi.model.ResultSet')
            %                 error('We should never get here. Debug code.');
            %                 %Ensure valid
            %                 [b,idx] = ismember(val, obj.Results); %#ok<MCSUP>
            %                 assert(all(b), 'invalid selection');
            %                 val = idx(b);
            %                 val = val(:).';
            %
            %             else
            %
            %                 %This type of view can't present whatever is being selected.
            %                 % COULD issue a warning and continue
            %                 % warning(['selection set to class ''', class(val), '''']);
            %                 %
            %                 % but better I think to just continue, with no change, silently
            %                 return;
            %
            %             end
            
            %Make a note
            obj.PlotTypeSelection = val;
            
        end
        
        function set.PlotQuantitySelection(obj, val)
            
            validateattributes(val, {'numeric'}, {'nonnegative'}, ...
                class(obj), 'Selection');
            
            %             %Selection stores as index into list of available results
            %             if isnumeric(val)
            %
            %                 %Ensure valid
            %                 assert(isempty(val) || all(ismember(val, 1:numel(obj.AvailableResultsQuantities))), 'invalid selection');  %#ok<MCSUP>
            %
            %             elseif isa(val, 'awi.model.ResultSet')
            %                 error('We should never get here. Debug code.');
            %                 %Ensure valid
            %                 [b,idx] = ismember(val, obj.Results); %#ok<MCSUP>
            %                 assert(all(b), 'invalid selection');
            %                 val = idx(b);
            %                 val = val(:).';
            %
            %             else
            %
            %                 %This type of view can't present whatever is being selected.
            %                 % COULD issue a warning and continue
            %                 % warning(['selection set to class ''', class(val), '''']);
            %                 %
            %                 % but better I think to just continue, with no change, silently
            %                 return;
            %
            %             end
            
            %Make a note
            obj.PlotQuantitySelection = val;
            
        end
        
        function val = get.CurrentStyle(obj)
            
            %MUST return a cell-array with at least this many values
            nval = obj.Rows * obj.Cols;
            
            %Start here
            val = obj.CurrentStyle;
            
            %Nothing yet ?
            if isempty(val)
                
                %Go with this by default
                val = repmat(obj.SupportedStyles(1), 1, nval);
                
            elseif numel(val) < nval
                
                %Pad with defaults
                val(end+1:nval) = repmat(obj.SupportedStyles(1), 1, nval - numel(val));
                
            end
            
        end
        
        function val = get.CurrentLayoutIndex(obj)
            
            val = find(strcmp(obj.CurrentLayout, obj.SupportedLayouts));
            
        end
        
        function val = get.CurrentLayout(obj)
            
            %Start here
            val = obj.CurrentLayout;
            
            %Nothing yet ?
            if isempty(val)
                
                %Go with this by default
                val = obj.SupportedLayouts{1};
                
            end
            
        end
        
        function val = get.Arrangement(obj)
            
            %Name
            val.Name = obj.CurrentArrangement;
            
            %How many rows and columns ?
            val.Rows = obj.Rows;
            val.Cols = obj.Cols;
            
        end
        
        function set.Arrangement(obj, val)
            
            %If specified as numeric
            if isnumeric(val)
                
                %Treat as index into list of supported arrangements
                val = obj.SupportedArrangements(val);
                
            end
            
            %Make a note of name
            obj.CurrentArrangement = val.Name;
            
            %Set details
            obj.NumAxes = val.Rows;
            obj.Cols = val.Cols;
            %obj.VariablesX = val.VariablesX;
            %obj.VariablesY = val.VariablesY;
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.VariablesYLeft(obj)
            val = obj.VariablesYLeft;
            if numel(val) < numel(obj.hAxes)
                val = [val, cell(1, numel(obj.hAxes) - numel(val))];
            end
        end
        
        function val = get.VariablesYRight(obj)
            val = obj.VariablesYRight;
            if numel(val) < numel(obj.hAxes)
                val = [val, cell(1, numel(obj.hAxes) - numel(val))];
            end
        end
        
        function val = get.ResultsLeft(obj)
            val = obj.ResultsLeft;
            if numel(val) < numel(obj.hAxes)
                val = [val, cell(1, numel(obj.hAxes) - numel(val))];
            end
            %Check for empties - Force empty matrix, not empty objects
            idx = cellfun(@(x) isempty(x), val);
            val(idx) = {[]};
        end
        
        function val = get.ResultsRight(obj)
            val = obj.ResultsRight;
            if numel(val) < numel(obj.hAxes)
                val = [val, cell(1, numel(obj.hAxes) - numel(val))];
            end
            %Check for empties - Force empty matrix, not empty objects
            idx = cellfun(@(x) isempty(x), val);
            val(idx) = {[]};
        end
        
        function val = get.Results(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no Results
                val = [];
                
            else
                
                %Pass it on
                val = findall(obj.Model,'isa','awi.model.ResultSet');
                
                %Remove any transient result objects
                idx = arrayfun(@(x) isa(x, 'awi.model.TransientResultSet'), val);
                val = val(~idx);
                
            end
            
        end
        
        function val = get.CurrentResults(obj)
            
            val = obj.Results(obj.Selection);
            
        end
        
        function val = get.AvailableResultsNames(obj)
            val = {obj.Results.Name};
        end
        
        function val = get.CurrentResultsNames(obj)
            val = {obj.CurrentResults.Name};
        end
        
        function val = get.CurrentBeamResults(obj)
            
            %Each results set can have multiple beam results (if it has
            %results for more than one beam!)
            allBeamResults = {obj.CurrentResults.BeamResults};
            
            %Simply index the list
            val = cellfun(@(x) x(obj.BeamSelection), allBeamResults, 'Unif', false);
            
            %Return an object array
            val = horzcat(val{:});
        end
        
        function val = get.AvailableBeams(obj)
            
            sel   = obj.Selection;
            beams = arrayfun(@(br) br.Beam, [obj.Results(sel).BeamResults], 'Unif', false);
            val   = horzcat(beams{:});
            
        end
        
        function val = get.CurrentBeam(obj)
            %Simply index the list
            val = obj.AvailableBeams(obj.BeamSelection);
        end
        
        function val = get.AvailableBeamNames(obj)
            
            %Start with the handles
            beams = obj.AvailableBeams;
            
            %If there is not selection then we may not return anything. In
            %which case an empty cell array is fine.
            if isempty(beams)
                val = {''};
                return
            end
            
            %Get the unique beams and their names
            %   TODO - Should be able to invoke 'unique' on the object
            %   array 'beams' and it should JUST WORK. It is failing in the
            %   Collectable class method 'sort'. Need Phil to fix this.
            names = {beams.Name};
            val   = unique(names, 'stable');
            
        end
        
        function val = get.CurrentBeamName(obj)
            val = '';
            if isempty(obj.CurrentBeam)
                return
            end
            val = obj.CurrentBeam.Name;
        end
        
        function val = get.AvailableResultsTypes(obj)
            
            %Start with an empty cell-array
            val = {''};
            
            %Grab all the beam results in the available selection
            allBeamResults = obj.CurrentBeamResults;
            
            %If nothing is selected this will return an empty so just quit
            %out early.
            if isempty(allBeamResults)
                return
            end
            
            %Grab the name of all BeamProperty types and the unique values
            allTypes = {allBeamResults.BeamPropertyTypes};
            uniqueTypes = unique(horzcat(allTypes{:}), 'stable');
            
            %Find the common types across all beam results and return this
            nTypes = cellfun(@(x) numel(intersect(uniqueTypes, x)), allTypes);
            
            %Select the cell with the lowest number of interestions
            [~, ind] = min(nTypes);
            val  = allTypes{ind};
            
        end
        
        function val = get.CurrentResultsType(obj)
            val = obj.AvailableResultsTypes(obj.PlotTypeSelection);
        end
        
        function val = get.AvailableResultsQuantities(obj)
            
            val = {};
            
            %All beam results for current selection
            BeamResult = obj.CurrentBeamResults;
            
            if isempty(BeamResult)
                return
            end
            
            %All results types belonging to this
            allResultTypes = obj.AvailableResultsTypes;
            
            %Particular results type that we are interested in
            resultType = allResultTypes{obj.PlotTypeSelection};
            
            %Grab the underlying 'awi.mixin.BeamProperty' objects
            BP = getBeamProperty(BeamResult(1), resultType);
            
            %Filter out any data that is just NaN
            idx = ~cellfun(@(x) any(any(isnan(x))), {BP.Value});
            
            %Return all the names of the results quantities in a cell array
            val = {BP(idx).Name};
            
        end
        
        function val = get.CurrentResultsQuantities(obj)
            val = obj.AvailableResultsQuantities(obj.PlotQuantitySelection);
        end
        
        function val = get.CurrentBeamProperties(obj)
            if isempty(obj.CurrentBeamResults)
                val = [];
                return
            end
            
            %Stack ALL of the beam properties
            bp  = vertcat(obj.CurrentBeamResults.BeamProperties);
            
            %Down-select to 'awi.model.BeamProperty' objects of type
            %'obj.CurrentResultsType' && 'obj.CurrentResultsQuantities'
            idx = and( ...
                ismember({bp(1, :).Type}, obj.CurrentResultsType{1}), ...
                ismember({bp(1, :).Name}, obj.CurrentResultsQuantities)); %Make sure we filter by 'Name' and not 'Quantity' as the list view (Results Quantities) is populated using the Names!
            val = bp(:,idx);
            
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
        
        function set.ExportPath(obj, val)
            
            validateattributes(val, {'char'}, {'scalartext'}, class(obj), ...
                'ExportPath');
            
            %Check the path is fully qualified
            [p, ~, ~] = fileparts(val);
            if isempty(p)
                val = fullfile(pwd, val);
            end
            
            obj.ExportPath = val;
            
        end
        
        function set.Res2Export(obj, val)
            
            if isempty(val)
                obj.Res2Export = val;
            else
                validateattributes(val, {'awi.model.ResultSet'}, {'vector'}, ...
                    class(obj), 'Res2Export');
                obj.Res2Export = val;
            end
            
        end
        
        function set.BeamNames2Export(obj, val)
            
            if ~isempty(val)
                assert(iscellstr(val), ['Expected ''BeamNames2Export'' ', ...
                    'to be a cell array of strings']);
            end
            
            %Assign
            obj.BeamNames2Export = val;
            
        end
        
        function set.ResTypes2Export(obj, val)
            
            if ~isempty(val)
                assert(iscellstr(val), ['Expected ''ResTypes2Export'' ', ...
                    'to be a cell array of strings']);
            end
            
            %Assign
            obj.ResTypes2Export = val;
            
        end
        
        function set.ExportAWI(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'ExportAWI');
            obj.ExportAWI = val;
        end
        
        function set.ExportCSV(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'ExportCSV');
            obj.ExportCSV = val;
        end
        
        function set.ExportExcel(obj, val)
            validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
                'ExportExcel');
            obj.ExportExcel = val;
        end
        
    end
    
    methods % construction
        
        function obj = BeamResultsViewer( model, varargin )
            
            %Caller may be supplying properties we need to grab at this level
            prp = {'PushButtons'};
            
            %Look for match in varargin
            idx = zeros(size(varargin));
            b = cellfun(@ischar, varargin);
            [~, idx(b)] = cellfun(@(x)ismember(x, prp), varargin(b));
            
            %Anything ?
            if any(idx)
                
                %Pull out corresponding values
                val = varargin(find(idx > 0) + 1);
                
                %Elliminate from varargin
                varargin(sort([find(idx > 0), find(idx > 0) + 1])) = [];
                
            end
            
            %So the properties we found in varargin are
            prp = prp(idx(idx > 0));
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Assign values at this level, if any were passed in
            if ~isempty(prp)
                pav = [prp; val];
                set(obj, pav{:});
            end
            
            %Start with an HBox
            hh = uiextras.HBoxFlex('Parent', obj.UIContainer, 'Spacing', 6, 'Padding', 6);
            
            %Add a panel in which the results are displayed
            obj.hPanel = uiextras.Panel('Parent', hh, 'Padding', 6);
            
            %Add a tab panel for the tool strip and the legend
            ht = uix.TabPanel('Parent', hh);
            
            %Add a scrolling panel to house the tool-strip
            hs = uix.ScrollingPanel('Parent', ht);
            
            %Add a vertical strip into which controls are embedded
            hv = uiextras.VBox('Parent', hs, 'Padding', 6, 'Spacing', 6);
            
            %Good width
            hh.Widths(end) = 200;
            
            %Good height for simple text labels ?
            hLab = 20;
            
            %A panel in which the results are selected
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Analysis Results');
            hv.Heights(end) = hLab;
            obj.hResults = uicontrol('style', 'listbox', ...
                'Parent',hv,...
                'min'   ,0, ...
                'max'   ,2, ...
                'Callback',@obj.onSelectionChange, ...
                'String'  , '-Select a Results Set-');
            hv.Heights(end) = 100;
            
            %A panel in which the beam results are selected
            uicontrol('Parent', hv, 'Style', 'text', 'String', 'Beam Results');
            hv.Heights(end) = hLab;
            obj.hBeamResults = uicontrol('style', 'listbox', ...
                'Parent'  , hv, ...
                'min'     , 0 , ...
                'max'     , 1 , ...
                'Callback', @obj.onSelectionChange, ...
                'String'  , '-Select a Results Set-');
            hv.Heights(end) = 75;
            
            %Add a selection for results quantity types
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Results Types');
            hv.Heights(end) = hLab;
            obj.hResultsType = uicontrol('Parent', hv, ...
                'Style'   , 'popup', ...
                'Value'   , 1, ...
                'String'  , '-Select a Results Set-', ...
                'Callback', @obj.onSelectionChange);
            hv.Heights(end) = 25;
            
            %A panel for down-selecting results quantities
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Results Quantities');
            hv.Heights(end) = hLab;
            obj.hResultsList = uicontrol('Parent', hv, ...
                'Style'   , 'listbox', ...
                'min'     , 0, ...
                'max'     , 2, ...
                'Callback', @obj.onSelectionChange, ...
                'String'  , '-Select a Results Set-');
            hv.Heights(end) = 100;
            
            %Add a popup for layout
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Layout');
            hv.Heights(end) = hLab;
            obj.hLayout = uicontrol('Parent', hv,...
                'Style'   , 'popup',...
                'Value'   , obj.CurrentLayoutIndex, ...
                'String'  , obj.SupportedLayouts, ...
                'Callback', @obj.cbLayout);
            hv.Heights(end) = 25;
            
            %Edit box for Rows
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'No. of axes');
            hv.Heights(end) = hLab;
            obj.hNumAxes = uicontrol('Parent', hv,...
                'Style'   , 'edit',...
                'String'  , num2str(obj.Rows), ...
                'Callback', @obj.cbNumAxes);
            hv.Heights(end) = 25;
            
            %Edit box for Columns
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Cols');
            hv.Heights(end) = hLab;
            obj.hCols = uicontrol('Parent', hv,...
                'Style'   , 'edit',...
                'String'  , num2str(obj.Cols), ...
                'Callback', @obj.cbCols);
            hv.Heights(end) = 25;
            
            %Add a pushbutton for exporting data
            uicontrol('Parent', hv, ...
                'Style'   , 'pushbutton'    , ...
                'String'  , 'Export Results...', ...
                'Callback', @obj.cbExportResults);
            hv.Heights(end) = 25;
            
            %Add a push button for calculating the envelope
            obj.hLoadsEnv = uicontrol('Parent', hv, ...
                'Style'   , 'pushbutton' , ...
                'String'  , 'Show Loads Envelope', ...
                'Callback', @obj.cbShowLoadsEnvelope, ...
                'Enable'  , 'off');
            hv.Heights(end) = 25;
            
            %And a container, to be used later to attach legend
            obj.hLegend = uicontainer('Parent', ht);
            
            %Update the tab panel names and selection
            ht.TabTitles = {'Select Results', 'Legend'};
            ht.Selection = 1;
            ht.TabWidth  = hh.Widths(end) * 0.45;
            
            %Can we find a default view arrangements file ?
            arrangement(obj, 'default');
            
            %Update the current height of the uix.VBox inside the
            %uix.ScrollingPanel --> Will ensure that the scroll bar appears
            hs.Heights        = -1;
            hs.MinimumHeights = sum(hv.Heights) + ...       %Account for heights of children
                ((numel(hv.Heights)-1) * hv.Padding) + ...  %Account for padding between each child
                (2 * hv.Spacing);                           %Account for spacing around all children
            
            %Store the model
            obj.Model = model;
            
        end
        
    end
    
    methods (Access = protected) % updating and displaying the data
        
        function update(obj)
            %update Updates the view.
            %
            % Performs the following actions:
            %   - Clears all axes
            %   - Refreshes the content in the right-hand panel based on
            %     the current selection in the list-boxes.
            %   - Creates the new axes based on the desired layout
            %   - Populates the axes with data based on the current
            %     selection and any additional data added via the context
            %     menu.
            %   - Adds a legend on the second tab panel of the right-hand
            %     panel.
            
            %Clean sheet
            delete(obj.hAxes);
            
            %Invoke get method once!
            Res = obj.Results;
            
            %Anything to do ?
            if isempty(Res)
                
                %Ensure listbox empty
                set(obj.hResults, 'String', {});
                
                %Ensure panel has no content
                delete(obj.hPanel.Children);
                
                %We're done
                return;
                
            end
            
            %Grab everything we need to populate the GUI with
            Data = getEverythingOnce(obj, Res);
            
            %Beam properties of type 'obj.CurrentResultsType' &&
            %'obj.CurrentResultsQuantities'
            BeamProperty = Data.BeamProperty;
            
            %% Update the list-boxes
            
            %List available results and show selection
            set(obj.hResults, ...
                'String', Data.ResultNames, ...
                'Value' , Data.ResultSelection);
            
            %List available 'awi.model.Beam' objects within the selected
            %results.
            set(obj.hBeamResults, ...
                'String', Data.AvailableBeamNames, ...
                'Value' , Data.BeamSelection);
            
            %List available results sets belonging to the current selected
            %beam result
            set(obj.hResultsType, ...
                'String', Data.ResultTypes, ...
                'Value' , Data.ResultTypeSelection);
            
            %List all the available results quantities that can be plotted
            set(obj.hResultsList, ...
                'String', Data.AllResultQuantities, ...
                'Value' , Data.ResultQuantitySelection);
            
            %% Update the edit-box and "Layout" drop-down menu
            
            %Layout and dimensions of grid
            set(obj.hLayout , 'String', obj.SupportedLayouts, 'Value', obj.CurrentLayoutIndex);
            set(obj.hNumAxes, 'String', num2str(obj.NumAxes));
            set(obj.hCols   , 'String', num2str(obj.Cols));
            
            %% Create the axes layout
            
            nAx = obj.NumAxes;
            
            %Results are viewed in this layout
            switch obj.CurrentLayout
                
                case 'grid'
                    
                    %Simple grid
                    hg = uiextras.Grid('Parent', obj.hPanel);
                    
                case 'tabs'
                    
                    %Make a note of current selection before overwriting
                    nsel = obj.CurrentTab;
                    
                    %Create the tab panel inside a scrolling panel
                    hScroll = uix.ScrollingPanel('Parent', obj.hPanel);
                    
                    %Separate tab for each result
                    hg = uiextras.TabPanel('Parent', hScroll, ...
                        'TabWidth', 100, ...
                        'SelectionChangedFcn', @obj.cbTabSelectionChange);
                    
                otherwise
                    error('bad layout');
                    
            end
            
            %Create one axes for each result - Store the axes in a
            %uicontainer to avoid the GLT bug with axes
            ha = arrayfun(@(x)axes('Parent', uicontainer('Parent', hg), ...
                'UserData', x, 'Box', 'on'), 1 : nAx, 'Unif', false);
            
            %Store in object (maybe needed later, e.g. dynamic update of visible state based on legend click)
            obj.hAxes = [ha{:}];
            
            %Results are viewed in this layout
            switch obj.CurrentLayout
                
                case 'grid'
                    
                    %Arranged accordingly
                    hg.Widths = -ones(1, obj.Cols);
                    
                case 'tabs'
                    
                    %Construct the tab title from a combination of the
                    %primary and additional beam properties...
                    titles = i_getPlotTitles(nAx, BeamProperty, ...
                        obj.VariablesYLeft, obj.VariablesYRight, 'Name');
                    
                    %Set the tab titles
                    hg.TabTitles = titles;
                    
                    %Currently selected tab ?
                    if nsel >= 1 && nsel <= numel(titles)
                        hg.Selection   = nsel;
                        obj.CurrentTab = nsel;
                    end
                    
                    %Update minimum widths of the scrolling panel
                    %   - Buffer of 50 pixels to account for extra spaces
                    %     between tabs
                    hScroll.MinimumWidths = (hg.TabWidth .* numel(hg.TabTitles)) + 50;
                    
                otherwise
            end
            
            %% Populate the axes with data
            
            nDataSets = size(BeamProperty, 2);
            
            %For each element in grid
            for iAx = 1 : nAx
                
                %Primary data from the selection
                if isempty(BeamProperty) || iAx > nDataSets
                    bpData = [];
                else
                    bpData = BeamProperty(:, iAx);
                end
                
                %User helper function to populate grid
                updateElement(obj, iAx, bpData, Data.CurrentResult);
                
            end
            
            %Use helper function to maintain legend
            %   - This is the 'global' legend. It just keeps track of the
            %     colour of the lines associated with the currently
            %     selected results.
            %   - If additional data is added to the axes using the context
            %     menu then a new ('local') legend is made in that axes.
            updateLegend(obj);
            
            function title = i_getPlotTitles(nAx, BP, BP_L, BP_R, tok)
                %i_getPlotTitles Construct the titles of the plot based on
                %the 'Name' of the beam properties on the primary (BP),
                %left (BP_L) and right (BP_R) axes. Can provide a different
                %property name in 'tok' such as 'Quantity' if the shorter
                %beam property name is desired.
                
                if nargin < 4
                    tok = 'Name';
                end
                
                %Primary
                if isempty(BP)           %Grab the beam property string
                    pBPNames = {};
                else
                    pBPNames = {BP(1, :).(tok)};
                end
                if numel(pBPNames) < nAx %Ensure correct dimensions
                    pBPNames = [pBPNames, cell(1, nAx - numel(pBPNames))];
                elseif numel(pBPNames) > nAx
                    pBPNames = pBPNames(1 : nAx);
                end
                
                %Left
                lBPNames = cell(1, nAx);
                idxL = ~cellfun(@isempty, BP_L);
                lBPNames(idxL) = cellfun(@(x) strjoin(unique({x.(tok)}), ', '), BP_L(idxL), 'Unif', false);
                
                %Right
                rBPNames = cell(1, nAx);
                idxR = ~cellfun(@isempty, BP_R);
                rBPNames(idxR) = cellfun(@(x) strjoin(unique({x.(tok)}), ', '), BP_R(idxR), 'Unif', false);
                
                %Combine them
                BPNames = [pBPNames ; lBPNames ; rBPNames];
                
                %Join all names
                title = cell(1, nAx);
                for i = 1 : nAx
                    idx = ~cellfun(@isempty, BPNames(:, i));
                    title{i} = strjoin(BPNames(idx, i), ', ');
                end
                
                %Replace any empties with '-Unassigned-'
                title(cellfun(@isempty, title)) = {'-Unassigned-'};
                
            end
            
        end
        
        function updateElement(obj, iAx, bpData, CurrentResult)
            
            %Grab Beam Properties and current results data from selection
            if nargin < 3 || nargin < 4
                Data = getEverythingOnce(obj);
                if nargin < 3
                    %Primary data from the selection
                    if isempty(Data.BeamProperty) || iAx > size(Data.BeamProperty, 2)
                        bpData = [];
                    else
                        bpData = Data.BeamProperty(:, iAx);
                    end
                end
                if nargin < 4
                    CurrentResult = Data.CurrentResult;
                end
            end
            
            %Grab additional Beam Properties & Results Sets from the
            %context menu (L)
            bpLeft  = obj.VariablesYLeft{iAx};
            resLeft = obj.ResultsLeft{iAx};
            
            %Grab additional Beam Properties & Results Sets from the
            %context  menu (R)
            bpRight  = obj.VariablesYRight{iAx};
            resRight = obj.ResultsRight{iAx};
            
            %Placeholders for lines, legends and messages
            hl  = {};
            
            %Force empty numeric, not empty 'awi.mixin.BeamProperty', and
            %index the additional beam properties to respect the current
            %selection.
            if isempty(bpData)
                bpData = [];
            end
            if isempty(bpLeft)
                bpLeft = [];
            else
                %Down-select so we only show data for the currently
                %selected results set
                idx     = ismember({resLeft.FullName}, {CurrentResult.FullName});
                bpLeft = bpLeft(idx);
                %Force an empty numeric if we have no data
                if nnz(idx) ==0, bpLeft = [];end
            end
            if isempty(bpRight)
                bpRight = [];
            else
                %Down-select so we only show data for the currently
                %selected results set
                idx     = ismember({resRight.FullName}, {CurrentResult.FullName});
                bpRight = bpRight(idx);
                %Force an empty numeric if we have no data
                if nnz(idx) ==0, bpRight = [];end
            end
            
            %Get rid of any children of the current axes
            delete(obj.hAxes(iAx).Children);
            if numel(obj.hAxes(iAx).YAxis) > 1
                yyaxis(obj.hAxes(iAx), 'right');
                delete(obj.hAxes(iAx).Children);
                yyaxis(obj.hAxes(iAx), 'left');
            end
            
            %Careful
            try
                
                %Context menu applied to all axes
                %   - Pass the index of the axes as we want the contextmenu
                %   to be specific to the axes
                hcm = uicontextmenu( ...
                    'Parent', ancestor(obj.hAxes(iAx),'Figure'), ...
                    'CallBack', {@obj.cbVariables, iAx});
                uimenu(hcm, 'Label', 'Add/Remove (y-axis left)');
                uimenu(hcm, 'Label', 'Add/Remove (y-axis right)');
                uimenu(hcm, 'Label', 'Clear additional data', ...
                    'Separator', 'on');
                set(obj.hAxes(iAx), 'UiContextMenu', hcm);
                
                %Define line colours, line style and markers
                nRes = numel(obj.Selection);    % TODO - What is nRes is empty
                %if nRes == 0, nRes = 1; end
                nDat = [size(bpData, 1) / nRes, size(bpLeft, 1) / nRes, size(bpRight, 1) / nRes];    %Number of data quantities to be plotted
                nDat(isnan(nDat)) = 0; %Replace NaN
                nDat(isinf(nDat)) = 0; %Replace Inf
                clr  = uniqueColours(nRes, {'w','k'});  %RGB colours which are distinguishable against a white or black background
                [lin, mkr] = allLineStyles(sum(nDat));
                
                %Create ColorSpec & LineSpec for data on the L/R axes
                nL = (nDat(1) + nDat(2));
                nR = nDat(3);
                
                %Colour
                clrLeft   = num2cell(repmat(clr, [nL, 1]), 2);
                clrRight  = num2cell(repmat(clr, [nR, 1]), 2);
                
                %Line Style
                linStyleL = lin(1 : nL)';
                linStyleR = lin(nL + 1 : end)';
                linStyleL = repmat(linStyleL, [nRes, 1]);  %Repeat by the number of results
                linStyleR = repmat(linStyleR, [nRes, 1]);
                linStyleL = linStyleL(:);                  %Vectorise
                linStyleR = linStyleR(:);
                
                %Marker
                mkrStyleL = mkr(1 : nL)';
                mkrStyleR = mkr(nL + 1 : end)';
                mkrStyleL = repmat(mkrStyleL, [nRes, 1]);  %Repeat by the number of results
                mkrStyleR = repmat(mkrStyleR, [nRes, 1]);
                mkrStyleL = mkrStyleL(:);                  %Vectorise
                mkrStyleR = mkrStyleR(:);
                
                %Use the BeamProperty 'Name' as the 'DisplayName'
                if ~isempty(bpData)
                    names = {bpData.Name}';
                else
                    names = {};
                end
                if ~isempty(bpLeft)
                    namLeft = {bpLeft.Name}';
                else
                    namLeft = {};
                end
                if ~isempty(bpRight)
                    namRight = {bpRight.Name}';
                else
                    namRight = {};
                end
                
                %Combine any extra data belonging to the left-hand y-axis
                %with the current beam property list from the selection
                bpLeft  = [bpData ; bpLeft(:)];
                namLeft = [names ; namLeft];
                
                %Check if we actually have any data to plot!
                %   - If not, ensure the axes is reset
                if isempty(bpLeft)
                    dataLeft = false;
                    %Clear the axes
                    yyaxis(obj.hAxes(iAx), 'left');
                    cla(obj.hAxes(iAx), 'reset');
                else
                    dataLeft = true;
                end
                if isempty(bpRight)
                    dataRight = false;
                    %Clear the axes
                    yyaxis(obj.hAxes(iAx), 'right');
                    cla(obj.hAxes(iAx), 'reset');
                else
                    dataRight = true;
                end
                
                %Bail out if there is no data to plot on left or right axes
                if ~dataLeft && ~dataRight
                    %Remind the user to select some data
                    xl = xlim(obj.hAxes(iAx));
                    yl = ylim(obj.hAxes(iAx));
                    text(obj.hAxes(iAx), mean(xl), mean(yl), ...
                        {'Select data using the GUI element in the ', ...
                        'toolstrip or right-click on the axes to see', ...
                        'a selection of options.'}, ...
                        'VerticalAlignment', 'middle', ...
                        'HorizontalAlignment', 'center', ...
                        'Color', [1, 0, 0]);
                    return
                end
                
                %Hold, and apply labels
                %   - Labels on the left and right y-axis are a
                %   concatenation of all quantities being displayed on
                %   those axes.
                hold(obj.hAxes(iAx), 'on');
                xlabel(obj.hAxes(iAx), 'Eta [-]');
                if dataLeft
                    labelstr = strjoin(unique({bpLeft.Label}, 'stable'), ', ');
                    ylabel(obj.hAxes(iAx), labelstr);
                end
                if dataRight
                    yyaxis(obj.hAxes(iAx), 'right')
                    labelstr = strjoin(unique({bpRight.Label}, 'stable'), ', ');
                    ylabel(obj.hAxes(iAx), labelstr);
                    yyaxis(obj.hAxes(iAx), 'left');
                end
                
                %Interpolate the data from each results set to a common eta
                %distribution.
                if dataLeft  %Data on yy-axis left
                    [x, y] = obj.stackBeamPropertyData(bpLeft);
                end
                if dataRight %Data on yy-axis right
                    [xR, yR] = obj.stackBeamPropertyData(bpRight);
                end
                
                %Present the data
                switch obj.CurrentStyle{iAx}
                    
                    case 'line'
                        
                        %Simple
                        if dataLeft  %Plot data on left yyaxis
                            
                            %Plot the data - TODO : Use 'line' instead of
                            %'plot'
                            hl{end+1} = plot(obj.hAxes(iAx), x, y);
                            
                            %Force the y-axis to be black
                            obj.hAxes(iAx).YColor = 'k';
                            
                            %Update line properties
                            set(hl{end}, ...
                                {'Color'}          , clrLeft  , ...
                                {'LineStyle'}      , linStyleL, ...
                                {'Marker'}         , mkrStyleL, ...
                                'MarkerEdgeColor'  , 'k'      , ...
                                {'MarkerFaceColor'}, clrLeft  , ...
                                {'DisplayName'}    , namLeft);
                            
                        end
                        if dataRight %Plot data on left yyaxis
                            
                            %Switch axes
                            yyaxis(obj.hAxes(iAx), 'right');
                            
                            %Plot the data - TODO : Use 'line' instead of
                            %'plot'
                            hl{end + 1} = plot(obj.hAxes(iAx), xR, yR);
                            
                            %Force the second y-axis to be black
                            obj.hAxes(iAx).YColor = 'k';
                            
                            %Switch back to default y-axis
                            yyaxis(obj.hAxes(iAx), 'left');
                            
                            %Update line properties
                            set(hl{end}, ...
                                {'Color'}          , clrRight , ...
                                {'LineStyle'}      , linStyleR, ...
                                {'Marker'}         , mkrStyleR, ...
                                'MarkerEdgeColor'  , 'k'      , ...
                                {'MarkerFaceColor'}, clrRight , ...
                                {'DisplayName'}    , namRight);
                            
                        end
                        
                    case 'scatter'
                        error('Update ''scatter'' plot for left and right y-axes');
                        %                         %Some care required
                        %                         if size(x,2) == 1
                        %                             htemp = arrayfun(@(k)scatter(obj.hAxes(iAx), x, y(:,k)), 1:size(y,2),'UniformOutput',false);
                        %                             hl{end+1} = vertcat(htemp{:});
                        %                         elseif size(x,2) == size(y,2)
                        %                             htemp = arrayfun(@(k)scatter(obj.hAxes(iAx), x(:,k), y(:,k)), 1:size(y,2),'UniformOutput',false);
                        %                             hl{end+1} = vertcat(htemp{:});
                        %                         else
                        %                             error('size mis-match');
                        %                         end
                        
                    otherwise
                        
                        %TODO
                        hl{end+1} = [];
                        
                end
                
                %Gather line data
                hl = vertcat(hl{:});
                
                %Extend legend string
                resNames = obj.CurrentResultsNames;
                
                %Decide whether we need a legend for this axes
                if numel(hl) == numel(names)
                    
                    %No additional data in the GUI, just name the lines
                    %based on the results sets.
                    set(hl, {'DisplayName'}, resNames');
                    
                else
                    
                    %Additional data on the GUI - Add a legend to this axis
                    n = numel(hl) / numel(resNames);
                    assert(mod(numel(hl), numel(resNames)) ==0, ...
                        ['Expected the number of lines to be an ', ...
                        'integer multiple of the number of results sets.']);
                    
                    %Combine the results name with the beam property names
                    prefix = strcat(repmat(resNames', [n, 1]), {' : '});
                    nam    = strcat(prefix, get(hl, 'DisplayName'));
                    
                    %Just deal with what we can at the moment...
                    set(hl, {'DisplayName'}, nam);
                    
                    %Make the legend - Let MATLAB choose...
                    legend(obj.hAxes(iAx), hl, nam  , ...
                        'Location'   , 'EastOutside', ...
                        'ItemHitFcn' , @obj.legendItemHit);
                    
                end
                
                %Force axis formatting
                obj.hAxes(iAx).Box = 'on';
                
                %Update the global legend
                updateLegend(obj);
                
            catch err
                
                %What went wrong ?
                xl = xlim(obj.hAxes(iAx));
                yl = ylim(obj.hAxes(iAx));
                text(obj.hAxes(iAx), mean(xl), mean(yl), ...
                    {'View creation FAILED with error:', err.message}, ...
                    'VerticalAlignment', 'middle', ...
                    'HorizontalAlignment', 'center');
                
            end
            
        end
        
        function updateLegend(obj)
            %Attempt to display a single legend somewhere more convenient
            %than embedded in each chart
            
            %Ensure clean sheet in container to which legend is attached
            delete(obj.hLegend.Children);
            
            %Find the common lines between all axes based on 'DisplayName'
            dispNames = cell(1, numel(obj.hAxes));
            hl        = cell(size(dispNames));
            for iAx = 1 : numel(obj.hAxes)
                
                ch = cell(1, 2);
                
                %Grab all children on the left and right axes
                if numel(obj.hAxes(iAx).YAxis) > 1
                    yyaxis(obj.hAxes(iAx), 'left');
                    ch{1} = obj.hAxes(iAx).Children;
                    yyaxis(obj.hAxes(iAx), 'right');
                    ch{2} = obj.hAxes(iAx).Children;
                else
                    ch{1} = obj.hAxes(iAx).Children;
                end
                
                %Combine data from left & right axes
                ch_   = vertcat(ch{:});
                
                %Remove any text - Text will be present if the view has
                %errored
                ch_(arrayfun(@(x) isa(x, 'matlab.graphics.primitive.Text'), ch_)) = [];
                if isempty(ch)
                    continue
                end
                
                %Store the handles for later...
                hl{iAx} = flipud(ch_);
                
                %Get display names and store them
                if isa(hl{iAx}, 'matlab.graphics.GraphicsPlaceholder')
                    continue
                end
                dn = get(hl{iAx}, 'DisplayName');
                if isempty(dn)
                    dn = '';
                end
                if ~iscell(dn)
                    dn = {dn};
                end
                dispNames{iAx} = dn;
            end
            
            %Find common data across all between the display names
            uniqueNames = unique(vertcat(dispNames{:}), 'stable'); 
            uniqueNames(cellfun(@isempty, uniqueNames)) = []; %Get rid of empties
            if isempty(uniqueNames ) %Escape route
                return
            end            
            nTypes      = cellfun(@(x) numel(intersect(uniqueNames, x)), dispNames);
            val         = min(nTypes(nTypes ~= 0)); %Select the cell with the lowest number of interestions
            ind         = find(nTypes == val, 1);
            
            if isempty(ind) %Escape route
                return
            end
            
            %Grab the display names and graphics handles
            names = dispNames{ind};
            hline = hl{ind}(ismember(dispNames{ind}, uniqueNames));
            
            %Create an axes inside the legend parent panel
            ax = axes('Parent', obj.hLegend, ...
                'Units'   , 'normalized', ...
                'Position', [0, 0, 1, 1], ...
                'Visible' , 'off');
            
            %Copy the line data onto the new axes
            hl_copy = copyobj(hline, ax);
            
            %Get rid of the line (x, y) data
            arrayfun(@(x) set(x, 'XData', nan, 'YData', nan), hl_copy);
            
            %Add legend to this axes, with custom hit callback
            legend(ax, hl_copy, names, ...
                'Interpreter', 'none', ...
                'ItemHitFcn' , @obj.legendItemHit);
            
        end
        
        function legendItemHit(obj, hc, evtdata)
            
            %Which line are we ?
            nam = evtdata.Peer.DisplayName;
            %disp([nam, ' hit...']);
            
            %Find line in the "dummy" plot used to create legend
            hl = findobj(evtdata.Peer.Parent, 'DisplayName', nam);
            
            if numel(hl) >  1
                hl = hl(1);
            end
            
            %Get visible state
            bVis = strcmp({hl.Visible}, 'on');
            
            %Toggle it
            vis = obj.bool2offon(~bVis);
            set(hl, 'Visible', vis);
            
            %Find corresponding line in all charts
            hl = findobj(obj.hAxes, 'DisplayName', nam);
            
            %Apply new visible state
            set(hl, 'Visible', vis);
            
        end
        
        function onModelChanged(obj, hc, ~ )
            
            %Might take a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Update content
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
            end
            
        end
        
        function onSelectionChange(obj, hc, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note of new selection (for every list/popup box)
                obj.Selection             = get(obj.hResults    , 'Value');
                obj.BeamSelection         = get(obj.hBeamResults, 'Value');
                obj.PlotTypeSelection     = get(obj.hResultsType, 'Value');
                obj.PlotQuantitySelection = get(obj.hResultsList, 'Value');
                
                %Decide whether to enable the loads envelope button
                switch obj.hResultsType.String{obj.PlotTypeSelection}
                    case 'Internal Loads'
                        obj.hLoadsEnv.Enable = 'on';
                    otherwise
                        obj.hLoadsEnv.Enable = 'off';
                end
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
    end
    
    methods (Access = protected) % callbacks
        
        function cbNumAxes(obj, hc, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note
                obj.NumAxes = str2num(get(obj.hNumAxes,'String')); %#ok<ST2NM>
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbCols(obj, hc, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note
                obj.Cols = str2num(get(obj.hCols,'String')); %#ok<ST2NM>
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbLayout(obj, hc, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note
                obj.CurrentLayout = obj.hLayout.String{obj.hLayout.Value};
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbVariables(obj, hc, ~, gdx)
            %cbVariables Builds the context menuo for
            %adding/removing/exporting data to & from the axes
            
            %Ensure the context menu is visible!
            hc.Visible = 'on';
            
            %Delete current child menus below the "Add"/"Remove" level
            %   - Cannot index into hm (e.g. hm(1)) in case findobj return
            %   an empty GraphicsPlaceholder
            hm = findobj(hc, 'Label', 'Add/Remove (y-axis left)');
            hm = [hm, findobj(hc, 'Label', 'Add/Remove (y-axis right)')];
            ch = get(hm, 'Children'); %Shouldn't be necessary!
            if ~isempty(ch)
                delete(vertcat(ch{:}));
            end
            
            %Get information about the results, selection, etc.
            Data = getEverythingOnce(obj);
            
            %Grab pointer to the results set
            Res = Data.CurrentResult;
            
            %Grab header names for the uimenu(s)
            allResultTypes   = obj.AvailableResultsTypes;
            
            %Index into all the beam property objects that have been
            %returned by 'currentBeam' and grabs the beam properties for
            %each results type.
            BeamProps = arrayfun(@(b) ...
                cellfun(@(x) getBeamProperty(b, x), Data.ResultTypes, 'Unif', false), ...
                Data.CurrentBeamResult, 'Unif', false);
            BeamProps = vertcat(BeamProps{:});
            
            %If 'BeamProps' is empty then we can't create the uimenu
            %   - The only alternative is to provide the user with an
            %   enhanced selection which would allow them to choose from
            %   EVERY results set for EVERY beam in those results sets.
            %   This is too cumbersome so just quit out, perhaps inform the
            %   user?
            if isempty(BeamProps)
                %Suppress the context menu
                set(hc, 'Visible', 'off');
                %Provide some guidance to the user
                uiwait(helpdlg(['Before adding data using the context ', ...
                    'menu you must first make a selection in the ', ...
                    'toolstrip in the right-hand panel.'], obj.Title));
                return
            end
            
            %Define the menu identifier for "Left" or "Right" axis
            menuLabel = {'VariablesYLeft', 'VariablesYRight'};
            
            %Define the tag for allocating results data to the "Left" or
            %"Right" axis
            resLabel  = {'ResultsLeft', 'ResultsRight'};
            
            %Grab the names of all the beam properties
            BeamPropNames = cellfun(@(x) {x.Name}, BeamProps(1, :), 'Unif', false);
            
            %Make a note of the current results quantity
            if numel(Data.CurrentResultQuantity) >= gdx
                currentResultQ = Data.CurrentResultQuantity{gdx};
            else
                currentResultQ = 'none';
            end
            
            %Populate the 'Add/Remove' uimenu(s)
            for i = 1 : numel(hm) % Loop through 'Left'/'Right' axes
                for ii = 1 : numel(allResultTypes) %Loop through results types
                    
                    %Grab the current beam property names
                    names = {BeamProps{1, ii}.Name};
                    
                    %First level of indexing to filter by type...
                    bp = vertcat(BeamProps{:, ii});
                    
                    %Determine whether any of the children of this menu
                    %are already plotted
                    %idxP = any(ismember(bp, currentBeamProps));
                    %idxL = any(ismember(bp, obj.(menuLabel{i}){gdx}));
                    %idx  = or(idxP, idxL);
                    idx = any(ismember(bp, obj.(menuLabel{i}){gdx}));
                    tf  = any(idx);
                    
                    %Create the first 'uimenu' for down-selecting by
                    %results type.
                    hM = uimenu(hm(i), 'Label', allResultTypes{ii}, ...
                        'Checked' , mvc.mixin.UiTools.bool2offon(tf));
                    
                    for j = 1 : numel(BeamPropNames{ii})
                        
                        %Determine whether the results are present in the
                        %current axes
                        if size(obj.(menuLabel{i}),1) < gdx
                            tf = false;
                        else
                            tf = any(ismember(bp(:, j), obj.(menuLabel{i}){gdx}));
                        end
                        
                        if tf
                            hM.Checked = mvc.mixin.UiTools.bool2offon(tf(1));
                        end
                        
                        %Determine whether we are on the primary selection
                        enabled = ~strcmp(currentResultQ, BeamPropNames{ii}{j});
                        
                        %If we are on the primary selection then we don't
                        %want to add another one so we must force 'Checked'
                        if ~enabled
                            tf = true;
                        end
                        
                        %Index to find the beam property of interest
                        idx = ismember(names, BeamPropNames{ii}{j});
                        
                        %Check for NaN terms!
                        idx_nan = cellfun(@(x) any(any(~isnan(x))), {bp(:, idx).Value});
                        
                        %Create the second 'uimenu' for down-selecting by
                        %results quantity BUT only do this if this results
                        %quantity has data in one of the results sets
                        if any(idx_nan)
                            uimenu(hM, ...
                                'Label'   , BeamPropNames{ii}{j}, ...
                                'Callback', {@obj.cbVariable, gdx, menuLabel{i}, bp(:, idx), ~tf, Res, resLabel{i}}, ...
                                'Checked' , mvc.mixin.UiTools.bool2offon(tf(1)), ...
                                'Enable'  ,  mvc.mixin.UiTools.bool2offon(enabled));
                        end
                        
                    end
                end
            end
            
            %Assign a callback to the 'Remove additional data' menu
            hm = findobj(hc, 'Label', 'Clear additional data');
            hm.Callback = {@obj.cbRemoveAxisData, gdx};
            
        end
        
        function cbVariable(obj, hc, ~, gdx, menuFlag, BeamProp, addIfTrue, Res, resFlag)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Check if this BeamProp is already on the list
                if size(obj.(menuFlag),1) < gdx
                    b = false;
                else
                    b = ismember(obj.(menuFlag){gdx, 1}, BeamProp);
                end
                
                %Store a reference to the data and a handle to the ResultsSet
                %that the data refers to
                if addIfTrue && ~any(b)
                    
                    %Store all beam properties that we want to display on the
                    %axes
                    if size(obj.(menuFlag),1) < gdx || isempty(obj.(menuFlag){gdx, 1})
                        obj.(menuFlag){gdx} = BeamProp(:);
                        obj.(resFlag){gdx}  = Res;
                    else
                        obj.(menuFlag){gdx}(end + 1 : end + numel(BeamProp)) = BeamProp(:);
                        obj.(resFlag){gdx}(end + 1 : end + numel(Res))       = Res;
                    end
                    
                elseif ~addIfTrue && any(b)
                    
                    %Remove from list
                    obj.(menuFlag){gdx}(b) = [];
                    obj.(resFlag){gdx}(b)  = [];
                    
                else
                    
                    %Ideally this should never happen
                    warning('unexpected inconsistency between caller arguments and current plot configuration');
                    
                end
                
                %Trigger an update of the GUI
                updateElement(obj, gdx);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbArrangements(obj, hc, ~)
            
            %Find the 'Arrangement' menu
            hm = findobj(hc, 'Label', 'Arrangement');
            
            %Clean sheet
            delete(hm.Children);
            
            %Any arrangements available ?
            if isempty(obj.SupportedArrangements)
                
                %No
                uimenu(hm, 'Label', 'none available', 'enable', 'off');
                
            else
                
                %Get name of arrangements
                nam = {obj.SupportedArrangements.Name};
                
                %With separator ?
                sep = cellfun(@(x)x(1) == '|', nam);
                nam(sep) = cellfun(@(x)x(2:end), nam(sep), 'UniformOutput', false);
                sep = arrayfun(@mvc.mixin.UiTools.bool2offon, sep, 'UniformOutput', false);
                
                %Add option to select each, with checkbox to indicate which one is currently selected
                cellfun(@(x,y)uimenu(hm, 'Label', x, ...
                    'Separator', y, ...
                    'Check', mvc.mixin.UiTools.bool2offon(strcmp(x, obj.CurrentArrangement)), ...
                    'Callback', {@obj.cbArrangement, 'apply', x}), ...
                    nam, sep);
                
            end
            
            %Additional options to manage arrangements
            hm = uimenu(hm, 'Label', 'Manage', 'Separator', 'on');
            if ~isempty(obj.CurrentArrangement)
                uimenu(hm, 'Label', ['Save arrangement ''', obj.CurrentArrangement, ''''], 'Callback', {@obj.cbArrangement, 'save'});
            end
            uimenu(hm, 'Label', 'New arrangement...', 'Callback', {@obj.cbArrangement, 'new'});
            uimenu(hm, 'Label', 'Remove arrangement(s)...', 'Callback', {@obj.cbArrangement, 'remove'});
            uimenu(hm, 'Separator', 'on', 'Label', 'Export arrangements...', 'Callback', {@obj.cbArrangement, 'export'});
            uimenu(hm, 'Label', 'Import arrangements...', 'Callback', {@obj.cbArrangement, 'import'});
            
        end
        
        function cbArrangement(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it on
                arrangement(obj, varargin{:});
                
                %Force an update
                update(obj)
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbStyles(obj, hc, gdx)
            
            %Find the 'Style' menu
            hm = findobj(hc, 'Label', 'Style');
            
            %Clean sheet (TODO: is this necessary ?)
            delete(hm.Children);
            
            %Add submenus
            hmm = cellfun(@(x)uimenu(hm, 'Label', x, 'Callback', {@obj.cbStyle, gdx, x}), ...
                obj.SupportedStyles, 'Unif', false);
            hmm = [hmm{:}];
            
            %Selected style ?
            sdx = strcmp(obj.CurrentStyle(gdx), obj.SupportedStyles);
            [hmm(sdx).Checked] = deal('on');
            
        end
        
        function cbStyle(obj, hc, ~, gdx, val)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note
                obj.CurrentStyle{gdx} = val;
                
                %Force an update
                update(obj)
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function varargout = arrangement(obj, action, name)
            
            %TODO - UPDATE THE DEFAULT VIEW ARRANGMENT FILE (.vaf)
            
            %Action not specified
            if nargin < 2 || isempty(action)
                
                %Good default choice
                action = 'apply';
                
            end
            
            %Name provided ?
            if nargin < 3
                
                %No
                name = [];
                
            end
            
            %Do what ?
            switch action
                
                case {'save', 'new'}
                    
                    %If 'save'
                    if strcmp(action, 'save')
                        
                        %Get name of current arrangenent
                        name = obj.CurrentArrangement;
                        
                    end
                    
                    %Get current arrangement
                    A = obj.Arrangement;
                    
                    %Name provided ?
                    if isempty(name)
                        
                        %Ask the user
                        val = inputdlg({'Name of this arrangement:'}, obj.Title, 1, {char(A.Name)});
                        
                        %Unpick
                        if isempty(val) || isempty(val{1})
                            return;
                        else
                            name = val{1};
                        end
                        
                    end
                    
                    %Assign
                    A.Name = name;
                    
                    %Add to store
                    if isempty(obj.SupportedArrangements)
                        obj.SupportedArrangements = A;
                    else
                        
                        %Check for existing name
                        [b, idx] = ismember(A.Name, {obj.SupportedArrangements.Name});
                        
                        %If 'save'
                        if strcmp(action, 'save')
                            
                            %Overwrite existing, or add new as appropriate
                            if ~b
                                idx = numel(obj.SupportedArrangements) + 1;
                            end
                            
                        else
                            
                            %Already exists ?
                            if b
                                
                                %Seek confirmation to overwrite
                                if ~isequal('Yes', questdlg('Overwrite existing arrangement ?', obj.Title, 'Yes', 'No', 'Yes'))
                                    
                                    %Bail out
                                    return;
                                    
                                end
                                
                            else
                                
                                %Next
                                idx = numel(obj.SupportedArrangements) + 1;
                                
                            end
                            
                        end
                        
                        %Add to store
                        obj.SupportedArrangements(idx) = A;
                        
                    end
                    
                case 'remove'
                    
                    %Remove what ?
                    idx = i_select(obj.SupportedArrangements, name, 'Choose view(s) to remove...', ...
                        'SelectionMode', 'multiple');
                    if isempty(idx)
                        return;
                    end
                    
                    %Remove it
                    obj.SupportedArrangements(idx) = [];
                    
                case 'apply'
                    
                    %Apply what ?
                    idx = i_select(obj.SupportedArrangements, name);
                    
                    %Implement it
                    obj.Arrangement = obj.SupportedArrangements(idx);
                    
                case {'import', 'default'}
                    
                    %Import from where ?
                    if strcmp(action, 'default')
                        
                        %Go with this file
                        name = fullfile(fileparts(which(class(obj))), 'defaults.raf');
                        
                        %If not found
                        if exist(name, 'file') ~= 2
                            
                            %Return silently
                            return;
                            
                        end
                        
                    end
                    
                    %What name has been passed in ?
                    if isempty(name)
                        
                        %Not specified - ask the user
                        [f,p] = uigetfile(obj.ArrangementFileMask, 'Import arrangements from...', ...
                            char(obj.ArrangementFile));
                        
                        %Cancelled ?
                        if isempty(f) || isequal(f, 0)
                            return;
                        end
                        
                        %So the file is
                        fn = fullfile(p,f);
                        
                        %Load it
                        TEMP = load('-mat', fn, 'SupportedArrangements');
                        
                        %Validate
                        assert(isfield(TEMP, 'SupportedArrangements'), 'file does not contain a variable named ''SupportedArrangements''');
                        assert(isstruct(TEMP.SupportedArrangements), 'file does not contain a data structure named ''SupportedArrangements''');
                        
                        %Make a note
                        obj.ArrangementFile = fn;
                        
                        %And store the detail
                        obj.SupportedArrangements = TEMP.SupportedArrangements;
                        
                    elseif ischar(name)
                        
                        %Then treat it as a file whose content includes the ViewArrangement
                        TEMP = load('-mat', name, 'SupportedArrangements');
                        
                        %If nothing found
                        if ~isfield(TEMP, 'SupportedArrangements') || isempty(TEMP.SupportedArrangements)
                            
                            %Bail out silently
                            return;
                            
                        end
                        
                        %Make a note
                        obj.ArrangementFile = name;
                        
                        %And store the detail
                        obj.SupportedArrangements = TEMP.SupportedArrangements;
                        
                    elseif isstruct(name)
                        
                        %Assume it defines the arrangements directly
                        obj.ArrangementFile = [];
                        
                        %And store the detail
                        obj.SupportedArrangements = name;
                        
                    else
                        error('bad input');
                    end
                    
                    %And apply the first in the list
                    obj.Arrangement = obj.SupportedArrangements(1);
                    
                case 'export'
                    
                    %Save to where ?
                    if isempty(name)
                        
                        %Ask the user
                        [f,p] = uiputfile(obj.ArrangementFileMask, 'Export arrangement to file...', char(obj.ArrangementFile));
                        
                        %Cancelled ?
                        if isempty(f) || isequal(f, 0)
                            return;
                        end
                        
                        %So the file is
                        fn = fullfile(p,f);
                        
                        %Save it
                        SupportedArrangements = obj.SupportedArrangements; %#ok<NASGU,PROPLC>
                        save('-mat', fn, 'SupportedArrangements');
                        
                        %Make a note of file
                        obj.ArrangementFile = fn;
                        
                    else
                        
                        %Send back the answer
                        varargout{1} = Column;
                        
                    end
                    
                otherwise
                    error(['bad action ''', action, '''']);
            end
            
            function idx = i_select(A, name, prompt, varargin)
                
                %Valid ?
                assert(~isempty(A), 'no arrangements available');
                
                %Default prompt
                if nargin < 3 || isempty(prompt)
                    prompt = 'Select arrangement...';
                end
                
                %Caller may specify arrangement by name, index or not at all (in which case, ask)
                if nargin < 2 || isempty(name)
                    
                    %Ask
                    idx = listdlg('ListString', {A.Name}, ...
                        'Name', obj.Title, ...
                        'SelectionMode', 'single', ...
                        'PromptString', prompt, ...
                        varargin{:});
                    
                elseif isnumeric(name)
                    
                    %Specified by index directly
                    idx = name;
                    
                    %Validity ?
                    assert(idx > 0 && idx <= numel(A), 'out of range');
                    
                elseif ischar(name)
                    
                    %Specified by name
                    [~, idx] = ismember(name, {A.Name});
                    assert(idx ~= 0, 'unknown arrangement');
                    
                else
                    error('bad input');
                end
            end
            
        end
        
        function cbTabSelectionChange(obj, hc, ~)
            
            %Make a note
            obj.CurrentTab = hc.Selection;
            
        end
        
        function cbRemoveAxisData(obj, ~, ~, gdx)
            
            %Tell the object to clear the data from this axis next time
            %around
            obj.VariablesYLeft{gdx}  = [];
            obj.VariablesYRight{gdx} = [];
            obj.ResultsLeft{gdx}     = [];
            obj.ResultsRight{gdx}    = [];
            
            %Trigger an update of this element of the GUI
            updateElement(obj, gdx);
            
        end
        
        function cbShowLoadsEnvelope(obj, ~, ~)
            %cbShowLoadsEnvelope Calculates the loads envelope and opens a
            %new figure showing it...
            
            %Get view state
            Data = getEverythingOnce(obj);
            
            %Grab the beam properties
            BP = Data.BeamProperty;
            
            %Verify they are of type 'Internal Loads'
            assert(nnz(~ismember({BP.Type}, 'Internal Loads')) == 0, ...
                'Expected all of the BeamProperties to be of type ''Internal Loads''.');
            
            %Create colours 
            clr  = uniqueColours(numel(obj.Selection), {'w','k'});
            
            %No. of axes to plot data on?
            lim = numel(obj.hAxes);
            if lim > size(BP, 2)
                lim = size(BP, 2);
            end
            
            %Loop through axes and create the envelope
            for iAx = 1 : lim
                
                %Create matrix of values at common eta positions but force
                %the data in cell array form so we can check for gust data
                [eta, val] = obj.stackBeamPropertyData(BP(:, iAx), true);
                
                %Find loads which have two rows
                %   - Assume these are the max/min loads during a gust
                nRow   = cellfun(@(x) size(x, 1), val);
                idx    = (nRow == 2);
                nExtra = nnz(idx);
                
                %Grab the data in the 2nd row and then remove it
                col2     = cellfun(@(x) x(2, :), val(idx), 'Unif', false);
                val(idx) = cellfun(@(x) x(1, :), val(idx), 'Unif', false);
                
                %Append the new data & update colour index
                val(end +1 : end + nExtra) = col2;
                clr_ = [clr ; clr(idx, :)];
                
                %Combine into matrix to allow vectorised operations
                val = vertcat(val{:})';
                
                %Calculate the max/min values
                [maxVal, maxInd] = max(val, [], 2);
                [minVal, minInd] = min(val, [], 2);
                
                %Create single 'flipped' vector for plotting envelope
                eta_ = [eta ; flipud(eta)];
                env  = [maxVal ; flipud(minVal)];
                ind  = [maxInd ; minInd];
                uInd = unique(ind);
                
                %Plot the data
                plot(obj.hAxes(iAx), eta_, env, 'k--');

                %Plot the indiviual (colour-coded) markers but group them
                %to reduce the number of graphics objects
                for ii = 1 : numel(uInd)
                    
                    %Logical indexing
                    idx = ismember(ind, uInd(ii));
                    
                    %Make the markers
                    line(obj.hAxes(iAx)  , ...
                        'XData'          , eta_(idx), ...
                        'YData'          , env(idx) , ...
                        'Marker'         , 'o'     , ...
                        'LineStyle'      , 'none'  , ...
                        'MarkerEdgeColor', clr_(uInd(ii), :), ...
                        'MarkerFaceColor', clr_(uInd(ii), :));
                    
                end
                
            end
            
        end
        
    end
    
    methods % exporting data from the views
        
        function cbExportResults(obj, ~, ~, varargin)
            %cbExportResults Opens an export manager which allows the user
            %to choose which results sets they want to be exported and in
            %what format.
            
            %Make our own export manager and allow it to populate the
            %properties that control the export process
            awi.view.ResultsExporter(obj);
            
            %Do the export!
            if obj.ExportAWI   % .mat
                %export_mat(obj);
            end
            if obj.ExportCSV   % .csv
                [b, msg] = export_csv(obj);
            end
            if obj.ExportExcel %.xlsx
                [b, msg] = export_xls(obj);
            end
            
            %Revert to defaults
            obj.ExportAWI   = false;
            obj.ExportCSV   = false;
            obj.ExportExcel = false;
            obj.ExportPath  = '';
            obj.Res2Export  = [];
            obj.BeamNames2Export = {};
            obj.ResTypes2Export  = {};
            
        end
        
        function data2Export = prepateDataForExport(obj)
            %prepateDataForExpor Sorts and groups the results sets ready
            %for writing the data to a file.
            
            %Preliminary data
            nBeam2Export    = numel(obj.BeamNames2Export);
            nResType2Export = numel(obj.ResTypes2Export);
            
            %Grab the actual beam results objects
            BeamRes = arrayfun(@(x) x.BeamResults, obj.Res2Export, 'Unif', false);
            BeamRes = horzcat(BeamRes{:});
            
            %Grab the names of the beams
            allBeamNames = get([BeamRes.Beam], 'Name');
            
            %Preallocate
            temp_       = cell(1, nResType2Export);
            data2Export = cell(size(obj.BeamNames2Export));
            [data2Export{:}] = deal(temp_);
            
            %Group by beam name, result type & result quantity
            for i = 1 : nBeam2Export
                
                %Index by beam type
                idxBeam  = ismember(allBeamNames, obj.BeamNames2Export{i});
                BeamData = BeamRes(idxBeam);
                
                for j = 1 : nResType2Export
                    
                    %Grab the (hidden) beam property objects
                    BP = {BeamData.BeamProperties};
                    
                    bp = cell(1, numel(BP));
                    
                    for k = 1 : numel(BP)
                        
                        %What are the BeamProperty types?
                        type = {BP{k}.Type};
                        
                        %Filter by requested type
                        idx = ismember(type, obj.ResTypes2Export{j});
                        
                        %Keep a soft copy of these beam properties for
                        %further processing...
                        bp{k} = BP{k}(idx);
                        
                    end
                    
                    %Group by beam property name...
                    
                    %What are the names of the beam properties?
                    propNames = cellfun(@(x) {x.Name}, bp, 'Unif', false);
                    
                    %Unique beam property names
                    uniqueNames = unique(horzcat(propNames{:}), 'stable');
                    
                    %Find the common types across all beam results and return this
                    nNames = cellfun(@(x) numel(intersect(uniqueNames, x)), propNames);
                    
                    %Select the cell with the lowest number of interestions
                    [~, ind] = min(nNames);
                    commonPropNames  = propNames{ind};
                    
                    %Create a cell for each beam property
                    %   - Each beam property should have a value for each
                    %     results set.
                    for k = 1 : numel(commonPropNames)
                        data = cell(1, numel(bp));
                        for ii = 1 : numel(bp)
                            %Grab all beam properties across the results
                            %sets which are of this Name/Quantity
                            data{ii} = bp{ii}(ismember(propNames{ii}, commonPropNames{k}));
                            
                            %Assign to the export data
                            data2Export{i}{j}{k} = [data{:}];
                        end
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods % custom exporters
        
        function export_mat(obj)
            %export_mat Exports the raw AWI Framework objects in .mat
            %format.
            
            %Make a temp. variable so that we can save the data
            %   - N.B. Matlab really needs to sort out the saving process!
            AWIResults = obj.Res2Export;  %#ok<NASGU>
            
            %Make a filename
            %    - Timestamp the file to ensure uniqueness
            filename = sprintf('AWI_Results_%s.mat', datestr(now, 30));
            
            %Save in .mat format
            save(fullfile(obj.ExportPath, filename), 'AWIResults');
            
        end
        
        function [b, msg] = export_xls(obj)
            %export_xls Exports the results data into an Excel spreadsheet
            % file.
            %
            % 1 x Excel spreadsheet file is created for each beam element
            % and results type that has been requested.
            %
            %   e.g. Suppose the user has requested the results for the
            %   starboard wing and the VTP (nBeam = 2) and has asked for
            %   the Internal Loads and the Aerodynamic Loads (nType = 2)
            %   for both of these elements.
            %
            %   The number of files that are written will therefore be
            %   equal to nBeam x nType.
            %
            %   The individual results quantities are exported to different
            %   sheets in the various Excel spreadsheets/workbooks. For
            %   this reason, exporting data to an Excel sheet is much more
            %   concise than exporting to .csv.
            %
            %   HOWEVER!! The in-built Matlab function 'xlswrite' is much
            %   slower than the bespoke 'csvwrite' method in the AWI
            %   Framework. One of the reasons for this is that when
            %   exporting data to the Excel workbook the exported data is
            %   APPENDED to a workbook (if it already exists). This
            %   requires the excel file to be parsed to find the extent of
            %   the existing data before the new data can be written to the
            %   file.
            %
            % If a results quantity does not have any data (this is
            % represented by the 'Value' property of the 'BeamProperty'
            % object having all NaN) then the writing of that results
            % quantity will be skipped.
            
            %Has the export path been set?
            if isempty(obj.ExportPath)
                obj.ExportPath = pwd;
            end
            
            %Prepate data
            data2Export = prepateDataForExport(obj);
            
            %Grab beam names, result types & result set names
            beamNames = obj.BeamNames2Export;
            resTypes  = obj.ResTypes2Export;
            resNames  = get(obj.Res2Export, 'Name');
            if ~iscell(resNames) %Force cell
                resNames = {resNames};
            end
            
            %Write the files
            [b, msg] = i_write2xls(data2Export, resNames, beamNames, resTypes, obj.ExportPath);
            
            %If caller did not want anything back
            if nargout == 0
                
                %OK ?
                if b
                    
                    %Show us
                    winopen(obj.ExportPath);
                    
                else
                    
                    %What went wrong ?
                    error(msg);
                    
                end
                
            end
            
            function [b, msg] = i_write2xls(data2Export, resNames, beamNames, resTypes, ExportPath)
                %i_write2xls Writes the various .xlsx files.
                
                %Suppress Matlab warning for adding a new sheet
                warning('off', 'MATLAB:xlswrite:AddSheet');
                
                %Start with some defaults...
                b   = true;
                msg = '';
                
                %How many files are we about to write?
                nFile = numel(resNames) * numel(beamNames);
                
                %Tell the user about it
                hwb = waitbar(0, 'Exporting results. Hit ''Cancel'' to interrupt', ...
                    'CreateCancelBtn', 'delete(gcbf)', ...
                    'Name', 'AWI Framework - Exporting in Progress');
                clu = onCleanup(@()delete(hwb(ishghandle(hwb))));
                counter = 0;
                
                for iB = 1 : numel(beamNames) %For each beam ...
                    for iT  = 1 : numel(resTypes) %For each type of result
                        
                        counter = counter + 1;
                        
                        %Update the waitbar
                        waitbar(counter / nFile, hwb);
                        
                        %Construct the file name
                        xlsName = [beamNames{iB}, '_', resTypes{iT}, '.xlsx'];
                        xlsFile = fullfile(ExportPath, xlsName);
                        
                        %Check if this excel file already exists
                        if exist(xlsFile, 'file') == 2
                            newFile = false;
                        else
                            newFile = true;
                        end
                        
                        %If the excel file already exists grab sheet names
                        if ~newFile
                            [~, sheets] = xlsfinfo(xlsFile);
                        end
                        
                        for iP = 1 : numel(data2Export{iB}{iT}) %For each beam property
                            
                            %What beam property are we exporting?
                            propName = data2Export{iB}{iT}{iP}(1).Name;
                            qName    = data2Export{iB}{iT}{iP}(1).Quantity; %Only used if 'propName' is an invalid Excel sheet name
                            
                            %Index the data
                            data = data2Export{iB}{iT}{iP};
                            
                            %Check that we actually have anything to plot!
                            %   - Check for any NaN data...
                            %   - ...cannnot write NaN data so skip.
                            %   - Need additional 'any' term to accound for
                            %     data with multiple rows
                            val = {data.Value};
                            idx = cellfun(@(x) any(any(isnan(x))), val); 
                            
                            %Down-select
                            data = data(~idx);
                            nam  = resNames(~idx);
                            
                            if isempty(data) %Anything left?
                                continue
                            end
                            
                            %Grab data
                            eta = {data.Distribution};
                            val = {data.Value};
                            
                            %val = arrayfun(@(bp) num2cell(bp.Value')       , data, 'Unif', false);
                            %eta = arrayfun(@(bp) num2cell(bp.Distribution'), data, 'Unif', false);
                            
                            %Check for multidimensionsal data sets
                            nRows = cellfun(@(x) size(x, 1), val);
                            
                            %Cannot export data with more than 2 rows...
                            %   - TODO : Understand what situation would
                            %   lead to there being more than 2 rows. Time
                            %   series data?
                            if any(nRows > 2)
                                continue
                            end
                            
                            %Construct header string
                            etaStr = cellfun(@(rName) sprintf('%s - Eta (%s-%s-%s)', ...
                                propName, beamNames{iB}, resTypes{iT}, rName), ...
                                nam, 'Unif', false);
                            valStr = cellfun(@(rName) sprintf('%s (%s-%s-%s)', ...
                                propName, beamNames{iB}, resTypes{iT}, rName), ...
                                nam, 'Unif', false);
                            
                            %Any data with two rows should be an extra
                            %column in the output table - Use cell notation
                            %to expand the data set 
                            idx = (nRows == 2);
                            eta(~idx) = cellfun(@(x) {x}, eta(~idx), 'Unif', false);
                            val(~idx) = cellfun(@(x) {x}, val(~idx), 'Unif', false);
                            eta(idx)  = cellfun(@(x) {x, x}            , eta(idx) , 'Unif', false);
                            val(idx)  = cellfun(@(x) {x(1, :), x(2, :)}, val(idx) , 'Unif', false);
                            eta = horzcat(eta{:});
                            val = horzcat(val{:});
                            
                            %Need to create extra headers for the new data
                            etaStr(~idx) = cellfun(@(x) {x}, etaStr(~idx), 'Unif', false);
                            valStr(~idx) = cellfun(@(x) {x}, valStr(~idx), 'Unif', false);
                            etaStr(idx)  = cellfun(@(x) {[x, '(MAX)'], [x, '(MIN)']}, etaStr(idx), 'Unif', false);
                            valStr(idx)  = cellfun(@(x) {[x, '(MAX)'], [x, '(MIN)']}, valStr(idx), 'Unif', false);
                            etaStr = horzcat(etaStr{:});
                            valStr = horzcat(valStr{:});
                            
                            %Return to cell notation for writing to .xlsx
                            eta = cellfun(@(x) num2cell(x'), eta, 'Unif', false);
                            val = cellfun(@(x) num2cell(x'), val, 'Unif', false);
                            
                            %Check data lengths
                            nVal   = cellfun(@(x) numel(x), val);
                            maxVal = max(nVal);
                            
                            %Anything need padding?
                            idx = nVal < maxVal;
                            ind = find(idx == true);
                            if any(idx)
                                for i = 1 : numel(ind)
                                    %Pad the other lines with NaN
                                    pad = repmat({[]}, [(maxVal - nVal(ind(i))), 1]);                                    
                                    eta{ind(i)} = [eta{ind(i)} ; pad];
                                    val{ind(i)} = [val{ind(i)} ; repmat(pad, [1, size(val{ind(i)}, 2)])];
                                end
                            end
                            
                            %Combine...
                            etaAndVal  = [eta ; val];
                            data2Write = horzcat(etaAndVal{:});
                            
                            %...and header data as well!
                            temp    = [etaStr ; valStr];
                            header  = temp(:)';
                            
                            %Check that the header and the data have the
                            %same number of columns
                            %   - Different number of columns is possible
                            %     because some of the data has multiple 
                            %     rows (e.g. max/min values from a gust) so
                            %     we need to account for this!
                            if size(header, 2) ~= size(data2Write, 2)
                                header = [header, cell(1, size(data2Write, 2) - size(header, 2))];
                            end
                            
                            %Combine header & numeric data
                            ExcelData = [header ; data2Write];
                            
                            %Determine the starting cell for the data
                            if newFile
                                xlrange = 'A1';
                            else
                                %Check the sheet already exists...
                                if any(ismember(sheets, propName)) || any(ismember(sheets, qName))
                                    % ... if it does, need to find the
                                    % extent of the current data in the
                                    % sheet.
                                    
                                    %What is the sheet name?
                                    if any(ismember(sheets, propName))
                                        sheetName = propName;
                                    elseif any(ismember(sheets, qName))
                                        sheetName = qName;
                                    end
                                    
                                    %Read the data from the sheet
                                    [d1, ~, ~] = xlsread(xlsFile, sheetName);
                                    
                                    %How many columns?
                                    colNum = size(d1, 2) + 1;
                                    
                                    %Convert to a valid Excel column name
                                    colName = num2xlsCol(colNum);
                                    
                                    %Append the row number - Assume row 1
                                    xlrange = [colName, '1'];
                                else
                                    % ... if not, just start at A1!
                                    xlrange = 'A1';
                                end
                            end
                            
                            %Write the data to the file
                            [b, message] = xlswrite(xlsFile, ExcelData, propName, xlrange);
                            
                            if ~b %Exported okay?
                                switch message.identifier
                                    case 'MATLAB:xlswrite:InvalidSheetName'
                                        %Try a different sheet name
                                        [b, message] = xlswrite(xlsFile, ExcelData, qName);
                                        if ~b %Exported okay?
                                            %Grab the messge
                                            msg = message.message;
                                            %Re-enable warning
                                            warning('on', 'MATLAB:xlswrite:AddSheet');
                                            return
                                        end
                                    otherwise
                                        %Grab the messge
                                        msg = message.message;
                                        %Re-enable warning
                                        warning('on', 'MATLAB:xlswrite:AddSheet');
                                        return
                                end
                            end
                            
                        end
                    end
                end
                
                delete(hwb);
                
                %Re-enable warning
                warning('on', 'MATLAB:xlswrite:AddSheet');
                
            end
            
            function col = num2xlsCol(num)
                %num2xlscol Converts an index number into a valid Excel
                %column name.
                
                %26 letters in the alphabet
                base = 26;
                
                %Algorithm ...
                n   = ceil(log(num)/log(base));  % estimate number of digits
                d   = cumsum(base.^(0:n+1));     % offset
                n   = find(num >= d, 1, 'last'); % actual number of digits
                d   = d(n:-1:1);                 % reverse and shorten
                r   = mod(floor((num-d)./base.^(n-1:-1:0)), base) + 1;  % modulus
                col = char(r+64);  % convert number to ASCII
                
            end
            
        end
        
        function [b, msg] = export_csv(obj, varargin)
            %export_csv Exports the results data into a .csv file.
            %
            % 1 x .csv file is created for each results property that has
            % been requested.
            %
            %   e.g. Suppose the user has requested the results for the
            %   starboard wing and the VTP (nBeam = 2) and has asked for
            %   the Internal Loads and the Aerodynamic Loads (nType = 2)
            %   for both of these elements.
            %
            %   The number of files that are written is dependent on how
            %   many quantities are defined for various results types.
            %   Let us assume that Internal Loads has 5 quantities and
            %   Aerodynamic Loads have 3. For this case the export routine
            %   will generate 32 .csv files. Internal Loads will have
            %   5 x 2 x 2 files and Aerodynamaic Loads will have 3 x 2 x 2
            %   files.
            %
            %   In the case of writing .csv files it is easy to see how
            %   a large number of files can be generated if the user is not
            %   careful.
            %
            % If a results quantity does not have any data (this is
            % represented by the 'Value' property of the 'BeamProperty'
            % object having all NaN) then the writing of that results
            % quantity will be skipped.
            
            %Has the export path been set?
            if isempty(obj.ExportPath)
                obj.ExportPath = pwd;
            end
            
            %Prepate data
            data2Export = prepateDataForExport(obj);
            
            %Grab beam names, result types & result set names
            beamNames = obj.BeamNames2Export;
            resTypes  = obj.ResTypes2Export;
            resNames  = get(obj.Res2Export, 'Name');
            
            %Write the files
            [b, msg] = i_write2csv(data2Export, resNames, beamNames, resTypes, obj.ExportPath);
            
            %If caller did not want anything back
            if nargout == 0
                
                %OK ?
                if b
                    
                    %Show us
                    winopen(obj.ExportPath);
                    
                else
                    
                    %What went wrong ?
                    error(msg);
                    
                end
                
            end
            
            function [b, msg] = i_write2csv(data2Export, resNames, beamNames, resTypes, ExportPath)
                %i_write2csv Writes the various .csv files.
                
                %Start with some defaults...
                b   = true;
                msg = '';
                
                %How many files are we about to write?
                temp = cellfun(@(x) cellfun(@(y) numel(y), x), data2Export, 'Unif', false);
                nFile = sum(cellfun(@(x) sum(x), temp));
                
                %Tell the user about it
                hwb = waitbar(0, 'Exporting results. Hit ''Cancel'' to interrupt', ...
                    'CreateCancelBtn', 'delete(gcbf)', ...
                    'Name', 'AWI Framework - Exporting in Progress');
                clu = onCleanup(@()delete(hwb(ishghandle(hwb))));
                counter = 0;
                
                for iB = 1 : numel(beamNames) %For each beam ...
                    for iT  = 1 : numel(resTypes) %For each type of result
                        for iP = 1 : numel(data2Export{iB}{iT}) %For each beam property
                            
                            counter = counter + 1;
                            
                            %Update the waitbar
                            waitbar(counter / nFile, hwb);
                            
                            %What beam property are we exporting?
                            propName = data2Export{iB}{iT}{iP}(1).Name;
                            
                            %Index the data
                            data = data2Export{iB}{iT}{iP};
                            
                            %Check that we actually have anything to plot!
                            %   - Check for any NaN data...
                            %   - ...cannnot write NaN data so skip.
                            val = {data.Value};
                            idx = cellfun(@(x) any(any(isnan(x))), val);
                            
                            %Down-select
                            data = data(~idx);
                            nam  = resNames(~idx);
                            
                            if isempty(data) %Anything left?
                                continue
                            end
                            
                            %Grab eta and value of each property in 'data'
                            eta = {data.Distribution};
                            val = {data.Value};
                            
                            %Check for multidimensionsal data sets
                            nRows = cellfun(@(x) size(x, 1), val);
                            
                            %Cannot export data with more than 2 rows...
                            %   - TODO : Understand what situation would
                            %   lead to there being more than 2 rows. Time
                            %   series data?
                            if any(nRows > 2)
                                continue
                            end
                            
                            %Construct header string
                            etaStr = cellfun(@(rName) sprintf('%s - Eta (%s-%s-%s)', ...
                                propName, beamNames{iB}, resTypes{iT}, rName), ...
                                nam, 'Unif', false);
                            valStr = cellfun(@(rName) sprintf('%s (%s-%s-%s)', ...
                                propName, beamNames{iB}, resTypes{iT}, rName), ...
                                nam, 'Unif', false);                            
                            
                            %Any data with two rows should be an extra
                            %column in the output table - Use cell notation
                            %to expand the data set 
                            idx = (nRows == 2);
                            eta(~idx) = cellfun(@(x) {x}, eta(~idx), 'Unif', false);
                            val(~idx) = cellfun(@(x) {x}, val(~idx), 'Unif', false);
                            eta(idx)  = cellfun(@(x) {x, x}            , eta(idx) , 'Unif', false);
                            val(idx)  = cellfun(@(x) {x(1, :), x(2, :)}, val(idx) , 'Unif', false);
                            eta = horzcat(eta{:});
                            val = horzcat(val{:});
                            
                            %Need to create extra headers for the new data
                            etaStr(~idx) = cellfun(@(x) {x}, etaStr(~idx), 'Unif', false);
                            valStr(~idx) = cellfun(@(x) {x}, valStr(~idx), 'Unif', false);
                            etaStr(idx)  = cellfun(@(x) {[x, '(MAX)'], [x, '(MIN)']}, etaStr(idx), 'Unif', false);
                            valStr(idx)  = cellfun(@(x) {[x, '(MAX)'], [x, '(MIN)']}, valStr(idx), 'Unif', false);
                            etaStr = horzcat(etaStr{:});
                            valStr = horzcat(valStr{:});
                            
                            %Concatenate any data with precisely 2 rows
                            %eta(idx) = cellfun(@(x) [x, x]            , eta(idx), 'Unif', false);
                            %val(idx) = cellfun(@(x) [x(1, :), x(2, :)], val(idx), 'Unif', false);
                            
                            %Convert data to strings ready for writing
                            eta = cellfun(@(x) cellstr(num2str(x')), eta, 'Unif', false);
                            val = cellfun(@(x) cellstr(num2str(x')), val, 'Unif', false);
                            
                            %Check data lengths
                            nVal   = cellfun(@(x) numel(x), val);
                            maxVal = max(nVal);
                            
                            %Anything need padding?
                            idx = nVal < maxVal;
                            ind = find(idx == true);
                            if any(idx)
                                for i = 1 : numel(ind)
                                    %Pad the other lines with empty cells
                                    pad = repmat({''}, [maxVal - nVal(ind(i)), 1]);
                                    eta{ind(i)} = [eta{ind(i)} ; pad];
                                    val{ind(i)} = [val{ind(i)} ; pad];
                                end
                            end
                            
                            %Combine...
                            etaAndVal  = [eta ; val];
                            etaAndVal  = etaAndVal(:)';
                            etaAndVal  = horzcat(etaAndVal{:});
                            data2Write = arrayfun(@(i) strjoin(etaAndVal(i, :), ','), 1 : maxVal, 'Unif', false);
                            
                            %...header as well!
                            temp    = [etaStr, valStr]';
                            header  = temp(:);
                            hString = strjoin(header, ',');
                            
                            %Construct the file name
                            csvName = [beamNames{iB}, '_', resTypes{iT}, '_', propName, '.csv'];
                            csvFile = fullfile(ExportPath, csvName);
                            
                            %Open the file & write the data
                            fid = fopen(csvFile, 'w');
                            if fid == -1 %Check the file was opened
                                b = false;
                                msg = sprintf(['Unable to open file ''%s'' ', ...
                                    'for writing data.'], csvFile);
                                return
                            end
                            fprintf(fid, '%s\r\n', hString);       %Header
                            fprintf(fid, '%s\r\n', data2Write{:}); %Data
                            fclose(fid);
                            
                        end
                    end
                end
                
                delete(hwb);
                
            end
            
        end
        
    end
    
    methods (Access = private) % helper functions
        function Data = getEverythingOnce(obj, Res)
            %getEverythingOnce Grabs quantities relating to the results,
            %their underlying beam elements and 'Beamable' results sets in
            %order to speed up the code. Also, validates the selection
            %values for the various uicontrol objects.
            
            %Start with a blank structure
            Data = struct();
            
            if nargin < 2
                Res = obj.Results;
            end
            
            %% Grab the results
            
            %Name of all results
            Data.ResultNames  = {Res.Name};
            
            %Current selection in the "Analysis Results" list-box.
            sel = obj.Selection;
            
            if isnumeric(sel)   %Validate
                %Ensure valid
                sel(sel < 1 |sel > numel(Res)) = [];
            elseif isa(sel, 'mvc.mixin.ResultSet')
                %Convert to numeric
                [b, sel] = ismember(sel, Res);
                sel = sel(b);
            else
                %Not valid
                sel = [];
            end
            
            %Currently selected results
            Data.CurrentResult     = Res(sel);
            Data.CurrentResultName = {Res(sel).Name};
            Data.ResultSelection   = sel;
            
            %% Underlying 'awi.model.Beam' objects...
            
            %Grab the handle to the 'awi.model.Beam' object
            allBeams = arrayfun(@(res) [res.BeamResults.Beam], Data.CurrentResult, 'Unif', false);
            
            %Grab the names & fullnames
            allBeamNames  = cellfun(@(x) {x.Name}    , allBeams, 'Unif', false);
            allBeamFNames = cellfun(@(x) {x.FullName}, allBeams, 'Unif', false);
            
            %Intersect the list of beams to find the common beam names
            uniqueBeamNames = unique(horzcat(allBeamFNames{:}), 'stable');
            
            %Find the common types across all beam results and return this
            nNames = cellfun(@(x) numel(intersect(uniqueBeamNames, x)), allBeamFNames);
            
            %Select the cell with the lowest number of interestions
            [~, ind] = min(nNames);
            
            %Check for empties
            if isempty(ind)
                Data.AvailableBeams     = allBeams;
                Data.AvailableBeamNames = {'-Select a Results Set-'};
            else
                Data.AvailableBeamNames = allBeamNames{ind};
                Data.AvailableBeams     = allBeams{ind};
            end
            
            %Current selection in the "Beam Results" list-box
            beamSel = obj.BeamSelection;
            
            if isnumeric(beamSel)   %Validate
                %Ensure valid
                beamSel(beamSel < 1 |beamSel > numel(Data.AvailableBeamNames)) = [];
            elseif ischar(beamSel) || iscellstr(beamSel)
                %Convert to numeric
                [b, beamSel] = ismember(beamSel, Data.AvailableBeamNames);
                beamSel = beamSel(b);
            else
                %Not valid
                beamSel = [];
            end
            
            %Currently selected beams
            Data.BeamSelection = beamSel;
            
            %Name of the currently selected beam
            if isempty(Data.AvailableBeams)
                Data.CurrentBeam         = [];
                Data.CurrentBeamName     = '';
                Data.CurrentBeamFullName = '';
            else
                Data.CurrentBeam         = Data.AvailableBeams(Data.BeamSelection);
                Data.CurrentBeamName     = Data.CurrentBeam.Name;
                Data.CurrentBeamFullName = Data.CurrentBeam.FullName;
            end
            
            %% Underlying 'awi.model.BeamResult' objects...
            
            %All 'awi.model.BeamResult' objects defined for this beam
            Data.AllBeamResults = {Data.CurrentResult.BeamResults};
            
            %The 'awi.model.BeamResult' objects that are currently selected
            %within the "Beam Results" list-box
            currentBeamRes = cell(size(Data.AllBeamResults));
            for i = 1 : numel(currentBeamRes)
                beamNames = get([Data.AllBeamResults{i}.Beam], 'FullName');
                if ~iscell(beamNames) %Force cell
                    beamNames = {beamNames};
                end
                idx = ismember(beamNames, Data.CurrentBeamFullName);
                currentBeamRes{i} = Data.AllBeamResults{i}(idx);
            end
            %Return an object array
            Data.CurrentBeamResult = horzcat(currentBeamRes{:});
            
            %% Results types belonging to the currently selected Beam
            
            %Grab the available results types
            if isempty(Data.CurrentBeamResult)
                
                Data.ResultTypes = {'-Select a Results Set-'};
                
            else
                
                %Grab the name of all BeamProperty types and the unique values
                allTypes    = {Data.CurrentBeamResult.BeamPropertyTypes};
                uniqueTypes = unique(horzcat(allTypes{:}), 'stable');
                
                %Find the common types across all beam results and return this
                nTypes = cellfun(@(x) numel(intersect(uniqueTypes, x)), allTypes);
                
                %Select the cell with the lowest number of interestions
                [~, ind] = min(nTypes);
                Data.ResultTypes = allTypes{ind};
                
            end
            
            %Current selection in the "Results Types" drop-down menu
            resTypeSel = obj.PlotTypeSelection;
            
            if isnumeric(resTypeSel)   %Validate
                %Ensure valid
                resTypeSel(resTypeSel < 1 |resTypeSel > numel(Data.ResultTypes)) = [];
            elseif ischar(resTypeSel) || iscellstr(resTypeSel)
                %Convert to numeric
                [b, resTypeSel] = ismember(resTypeSel, Data.AvailableBeamNames);
                resTypeSel = resTypeSel(b);
            else
                %Not valid
                resTypeSel = [];
            end
            
            %What is the currently selected results type?
            Data.CurrentResultType   = Data.ResultTypes(resTypeSel);
            Data.ResultTypeSelection = resTypeSel;
            
            %% Results quantities belonging to the current results type
            
            %Grab results quantities
            if isempty(Data.CurrentBeamResult)
                
                Data.AllResultQuantities     = {'-Select a Results Set-'};
                Data.CurrentResultQuantity   = {''};
                Data.ResultQuantitySelection = 1;
                
            else
                
                %Grab the underlying 'awi.mixin.BeamProperty' objects
                BP = getBeamProperty(Data.CurrentBeamResult(1), Data.CurrentResultType);
                
                %Filter out any data that is just NaN
                idx = ~cellfun(@(x) any(any(isnan(x))), {BP.Value});
                
                %Check to see if there is any data left!
                if ~any(idx)
                    
                    %Tell the user...
                    Data.AllResultQuantities     = {'-No Data-'};
                    Data.CurrentResultQuantity   = {'-No Data-'};
                    Data.ResultQuantitySelection = 1;
                    
                else
                    
                    %Return all the names of the results quantities in a cell array
                    Data.AllResultQuantities = {BP(idx).Name};
                    
                    %Current selection in the "Results Quantities" list-box
                    resQuantSel = obj.PlotQuantitySelection;
                    
                    if isnumeric(resQuantSel)   %Validate
                        %Ensure valid
                        resQuantSel(resQuantSel < 1 |resQuantSel > numel(Data.AllResultQuantities)) = [];
                    elseif ischar(resQuantSel) || iscellstr(resQuantSel)
                        %Convert to numeric
                        [b, resQuantSel] = ismember(resQuantSel, Data.AvailableBeamNames);
                        resQuantSel = resQuantSel(b);
                    else
                        %Not valid
                        resQuantSel = [];
                    end
                    
                    %Grab the name of the currently selected beam property
                    Data.CurrentResultQuantity   = Data.AllResultQuantities(resQuantSel);
                    Data.ResultQuantitySelection = resQuantSel;
                    
                end
                
            end
            
            %% Underlying 'awi.mixin.BeamProperty' objects
            
            %Get the beam properties
            if isempty(Data.CurrentBeamResult) || isempty(Data.CurrentResultQuantity)
                
                Data.BeamProperty = [];
                
            else
                
                %Stack ALL of the beam properties
                bp  = vertcat(Data.CurrentBeamResult.BeamProperties);
                
                %Down-select to 'awi.model.BeamProperty' objects of type
                %'obj.CurrentResultsType' && 'obj.CurrentResultsQuantities'
                idx = and( ...
                    ismember({bp(1, :).Type}, Data.CurrentResultType{1}), ...
                    ismember({bp(1, :).Name}, Data.CurrentResultQuantity)); %Make sure we filter by 'Name' and not 'Quantity' as the list view (Results Quantities) is populated using the Names!
                Data.BeamProperty = bp(:,idx);
                
            end
            
            %% Format the selection values
            %   - Cannot have any empty selection values as this will cause
            %     the uicontrol objects to dissapear
            %   - Force a selection value of 1.
            
            if isempty(Data.BeamSelection)
                Data.BeamSelection = 1;
            end
            if isempty(Data.ResultTypeSelection)
                Data.ResultTypeSelection = 1;
            end
            if isempty(Data.ResultQuantitySelection)
                Data.ResultQuantitySelection = 1;
            end
            
        end
    end
    
    methods (Static) % helper functions
        function [x, y] = stackBeamPropertyData(BeamProperty, forceCell)
            %stackBeamPropertyData Interpolates the data for each
            %BeamProperty to a common set of eta positions and returns
            %the data in a format suitable for plotting or performing
            %vectorised operations on.
            %
            % If a BeamProperty has more than 2 rows in its 'Value'
            % property then a row vector is formed. An error is output
            % if more than 2 rows are present.
            
            if nargin < 2
                forceCell = false;
            end
            
            %Interpolate the data
            eta  = getEta(BeamProperty);
            data = arrayfun(@(x) getBPV(x, x.Quantity, eta), ...
                BeamProperty, 'Unif', false);
            
            %Force cell output of half-processed data
            %   - Useful for calculating the envelope correctly.
            if forceCell
                x = eta';
                y = data;
                return
            end
            
            %Check for data that has multiple rows
            nRows = cellfun(@(x) size(x, 1), data);
            
            if any(nRows > 2) %Handle the case where transient data has been accidentally provided...
                error('The ''BeamResultsViewer'' class cannot handle data with more than 2 rows.');
            end
            
            %Assume any data with 2 rows is the max/min loads for the
            %gust envelope
            %   - TODO : Is this acceptable? How will this effect the
            %   exporting of results? Surely this is not robust enough!
            ind = nRows == 2;
            if any(ind)
                
                %To continue to enable vectorised plotting we must
                %reflect ALL data, not just the ones with 2 rows.
                data(ind)  = cellfun(@(x) [x(1, :), fliplr(x(2, :))], data(ind), 'Unif', false);
                data(~ind) = cellfun(@(x) [x, nan(size(x))], data(~ind), 'Unif', false);
                
                %Eta vector must be extended
                eta = [eta, fliplr(eta)];
                
            end
            
            %Plotting code expects the (x,y) data to be in column
            %format -> size = [nEta, nBeam]
            x = eta';
            y = vertcat(data{:})';
            
        end        
    end
    
end


