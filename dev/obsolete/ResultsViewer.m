classdef ResultsViewer < mvc.view.Container & mvc.mixin.Exportable
    %
    % Generic awi results viewer
    
    properties
        
        %The model to which we are attached
        Model;
        
        %Function to call on selection change (not used, but required for compatibility with ViewManager)
        SelectionChangeFcn;
        
        %Track currently selected item(s) (not used, but required for compatibility with ViewManager)
        Selection;
        
        %Numbers of rows and columns into which results are arranged
        Rows = 1;
        Cols = 1;
        
        %List of all variables that can be viewed
        AllVariables = {'Eta', 'Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz'};
        
        %Lists of variables that are actually being viewed
        VariablesX = {};
        VariablesY = {};
        
        %What arrangements of results are supported by this viewer ?
        SupportedArrangements;
        
        %Which arrangement is currently selected ?
        CurrentArrangement;
        
        %What styles are supported ?
        SupportedStyles = {'line', 'scatter', 'table'};
        
        %What style is currently selected ?
        CurrentStyle;
        
        %What layouts are supported ?
        SupportedLayouts = {'grid', 'tabs'};
        
        %What layout is currently selected ?
        CurrentLayout;
        
        %What tab is currently selected (if any) ?
        CurrentTab = 1;
        
    end
    
    properties (Dependent)
        
        Aircraft;
        LoadCases;
        Results;
        
        Arrangement;
        CurrentLayoutIndex;
        
    end
    
    properties (Access = protected) % To help persist view arrangements
        
        ArrangementFileMask = {'*.raf', 'Results viewer arrangement files (*.raf)'; ...
            '*.mat', 'MATLAB data files (*.mat)'; ...
            '*.*', 'All files (*.*)'};
        ArrangementFile;
        
    end
    
    properties (Access = protected)
        
        %List of pushbutton labels and functions, to be displayed in view
        Pushbuttons = {};
        
        %Results selection
        hResults;
        
        %The panel in which we display the results
        hPanel;
        
        %Axes embedded within panel
        hAxes;
        
        %Pushbutton(s) to control the analysis (may or may not be relevant, depending on analyis)
        hButtons;
        
        %Placeholder for listeners
        Listeners;
        
        %Controls displaying numbers of rows / columns, layout
        hRows;
        hCols;
        hLayout;
        
        %Somewhere to attach a single legend
        % (avoiding duplication and freeing up real-estate)
        hLegend;
        
    end
    
    methods % get/set
        
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
            
            %Make a note of variables
            val.VariablesX = obj.VariablesX;
            val.VariablesY = obj.VariablesY;
            
        end
        
        function set.Selection(obj, val)
            
            %Selection stores as index into list of available results
            if isnumeric(val)
                
                %Ensure valid
                assert(isempty(val) || all(ismember(val, 1:numel(obj.Results))), 'invalid selection'); %#ok<MCSUP>
                
            elseif isa(val, 'awi.model.ResultSet')
                
                %Ensure valid
                [b,idx] = ismember(val, obj.Results); %#ok<MCSUP>
                assert(all(b), 'invalid selection');
                val = idx(b);
                val = val(:).';
                
            else
                
                %This type of view can't present whatever is being selected.
                % COULD issue a warning and continue
                % warning(['selection set to class ''', class(val), '''']);
                %
                % but better I think to just continue, with no change, silently
                return;
                
            end
            
            %Make a note
            obj.Selection = val;
            
            %Force an update
            update(obj);
            
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
            obj.Rows = val.Rows;
            obj.Cols = val.Cols;
            obj.VariablesX = val.VariablesX;
            obj.VariablesY = val.VariablesY;
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.VariablesX(obj)
            
            %Start here
            val = obj.VariablesX;
            
            %But if number of elements is less than expected from rows/cols
            if numel(val) < obj.Rows * obj.Cols
                
                %Pad with empties
                val(end+1:end+(obj.Rows*obj.Cols)-numel(val)) = cell(1,obj.Rows*obj.Cols-numel(val));
                
            end
            
        end
        
        function val = get.VariablesY(obj)
            
            %Start here
            val = obj.VariablesY;
            
            %But if number of elements is less than expected from rows/cols
            if numel(val) < obj.Rows * obj.Cols
                
                %Pad with empties
                val(end+1:end+(obj.Rows*obj.Cols)-numel(val)) = cell(1,obj.Rows*obj.Cols-numel(val));
                
            end
            
        end
        
        function val = get.LoadCases(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no LoadCases
                val = [];
                
            else
                
                %What have we got ?
                val = findall(obj.Model, 'Type', 'LoadCase');
                
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
        
        function val = get.Results(obj)
            
            %No model ?
            if isempty(obj.Model)
                
                %Then no Results
                val = [];
                
            else
                
                %Pass it on
                val = findall(obj.Model,'isa','awi.model.ResultSet');
                
            end
            
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
        
        function obj = ResultsViewer( model, varargin )
            
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
            
            %Initialise Exportable properties
            obj.ExportSheetTag = [];
            obj.ExportSpec(end+1,:) = {'*.csv', 'Comma-separated variable files (*.csv)', @export_csv};
            obj.sortExportSpec({'xls', 'csv', 'mat'});
            
            %Assign values at this level, if any were passed in
            if ~isempty(prp)
                pav = [prp; val];
                set(obj, pav{:});
            end
            
            %Start with an HBox
            hh = uiextras.HBoxFlex('Parent', obj.UIContainer, 'Spacing', 6, 'Padding', 6);
            
            %Add a panel in which the results are displayed
            obj.hPanel = uiextras.Panel('Parent', hh, 'Padding', 6);
            
            %Add a vertical strip into which controls are embedded
            hv = uiextras.VBox('Parent', hh, 'Padding', 6, 'Spacing', 6);
            
            %Good width
            hh.Widths(end) = 100;
            
            %Good height for simple text labels ?
            hLab = 20;
            
            %A panel in which the results are selected
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Results');
            hv.Heights(end) = hLab;
            obj.hResults = uicontrol('style', 'listbox', ...
                'Parent',hv,...
                'min',0,'max',2,...
                'Callback',@obj.onSelectionChange);
            hv.Heights(end) = 100;
            
            %Add a popup for layout
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Layout');
            hv.Heights(end) = hLab;
            obj.hLayout = uicontrol('Parent', hv,...
                'Style', 'popup',...
                'Value', obj.CurrentLayoutIndex, ...
                'String', obj.SupportedLayouts, ...
                'Callback', @obj.cbLayout);
            hv.Heights(end) = 25;
            
            %Edit box for Rows
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Rows');
            hv.Heights(end) = hLab;
            obj.hRows = uicontrol('Parent', hv,...
                'Style', 'edit',...
                'String', num2str(obj.Rows), ...
                'Callback', @obj.cbRows);
            hv.Heights(end) = 25;
            
            %Edit box for Columns
            uicontrol('Parent', hv, 'Style', 'Text', 'String', 'Cols');
            hv.Heights(end) = hLab;
            obj.hCols = uicontrol('Parent', hv,...
                'Style', 'edit',...
                'String', num2str(obj.Cols), ...
                'Callback', @obj.cbCols);
            hv.Heights(end) = 25;
            
            %If pushbuttons also required
            if ~isempty(obj.Pushbuttons)
                
                %Add a button box
                % TODO - Cater for controls of different sizes
                obj.hButtons = uiextras.VButtonBox('Parent', hh, 'Spacing', 6);
                
                %Add the buttons
                for i = 1:2:numel(obj.Pushbuttons)
                    
                    %Add the button
                    uicontrol('Parent', obj.hButtons, 'Style', 'Pushbutton', ...
                        'String', obj.Pushbuttons{i}, ...
                        'Callback', {@obj.onAnalyse, obj.Pushbuttons{i + 1}});
                    
                end
                
            end
            
            %And a container, to be used later to attach legend
            obj.hLegend = uicontainer('Parent', hv); % axes('Parent', hv, 'Visible', 'off');
            
            %Can we find a default view arrangements file ?
            arrangement(obj, 'default');
            
            %Store the model
            obj.Model = model;
            
        end
        
    end
    
    methods ( Access = protected )
        
        function update(obj)
            
            %Clean sheet
            delete(obj.hAxes);
            
            %Anything to do ?
            if isempty(obj.Results)
                
                %Ensure listbox empty
                set(obj.hResults, 'String', {});
                
                %Ensure panel has no content
                delete(obj.hPanel.Children);
                
                %We're done
                return;
                
            end
            
            %Selection must match available results
            sel = obj.Selection;
            if isnumeric(sel)
                
                %Ensure valid
                sel(sel < 1 | sel > numel(obj.Results)) = [];
                
            elseif isa(sel, 'mvc.mixin.ResultSet')
                
                %Convert to numeric
                [b, sel] = ismember(sel, obj.Results);
                sel = sel(b);
                
            else
                
                %Not valid
                sel = [];
                
            end
            
            %List available results and show selection
            set(obj.hResults, 'String', {obj.Results.Name}, 'Value', sel);
            
            %Layout and dimensions of grid
            set(obj.hLayout, 'String', obj.SupportedLayouts, 'Value', obj.CurrentLayoutIndex);
            set(obj.hRows, 'String', num2str(obj.Rows));
            set(obj.hCols, 'String', num2str(obj.Cols));
            
            %Results are viewed in this layout
            switch obj.CurrentLayout
                
                case 'grid'
                    
                    %Simple grid
                    hg = uiextras.Grid('Parent', obj.hPanel);
                    
                case 'tabs'
                    
                    %Make a note of current selection before overwriting
                    nsel = obj.CurrentTab;
                    
                    %Separate tab for each result
                    hg = uiextras.TabPanel('Parent', obj.hPanel, 'TabWidth', 100, 'SelectionChangedFcn', @obj.cbTabSelectionChange);
                    
                otherwise
                    error('bad layout');
            end
            
            %Create one axes for each result
            ha = arrayfun(@(x)axes( ...
                'Parent', hg, ... % uicontainer('Parent', hg), ... no need to insert uicontainer, if legend handled separately
                'UserData', x), ...
                1:obj.Rows*obj.Cols, 'Unif', false);
            
            %Store in object (maybe needed later, e.g. dynamic update of visible state based on legend click)
            obj.hAxes = [ha{:}];
            
            %Results are viewed in this layout
            switch obj.CurrentLayout
                
                case 'grid'
                    
                    %Arranged accordingly
                    hg.Widths = -ones(1,obj.Cols);
                    
                case 'tabs'
                    
                    %Good tab titles (TODO: generalise this)
                    ttit = cell(1,obj.Rows*obj.Cols);
                    for i = 1:obj.Rows*obj.Cols
                        ttit{i} = [obj.VariablesY{i}, '(', obj.VariablesX{i}, ')'];
                    end
                    hg.TabTitles = ttit;
                    
                    %Currently selected tab ?
                    if nsel >= 1 && nsel <= numel(ttit)
                        hg.Selection = nsel;
                        obj.CurrentTab = nsel;
                    end
                    
                otherwise
            end
            
            %For each element in grid
            for i = 1:obj.Rows*obj.Cols
                
                %User helper function to populate grid
                leg = updateElement(obj, i, sel);
                
            end

            %Use helper function to maintain legend
            updateLegend(obj, leg);
            
        end
        
        function [leg, msg] = updateElement(obj, i, sel)
            
            %Placeholders for lines, legends and messages
            hl = {};
            leg = {};
            msg = {};
            
            %Careful
            try
                
                %Context menu applied to all results
                hcm = uicontextmenu('Parent', ancestor(obj.hAxes(i),'Figure'), 'CallBack', {@obj.cbVariables, i});
                uimenu(hcm, 'Label', 'X');
                uimenu(hcm, 'Label', 'Y');
                uimenu(hcm, 'Label', 'Style');
                uimenu(hcm, 'Separator', 'on', 'Label', 'Arrangement');
                uimenu(hcm, 'Separator', 'on', 'Label', 'Export...', 'Callback', @obj.onExport);
                set(obj.hAxes(i), 'UiContextMenu', hcm);
                
                %Hold, and apply labels
                hold(obj.hAxes(i), 'on');
                xlabel(obj.hAxes(i), obj.VariablesX{i});
                ylabel(obj.hAxes(i), obj.VariablesY{i});
                
                %Do we know what to plot ?
                if isempty(obj.VariablesX{i}) || isempty(obj.VariablesY{i})
                    
                    %No
                    return;
                    
                end
                
                %For each selected result
                for j = sel
                    
                    %Extract x and y data, with care
                    [x, y, msg{end+1}] = getVariables(obj.Results(j), obj.VariablesX{i}, obj.VariablesY{i});
                    
                    %Problem ?
                    if ~isempty(msg{end})
                        
                        %Move on
                        continue;
                        
                    end
                    
                    %Present the data
                    switch obj.CurrentStyle{i}
                        
                        case 'line'
                            
                            %Simple
                            hl{end+1} = plot(obj.hAxes(i), x, y);
                            
                        case 'scatter'
                            
                            %Some care required
                            if size(x,2) == 1
                                htemp = arrayfun(@(k)scatter(obj.hAxes(i), x, y(:,k)), 1:size(y,2),'UniformOutput',false);
                                hl{end+1} = vertcat(htemp{:});
                            elseif size(x,2) == size(y,2)
                                htemp = arrayfun(@(k)scatter(obj.hAxes(i), x(:,k), y(:,k)), 1:size(y,2),'UniformOutput',false);
                                hl{end+1} = vertcat(htemp{:});
                            else
                                error('size mis-match');
                            end
                            
                        otherwise
                            
                            %TODO
                            hl{end+1} = [];
                            
                    end
                    
                    %Extend legend string (TODO properly)
                    leg{end+1} = obj.Results(j).Legend;
                    
                end
                
                %Expand lines and legend
                hl = vertcat(hl{:});
                leg = [leg{:}];
                
                %NO - legend on each individual axes takes up too much real estate
                %if ~isempty(leg)
                %    hleg = legend(obj.hAxes(i), leg, 'Interpreter', 'none');
                %end
                %
                if iscell(leg)
                %Instead, assign legend strings as display names (used later)
                [hl.DisplayName] = deal(leg{:});
                end
                
                %Any problems to report ?
                b = ~cellfun(@isempty, msg);
                if any(b)
                    
                    %Consolidate
                    msg = arrayfun(@(x)[obj.Results(x).Name, ': ', msg{x}], find(b), 'UniformOutput', false);
                    
                    %Add to plot
                    xl = xlim(obj.hAxes(i));
                    yl = ylim(obj.hAxes(i));
                    text(obj.hAxes(i), mean(xl), mean(yl), msg, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle');
                    
                end
                
            catch err
                
                %What went wrong ?
                xl = xlim(obj.hAxes(i));
                yl = ylim(obj.hAxes(i));
                text(obj.hAxes(i), mean(xl), mean(yl), ...
                    {'View creation FAILED with error:', err.message}, ...
                    'VerticalAlignment', 'middle', ...
                    'HorizontalAlignment', 'center');
                
            end
            
        end
        
        function updateLegend(obj, leg)
            %Attempt to display a single legend somewhere more convenient than embedded in each chart
            
            %Ensure clean sheet in container to which legend is attached
            delete(obj.hLegend.Children);
        
            %Anything else to do ?
            if isempty(leg)
                
                %No
                return;
                
            end
            
            %Copy last axes created into the hLegend container
            ha_copy = copyobj(obj.hAxes(end), obj.hLegend);
            
            %But need to tweak units and position, and hide, the helper axes
            set(ha_copy, 'Units', 'normalized', 'Position', [0, 0, 1, 1], 'Visible', 'off');
            
            %Get content
            hc_copy = ha_copy.Children;
            
            %Eliminate any text (wehich will be error messages, and have no legend associated)
            hc_copy(arrayfun(@(x)isa(x, 'matlab.graphics.primitive.Text'), hc_copy)) = [];
            
            %Hide the rest, but without setting visible off
            % (because that causes the legend to respond)
            arrayfun(@(x)set(x, 'XData', nan, 'YData', nan), hc_copy);
            
            %Add legend to this copy, with custom hit callback
            hleg = legend(ha_copy, leg, ...
                'Interpreter', 'none', ...
                'ItemHitFcn', @obj.legendItemHit); %% {@i_legendItemHit, obj.hAxes});
                
        end
            
        function legendItemHit(obj, hc, evtdata)
            
            %Which line are we ?
            nam = evtdata.Peer.DisplayName;
            disp([nam, ' hit...']);
            
            %Find line in the "dummy" plot used to create legend
            hl = findobj(hc.Parent, 'DisplayName', nam);
            
            %Get visible state
            bVis = strcmp(hl.Visible, 'on');
            
            %Toggle it
            vis = obj.bool2offon(~bVis);
            hl.Visible = vis;
            
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
        
        function onLoadCaseChange(obj, hc, ~)
            
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
        
        function onAnalyse(obj, hc, ~, fcn)
            
            %Might take a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Do it
                fcn(obj);
                
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
                
                %Make a note of new selection
                obj.Selection = get(obj.hResults, 'Value');
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbRows(obj, hc, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note
                obj.Rows = str2num(get(obj.hRows,'String'));
                
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
                obj.Cols = str2num(get(obj.hCols,'String'));
                
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
                obj.CurrentLayout = popupstr(obj.hLayout);
                
                %Force an update
                update(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function cbVariables(obj, hc, ~, gdx)
            
            %List of all variables
            vars = obj.AllVariables;
            
            %Choice of X variable
            var = obj.VariablesX{gdx};
            
            %Displayed here
            hx = findobj(hc, 'Label', 'X');
            
            %Clean sheet
            delete(hx.Children);
            
            %Rebuild list
            cellfun(@(x)uimenu(hx, 'Label', x, ...
                'Check', mvc.mixin.UiTools.bool2offon(strcmp(var, x)), ...
                'CallBack', {@obj.cbVariable, 'VariablesX', gdx, x}), vars);
            
            %Choice of Y variable
            var = obj.VariablesY{gdx};
            
            %Displayed here
            hy = findobj(hc, 'Label', 'Y');
            
            %Clean sheet
            delete(hy.Children);
            
            %Rebuild list
            cellfun(@(x)uimenu(hy, 'Label', x, ...
                'Check', mvc.mixin.UiTools.bool2offon(strcmp(var, x)), ...
                'CallBack', {@obj.cbVariable, 'VariablesY', gdx, x}), vars);
            
            %Styles handled by separate callback
            cbStyles(obj, hc, gdx);
            
            %Arrangements handled by separate callback
            cbArrangements(obj, hc);
            
        end
        
        function cbVariable(obj, hc, ~, fld, gdx, var)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Make a note of new value
                obj.(fld){gdx} = var;
                
                %Force an update
                update(obj)
                
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
        
        function onExport(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it on
                export(obj)
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function tab = tabulate(obj, varargin)
            
            %Results to be exported
            res = obj.Results(obj.Selection);
            
            %Placeholder for results sand messages
            dat = {};
            hdg = {};
            nam = {};
            msg = {};
            
            %For each selected result-set
            for i = 1:numel(res)
                
                %For each element in grid
                for j = 1:obj.Rows*obj.Cols
                    
                    %Extract x and y data, with care
                    [x, y, msg{end+1}] = getVariables(res(i), obj.VariablesX{j}, obj.VariablesY{j});
                    
                    %Problem ?
                    if ~isempty(msg{end})
                        
                        %Move on
                        continue;
                        
                    end
                    
                    %Add to store
                    dat(1:size(x,1), end + (1:size(x,2))) = num2cell(x);
                    dat(1:size(y,1), end + (1:size(y,2))) = num2cell(y);
                    
                    %Headings
                    hdg{end+1} = obj.VariablesX{j};
                    if size(y,2) == 1
                        hdg{end+1} = obj.VariablesY{j};
                    else
                        hdg(end+1:end+size(y,2)) = arrayfun(@(x)[obj.VariablesY{j}, ' (', num2str(x), ')'], 1:size(y,2), 'UniformOutput', false);
                    end
                    
                    %Names
                    nam{end+1} = res(i).Name;
                    nam(end+1:end+size(y,2)) = cell(1,size(y,2));
                    
                end
                
            end
            
            %Splurge
            tab = [nam; hdg; dat];
            
        end
        
        function cbTabSelectionChange(obj, hc, ~)
            
            %Make a note
            obj.CurrentTab = hc.Selection;
            
        end
        
    end
    
    methods
        
        function export_mat(obj, varargin)
            
            %Tabulate data
            tab = tabulate(obj); %#ok<NASGU>
            
            %Save to file
            save(obj.ExportFullFile, 'tab');
            
        end
        
        function [b, msg] = export_xls(obj, varargin)
            
            %Tabulate data
            tab = tabulate(obj);
            
            %Write to Excel
            [b, msg] = obj.xlswrite(obj.ExportFullFile, tab, obj.ExportSheet);
            
            %If caller did not want anything back
            if nargout == 0
                
                %OK ?
                if b
                    
                    %Show us
                    winopen(obj.ExportFullFile);
                    
                else
                    
                    %What went wrong ?
                    error(msg);
                    
                end
                
            end
            
        end
        
        function [b, msg] = export_csv(obj, varargin)
            
            %Tabulate data
            tab = tabulate(obj);
            
            %Write to CSV using generic writer (TODO: fine-tune for Airbus requirements)
            [b, msg] = obj.csvwrite(obj.ExportFullFile, tab);
            
            %If caller did not want anything back
            if nargout == 0
                
                %OK ?
                if b
                    
                    %Show us
                    winopen(obj.ExportFullFile);
                    
                else
                    
                    %What went wrong ?
                    error(msg);
                    
                end
                
            end
            
        end
        
    end
    
end
