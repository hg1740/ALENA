classdef SimulationResult < mvc.view.Container
    
    properties (AbortSet, SetObservable)
        
        %Laid out how ?
        Layout = [1, 1];
        
        %Which variable(s) to plot ?
        VariablesX = {};
        VariablesY = {};
        VariablesZ = {};
        
        %What style(s) of plot ?
        Style = {};
        
        %What filter(s) to apply ?
        Filter = {};
        
    end
    
    properties (Constant)
        
        Layouts = {[1, 1], [2, 1], [3, 1], [4, 1], [1, 2]};
        Styles = {'timeseries', 'waterfall', 'image', 'scatter', 'histogram', 'barh'};
        
    end
    
    properties
        
        %The model to which we are attached
        Model;
        
        %The field to display
        Field = 'Name';
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %Track currently selected item in model
        % (required if model is a collection, because this view can present only one item)
        Selection;
        
    end
    
    properties (Access = protected)
        
        %Selection popup
        hSelection;
        
        %The panel in which we display results
        hPanel;
        
        %Placeholder for context menu
        hContextMenu;
        
        %Placeholder for listeners
        Listeners;
        
    end
    
    methods % get/set
        
        function set.Layout(obj, val)
            
            %Make a note
            obj.Layout = val;
            
            %Update
            update(obj);
            
        end
        
        function set.VariablesX(obj, val)
            
            %Make a note
            obj.VariablesX = val;
            
            %Update
            update(obj);
            
        end
        
        function set.VariablesY(obj, val)
            
            %Make a note
            obj.VariablesY = val;
            
            %Update
            update(obj);
            
        end
        
        function set.VariablesZ(obj, val)
            
            %Make a note
            obj.VariablesZ = val;
            
            %Update
            update(obj);
            
        end
        
        function set.Style(obj, val)
            
            %Make a note
            obj.Style = val;
            
            %Update
            update(obj);
            
        end
        
        function val = get.Selection(obj)
            
            %Get from uicontrol
            val = get(obj.hSelection, 'Value');
            
            %Return the actual selection, rather than the index of the selection
            val = obj.Model(val);
            
        end
        
        function set.Selection(obj, val)
            
            %Make a note
            obj.Selection = val;
            
            %No selection is valid from user's point of view
            if isempty(val)
                
                %But we need something to put into popup
                idx = 1;
                
            elseif isnumeric(val)
                
                %TODO: ensure within range ?
                idx = val;
                
            elseif isa(val, class(obj.Model)) %#ok<MCSUP>
                
                %Allow selection to be set with an actual member of the list
                %
                %As long as it is in the list
                % Would like to do this, but it fails
                %[b, idx] = ismember(val, obj.Model);
                %assert(all(b), 'invalid selection');
                %
                %This doesn't work either, under all circumstances, because of heterogeneous collection
                %idx = arrayfun(@(x)find(x == obj.Model), val);
                %assert(all(idx > 0), 'invalid selection');
                
                %Instead
                mdl = obj.Model; %#ok<MCSUP>
                b = arrayfun(@(x)handle(val) == handle(mdl(x)), 1:numel(mdl));
                
                %TODO: this test not quite adequate ?
                assert(any(b), 'invalid selection');
                
                %Make a note
                idx = find(b);
                
            end
            
            %Assign in uicontrol
            set(obj.hSelection, 'Value', idx); %#ok<MCSUP>
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.Model(obj)
            
            %Get the underlying value
            val = obj.Model;
            
            %If Collectable
            if isa(val, 'mvc.mixin.Collectable')
                
                %Flatten it
                val = flatlist(val, 2);
                
            end
            
        end
        
        function set.Model(obj, val)
            
            %Make a note
            obj.Model = val;
            
            %Update content
            update(obj);
            
        end
        
    end
    
    methods % construction
        
        function obj = SimulationResult( model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.SimulationResult;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Start with a VBox
            hv = uiextras.VBox('Parent', obj.UIContainer);
            
            %Containing a grid
            hg = uiextras.Grid('Parent', hv);
            
            %Containing a label, and popup (to control selection)
            uicontrol('Parent', hg, 'Style', 'Text', 'String', 'Selection:')
            obj.hSelection = uicontrol('Parent', hg, 'Style', 'popup', 'callback', @obj.onSelectionChange);
            hg.Widths = [75, -1];
            
            %And a panel
            obj.hPanel = uiextras.Panel('Parent', hv, 'Padding', 6);
            
            %Set heights
            hv.Heights = [25, -1];
            
            %Store the model
            obj.Model = model;
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
            %             %Context menus
            %             hcm = uicontextmenu('Parent', ancestor(obj.hPanel, 'Figure'));
            %             uimenu(hcm, 'Label', 'Layout...', 'Callback', {@obj.onCallback, 'layout'});
            %             uimenu(hcm, 'Label', 'Variables...', 'Callback', {@obj.onCallback, 'variables'});
            %             uimenu(hcm, 'Label', 'Filter...', 'Callback', {@obj.onCallback, 'filter'});
            %             uimenu(hcm, 'Label', 'Style...', 'Callback', {@obj.onCallback, 'style'});
            %             set(obj.hPanel, 'UiContextMenu', hcm);
            %             obj.hContextMenu = hcm;
            
        end
        
        function delete(obj)
            
            %Suicide pact
            delete(obj.hPanel);
            
        end
        
    end
    
    methods
        
        function layout(obj, varargin)
            
            %Options ?
            str = cellfun(@mat2str, obj.Layouts, 'UniformOutput', false);
            
            %Current selection ?
            sel = find(cellfun(@(x)isequal(obj.Layout, x), obj.Layouts));
            
            %Offer to user
            sel = listdlg(obj, str, 'Choose layout...', ...
                'SelectionMode', 'single', ...
                'InitialValue', sel);
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Unpick
            obj.Layout = obj.Layouts{sel};
            
        end
        
        function style(obj, n)
            
            %What have we got so far ?
            [b, sel] = ismember(obj.getStyle(n), obj.Styles);
            
            %What have we got at the moment ?
            sel = listdlg(obj, obj.Styles, 'Choose style...', ...
                'SelectionMode', 'single', ...
                'InitialValue', sel(b));
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Make a note
            obj.setStyle(n, sel);
            
        end
        
        function hilite(obj, n, varargin)
            
            %What is available ?
            var = obj.Selection.Variables;
            
            %What have we got so far ?
            sel = obj.getVariables(n, varargin{:});
            sel(sel > numel(var)) = [];
            
            %If more than one
            if numel(sel) > 1
                
                %Hilite what, exactly ?
                subsel = listdlg(obj, var(sel), 'Choose variable(s)...', ...
                    'ListSize', [600, 300], ...
                    'SelectionMode', 'single');
                
                %Cancelled ?
                if isempty(subsel)
                    return;
                end
                
                %Down-select accordingly
                sel = sel(subsel);
                
            end
            
            %Find the corresponding element in model
            nam = var{sel};
            tok = strsplit(nam, '.');
            el = obj.Selection.getElement(tok{1});
            pth = el.BlockPath.convertToCell;
            blk = pth{1};
            
            %How tedious is this ?
            for i = 2:numel(tok)
                
                %Search at this level
                blk2 = find_system(blk, 'findall', 'on', 'LookUnderMasks', 'on', 'Type', 'line', 'Name', tok{i});
                
                %Anything ?
                if numel(blk2) == 1
                    
                    %Yes
                    blk = blk2;
                    
                end
                
            end
            
            %And show us
            hilite_system(blk);
            
        end
        
        function variables(obj, n, varargin)
            
            %What is available ?
            var = obj.Selection.Variables;
            
            %What have we got so far ?
            sel = obj.getVariables(n, varargin{:});
            sel(sel > numel(var)) = [];
            
            %Offer to user
            sel = listdlg(obj, var, 'Choose variable(s)...', ...
                'ListSize', [600, 300], ...
                'InitialValue', sel);
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Make a note
            obj.setVariables(n, varargin{:}, sel);
            
        end
        
        function filter(obj, n)
            
        end
        
    end
    
    methods (Access = private)
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onSelectionChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
            %Anything else to do ?
            if ~isempty(obj.SelectionChangeFcn)
                
                %Pass it on
                obj.SelectionChangeFcn(obj.Selection);
                
            end
            
        end
        
        function update(obj)
            
            %Get the model
            mdl = obj.Model;
            
            %Anything to do ?
            if isempty(mdl)
                
                %Bail
                return;
                
            end
            
            %Get content to be displayed
            str = {mdl.(obj.Field)};
            
            %If the underlying list is actually hierarchical
            if isa(mdl, 'mvc.mixin.Collectable')
                
                %Then indent each member of string by depth of corresponding element in tree structure
                n = [mdl.DepthInTree];
                str = arrayfun(@(i)[repmat(' ', 1, n(i)), str{i}], 1:numel(str), 'UniformOutput', false);
                
            end
            
            %Get current selection
            sel = get(obj.hSelection, 'Value');
            
            %In case multiple selection has been specified (from external source)
            sel(2:end) = [];
            
            %Ensure selection popup remains valid
            sel = min(sel, numel(str));
            
            %In case of no selection
            if isempty(sel)
                set(obj.hSelection, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
            else
                set(obj.hSelection, 'String', str, 'Value', sel, 'Enable', 'on');
            end
            
            %Does selected item support context menus ?
            if isa(mdl(sel), 'mvc.mixin.Contextable')
                
                %Yes - pass it on
                context(mdl(sel), obj.hSelection);
                
            end
            
            %Redraw view of selected signals
            delete(allchild(obj.hPanel)); % expensive ?
            
            %Can we draw anything ?
            if ~any(cellfun(@(x)isa(obj.Selection, x), {'mvc.model.SimulationResult', 'mvc.model.SimulationSweep'}))
                
                %No
                return;
                
            end
            
            %Construct view(s) in a grid
            hg = uiextras.GridFlex('Parent', obj.hPanel, 'Padding', 6, 'Spacing', 6);
            
            %Containing axes (with intermediate container so legend etc don't mess up)
            ha = arrayfun(@(x)axes('Parent', uicontainer('Parent', hg)), ...
                1:prod(obj.Layout));
            
            %Set layout
            hg.Widths = -1 .* ones(1, obj.Layout(2));
            
            %For each element of grid
            for i = 1:prod(obj.Layout)
                
                %Careful
                try
                    
                    %Pass it on
                    plot(obj, i, ha(i));
                    
                catch err
                    
                    %Display what went wrong in the axes
                    text(ha(i), 0, 0.5, ['FAILED with error:', 10, 10, err.message], ...
                        'HorizontalAlignment', 'left');
                    
                end
                
                %Update context menus
                context(obj, ha(i), i);
                
            end
            
            %Sort out any labelling requirement
            obj.axesLabels(hg);
            
            %Ensure labels are adjusted as required, on size change
            set(hg, 'SizeChangedFcn', @obj.axesLabels);
            
            %And on zoom
            
        end
        
    end
    
    methods % helpers
        
        function setStyle(obj, n, val)
            
            %Valid ?
            if ischar(val)
                assert(ismember(val, obj.Styles), 'invalid style');
                val = {val};
            elseif iscellstr(val)
                assert(all(ismember(val, obj.Styles)), 'invalid style');
            elseif isnumeric(val)
                assert(val >= 1 && val <= numel(obj.Styles), 'invalid style');
                val = obj.Styles(val);
            end
            
            %Assign in cell
            obj.Style(n) = val;
            
        end
        
        function val = getStyle(obj, n)
            
            %Enough ?
            if numel(obj.Style) < n || isempty(obj.Style{n})
                
                %No - go with this by default
                val = obj.Styles{1};
                
            else
                
                %Get it
                val = obj.Style{n};
                
            end
            
        end
        
        function setVariables(obj, n, typ, varargin)
            
            %What variable type ?
            if isnumeric(typ)
                
                %Include typ in varargin
                varargin = [{typ}, varargin];
                
                %Type not specified - go with Y-axis by default
                typ ='y';
                
            end
            
            %Field of interest
            foi = ['Variables', upper(typ)];
            
            %Assign in cell
            obj.(foi)(n) = varargin;
            
        end
        
        function val = getVariables(obj, n, typ)
            
            %What variable type ?
            if nargin < 3
                
                %Not specified - go with Y-axis by default
                typ ='y';
                
            end
            
            %Field of interest
            foi = ['Variables', upper(typ)];
            
            %Enough ?
            if numel(obj.(foi)) < n
                
                %No - go with this by default
                val = 1;
                
            else
                
                %Get it
                val = obj.(foi){n};
                
            end
            
        end
        
        function [val, nam, lab] = getPlotData(obj, n, varargin)
            
            %Corresponding variables ?
            var = obj.getVariables(n, varargin{:});
            
            %Get the data
            [val, nam] = obj.Selection.getVariable(var);
            
            %Get labels, if required
            if nargout > 2
                lab = obj.Selection.getLabels(val, nam);
            end
            
        end
        
        function ho = plot(obj, n, ha, varargin)
            
            %In what style ?
            style = obj.getStyle(n);
            
            %Get the data to be plotted
            [dat, tit, lab] = getPlotData(obj, n);
            
            %Offset for multiple elements
            dy = 0;
            
            %Placeholder for label detail
            labels.mode = 'none';
            labels.labels = [];
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d, l{i}] = obj.flatten(dat(i), lab{i});
                
                %Suitable y-values
                y = 1:size(d,2);
                
                %Handle according to style
                switch style
                    
                    case 'timeseries'
                        
                        %Plot as line(s)
                        ho{i} = plot(ha, t, d, ...
                            'HitTest', 'off', ... % so right click falls through
                            varargin{:});
                        
                        %Apply labels as ?
                        labels.mode = 'Legend';
                        
                    case 'waterfall'
                        
                        %Adding an offset to enuerated types is tricky
                        if isa(d, 'Simulink.IntEnumType')
                            d = int32(d);
                        end
                        
                        %And plot it - NOT using standard 'waterfall', which is not (I think) what we want.
                        % Instead, work out a suitable offset between rows
                        rng = max(range(d, 1));
                        
                        %Build an array of offsets
                        offset = rng .* ((1:size(d,2)).' * ones(1,size(d,1))).';
                        %offset = rng .* (cumsum(ones(1,size(d,2), 'like', d)).' * ones(1,size(d,1), 'like', d)).';
                        
                        %And the just plot as normal
                        ho{i} = plot(ha, t, d + offset + dy, ...
                            'HitTest', 'off', ... % so right click falls through
                            varargin{:});
                        
                        %Best if y-axes point downwards, as they would with image
                        set(ha, 'YDir', 'rev');
                        
                        %Offset for next
                        dy = dy + max(offset(:));
                        
                        %Apply labels as ?
                        labels.mode = 'YTickLabel';
                        
                    case 'image'
                        
                        %Plot as image
                        ho{i} = imagesc(ha, t, y + dy, d.', ...
                            'HitTest', 'off', ... % so right click falls through
                            varargin{:});
                        
                        %Offset for next
                        dy = dy + y(end);
                        
                        %Apply labels as ?
                        labels.mode = 'YTickLabel';
                        
                    otherwise
                        error(['unknown style ''', style, '''']);
                end
                
                %Fill the available real-estate
                axis(ha, 'tight');
                
                %Hold for next
                hold(ha, 'on');
                
            end
            
            %Why is this necessary ?
            %             if dy > 0
            %                 ylim(ha, [0 dy] + 0.5);
            %             end
            
            %Labels
            title(ha, tit, 'Interpreter', 'none');
            
            %Make a note of the label details, the parent SizeChangeFcn
            % will take responsibility for sorting it all out
            labels.labels = [l{:}];
            setappdata(ha, 'Labels', labels);
            
        end
        
        function scatter(obj, n, ha, varargin)
            
            %Get the data to be plotted (x-axis)
            [datx, titx] = getPlotData(obj, n, 'x');
            [tx, dx] = obj.flatten(datx);
            
            %Get the data to be plotted (y-axis)
            [dat, tit] = getPlotData(obj, n, 'y');
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d] = obj.flatten(dat(i));
                
                %And plot it
                scatter(ha, dx, d, varargin{:});
                
                %Hold it
                hold(ha, 'on');
                
            end
            
            %Labels
            xlabel(ha, titx, 'Interpreter', 'none');
            ylabel(ha, tit, 'Interpreter', 'none');
            
        end
        
        function timeseries(obj, n, ha, varargin)
            
            %Get the data to be plotted
            [dat, tit] = getPlotData(obj, n);
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d] = obj.flatten(dat(i));
                
                %And plot it
                plot(ha, t, d, varargin{:});
                
                %Hold it
                hold(ha, 'on');
                
            end
            
            %Labels
            title(ha, tit, 'Interpreter', 'none');
            
        end
        
        function waterfall(obj, n, ha, varargin)
            
            %Offset for multiple elements
            dy = 0;
            
            %Get the data to be plotted
            [dat, tit] = getPlotData(obj, n);
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d] = obj.flatten(dat(i));
                
                %And plot it - NOT using standard 'waterfall', which is not (I think) what we want.
                % Instead, work out a suitable offset between rows
                rng = max(range(d, 1));
                
                %Build an array of offsets
                offset = rng .* ((1:size(d,2)).' * ones(1,size(d,1))).';
                
                %And the just plot as normal
                plot(ha, t, d + offset + dy, varargin{:});
                
                %Hold for next
                hold(ha, 'on');
                
                %Offset for next
                dy = dy + max(offset(:));
                
            end
            
            %Labels
            title(ha, tit, 'Interpreter', 'none');
            
        end
        
        function image(obj, n, ha, varargin)
            
            %Get the data to be plotted
            [dat, tit, lab] = getPlotData(obj, n);
            
            %Offset for multiple elements
            dy = 0;
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d, l{i}] = obj.flatten(dat(i), lab{i});
                
                %Suitable y-values
                y = 1:size(d,2);
                
                %And plot it
                hm(i) = imagesc(ha, t, y + dy, d.', ...
                    'HitTest', 'off', ... % so right click falls through
                    varargin{:});
                
                %Hold for next
                hold(ha, 'on');
                
                %Offset for next
                dy = dy + y(end);
                
            end
            
            %Why is this necessary ?
            ylim(ha, [0 dy] + 0.5);
            
            %Labels
            title(ha, tit, 'Interpreter', 'none');
            
            %For tick labels, can't just do this
            % set(ha, 'YTickLabels', lab);
            %
            %Because we can't assume that the axes scaling would show every tick,
            % so intead just make a note of the label details, the parent SizeChangeFcn
            % will take responsibility for sorting it all out
            setappdata(ha, 'YTickLabel', [l{:}]);
            
        end
        
        function barh(obj, n, ha, varargin)
            
            %Get the data to be plotted
            [dat, tit, lab] = getPlotData(obj, n);
            
            %Offset for multiple elements
            dy = 0;
            
            %For each in turn
            for i = 1:numel(dat)
                
                %Ensure plottable
                [t, d, l{i}] = obj.flatten(dat(i), lab{i});
                
                %Suitable y-values
                y = 1:size(d,2);
                
                %Intervals between time points ?
                dt = diff(t);
                
                %Hence up-time
                d = sum(dt .* d(1:end-1,:), 1);
                
                %And plot it
                hm(i) = barh(y + dy, d, ...
                    'HitTest', 'off', ... % so right click falls through
                    'Parent', ha, ...
                    varargin{:});
                
                %Hold for next
                hold(ha, 'on');
                
                %Offset for next
                dy = dy + y(end);
                
            end
            
            %Why is this necessary ?
            ylim(ha, [0 dy] + 0.5);
            
            %Labels
            title(ha, tit, 'Interpreter', 'none');
            
            %For tick labels, can't just do this
            % set(ha, 'YTickLabels', lab);
            %
            %Because we can't assume that the axes scaling would show every tick,
            % so intead just make a note of the label details, the parent SizeChangeFcn
            % will take responsibility for sorting it all out
            setappdata(ha, 'YTickLabel', [l{:}]);
            
        end
        
        function context(obj, hh, n)
            
            %Get containing figure
            hf = unique(arrayfun(@(x)ancestor(x, 'Figure'), hh));
            assert(numel(hf) == 1, 'non-unique figure handles');
            
            %Apply standard set of context menus to handle(s) h
            hcm = uicontextmenu('Parent', hf);
            uimenu(hcm, 'Label', 'Layout...', 'Callback', {@obj.onCallback, 'layout'});
            
            %Options for choice of variables is style-dependent
            switch obj.getStyle(n)
                
                case 'scatter'
                    
                    %Need separate variables for x and y (and z ?)
                    hmm = uimenu(hcm, 'Label', 'Variables');
                    uimenu(hmm, 'Label', 'X', 'Callback', {@obj.onCallback, 'variables', n, 'x'});
                    uimenu(hmm, 'Label', 'Y', 'Callback', {@obj.onCallback, 'variables', n, 'y'});
                    uimenu(hmm, 'Label', 'Z', 'Callback', {@obj.onCallback, 'variables', n, 'z'});
                    
                otherwise
                    
                    %Just need dependent variable
                    uimenu(hcm, 'Label', 'Variables...', 'Callback', {@obj.onCallback, 'variables', n});
                    
            end
            
            %Filter and ability to change style common to all
            uimenu(hcm, 'Label', 'Filter...', 'Callback', {@obj.onCallback, 'filter', n});
            uimenu(hcm, 'Label', 'Style...', 'Callback', {@obj.onCallback, 'style', n});
            uimenu(hcm, 'Label', 'Hilite...', 'Callback', {@obj.onCallback, 'hilite', n});
            set(hh, 'UiContextMenu', hcm);
            
        end
        
    end
    
    methods (Access = private, Static)
        
        function axesLabels(hp, ~, ~)
            
            %For all axes contained by parent
            ha = findobj(hp, 'Type', 'Axes');
            for i = 1:numel(ha)
                
                %Look for relevant content
                if ~isappdata(ha(i), 'Labels')
                    continue;
                end
                
                %Get the details
                labels = getappdata(ha(i), 'Labels');
                
                %How are they applied ?
                switch lower(labels.mode)
                    
                    case 'none'
                        
                        %Do nothing
                        
                    case 'legend'
                        
                        %How many entries is reasonable to show ?
                        if numel(labels.labels) <= 10
                            
                            %Display as legend
                            legend(ha(i), labels.labels);
                        
                        end
                        
                    case 'yticklabel'
                        
                        %Do not simply display all - better to match against tick marks
                        yt = get(ha(i), 'YTick');
                        str = cell(size(yt));
                        
                        %Only interested in integers
                        idx = find(rem(yt, 1) == 0 & yt > 0 & yt <= numel(labels.labels));
                        str(idx) = labels.labels(yt(idx));
                        
                        %Apply
                        set(ha(i), 'YTickLabel', str);
                        
                    otherwise
                        error(['bad label mode ''', labels.mode, '''']);
                end
                
            end
            
        end
        
        function [t, d, l] = flatten(ts, lab)
            
            %Unpick the time-series
            t = ts.Time;
            d = ts.Data;
            
            %Ensure the dimensions are ok
            sz = size(d);
            switch numel(sz)
                
                case 3
                    
                    %Permute accordingly
                    d = permute(d, [3 1 2]);
                    
                    %And flatten anything left
                    d = reshape(d, [sz(3), prod(sz(1:2))]);
                    
                    %Ditto labels, if required
                    if nargin > 1 && nargout > 2
                        
                        %Have we got anything ?
                        if all(cellfun(@isempty, lab))
                            
                            %No
                            l = {};
                            
                        else
                            
                            %TODO: check dimensions are consistent
                            
                            %Map the labels
                            ff = fullfact(sz(1:2));
                            l = arrayfun(@(x){lab{1}{ff(x,1)}, lab{2}{ff(x,2)}}, ...
                                1:size(ff,1), 'UniformOutput', false);
                            l = cellfun(@(x)x(~cellfun(@isempty, x)), l, 'UniformOutput', false);
                            l = cellfun(@(x)strjoin(x, ', '), l, 'UniformOutput', false);
                            
                        end
                        
                    end
                    
                case {1, 2}
                    
                    %Nothing to do to data
                    
                    %Build list of labels, if required
                    if nargin > 1 && nargout > 2
                        l = lab{1};
                    end
                    
                otherwise
                    error('unsupport data size');
            end
            
        end
        
    end
    
end
