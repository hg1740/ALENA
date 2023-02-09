classdef LayeredDrawing < mvc.view.Drawing & mvc.mixin.Importable & mvc.mixin.Exportable
    %
    % Extends Drawing with the concept of "layers".
    %
    %  Graphical content tagged with a "layer" name can be selectively shown / hidden.

    properties

        %Draw what ?
        DrawOption;
        
        %Draw which layers ?
        DrawLayers = {};
        
    end

    properties (SetAccess = protected)
        
        %Draw what, from the following
        DrawOptions = {'All', 'Selection and children', 'Selection only'};
        
        %Groups of layers
        LayerGroups;
        
    end
    
    properties (Access = protected, Dependent)
        
        %Handy helper
        DrawOptionValue;
        
    end

    properties (Access = protected, Transient)

        %An hbox, to allow for additional controls to be added
        hHbox;

        %A grid by which layers checkboxes are parented (and dynamically updated as required)
        hGrid;
        
        %Allow user control over what is drawn
        hDrawOptions;
        hDrawLayers;

        %Fine-control over layers
        FirstTimeThrough = true;
        
    end
    
    methods % set / get
    
        function val = get.DrawOptionValue(obj)
            
            %Look for match
            val = find(strcmp(obj.DrawOption, obj.DrawOptions));
            
        end
        
        function set.DrawOption(obj, val)
            
            %If numeric passed in
            if isnumeric(val)
                
                %Treat as index
                val = obj.DrawOptions{val}; %#ok<MCSUP>
                
            end
            
            %Make a note
            obj.DrawOption = val;
            
            %Force an update
            update(obj);
            
        end
        
        function val = get.DrawOption(obj)
            
            %Start here
            val = obj.DrawOption;
            
            %Nothing set ?
            if isempty(val)
                
                %Go with first on list
                val = obj.DrawOptions{1};
                
            end
            
        end
    
        function set.DrawLayers(obj, val)
            
            %Make a note
            obj.DrawLayers = val;
            
            %Update content
            update(obj);            
            
        end
        
    end
    
    methods % construction
        
        function obj = LayeredDrawing( varargin )
            
            %Call superclass constructor
            obj@mvc.view.Drawing(varargin{:});
             
            %Initialise list of Layer Groups
            obj.LayerGroups = struct('Name', {}, 'Layers', {});
            
            %Icons placed in subfolder of classdef
            dn = fullfile(fileparts(which(class(obj))), 'ico');
                        
            %Add toolbar buttons
            %   - Toggle all drawable objects 
            addToggleTool(obj, '-DisplayAll', ...
                {@obj.cbSetDrawLayers, '-all'}  , ...
                fullfile(dn, 'aircraft_16.png'), ...
                'View all aircraft elements');
            %   - Toggle lifting surface view
            addToggleTool(obj, 'DisplayAeroPanels', ...
                @obj.cbSetDrawLayers, fullfile(dn, 'lifting_surfaces_16.png'), ...
                'View all lifting surfaces', 'UserData', {'Aero Panel'});
            %   - Toggle control surface view
            addToggleTool(obj, 'DisplayAeroDevices', ...
                @obj.cbSetDrawLayers, fullfile(dn, 'devices_16.png'), ...
                'View all aero devices', ...
                'UserData', {'Ailerons', 'Slats', 'Flaps', 'Spoilers', ...
                'Generic Control Surfaces'});
            %   - Toggle all bluff-body cross-sections
            addToggleTool(obj, 'DisplayCrossSection', ...
                @obj.cbSetDrawLayers, fullfile(dn, 'cross-sections_16.png'), ...
                'View all bluff-bodies', 'UserData', {'Cross-Section'});
            
            %Insert VBox into hierarchy
            hv = uiextras.VBox('Parent', obj.UIContainer, 'HitTest', 'off');
            
            %Send hittest upstream
            hv.Parent.HitTest = 'off';
            
            %Containing a grid
            hg = uiextras.Grid('Parent', hv);
            
            %Containing a label, and popup (to control selection)
            uicontrol('Parent', hg, 'Style', 'Text', ...
                'String'           , 'Selection:')
            obj.hSelection = uicontrol('Parent', hg, 'Style', 'popup', ...
                'String', {'-no selection-'}, 'Value', 1, 'Callback', @obj.onSelectionChange);
            hg.Widths = [75, -1];
            
            %Add a checkbox for the 'SelectionLinked' flag
            obj.hSelectionLinked = uicontrol('Parent', hg, ...
                'Style', 'Checkbox', ...
                'String', 'linked', ...
                'TooltipString', 'Selection in this View is linked to selection in Application', ...
                'Callback', @obj.onSelectionLinked);
            
            %Containing an HBox, which reparents the original panel
            obj.hPanel.Parent = uiextras.HBoxFlex('Parent', hv, ...
                'Padding', 6, 'Spacing', 6);
            
            %Set heights
            hv.Heights = [25, -1];
            
            %Somewhere to put additional controls
            obj.hGrid = uiextras.Grid('Parent', obj.hPanel.Parent);
            obj.hPanel.Parent.Widths(end) = 100;
            
            %Draw option ?
            uicontrol('Parent', obj.hGrid, ...
                'Style' , 'Text'         , ...
                'String', [10, 'Draw...'], ...
                'HorizontalAlignment', 'left');
            uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'popup', ...
                'String'  , obj.DrawOptions, ...
                'Value'   , obj.DrawOptionValue, ...
                'Callback', @obj.onDrawOption, ...
                'HorizontalAlignment', 'left');  
            
            %Other handy options
            uicontrol('Parent', obj.hGrid, ...
                'Style', 'Text', ...
                'String', [10, 'Axes...'], ...
                'HorizontalAlignment', 'left');
            uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'popup', ...
                'String'  , {'normal', 'equal', 'tight'}, ...
                'Value'   , 1    , ...
                'UserData', @axis, ...
                'Callback', @obj.onAxesOption);
            uicontrol('Parent', obj.hGrid, ...
                'Style'  , 'checkbox', ...
                'String' , 'Visible' , ...
                'Callback', @obj.onAxesOption);
            uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'checkbox', ...
                'String'  , 'Box'     , ...
                'Callback', @obj.onAxesOption);
            uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'checkbox', ...
                'String'  , 'Grid'    , ...
                'UserData', @grid     , ...
                'Callback', @obj.onAxesOption);
            uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'checkbox', ...
                'String'  , 'Coordinate System', ...
                'Callback', @obj.cbToggleCoordSys);
            
            %Placeholder for layers (one listbox better than individual checkboxes)
            uicontrol('Parent', obj.hGrid, ...
                'Style', 'Text', ...
                'String', [10, 'Layers...'], ...
                'HorizontalAlignment', 'left');
            obj.hDrawLayers = uicontrol('Parent', obj.hGrid, ...
                'Style', 'Listbox', ...
                'Min', 0, 'Max', 2, ...
                'String', {}, ...
                'TooltipString', 'Selected layers will be drawn in the view', ...
                'Callback', @obj.onDrawLayers);
            
            %Arranged in a column
            obj.hGrid.Widths = -1;
            obj.hGrid.Heights = repmat(25, size(obj.hGrid.Heights));
            obj.hGrid.Heights(end) = 150;
            
            %Update list of available layers
            updateLayerListbox(obj);
            
            %Access to "groups" of layers
            hcm = uicontextmenu(ancestor(obj.hGrid, 'figure'), 'callback', {@obj.onDrawLayers, 'callback'});
            uimenu(hcm, 'Label', 'All', 'Callback', {@obj.onDrawLayers, 'all'});
            uimenu(hcm, 'Label', 'Groups');
            hm = uimenu(hcm, 'Label', 'Manage', 'Separator', 'on');
            uimenu(hm, 'Label', 'Save'     , 'Callback', {@obj.onDrawLayers, 'save'});
            uimenu(hm, 'Label', 'New...'   , 'Callback', {@obj.onDrawLayers, 'new'});
            uimenu(hm, 'Label', 'Remove...', 'Callback', {@obj.onDrawLayers, 'remove'});
            uimenu(hm, 'Label', 'Import...', 'Callback', {@obj.onDrawLayers, 'import'}, 'Separator', 'on');
            uimenu(hm, 'Label', 'Export...', 'Callback', {@obj.onDrawLayers, 'export'});
            set(obj.hDrawLayers, 'UiContextMenu', hcm);
            
            %Initialise Importable/Exportable properties
            obj.ExportSheetTag = [];
            obj.ExportSpec = {'*.gaf', 'Layer Group Arrangement File (*.gaf)', @export_gaf};
            obj.ImportSpec = {'*.gaf', 'Layer Group Arrangement File (*.gaf)', @import_gaf};
            obj.ImportShowProgress = false;
            
            %And look for default list of Layer Groups
            obj.onDrawLayers([], [], 'default');
            
        end
        
    end
    
    methods (Access = protected) % import/export mask files
        
        function export_gaf(obj)
            
            %Save details
            LG = obj.LayerGroups;
            save(obj.ExportFullFile, '-mat', 'LG');
            
        end
        
        function import_gaf(obj)
            
            %Load details
            TEMP = load(obj.ImportFullFile, '-mat', 'LG');
            
            %Make a note
            obj.LayerGroups = TEMP.LG;
            
        end
        
    end
    
    methods (Access = protected) % callbacks
        
        function updateLayerListbox(obj)
            %
            % Alternative approach to checkboxes for each layer
            %  is to add a single listbox, perhaps more scaleable
            
            %What layers have we got ?
            if isempty(obj.hPanel.Children)
                
                %None
                layers = {};
                
            else
                
                %Get from tagged content
                layers = unique(get(findobj(obj.hPanel.Children.Children), 'Tag'));
                layers(cellfun(@isempty, layers)) = [];
                
            end
            
            %Get rid of the '-CoordSys' token
            layers(ismember(layers, '-CoordSys')) = [];            
            
            %Update listbox contents
            set(obj.hDrawLayers, ...
                'String', layers, ...
                'Value' , 1:numel(layers));
            
            %If first time through
            if obj.FirstTimeThrough && ~isempty(layers)
                
                %Make a note
                obj.DrawLayers = layers;
                
                %And lower the first-time flag
                obj.FirstTimeThrough = false;
                
            end
            
        end
                    
        function [L, ha] = update(obj)
            
            %Start with base-class
            [L, ha] = update@mvc.view.Drawing(obj);
           
            %Draw what ?
            switch lower(obj.DrawOption)
                
                case 'all'
                    
                    %Everything is visible
                    bVis = true(size(L));
                    setProperty(L, 'Visible', 'on', ha);
                
                case 'selection only'
                    
                    %Find selection in list of drawables
                    bVis = ismember(L, obj.Selection);
                    
                case 'selection and children'
                    
                    %Find selection in list of drawables
                    bVis = ismember(L, obj.Selection.findall);
                    
                otherwise
                    error('bad draw option');
            end
            
            %Apply
            setProperty(L(~bVis), 'Visible', 'off', ha);
            setProperty(L(bVis), 'Visible', 'on', ha);
            
            %Find content that is currently visible
            hc = findobj(ha, 'Visible', 'on');
             
            %Make the view interactive (so you can click on a
            % graphic feature and select the corresponding object in tree)
            set(hc, 'ButtonDownFcn', @obj.onButtonDown);
            
            %Get tagged content
            tag = get(hc, {'Tag'});
            
            %Only interested if tag non-empty
            b = cellfun(@isempty, tag);
            tag(b) = [];
            hc(b) = [];
            
            %Identify tags that do NOT match what we're looking for
            b = ~ismember(tag, obj.DrawLayers);
            
            %Hide that content
            set(hc(b), 'Visible', 'off');
            
        end
        
        function onButtonDown(obj, hc, ~)
            
            %Might take some time
            clu = mvc.mixin.UiTools.pointer(obj.hPanel); %#ok<NASGU>
            
            %Careful
            try
                
                %Work upstream till we find a transform
                ht = ancestor(hc, 'HGTransform');
                
                %Whose user-data is set to a collectable
                ud = get(ht, 'UserData');
                if isa(ud, 'mvc.mixin.Collectable')
                    
                    %Select it. Unfortunately can not just rely on this
                    % obj.Selection = ud;
                    % because of the way the "SelectionLinked" flag works
                    % to, if set false, effectively ignore changes in selection from
                    % what it thinks is the outside world (but in this
                    % case the change is coming from a sub-classed instance
                    % rather than the outside world, but the SelectionLinked
                    % flag doesn't know that).
                    %
                    %Do it this way instead, so the view manager can't tell
                    % that the user didn't just change the selection in popup
                    [b, sel] = ismember(ud, obj.Drawables);
                    if b
                        
                        %Pass it on
                        set(obj.hSelection, 'Value', sel);
                        onSelectionChange(obj);
        
                        %Nice touch
                        context(ud, hc);
                        
                    end
                    
                end
                
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
            end
            
        end
        
        function onModelChanged(obj, ~, ~ )
            
            %Might take some time
            clu = mvc.mixin.UiTools.pointer(obj.hPanel); %#ok<NASGU>
            
            %Careful
            try
                
                %Redraw view
                draw(obj.Model, obj.hPanel);
                
                %Refresh the list of layers presented in the view
                %obj.DrawLayers = addLayerCheckboxes(obj);
                updateLayerListbox(obj);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
            end
            
        end
                
        function onDrawLayers(obj, hc, ~, action, nam)
            
            %Get selected layers, by name, from uicontrol
            str = obj.hDrawLayers.String;
            val = obj.hDrawLayers.Value; % {obj.hDrawLayers(logical([obj.hDrawLayers.Value])).String};
            str = str(val);
            
            %Get rid of the '-CoordSys' token
            str(ismember(str, '-CoordSys')) = [];            
            obj.DrawLayers = str;
            
            %Any special action required ?
            if nargin > 3
            
                %Careful
                try
                    
                    %Do what ?
                    switch action
                        
                        case 'export'
                            
                            %Pass it on
                            export(obj);
                            
                            %Nothing else required
                            return;
                
                        case 'default'
                            
                            %Go with this file
                            fn = fullfile(fileparts(which(class(obj.Model))), 'defaults.gaf');
                            
                            %If found
                            if exist(fn, 'file') == 2
                                
                                %Import it
                                import(obj, fn);
                                
                            end
                            
                            %Nothing else required
                            return;

                        case 'import'
                            
                            %Pass it on
                            import(obj);
                            
                            %Nothing else required
                            return;
                            
                        case 'remove'
                            
                            %Remove what ?
                            sel = listdlg(obj, {obj.LayerGroups.Name}, 'Remove group(s)...');
                            
                            %Cancelled ?
                            if isempty(sel)
                                return;
                            end
                            
                            %Do it
                            obj.LayerGroups(sel) = [];
                            
                            %Nothing else required
                            return;
                            
                        case 'new'
                            
                            %Name of new group
                            nam = inputdlg(obj, matlab.lang.makeUniqueStrings('Group', {obj.LayerGroups.Name}));
                            
                            %Cancelled ?
                            if isempty(nam)
                                return;
                            end
                            
                            %Make a note
                            obj.LayerGroups(end+1).Name = nam;
                            obj.LayerGroups(end).Layers = str;
                            
                            %Nothing else to do
                            return;
                            
                        case 'all'
                            
                            %Caller is asking for all layers to be drawn
                            str = obj.hDrawLayers.String;
                            
                            %Update UI too
                            obj.hDrawLayers.Value = 1:numel(str);
                            
                        case 'callback'
                            
                            %Locate 'Groups' uimenu
                            hm = findobj(hc, 'Label', 'Groups');
                            
                            %Clean sheet
                            delete(hm.Children);
                            
                            %Populate uimenu with list of available Layer Groups
                            nam = {obj.LayerGroups.Name};
                            cellfun(@(x)uimenu(hm, 'Label', x, 'Callback', {@obj.onDrawLayers, 'apply', x}), nam);
                            
                            %And return
                            return;
                            
                        case 'apply'
                            
                            %Look for match with list of known Layer Groups
                            idx = find(strcmp(nam, {obj.LayerGroups.Name}));
                            assert(~isempty(idx), ['Group ''', nam, ''' not found']);
                            
                            %Hence layers in this group
                            str = obj.LayerGroups(idx).Layers;
                            
                            %Update UI too
                            [b, val] = ismember(str, obj.hDrawLayers.String);
                            obj.hDrawLayers.Value = val(b);
                            
                        otherwise
                            
                            %Not acceptable
                            error('bad action');
                            
                    end
                
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj, err));
                    
                end
                
            end
            
            %Should any of the toggle buttons be updated?
            %   - TODO: Is this the right place to update this?
            tag  = {obj.hUIToggleTool.UserData};
            idx  = false(size(tag));
            idx_ = cellfun(@isempty, tag);
            idx(~idx_) = cellfun(@(x) all(ismember(x, str)), tag(~idx_));
            
            %Update it!
            set(obj.hUIToggleTool(idx) , 'State', 'on');
            set(obj.hUIToggleTool(~idx), 'State', 'off');
            
            %Pass it on
            obj.DrawLayers = str;
            
        end
        
        function onDrawOption(obj, hc, ~)
            
            %Get value from control and pass it on
            obj.DrawOption = hc.Value;
            
        end
        
        function onAxesOption(obj, hc, ~)
            
            %Find axes
            ha = obj.hAxes;

            %What have we got ?
            switch get(hc, 'Style')
                
                case 'popupmenu'
                    
                    %Value ?
                    val = popupstr(hc);
                    
                case 'checkbox'
                    
                    %What property are we setting ?
                    prp = get(hc, 'String');
                    
                    %Value ?
                    val = mvc.mixin.UiTools.bool2offon(get(hc, 'Value'));
            
                otherwise
                    
                    %Not (yet) handled
                    error(['uicontrol style ''', get(hc, 'Style'), ''' not (yet) handled']);
                    
            end
            
            %Get the control user-data
            ud = get(hc, 'UserData');
            if isempty(ud)
            
                %No - just apply the property / value pair to the axes
                set(ha, prp, val);
    
            elseif isa(ud, 'function_handle')
            
                %Call it
                arrayfun(@(ax) ud(ax, val), ha);

            else
                error(['user data of class ''', class(ud), ''' not (yet) handled']);
            end
            
        end
        
        function cbToggleCoordSys(obj, hc, ~)
            
            %Don't toggle the value if there is no data visible for the
            %model
            if isempty(obj.Model) || ~obj.Model.VisibleDescendants
                return
            end
            
            %Grab axes and their respective tags
            ha  = obj.hAxes;
            idx = ismember({ha.Tag}, '-CoordSys');
            
            if nnz(idx) == 0 %Escape route
                return
            end
            
            %Grab children            
            ch  = ha(idx).Children;
            
            %Are they visible
            vis = {ch.Visible};
            tok = unique(vis);
            if numel(tok) > 1
                error(['Disambiguous graphics object state. Expected ', ...
                    'all the objects to have their ''Visible'' state ', ...
                    'to ''off'' or ''on''.']);
            end
            
            switch tok{1}
                case 'on'
                    set(ch, 'Visible', 'off');
                    hc.Value = 0;
                case 'off'
                    set(ch, 'Visible', 'on');
                    hc.Value = 1;
            end
            
        end
        
        function cbSetDrawLayers(obj, ~, ~, varargin)
            %cbSetDrawLayers Sets the current value of 'DrawLayers' based 
            %on whether the 'uipushtools' in the toolbar are active or not.
            %Also, sets the selection in the 'hDrawLayers' uicontrol object 
            %to reflect this.
            
            %Look for special tokens
            tok = '-all';
            idx = arrayfun(@(x) ismember(x, tok), varargin);
            
            %What push tools are currently 'on'
            %   - Discount the '-DisplayAll' push button
            hpt   = findobj(obj.hUIToggleTool, '-not', 'Tag', '-DisplayAll', 'State', 'on');
            
            %Get the handle of the special '-DisplayAll' push tool
            hpt_all = findobj(obj.hUIToggleTool, 'Tag', '-DisplayAll');
            
            %Choose next action based on what is currently selected
            if isempty(hpt) || any(idx)
                %Draw everything - all other push tools are set to 'off'
                hpt = findobj(obj.hUIToggleTool, '-not', 'Tag', '-DisplayAll');
                set(hpt, 'State', 'off');
                %Draw everything
                hpt_all.State = 'on';
                obj.hDrawLayers.Value = 1 : numel(obj.hDrawLayers.String);
                obj.DrawLayers = obj.hDrawLayers.String;
                return
            else
                hpt_all.State = 'off';
            end
            
            %Layer names are stored as user data in the hgobjects
            layers = {hpt.UserData};
            idx    = cellfun(@(x) iscellstr(x), layers); % only cellstr
            layers = unique(horzcat(layers{idx}));
            
            if isempty(layers) %Escape route
                return
            end
            
            %Update the selection
            allLayers = obj.hDrawLayers.String;
            index  = cellfun(@(x) find(ismember(allLayers, x), 1), layers, 'Unif', false);
            obj.hDrawLayers.Value = vertcat(index{:});
            
            %Pass it on            
            obj.DrawLayers = layers;
            
        end
        
    end
    
end
