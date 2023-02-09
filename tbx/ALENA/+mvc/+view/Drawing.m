classdef Drawing < mvc.view.Container
    
    properties
        
        %The model to which we are attached
        Model;
        
        %The field to display in drawing
        Field = 'Name';
        
        %The selection (which should be part of the model, and highlighted when drawn)
        Selection;
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %Control whether we respond to externally applied changes in selection
        SelectionLinked = true;
        
    end
    
    properties (Access = protected)
        
        %Selection popup
        hSelection;
        
        %SelectionLinked checkbox
        hSelectionLinked;
        
        %The panel in which we display the actual drawing
        hPanel;
        
        %Placeholder for listeners
        Listeners;
        
    end
    
    properties (Dependent, Access = protected)
        
        %List of drawable objects within the model (for internal use)
        Drawables;
        
        %Handle to the 'matlab.graphics.axis.Axes' object inside the panel
        hAxes
        
        %Absolute position of the axes in pixles
        AxesPixelPosition
        
    end
    
    methods % set / get
        
        function val = get.Selection(obj)
            
            %If control is not drawn
            if isempty(obj.hSelection)
                
                %Get directly
                val = obj.Selection;
                
            else
                
                %Get from uicontrol
                val = get(obj.hSelection, 'Value');
                
                %Return the actual selection, rather than the index of the selection
                if val <= numel(obj.Drawables)
                    val = obj.Drawables(val);
                else
                    val = [];
                end
                
            end
            
        end
        
        function set.Selection(obj, val)
            
            %Interested in responding to externally applied change in selection ?
            if ~obj.SelectionLinked %#ok<MCSUP>
                
                %No
                return;
                
            end
            
            %Make a note
            obj.Selection = val;
            
            %If control is drawn
            if ~isempty(obj.hSelection) %#ok<MCSUP>
                
                %No selection is valid from user's point of view
                if isempty(val)
                    
                    %But we need something to put into popup
                    idx = 1;
                    
                elseif isnumeric(val)
                    
                    %TODO: ensure within range ?
                    idx = val;
                    
                else
                    
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
                    mdl = obj.Drawables; %#ok<MCSUP>
                    b = arrayfun(@(x)handle(val) == handle(mdl(x)), 1:numel(mdl));
                    
                    %If we have a selection
                    if any(b)
                        
                        %Make a note
                        idx = find(b);
                        
                    else
                        
                        %No change
                        idx = get(obj.hSelection, 'Value');
                        
                    end
                    
                end
                
                %Assign in uicontrol
                set(obj.hSelection, 'Value', idx); %#ok<MCSUP>
                
            end
            
            %Update content
            update(obj);
            
        end
        
        function val = get.Drawables(obj)
            
            %Anything ?
            if isempty(obj.Model)
                
                %No
                val = [];
                
            else
                
                %Flatten the list
                val = obj.Model.flatlist(2);
                
                %Reduce to just those that are drawable
                if ~isa(val, 'mvc.mixin.Drawable')
                    
                    %Down-select accordingly
                    b = arrayfun(@(x)isa(x, 'mvc.mixin.Drawable'), val);
                    val = val(b);
                    
                end
                
                %You might think we'd want to return only drawables that are also flagged as Visible,
                % but not quite that simple - have to return items that are either visible themselves,
                % of have visible descendents
                val(~[val.VisibleDescendants]) = [];
                
            end
            
        end
        
        function val = get.hAxes(obj)
            %get.hAxes Retrieves the handle to any axes contained within
            %the panel excluding the special '-CoordSys' axes containing
            %the linked coordinate system.
            
            val = findobj(obj.hPanel, 'Type', 'axes', '-not', 'Tag', '-CoordSys');
            
        end
        
        function val = get.Model(obj)
            
            %Get the underlying value
            val = obj.Model;
            
        end
        
        function set.Model(obj, val)
            
            %Make a note
            obj.Model = val;
            
            %TODO: Why not just Update content ??
            %update(obj);
            
            %Redraw content
            if numel(obj.Model) == 1
                draw(obj.Model, obj.hPanel); %#ok<MCSUP>
            end
            
        end
        
    end
    
    methods % construction
        
        function obj = Drawing(model, varargin)
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container(varargin{:});
            
            %Set the default UI mode
            obj.DefaultUIMode = 'rotate';
            
            %Update the context
            obj.addContext(2, '|Set default mouse behaviour>On mouse click>Do nothing' , {@setDefaultUIMode, []});
            obj.addContext(2, '|Set default mouse behaviour>On mouse click>Rotate'     , {@setDefaultUIMode, 'Rotate'});
            obj.addContext(2, '|Set default mouse behaviour>On mouse click>Pan'        , {@setDefaultUIMode, 'Pan'});
            obj.addContext(2, '|Set default mouse behaviour>On mouse scroll>Do nothing', {@setDefaultUIMode, []});
            obj.addContext(2, '|Set default mouse behaviour>On mouse scroll>Zoom'      , {@setDefaultUIMode, 'Zoom'});
            
            %Initialise panel
            obj.hPanel = uiextras.Panel('Parent', obj.UIContainer, ...
                'Padding', 6, 'HitTest', 'off');
            
            %Send hittest upstream
            obj.hPanel.Parent.HitTest = 'off';
            
            %Store the model
            obj.Model = model;
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
        end
        
    end
    
    methods (Access = protected) % updating the view
        
        function [L, ha] = update(obj)
            
            %             %Redraw content
            %             if numel(obj.Model) == 1
            %                 draw(obj.Model, obj.hPanel);
            %             end
            
            %Get the drawable content
            L = obj.Drawables;
            
            %If selection listbox is drawn
            if ~isempty(obj.hSelection)
                
                %Get content to be displayed
                str = {L.(obj.Field)};
                
                %If the underlying list is actually hierarchical
                if isa(L, 'mvc.mixin.Collectable')
                    
                    %Then indent each member of string by depth of corresponding element in tree structure
                    n = [L.DepthInTree];
                    str = arrayfun(@(i)[repmat(' ', 1, n(i)), str{i}], 1:numel(str), 'UniformOutput', false);
                    
                end
                
                %Get current selection
                sel = get(obj.hSelection, 'Value');
                
                %In case multiple selection has been specified (from external source)
                sel(2:end) = [];
                
                %Ensure selection popup remains valid
                sel = min(sel, numel(str));
                
                %In case of no selection
                if isempty(sel) || sel == 0
                    set(obj.hSelection, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
                else
                    set(obj.hSelection, 'String', str, 'Value', sel, 'Enable', 'on');
                end
                
                %Linked option
                set(obj.hSelectionLinked, 'Value', obj.SelectionLinked);
                
            end
            
            %Apply this state only to content being hosted by this view
            ha = findobj(obj.hPanel, 'Type', 'Axes');
            if ~isempty(ha)
                
                if numel(ha) > 1
                    idx = ismember({ha.Tag}, '-CoordSys');
                    ha  = ha(~idx);
                end
                
                %Find selection in list of drawables
                bSel = ismember(L, obj.Selection);
                
                %Pass it on
                setProperty(L(~bSel), 'Selected', 'off', ha);
                setProperty(L(bSel), 'Selected', 'on', ha);
                
                %Good title ?
                if isempty(obj.Selection)
                    tit = '[no selection]';
                else
                    tit = dlgtitle(obj.Selection, true);
                end
                title(ha, tit, 'Interpreter', 'none');
                
            end
            
        end
        
    end
    
    methods (Access = protected) % specific callbacks for model & selection changed
        
        function onModelChanged(obj, ~, ~ )
            
            %Might take some time
            clu = mvc.mixin.UiTools.pointer(obj.hPanel); %#ok<NASGU>
            
            %Careful
            try
                
                %Redraw view
                draw(obj.Model, obj.hPanel);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
            end
            
        end
        
        function onSelectionChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
            %Anything else to do ?
            if ~obj.SelectionLinked
                
                %No
                
            elseif ~isempty(obj.SelectionChangeFcn)
                
                %Pass it on
                obj.SelectionChangeFcn(obj.Selection);
                
            end
            
        end
        
        function onSelectionLinked(obj, ~, ~)
            
            %Make a note
            obj.SelectionLinked = boolean(get(obj.hSelectionLinked, 'Value'));
            
        end
        
    end
    
    methods (Access = protected) % view specific callbacks for mouse interaction
        
        function cbViewButtonDownFcn(obj, hFig, evt)
            %cbViewButtonDownFcn Callback for when a mouse button is
            %pressed on the 'Drawing' view. 
            %
            % Actions:
            %   1. Get the figure manager
            %   2. Turn off the current mode.
            
            %Add to the context structure based on the current mode
            
            %Pass it on
            cbViewButtonDownFcn@mvc.view.Container(obj, hFig, evt);
            
            %Set the menu appearence to reflect the current ui mode
            menus = findall(obj.UIContextMenu, 'Type', 'uimenu');
            if isempty(menus) %Escape route
                return
            end
            idx   = strcmpi({menus.Label}, obj.DefaultUIMode);
            if any(idx)
                menus(idx).Checked = 'on';
            end
            
            %Turn off the current mode
            obj.switchOffUIMode(hFig);
            
        end
        
        function cbViewButtonUpFcn(obj, hFig, evt)
            %cbViewButtonUpFcn Callback for when a mouse button is released
            %on the 'Drawing' view. 
            %
            % Actions:
            %   1. Turn off the current mode.
            
            %Turn off the current mode
            obj.switchOffUIMode(hFig)
            
        end
        
        function cbViewMouseMoveFcn(obj, hFig, evt)
            %cbViewMouseMoveFcn Callback for when the mouse moves over the
            %'Drawing' view.
            %
            % Actions:
            %   1. Query whether there is a uimode set to be active when
            %      the mouse hovers over the axes.
            %   2. Check if we are over the axes
            %   3. If we are over the axes then set the uimode to the users
            %      preferred mode BUT only do this if the axes contains
            %      hgobjects that are appropriate for rotation/panning/etc.
            
            %Has the mouse behaviour been set?
            if isempty(obj.DefaultUIMode)
                return
            end
            
            %Get current point in figure units
            curr_units = evt.Point;
            set(hFig, 'CurrentPoint', curr_units);
            
            %Check to see if we are over an axes
            hAx = obj.localFindAxes(hFig, evt);
            
            %Get the 'ModeManager'
            hManager = uigetmodemanager(hFig);
            
            %Get the function handle for the 'uimode'
            if isempty(hManager.CurrentMode)
                %What is the default UI mode?
                switch lower(obj.DefaultUIMode)
                    case 'rotate'
                        modeFcn = @rotate3d;
                    case 'pan'
                        modeFcn = @pan;
                    case 'zoom'
                        modeFcn = @zoom;
                    otherwise
                        return
                end  
            else
                %Has a UI mode been instigated outside of the view
                %environment - i.e. From the parent figure toolstrip...
                switch hManager.CurrentMode.Name
                    case 'Exploration.Rotate3d'
                        modeFcn = @rotate3d;
                    case 'Exploration.Zoom'
                        modeFcn = @zoom;
                    case 'Exploration.Pan'
                        modeFcn = @pan;
                    otherwise
                        return
                end
            end
            
            %Perform actions based on whether we are over the axes
            if ~isempty(hAx)
                %Does the axes have any children that are suitable for
                %rotating, panning, zooming, etc.
                %   - That is, if the axes only has text then do nothing!
                ch = findobj(hAx.Children, '-not', 'Type', 'text');
                if isempty(ch) %Escape route
                    return
                end
                %Select any 'Transform' objects
                idx = arrayfun(@(g) isa(g, 'matlab.graphics.primitive.Transform'), ch);
                if any(idx)  %Any transforms?
                    idx = arrayfun(@(t) ~isempty(t.Children), ch(idx));
                    if ~any(idx) %Any populated transforms?
                        return
                    end
                end                
                %Initialise the mode object
                h = modeFcn(hFig);
                %Enable the mode
                set(h, 'Enable', 'on');
                %Set behaviour at the start & end of the mouse action
%                 h.ActionPreCallback  = @obj.cbViewButtonDownFcn;
                h.ActionPostCallback = @obj.cbViewButtonUpFcn;
                %Allow the mode to be interrupted
                %  - https://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
                [hManager.WindowListenerHandles.Enabled] = deal(false);  
            else
                %Turn off the current uimodes
                h = modeFcn(hFig);
                set(h, 'Enable', 'off');
                %Restore the pointer
                setptr(hFig, 'arrow');
            end
            
        end
        
    end
    
    methods (Static, Hidden) % helper functions for managing view content
        
        function ax = localFindAxes(fig, evd)
            % Return the axes that the mouse is currently over
            % Return empty if no axes found (i.e. axes has hidden handle)
            %
            % This function is stripped from the 'rotate3D' code.
            
            if ~any(ishghandle(fig))
                return;
            end
            
            % Return all axes under the current mouse point
            allAxes = matlab.graphics.interaction.internal.hitAxes(fig, evd);
            ax = mvc.view.Drawing.locValidateAxes(allAxes);
            
        end
        
        function ax = locValidateAxes(allAxes)
            %locValidateAxes Validates the axes handles in 'allAxes'.
            %
            % This function is stripped from the 'rotate3D' code.
            
            ax = [];
            
            for i=1:length(allAxes)
                candidate_ax=allAxes(i);
                if strcmpi(get(candidate_ax,'HandleVisibility'),'off')
                    % ignore this axes
                    continue;
                end
                b = hggetbehavior(candidate_ax,'Rotate3d','-peek');
                % b can be an object (MCOS) or a handle (UDD)
                if  isobject(b) && ~get(b,'Enable')
                    % ignore this axes
                    
                    % 'NonDataObject' & 'unrotatable' are legacy flags
                elseif ~isappdata(candidate_ax,'unrotatable') ...
                        && ~isappdata(candidate_ax,'NonDataObject')
                    ax = candidate_ax;
                    break;
                end
            end
        end
        
        function switchOffUIMode(hFig)
            %switchOffUIMode Finds the 'mode-manager' of the figure 'hFig'
            %and switches off any 'uimodes' that are active.
            %
            % Actions:
            %   1. Get the figure manager
            %   2. Allow the current mode to be interrupted
            %   3. Turn off the current mode.
            
            %Get the 'ModeManager'
            hManager = uigetmodemanager(hFig);
            
            %Is a mode active?
            if ~isempty(hManager.CurrentMode)
                %If a current mode is active then prevent from blocking
                %callbacks
                hManager.CurrentMode.Blocking = false;
                %Turn off the current mode
                hManager.CurrentMode = [];
                %Quit out
                return
            end
            
        end
        
    end
    
end
