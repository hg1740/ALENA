classdef Container < matlab.mixin.SetGet ...
        & matlab.mixin.Heterogeneous ...
        & mvc.mixin.UiTools ...
        & mvc.mixin.Contextable
    %Container  Peer to uicontainer. Allows user to reparent and reposition
    %the container, while ensuring that the lifecycles of the objects are
    %tied together.
    %
    %   Copyright 2013 The MathWorks, Inc.
    %   $Date: 2016-06-17 15:40:31 +0100 (Fri, 17 Jun 2016) $
    
    properties %Title
        %Title of the container
        Title; 
    end
    
    properties( Access = protected )
        %Handle to the underlying 'uicontainer'
        UIContainer         
        %Default 'uimode' for when the mouse hovers over the view.
        DefaultUIMode = [];
        %Current 'uimode'
        CurrentUIMode = [];        
    end
    
    properties (SetAccess = private, GetAccess = protected) % ui-tools
        %Handle to any user-specified toolbar buttons including
        %'uipushtool', 'uitoggletool' and 'uisplittool'.
        hToolBarButtons
        %Structure containing data about the toolbar buttons
        %   - Name  : Name of the tool
        %   - Type  : Type of gobject
        %   - Fcn   : Callback function
        %   - ToolTipString : String to be displayed when mouse hovers over
        %   - Index : Order of the tool in the toolbar
        %   - Group : Name of the group that this tool belongs too
        %   - GroupIndex    : Order of the tool within the group
        %   - GroupPriority : Ordering of the groups
        ToolBarButtonData = struct('Name', {}, 'Type', {}, ...
            'Fcn', {}, 'ToolTipString', {}, 'Index', {}, 'Group', {}, ...
            'GroupIndex', {}, 'GroupPriority', {}, 'Extras', {});
        ToolBarButtonGroups = {};
        %Handles to any user-specified 'uipushtool' objects
        hUIPushTool
        %Handles to any user-specified 'uitoggletool' objects
        hUIToggleTool
    end
    
    properties(Dependent)
        %Parent of the container
        Parent
        %Units of measure for the position of the container
        Units 
        %Inner-position of the container
        Position 
        %Logical flag for indicating whether the container is visible
        Visible 
        %Handle to the underlying context menu
        UIContextMenu
    end
    
    properties (Dependent, Hidden = true)
        %Handle to the parent figure
        hFigure_
        %Handle to the toolbar of the parent figure
        hToolBar_
    end
    
    methods % construction / destruction
        
        function obj = Container( varargin )
            %Container  Peer to uicontainer.
            
            % Create container
            hContainer = uicontainer( varargin{:}, ...
                'DeleteFcn', @obj.onDeleted );
            
            % Set properties
            obj.UIContainer = hContainer;
            
            %Update the context
            obj.addContext(1, 'Help', 'openHelp')
            
        end % constructor
        
        function delete( obj )
            %delete  Destructor

                delete( obj.UIContainer )
            
        end % destructor
        
    end 
    
    methods % set / get
        
        function val = get.Title(obj)
        
            %What have got ?
            val = obj.Title;
            
            %Nothing yet ?
            if isempty(val)
                
                %Go with classname
                tok = strsplit(class(obj), '.');
                val = tok{end};
                
            end
            
        end
        
        function value = get.Parent( obj )
            
            value = get( obj.UIContainer, 'Parent' );
            
        end % get.Parent
        
        function set.Parent( obj, value )
            
            set( obj.UIContainer, 'Parent', value );
            
        end % set.Parent
        
        function value = get.Units( obj )
            
            value = get( obj.UIContainer, 'Units' );
            
        end % get.Units
        
        function set.Units( obj, value )
            
            set( obj.UIContainer, 'Units', value );
            
        end % set.Units
        
        function value = get.Position( obj )
            
            value = get( obj.UIContainer, 'Position' );
            
        end % get.Position
        
        function set.Position( obj, value )
            
            set( obj.UIContainer, 'Position', value );
            
        end % set.Position
        
        function value = get.Visible( obj )
            
            value = get( obj.UIContainer, 'Visible' );
            
        end % get.Visible
        
        function set.Visible( obj, value )
            
            set( obj.UIContainer, 'Visible', value );
            
        end % set.Visible

        function value = get.UIContextMenu(obj) %get.UIContextMenu
            %get.UIContextMenu Get method for the dependent property
            %'UIContextMenu'.
            
            value = obj.UIContainer.UIContextMenu;
            
        end
        
        function set.DefaultUIMode(obj, value) %set.DefaultUIMode
           %set.DefaultUIMode Set method for the property 'DefaultUIMode'.
           %
           % The 'DefaultUIMode' can either be empty of one of the valid
           % strings:
           %    * 'rotate' - The built-in 'rotate3d' mode.
           %    * 'pan'    - The built-in 'pan' mode.
           %    * 'zoom'   - The built-in 'zoom' mode.
           
           if isempty(value) %Always okay
                obj.DefaultUIMode = [];
                return
           end
           
           %Allowable values for the 'DefaultUIMode'
           tok = {'Rotate', 'Pan', 'Zoom'};
           
           %Validate the entry
           validatestring(lower(value), tok, class(obj), 'DefaultUIMode');

           %Ensure case-sensitive match with 'uimenu' label
           idx = strcmpi(tok, value);
           
           %Set it
           obj.DefaultUIMode = tok{idx};
            
        end
        
        function value = get.hFigure_(obj) %get.hFigure_
           %get.hFigure_ Get method for the property 'hFigure_'.
           
           %Simple!
           value = ancestor(obj.UIContainer, 'matlab.ui.Figure');
           
           if isempty(value) %This should never happen...
               value = [];
           end
            
        end
            
        function value = get.hToolBar_(obj) %get.hToolBar_
            %get.hToolBar_ Get method for the dependet property
            %'hToolBar_'.
            
            %Sensible default
            value = [];
            
            %Invoke the 'get' method once!
            hFig = obj.hFigure_;
            
            %Check if the container is parented by a figure
            if isempty(hFig)
                return
            end
            
            %Make sure all handles are visible
            r   = groot;
            val = r.ShowHiddenHandles;
            r.ShowHiddenHandles = 'on';
            
            %Find it!
            value = findobj(hFig, 'Type', 'uitoolbar');
            
            %Return to original state
            r.ShowHiddenHandles = val;
            
            %If the figure doesn't have on then create one
            if isempty(value)
                value = uitoolbar('Parent', hFig);
            end
            
        end
        
    end 
    
    methods (Access = private)
        
        function onDeleted( obj, ~, ~ )
            
            % Call destructor
            obj.delete()
            
        end % onDeleted
        
    end % event handlers
    
    methods (Access = protected) % view specific callbacks
        
        function cbViewButtonDownFcn(obj, hFig, evt, varargin)
            %cbViewButtonDownFcn Generic callback for handling the case
            %where a mouse-button is pressed on the view. 
            %
            % If the user right-clicks on the view then a context menu will
            % open providing the user with view-specific context options.
            
            %Dependent on mouse button
            switch hFig.SelectionType
                case 'alt'
                    %Open up context menu - Use 'Contextable' for this.
                    context(obj, obj.UIContainer, obj.Context); 
                    %Get the handle to the context menu
                    hcm = obj.UIContextMenu;
                    %Set the position
                    hcm.Position = evt.Point;
                    %Make sure it is visible
                    hcm.Visible = 'on';
                otherwise 
            end
            
        end
        
        function cbViewButtonUpFcn(obj, src, evt)
            %cbViewButtonDownFcn Generic callback for handling the case
            %where a mouse-button is released on the view. 
            %
            % N.B. This is just a placeholder method. It is expected that
            % this method will be over-loaded at the subclass level.
        end
        
        function cbViewMouseMoveFcn(obj, src, evt)
            %cbViewButtonDownFcn Generic callback for handling the case
            %where the mouse is moved on the view.
            %
            % N.B. This is just a placeholder method. It is expected that
            % this method will be over-loaded at the subclass level.
        end
        
    end
    
    methods % callbacks with public access
        
        function openHelp(obj)
            %openHelp Opens the help page for the object 'obj'. 
            %
            % TODO - Expand this so that it checks for user-defined
            % documentation first and accesses that.
            
            %Open up the documentation centre
            eval(['doc ', class(obj)]);
            
        end
        
        function setDefaultUIMode(obj, val)
            %setDefaultUIMode Thin wrapper around set-method for the
            %property 'DefaultUIMode' - allows the property to be set in
            %callbacks.
            
            obj.DefaultUIMode = val;
            
        end
    end
    
    methods (Access = protected) % configuring the figure toolstrip
        
        function varargout = addToolBarButton(obj, htbb, nam, fn, cdata, desc, idx, grp, grp_idx, grp_pr, varargin)
            %addToolBarButton Helper function for adding a 'ui-button' to
            %the figure's toolbar. 
            %
            % The 'ui-button' objects can be arranged in 'groups'. If there
            % are any groups in the 'ToolBarButtonData' then these groups
            % are arranged starting on the left of the toolbar. All other
            % individual buttons are placed to the right of the grouped
            % buttons.
            
            %Grab the toolbar handle 
            htb = obj.hToolBar_;
            
            if isempty(htb) %Escape route
                return
            end
            
            %Do we already have a toolbar button with this name?
            buttons = {obj.ToolBarButtonData.Name};
            if any(ismember(buttons, nam))
                %If we already have a button with this name then quit out
                varargout{1} = [];
                return
            end
            
            %What have we got?
            if ishghandle(htbb)
                %Can start from the object - Ensure correct parent
                set(htbb, 'Parent', htb);
            elseif isa(htbb, 'function_handle')
                %Can start from the function handle
                htbb = htbb(htb);
            elseif ischar(htbb)
                %Can start from a character
                fn   = str2func(htbb);
                htbb = fn(htb);
            end
            assert(ishghandle(htbb), 'Expected the toolbar button to be a handle-graphics object');
                        
            %Function can be specified as a string
            if ischar(fn)
                fn = str2func(fn);
            end
            assert(isa(fn, 'function_handle'), 'Expected the callback function to be a function handle');
            
            %Icon can be specified as a file
            if ischar(cdata)
                if ~exist(cdata, 'file')
                    cdata = [];
                else
                    %Import the image
                    [img, map] = imread(cdata);
                    if isa(img, 'uint8')
                        %Continue
                        cdata = img;
                        %                     %Replace any white space with the default MATLAB colour
                        %                     clr = uint8( obj.hFigure_.Color .* 256);
                        %                     idx = and(cdata(:, :, 1) == 255, range(cdata, 3) == 0);
                        %                     for i = 1 : 3
                        %                         cdata(idx, i) = clr(1);
                        %                     end
                        %                     disp('Wait');
                    else
                        % Convert image from indexed to truecolor
                        if isempty(map)
                            cdata = ind2rgb(img);
                        else
                            cdata = ind2rgb(img, map);
                        end
                    end
                end
            else
                %Check size
                assert(isequal(size(cdata), [16, 16, 3]), 'Expected cdata to be a 16 x 16 RGB array');
            end
            
            %Configure the ordering of the buttons
            if nargin < 7 || isempty(idx) %Index specified?
                idx = numel(obj.ToolBarButtonData) + 1;
            end
            if nargin < 8  %Group specified?
                grp = '';
            end
            if nargin < 9  %Group index specified?
                grp_idx = [];
            end
            if nargin < 10 %Group priority specified?
                %Group already specified?
                idx = ismember({obj.ToolBarButtonData.Group}, grp);
                if any(idx)
                    %Use exisiting group priority
                    ind    = find(idx == true, 1, 'first');
                    grp_pr = obj.ToolBarButtonData(ind).GroupPriority;
                else
                    %Increment the group...
                    nGrp   = numel(unique({obj.ToolBarButtonData.Group}));
                    grp_pr = nGrp + 1;
                end
            end

            %Pull any extras out of varargin
            if isa(htbb, 'uisplittool')
                
            end
            
            %Configure the group 
            %   - If the user has specified a group, use this to infer the
            %     'grp_idx' if it hasn't been specified
            if ~isempty(grp)
                %Check if the group already exists
                idx = ismember({obj.ToolBarButtonData.Group}, grp);
                if any(idx)
                    if isempty(grp_idx) %Group index specified?
                        %Just add it onto the end
                        grp_idx = nnz(idx) + 1;
                    end
                else
                    %If this is the first time the group has been defined
                    %then the index within that group must be 1!
                    grp_idx = 1;
                end
            end
%             
%             %Set the 'GroupPriority'
%             if isempty(grp_pr)
%                 %What groups do we have?
%                 grps = {obj.ToolBarButtonData.Group};
%                 %Just tag it on to the end
%                 grp_pr = numel(unique(grps)) + 1;
%             end
            
            %Update the structure
            obj.ToolBarButtonData(end + 1).Name      = nam;
            obj.ToolBarButtonData(end).Type          = htbb.Type;
            obj.ToolBarButtonData(end).Fcn           = fn;
            obj.ToolBarButtonData(end).ToolTipString = desc;
            obj.ToolBarButtonData(end).Index         = idx;
            obj.ToolBarButtonData(end).Group         = grp;
            obj.ToolBarButtonData(end).GroupIndex    = grp_idx;
            obj.ToolBarButtonData(end).GroupPriority = grp_pr;
            obj.ToolBarButtonData(end).Extras        = varargin;            
            
            %Set all the properties of the button object
            set(htbb, ...
                'Tag'            , nam  , ...
                'ClickedCallback', fn   , ...
                'CData'          , cdata, ...
                'ToolTipString'  , desc);
            
            %Configure the order ...
            grps = {obj.ToolBarButtonData.Group};
            
            %How many in each group?
            [ugrps, idx] = unique(grps);
            grp_ntb = cellfun(@(x) nnz(ismember(grps, x)), ugrps);
            grp_prs = [obj.ToolBarButtonData(idx).GroupPriority];
            
            %Grab any buttons that aren't grouped and make a note of their
            %current order
            idx   = cellfun(@isempty, grps);            
            index = [obj.ToolBarButtonData(idx).Index];    
            
            %Does each group have the same priority?
            grp_pr_i = cellfun(@(x) range( ...
                [obj.ToolBarButtonData(ismember(grps, x)).GroupPriority]), ugrps);
            assert(all(grp_pr_i == 0), ['All tool bar button groups ', ...
                'must have the same priority']);
            
            %Sort the groups and assign index numbers
            [grp_prs, ind] = sort(grp_prs, 'ascend');
            grp_ntb        = grp_ntb(ind);
            
            %Update the java container of the MATLAB toolbar to refresh
            %contents
            
            %User wants it back?
            if nargout > 0
                varargout{1} = htbb;
            end
            
        end
        
        function varargout = addToggleTool(obj, nam, fn, cdata, desc, varargin)
            %addToggleTool Helper function for adding a 'uitoggletool' 
            %object to the figure toolbar.
            
            if nargin < 4
                cdata = zeros(16, 16, 3);
            end
            if nargin < 5
                desc = '';
            end
            
            %Pull out any user-data
            ind = find(ismember(varargin(1 : 2 : end), 'UserData'));
            if isempty(ind)
                ud =[];
            else
                ud = varargin{ind * 2};
            end
            
            %Function can be specified as a string
            if ischar(fn) 
                fn = str2func(fn);
            end
            
            %Icon can be specified as a file
            if ischar(cdata)
                if ~exist(cdata, 'file')
                    cdata = [];
                else
                %Import the image
                [img, map] = imread(cdata);
                if isa(img, 'uint8')
                    %Continue
                    cdata = img;
%                     %Replace any white space with the default MATLAB colour
%                     clr = uint8( obj.hFigure_.Color .* 256);
%                     idx = and(cdata(:, :, 1) == 255, range(cdata, 3) == 0);
%                     for i = 1 : 3
%                         cdata(idx, i) = clr(1);
%                     end
%                     disp('Wait');
                else                    
                    % Convert image from indexed to truecolor
                    if isempty(map)
                        cdata = ind2rgb(img);
                    else
                        cdata = ind2rgb(img,map);
                    end
                end
                end
            end
            
            %Grab the toolbar handle 
            htb = obj.hToolBar_;
            
            %Make the tool
            ht = uitoggletool(htb, ....
                'Tag'            , nam  , ... % meaningful name
                'ClickedCallback', fn   , ... % callback function
                'CData'          , cdata, ... % icon data
                'TooltipString'  , desc , ... % helpful text
                'UserData'       , ud);       % useful data
            
            %Store a reference in this object
            if isempty(obj.hUIToggleTool)
                obj.hUIToggleTool = ht;
            else
                obj.hUIToggleTool(end + 1) = ht;
            end
            
            %Give it back?
            if nargout == 1
                varargout{1} = ht;
            end            
            
        end
        
        function varargout = addPushTool(obj, nam, fn, cdata, desc, varargin)
            %addToggleTool Helper function for adding a 'uitoggletool' 
            %object to the figure toolbar.
            
            if nargin < 4
                cdata = zeros(16, 16, 3);
            end
            if nargin < 5
                desc = '';
            end
            
            %Pull out any user-data
            ind = find(ismember(varargin(1 : 2 : end), 'UserData'));
            if isempty(ind)
                ud =[];
            else
                ud = varargin{ind * 2};
            end
            
            %Function can be specified as a string
            if ischar(fn) 
                fn = str2func(fn);
            end
            
            %Icon can be specified as a file
            if ischar(cdata)
                if ~exist(cdata, 'file')
                    cdata = [];
                else
                %Import the image
                [img, map] = imread(cdata);
                if isa(img, 'uint8')
                    %Continue
                    cdata = img;
%                     %Replace any white space with the default MATLAB colour
%                     clr = uint8( obj.hFigure_.Color .* 256);
%                     idx = and(cdata(:, :, 1) == 255, range(cdata, 3) == 0);
%                     for i = 1 : 3
%                         cdata(idx, i) = clr(1);
%                     end
%                     disp('Wait');
                else                    
                    % Convert image from indexed to truecolor
                    if isempty(map)
                        cdata = ind2rgb(img);
                    else
                        cdata = ind2rgb(img,map);
                    end
                end
                end
            end
            
            %Grab the toolbar handle 
            htb = obj.hToolBar_;
            
            %Make the tool
            ht = uipushtool(htb, ....
                'Tag'            , nam  , ... % meaningful name
                'ClickedCallback', fn   , ... % callback function
                'CData'          , cdata, ... % icon data
                'TooltipString'  , desc , ... % helpful text
                'UserData'       , ud);       % useful data
            
            %Store a reference in this object
            if isempty(obj.hUIPushTool)
                obj.hUIPushTool = ht;
            else
                obj.hUIPushTool(end + 1) = ht;
            end
            
            %Give it back?
            if nargout == 1
                varargout{1} = ht;
            end            
            
        end
        
    end
    
    methods
        
        function varargout = onCallback(obj, hc, evtdata, action, varargin)
            %
            % Generic callback handler, allows specific callbacks to not worry
            % about setting (and restoring) pointer, or try-catch constructs
            
            %May take a while
            clu = obj.pointer; %#ok<NASGU>
            
            %Careful
            try
                
                %Do what ?
                fcn = str2func(action);
                
                %Do it
                [varargout{1:nargout}] = fcn(obj, varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err));
                
            end
            
        end
        
        function str = dlgtitle(obj, varargin)
   
            %Simple
            str = obj.Title;
            
        end
        
    end
    
    methods(Sealed) %helper functions for heterogeneous arrays
        function varargout = set(obj,varargin)
            [varargout{1:nargout}] = set@matlab.mixin.SetGet(obj,varargin{:});
        end
        function varargout = get(obj,varargin)
            [varargout{1:nargout}] = get@matlab.mixin.SetGet(obj,varargin{:});
        end
    end
    
end % classdef