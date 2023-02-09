classdef ViewManager < mvc.view.Container & mvc.mixin.UiTools
    %
    % ViewManager provides a container within which one or more view are managed.
    
    properties
        
        %The model to which we are attached
        Model;
        
        %Track currently selected item in model
        Selection;
        
        %What views are supported by this manager ?
        SupportedViews = {};
        
        %Offer option to use to "browse" for additional view(s) ?
        OfferBrowseForNewView = true;
        
        %What arrangements of views are supported by this manager ?
        SupportedArrangements;
        
        %Which arrangement is currently selected ?
        CurrentArrangement;
                
        %Useful when prompting for input
        GuiName = 'View Manager';
        
    end
        
    properties (SetAccess = private)
        
        %The list of actual views
        Views;
        
        %The position of the views (in pixels) with respect to the figure
        ViewPixelPosition
        
        %What figure do we live inside ?
        hFigure;
        
        %Where do the view(s) live ?
        hTabPanel;
         
        %Where do panel(s) live ?
        hPanels;
        
        %Placeholder for listeners
        Listeners;  
               
        %List of temp files, awaiting deletion
        TempFiles = {};
        
    end
    
    properties (Dependent)
    
        %How many rows of panels have we got, in each column of panels ?
        PanelRowsPerColumn;
        
        %Control over widths of separate panels
        PanelWidths;
        
        %Description of how view(s) are currently arranged
        ViewArrangement;
        
    end
    
    properties (Access = protected) % To help persist view arrangements
        
        ViewArrangementFileMask = {'*.vaf', 'View arrangement files (*.vaf)'; ...
            '*.mat', 'MATLAB data files (*.mat)'; ...
            '*.*', 'All files (*.*)'};
        ViewArrangementFile;

    end
    
    methods % get/set
        
        function val = get.Views(obj)
        
            %Start here
            val = obj.Views;
            
            %Anything ?
            if isempty(val)
                return;
            end
            
            %Check for validity
            b = ~isvalid(val);
            if any(b)
                
                %Throw away any duffers
                val(b) = [];
                
                %And do same for property
                obj.Views(b) = [];
                
            end
            
        end
        
        function set.ViewArrangement(obj, val)
            
            %If caller passing in a filename
            if ischar(val)
            
                %Look for relevant content
                TEMP = load('-mat', val, 'ViewArrangement');
                if isempty(TEMP) || ~isfield(TEMP, 'ViewArrangement') || isempty(TEMP.ViewArrangement)
                    
                    %Nothing to do - return silently
                    return;
                    
                end
                
                %Get the arrangement and carry on
                val = TEMP.ViewArrangement;
                
            end
            
            %Remove all existing views
            obj.removeAllPanels;
 
            %Panel counter
            n = 1;
            
            %For each column
            for i = 1:numel(val.ColumnWidths)
                
                %Add the panel
                obj.addPanel('right');
                
                %For each row
                for j = 1:numel(val.RowHeights{i})
                
                    %For second and subsequent
                    if j > 1
                        
                        %Add the panel
                        obj.addPanel('below');
                        
                    end
                    
                    %Add the view(s)
                    cellfun(@(x)obj.addView(x, [], 'ignore'), val.Views{n}, 'UniformOutput', false);
                    
                    %Next
                    n = n + 1;
                    
                end

                %Set heights (carefully - allowing for the possibility of annotation in view)
                ht = val.RowHeights{i};
                obj.hPanels.Contents(i).Heights(end-numel(ht)+1:end) = ht;
                
            end

            %Set widths
            obj.hPanels.Widths = val.ColumnWidths;
            
            %Make a note of arrangement name
            obj.CurrentArrangement = val.Name;
            
        end
        
        function val = get.ViewArrangement(obj)
        
            %Placeholder for name
            val.Name = [];
            
            %How many panel rows per column ?
            val.PanelRowsPerColumn = obj.PanelRowsPerColumn;
            
            %Column widths ?
            val.ColumnWidths = obj.hPanels.Widths;
            
            %Panel counter
            n = 1;
                
            %Get columns
            hc = obj.hPanels.Contents;
            
            %For each column
            for i = 1:numel(hc)
            
                %Row heights ?
                val.RowHeights{i} = hc(i).Heights;
                
                %Get rows
                hr = hc(i).Contents;
                
                %For each row
                for j = 1:numel(hr)
                    
                    %Allow for possibility of non-tab panel in view
                    if isa(hr(j), 'uiextras.TabPanel')
                        
                        %Get the view(s) by name
                        val.Views{n} = hr(j).TabNames;
                    
                        %Next
                        n = n + 1;
                
                    end
                    
                end
                
            end
            
        end
        
        function val = get.PanelRowsPerColumn(obj)
            %
            %Return a vector, one element for each column of panels,
            % where the value of the element is the number of rows in that column
            val = arrayfun(@(x)numel(x.Children), obj.hPanels.Children);
            
        end

        function set.GuiName(obj, val)
            
            %Make a note
            obj.GuiName = val;
            
            %Update content
            update(obj);
            
        end
        
        function val = get.PanelWidths(obj)
            
            %Pass it on
            val = obj.hPanels.Widths;
            
        end

        function set.PanelWidths(obj, val)
            
            %Pass it on
            obj.hPanels.Widths = val;
            
        end
        
        function set.Selection(obj, val)
            
            %Make a note
            obj.Selection = val;
            
            %Pass it on to any views
            for i = 1:numel(obj.Views) %#ok<MCSUP>
                obj.Views(i).Selection = val; %#ok<MCSUP>
            end
            
            %Update content
            update(obj);
            
        end
        
        function set.Model(obj, val)
            
            %Make a note
            obj.Model = val;
            
            %Pass it on to any views
            for i = 1:numel(obj.Views) %#ok<MCSUP>
                obj.Views(i).Model = val; %#ok<MCSUP>
            end
            
            %Update content
            update(obj);
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = ViewManager( model, varargin)
            
            %No point carrying on without GLT
            v = ver('layout');
            assert(~isempty(v), ['GUI Layout Toolbox (GLT) not found.', 10, ...
                'Please download GLT from MATLAB File Exchange']);
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = [];
                
            end
            
            %If parent not specified in varargin
            if ~any(strcmpi('Parent', varargin(1:2:end)))
                
                %Go with a new figure
                hf = figure('NumberTitle', 'off', ...
                    'IntegerHandle'      , 'off', ...
                    'HandleVisibility'   , 'off');
                
                %Make a note of figure handle, to be passed on to container
                extras = {'Parent', hf};
             
            else
                
                %Needed later
                hf = ancestor(varargin{find(strcmpi('Parent', varargin(1:2:end))) * 2}, 'Figure');
                extras = {};
                
            end

            %Start with base class
            obj@mvc.view.Container(extras{:}, varargin{:});
            %Views are presented in a series of tab pages, displayed in panels of an HBoxFlex
            obj.hPanels = uiextras.HBoxFlex('Parent', obj.UIContainer, 'Padding', 6, 'Spacing', 6);
            
            %Icons placed in subfolder of classdef
            dn = fullfile('+mvc', '+view', 'ico');
%             dn = fullfile(fileparts(which(class(obj))), 'ico');

            %Add toolbar buttons
            if isa(model, 'mvc.mixin.Importable')
               
                hT(1) = addToolBarButton(obj, 'uitoggletool', 'Open', ...
                    @obj.cbImport, zeros(16, 16, 3), ...
                    'Open a model', [], 'File', 1);
                hT(2) = addToolBarButton(obj, 'uitoggletool', 'Save', ...
                    @obj.cbImport, zeros(16, 16, 3), ...
                    'Save a session', [], 'File', 2);                
                hT(3) = addToolBarButton(obj, 'uitoggletool', 'Import', ...
                    @obj.cbImport, fullfile(dn, 'folder-import_16.png'), ...
                    'Import a model', [], 'File', inf);
                %Add import button
                hT = addPushTool(obj, 'Import', @obj.cbImport, ...
                    fullfile(dn, 'folder-import_16.png'), 'Import a model');
                
                %TODO - Implement uitoggletool etc. as generic tool buttons
                %uisplittool, uitogglesplittool, 
%                 %Grab the recent filenames from the model
%                 fNames = model.ImportRecent;
%                 
%                 %Use magic Java dust to add the drop-down menu
%                 %   - https://undocumentedmatlab.com/blog/customizing-standard-figure-toolbar-menubar
%                 %   - TODO : Update this list during the view managers
%                 %     'update' method.
%                 jT = get(hT, 'JavaContainer');
%                 jTMenu = jT.getMenuComponent;
%                 for iF = 1 : numel(fNames)
%                     jMenuItem = handle(jTMenu.add(fNames{iF}),'CallbackProperties');
%                     set(jMenuItem,'ActionPerformedCallback',{@obj.Model.import, fNames{iF}});
%                 end
                                
            end
            
            %Create panel using a helper function (in order that user can subsequently add/remove more panels)
            obj.addPanel(obj.hPanels);

            %What views do we support ?  Depends on the application, but can infer a reasonable starting point
            if isa(model, 'mvc.mixin.Collectable')
                obj.SupportedViews(end+1,:) = {'Tree', @mvc.view.Tree};
                obj.SupportedViews(end+1,:) = {'List', @mvc.view.List};
            end
            if isa(model, 'mvc.mixin.Nameable')
                obj.SupportedViews(end+1,:) = {'Properties', @mvc.view.Properties};
            end
            if isa(model, 'mvc.mixin.Drawable')
                if isa(model, 'mvc.mixin.Collectable')
                    obj.SupportedViews(end+1,:) = {'Drawing', @mvc.view.LayeredDrawing};
                else
                    obj.SupportedViews(end+1,:) = {'Drawing', @mvc.view.Drawing};
                end
            end
            if isa(model, 'mvc.mixin.Simulatable')
                obj.SupportedViews(end+1,:) = {'Simulation Result', @mvc.view.SimulationResult};
            end
            if isa(model, 'mvc.mixin.Auditable')
                obj.SupportedViews(end+1,:) = {'Audit Trail', @mvc.view.AuditTrail};
            end
            
            %If the model is also Contextable
            if isa(model, 'mvc.mixin.Contextable')
                
                %Pass it on
                context(model, hf);
                
            end
            
            %Store the model
            obj.Model = model;

            %Listen for changes
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
            %Initialise the figure window in which the view manager lives
            bArrange = initialiseFigure(obj, hf);
            
            %Can we find a default view arrangements file ?
            arrangeViews(obj, 'default');

            %If caller does not want the manager back
            if nargout == 0 && bArrange
                
                %If content is collectable
                if isa(obj.Model, 'mvc.mixin.Collectable')
                    
                    %It would be reasonable to start with a Tree View on the left
                    obj.addView('Tree');
                    
                    %Tree views don't play well with other types
                    obj.addPanel;
                    obj.PanelWidths = [-1 -4];

                end
                
                %If the content is nameable
                if isa(obj.Model, 'mvc.mixin.Nameable')
                    
                    %And a Properties view
                    obj.addView('Properties');
                    
                end
                
                %If the content is drawable
                if isa(obj.Model, 'mvc.mixin.Drawable')
                    
                    %And a Drawing view
                    obj.addView('Drawing');
                    
                end
                    
            end
            
            %Store a reference to the 'ViewManager' in the figure so that
            %views can be added/removed/etc. at the view level.
            obj.hFigure.UserData = [obj.hFigure.UserData, {obj}];
            
            %Anything required in resize?
            obj.hFigure.ResizeFcn = @obj.cbOnFigResize;
            
            %Assign callbacks for enabling interactive figure 
            %   'ButtonMotionFcn', 'ButtonDownFcn', 'ButtonUpFcn', etc
            set(obj.hFigure, 'WindowButtonMotionFcn', ...
                {@obj.cbViewManagerMouseAction, @cbViewMouseMoveFcn});
            set(obj.hFigure, 'WindowButtonDownFcn', ...
                {@obj.cbViewManagerMouseAction, @cbViewButtonDownFcn});
            set(obj.hFigure, 'WindowButtonUpFcn', ...
                {@obj.cbViewManagerMouseAction, @cbViewButtonUpFcn});
            
        end
        
        function delete(obj)
        
            %Suicide pact ?? NOT sure that a full-blown delete is a good idea
            % (because might have multiple views open of the same underlying model)
            % delete(obj.Model);
            
        end
        
        function close(obj, varargin)

            %If model is Serializable
            if isa(obj.Model, 'mvc.mixin.Serializable')
                
                %Ensure foregrounded
                figure(obj.hFigure);
                
                %Pass it on
                if ~canClose(obj.Model, varargin{:})
                    
                    %Bail out
                    return;
                    
                end
                
                %Make a note of position
                obj.Model.LastPosition = get(obj.hFigure, 'Position');
                
            end    
            
            %OK to actually close
            delete(obj.hFigure);
            
        end
        
    end
    
    methods % Manipulating the View(s)
         
        function hp = addPanel(obj, pos, par, vbox)
            
            %Check if handle to 'VBoxFlex' has been passed in
            if nargin < 4
                vbox = [];
            end
            
            %Unpick inputs
            if nargin == 1
                
                %Default position
                pos = 'right';
                
                %Default parent
                par = obj.UIContainer.Children;
                
            elseif nargin == 2 && ischar(pos)
                
                %Default parent
                if strcmp(pos, 'right')
                    par = obj.UIContainer.Children(1);
                elseif strcmp(pos, 'below')
                    par = obj.UIContainer.Children(1).Children(1);
                else
                    error('bad position: must be either ''right'' or ''below''');
                end

            elseif nargin == 2 && ishghandle(pos)
                
                %Parent is passed in
                par = pos;
                
                %Default position
                pos = 'right';
                
            elseif ischar(pos) && ishghandle(par)
                
                %We're good to go
                
            elseif ishghandle(pos) && ischar(par)
                
                %Wrong way round
                [pos, par] = deal(par, pos);
                
            else
                %Nothing else will do
                error('bad inputs');
            end
            
            %A new panel to the right or left
            if strcmpi(pos, 'right') || strcmpi(pos, 'left')
                
                %Parent is a 'uix.HBoxFlex' - Stash a reference as we are
                %about to over-write it
                hbf = par;
                
                %Means start with a VBoxFlex (so that further panels can be added below this one)
                par = uiextras.VBoxFlex('Parent', par, 'Padding', 6, 'Spacing', 6);
                
                %Need to switch the order of the children to account for
                %'left'/'right' -> If we don't have a reference to the
                %'VBoxFlex' then just append the new panel
                if isa(vbox, 'uix.VBoxFlex')
                    
                    %Keep a list of the 'VBoxFlex' indices
                    index = 1 : numel(hbf.Children);
                                        
                    %Which number 'VBoxFlex' was the context menu invoked
                    %from? N.B. The order is 'reversed', first is last and
                    %last is first.
                    ind   = find(hbf.Children == vbox);
                    ind_  = find(hbf.Children == par);
                   
                    %Remove new panel from the index list (for now...)
                    index(ind_) = [];
                    
                    %Split into left and right views                    
                    switch pos
                        case 'right'
                            r = index(1 : find(index == ind) - 1);
                            l = index(find(index == ind) : end);
                        case 'left'
                            r = index(1 : find(index == ind));
                            l = index(find(index == ind) + 1 : end);
                    end

                    %Insert new panel index
                    index = [r, ind_, l];

                    %Change the order of children
                    hbf.Children = hbf.Children(index);
                    
                end
                
            end
            
            %TODO: cater for user-friendly, helpful message in an empty panel
            %
            %             %Start with an annotation showing some helpful text
            %             ha = annotation(par, 'textbox', ...
            %                 'String', 'Right-click for options', ...
            %                 'HorizontalAlignment', 'center', ...
            %                 'VerticalAlignment', 'middle');
            
            %Create a new tab panel
            hp = uiextras.TabPanel('Parent', par, ...
                'TabWidth', 100, ...
                'Padding', 6, ...
                'SelectionChangedFcn', @i_selectionChange);            
            
            %Add context
            hcm = uicontextmenu('Parent', ancestor(obj.Parent, 'Figure'), 'Callback', {@obj.cbContext, hp});
            set(hp, 'UiContextMenu', hcm);
            %set(ha, 'UiContextMenu', hcm);
            
            % This issues a prompt to select view from list
            %uimenu(hcm, 'Label', 'Add view...', 'Callback', {@obj.cbAddView, [], hp});
            % This adds available views to uimenu on the fly - better I think
            uimenu(hcm, 'Label', 'Add view');
            
            %Manipulate existing views (remove, move, etc)
            uimenu(hcm, 'Label', 'Remove view...', 'Callback', {@obj.cbRemoveView, [], hp});
            uimenu(hcm, 'Label', 'Remove views...', 'Callback', {@obj.cbRemoveView, 'ask', hp});
            uimenu(hcm, 'Label', 'Move view to');
            hmm = uimenu(hcm, 'Label', 'Copy view to');
            uimenu(hmm, 'Label', 'Figure', 'Callback', {@obj.cbCopyView, hp, 'figure'});
            uimenu(hmm, 'Label', 'File', 'Callback', {@obj.cbCopyView, hp, 'file'});
            uimenu(hmm, 'Label', 'Clipboard', 'Callback', {@obj.cbCopyView, hp, 'clipboard'});
            uimenu(hmm, 'Label', 'Word', 'Callback', {@obj.cbCopyView, hp, 'word'});
            uimenu(hmm, 'Label', 'Excel', 'Callback', {@obj.cbCopyView, hp, 'excel'});
            hmm = uimenu(hcm, 'Separator', 'on', 'Label', 'Add panel');
            uimenu(hmm, 'Label', 'Right...' , 'Callback', {@obj.cbAddPanel, 'right', ancestor(par, 'uiextras.HBoxFlex'), hp.Parent});
            uimenu(hmm, 'Label', 'Left...'  , 'Callback', {@obj.cbAddPanel, 'left' , ancestor(par, 'uiextras.HBoxFlex'), hp.Parent});
            uimenu(hmm, 'Label', 'Below...' , 'Callback', {@obj.cbAddPanel, 'below', ancestor(par, 'uiextras.VBoxFlex')});
            uimenu(hcm, 'Label', 'Remove panel...', 'Callback', {@obj.cbRemovePanel, hp});
            
            %Add / remove / apply / export / import arrangements of panels
            uimenu(hcm, 'Separator', 'on', 'Label', 'Arrangement');                        
            
            %Add to list of tab panels
            %obj.hTabPanel(end+1) = hp; EDIT. C.Szczyglowski 26/10/18
            obj.hTabPanel = [obj.hTabPanel, hp];
            
            %Ensure everything else in sync
            obj.update;
            
            function i_selectionChange(hc, ~)
            
                %Locate any Tree views
                idx = find(arrayfun(@(x)isa(x, 'mvc.view.Tree'), obj.Views));
                if isempty(idx)
                    
                    %Nothing to do
                    return;
                    
                end
                
                %Are any of these in the current tab panel ?
                [b, jdx] = ismember([obj.Views(idx).UIContainer], flip(hc.Children));
                if ~any(b)
                    
                    %Nothing to do
                    return;
                    
                end
                
                %If there is only one tab in the panel
                if numel(hc.Children) == 1
                    
                    %Short pause (not too sure why necessary - some sort of race condition ?)
                    pause(0.1);
                    
                end
                
                %Show/hide trees depending on whether they are im the selected tab
                [obj.Views(idx(b)).TreeVisible] = jdx(b) == hc.Selection;
                
            end
            
        end
        
        function removeAllPanels(obj)
            
            %Pass it on, working backwards from end
            arrayfun(@(x)removePanel(obj, x), numel(obj.hTabPanel):-1:1);
            
        end
        
        function removePanel(obj, hp)
            
            %Remove what ?
            if isnumeric(hp)
                
                %Treat as index directly
                idx = hp;
                
            elseif ishandle(hp)
                
                %Search in list
                idx = find(hp == obj.hTabPanel);
                assert(~isempty(idx), 'panel not found');
            
                %If the panel has no content
                if isempty(hp.Children)
                    
                    %Just crack on - deleting it wont hurt anyone
                    
                elseif ~obj.confirm(['Remove panel AND the ', num2str(numel(hp.Children)), ' view(s) it contains'])
                    
                    %Bail out
                    return;
                    
                end
                
            else
                error('bad input');
            end
            
            %If index specified by a two-element numeric
            if numel(idx) == 2
                
                %Treat it as directly [column, row] index into panels
                par = obj.hPanels.Children(idx(1)).Children(idx(2));
                
                %Delete it
                delete(par);
                
            else
                
                %Treat it as index into list of tab panels - get the corresponding parent
                par = get(obj.hTabPanel(idx), 'Parent');
                
                %Do it
                delete(obj.hTabPanel(idx));
                obj.hTabPanel(idx) = [];
                
                %If the parent now has no children
                if isempty(par.Children)
                    
                    %Then delete the parent too
                    delete(par);
                    
                end
                
            end
            
            %Ensure everything else in sync
            obj.update;
            
        end
        
        function vw = addView(obj, nam, hp, actionIfNotListed)
            
            %What to do if supplied view does not match anything we know about ?
            if nargin < 4
                actionIfNotListed = 'error';
            end
            
            %Add what ?
            if nargin < 2 || isempty(nam)
                
                %Look for match from this list
                sv = obj.SupportedViews(:,1);
                
                %Optionally allow user to browse for a new view type on the fly
                if obj.OfferBrowseForNewView
                    sv{end+1} = ' browse...';
                end
                
                %Ask the user
                [sel, bOK] = listdlg('ListString', sv, ...
                    'SelectionMode', 'single');
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    vw = [];
                    return;
                
                elseif sel == size(obj.SupportedViews,1) + 1
                    
                    %Browse for new view
                    sel = i_browse(obj);
                    if isempty(sel)
                        
                        %Bail out
                        vw = [];
                        return;
                        
                    end
                                       
                end
            
            elseif isnumeric(nam)
                
                %Treat as index directly
                sel = nam;
             
            elseif isa(nam, 'function_handle')
                
                %Caller is passing in view function directly
                sel = {nam};
            
            elseif strcmp(nam, 'browse')
                
                %Pass it on
                sel = i_browse(obj);
                if isempty(sel)
                    
                    %Bail out
                    vw = [];
                    return;
                    
                end

            else
                
                %Look up in list
                [b, sel] = ismember(nam, obj.SupportedViews(:,1));
                
                %Allow specific class to be specified by name too
                if ~any(b)
                    [b, sel] = ismember(nam, cellfun(@func2str, obj.SupportedViews(:,2), 'UniformOutput', false));
                end
                
                %Still nothing ?
                if ~all(b)
                    
                    %What to do ?
                    switch actionIfNotListed
                        
                        case 'allow'
                            
                            %Hope for the best
                            sel = {str2func(nam)};
                            
                        case {'skip', 'ignore'}
                            
                            %Do not add the view
                            vw = [];
                            return;
                            
                        otherwise
                            
                            %Treat as error
                            error(['no match with view ''', nam, '''']);
                            
                    end
                    
                end
                
            end
            
            %Add to which panel ?
            if nargin < 3 || isempty(hp)
                
                %Go with the last one in list
                hp = obj.hTabPanel(end);
                
            elseif ishandle(hp)
                
                %Just go with it
                
            elseif isnumeric(hp)
                
                %Treat as index into list
                hp = obj.hTabPanel(hp);
                
            else
                error('bad input');
            end
            
            %Careful
            err = [];
            try
                
                %Allow for multiple views to be created
                for i = 1:numel(sel)
                    
                    %Get view create function
                    if iscell(sel) && isa(sel{i}, 'function_handle')
                        
                        %Get the create function
                        fcn = sel{i};
                        
                        %Get title, if known
                        b = cellfun(@(x)isequal(fcn, x), obj.SupportedViews(:,2));
                        if any(b)
                            tit = obj.SupportedViews{b,1};
                        else
                            tit = [];
                        end
                        
                    elseif isnumeric(sel(i))
                        
                        %Get the create function and title
                        fcn = obj.SupportedViews{sel(i),2};
                        tit = obj.SupportedViews{sel(i),1};
                        
                    end
                    
                    %Create the view
                    vw(i) = fcn(obj.Model, 'Parent', hp);
                    
                    %Name it
                    vw(i).Title = tit;
                    
                    %Add to list
                    if isempty(obj.Views)
                        obj.Views = vw(i);
                    else
                        obj.Views(end+1) = vw(i);
                    end
                    
                    %Make it the selected view in this panel
                    set(hp, 'Selection', numel(get(hp, 'Children')));
                    
                    %What to do on selection change ?
                    vw(i).SelectionChangeFcn = @obj.onSelectionChange;
                    
                    %Set current section
                    vw(i).Selection = obj.Selection;
                    
                end
            
            catch err
            
                %For time being, move on
                
            end
            
            %Force an update
            obj.update;

            %If anything failed
            if ~isempty(err)
                
                %Rethrow
                rethrow(err);
                
            end
            
            %Helper function
            function fcn = i_browse(obj)
                
                %Persist in case we do this a number of times
                persistent FN;
                
                %User selected the 'browse...' option
                [f,p] = uigetfile({'*.m;*.p', 'MATLAB code files (*.m,*.p)'}, ...
                    'Select file...', char(FN));
                
                %Cancelled ?
                if isempty(f) || isequal(f, 0)
                    
                    %Bail out
                    fcn = [];
                    return;
                    
                else
                    
                    %Full filename
                    fn = fullfile(p,f);
                    
                    %Strip file extension
                    [~,f] = fileparts(f);
                    
                    %Convert fully-qualified filename to package name
                    tok = strsplit(fn, filesep);
                    pdx = find(cellfun(@(x)x(1) == '+', tok));
                    pdx(1:find(diff(pdx) > 1, 1, 'last')) = [];
                    tok = tok(pdx);
                    tok = cellfun(@(x)x(2:end), tok, 'UniformOutput', false);
                    tok{end+1} = f;
                    n = strjoin(tok, '.');
                    
                    %Ensure that the selected file is on the path
                    assert(strcmp(fn, which(n)), 'file must be on path');
                    
                    %And the corresponding function
                    fcn = {str2func(n)};
                    
                    %Make a note for next time
                    FN = fn;
                    
                    %Optionally add this view to the list ?
                    if confirm(obj, 'Remember this type of view ?')
                        
                        %Yes - add to list
                        obj.SupportedViews(end+1,:) = [{f}, fcn];
                        
                    end
                    
                end
                
            end
            
        end
        
        function vw = removeView(obj, nam, hp)
            
            %Current list of views ?
            str = {obj.Views.Title};
            
            %Remove what ?
            if nargin < 2 || isempty(nam) || isequal(nam, 'ask')
                
                %Infer from current view in current panel
                hc = flip(hp.Children); % because of the way hg orders children
                hc = hc(hp.Selection);
                sel = find([obj.Views.UIContainer] == hc);
                
                %If nothing
                if isempty(sel) || isequal(nam, 'ask')
                    
                    %Ask the user
                    [sel, bOK] = listdlg('ListString', str, ...
                        'SelectionMode', 'multiple', ...
                        'InitialValue', sel, ...
                        'PromptString', 'Remove the following view(s)...');
                    
                    %Cancelled ?
                    if isempty(sel) || ~bOK
                        
                        %Bail out
                        vw = [];
                        return;
                        
                    end
                
                end
                
                %Seek confirmation
                if ~obj.confirm(['Remove view(s) ', strjoin(str(sel), ',')])
                    vw = [];
                    return;
                end
                
            elseif isnumeric(nam)
                
                %Treat as index directly
                sel = nam;
                
            elseif isa(nam, 'mvc.view.Container')
                
                %Has the view been passed in?
                %   - Just use its title...
                nam = nam.Title;
                
                %Look for match
                [b, sel] = ismember(nam, str);
                assert(all(b), 'no match');
                
            else
                
                %Look for match
                [b, sel] = ismember(nam, str);
                assert(all(b), 'no match');
                
            end
            
            %Do it
            vw = obj.Views(sel);
            obj.Views(sel) = [];
            
            %Maintain views in sync (this shouldn't be necessary ?)
            par = unique([vw.Parent]);
            delete(setdiff([par.Contents], [obj.Views.UIContainer]));

            %Ensure all up to date
            obj.update;
            
        end
        
        function varargout = copyView(obj, hv, dest, varargin)
            
            %Caller may pass in a handle to the view to be copied
            if ishghandle(hv)
                
                %OK
            
            elseif isempty(hv)
                
                %No view passed in is also OK - used to trigger clean-up,
                % for example on completion of a sequence of copies to Word doc
                
            elseif isa(hv, 'mvc.view.Container')
                
                %Caller has passed in the actual view to be copied
                hv = hv.UIContainer;
                
            elseif isnumeric(hv)
                
                %Caller passing in a view by index into list
                hv = obj.Views(hv).UIContainer;
                
            else
                
                %Nothing else will do
                error('bad input');
                
            end
            
            %If we get here but have no view
            if isempty(hv)
                
                %Nothing to do at this point
                hf = [];
                
            else
                
                %Best way to capture the view ? Get the position of the view
                rect = hv.Position;
                
                %Create a new figure, with same position (so same aspect ratio)
                hf = figure('Position', rect, ...
                    'Visible', 'off');
                
                %Copy the content
                copyobj(hv, hf);
            
            end
            
            %Copy to where ?
            switch lower(dest)
                
                case 'figure'
                    
                    %Easy, just show us
                    figure(hf);
                    
                    %Caller may want it back
                    if nargout > 0
                        varargout{1} = hf;
                    end
                    
                    %Don't delete it before exit
                    hf = [];
                    
                case 'file'
                    
                    %Save to where ?
                    persistent FN;
                    [f,p] = uiputfile({'*.fig', 'MATLAB figure file (*.fig)'}, 'Save to...', char(FN));
                    if isempty(f) || isequal(f,0)
                        
                        %Cancelled
                        
                    else
                        
                        %Make a note for next time
                        FN = fullfile(p,f);
                        
                        %Save to fig file
                        hgsave(hf, FN);
                        
                    end
                    
                    %Caller may want it back
                    if nargout > 0
                        varargout{1} = FN;
                    end
                    
                case 'clipboard'
                    
                    %Easy
                    print(hf, '-dmeta');
                    
                case 'word'
                    
                    %Anything passed in ?
                    if isempty(varargin)
                        
                        %No - get help from Report Generator to create a new file,
                        % allowing for possibility of multiple docs already open
                        sfx = 1;
                        while true
                            try
                                doc = mlreportgen.dom.Document(['untitled-', num2str(sfx)], 'docx');
                                doc.open;
                                break
                            catch
                                sfx = sfx + 1;
                            end
                        end
                    
                    elseif ischar(varargin{1})
                        
                        %Name of document explicitly passed in
                        doc = mlreportgen.dom.Document(varargin{1}, 'docx');
                        doc.open;
                        
                    elseif isa(varargin{1}, 'mlreportgen.dom.Document')
                        
                        %Handle to document explicitly passed in
                        doc = varargin{1};
                        
                    else
                        error('bad input');
                    end
                    
                    %Any other inputs ?
                    for i = 2:numel(varargin)
                        
                        %Treat as text to be written to document before the view
                        doc.append(mlreportgen.dom.Paragraph(varargin{i}));
                        
                    end
                     
                    %If we have no figure at this point
                    if isempty(hf)
                        
                        %Caller is asking us to close up handle to Word
                        doc.close;
                        
                        %Delete any pending temp files
                        cellfun(@delete, obj.TempFiles);
                        obj.TempFiles = {};
                        
                        %Caller may want doc name back
                        if nargout > 0
                            
                            %Send it back
                            varargout{1} = doc.OutputPath;
                            
                        else
                            
                            %Show us
                            rptview(doc.OutputPath);
                            
                        end
                        
                    else
                        
                        %Is this really necessary ?
                        tn = tempname;
                        print(hf, [tn, '.png'], '-dpng');
                        
                        %Add the image
                        im = mlreportgen.dom.Image([tn, '.png']);
                        im.Style = {mlreportgen.dom.ScaleToFit};
                        doc.append(im);
                        
                        %Caller may want doc handle back
                        if nargout > 0
                            
                            %Can't delete the temp file yet, so add to pending list
                            obj.TempFiles{end+1} = [tn, '.png'];
                            
                            %Send it back
                            varargout{1} = doc;
                            
                        else
                            
                            %Close the document
                            doc.close;
                            
                            %Delete the temp file
                            delete([tn, '.png']);
                            
                            %Show us
                            rptview(doc.OutputPath);
                            
                        end
                    
                    end
                    
                case 'excel'
                    
                    %Copy to Excel handled by each type of view individually
                    
                otherwise
                    error('bad destination');
            end
            
            %We done with the figure
            delete(hf);
            
        end
        
        function vw = moveView(obj, hv, pdx)
            
            %Caller may pass in a handle to the view to be moved
            if ishghandle(hv)
                
                %OK
                
            elseif isa(hv, 'mvc.view.Container')
                
                %Caller has passed in the actual view to be moved
                hv = hv.UIContainer;
                
            elseif isnumeric(hv)
                
                %Caller passing in a view by index into list
                hv = obj.Views(hv).UIContainer;
                
            else
                
                %Nothing else will do
                error('bad input');
                
            end
            
            %Moving to a specific panel ?
            if nargin < 3 || isempty(pdx)
                
                %No panel specified, so create a new one
                ht = obj.addPanel;
            
            elseif ischar(pdx) && strcmpi(pdx, 'figure')
                
                %Reparent the view, adjust units and position while we're at it
                set(hv, 'Parent', figure, 'Units', 'normalized', 'Position', [5, 5, 90, 90] ./ 100);
                
                %Nothing else required
                return;
            
            elseif ischar(pdx) && strcmpi(pdx, 'right')
                
                %Pass it on
                ht = obj.addPanel(pdx, hv.Parent.Parent.Parent);
            
            elseif ischar(pdx) && strcmpi(pdx, 'below')
                
                %Pass it on
                ht = obj.addPanel(pdx, hv.Parent.Parent);
                
            else
                
                %Get handle to destination tab panel - pdx(1) identifies which panel horizontally
                hb = ancestor(hv, 'uiextras.HBoxFlex');
                hc = flip(hb.Children); % because of the order in which children are stored
                ht = hc(pdx(1));
                
                %If supplied, pdx(2) identifies which panel vertically
                if numel(pdx) == 1
                    pdx(2) = 1;
                end
                hc = flip(ht.Children); % because of the order in which children are stored
                ht = hc(pdx(2));
                
            end
            
            %And move the view from source to destination panel
            hv.Parent = ht;
            
            %Flush pending graphics before forcing selection
            drawnow;
            
            %Select it
            ht.Selection = numel(ht.Children);
            
            %Force an update
            obj.update;            
            
        end        
        
        function str = dlgtitle(obj, varargin)
   
            %Simple
            str = obj.GuiName;
            
        end
        
        function varargout = arrangeViews(obj, action, name)
        
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
                    A = obj.ViewArrangement;
                    
                    %Name provided ?
                    if isempty(name)
                        
                        %Ask the user
                        val = inputdlg({'Name of this arrangement:'}, obj.GuiName, 1, {char(A.Name)});
                        
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
                                if ~isequal('Yes', questdlg('Overwrite existing arrangement ?', obj.GuiName, 'Yes', 'No', 'Yes'))
                                    
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
                    obj.ViewArrangement = obj.SupportedArrangements(idx);

                case {'import', 'default'}
                    
                    %Import from where ?
                    if strcmp(action, 'default')
                        
                        %Go with this file
                        name = fullfile(fileparts(which(class(obj))), 'defaults.vaf');
                        
                        %If not found
                        if exist(name, 'file') ~= 2
                            
                            %Return silently
                            return;
                            
                        end
                    
                    end
                    
                    %What name has been passed in ?
                    if isempty(name)
                        
                        %Not specified - ask the user
                        [f,p] = uigetfile(obj.ViewArrangementFileMask, 'Import arrangements from...', ...
                            char(obj.ViewArrangementFile));
                        
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
                        obj.ViewArrangementFile = fn;
                        
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
                        obj.ViewArrangementFile = name;
                        
                        %And store the detail
                        obj.SupportedArrangements = TEMP.SupportedArrangements;
                        
                    elseif isstruct(name)
                        
                        %Assume it defines the arrangements directly
                        obj.ViewArrangementFile = [];
                        
                        %And store the detail
                        obj.SupportedArrangements = name;
                        
                    else
                        error('bad input');
                    end
                    
                    %And apply the first in the list
                    obj.ViewArrangement = obj.SupportedArrangements(1);
                    
                case 'export'
                    
                    %Save to where ?
                    if isempty(name)
                        
                        %Ask the user
                        [f,p] = uiputfile(obj.ViewArrangementFileMask, 'Export arrangement to file...', char(obj.ViewArrangementFile));
                        
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
                        obj.ViewArrangementFile = fn;
                        
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
                        'Name', obj.GuiName, ...
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
        
        function updateArrangementMenu(obj, hm, varargin)
            
            %Anything to do ?
            if isempty(hm) || ~isvalid(hm)
                
                %No
                return;
                
            end
                        
            %Clean sheet
            delete(hm.Children)
            
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
                    'Callback', {@obj.cbArrangeViews, 'apply', x}), ...
                    nam, sep);
                
            end
            
            %Additional options to manage arrangements
            hm = uimenu(hm, 'Label', 'Manage', 'Separator', 'on');
            if ~isempty(obj.CurrentArrangement)
                uimenu(hm, 'Label', ['Save arrangement ''', obj.CurrentArrangement, ''''], 'Callback', {@obj.cbArrangeViews, 'save'});
            end
            uimenu(hm, 'Label', 'New arrangement...', 'Callback', {@obj.cbArrangeViews, 'new'});
            uimenu(hm, 'Label', 'Remove arrangement(s)...', 'Callback', {@obj.cbArrangeViews, 'remove'});
            uimenu(hm, 'Separator', 'on', 'Label', 'Export arrangements...', 'Callback', {@obj.cbArrangeViews, 'export'});
            uimenu(hm, 'Label', 'Import arrangements...', 'Callback', {@obj.cbArrangeViews, 'import'});
            
        end
        
        function updateViewMenu(obj, hm, varargin)
            
            %Anything to do ?
            if isempty(hm) || ~isvalid(hm)
                
                %No
                return;
                
            end
            
            %Clean sheet
            delete(hm.Children)
            
            %Any views available ?
            if isempty(obj.SupportedViews)
                
                %No
                uimenu(hm, 'Label', 'none available', 'enable', 'off');
                
            else
                
                %Get name of views
                nam = obj.SupportedViews(:,1);
                
                %With separator ?
                sep = cellfun(@(x)x(1) == '|', nam);
                nam(sep) = cellfun(@(x)x(2:end), nam(sep), 'UniformOutput', false);
                sep = arrayfun(@mvc.mixin.UiTools.bool2offon, sep, 'UniformOutput', false);
                
                %Add quick access to additional views
                cellfun(@(x,y)uimenu(hm, 'Label', x, 'Separator', y, 'Callback', {@obj.cbAddView, x, varargin{:}}), nam, sep);
            
            end
            
            %Offer option to "browse" for a new view ?
            if obj.OfferBrowseForNewView
                
                %Add to menu
                uimenu(hm, 'Separator', 'on', 'Label', ' browse...', 'Callback', {@obj.cbAddView, 'browse', varargin{:}});
                
            end                

        end
                        
    end
    
    methods (Access = protected) % Callbacks
        
        function bArrange = initialiseFigure(obj, hf)
            
            %Initialise return arg
            bArrange = true;
                    
            %Make a note of figure
            obj.hFigure = hf;
            
            %There might be nothing to do
            if isempty(hf)
                return;
            end
            
            %Action on close ?
            set(hf, 'CloseRequestFcn', @obj.close);
            
            %Set name accordingly
            set(hf, 'Name', obj.GuiName);
            
            %Store the view manager itself in the figure appdata
            setappdata(hf, matlab.lang.makeValidName(class(obj)), obj);
            
            %Special handling if caller has created this view using the "new" method
            hf_caller = obj.currentFigure(gcbo);
            if ~isempty(hf_caller)
                
                %IF the current callback figure is ALSO associated with same view
                obj_caller = getappdata(hf_caller, matlab.lang.makeValidName(class(obj)));
                if isa(obj_caller, class(obj))
                    
                    %Initialise position accordingly
                    pos = get(hf_caller, 'Position');
                    pos(1:2) = pos(1:2) + [25, -25];
                    set(hf, 'Position', pos);
                    
                    %Copy the arrangement of views from the calling object
                    obj.ViewArrangement = obj_caller.ViewArrangement;
                    
                    %And no need to make up an arrangement for ourselves now
                    bArrange = false;
                    
                end
                
            end
            
            %If we get to here and the arrange flag remains true
            if bArrange && isa(obj.Model, 'mvc.mixin.Serializable')
                
                %Then see if we can find a last position
                % (TODO: is the Model the right place to store this ?)
                pos = obj.Model.LastPosition;
                if ~isempty(pos)
                    
                    %To allow for laptop / docking station combos, need to check
                    % that whatever lastpos was it is still visible on any available screen.
                    %
                    %Start with polygon for position
                    [xq, yq] = i_pos2poly(pos);
                    
                    %Then for each available monitor position
                    mpos = get(0, 'MonitorPositions');
                    for i = 1:size(mpos,1)
                        
                        %Create polygon for monitor
                        [xm, ym] = i_pos2poly(mpos(i,:));
                        
                        %Are we in it ?
                        if any(inpolygon(xq, yq, xm, ym))
                            
                            %Yes - apply position
                            set(hf, 'Position', pos);
                            break;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            %Where the helper function is
            function [x, y] = i_pos2poly(pos)
                x = [pos(1), sum(pos([1 3])) .* [1 1], pos(1) .* [1 1]];
                y = [pos(2) .* [1 1], sum(pos([2 4])).* [1 1], pos(2)];
            end
            
        end
        
        function update(obj, varargin)
             
            %Refresh figure name
            hf = ancestor(obj.UIContainer, 'Figure');
            nam = obj.GuiName;
            if isa(obj.Model, 'mvc.mixin.Serializable')
                
                %Append relevant details
                if obj.Model.ChangedSinceLastSave
                    nam = [nam, '*'];
                end
                if ~isempty(obj.Model.SerializeFile)
                    nam = [nam, ' (', obj.Model.SerializeFile, ' in ', obj.Model.SerializePath, ')'];
                end
                
            end
            set(hf, 'Name', nam);
           
            %Refresh list of tab-titles
            if isempty(obj.Views)
                
                %Nothing to do
                return;
                
            end
            
            %Safety net - should not be necessary
            obj.Views(~isvalid(obj.Views)) = [];
            
            %Who belongs to whom ?
            for i = 1:numel(obj.hTabPanel)
                
                %Look for an annotation associated with this panel
                ha = findobj(allchild(get(obj.hTabPanel(i), 'Parent')), 'Type', 'annotationpane');
                
                %Look for match
                [b, idx] = ismember(flip(get(obj.hTabPanel(i), 'Children')), [obj.Views.UIContainer]);
            
                %Anything ?
                if ~any(b)
                    
                    %Need annotation
                    if ~isempty(ha)
                        ha.Parent.Heights(1) = -1;
                        ha.Parent.Heights(2:end) = 0;
                        ha.Visible = 'on';
                    end
                    
                    %No
                    continue;
                    
                end
                    
                %Don't need annotation
                if ~isempty(ha)
                    ha.Parent.Heights(1) = 0;
                    ha.Parent.Heights(2:end) = -1;
                    ha.Visible = 'off';
                end
                
                %Update tab titles accordingly (carefully, so as to avoid confusing error msg
                % if view creation has already failed for other reasons)
                ttCurrent = get(obj.hTabPanel(i), 'TabTitles');
                ttExpected = {obj.Views(idx(b)).Title};
                n = min([numel(ttCurrent), numel(ttExpected)]);
                ttCurrent(1:n) = ttExpected(1:n);
                set(obj.hTabPanel(i), 'TabTitles', ttCurrent);

                %Ensure consistency across context menus
                hcm = get(obj.hTabPanel(i), 'UIContextMenu');
                set(obj.hTabPanel(i), 'TabContextMenus', repmat({hcm}, 1, numel(ttCurrent)));                

            end
            
            %Stash the pixel positions of the views within the figure
            obj.ViewPixelPosition = getViewPixelPosition(obj);
            
        end
        
        function cbContext(obj, hc, ~, hp)
            %
            % Manage enable status of context menus
             
            %Careful
            try
                
                %Refresh the list of available views
                obj.updateViewMenu(findobj(hc, 'Label', 'Add view'), hp);
                
                %Refresh the list of available arrangements (and other related options)
                obj.updateArrangementMenu(findobj(hc, 'Label', 'Arrangement'));

                %Each panel sits in a vbox
                hv = hp.Parent;
                
                %Where are we in child list ?
                %pdy = find(hp == flip(hv.Children));
                
                %Each vbox sits in an hbox
                hh = hv.Parent;
                
                %Where are we in child list ?
                %pdx = find(hv == flip(hh.Children));
                
                %Find the 'Remove panel' submenu
                hv = findobj(hc, 'Label', 'Remove panel...');
                
                %Enable option to remove panel ONLY if this is not the last one
                hv.Enable = mvc.mixin.UiTools.bool2offon(sum(obj.PanelRowsPerColumn) > 1);
                
                %The currently selected tab of the current panel
                cdx = hp.Selection;
                
                %Find the 'Remove view' and 'Move to' submenu
                hr = findobj(hc, 'Label', 'Remove view...');
                hm = findobj(hc, 'Label', 'Move view to');
                delete(hm.Children);
                
                %Anything ?
                if isempty(cdx) || cdx == 0
                    
                    %No - disable options to remove and move view
                    hr.Enable = 'off';
                    hm.Enable = 'off';
                    
                else
                    
                    %Enable options to remove and move view
                    hr.Enable = 'on';
                    hm.Enable = 'on';
                    
                    %So the handle of the view to be moved
                    hc = flip(hp.Children); % because of the way children are ordered in reverse
                    hv = hc(cdx);
                    
                    %Placeholders
                    lab = {};
                    arg = {};
                    
                    %What panels do we have ?
                    for i = 1:numel(hh.Children)
                        for j = 1:numel(hh.Children(i).Children)
                            
                            %No point sending a view to the panel that alredy owns it
                            if hp == hh.Children(i).Children(j)
                                continue;
                            end
                            
                            %Work out position of destination panel
                            colrow = [numel(hh.Children) - i, numel(hh.Children(i).Children) - j] + 1;
                            
                            %Suitable label
                            if numel(hh.Children(i).Children) == 1
                                lab = ['Panel ', num2str(colrow(1))];
                            else
                                lab = ['Panel [', num2str(colrow(1)), ', ', num2str(colrow(2)), ']'];
                            end
                            
                            %Add to menu
                            uimenu(hm, 'Label', lab, 'Callback', {@obj.cbMoveView, hv, colrow});
                            
                        end
                    end
                    
                    %OR we could move the view to a new panel
                    hmm = uimenu(hm, 'Label', 'New panel');
                    uimenu(hmm, 'Label', 'Right', 'Callback', {@obj.cbMoveView, hv, 'right'});
                    uimenu(hmm, 'Label', 'Below', 'Callback', {@obj.cbMoveView, hv, 'below'});
                
                    %OR we could move the view to a separate figure window
                    uimenu(hm, 'Label', 'New figure', 'Callback', {@obj.cbMoveView, hv, 'figure'});
                    
                end
                               
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
                       
        end
        
        function cbAddPanel(obj, hc, ~, pos, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.addPanel(pos, varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end             
        
        function cbRemovePanel(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.removePanel(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbAddView(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.addView(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbRemoveView(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.removeView(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbCopyView(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.copyView(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbMoveView(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.moveView(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbArrangeViews(obj, hc, ~, varargin)
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pass it in
                obj.arrangeViews(varargin{:});
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.GuiName, 'modal'));
                
            end
            
        end
        
        function cbOnFigResize(obj, ~, ~)
            %cbOnFigResize Callback for whenever the figure window is
            %resized.
            %
            % Need to update the stored value for 'ViewPixelPosition' as we
            % have changed the size of the layout!
            
            %Stash the pixel positions of the views within the figure
            obj.ViewPixelPosition = getViewPixelPosition(obj);
            
        end
        
        function cbImport(obj, ~, ~)
            
            %Pass it on to the model
            if ~isempty(obj.Model)
                import(obj.Model)
            end
            
        end
        
    end
        
    methods ( Access = private ) % Callbacks for model and selection changed
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onSelectionChange(obj, sel)
            
            %Make a note of new selection
            obj.Selection = sel;
                
        end        
       
    end
    
    methods ( Access = protected ) % view specific callbacks for mouse interaction
        
        function cbViewManagerMouseAction(obj, hFig, evt, fn, varargin)
            %cbViewManagerMouseAction Handles all callbacks for the view
            %manager mouse actions (MouseMotion, ButtonDown, ButtonUp, etc)
            %
            % Actions:
            %   - Get current mouse position
            %   - Get current view that the mouse is hovering over
            %   - Pass the method onto the view-specific callback

             %Get the current location of the mouse
            %   - Use 'evt.Point' as it is with respect to the parent
            %     figure whereas 'get(0, 'PointerLocation')' is with 
            %     respect to the the computer monitor(s).
            mousePos = evt.Point;  
            
            %What view have we clicked on?
            cvw = obj.getCurrentView(mousePos);
            
            if isempty(cvw) %Escape route
                return
            end
                        
            %Pass it on
            try
                fn(cvw, hFig, evt, varargin{:});
            catch ME
                rethrow(ME);
            end
            
        end
        
        function cvw = getCurrentView(obj, mousePos)
            %getCurrentView Returns the handle to the current view that the
            %mouse is positioned over.
            %
            % TODO - Account for views that are hidden in tab panels.
            
            %Sensible default
            cvw = [];
            
            %Gather the positions of each view
            viewPos = obj.ViewPixelPosition;
            
            if isempty(viewPos) %Escape route
                return
            end
            
            %Check to see if we are inside one of the views
            onViewX = and(mousePos(1) > viewPos(:, 1), mousePos(1) < viewPos(:, 1) + viewPos(:, 3));
            onViewY = and(mousePos(2) > viewPos(:, 2), mousePos(2) < viewPos(:, 2) + viewPos(:, 4));
            onView  = and(onViewX, onViewY);
            
            if ~any(onView) %Escape route
                return
            end
            
            %What view are we on?
            cvw = obj.Views(onView);
            
            if numel(cvw) > 1
                warning('Cannot have two active views at once!');
                cvw = [];
            end
            
        end        
        
        function viewPixelPos = getViewPixelPosition(obj)
            %getViewPixelPosition Retrieves the relative pixel position of
            %each view in the 'ViewManager'.
            %
            % For each view, we traverse up the object hierachy and grab
            % the pixel position until we reach the parent figure.
            %
            % TODO - Account for the case where views are in a tab panel
            % and can be overlaid on top of one another.
            
            %Sensible default
            viewPixelPos = [];
            
            if isempty(obj.Views) %Escape route
                return
            end
            
            %Preallocate
            viewPixelPos = zeros(numel(obj.Views), 4);
            
            %For each view in turn ...
            for ivw = 1 : numel(obj.Views)
                
                %Start with the current object
                gObj = obj.Views(ivw);
                
                %Get its position
                pos = i_getPixelPos(gObj);
                
                %Traverse up the tree
                while ~isempty(gObj.Parent)                    
                    %New object
                    gObj = gObj.Parent;
                    %Stop once we get to the parent of the ViewManager
                    if gObj == obj.Parent
                        break
                    end
                    %Get its position
                    pos(end + 1, :) = i_getPixelPos(gObj);                    
                end
                
                %Calculate view position [x, y, w, h]
                x = sum(pos(:, 1), 1);
                y = sum(pos(:, 2), 1);
                viewPixelPos(ivw, :) = [x, y, pos(1, 3), pos(1, 4)];
                %viewPixelPos(ivw, :) = [x, y, x + pos(1, 3), y + pos(1,
                %4)];%--> Hot-spot
                
            end
            
            function pos = i_getPixelPos(gObj)
                %i_getPixelPos Retrieves the pixel position of the view.
                
                %Get current units
                cu = gObj.Units;
                
                %Enforce pixels
                gObj.Units = 'pixels';
                
                %Grab position
                pos = gObj.Position;
                
                %Return to original state
                gObj.Units = cu;
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function obj = defaultView(varargin)
            
            %Create view manager
            obj = mvc.view.ViewManager(varargin{:});
            
            %Add a Tree view on the left
            obj.addView('Tree');
            
            %Add Properties and Drawing views on the right
            obj.addPanel;
            obj.addView('Properties');
            obj.addView('Drawing');
            
        end
        
    end
    
end