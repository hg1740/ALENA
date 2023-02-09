classdef (ConstructOnLoad) Nameable < matlab.mixin.SetGet ...
        & matlab.mixin.Heterogeneous ...
        & matlab.mixin.CustomDisplay ...
        & matlab.mixin.Copyable
    %
    % The Nameable class provides basic functionality required of any object that has metadata
    %  such as Name, Type, Description, along with the "proper" data associated with that object.
    %
    % On the assumption that this class will be subclassed by application-specific objects
    %  that manage potentially large numbers of properties, Nameable implements a scalable
    %  property viewer / editor, in which properties are presented to the user in groups,
    %  each group presented in a separate tab-page of a dialog.
    %
    % Examples: where fcn constructs a class that implements Nameable
    %
    %  obj = fcn('Name', 'Peter', 'Type', 'Employee', 'Description', {'Peter is an employee'});
    %
    % Basic get/set functionality:
    %
    %   nam = get(obj, 'Name');
    %   nam = obj.Name;
    %   prp = get(obj, {'Name', 'Type'});
    %
    % Property view/edit:
    %
    %  propedit(obj)
    
    properties (AbortSet, SetObservable)
        
        %Everything has a name, type and description
        Name % char; %% property data typing not supported in R2012b
        Type % char; %% property data typing not supported in R2012b
        Description = {''};
        
        %Are properties locked (i.e. do not permit editing of values)
        PropertiesLocked = false;
        
    end
    
    properties (Dependent)
        
        %Saves concatenation in various places
        NameAndType;
        
    end
    
    properties (SetAccess = protected, Transient)
        
        %Does user have direct access to the state of the "...Locked" flags ?
        PropertiesLockable = true;
        
        %Some properties should be always/never locked, regardless of value of the PropertiesLocked flag
        PropertiesNeverLocked = {'Name', 'Description', 'PropertiesLocked'};
        PropertiesAlwaysLocked = {'Type'};
        
        %How should properties be grouped (e.g. when presented by "propedit") ?
        PropertyGroups = struct('Title', 'General', ...
            'Visible', true, ...
            'Enable', true, ...
            ... % 'Export', true, ...
            'Properties', {{'Name', 'Type', 'Description'}}, ...
            'Description', {{'',    '',     ''}}, ...
            'TooltipString', {{'',  '',     ''}}, ...
            'PopupOptions',  {{[],  [],     []}}, ...
            'CreateFcn', []);
        
    end
    
    properties (Hidden, Transient)
        
        %If Name not (yet) explicitly set, how should it be derived ?
        DefaultName;
        DefaultField; % optionally use another field to derive name
        
    end
    
    methods % construction / destruction
        
        function obj = Nameable(varargin)
            
            %Anything to pass on ?
            if ~isempty(varargin)
                
                %Pull out any special tokens
                %   - TODO : Switch to 'contains' once we are sure of ver.
                ind = find(cellfun(@(x) ~isempty(strfind(x, '-')), ...
                    varargin(1 : 2 : end - 1)));
                varargin([ind + 1, 2 * ind]) = [];
                
                %Yes
                if numel(varargin) > 1
                    set(obj, varargin{:});
                end
                
            end
            
            %Grant user access to the "...Locked" flags ?
            if obj.PropertiesLockable
                
                %Extend property groups
                obj.addPropertyGroup('General', ...
                    'PropertiesLocked', 'PropertiesLocked');

            end
            
            %Extend context, if applicable
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, 'Edit>Properties...', 'propedit');
            end
            
        end
        
    end
    
    methods % setters / getters
        
        function str = get.NameAndType(obj)
            
            %Allow for Type to be unset, or (commonly the case) to be same as Name
            if isempty(obj.Type) || isequal(obj.Name, obj.Type)
                
                %Just go with name
                str = obj.Name;
                
            else
                
                %Simple concatenation
                str = [obj.Name, ' (', obj.Type, ')'];
        
            end
            
        end
        
        function set.Name(obj, val)
            
            %Validate
            validateattributes(val, {'char'}, {});
            
            %Make a note
            obj.Name = val;
            
        end
        
        function str = get.DefaultName(obj)
            
            %Start from whatever we've got
            str = obj.DefaultName;
            
            %If nothing yet
            if isempty(str)
                
                %Derived from classname
                tok = strsplit(class(obj), '.');
                str = tok{end};
                
            end
            
        end
            
        function str = get.Name(obj)
            
            %Start from whatever we've got
            str = obj.Name;
            
            %If set
            if ~isempty(str)
                
                %Nothing else to do
                
            elseif ~isempty(obj.DefaultField) && ~isempty(obj.(obj.DefaultField))
                
                %Go with it
                str = obj.(obj.DefaultField);
                
            elseif ~isempty(obj.DefaultName)
                
                %Go with it
                str = obj.DefaultName;
                
            else
                
                %Something is better than nothing
                str = 'unnamed';
                
            end
            
        end
        
        function set.Type(obj, val)
            
            %Validate
            validateattributes(val, {'char'}, {});
            
            %Make a note
            obj.Type = val;
            
        end
        
        function str = get.Type(obj)
            
            %Start from whatever we've got
            str = obj.Type;
            
            %If not set
            if isempty(str)
                
                %Go with classname instead (exluding package detail ?)
                tok = strsplit(class(obj), '.');
                str = tok{end};
                
            end
            
        end
        
        function set.Description(obj, val)
            
            %Caller can provide a string
            if ischar(val)
                
                %But it is always stored as a cellstr
                val = {val};
                
            elseif isempty(val)
                
                %Nothing is ok too, but ensure it is a cell array of nothing
                val = {};
                
            end
            
            %Validate
            validateattributes(val, {'cell'}, {});
            
            %Make a note
            obj.Description = val;
            
        end
        
    end
    
    methods (Access = protected) % to make collection copyable
        
        function cpy = copyElement(obj)
            
            %If object implements dynamicable
            if isa(obj, 'mvc.mixin.Dynamicable')
                
                %Then do NOT call base-class copyElement - it has already been called
                cpy = obj;
                
            else
                
                %Start with base-class
                cpy = copyElement@matlab.mixin.Copyable(obj);
                
            end
            
            %Check property groups
            for i = 1:numel(cpy.PropertyGroups)

                %Pass it on
                cpy.PropertyGroups(i) = i_updateFunctionHandles(cpy, cpy.PropertyGroups(i));
                
            end
            
        end
        
        function S = i_updateFunctionHandles(obj, S)
            
            %Check each field in turn
            fld = fieldnames(S);
            for i = 1:numel(fld)
                
                %What have we got ?
                if isstruct(S.(fld{i}))
                   
                    %Recurse
                    S.(fld{i}) = i_updateFunctionHandles(obj, S.(fld{i}));
                
                elseif iscell(S.(fld{i}))
                    
                    %Recurse
                    S.(fld{i}) = cellfun(@(x)i_updateFunctionHandle(obj, x), S.(fld{i}), 'UniformOutput', false);
                    
                elseif isa(S.(fld{i}), 'function_handle')
                    
                    %Update it
                    S.(fld{i}) = i_updateFunctionHandle(obj, S.(fld{i}));
                    
                    %                 elseif iscell(S.(fld{i})) && numel(S.(fld{i})) > 0 && isa(S.(fld{i}){1}, 'function_handle')
                    %
                    %                     %Update it
                    %                     S.(fld{i}){1} = i_updateFunctionHandle(obj, S.(fld{i}){1});
                    
                end
                
            end
            
        end
        
        function fcn = i_updateFunctionHandle(obj, fcn) %#ok<INUSL>
        
            %Anything to do ?
            if ~isa(fcn, 'function_handle')
                
                %No
                return;
                
            end
            
            %Update the function call by casting the handle to a string
            str = func2str(fcn);
            
            %Then back to a function handle
            % (use eval, rather than str2func, because the latter has problems if str contains patterns like .FieldName{i}
            fcn = eval(str);
            
        end
        
    end
        
    methods (Access = protected) % custom display
        
        function str = displayEmptyObject(obj)
            
            %What have we got ?
            str = ['empty ', class(obj)];
            
            %If caller did not want the string back
            if nargout == 0
                
                %Display it
                disp(str);
                
            end
            
        end
        
        function str = displayScalarObject(obj, bShort)
            
            %Go with basic properties
            str = obj.NameAndType;
            
            %Caller can ask for a short string only
            if nargin > 1 && bShort
                
                %Nothing to add
                
            else
                
                %Anything to add from description ?
                b = ~cellfun(@isempty, obj.Description);
                if any(b)
                    
                    %Yes - append
                    str = [str, ': ', strjoin(obj.Description(b), ', ')];
    
                end
                
                %Maybe needed later
                bShort = false;
                
            end
            
            %And more (TODO: is this the right place for this ?)
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %What have we got ?
                str = [str, ': ', displayNonScalarObject(obj.Children, bShort)];
                
            end
            
            %If caller did not want the string back
            if nargout == 0
                
                %Display it
                disp(str);
                
            end
            
        end
        
    end
    
    methods (Access = protected, Sealed) % custom display
        
        function str = displayNonScalarObject(obj, bShort)
            
            %How many items to list, explicitly ?
            nList = min([2, numel(obj)]);
            
            %Caller can ask for a short string only
            if nargin > 1 && bShort
                
                %Just go with
                str = strjoin({obj(1:nList).Name}, ', ');

                
            else
                
                %What have we got ?
                str = [obj.convertDimensionsToString(obj), ' ', class(obj)];
                
                %Add a little more
                str = [str, ': ', strjoin({obj(1:nList).NameAndType}, ', ')];
                
            end
            
            %And the extras
            if numel(obj) > nList
                str = [str, '...and ', num2str(numel(obj) - nList), ' more'];
            end
            
            %If caller did not want the string back
            if nargout == 0
                
                %Display it
                disp(str);
                
            end
            
        end
        
    end
    
    methods (Sealed) % custom display functions needs to be accessible to peers of this class
        
        function str = dlgtitle(obj, varargin)
            %
            % Returns a simple, one-line, string that can be assigned to a figure's 'Name' property.
            
            %What have we got ?
            if numel(obj) == 0
                
                %Nothing
                str = displayEmptyObject(obj, varargin{:});
                
            elseif numel(obj) > 1
                
                %More than one
                str = displayNonScalarObject(obj, varargin{:});
                
            else
                
                %Just the one
                str = displayScalarObject(obj, varargin{:});
                
            end
            
        end
        
    end
    
    methods (Sealed) % Property view / edit
        
        function sortPropertyGroup(obj, action, varargin)
            
            %Do what ?
            switch action
                
                case 'move to start'
                    
                    %Locate the specified property groups
                    b = ismember({obj.PropertyGroups.Title}, varargin);
                    
                    %So the new sort order is
                    sdx = [find(b), find(~b)];
                    
                case 'move to end'
                    
                    %Locate the specified property groups
                    b = ismember({obj.PropertyGroups.Title}, varargin);
                    
                    %So the new sort order is
                    sdx = [find(~b), find(b)];
                    
                otherwise
                    error('bad action');
            end
            
            %Apply change
            obj.PropertyGroups = obj.PropertyGroups(sdx);
            
        end
        
        function setPropertyGroup(obj, tit, varargin)
            
            %Does this property group already exist ?
            tdx = find(strcmp(tit, {obj.PropertyGroups.Title}));
            
            %To set anything, it MUST
            assert(~isempty(tdx), ['property group ''', tit, ''' not found']);
            
            %Set whatever
            for i = 1:2:numel(varargin)
                obj.PropertyGroups(tdx).(varargin{i}) = varargin{i+1};
            end
            
        end
        
        function b = getPropertyAttributes(obj, tit, nam, def)
            
            %Caller may pass in an index directly
            if isnumeric(tit)
                
                %Just go with it
                tdx = tit;
                
            else
                
                %Does this property group exist ?
                tdx = find(strcmp(tit, {obj.PropertyGroups.Title}));
                
                %To get anything, it MUST
                assert(~isempty(tdx), ['property group ''', tit, ''' not found']);
            
            end
            
            %             %Does this property name exist within this group ?
            %             ndx = find(strcmp(nam, obj.PropertyGroups(tdx).Properties));
            %
            %             %To get anything, it MUST
            %             assert(~isempty(ndx), ['property ''', nam, ''' in group ''', tit, ''' not found']);
            
            %If we have the specified attribute
            if isfield(obj.PropertyGroups(tdx), ['Prop', nam])
                
                %Get the values
                b = obj.PropertyGroups(tdx).(['Prop', nam]);
                
                %Empties mapped to default value
                [b{cellfun(@isempty, b)}] = deal(def);
                
                %Allowing for any of them to be function handles rather than logicals
                bFH = cellfun(@(x)isa(x, 'function_handle'), b);
                b(bFH) = cellfun(@feval, b(bFH), 'UniformOutput', false);
                b = [b{:}];
                
            else
                
                %For backward compat
                b = repmat(def, size(obj.PropertyGroups(tdx).Properties));
                
            end                
            
        end
        
        function setPropertyAttribute(obj, tit, nam, varargin)
            
            %Special cases (TO CONSIDER: good idea ?? maybe need separate method)
            if nargin == 3 && strcmpi(tit, 'NeverLocked')
                
                %Add to list
                obj.PropertiesNeverLocked{end+1} = nam;
                
            elseif nargin == 3 && strcmpi(tit, 'AlwaysLocked')
                
                %Add to list
                obj.PropertiesAlwaysLocked{end+1} = nam;
                
            else
                
                %Does this property group already exist ?
                tdx = find(strcmp(tit, {obj.PropertyGroups.Title}));
                
                %To set anything, it MUST
                assert(~isempty(tdx), ['property group ''', tit, ''' not found']);
                
                %Caller can supply a single name, a list of names, or nothing
                % (which implies all names)
                if isempty(nam)
                    
                    %Go with all
                    ndx = 1:numel(obj.PropertyGroups(tdx).Properties);
                    
                else
                    
                    %Does this property name already exist within this group ?
                    [b, ndx] = ismember(nam, obj.PropertyGroups(tdx).Properties);
                    
                    %To set anything, it MUST
                    if ischar(nam)
                        assert(b, ['property ''', nam, ''' in group ''', tit, ''' not found']);
                    else
                        assert(all(b), ['properties ''', strjoin(nam(b), ''', '''), ''' in group ''', tit, ''' not found']);
                    end
                    
                end
                
                %Set whatever
                for i = 1:2:numel(varargin)
                    [obj.PropertyGroups(tdx).(['Prop', varargin{1}]){ndx}] = deal(varargin{2});
                end
            
            end
            
        end
        
        function addPropertyGroup(obj, tit, varargin)
            
            %Does this property group already exist ?
            tdx = find(strcmp(tit, {obj.PropertyGroups.Title}));
            if isempty(tdx)
                
                %No - extend with new group
                tdx = numel(obj.PropertyGroups) + 1;
                obj.PropertyGroups(tdx).Title = tit;
                
                %By default it is Visible and Enabled (TODO: and Exportable ??)
                obj.PropertyGroups(tdx).Visible = true;
                obj.PropertyGroups(tdx).Enable = true;
                % obj.PropertyGroups(tdx).Export = true;
                
            end
            
            %Caller may optionally provide a function handle, used to modify creation of the view of this group
            if ~isempty(varargin) && isa(varargin{1}, 'function_handle')
                
                %Make a note
                obj.PropertyGroups(tdx).CreateFcn = varargin{1};
                
                %Remove from list
                varargin(1) = [];
                
            end
            
            %Caller may specify properties in short form, with just property name and description,
            % or long form with tooltip string and popup options too.
            %
            %Can unambiguously tell one from the other because popup options will either be [],
            % if the property isn't a popup, or a cellstr of options, or a function handle
            % using which options are evaluated - none of which are char, so
            if iscellstr(varargin)
                
                %Short form.  %How many properties are we adding ?
                assert(rem(numel(varargin), 2) == 0, 'arguments must be supplied in groups of TWO');
                n = numel(varargin) / 2;
                
                %Allow for possibility that we've ALREADY added these details
                if ~isempty(obj.PropertyGroups(tdx).Properties)
                    
                    %Look for duplicates
                    [b, ndx] = ismember(varargin(1:2:end), obj.PropertyGroups(tdx).Properties);
                    if any(b)
                        
                        %Check that other details are same too
                        [bb, ddx] = ismember(varargin(2:2:end), obj.PropertyGroups(tdx).Description);
                        
                        %If all the same
                        if isequal(b, bb) && isequal(ndx, ddx)
                            
                            %COULD be a mistake by the user, or a result of multiple inheritance,
                            % TODO: consider whether a warning should be posted ?
                            
                            %Eliminate duplicate(s)
                            varargin([ndx, ndx + 1]) = [];
                            n = n - numel(ndx);
                            
                        end
                        
                    end
                
                end
                
                %Treat each pair of further arguments as Name and Description
                obj.PropertyGroups(tdx).Properties(end+1:end+n) = varargin(1:2:end);
                obj.PropertyGroups(tdx).Description(end+1:end+n) = varargin(2:2:end);
                
                %Nothing extras for tooltip or popup options
                obj.PropertyGroups(tdx).TooltipString(end+1:end+n) = repmat({''}, 1, n);
                obj.PropertyGroups(tdx).PopupOptions(end+1:end+n) = cell(1, n);
            
            else
                
                %Long form.  %How many properties are we adding ?
                assert(rem(numel(varargin), 4) == 0, 'arguments must be supplied in groups of FOUR');
                n = numel(varargin) / 4;
                
                %Allow for possibility that we've ALREADY added these details
                if ~isempty(obj.PropertyGroups(tdx).Properties)
                    
                    %Look for duplicates
                    [b, ndx] = ismember(varargin(1:4:end), obj.PropertyGroups(tdx).Properties);
                    if any(b)
                        
                        %Check that other details are same too
                        [bb, ddx] = ismember(varargin(2:4:end), obj.PropertyGroups(tdx).Description);
                        
                        %TODO: compare tooltip string etc too ??
                        
                        %If all the same
                        if isequal(b, bb) && isequal(ndx, ddx)
                            
                            %COULD be a mistake by the user, or a result of multiple inheritance,
                            % TODO: consider whether a warning should be posted ?
                            
                            %Eliminate duplicate(s)
                            varargin([ndx, ndx + 1, ndx + 2, ndx + 3]) = [];
                            n = n - numel(ndx);
                            
                        end
                        
                    end
                    
                end
                
                %Treat each group of 4 further arguments as Name, Description, Tooltip and Popup
                obj.PropertyGroups(tdx).Properties(end+1:end+n) = varargin(1:4:end);
                obj.PropertyGroups(tdx).Description(end+1:end+n) = varargin(2:4:end);
                obj.PropertyGroups(tdx).TooltipString(end+1:end+n) = varargin(3:4:end);
                obj.PropertyGroups(tdx).PopupOptions(end+1:end+n) = varargin(4:4:end);
            
            end
            
            %Placeholders for fine control of enable, visible and export states, of individual properties
            obj.PropertyGroups(tdx).PropEnable = cell(size(obj.PropertyGroups(tdx).Properties));
            obj.PropertyGroups(tdx).PropVisible = cell(size(obj.PropertyGroups(tdx).Properties));
            obj.PropertyGroups(tdx).PropExport = cell(size(obj.PropertyGroups(tdx).Properties));
            
            %And any extras
            obj.PropertyGroups(tdx).PropExtras = cell(size(obj.PropertyGroups(tdx).Properties));
            
        end
        
        function ht = propedit(obj, hPanel, varargin)
            
            %Only applicable to single instance (or nothing)
            assert(numel(obj) <= 1, 'only applicable to a single instance');
            
            %Has panel been provided into which properties are displayed ?
            if nargin < 2 || isempty(hPanel)
                
                %No - create from scratch
                [hPanel, hButton] = i_panel(obj);
                
                %Treat as modal
                bWait = true; 
            else
                
                %No need to wait
                bWait = false;
                hButton = [];
                
            end
            
            %Set title of panel - use a helper function
            i_refreshTitle;
            
            %Where the helper function is defined as follows
            function i_refreshTitle(varargin)
                
                %Anything to do ?
                if ishandle(hPanel)
                    
                    %Refresh with dialog title (with 'short' option)
                    set(hPanel, 'Title', dlgtitle(obj, true));
                    
                end
                
            end
            
            %Which makes it easy to listen for relevant changes
            % (TODO: Any way of not having to know which properties are used in dlgtitle ??)
            addlisteners(obj, {'Name', 'Type', 'Description'}, @i_refreshTitle);
            
            %Storage for new values of properties (only used if bWait == true)
            prop_and_value = {};
            
            %Background colours of uicontrol to indicate NoChange, and Change
            colNoChange = [1, 1, 1];
            colChanged = [0.2, 1, 1];
            
            %Apply and OK buttons (if present) are, initially, disabled
            set(hButton, 'Enable', 'off');
            
            %Display properties
            ht = i_propedit(obj, hPanel);
            
            %Add context, if applicable
            if isa(obj, 'mvc.mixin.Contextable')
                context(obj, ht.Children);
            end
            
            %Wait till close ?
            if bWait
                
                %The panel was created here, so treat it as modal and wait till user closes
                uiwait(ancestor(hPanel, 'Figure'));
                
            end
            
            %Where helper functions are as follows
            function hTab = i_propedit(obj, hPanel)
                
                %Properties themselves are displayed in tab pages
                if isempty(hPanel.Children)
                    
                    %Insert a scrolling panel, to help manage scale
                    hScrollingPanel = uix.ScrollingPanel('Parent', hPanel);

                    %Create
                    hTab = uiextras.TabPanel('Parent', hScrollingPanel, ...
                        'Padding', 3, ...
                        'TabWidth', 100);
                    
                    %Set width of scroll panel, using helper function
                    i_panelSizeChanged(hTab, [], hScrollingPanel, 'Width');
                    
                    %And scrolling panel needs to respond to change in size
                    hTab.SizeChangedFcn = {@i_panelSizeChanged, hScrollingPanel, 'Width'};
                    
                    %Default selection
                    sel = 1;
                    
                else
                    
                    %Reuse
                    hTab = hPanel.Children.Children;
                    
                    %Make a note of selection
                    sel = hTab.Selection;
                    
                    %But rebuild from scratch (expensive ?)
                    delete(hTab.Children);
                    
                end
                
                %In case called with no object
                if isempty(obj)
                    
                    %That's all we can do
                    return;
                    
                end
                
                %There could be nothing else to do
                if isempty(obj.PropertyGroups)
                    return;
                end
                
                %Get Visible states
                bVis = {obj.PropertyGroups.Visible};
                
                %Allowing for any of them to be function handles rather than logicals
                bFH = cellfun(@(x)isa(x, 'function_handle'), bVis); 
                bVis(bFH) = cellfun(@feval, bVis(bFH), 'UniformOutput', false);
                bVis = [bVis{:}];
                if ~any(bVis)
                    
                    %Nothing to do
                    return;
                    
                end
                
                %Ditto enables
                bEn = {obj.PropertyGroups.Enable};
                bFH = cellfun(@(x)isa(x, 'function_handle'), bEn); 
                bEn(bFH) = cellfun(@feval, bEn(bFH), 'UniformOutput', false);
                bEn = [bEn{:}];
                
                %For each (visible) property group
                for i = find(bVis)
                    
                    %Check custom create field
                    if isempty(obj.PropertyGroups(i).CreateFcn)
                        
                        %Nothing specified - so just crack on, displaying this property group directly in hTab
                        hPar = hTab;
                        
                    elseif isa(obj.PropertyGroups(i).CreateFcn, 'function_handle')
                        
                        %Call it - allowing the callee to modify where this property group is displayed
                        hPar = obj.PropertyGroups(i).CreateFcn(obj, hTab);
                        
                        %If custom function returns nothing
                        if isempty(hPar)
                            
                            %Then we're done, for this property group
                            continue;
                            
                        end
                        
                    else
                        error('bad value for custom create function');
                    end
                    
                    %Allow fine control over visible state of each property this group
                    bVisP = getPropertyAttributes(obj, i, 'Visible', true); % i_propertyState(obj, i, 'Visible', true);
                    
                    %Anything to do ?
                    if ~any(bVisP) || i_omitPropertyGroup(obj, i)
                        
                        %No
                        continue;
                        
                    elseif i_omitPropertyNames(obj, i)
                        
                        %Add a simple panel
                        hGroup = uiextras.Panel('Parent', hPar, 'Padding', 3);
                        
                        %No layout required
                        bLayout = false;
                        
                    else

                        %Add a scrolling panel in which to display content
                        hScrollPanel = uix.ScrollingPanel('Parent', hPar);
                        
                        %Add a panel, in the form of a grid, in which to display multiple properties
                        hGroup = uiextras.GridFlex('Parent', hScrollPanel, 'Spacing', 3, 'Padding', 3);
                        
                        %Use property names as labels
                        nam = obj.PropertyGroups(i).Properties;
                        
                        %Unless description is set
                        b = ~cellfun(@isempty, obj.PropertyGroups(i).Description);
                        nam(b) = obj.PropertyGroups(i).Description(b);
                        
                        %Down-select i.a.w. individual visible flags
                        nam = nam(bVisP);
                        
                        %Label controls
                        cellfun(@(x)uicontrol('Parent', hGroup, 'Style', 'Text', 'HitTest', 'off', 'String', x), nam);
                        
                        %Layout required
                        bLayout = true;
                        
                    end
                    
                    %Set tab title and enable state
                    hTab.TabTitles{end} = obj.PropertyGroups(i).Title;
                    hTab.TabEnables{end} = mvc.mixin.UiTools.bool2offon(bEn(i));
                    
                    %Allow fine control over enable state of each property this group
                    bEnP = getPropertyAttributes(obj, i, 'Enable', true); % i_propertyState(obj, i, 'Enable', true);
                    
                    %Add property values
                    [~, dy] = cellfun(@(x,y,z,b,ext)i_uicontrol(x, y, z, b, ext, 'Parent', hGroup), ...
                        obj.PropertyGroups(i).Properties(bVisP), ...
                        obj.PropertyGroups(i).TooltipString(bVisP), ...
                        obj.PropertyGroups(i).PopupOptions(bVisP), ...
                        mvc.mixin.UiTools.bool2offon({bEnP(bVisP)}), ...
                        obj.PropertyGroups(i).PropExtras(bVisP), ...
                        'UniformOutput', false);
                    
                    %Set layout, if required
                    if bLayout
                        
                        %Set widths and heights of controls in the grid
                        hGroup.Widths = [-1, -4];
                        hGroup.Heights = [dy{:}];
                        
                        %Set height of scroll panel, using helper function
                        i_panelSizeChanged(hPanel, [], hScrollPanel, 'Height');
                        
                        %And dynamically link the scroll panel height to the hPanel height
                        hPanel.SizeChangedFcn = @i_panelSizeChanged;
                        
                    end
                    
                end
                
                %Selected tab ?
                sel = max(1, min(sel, numel(hTab.Children)));
                hTab.Selection = sel;
                
                %In case number of tabs exceeds availabnle width
                hcm = uicontextmenu('Parent', ancestor(hTab, 'Figure'));
                arrayfun(@(x)uimenu(hcm, 'Label', hTab.TabTitles{x}, 'Callback', {@i_tabSelect, hTab, x}), 1:numel(hTab.TabTitles));
                function i_tabSelect(~, ~, hTab, i)
                    hTab.Selection = i;
                end
                hTab.TabContextMenus = repmat({hcm}, 1, numel(hTab.TabTitles));
                
            end
                
            %Where the scroll-panel size is maintained thus
            function i_panelSizeChanged(hPanel, ~, hScrollPanel, wh)
                
                %Enough inputs ?
                if nargin < 3
                    
                    %No
                    return;
                    
                end
                
                %Width or height ?
                switch lower(wh)
                    
                    case 'width'
                        
                        %Set width of scroll panel - Edit 15/06/2018 - C.Szczyglowski
                        w(1) = hPanel.Position(3);
                        if  isa(hScrollPanel.Children, 'uix.TabPanel')
                            hTab = hScrollPanel.Children;
                            w(2) =((numel(hTab.TabTitles) + 1) * hTab.TabWidth) + (2 * hTab.Padding);
                        else
                            w(2) = 0;
                        end
                        
                        %Take the maximum of the panel width or the tab
                        %panel width
                        hScrollPanel.Widths = max(w);
                        
                    case 'height'
                        
                        %Set height of scroll panel
                        hScrollPanel.Heights = hPanel.Position(4);
                        
                    otherwise
                        error('bad option');
                end
                
            end
            
            function b = i_propertyState(obj, i, nam, def)
                
                %If the named field exists (prefixed with ''Prop'', as it applies to each property individually)
                if isfield(obj.PropertyGroups(i), ['Prop', nam])
                    
                    %Get the values
                    b = obj.PropertyGroups(i).(['Prop', nam]);
                    
                    %Empties mapped to default value
                    [b{cellfun(@isempty, b)}] = deal(def);
                    
                    %Allowing for any of them to be function handles rather than logicals
                    bFH = cellfun(@(x)isa(x, 'function_handle'), b);
                    b(bFH) = cellfun(@feval, b(bFH), 'UniformOutput', false);
                    b = [b{:}];
                    
                else
                    
                    %For backward compat
                    b = rempat(def, size(obj.PropertyGroups(i).Properties));
                    
                end
                
            end
            
            function b = i_omitPropertyGroup(obj, i)
            
                %Assume that we do NOT want to omit this group
                b = false;
                
                %For each property
                for j = 1:numel(obj.PropertyGroups(i).Properties)
                    
                    %If property if not empty
                    if ~isempty(obj.(obj.PropertyGroups(i).Properties{j}))
                        
                        %Then we have something to display
                        return;

                    elseif isEditable(obj, obj.PropertyGroups(i).Properties{j})
                        
                        %Then we have nothing at the moment, but user could edit
                        return;
                        
                    end
                    
                end
                
                %IF we get this far, there's nothing to display so we might as well skip
                b = true;
                
            end
            
            function b = i_omitPropertyNames(obj, i)
                %
                % Determine whether it is sensible to include property names, as well as values, in dialog.
                
                %How many properties are we displaying ?
                if numel(obj.PropertyGroups(i).Properties) > 1
                    
                    %More than one, so it makes sense to NOT omit names
                    b = false;
                    
                else
                   
                    %What have we got ?
                    val = obj.(obj.PropertyGroups(i).Properties{1});
                    
                    %If a table, a time-table or a cellstr, it makes sense to omit names
                    b = iscell(val) || istable(val) || (~verLessThan('matlab', '9.3') && istimetable(val));
                    
                end
            end
            
            function [hPanel, hButton] = i_panel(obj)
                
                %Create helper dialog
                hf = dialog('Visible', 'off', 'Name', dlgtitle(obj), ...
                    'CloseRequestFcn', @i_close, ...
                    'Resize', 'on'); %Edit 15/06/2018 - C.Szczyglowski (Airbus wanted to be able to resize the window so they could view all tabs)
                
                %Split vertically
                hv = uiextras.VBox('Parent', hf);
                
                %Display properties in a panel
                hPanel = uiextras.Panel('Parent', hv, 'Padding', 3);
                                               
                %Make room for Apply, OK and cancel buttons
                hb = uiextras.HButtonBox('Parent', hv, 'Spacing', 3, 'Padding', 3);
                hButton(1) = uicontrol('Parent', hb, 'String', 'Apply', 'callback', @i_apply, 'TooltipString', 'Applies pending changes to property values');
                hButton(2) = uicontrol('Parent', hb, 'String', 'OK', 'callback', @i_ok, 'TooltipString', 'Close dialog, after applying pending changes to property values');
                uicontrol('Parent', hb, 'String', 'Close', 'callback', @i_close, 'TooltipString', 'Close dialog, without applying pending changes to property values');
                
                %Set layout
                hv.Heights = [-1, 50];
                
            end
            
            function i_ok(varargin)
                
                %Apply any changes
                i_apply;
                
                %Close
                i_close;
                
            end
            
            function i_apply(varargin)
                
                %Anything to apply ?
                if isempty(prop_and_value)
                    
                    %No
                    return;
                    
                end
                
                %Apply changes
                set(obj, prop_and_value{:});
                
                %Clear local store
                prop_and_value = {};
                
                %There's now nothing to apply
                set(hButton, 'Enable', 'off');
                
                %And something has changed
                i_somethingHasChanged;
                
            end
            
            function i_somethingHasChanged(varargin)
                                
                %If serializable (either directly in this object, or indirectly via parent(s))
                if isa(obj, 'mvc.mixin.Serializable')
                    
                    %Make a note that something has changed
                    obj.ChangedSinceLastSave = true;
                    
                elseif isa(obj, 'mvc.mixin.Collectable') && isa(obj.Root, 'mvc.mixin.Serializable')
                    
                    %Make a note that something has changed
                    obj.Root.ChangedSinceLastSave = true;
                    
                end
                
            end
            
            function i_close(varargin)
                
                %Simply close
                delete(ancestor(hPanel, 'Figure'));
                
            end
            
            function [hc, dy] = i_uicontrol(prp, tip, opt, en, extras, varargin)
                %
                % Creates a uicontrol suitable for display of property prp
                
                %Anything to add ?
                if isempty(extras)
                    
                    %No
                    
                elseif iscell(extras)
                    
                    %Append to varargin
                    varargin(end+1:end+numel(extras)) = extras;
                    
                else
                    error('bad input');
                end
                
                %Modify default enable state of this control ?
                if ismember(prp, obj.PropertiesNeverLocked)
                    
                    %Ignore state of the "PropertiesLocked" flag
                    en = 'on';
                    
                elseif ismember(prp, obj.PropertiesAlwaysLocked)
                    
                    %Ignore state of the "PropertiesLocked" flag
                    en = 'off';
                    
                elseif obj.PropertiesLocked
                    
                    %Set to off, regardless
                    en = 'off'; % mvc.mixin.UiTools.bool2offon(~obj.PropertiesLocked);
                    
                end
                
                %What have we got ?
                [val, dy, args] = i_val2str(obj.(prp), en, tip, opt);
                
                %Tables are 'special' TODO: rationalise this, e.g. merge into i_val2str
                % Need version check because 2015b doesn't support timetable
                if istable(val) || (~verLessThan('matlab', '9.3') && istimetable(val))
                    
                    %Make a note of timestamps
                    if (~verLessThan('matlab', '9.3') && istimetable(val))
                        t = val.when;
                    else
                        t = [];
                    end
                    
                    %Column format may be specified in userdata
                    if iscell(val.Properties.UserData)
                        args = {'ColumnFormat', val.Properties.UserData};
                    else
                        args = {};
                    end
                    
                    %Can't display directly
                    nam = val.Properties.VariableNames;
                    val = table2cell(val);
                    val = obj.tabulatable(val);
                    
                    %Display in uitable
                    hc = uitable('Data', val, ...
                        'ColumnName', nam, ...
                        'Tag', prp, ...
                        'Enable', en, ...
                        'TooltipString', char(tip), ...
                        args{:}, varargin{:});

                    %If it was a timetable
                    if ~isempty(t)
                    
                        %Set row-names
                        hc.RowName = char(t);
                        
                    end
                
                elseif ~strcmpi('popup', args{find(strcmpi('style', args(1:2:end)))*2}) && iscellstr(val) && size(val,2) > 1
                    
                    %Can we display in tabular form ?
                    hc = uitable('Data', val, ...
                        'Tag', prp, ...
                        'Enable', en, ...
                        'TooltipString', tip, ...
                        varargin{:});
                    
                else
                    
                    %This is a VERY special case that should be catered for more generically
                    if strcmp(prp, 'StructureLocked') && isa(obj, 'mvc.mixin.Collectable') && obj.StructureLockedByParent
                        
                        %Treat as NOT editable
                        varargin(end+1:end+4) = {'Value', true, 'Enable', 'off'};
                        
                    elseif isEditable(obj, prp)
                        
                        %Set action on change
                        varargin(end+1:end+2) = {'Callback', {@i_edit, prp}};
                        
                    else
                        
                        %No
                        varargin(end+1:end+2) = {'Enable', 'off'};
                        
                    end
                    
                    %Tag the control with name of property it is displaying,
                    % and set the initial background colour of the control
                    args(end+1:end+4) = {'Tag', prp, 'BackgroundColor', colNoChange};
                    
                    %Do it
                    hc = uicontrol(args{:}, varargin{:});
                    
                end
                
                %Listen for changes made elsewhere
                addlisteners(obj, prp, @i_refreshContent);
                
            end
            
            function [str, dy, args] = i_val2str(val, en, tip, opt)
                
                %Good height ?  By default, go with good height (pixels) for a single character
                dy = 25;
                
                %May need multiline
                bMultiLine = false;
                
                %Default enable state ?
                if nargin < 2
                    en = 'on';
                end
                
                %Added later
                if nargin < 3 || isempty(tip)
                    argx = {};
                else
                    argx = {'TooltipString', tip};
                end
                    
                %Default style is
                if nargin < 4 || isempty(opt)

                    %Edit box
                    sty = 'edit';
                    
                else
                    
                    %Unless popup options have been specified
                    sty = 'popup';
                    
                    %And the popup strings are
                    if iscellstr(opt)
                        
                        %Specified explicitly
                        str = opt;
                        
                    elseif isa(opt, 'function_handle')
                        
                        %Specified indirectly
                        str = opt();
                        
                    else
                        error('bad input');
                    end
                    
                    %The value must be one of the options
                    vdx = find(strcmp(val, str));
                    
                    %Safety net
                    % assert(~isempty(vdx), 'value does not match any of the options');
                    if isempty(vdx)
                        
                        %Don't bomb, might be recoverable
                        warning(['value ''', val, ''' does not match any of the available options ''', strjoin(str, ''', ''')]);
    
                        %Better than nothing
                        val = 1;
                        
                    else
                        val = vdx;
                    end

                end
                
                %Is it easily displayable ?
                if isnumeric(val) && strcmp(sty, 'popup')
                    
                    %Include value in additional args
                    argx(end+1:end+2) = {'Value', val};
                    
                elseif ischar(val)
                    
                    %No conversion required
                    str = val;
                    
                    %Allow for multi-line text
                    if size(val,1) > 1
                        
                        %Adjust properties of uicontrol
                        bMultiLine = true;
                        
                        %Allow height to stretch
                        dy = -1;
                        
                    end
                
                elseif isstring(val)
                    
                    %A bit experimental, needs more work
                    str = strjoin(val, ', ');
                    
                    %So turn editing off
                    en = 'off';
                    
                elseif iscellstr(val)
                    
                    %No conversion required
                    str = val;
                    
                    %Adjust properties of uicontrol
                    bMultiLine = true;
                    
                    %Allow height to stretch
                    dy = -1;
                
                elseif islogical(val)
                    
                    %Adjust style accordingly
                    sty = 'checkbox';
                    
                    %String not required (because label created elsewhere)
                    str = '';
                    
                    %Include value in additional args
                    argx(end+1:end+2) = {'Value', val};
                    
                elseif isnumeric(val)
                    
                    %Care may be required over how to show the value
                    sz = size(val);
                    if numel(sz) > 2
                        
                        %Far too difficult !
                        str = ['-', class(val), ' with dimensions ', mat2str(sz), ' -'];
                        
                    elseif any(sz == 0)
                        
                        %No conversion required
                        str = val;
                        
                    elseif any(sz == 1)
                        
                        %Single-line is fine
                        str = mat2str(val);
                        
                    else
                        
                        %No conversion required
                        str = val;
                        
                        %Go with multi-line
                        bMultiLine = true;
                        
                        %Allow height to stretch
                        dy = dy .* (size(val,1) + 1);
                        
                    end
                    
                elseif isdatetime(val)
                    
                    %Display as string
                    str = char(val);
                    
                elseif istable(val) || (~verLessThan('matlab', '9.3') && istimetable(val))
                    
                    %No conversion required here (handled differently elsewhere)
                    str = val;
                    
                    %Adjust properties of uicontrol
                    bMultiLine = true;
                    
                    %Allow height to stretch
                    dy = -1;
                
                elseif iscell(val)
                    
                    %Recurse to create a cellstr
                    str = cellfun(@i_val2str, val, 'UniformOutput', false);
                    
                    %Adjust properties of uicontrol
                    bMultiLine = true;
                    
                    %Allow height to stretch
                    dy = -1;
                
                elseif isa(val, 'mvc.mixin.Nameable')
                    
                    %Not editable, but displayable
                    str = dlgtitle(val);
                    en = 'off';
                    
                else
                    
                    %Not editable or easily displayable
                    str = ['-', class(val), '-'];
                    en = 'off';
                    
                end
                
                %Extras
                args = {'Style', sty, 'HorizontalAlignment', 'left', 'String', str, 'Enable', en};
                
                %Allow for multi-line if necessary
                if bMultiLine
                    args(end+1:end+4) = {'Min', 0, 'Max', 2};
                end
                
                %And any extras
                args(end+1:end+numel(argx)) = argx;
                
            end
            
            function i_refreshContent(hc, ~)
                
                %Which property has changed ?
                if isa(hc, 'meta.property')
                    prp = hc.Name;
                else
                    prp = hc.Source.Name;
                end
                
                %Which uicontrol needs to be updated ?
                hu = findobj(hPanel, 'Tag', prp);
                
                %New value (as string)
                val = i_val2str(obj.(prp));
                
                %Update content
                if isa(hu, 'matlab.ui.control.UIControl')
                    
                    %Allow for multiple controls (good idea ??)
                    for i = 1:numel(hu)
                        if strcmpi('popupmenu', hu(i).Style)
                            set(hu(i), 'Value', find(strcmp(val, hu(i).String)), 'BackgroundColor', colNoChange);
                        else
                            set(hu(i), 'String', val, 'BackgroundColor', colNoChange);
                        end
                    end
                    
                elseif isa(hu, 'matlab.ui.control.Table')
                    
                    %If original data are in table form
                    if istable(val)
                        
                        %Can't display directly in uitable
                        nam = val.Properties.VariableNames;
                        val = table2cell(val);
                        val = obj.tabulatable(val);
                        set(hu, 'Data', val, 'ColumnName', nam);
                    
                    elseif (~verLessThan('matlab', '9.3') && istimetable(val))
                        
                        %Can't display directly in uitable
                        nam = val.Properties.VariableNames;
                        val = table2cell(val);
                        val = obj.tabulatable(val);
                        set(hu, 'Data', val, 'ColumnName', nam);
                        
                    else
                        
                        %Just apply
                        set(hu, 'Data', val);
                        
                    end
                    
                end
                
            end
            
            function i_edit(hc, ~, prp)
                
                %Careful
                try
                    
                    %Is it legal to change content ?
                    if ismember(prp, obj.PropertiesNeverLocked)
                        
                        %Yes
                        
                    elseif ismember(prp, obj.PropertiesAlwaysLocked)
                        
                        %No
                        error(['property ''', prp, ''' of this object is LOCKED for editing']);
                        
                    elseif obj.PropertiesLocked
                        
                        %No
                        error('properties of this object are LOCKED for editing');
                        
                    end

                    %Check for specific data types
                    if strcmpi('popupmenu', get(hc, 'Style'))
                        
                        %Get new value, as popupstr
                        newval = popupstr(hc);
                        
                    elseif islogical(obj.(prp))
                        
                        %Get new value, as logical
                        newval = logical(get(hc, 'Value'));
                        
                    elseif isnumeric(obj.(prp))
                        
                        %Get new value, as numeric
                        newval = str2num(get(hc, 'String'));
                        
                        %Conserve data type
                        fcn = str2func(class(obj.(prp)));
                        newval = fcn(newval);
                    
                    else
                        
                        %Just get new value as string
                        newval = get(hc, 'String');
                        
                    end
                    
                    %Are we waiting for user to hit OK, or Apply ?
                    if bWait
                        
                        %Locate property name in existing store of new values
                        pdx = find(strcmp(prp, prop_and_value(1:2:end))) * 2 - 1;
                        if isempty(pdx)
                            pdx = numel(prop_and_value) + 1;
                            prop_and_value{pdx} = prp;
                        end
                        
                        %Make a note of new values (without immediate apply)
                        prop_and_value{pdx + 1} = newval;
                        
                        %Tweak the appearance of the uicontrol
                        set(hc, 'BackgroundColor', colChanged);
                        
                        %There's now something to apply
                        set(hButton, 'Enable', 'on');
                        
                    else
                        
                        %Safety net
                        
                        %Just get on with it
                        obj.(prp) = newval;
                        
                        %And something has changed
                        i_somethingHasChanged;
                        
                    end
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj, err, 'modal'));
                    
                    %Revert content
                    switch get(hc, 'Style')
                        
                        case 'popupmenu'
                            
                            %Care required
                            idx = find(strcmp(obj.(prp), get(hc, 'String')));
                            if ~isempty(idx)
                                set(hc, 'Value', idx);
                            else
                                warning('unable to match whilst trying to revert popup value');
                            end
                            
                        otherwise
                            
                            %Easy, I hope
                            set(hc, 'String', i_val2str(obj.(prp)));
        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    methods (Sealed) % Helper functions acting on a heterogrenous array of Nameables
        
        function makeUniqueNames(obj)
        
            %What have we got now ?            
            nam = {obj.Name};
            
            %Exclude any default names
            def = {obj.DefaultName};
            
            %Make them all unique
            nam = matlab.lang.makeUniqueStrings(nam, def);
            
            %Putback
            [obj.Name] = deal(nam{:});
            
        end
        
        function strDiff = compare(obj, ref, varargin)
            %
            % Compares OBJ with REF, a reference object
            
            %Create a list of differences (descriptive)
            strDiff = {};
            
            %What function to use to compare ?
            fcn = @isequal;

            %Caller MAY override
            b = cellfun(@(x)ischar(x) && strcmp(x, '-ignorenan'), varargin);
            if any(b)
                fcn = @isequaln;
                varargin(b) = [];
            end
            
            %Start with some simple tests
            if ~strcmp(class(obj), class(ref))
                strDiff{end+1} = 'classes differ';
            elseif numel(obj) ~= numel(ref)
                strDiff{end+1} = 'numbers differ';
            elseif ~isequal(size(obj), size(ref))
                strDiff{end+1} = 'sizes differ';
            elseif all(arrayfun(@eq, obj, ref)) % can't just do "obj == ref" because eq not sealed
                
                %No point going any further, there can't be any differences
                strDiff = {};
                
            elseif fcn(obj, ref)
                
                %No point going any further, there aren't be any differences
                strDiff = {};
              
            else
                
                %Consider each object in turn
                for i = 1:numel(obj)
                    
                    %Get list of properties we need to test
                    prp = properties(obj(i));
                    
                    %If obj is Collectable
                    if isa(obj(i), 'mvc.mixin.Collectable')
                        
                        %Then do NOT bother with Root and Parent properties (else we recurse endlessly),
                        % and StructureLocked (incl by Parent) (because not that bothered (??))
                        prp = setdiff(prp, {'Root', 'Parent', 'StructureLocked', 'StructureLockedByParent'}, 'stable');
                        
                    end
                    
                    %Get values from this and that
                    valObj = get(obj(i), prp);
                    valRef = get(ref(i), prp);
                    
                    %Quick check property by property
                    bSame = cellfun(fcn, valObj, valRef);
                    
                    %Allow for comparisons of object arrays
                    bSame = all(bSame, 1);
                    
                    %So the answer is
                    strDiff = cellfun(@(x)[x, ' differ'], prp(~bSame), 'UniformOutput', false);
                    
                    %BUT were any of the properties that differed also Nameables ?
                    if any(~bSame)
                        
                        %If so, need to recursively call this function
                        bRecurse = cellfun(@(x,y)isa(x, 'mvc.mixin.Nameable') && isa(y, 'mvc.mixin.Nameable'), valObj(~bSame), valRef(~bSame));
                        if any(bRecurse)
                            
                            %Pass it on
                            idx = find(~bSame);
                            strRecurse = cell(numel(idx),1);
                            strRecurse(bRecurse) = cellfun(@compare, valObj(idx(bRecurse)), valRef(idx(bRecurse)), 'UniformOutput', false);
                            
                            %Some of these may have turned out NOT to differ
                            jdx = find(bRecurse(:) & cellfun(@isempty, strRecurse));
                            if ~isempty(jdx)
                                
                                %Backtrack
                                strDiff(jdx) = [];
                                bRecurse(jdx) = []; %% false;
                                strRecurse(jdx) = [];
                                kdx = find(~bSame);
                                bSame(kdx(jdx)) = true;
                                
                            end
                            
                            %Append details
                            if any(bRecurse)
                                if numel(obj) == 1
                                    prefix = ': ';
                                else
                                    prefix = ['(', num2str(i), '): '];
                                end
                                strDiff(bRecurse) = cellfun(@(x,y)[x, prefix, strjoin(y, ', ')], strDiff(bRecurse), strRecurse(bRecurse), 'UniformOutput', false);
                            end
                            
                        end
                        
                    end
                
                end
                
            end
            
            %If caller does not want the answer back
            if nargout == 0
                
                %Tell us about it
                if isempty(strDiff)
                    disp('no differences');
                else
                    cellfun(@disp, strDiff);
                end
                
            end
            
        end
        
        function bIsEditable = isEditable(obj, varargin)
            %
            % Returns array bIsProp of size numel(obj) by numel(varargin),
            %  indicating whether each name listed in varargin is an editable property of each element of object array
            
            %May optionally consider settable, but un-gettable, Dependent properties, as uneditable
            [bExcludeUngettables, varargin] = checkOption('-excludeUngettables', false, varargin{:});
            
            %Start here
            bIsEditable = isProperty(obj, varargin{:});
            
            %For valid properties, further test of editability
            for i = 1:numel(obj)
                
                %Get help from metaclass (TODO: do this once only for each class in obj array)
                mc = metaclass(obj(i));
                pl = mc.PropertyList;
                
                %Locate required properties in metadata
                [b, pdx] = ismember(varargin(bIsEditable(i,:)), {pl.Name});
                
                %TODO - Cater for dynamic properties as these do not appear
                %in the PropertyList field of the MetaClass object.
                bNotEd = false(size(b));

                %Handle Dependents differently, depending on option to excludeUngettables
                if bExcludeUngettables
                    bDepNotEd = [pl(pdx(b)).Dependent] & (cellfun(@isempty, {pl(pdx(b)).SetMethod}) | cellfun(@isempty, {pl(pdx(b)).GetMethod}));
                else
                    bDepNotEd = [pl(pdx(b)).Dependent] & cellfun(@isempty, {pl(pdx(b)).SetMethod});
                end
                                
                %NOT editable also if any of the following are satisfied
                bOtherNotEd = [pl(pdx(b)).Constant] | [pl(pdx(b)).Hidden] | ~strcmp('public', {pl(pdx(b)).SetAccess});
                
                %Putback
                bIsEditable(i,b) = ~bDepNotEd & ~bOtherNotEd;
                
            end
            
        end
        
        function bIsProperty = isProperty(obj, varargin)
            %
            % Returns array bIsProp of size numel(obj) by numel(varargin),
            %  indicating whether each name listed in varargin is a property of each element of object array
            
            %Annoyingly, need a special case
            if isempty(varargin)
                
                %Because arrayfun would not return an [n by 0] logical
                bIsProperty = false(numel(obj), 0);
                
            else
                
                %Care required, to allow for obj as heterogenous array
                bIsProperty = arrayfun(@(x)cellfun(@(y)isprop(x,y), varargin), obj, 'UniformOutput', false);
                
                %Expand
                bIsProperty = vertcat(bIsProperty{:});
            
            end
            
        end
        
        function prp = allProperties(obj, varargin)
            %
            % Returns PRP, a list of all accessible properties associated with object array,
            %  derived from the internally stored 'PropertyGroups' detail.
            
            %For each object in turn
            prp = arrayfun(@(x)i_allProperties(x, varargin{:}), obj, 'UniformOutput', false);
            
            %Return a flat list of unique entries (preserving order)
            if ~isempty(prp)
                prp = unique([prp{:}], 'stable');
            end
            
            function prp = i_allProperties(obj, varargin)
                
                %If the 'force' flag is specified
                if nargin == 2 && islogical(varargin{1}) && varargin{1}                    
                    
                    %Get help from metaclass (TODO: do this once only for each class in obj array)
                    mc = metaclass(obj);
                    pl = mc.PropertyList;
                    
                    %NO do not even strip those that are not public - allProperties should mean exactly that
                    %                     %Omit anything that is not publicly settable
                    %                     pl(~strcmp('public', {pl.SetAccess})) = [];
                    %
                    %                     %Other tests may be required to determine what to return to caller
                    %                     %bRet = ~[pl(pdx).Dependent] | [pl(pdx).Constant] | [pl(pdx).Hidden] | ;
                
                    %So the answer is
                    prp = {pl.Name};
                    
                    %If we also have dynamic properties
                    if isa(obj, 'mvc.mixin.Dynamicable')
                        
                        %Include their names in the list
                        prp = [prp, obj.DynamicPropertyNames];
                        
                    end
                    
                else
                    
                    %Return only those properties detailed in property groups
                    prp = arrayfun(@(x)x.Properties, obj.PropertyGroups, 'UniformOutput', false);
                    
                    %Need ability to down-select to those with particular attributes
                    if nargin > 1
                        
                        %Get masks
                        b = arrayfun(@(x)getPropertyAttributes(obj, x, 'Export', true), 1:numel(obj.PropertyGroups), 'UniformOutput', false);

                        %Apply
                        prp = cellfun(@(x,y)x(y), prp, b, 'UniformOutput', false);
                        
                    end
                
                    %Return as flat list
                    prp = [prp{:}];
            
                    %But ensure we always return a cell array
                    if isempty(prp)
                        prp = {};
                    end
                    
                end
                
            end
            
        end
        
        function [val, bIsProperty] = obj2table(obj, varargin)
            %
            % Returns a table VAL containing values of properties PRP of object array OBJ.
            % If PRP not explicitly specified, it defaults to all properties of OBJ.
            %
            % For more details, see obj2cell.
            
            %Start with obj2cell
            [val, prp, bIsProperty] = obj2cell(obj, varargin{:});
            
            %Convert to table
            val = cell2table(val, 'VariableNames', prp);
            
        end
        
        function [val, prp, bIsProperty] = obj2cell(obj, prp, def, varargin)
            %
            % Returns a cell-array VAL containing values of properties PRP of object array OBJ.
            % If PRP not explicitly specified, it defaults to all properties of OBJ.
            %
            % OBJ may be heterogenous, so all of PRP may not be relevant to all of OBJ.
            %  bISPROP returns a map shows which of PRP applies to which of OBJ.
            %  VAL is returned with empties where bISPROP is false, UNLESS default value DEF
            %  is supplied, in which case DEF is assigned to VAL(~bISPROP).
            
            %Which properties ?
            if nargin < 2 || isempty(prp)
                
                %Go with these, by default
                prp = allProperties(obj);
                
                %Ensure a row
                prp = prp(:).';
                
            end
            
            %Start by checking what is actually a property
            bIsProperty = isProperty(obj, prp{:});
            
            %Corresponding values ?  Can't use 'get' directly, because it's not sealed
            % val = get(obj, prp);
            %You might think that this would work, but is not ideal because not all members of obj may have fields in prp
            % val = cellfun(@(x){obj.(x)}, prp, 'UniformOutput', false);
            %This is more flexible
            val = arrayfun(@(x)i_get(x, prp), obj, 'UniformOutput', false);
            
            %Where the helper get function is as follows
            function val = i_get(obj, prp)
                
                %Placeholder
                val = cell(size(prp));
                
                %What is applicable, for this object
                bIsProp = cellfun(@(x)isprop(obj, x), prp);
                
                %Only get what is applicable
                val(bIsProp) = get(obj, prp(bIsProp));
                
            end
            
            %Expand
            val = vertcat(val{:});
            
            %If a default value is provided for content that is not applicable
            if nargin > 2
                
                %Assign as appropriate
                [val{~bIsProperty}] = deal(def);
                
            end
            
        end
        
    end
    
    methods (Sealed) % Implementing Set/Get for heterogeneous arrays
       
        function varargout = set(obj, varargin)
            [varargout{1 : nargout}] = set@matlab.mixin.SetGet(obj, varargin{:});
        end
        
        function varargout = get(obj, varargin)
            [varargout{1 : nargout}] = get@matlab.mixin.SetGet(obj, varargin{:});
        end
        
    end
    
    methods (Static)
        
        function val = val2datetime(val)
            %
            % Whatever we get, cast it to an equivalent datetime
            
            %What have we got ?
            if isdatetime(val)
                
                %No worries
                
            elseif ischar(val)
                
                %Special case
                if strcmp(val, 'now')
                    
                    %Convert
                    val = datetime;
                    
                else
                    
                    %Attempt to convert
                    val = datetime(datevec(val));
                    
                end
                
            elseif isnumeric(val) && isscalar(val)
                
                %Attempt to convert
                val = datetime(datevec(val));
                
            elseif isnumeric(val) && any(numel(val) == [3, 6])
                
                %Attempt to convert
                val = datetime(val);
                
            else
                error(['value of class ''', class(val), ''' cannot be converted to datetime']);
            end
            
        end
        
        function [val, bInPlaceEditable] = tabulatable(val)
            %
            % Content of cell-array are converted into a form that can be written as 'Data' attribute of a uitable
            
            %Look for valid content
            bEmpty = cellfun(@isempty, val);
            bNumericEmpty = cellfun(@(x)isnumeric(x) && numel(x) == 0, val);
            bNumericScalar = cellfun(@(x)isnumeric(x) && numel(x) == 1, val);
            bLogicalScalar = cellfun(@(x)islogical(x) && numel(x) == 1, val);
            bChar = cellfun(@ischar, val);
            bDatetime = cellfun(@isdatetime, val);
            bCellstr = cellfun(@iscellstr, val);
            n = cellfun(@numel, val);
           
            %Numeric arrays treated as a special case
            bNumericArray = cellfun(@(x)isnumeric(x) && numel(x) > 0, val);
            
            %Convert numeric arrays to strings, for display purposes
            val(bNumericArray) = cellfun(@mat2str, val(bNumericArray), 'UniformOutput', false);
            
            %Datetimes can be displayed as strings
            val(bDatetime) = cellfun(@char, val(bDatetime), 'UniformOutput', false);
            
            %What is not valid content ?
            bNotValid = ~(bNumericScalar | bNumericEmpty | bLogicalScalar | bChar | bDatetime | bNumericArray | bCellstr);
            
            %Convert content accordingly
            val(bNotValid) = cellfun(@(x)['[', class(x), ']'], val(bNotValid), 'UniformOutput', false);
            
            %Which columns are in-place editable ?
            bInPlaceEditable = all(bNumericScalar | bLogicalScalar | bChar, 1);
            
            %Cellstr can be treated as a special case - handling them AFTER everything else,
            % because this might result in the cell-array we return being a different size from original table
            
            %Set empty cells to empty strings
            [val{bCellstr & n == 0}] = deal('');

            %Single-line cellstr map to just the single string
            val(bCellstr & n == 1) = cellfun(@(x)x{1}, val(bCellstr & n == 1), 'UniformOutput', false);
            
            %Multi-line cellstr ?
            [r,c] = find(bCellstr & n > 1);
            if numel(unique(c)) == 1
                
                %Expand as appropriate (working backwards, to avoid ovewrite during expansion)
                for i = numel(r):-1:1
                    
                    %Pull out content from line 2 onwards
                    rem = val{r(i), c(i)}(2:end,:);
                    
                    %Retain just first line in table
                    val{r(i), c(i)} = val{r(i), c(i)}{1,1};
                    
                    %Insert
                    val = [val(1:r(i),:); cell(size(rem,1),c(i)-1), rem, cell(size(rem,1), size(rem,2) - c(i)); val(r(i)+1:end,:)];
                    
                end
                
            else
                
                %Display first line(s) only - better than nothing (??)
                val(bCellstr & n > 1) = cellfun(@(x)[x{1}, '...'], val(bCellstr & n > 1), 'UniformOutput', false);
                
            end
            
        end
        
    end
    
end
