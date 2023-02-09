classdef (ConstructOnLoad) Collectable < matlab.mixin.SetGet & matlab.mixin.Copyable
    %
    % Collectable provides functionality associated with being part of a list, or collection.

    %EDIT : C.Szczyglowski 24/07/2018
    properties
        %Determines whether the object is currently being imported as part
        %of a larger import process
        ObjBeingImported = false;
    end
    
    properties (AbortSet) % (Access = protected)
        
        %The parent of this object
        Parent;
        
        %The children of this object
        Children;
        
        %Does this class represent a LEAF node in the collection ?
        IsLeafNode = false;
        
        %Are all items uniquely named ?
        AllNamesUnique = true;
        
        %Is this node expanded in Tree ? (TODO: store this in each relevant view ?)
        IsExpanded = true;
        
    end
   
    properties (AbortSet, SetObservable)
        
        %Is structure locked (i.e. do not permit add/remove/move/sort children)
        StructureLocked = false;
        
    end
   
    properties (SetAccess = protected)
        
        %What thing(s) are collectable ?
        CollectionSpec = struct('CreateFcn', {}, ... % Function used to create a new instance
            'Description', {}, ... % descriptive text
            'Class', {}, ... % class of new instance
            'MaxNumber', {}, ... % maximum number of instances we're prepared to allow
            'Hidden', {}, ... % if hidden, do not advertise to user the ability to parent this item
            'UniqueNames', {}, ... % enforce uniqueness of names to all instances
            'Group', {}, ... % Name of group, if required
            'DrawChildrenInTree', {}); % Control rendering in tree view
        
        %Classname of Collector nodes
        CollectorClass = 'mvc.model.Collector';
        
        %Are children of this item drawn in tree ?
        DrawChildrenInTree = true;
        
    end
    
    properties (SetAccess = protected, Transient)
        
        %Listeners for the Observable properties
        Listeners;
        
    end
    
    properties (Dependent)
        
        Root; % TODO - MAY need to be hidden else overloaded functions like "isequal" fall over on recursion
        DepthInTree;
        HasChildren;
        HasParent;
        IconFile;
        StructureLockedByParent;
        NumChildren;
        NumDescendants;
        FullName;
        FullNameRev;
        HasGrandChildren;
        GrandChildren;
        
    end
    
    events ( NotifyAccess = private )
        ModelChanging;                      
        ModelChanged;                      
    end

    methods % get/set

        function val = get.HasGrandChildren(obj)
        
            %Easy
            val = obj.HasChildren && any(arrayfun(@(x)x.HasChildren, obj.Children));
            
        end
        
        function val = get.GrandChildren(obj)
        
            %Anything ?
            if obj.HasChildren
                
                %Yes
                val = vertcat(obj.Children.Children);
    
            else
                
                %No
                val = [];
                
            end
            
        end
        
        function val = get.FullName(obj)
            %
            % Returns the name of this object, fully-qualified by
            %  prefixing name(s) of parent(s) all the way up to root node.
        
            %Start here
            val{1} = obj.Name;
            par = obj.Parent;
            
            %Work upwards
            while ~isempty(par)
                
                %Add to list
                val{end+1} = par.Name;
                
                %Next
                par = par.Parent;
                
            end
            
            %Flip the list
            val = flip(val);
            
            %And expand
            val = strjoin(val, filesep);
            
        end

        function val = get.FullNameRev(obj)
            %
            % Returns the name of this object, fully-qualified by
            %  suffixing name(s) of parent(s) all the way up to root node.
        
            %Start here
            val{1} = obj.Name;
            par = obj.Parent;
            
            %Work upwards
            while ~isempty(par)
                
                %Add to list
                val{end+1} = par.Name;
                
                %Next
                par = par.Parent;
                
            end
            
            %And expand
            val = strjoin(val, filesep);
            
        end
        
        function val = get.NumChildren(obj)
            
            %Easy
            val = numel(obj.Children);
            
        end
        
        function val = get.NumDescendants(obj)
        
            %Start with Children
            val = obj.NumChildren;
            
            %If anything
            if val > 0
            
                %Recurse
                val = val + sum(arrayfun(@(x)x.NumDescendants, obj.Children));
                
            end
            
        end
        
        function val = get.StructureLockedByParent(obj)
            
            %Initially false
            val = false;
            
            %Check upstream
            par = obj.Parent;
            while ~val && ~isempty(par)
                val = par.StructureLocked;
                par = par.Parent;
            end
            
        end
        
        function val = get.Root(obj)
            
            %This doesn't work, but it should (geckable ?)
            %             %No parent ?
            %             if isempty(obj.Parent)
            %
            %                 %This is the root
            %                 val = obj;
            %
            %             else
            %
            %                 %Pass it on
            %                 val = obj.Parent.Root;
            %
            %             end
            %
            %Try this instead
            val = obj;
            while ~isempty(val.Parent)
                val = val.Parent;
            end
                        
        end
        
        function n = get.DepthInTree(obj)
            
            %How far down the tree hierarchy is this object ?
            if isempty(obj.Parent)
                
                %We are at the top
                n = 0;
                
            else
                
                %Pass it on
                n = 1 + get(obj.Parent, 'DepthInTree');
                
            end
            
        end
        
        function b = get.HasChildren(obj)
            
            %Simple
            b = ~isempty(obj.Children);
            
        end
        
        function b = get.HasParent(obj)
            
            %Simple
            b = ~isempty(obj.Parent);
            
        end
        
        function val = get.Children(obj)
            
            %Start here
            val = obj.Children;
            
            %Anything ?
            if isempty(val)
                
                %No
                
            else
                
                %Any invalid content ?  TODO: This should never happen
                b = ~isvalid(val);
                if any(b)
                    
                    %Nuke em
                    val(b) = [];
                    
                end
                
            end
        
        end
        
        function set.Children(obj, val)
            
            %Safety-net (this SHOULD get picked up before we get here)
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED (!)']); %#ok<MCSUP>
%             assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent (!)']); %#ok<MCSUP>

            %Setting children to empty matrix
            if isequal([], val)
                
                %Is always fine
                obj.Children = [];
            
            else
                
                %Might be needed later
                prefix = @()[obj.Name, ': attempt to assign ', class(val), ' as child of class ''', class(obj), ''' '];
                
                %Check this object is not a leaf node
                assert(~obj.IsLeafNode, [prefix(), ' :leaf node cannot extend collection']); %#ok<MCSUP>
                
                %Assert that whatever is being added is itself collectable
                b = arrayfun(@(x)isa(x, 'mvc.mixin.Collectable'), val);
                assert(all(b), [prefix(), ' that is not collectable']);
                
                %Any other conditions on what can be collected ?
                cls = {obj.CollectionSpec.Class}; %#ok<MCSUP>
                if ~isempty(cls)
                    
                    %Verify that what is being set is acceptable
                    bCls = arrayfun(@(x)any(cellfun(@(y)isa(x, y), cls)), val);
                    assert(all(bCls), [prefix(), 'object(s) ', strjoin({val(~bCls).Name}, ', '), ' are not of class(es) ', strjoin(cls, ', ')]);
                    
                    %Do we care about the numbers of each (or any) members of the collection ?
                    nMax = [obj.CollectionSpec.MaxNumber]; %#ok<MCSUP>
                    if any(~isinf(nMax))
                        
                        %How many of each ?
                        n = cellfun(@(x)sum(arrayfun(@(y)isa(y, x), val)), cls);
                        assert(all(n <= nMax), [prefix(), ' that contains too many instances of class(es) ', strjoin(cls(n > nMax), ', ')]);
                        
                    end
                    
                    %Do we care about uniqueness of names within each class ?
                    %  Make sure we consider each item in val once only as we
                    %  iterate through elements of collection spec
                    mask = true(size(val));
                    
                    %For each type of collectable
                    for i = 1:numel(obj.CollectionSpec) %#ok<MCSUP>
                        
                        %All done ?
                        if ~any(mask)
                            break;
                        end
                        
                        %Where are instances of this class (that we have not already considered) ?
                        % We COULD identify based on "isa"
                        % b = mask & arrayfun(@(x)isa(x, obj.CollectionSpec(i).Class), val); %#ok<MCSUP>
                        % but that introduces problems with super-classes being matched, so
                        b = mask & arrayfun(@(x)strcmp(class(x), obj.CollectionSpec(i).Class), val); %#ok<MCSUP>
                        
                        %Do not consider these members of val again
                        mask(b) = false;
                        
                        %Anything to worry about ?
                        if ~any(b)
                            
                            %No
                            continue;
                            
                        end
                        
                        %If the spec limits this to a singleton instance
                        if obj.CollectionSpec(i).MaxNumber == 1 %#ok<MCSUP>
                            
                            %Do nothing
                            continue;
                            
                        end
                        
                        %Ensure names are unique
                        makeUniqueNames(val(b));
                        
                    end
                    
                end
                
                %Check if anything being set already has a parent that is NOT this
                bHasDifferentParent = arrayfun(@(x)x.HasParent && (x.Parent ~= obj), val);
                if any(bHasDifferentParent)
                    
                    %Need to detatch gracefully
                    val(bHasDifferentParent).detach;
                    
                end
                
                %Are all names unique, for all children regardless of class ?
                if obj.AllNamesUnique %#ok<MCSUP>
                    
                    %Ensure names are unique
                    makeUniqueNames(val);
                    
                end
                
                %Make a note
                obj.Children = val;
                
                %Ensure children know who they belong to
                [obj.Children.Parent] = deal(obj);
                
            end
            
            %Check to see if the object is being imported as part of a
            %wider import process
            
            %Collection has changed - send message to rebuild
            rebuild(obj);
            
            %If serializable
            if isa(obj.Root, 'mvc.mixin.Serializable') %#ok<MCSUP>
                
                %Make a note that something has changed
                obj.Root.ChangedSinceLastSave = true; %#ok<MCSUP>
                
            end
            
        end
        
        function val = get.IconFile(obj)
            
            %Assume the worst
            val = [];
            
            %Who are we ?
            cls = class(obj);
            
            %Icons placed in subfolder of classdef
            dn = fullfile(fileparts(which(cls)), 'ico');
            
            %Folder exists ?
            if ~isdir(dn)
                
                %No
                return;
                
            end
            
            %Name of icon file must match class name
            cls = strsplit(cls, '.');
            fn = dir(fullfile(dn, [cls{end}, '.*']));
            
            %Anything ?
            if isempty(fn)
                
                %No
                return;
                
            end
            
            %So the answer is
            val = fullfile(dn, fn(1).name);
            
        end
        
    end
    
    methods % construction / destruction
    
        function obj = Collectable(varargin)
            
            %If properties are managed by Nameable, and user has access to "...Locked" flags,
            % AND if this object is NOT inherrently a leaf node
            if isa(obj, 'mvc.mixin.Nameable') && obj.PropertiesLockable && ~obj.IsLeafNode
                
                %Extend property groups
                obj.addPropertyGroup('General', ...
                    'StructureLocked', 'StructureLocked');

                %Ensure this property is never locked
                obj.setPropertyAttribute('NeverLocked', 'StructureLocked');
                
            end
            
            %Extend context, if applicable
            if isa(obj, 'mvc.mixin.Contextable')
                
                %Add context functions to interact with collection
                addContext(obj, 1, ...
                    'Edit', 'editmenu', ...
                    'Edit>Add', [], ... % better than 'Edit>Add...', 'add', ...
                    'Edit>Remove>This...', 'removeThis', ...
                    'Edit>Remove>Children...', 'remove', ...
                    'Edit>Duplicate>This...', 'duplicateAndMoveThis', ...
                    'Edit>Duplicate>Children...', 'duplicateAndMove', ...
                    'Edit>Move>This...', 'moveThis', ...
                    'Edit>Move>Children...', 'move', ...
                    'Edit>Sort>This...', 'sortThis', ...
                    'Edit>Sort>Children...', 'sort', ...
                    'Edit>Properties>This...', 'propedit', ...
                    'Edit>Properties>Children...', 'propeditc');
                
                %And remove previous access to propedit (superceded by Edit>Properties>This...)
                rmContext(obj, 'Edit>Properties...');
                
                %If debuggable too
                if isa(obj, 'mvc.mixin.Debugable')
                
                    %Add an option to manually trigger a rebuild of the collection in any views listening to it
                    addContext(obj, '|Debug>Rebuild', 'rebuild');
                    
                end
                
            end
            
            %Listen for changes to content
            obj.listeners('add');
            
        end
        
        function delete(obj)
            
            %Anything to worry about ?
            if ~isempty(obj.Children)
                
                %Delete any content
                delete(obj.Children);
                
            end
            
        end
        
    end
    
    methods % listener management
        
        function varargout = listeners(obj, action, fcn, varargin)
        
            %Default action
            if nargin < 2
                action = 'get';
            end
            
            %Do what ?
            switch lower(action)
                
                case 'get'
                    
                    %Send them back
                    varargout{1} = [obj.Listeners];
                    
                case {'set', 'add'} % add is legacy - the term 'set' is far more appropriate
                    
                    %Function not specified ?
                    if nargin < 3
                        fcn = @(~,~)rebuild(obj);
                    end
                    
                    %Do it
                    L = arrayfun(@(x)addlisteners(x, fcn), obj, 'UniformOutput', false);
                    [obj.Listeners] = deal(L{:});
                    
                case {'remove', 'delete', 'clear'}
                    
                    %Do it
                    delete([obj.Listeners]);
                    [obj.Listeners] = deal([]);
                    
                case {'pause', 'disable'}
                    
                    %Do it
                    if ~isempty([obj.Listeners])
                        [obj.Listeners.Enabled] = deal(false);
                    end
                    
                case {'resume', 'enable'}
                    
                    %Do it
                    if ~isempty([obj.Listeners])
                        [obj.Listeners.Enabled] = deal(true);
                    end
                    
                otherwise
                    error('bad action');
            end
                   
        end
        
        function rebuild(obj)
        
            %Don't actually need to do anything - just raise the event
            % that triggers any viewers listening to this model,
            % but ensure event is raised from the top
            if isempty(obj.Parent)
                
                %To allow fine control over execution order, raise a "changing" event
                notify(obj, 'ModelChanging');
                
                %Followed by a "changed" event
                %   - BUT only if the object is not being imported (EDIT :
                %   C.Szczyglowski 24/07/2018) Want to stop the GUI
                %   updating all the time during an import process.
                %if ~obj.ObjBeingImported
                    notify(obj, 'ModelChanged');
                %end
                
            else
                
                %Pass it on
                rebuild(obj.Parent);
                
            end
           
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
            
            %But the copy does not automatically have the same parent
            [cpy.Parent] = deal([]);
            
            %Each copy has its own listeners
            cpy.listeners('set');
            
            %And any children need to recursively call their own custom copy
            bHasChildren = [obj.HasChildren];

            %Allow for vectorisation
            for i = find(bHasChildren)     
                
                %Temporarily lower the 'locked' flag
                bLocked = cpy(i).StructureLocked;
                if bLocked
                    cpy(i).StructureLocked = false;
                end
                
                %Do the copy
                cpy(i).Children = copy(obj(i).Children);
                
                %Restore the locked flag
                cpy(i).StructureLocked = bLocked;
                
            end
            
        end
        
    end
    
    methods (Sealed) % collection management
        
        function editmenu(obj, varargin)
            
            %Need current menu - could pass it in, or just go with this
            hc = gcbo;
            
            %Anything to worry about ?
            if isempty(hc)
                
                %No
                return;
                
            end

            %Assume, initially, that we are interested in THIS object
            sel = obj;
            
            %Which object are we really interested in ?
            if ishghandle(hc) && isa(obj, 'mvc.mixin.Serializable')
                
                %Look for relevant vciew manager
                mgr = viewManager(obj, ancestor(hc, 'Figure'));
                
                %Anything ?
                if ~isempty(mgr) && ~isempty(mgr.Selection)
                    
                    %Get the current selection
                    sel = mgr.Selection;
                    
                    %Rebuild 'Edit' submenus on the fly
                    context(sel, hc, 'Edit');
                    
                end
                
            end
            
            %Used in various places - best to get this only once
            bLockedByParent = sel.StructureLockedByParent;
            
            %In which case
            bLocked = sel.StructureLocked || bLockedByParent;
            
            %What actions are supported ?
            bCanAdd = ~bLocked && ~sel.IsLeafNode;
            bCanRemove = sel.HasParent && ~bLockedByParent;
            bCanDuplicate = ~bLockedByParent && ~isa(sel, 'mvc.model.Collector');
            bCanMove = ~bLockedByParent;
            bCanSort = ~bLocked && numel(sel.Children) > 1;
            
            %Dynamically extend menu with Addables (if found)
            i_addables(hc, bCanAdd, sel.CollectionSpec, sel);
            
            %Need to cater separately for HG and Java objects
            if ishghandle(hc)
                
                %Apply can-flags
                set(findobj(hc, 'Label', 'Add...'), 'Enable', mvc.mixin.UiTools.bool2offon(bCanAdd));
                set(findobj(findobj(hc, 'Label', 'Remove'), 'Label', 'This...'), 'Enable', mvc.mixin.UiTools.bool2offon(bCanRemove));
                set(findobj(findobj(hc, 'Label', 'Duplicate'), 'Label', 'This...'), 'Enable', mvc.mixin.UiTools.bool2offon(bCanDuplicate));
                set(findobj(findobj(hc, 'Label', 'Move'), 'Label', 'This...'), 'Enable', mvc.mixin.UiTools.bool2offon(bCanMove));
                set(findobj(hc, 'Label', 'Sort...'), 'Enable', mvc.mixin.UiTools.bool2offon(bCanSort));
                
                %Ability to manipulate Children determined by this or parent status,
                % and whether this object actually has any children anyway
                set(findobj(hc, 'Label', 'Children...'), 'Enable', mvc.mixin.UiTools.bool2offon(~bLocked && obj.HasChildren));
                
            elseif isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
            
                %Apply can-flags
                i_enable(hc, 'Add...', bCanAdd);
                i_enable(hc, 'Remove>This...', bCanRemove);
                i_enable(hc, 'Duplicate>This...', bCanDuplicate);
                i_enable(hc, 'Move>This...', bCanMove);
                i_enable(hc, 'Sort...', bCanSort);
                
                %Ability to manipulate Children determined by this or parent status,
                % and whether this object actually has any children anyway
                i_enable(hc, '*>Children...', ~bLocked && obj.HasChildren);
                
            else
                warning(['unhandled class ''', class(hc), '''']);
            end
            
            function hc = i_enable(hc, nam, val)
            
                %Tokenise
                tok = strsplit(nam, '>');
                
                %For each token
                for i = 1:numel(tok)
                    
                    %If it's a wild-card
                    if strcmp('*', tok{i})
                        
                        %Apply to all submenus
                        for j = 0:hc.getItemCount-1
                            hitem = hc.getItem(j);
                            if isa(hitem, 'javax.swing.JMenu')
                                i_enable(hitem, strjoin(tok(i+1:end)), val);
                            end
                        end
                        
                    else
                        
                        %Get sub-menu names
                        lab = arrayfun(@(x)get(hc.getItem(x), 'Text'), 0:hc.getItemCount-1, 'UniformOUtput', false);
                        
                        %Look for match
                        idx = find(strcmp(tok{i}, lab)) - 1;
                        
                        %Nothing ?
                        if isempty(idx)
                            
                            %Bail out
                            hc = [];
                            return;
                            
                        end
                        
                        %Get the item
                        hc = hc.getItem(idx);
                    
                    end
                    
                end
                
                %Apply flag
                hc.setEnabled(val);
                
            end
            
            function i_addables(hc, bCanAdd, spec, obj)
                
                %Need to cater separately for HG and Java objects
                if ishghandle(hc)
                    
                    %Look for the 'Add' menu
                    hc = findobj(hc, 'Label', 'Add');
                    
                    %Anything to worry about ?
                    if isempty(hc)
                        
                        %No
                        return;
                        
                    end
                    
                    %Is 'Add' explicitly disabled, or is there simply nothing addable ?
                    if ~bCanAdd || isempty(spec) || all([spec.Hidden])
                        
                        %No
                        hc.Enable = 'off';
                        return;
                        
                    end
                    
                elseif isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
                    
                    %Apply can-flags
                    hc = i_enable(hc, 'Add', bCanAdd && ~isempty(spec) && any(~[spec.Hidden]));
                    
                    %Anything to worry about ?
                    if isempty(hc)
                        
                        %No
                        return;
                        
                    end
                    
                end
                
                %So what can be added ?
                idx = find(~[spec.Hidden]);
                str = {spec(idx).Description};
                b = cellfun(@isempty, str);
                str(b) = arrayfun(@(x)func2str(spec(x).CreateFcn), idx(b), 'UniformOutput', false);
                
                %Create suitable sub-menus - borrowing from Contextable to build the menus
                ctx = struct('Label', str, ...
                    'Callback', arrayfun(@(x){'add', spec(x).CreateFcn}, idx, 'UniformOutput', false), ...
                    'Children', [], ...
                    'Priority', 1);
                
                %So we can just pass it on
                obj.context(hc, ctx);
                
            end
            
        end
        
        function obj = detach(obj)
            %
            % Disconnects obj from its parent gracefully, by first removing obj
            %  from list of children in parent, and then removing obj's parent
            
            %For each supplied object
            for i = 1:numel(obj)
                
                %Get parent
                par = obj(i).Parent;
                
                %Anything ?
                if isempty(par) || isempty(par.Children)
                    continue;
                end
                
                %Find this object in parent's children
                [b,idx] = ismember(obj(i), par.Children);
                
                %Anything ?
                if ~b
                    warning('unexpected connectivity error'); % or maybe error ?
                    continue;
                end
                
                %Break link
                par.Children(idx) = [];
                obj(i).Parent = [];
                
            end
            
        end
        
        function rmCollectionSpec(obj, varargin)
        
            %Placeholder
            idx = [];
            
            %For each input
            for i = 1:numel(varargin)
                
                %What have we got ?
                if isa(varargin{i}, 'function_handle')
                    
                    %Look for match
                    b = cellfun(@(x)isequal(x, varargin{i}), {obj.CollectionSpec.CreateFcn});
                    
                elseif ischar(varargin{i})
                    
                    %Look for match in description OR class
                    b = strcmp(varargin{i}, {obj.CollectionSpec.Description}) | ...
                        strcmp(varargin{i}, {obj.CollectionSpec.Class});
                    
                else
                    error('bad input');
                end
                
                %Anything ?
                if any(b)
                    
                    %Append to list
                    idx = [idx, find(b)];
                    
                end
                
            end
            
            %Do it
            obj.CollectionSpec(idx) = [];
            
        end
        
        function varargout = addCollectionSpec(obj, fcn, desc, cls, nMax, bHidden, bUnique, grp, bDrawInTree)
            
            %Special case of no extra input args
            if nargin == 1
                
                %Caller is asking for the current collection spec (e.g. for test execution)
                varargout{1} = obj.CollectionSpec;
                return;
            
            elseif nargin == 2 && isstruct(fcn)
                
                %Caller is passing in a ready-made collection spec structure
                obj.CollectionSpec(end+1) = fcn;
                
                %Needed later
                grp = fcn.Group;
                
            else
            
                %Extend details in store with create function
                obj.CollectionSpec(end+1).CreateFcn = fcn;
                
                %Descriptive text is optional
                if nargin > 2
                    obj.CollectionSpec(end).Description = desc;
                end
                
                %The class of object returned by the create function can also be specified
                if nargin < 4 || isempty(cls)
                    
                    %But if not, defaults as follows
                    cls = func2str(fcn);
                    
                end
                obj.CollectionSpec(end).Class = cls;
                
                %Max number of instances permitted
                if nargin < 5 || isempty(nMax)
                    nMax = Inf;
                end
                obj.CollectionSpec(end).MaxNumber = nMax;
                
                %Hidden ?
                if nargin < 6 || isempty(bHidden)
                    bHidden = false;
                end
                obj.CollectionSpec(end).Hidden = bHidden;
                
                %Unique names ?
                if nargin < 7 || isempty(bUnique)
                    bUnique = true;
                end
                obj.CollectionSpec(end).UniqueNames = bUnique;
                
                %Objects of this type may be grouped under a single node named as follows
                if nargin < 8
                    grp = [];
                end
                obj.CollectionSpec(end).Group = grp;
                
                %If grouped in a collector, they may optionally be not drawn in tree view
                if nargin < 9
                    bDrawInTree = true;
                end
                obj.CollectionSpec(end).DrawChildrenInTree = bDrawInTree;
                
            end
            
            %If the Group attribute is set
            if ~isempty(grp)
                
                %Need to ensure that this collection will permit Collectors to be added
                if ~ismember(obj.CollectorClass, {obj.CollectionSpec.Class})
                    addCollectionSpec(obj, str2func(obj.CollectorClass), [], [], [], true); % Hidden
                end
                
            end
            
        end
        
        function [b, strWhyNot] = canAdd(obj, item, bExcludeHidden, bExcludeGrouped)
            %
            % Checks whether ITEM can be added as child of OBJ
            
            %How best to describe this object, and the item(s) that might be added
            strThis = obj.Name;
            strItem = class(item);
            
            %Assume the worst
            b = false;
            
            %But hope for the best
            strWhyNot = '';
            
            %Easy checks
            if obj.IsLeafNode
                strWhyNot = [strThis, ' is a leaf node and cannot be extended'];
                return;
            elseif ~isa(item, 'mvc.mixin.Collectable')
                strWhyNot = [strItem, ' is not collectable'];
                return;
            elseif ismember(obj, item)
                strWhyNot = 'an object cannot be made a child of itself';
                return;
            elseif isempty(obj.CollectionSpec)
                %No other conditions on what can be collected, so we're good to go
                b = true;
                return;
            end
            
            %Need to consider the collection spec in detail
            spec = obj.CollectionSpec;
            
            %Exclude items in collection spec flagged as hidden ?
            if nargin < 3 || isempty(bExcludeHidden)
                
                %Check if this object is itself a Collection node
                if isa(obj, 'mvc.model.Collector') && ~isempty(obj.Parent)
                    
                    %Where is this group in parent collection spec
                    idx = find(strcmp(strThis, {obj.Parent.CollectionSpec.Group}));
                    
                    %Only hide in THIS if NOT hidden in parent
                    % (might seem a bit random, but it gives the expected behaviour)
                    bExcludeHidden = ~isempty(idx) && ~obj.Parent.CollectionSpec(idx(1)).Hidden;
                    
                else
                    
                    %Simple - yes we do exclude hiden itens by default
                    bExcludeHidden = true;
                    
                end
                
                
            end
            
            %Exclude items in collection spec flagged as grouped ?
            if nargin < 4 || isempty(bExcludeGrouped)
                bExcludeGrouped = true;
            end
            
            %Optionally excluding hidden items
            if bExcludeHidden
                spec([spec.Hidden]) = [];
            end

            
            %Optionally excluding grouped items
            if bExcludeGrouped
                spec(~cellfun(@isempty, {spec.Group})) = [];
            end
            
            %Anything left ?
            if isempty(spec)
                strWhyNot = 'no valid collection spec';
                return;
            end
            
            %Verify that what is being set is acceptable
            cls = {spec.Class};
            bCls = arrayfun(@(x)any(cellfun(@(y)isa(x, y), cls)), item);
            if ~all(bCls)
                strWhyNot = ['item(s) ', strjoin({item(~bCls).Name}, ', '), ' are not of class(es) ', strjoin(cls, ', ')];
                return;
            end
            
            %Do we care about the numbers of each (or any) members of the collection ?
            nMax = [spec.MaxNumber];
            if any(~isinf(nMax))
                
                %How many of each have we already got ?
                n(1,:) = cellfun(@(x)sum(arrayfun(@(y)isa(y, x), obj.Children)), cls);
                
                %How many of each are we thinking about adding ?
                n(2,:) = cellfun(@(x)sum(arrayfun(@(y)isa(y, x), item)), cls);
                
                %OK ?
                if any(sum(n,1) > nMax)
                    strWhyNot = [strThis, ' would then contain too many instances of class(es) ', strjoin(cls(sum(n,1) > nMax), ', ')];
                    return;
                end
                
            end
            
            %If we get this far
            b = true;
            
        end
        
        function item = add(obj, item, varargin)
            
            %Caller can optionally override the Locked state
            % (required, for example, when importing content from an XML file, and the file
            %  itself contains values of the locked flags setting them true.
            %  Without this override facility, such files become unimportable)
            [bIgnoreLockedFlags, varargin] = checkOption('-ignoreLockedFlags', false, varargin{:});
            %             bIgnoreLockedFlags = false;
            %             if ~isempty(varargin)
            %                 idx = find(cellfun(@(x)ischar(x) && strcmp(x, '-ignoreLockedFlags'), varargin));
            %                 if ~isempty(idx)
            %                     bIgnoreLockedFlags = true;
            %                     varargin(idx) = [];
            %                 end
            %             end
            
            %Safety net
            if ~bIgnoreLockedFlags
                assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
                assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
            end
            
            %Add what ?
            if nargin > 1
                
                %Caller has provided the thing to be added
                if isa(item, 'function_handle')
                    
                    %In the form of a function handle to be evaluated
                    item = item(varargin{:});
                    
                    %Group name, if applicable
                    [grp, gdx] = i_findGroup(item, obj.CollectionSpec);
                
                elseif ischar(item)
                    
                    %Caller providing a name - look it up in descriptions
                    desc = {obj.CollectionSpec.Description};
                    
                    %Replace any empties with constructor name
                    desc(cellfun(@isempty, desc)) = cellfun(@func2str, ...
                        {obj.CollectionSpec(cellfun(@isempty, desc)).CreateFcn}, 'UniformOutput', false);
                    
                    %Look for match, this would be the easy way
                    % b = contains(desc, item);
                    %
                    % but for backward compat with R2015b
                    b = ~cellfun(@(x)isempty(strfind(x, item)), desc);
                    
                    %Must be unique
                    assert(any(b), 'no match');
                    assert(sum(b) == 1, 'ambiguous match');
                    
                    %So the answer is
                    item = obj.CollectionSpec(b).CreateFcn(varargin{:});                    
                    grp = obj.CollectionSpec(b).Group;
                    gdx = find(b);
                    
                else
                    
                    %Any extras ?
                    if ~isempty(varargin)
                    
                        %Apply to item before adding
                        set(item, varargin{:});
                    
                    end
                    
                    %Group name, if applicable
                    [grp, gdx] = i_findGroup(item, obj.CollectionSpec);
                
                end
             
            elseif isempty(obj.CollectionSpec)
                
                %Add another one of ourselves (?)
                item = feval(class(obj));
                
                %Group name, if applicable
                [grp, gdx] = i_findGroup(item, obj.CollectionSpec);
                
            else
                
                %Only offer items detailed in spec that are NOT flagged as hidden
                bdx = find(~[obj.CollectionSpec.Hidden]);
                
                %Nothing ?
                if isempty(bdx)
                    
                    %In a well design interface, this would never happen
                    error('nothing addable by user');
                    
                end
                
                %Get descriptions
                desc = {obj.CollectionSpec(bdx).Description};
                
                %Replace any empties with constructor name
                desc(cellfun(@isempty, desc)) = cellfun(@func2str, ...
                    {obj.CollectionSpec(bdx(cellfun(@isempty, desc))).CreateFcn}, 'UniformOutput', false);
                
                %Ask the user
                [~, idx] = choose(obj, 'Add what...', desc{:});
                
                %Cancelled ?
                if isempty(idx)
                    
                    %Bail out
                    item = [];
                    return;
                    
                end

                %Now do it
                item = feval(obj.CollectionSpec(bdx(idx)).CreateFcn, varargin{:});
                
                %May need to add to a group node
                grp = obj.CollectionSpec(bdx(idx)).Group;
                gdx = find(bdx(idx));
                
            end
            
            %Is this item part of a group ?
            if isempty(grp)
                
                %What a pain
                if bIgnoreLockedFlags
                    bStructureLocked = obj.StructureLocked;
                    obj.StructureLocked = false;
                end
                
                %No - add item directly to collection
                if isempty(obj.Children)
                    
                    %Assign - always as column (??)
                    obj.Children = item(:);
                    
                else
                    
                    %Extend - always as column (??)
                    obj.Children(end+1:end+numel(item),1) = item(:);
                    
                end

                %Restore locked status, if necessary
                if bIgnoreLockedFlags && bStructureLocked
                    obj.StructureLocked = true;
                end
                
            else
                
                %Look for an existing group of this name
                if isempty(obj.Children)
                    
                    %Nothing, obviously
                    par = [];
                    
                else
                    
                    %Search existing
                    par = find(obj.Children, 'isa', item(1).CollectorClass, 'Name', grp);
                    
                end
                
                %What a pain
                if bIgnoreLockedFlags
                    ignoreArg = {'-ignoreLockedFlags'};
                else
                    ignoreArg = {};
                end  
                
                %Nothing yet ?
                if isempty(par)
                    
                    %Create the group node
                    par = obj.add(str2func(item(1).CollectorClass), 'Name', grp);
                    
                    %Configure the group node to accept items of the relevant class
                    % Is it as simple as this ?
                    gspec = obj.CollectionSpec(gdx);
                    gspec.Group = [];
                    par.addCollectionSpec(gspec);
                    
                    %Does this node draw its children in tree view ?
                    par.DrawChildrenInTree = gspec.DrawChildrenInTree;
                    
                end
                
                %Add the item to the group
                par.add(item, ignoreArg{:});
                
            end
                                
            function [grp,gdx] = i_findGroup(item, CollectionSpec)
               
                %If there is no spec
                if isempty(CollectionSpec)
                    
                    %Then there is no group
                    grp = [];
                    gdx = [];
                    
                else
                    
                    %Compare class of item with classes associated with collection
                    gdx = find(strcmp(class(item), {CollectionSpec.Class}));
                    if isempty(gdx)
                        
                        %Open up the search, using 'isa' instead
                        gdx = find(cellfun(@(x)isa(item, x), {CollectionSpec.Class}));
    
                    end
                    
                    %Nothing (? is this ok ?)
                    if isempty(gdx)
                        
                        %Then there is no group
                        grp = [];
                        
                    else
                        
                        %So candidate group name(s)
                        grps = {CollectionSpec(gdx).Group};
                        
                        %Eliminate any empties
                        gdx(cellfun(@isempty, grps)) = [];
                        grps(cellfun(@isempty, grps)) = [];
                        
                        %Anything left ?
                        if isempty(grps)
                            
                            %No - there is no group
                            grp = [];
                            
                        else
                            
                            %Eliminate any duplicates
                            [grps, udx] = unique(grps);
                            gdx = gdx(udx);
                            
                            %Safety net - may have to deal with this more gracefully
                            assert(numel(grps) == 1, 'ambiguous match');
                            
                            %So the group name is
                            grp = grps{1};
                            gdx = gdx(1);
                        
                        end
                        
                    end
                    
                end
                
            end                

        end
        
        function item = remove(obj, item, bForce, varargin)
            
            %Safety net
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
            assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
                        
            %Remove what ?
            if nargin > 1
                
                %Caller has provided item(s) to be removed,
                % locate them in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %User selects item(s) to be removed
                [sel, bOK] = select(obj);
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    item = [];
                    return;
                    
                end
                
            end
            
            %Caller can force, otherwise seek confirmation
            if nargin > 2 && islogical(bForce) && bForce
                
                %Just crack on
                
            elseif ~confirm(obj.Children(sel), 'Irreversibly remove object(s) and ALL associated children ?')
                
                %Bail
                item = [];
                return;
                
            end

            %If caller wants them back
            if nargout > 0
                
                %Pull 'em out
                item = obj.Children(sel);
                
            end
            
            %Remove them from collection
            obj.Children(sel) = [];            
                
        end
        
        function varargout = removeThis(obj, varargin)
            
            %Safety net
            assert(obj.HasParent, [obj.Name, ': no parent from whom to remove']);
            
            %Pass it on to parent
            [varargout{1:nargout}] = remove(obj.Parent, obj, varargin{:});
            
        end
        
        function [sel, bOK] = sort(obj, item, varargin)
            
            %Edit C.Szczyglowski 31/01/2019
            %   - If an object array is passed in then assume 'sort' is
            %   being called from unique. In which case just return the
            %   same data.
            if numel(obj) > 1
                sel = obj;
                bOK = 1 : numel(obj);
                return
            end
            
            %Safety net
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
            assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
                        
            %Sort what ?
            if nargin > 1
                                
                %Caller has provided item(s) to be sorted,
                % locate them in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %No initial selection
                sel = [];
                
            end
            
            %Ask the user to sort children
            [sel, bOK] = sortdlg(obj, arrayfun(@dlgtitle, obj.Children, 'UniformOutput', false), ...
                [], 'InitialValue', sel);
            
            %Cancelled ?
            if isempty(sel) || ~bOK
                
                %Bail out
                return;
                
            end
            
            %Apply sort
            obj.Children = obj.Children(sel);
            
        end
        
        function varargout = sortThis(obj, varargin)
            
            %Safety net
            assert(obj.HasParent, [obj.Name, ': no parent from whom to sort']);

            %Pass it on to parent
            [varargout{1:nargout}] = sort(obj.Parent, obj, varargin{:});
            
        end
        
        function item = move(obj, item, dest, varargin)
            
            %Safety net
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
            assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
                        
            %Move what ?
            if nargin > 1
                
                %Caller has provided item(s) to be moved,
                % locate them in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %User selects item(s) to be moved
                [sel, bOK] = select(obj);
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    item = [];
                    return;
                    
                end
                
            end
            
            %Pull 'em out
            item = obj.Children(sel);
            
            %Move to where ?  Could be anywhere within network
            if nargin > 2
                
                %Caller has provided destination
                assert(isa(dest, 'mvc.mixin.Collectable'), 'invalid destination');
                
            else
                
                %To where could the item(s) be moved ?
                lst = obj.Root.flatlist(2);
                b = arrayfun(@(x)x.canAdd(item), lst);
                assert(sum(b) > 1, 'no available destinations for selected item(s)');

                %Down-select
                lst = lst(b);

                %Any choice ?
                if numel(lst) == 1
                    
                    %No
                    dest = lst;
                    
                else
                    
                    %Ask the user, initialising selection with current parent
                    [b, val] = ismember(item.Parent, lst);
                    val = val(b);
                    
                    %User selects destination
                    dest = obj.selectFrom(lst, 'Move to...', ...
                        'InitialValue', val, ...
                        'SelectionMode', 'single');
                    
                    %Cancelled ?
                    if isempty(dest)
                        
                        %Bail out
                        item = [];
                        return;
                        
                    end
                
                end
                
            end
            
            %Move them
            dest.add(item);
                            
        end
        
        function varargout = moveThis(obj, varargin)
            
            %Safety net
            assert(obj.HasParent, [obj.Name, ': no parent from whom to move']);

            %Pass it on to parent
            [varargout{1:nargout}] = move(obj.Parent, obj, varargin{:});
            
        end
        
        function item = duplicate(obj, item, varargin)
            
            %Safety net
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
            assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
            
            %Duplicate what ?
            if nargin > 1
                
                %Caller has provided item(s) to be removed,
                % locate them in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %User selects item(s) to be duplicated
                [sel, bOK] = select(obj);
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    item = [];
                    return;
                    
                end
                
            end

            %If caller wants them back
            if nargout > 0
                
                %Pull 'em out
                item = obj.Children(sel);
                
            end
            
            %Duplicate them within collection
            obj.Children(end+1:end+numel(sel),:) = copy(obj.Children(sel));
                
        end
        
        function varargout = duplicateThis(obj, varargin)
            
            %Safety net
            assert(obj.HasParent, [obj.Name, ': no parent from whom to duplicate']);
            
            %Pass it on to parent
            [varargout{1:nargout}] = duplicate(obj.Parent, obj, varargin{:});
            
        end
        
        function item = duplicateAndMove(obj, item, dest, varargin)
            
            %Safety net
            assert(~obj.StructureLocked, [obj.Name, ': structure LOCKED']);
            assert(~obj.StructureLockedByParent, [obj.Name, ': structure LOCKED by parent']);
            
            %Duplicate what ?
            if nargin > 1
                
                %Caller has provided item(s) to be removed,
                % locate them in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %User selects item(s) to be duplicated
                [sel, bOK] = select(obj);
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    item = [];
                    return;
                    
                end
                
                %Pull 'em out
                item = obj.Children(sel);
                
            end
            
            %Move to where ?  Could be anywhere within network
            if nargin > 2
                
                %Caller has provided destination
                assert(isa(dest, 'mvc.mixin.Collectable'), 'invalid destination');
                
            else
                
                %To where could the item(s) be moved ?
                lst = obj.Root.flatlist(2);
                b = arrayfun(@(x)x.canAdd(item), lst);
                
                %If there are no feasible destinations, but we could
                % assign the duplicate to a new instance of application,
                % AND there are no return arguments
                if ~any(b) && isa(obj, 'mvc.model.Application')
                    
                    %Make a new instance of application, as destination
                    dest = obj.new;
                    
                    %No return arguments ?
                    if nargout == 0
                        
                        %Open it
                        view(dest);
                        
                    end
                    
                else
                    
                    %Valid to proceed ?
                    assert(any(b), 'no feasible destinations for selected item(s)');
                    
                    %Down-select
                    lst = lst(b);
                    
                    %Any choice ?
                    if numel(lst) == 1
                        
                        %No
                        dest = lst;
                        
                    else
                        
                        %Ask the user, initialising selection with current parent
                        [b, val] = ismember([item.Parent], lst);
                        val = unique(val(b));
                        
                        %User selects destination
                        dest = obj.selectFrom(lst, 'Move duplicate(s) to...', ...
                            'InitialValue', val(1), ...
                            'SelectionMode', 'single');
                        
                        %Cancelled ?
                        if isempty(dest)
                            
                            %Bail out
                            item = [];
                            return;
                            
                        end
                        
                    end
                    
                end
            
            end
            
            %Duplicate and move them
            dest.add(copy(obj.Children(sel)));
                
        end
        
        function varargout = duplicateAndMoveThis(obj, varargin)
            
            %If being called on an entire Application, with no return arguments
            if isa(obj, 'mvc.model.Application')
                
                %Create duplicate
                cpy = copy(obj, varargin{:});
                
                %TODO: should persistent file name be included in copy ?
                
                %No return arguments ?
                if nargout == 0
        
                    %Open it
                    view(cpy);
                
                else
                    
                    %Send it back
                    varargout{1} = cpy;
                    
                end
                
            else
                
                %Safety net
                assert(obj.HasParent, [obj.Name, ': no parent from whom to duplicate and move']);
                
                %Pass it on to parent
                [varargout{1:nargout}] = duplicateAndMove(obj.Parent, obj, varargin{:});
            
            end
            
        end
        
        function item = propeditc(obj, item, varargin)
            
            %Properties of what ?
            if nargin > 1
                
                %Caller has provided item - locate in collection
                [b, sel] = ismember(item, obj.Children);
                
                %They MUST all be present
                assert(all(b), 'invalid input');
                
            else
                
                %User selects item(s) to be duplicated
                [sel, bOK] = select(obj, [], 'SelectionMode', 'single');
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    item = [];
                    return;
                    
                end
                
                %Pull 'em out
                item = obj.Children(sel);
                
            end

            %Now pass it on
            propedit(item, varargin{:});
                            
        end
                
        function [sel, bOK] = selectFrom(obj, lst, varargin)
            
            %Caller may be telling us
            if numel(varargin) == 1 && (isnumeric(varargin{1}) || islogical(varargin{1}))
                
                %Index of selection being passed in directly
                sel = varargin{1};
                bOK = true;
                
            elseif numel(varargin) == 1 && isa(varargin{1}, 'mvc.mixin.Collectable')
                
                %The actual selection being passed in directly
                [b, sel] = ismember(varargin{1}, lst);
                assert(all(b), 'invalid selection');
                bOK = true;

            else
                
                %Good names to offer
                str = arrayfun(@dlgtitle, lst, 'UniformOutput', false);
                
                %Indentation i.a.w. hierarchy is good
                n = arrayfun(@(x)x.DepthInTree, lst);
                str = arrayfun(@(x)[repmat(' ', 1, n(x)), str{x}], 1:numel(str), 'UniformOutput', false);
                
                %Make a selection
                [sel, bOK] = listdlg(obj, str, varargin{:});
            
            end
            
            %Return selected object(s) (because index into flattened list may have little meaning)
            sel = lst(sel);
            
        end
        
        function [sel, bOK] = selectFromAll(obj, varargin)
            
            %Start from the top
            top = obj.Root;
            
            %Get all in network
            all = flatlist(top, 2);
            
            %Now pass it on
            [sel, bOK] = selectFrom(obj, all, varargin{:});
            
        end
        
        function [sel, bOK] = select(obj, varargin)
            
            %There must be something to select
            assert(~isempty(obj.Children), 'nothing to select');
            
            %Present choice using dlgtitle, or FullNameRev ?  The latter is better, I think
            %str = arrayfun(@dlgtitle, obj.Children, 'UniformOutput', false);
            str = {obj.Children.FullNameRev};
            
            %Ask the user to select
            [sel, bOK] = listdlg(obj, str, varargin{:});
            
        end
        
        function lst = flatlist(obj, mode)
        
            %Different ways of returning the flattened list
            if nargin < 2 || mode == 1
                
                %This looks neat, in terms of the code, but it returns the list "inside out",
                % with ordering of elements being sorted by depth on tree
                
                %Any children to worry about ?
                if isempty(vertcat(obj.Children))
                    
                    %No
                    lst = obj;
                    
                else
                    
                    %Return this object and any children in flat list
                    lst = vertcat(obj(:), flatlist(vertcat(obj.Children), 1));
                    
                end
                
                %Return a column
                lst = lst(:);
            
            elseif mode == 3
            
                %This looks clunky, but gives a more intuitive result
                % where the order in the flat list is the same as you'd see
                % if the application structure was drawn in a tree
                lst = {};
                for i = 1:numel(obj)
                    
                    %Start with this object
                    lst{end+1} = obj(i);
                    
                    %Safety net - should not be necessary
                    b = arrayfun(@isvalid, obj(i).Children);
                    
                    %Followed immediately by its children
                    lst(end+1:end+numel(obj(i).Children(b))) = arrayfun(@(x)flatlist(x, 2), obj(i).Children(b), 'UniformOutput', false);
                    
                end
                
                %And expand, and return as column
                lst = horzcat(lst{:});
                
            elseif mode == 2
            
                %Same as mode 3, but with pre-allocation (profile shows mode 3 is slow)
                lst = arrayfun(@i_flatten, obj, 'UniformOutput', false);
                
                %And expand
                lst = horzcat(lst{:});
                
            else
                error('bad mode');
            end
            
            function lst = i_flatten(obj)
                
                %Safety net - should not be necessary
                b = arrayfun(@isvalid, obj.Children);
                
                %Preallocate
                lst = cell(1, 1 + sum(b));
                
                %Start with this object
                lst{1} = obj;
                
                %Followed by its children
                lst(2:1+sum(b)) = arrayfun(@(x)flatlist(x, 3), obj.Children(b), 'UniformOutput', false);

                %And expand
                lst = horzcat(lst{:});
                
            end
            
        end
                
        function [hNode, nItems] = uitreenode(obj, reloadfcn, fld, hPar, nItems)
                        
            %Field to display in tree ?
            if nargin < 3 || isempty(fld)
                fld = 'Name';
            end
            
            %Useful to count the total number of items in tree
            if nargin < 5
                nItems = 0;
            end
            
            %Good name for this node ?
            nam = obj.(fld);
            
            %Append with indication of number of children
            if ~obj.DrawChildrenInTree && ~isempty(obj.Children)
                nam = [nam, ' (', num2str(numel(obj.Children)), ')'];
            end
            
            %Create a tree node
            hNode = uitreenode('v0', 0, nam, obj.IconFile, isempty(obj.Children));
            
            %Increment counter
            nItems = nItems + 1;
            
            %Add object as user data to its own node
            set(hNode, 'UserData', obj);
            
            %Listen for changes
            addlistener(obj, 'Name', 'PostSet', @i_rename);

            %Add to parent, if supplied
            if nargin > 3
                hPar.add(hNode);
            end
            
            %Draw children in tree ?
            if obj.DrawChildrenInTree
                
                %For each child
                for i = 1:numel(obj.Children)
                    
                    %Pass it on
                    [~, nItems] = uitreenode(obj.Children(i), reloadfcn, fld, hNode, nItems);
                    
                end
            
            end
            
            function i_rename(~, evtdata)
                
                %Rename the node
                hNode.Name = evtdata.AffectedObject.(fld);
                
                %Refresh the tree
                reloadfcn(hNode);
                
            end
            
        end
        
    end
    
end
