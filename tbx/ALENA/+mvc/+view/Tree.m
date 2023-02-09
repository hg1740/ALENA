classdef Tree < mvc.view.Container

    properties (AbortSet)
        
        %The model to which we are attached
        Model;
        
        %The field to display in tree
        Field = 'Name';
        
        %Track currently selected item in model
        Selection;
        
        %Function to call on selection change
        SelectionChangeFcn;
        
    end
    
    properties (Access = protected)
        
        %The tree control in which we display the model
        hTree;
        
        %Needed to prevent race condition
        bFindNode = true;
        
        %Placeholder for listeners
        Listeners;
        
        %Support for Drag-and-Drop
        EnableDragAndDrop = true;
        EnableFileDrop = true;
        
    end
    
    properties (Access = private)
        
        %Helper for file drop
        FileDropControl;
        
    end
    
    properties (Hidden, Dependent)
        
        TreeVisible;
        
    end
    
    methods % get/set
        
        function set.TreeVisible(obj, val)
            
            %Allow for logical, numeric or char input
            if islogical(val)
                
                %OK
                
            elseif isnumeric(val)
                
                %Cast
                val = logical(val);
                
            elseif ischar(val)
                
                %Convert
                if strcmpi(val, 'on')
                    val = true;
                elseif strcmpi(val, 'off')
                    val = false;
                else
                    error('bad value');
                end
            else
                error('bad input');
            end
            
            %Pass it on to underlying tree control
            obj.hTree.Visible = val;
            
        end
        
        function val = get.TreeVisible(obj)
            
            %Pass it on to underlying tree control
            val = obj.hTree.Visible;
            
        end
        
        function set.Selection(obj, val)
            
            %Allow selection to be set with an index into the list
            if isnumeric(val)
                
                %Get from (flattened) list
                L = flatten(obj.Model, 2); %#ok<MCSUP>
                val = L(val);
                
            end            
            
            %Anything to do ?
            if obj.Selection == val
                
                %No
                return;
                
            end
            
            %Make a note
            obj.Selection = val;
        
            %Ensure consistency with tree, if necessary
            if obj.bFindNode %#ok<MCSUP>
                findnode(obj, val);
            end
            
        end
        
        function set.Model(obj, val)
            
            %Allow model to be nothing at al
            if isequal(val, [])
                
            else
                
                %It MUST be Collectable
                assert(isa(val, 'mvc.mixin.Collectable'), 'model MUST be Collectable');
        
            end
            
            %Make a note
            obj.Model = val;
            
            %Update content
            update(obj);
            
        end
        
    end
    
    methods % construction
        
        function obj = Tree(model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Store the model
            obj.Model = model;
            
            %Create tree from model
            update(obj);
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
        end
        
        function delete(obj)
            
            %Make sure the underlying tree control goes with us
            try
                delete(obj.hTree);
            end
            
        end
        
    end
    
    methods (Access = private) % event handlers
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onFileDropped(obj, hc, evtdata)
        
            %Are we bothered ?
            if isempty(obj.EnableFileDrop)
                
                %No
                return;
                
            end
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                 
                %Generic file drop could be treated either as File->Open,
                % or as File->Import, not immediately obvious which should
                % take priority.  But Open seems like a better bet (? maybe ?)
                if isa(obj.Model, 'mvc.mixin.Serializable') && obj.Model.isOpenable(evtdata.Data{:})
                    
                    %Pass it on
                    obj.Model.open(evtdata.Data{:});
                    
                elseif isa(obj.Model, 'mvc.mixin.Importable') && obj.Model.isImportable(evtdata.Data{:})
                    
                    %Pass it on
                    obj.Model.import(evtdata.Data{:});
                
                else
                    warning('onFileDrop: don''t know what to do');
                end
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
               
        end
        
        function onNodeDropped(obj, hc, evtdata)
        
            %Are we bothered ?
            if ~obj.EnableDragAndDrop
                
                %No
                return;
                
            end
            
            %In case it takes a while
            clu = obj.pointer(hc); %#ok<NASGU>
            
            %Careful
            try
                
                %Pull out source and target nodes
                src = get(handle(evtdata.getSourceNode), 'UserData');
                tgt = get(handle(evtdata.getTargetNode), 'UserData');
                
                %Try to move
                src.moveThis(tgt);
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(err.message, obj.Title, 'modal'));
                
            end
            
        end
        
        function onExpanded(obj, hc, evtdata, flag)
            
            %Get underlying object
            sel = get(get(evtdata, 'CurrentNode'), 'UserData');
            
            %Make a note
            if isvalid(sel)
                sel.IsExpanded = flag;
            end
            
        end
        
        function onSelectionChange(obj, hc, evtdata)
            
            %Might take some time
            clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                        
            %Needed to prevent race condition
            obj.findnode(false);
            clu = onCleanup(@()obj.findnode(true));
            
            %Need to cater for possibility that the change is to NO selection
            if isempty(hc.getSelectedNodes)
                
                %Assign nothing
                obj.Selection = [];
                
            else
             
                %Make a note
                obj.Selection = get(get(evtdata, 'CurrentNode'), 'UserData');
    
            end
            
            % Does selected item support context menus ?
            if isa(obj.Selection, 'mvc.mixin.Contextable')
                
                %Yes - pass it on
                context(obj.Selection, obj.hTree);
                %context(obj.Selection, obj.hTree.Tree);
                
            end
                
            %Anything else to do ?
            if ~isempty(obj.SelectionChangeFcn)
                
                %Pass it on
                obj.SelectionChangeFcn(obj.Selection);
                
            end
            
        end
        
        function update(obj)
            
            %Can we make a note of what is currently expanded ?
            %TODO: not quite working yet, so revert to ALWAYS expand all
            if true || isempty(obj.hTree)
                
                %No
                bExpanded = [];
                
            else
                
                %Yes
                bExpanded = arrayfun(@(x)obj.hTree.Tree.isExpanded(x-1), 1:obj.hTree.Tree.getRowCount);
                
            end
            
            %Discard any existing tree
            delete(obj.hTree);
            
            %And filedrop helper, if it exists
            delete(obj.FileDropControl);
            
            %Allow for an empty model
            if isempty(obj.Model)
                
                %Nothing to show
                hRoot = [];
                nItems = 0;
                
            else
                
                %Convert model to corresponding treenode(s)
                [hRoot, nItems] = uitreenode(obj.Model, @i_reloadNode, obj.Field);
            
            end
            
            % Tree create does not honour the 'Parent' property,
            % so have to do it the hard way
            hf = ancestor(obj.UIContainer, 'Figure');
            cf = get(0, 'CurrentFigure');
            clu = onCleanup(@()set(0, 'CurrentFigure', cf(ishghandle(cf))));
            set(0, 'CurrentFigure', hf);
            
            % Create tree control
            obj.hTree = uitree('v0', 'Parent', hf, 'Root', hRoot);
            
            %Adjust position to fill container, using suitable helper
            set(obj.hTree, 'Units', get(obj.UIContainer, 'Units'));
            i_sizeChangeFcn;

            %Synchronise tree control position to that of container
            set(obj.UIContainer, 'SizeChangedFcn', @i_sizeChangeFcn);
            
            %Callback on node selection
            set(obj.hTree, 'NodeSelectedCallback', @obj.onSelectionChange);
            
            %Callbacks to manage Expanded property
            set(obj.hTree, 'NodeExpandedCallback', {@obj.onExpanded, true});
            set(obj.hTree, 'NodeCollapsedCallback', {@obj.onExpanded, false});
            
            %Support edit by drag and drop ?
            if obj.EnableDragAndDrop
                
                %Set callback accordingly
                set(obj.hTree, 'NodeDroppedCallback', @obj.onNodeDropped);
                
            end
            
            %Support handling of drop-files ?
            if ~isempty(obj.EnableFileDrop)
                
                %Careful - this could fail if Java path not set correctly or code not found
                try
                    
                    %Is FileDrop support enabled ?
                    if ~mvc.view.internal.dndcontrol.isInitialized
                        
                        %Initialise
                        mvc.view.internal.dndcontrol.initJava;
                        
                    end
                    
                    %Get help
                    obj.FileDropControl = mvc.view.internal.dndcontrol(obj.hTree.tree, @obj.onFileDropped, @obj.onStringDropped);
        
                catch err
                    
                    %Tell us about it
                    warning(['FAILED to enable file drop, with error: ', err.message]);
                    
                end
                
            end
            
            %Reinstate expanded state
            i_expand(obj.hTree, obj.Model, hRoot);
            
            function i_expand(hTree, mdl, hNode)
                
                %Anything to do ?
                if ~mdl.IsExpanded
                    
                    %No
                    return;
                    
                end
                
                %Expand this object
                hTree.expand(hNode);
                
                %Need to worry about children ?
                if mdl.DrawChildrenInTree
                    
                    %For each in turn
                    for i = 1:numel(mdl.Children)
                        
                        %Get corresponding node
                        if i == 1
                            hChild = hNode.getFirstChild;
                        else
                            hChild = hChild.getNextSibling;
                        end
                        
                        %Pass it on
                        i_expand(hTree, mdl.Children(i), hChild);
                        
                    end
                    
                end
                
            end
            
            function i_sizeChangeFcn(~,~)
                
                %Guard against unexpected change of units
                un = 'pixels';
                set(obj.hTree, 'Units', un);
                
                % Tree control 'Parent' is always the figure, not the container,
                % so we have to work out the absolute position (re figure)
                % by working our way up the parent hierarchy
                par = obj.Parent;
                parpos = i_position(par, un);
                while ~isa(par, 'matlab.ui.Figure')
                    par = par.Parent;
                    parpos(end+1,:) = i_position(par, un);
                end
                
                %Add them up, excluding the last 
                parpos = sum(parpos(1:end-1,:),1);
                
                %Get position of the container
                pos = i_position(obj.UIContainer, un);
                
                %Convert from relative to absolute
                pos(1:2) = pos(1:2) + parpos(1:2);
                
                %Apply to tree
                set(obj.hTree, 'Position', pos);
                
                function pos = i_position(hc, un)
                
                    %Ensure consistent units
                    un_orig = hc.Units;
                    hc.Units = un;
                    
                    %And get position
                    pos = hc.Position;
                    
                    %Restore if necessary
                    hc.Units = un_orig;
                    
                end
                
            end
            
            %Helper function returns tree to caller
            function i_reloadNode(hNode)
                
                %Ensure valid (this wouldn't be necessary if listener lifecycles were all managed correctly)
                if ~isvalid(obj)
                    return;
                end
                
                %Pass it on to tree
                obj.hTree.reloadNode(hNode);
                
                %Reinstate expand state
                %i_expand(obj.hTree, obj);
                
            end
            
        end
        
        function hNode = findnode(obj, varargin)
        
            %If called with a single boolean
            if numel(varargin) == 1 && islogical(varargin{1})
                
                %We are setting the 'findnode' mode
                obj.bFindNode = varargin{1};
                return;
                
            end

            %Check currently selected node, before bothering to traverse the entire tree
            if i_checknode(getSelectedNodes(obj.hTree), varargin{:})
                
                %We're there already
                hNode = getSelectedNodes(obj.hTree);
                hNode = hNode(1);
                
            else
                
                %Start from tree root
                hRoot = getRoot(obj.hTree);
                
                %Check this node
                hNode = i_findnode(hRoot, varargin{:});
                
                %Anything found, but caller did not want it back ?
                if nargout == 0 && ~isempty(hNode)
                    
                    %Make this the selected node in tree
                    setSelectedNode(obj.hTree, hNode);
                    
                end
                
            end
            
            function hMatch = i_findnode(hNode, varargin)
                
                %Loop until we find a match, or run out of steam
                while true
                    
                    %Check this node
                    if i_checknode(hNode, varargin{:})
                        
                        %We're done
                        hMatch = hNode;
                        return;
                        
                    end
                    
                    %Next
                    hNode = hNode.getNextNode;                
                    
                end
                 
                %If all else fails
                hNode = [];
                
            end
               
            function bMatch = i_checknode(hNode, varargin)
                
                %Assume the worst
                bMatch = false;
                
                %If nothing supplied
                if isempty(hNode)
                    
                    %Nothing else to do
                    return;
                    
                end
                
                %If single collectable argument supplied
                if numel(varargin) == 1 && isa(varargin{1}, 'mvc.mixin.Collectable')
                    
                    %Take a short-cut
                    if get(hNode(1), 'UserData') == varargin{1}
                        
                        %It's a match
                        bMatch = true;
                        
                    end
                    
                else
                    
                    %Get whatever details from current node
                    val = get(hNode(1), varargin(1:2:end));
                    
                    %Match ?
                    if all(cellfun(@(x,y)isequaln(x,y), val, varargin(2:2:end)))
                        
                        %It's a match
                        bMatch = true;
                        
                    end
    
                end
            
            end
            
        end
        
    end
    
end
