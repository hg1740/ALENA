classdef (ConstructOnLoad) Serializable < matlab.mixin.SetGet
    %
    %Serializable implements a number of features useful to an
    % application that needs to be capable of serializing a "session" to file:
    %
    %  Management of the name and path of the file to which the application is serialized
    %  Presentation of an "MRU" list (most-recently used files) to allow quick load from file
    %  Maintenance of a "ChangedSinceLastSaved" flag, so the application knows
    %   whether to prompt the user before closing a session down
    %  Management of preferences (as per setpref/getpref)
    
    properties (AbortSet, SetObservable, Transient) % Observable risks endless recursion when interacting with Auditable
        
        %Need to know whether content has changed since last saved to file
        ChangedSinceLastSave = false;
        
    end
    
    properties (AbortSet, SetObservable, SetAccess = protected, Transient)
        
        %File and path to/from which content is serialized
        SerializeFile;
        SerializePath;
        
    end
    
    properties (Access = protected, Transient)
        
        %To fine tune the serialization process
        SerializeSpec = { ...
            '*.mat',               'MATLAB data files (*.mat)', [], []; ...
            }; % '*.*',                 'All files (*.*)',           [], []};
        SerializeIndex = [];
        
        %A name is required to help distinguish MRU associated with this application, from any others
        SerializeName;
        
        %For managing MRU
        SerializeRecentMax = 5;
        
        %Whether to include the ViewArrangement in serialised file
        SerializeViewArrangement = false;
        
    end
    
    properties (Hidden, Transient)
        
        %List of properties to be copied from one object to the next, when 'New' is invoked
        CopyOnNew = {};
        
    end
    
    properties (Dependent)
        
        %For convenience
        SerializeFullFile;
        SerializeFullFileMRU;
        
        %For managing MRU
        SerializeRecent;
    
        %For position new windows when necessary
        LastPosition;
        
    end
    
    methods % get/set
        
        function val = get.LastPosition(obj)
            
            %Pass it on
            val = getpref(obj, 'Last Position', []);
            
        end
        
        function set.LastPosition(obj, val)
            
            %Pass it on
            setpref(obj, 'Last Position', val);
            
        end
        
        function val = get.SerializeRecent(obj)
            
            %Pass it on
            val = getpref(obj, 'Recent Files', {});
            
            %Ensure cellstr
            if isempty(val)
                val = {};
            end
            
        end
        
        function set.SerializeRecent(obj, val)
            
            %Ensure cellstr
            assert(iscellstr(val), 'must be a list of strings');
            
            %No point have anything appear twice on the list
            val = unique(val, 'stable');
            
            %Truncate, if necessary
            val(obj.SerializeRecentMax+1:end) = [];
            
            %Pass it on
            setpref(obj, 'Recent Files', val);
            
        end
        
        function val = get.SerializeName(obj)
            
            %Start here
            val = obj.SerializeName;
            
            %But ensure we send something back
            if isempty(val)
                val = class(obj);
            end
            
            %Ensure whatever we send back is a valid variable name
            val = matlab.lang.makeValidName(val);
            
        end
        
        function val = get.SerializeFile(obj)
            
            %If nothing
            if isempty(obj.SerializeFile)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.SerializeFile;
                
            end
            
        end
        
        function val = get.SerializePath(obj)
            
            %If nothing
            if isempty(obj.SerializePath)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.SerializePath;
                
            end
            
        end
        
        function val = get.SerializeFullFile(obj)
            
            %Combine and force to char
            val = fullfile(obj.SerializePath, obj.SerializeFile);
            
        end
        
        function set.SerializeFullFile(obj, val)
            
            %Split
            [p, n, e] = fileparts(val);
            
            %Make a note
            obj.SerializePath = p;
            obj.SerializeFile = [n, e];
            
        end
        
        function val = get.SerializeFullFileMRU(obj)
            
            %Same as other method
            val = obj.SerializeFullFile;
            
        end
        
        function set.SerializeFullFileMRU(obj, val)
            
            %Start here
            obj.SerializeFullFile = val;
            
            %Update MRU list with most recent file FIRST
            mru = obj.SerializeRecent;
            mru = [{val}, mru];
            obj.SerializeRecent = mru;
            
        end
        
    end
    
    methods % construction/destruction
        
        function obj = Serializable
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, -Inf, ... % specify a super-high priority, so this appears at start
                    'File', 'filemenu', ...
                    'File>New...!N', 'new', ...
                    'File>Open...!O', 'open', ...
                    'File>Recent', [], ...
                    'File>Close!W', 'close', ...
                    'File>Close All', 'closeall', ...
                    'File>Save!S', 'save', ...
                    'File>Save As...', 'saveAs');
            end
            
        end
        
        function delete(obj)
        
            %Treat same as close, with the 'force' flag raised
            close(obj, true);
            
        end
        
    end
    
    methods % preference management
        
        function varargout = getpref(obj, varargin)
            
            %Ensure preference name(s) are valid
            if nargin > 1
                varargin{1} = matlab.lang.makeValidName(varargin{1});
            end
            
            %Pass it on
            [varargout{1:nargout}] = getpref(obj.SerializeName, varargin{:});
            
        end
        
        function varargout = setpref(obj, varargin)
            
            %Ensure preference name(s) are valid
            if nargin > 1
                varargin{1} = matlab.lang.makeValidName(varargin{1});
            end
            
            %Pass it on
            [varargout{1:nargout}] = setpref(obj.SerializeName, varargin{:});
            
        end
        
    end
    
    methods % menu management
        
        function filemenu(obj, varargin)
            
            %Need the calling handle (is this always safe ?)
            hc = gcbo;
            
            %Anything to worry about ?
            if isempty(hc)
                
                %No
                return;
                
            end
            
            %Maintain MRU, the list of recently used files
            i_recent(hc, 'Recent', obj.SerializeRecent, @i_open);
            
            %Ditto recently imported files (if applicable)
            if isa(obj, 'mvc.mixin.Importable')
                i_recent(hc, 'Reimport', obj.ImportRecent, @i_import);                
            end
            
            function i_recent(hc, lab, nam, fcn)
                
                %Maintain MRU, the list of recently used files
                % Need to cater separately for HG and Java objects
                if ishghandle(hc)
                    
                    %Look for a sub-menu with specified label
                    hm = findobj(hc, 'Label', lab);
                    
                    %Anything to do ?
                    if isempty(hm)
                        
                        %No
                        return;
                        
                    end
                    
                    %Clean sheet
                    delete(allchild(hm));
                    
                    %Anything ?
                    if isempty(nam)
                        
                        %No
                        uimenu(hm, 'Label', '[none]', 'Enable', 'off');
                        
                    else
                        
                        %Ensure enabled
                        set(hm, 'Enable', 'on');
                        
                        %Add to menu
                        for i = 1:numel(nam)
                            uimenu(hm, 'Label', nam{i}, ...
                                'Callback', {fcn, nam{i}});
                        end
                        
                    end
                    
                elseif isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
                    
                    %Look for a sub-menu for recent files
                    labs = arrayfun(@(x)get(hc.getItem(x), 'Text'), 0:hc.getItemCount-1, 'UniformOUtput', false);
                    idx = find(strcmp(lab, labs)) - 1;
                    
                    %Anything to do ?
                    if isempty(idx)
                        
                        %No
                        return;
                        
                    end
                    
                    %Clean sheet
                    jmenu = hc.getItem(idx);
                    if isa(jmenu, 'javax.swing.JMenu')
                        
                        %Clean sheet
                        for i = jmenu.getItemCount()-1:-1:0
                            jmenu.remove(i);
                        end
                        
                    end
                    
                    %For each recent file
                    for i = 1:numel(nam)
                        
                        %Create a menu item with action
                        jitem = javax.swing.JMenuItem(nam{i});
                        
                        %Assign action
                        set(jitem, 'ActionPerformedCallback', {fcn, nam{i}});
                        
                        %And add to menu
                        jmenu.add(jitem);
                        
                    end
                    
                else
                    warning(['unhandled class ''', class(hc), '''']);
                end
                
            end
            
            function i_open(hc, ~, nam)
                
                %Might take some time
                clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                
                %Careful
                try
                    
                    %Pass it on
                    obj.open(nam);
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj, err, 'modal'));
                    
                end
                
            end
            
            function i_import(hc, ~, nam)
                
                %Might take some time
                clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                
                %Careful
                try
                    
                    %Pass it on
                    obj.import(nam);
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj, err, 'modal'));
                    
                end
                
            end
            
        end
        
    end
    
    methods % session management
        
        function varargout = new(obj, varargin)
            
            %Make a new instance of whatever this is
            fcn = str2func(class(obj));
            
            %Any properties to be copied from old to new ?
            if isempty(obj.CopyOnNew)
                
                %No
                pav = {};
            
            else
                
                %Get from old
                val = get(obj, obj.CopyOnNew);
                
                %Pass them on
                pav = [obj.CopyOnNew; val];
                
            end
            
            %Using any supplied arguments
            [varargout{1:nargout}] = fcn(pav{:}, varargin{:});
            
        end
        
        function [bOpenable, strWhyNot] = isOpenable(obj, fn)
            %
            % Helper function returns true of file FN exists,
            %  and has an extension that matches one defined in the serialize spec
            
            %Hope for the best
            bOpenable = true;
            strWhyNot = [];
            
            %File exists ?
            if exist(fn, 'file') ~= 2
                
                %Simple
                bOpenable = false;
                strWhyNot = 'file does not exist';
                
            else
            
                %Compare file extension with file spec
                [~, ~, e] = fileparts(fn);
                idxx = find(cellfun(@(x)~isempty(strfind(x,e)), obj.SerializeSpec(:,1)));
                
                %Anything ?
                if isempty(idxx)
                    
                    %No
                    bOpenable = false;
                    strWhyNot = 'file extension does not match any listed in serialization spec';
                    
                else
                    
                    %Is this guilding the lily ?  YES, because whos('-file', fn) loads objects,
                    % which takes time and causes problems (around ChangedSinceLastSaved)
                    return;
                    
                    %But if we did want to go into more detail, it could look like this
                    info = whos('-file', fn);
                    b = strcmp('obj', {info.name});
                    if ~any(b)
                        
                        %Not loadable
                        bOpenable = false;
                        strWhyNot = 'file does not contain a variable named ''obj''';
                        
                    elseif ~strcmp(class(obj), info(b).class)
                        
                        %Not loadable
                        bOpenable = false;
                        strWhyNot = ['file does not contain a variable named ''obj'' of class ''', class(obj), ''''];
                        
                    end
                    
                end
                
            end
            
        end
        
        function newobj = open(obj, fn, bForce, varargin)
            
            %Do we need to prompt for save-file ?
            if nargin < 2
                
                %Yes
                bPrompt = true;
                fn = [];
                
            elseif nargin < 3
                
                %Maybe
                bPrompt = isempty(fn);
                
            else
                
                %Maybe
                bPrompt = ~bForce;
                
            end
            
            %Check for mru option
            [bNoMRU, varargin] = checkOption('-nomru', false, varargin{:});
            
            %If caller specifies file as numeric
            if ~isempty(fn) && isnumeric(fn)
                
                %Assume it to be an index into the Recent File list
                fn = obj.SerializeRecent{fn};
                
            end
            
            %Helper function shared across classes and methods
            [fn, idx] = filepicker(fn, obj.SerializeSpec(:,1:2), bPrompt, false, dlgtitle(obj));
            if isempty(fn)
                
                %Bail out
                newobj = [];
                return;
                
            end
            
            %Check that file exists
            if exist(fn, 'file') ~= 2
                
                %Ensure this file is not on recent list
                mru = obj.SerializeRecent;
                mru(strcmpi(fn, mru)) = [];
                obj.SerializeRecent = mru;
                
                %Throw
                error(['file ''', fn, ''' does not exist']);
                
            end
            
            %Pull index-specific load function from serialzation spec
            fcn = obj.SerializeSpec{idx,4};
            
            %Anything ?
            if isempty(fcn)
                
                %No - just use standard load function, looking for content named 'obj'
                TEMP = load(fn, '-mat', 'obj');
                
                %Validate content
                assert(isfield(TEMP, 'obj'), 'selected file does not contain a variable named ''obj''');
                assert(isa(TEMP.obj, class(obj)), ['selected file does not contain a variable named ''obj'' of class ', class(obj)]);
                
                %In which case
                newobj = TEMP.obj;
                
            else
                
                %Pass it on, assuming open function wants filename as arg 1
                newobj = fcn(fn, varargin{:});
                
            end
            
            %However newobj was arrived at, ensure consistency with any 'CopyOnNew' properties
            if ~isempty(obj.CopyOnNew)
                
                %Get from old
                val = get(obj, obj.CopyOnNew);
                
                %Pass them on
                pav = [obj.CopyOnNew; val];
                set(newobj, pav{:});
                
            end
            
            %If caller does NOT want new object back
            if nargout == 0
                
                %If nothing has changed since last save
                if ~obj.ChangedSinceLastSave
                    
                    %No worries - crack on
                    
                elseif nargin > 2 && bForce
                    
                    %Caller is specifically saying force open regardless of unsaved changes
                    
                elseif ~obj.confirm('Continue WITHOUT saving changes to the current model ?')
                    
                    %Bail out
                    return;
                    
                end
                
                %Overwrite THIS object with new, but can't just do obj = TEMP.obj,
                % instead get whatever is gettable
                prp = allProperties(obj, true);
                
                %If collectable
                if isa(obj, 'mvc.mixin.Collectable')
                    
                    %Not interested in Parent, and Children are handled separately
                    prp = setdiff(prp, {'Parent', 'Children'});
                    
                    %Any children ?
                    if newobj.HasChildren
                        
                        %Pull them out before doing anything else
                        C = detach(newobj.Children);
    
                    else
                        
                        %Nothing to worry about
                        C = [];
                        
                    end
                    
                end
                
                %Limit to only those things that are also settable
                prp(~isEditable(obj, prp{:}, '-excludeUngettables')) = [];
                
                %Get corresponding values from object loaded from file
                val = get(newobj, prp);
                
                %Backward compatibility with R2015b
                if verLessThan('matlab', '9.3')
                    
                    %Assigning a value that is an empty object array throws error
                    %
                    % Cannot call method 'set' because one or more inputs
                    % of class 'awi.model.Entity' are heterogeneous and
                    % 'set' is not sealed.  For more details please see the
                    % method dispatching rules for heterogeneous arrays.
                    %
                    %Not a problem in 2017b - looks like a legacy bug in mcos
                    b = cellfun(@(x)isempty(x) && isa(x, 'mvc.mixin.Collectable'), val);
                    [val{b}] = deal([]);
                    
                end
                                
                %Assign in this object
                pav = [prp; val];
                set(obj, pav{:});
                
                %Now handle Children, gracefully
                if isa(obj, 'mvc.mixin.Collectable') && ~isempty(C)
                    obj.Children = C;
                end
                
                %If we got this far, make a note of file details
                if bNoMRU
                    obj.SerializeFullFile = fn;
                else
                    obj.SerializeFullFileMRU = fn;
                end
                obj.SerializeIndex = idx;
                
                %And this must be true
                obj.ChangedSinceLastSave = false;
                
                %Might seem odd that this is necessary, but without it when newobj goes
                % out of scope it would otherwise prompt to save changes
                newobj.ChangedSinceLastSave = false;
                
                %Are we bothered about view arrangement too ?
                if obj.SerializeViewArrangement
                    
                    %Find a view manager associated with this object
                    mgr = viewManager(obj);
                    
                    %Anything ?
                    if ~isempty(mgr)
                        
                        %Pass it on
                        mgr.ViewArrangement = fn;
                        
                    end
                
                end
                
            else
                
                %Ensure content all up to date
                if bNoMRU
                    newobj.SerializeFullFile = fn;
                else
                    newobj.SerializeFullFileMRU = fn;
                end
                newobj.SerializeIndex = idx;
                
                %And this must be true
                newobj.ChangedSinceLastSave = false;
                
            end
            
        end
        
        function varargout = close(obj, varargin)
            
            %Pass it on to view manager
            mgr = viewManager(obj);
            [varargout{1:nargout}] = close(mgr, varargin{:});
            
        end
        
        function closeall(obj, varargin)
            
            %Find ALL relevant peer objects
            mgr = viewManager(obj, 'all');
            
            %For each in turn
            for i = 1:numel(mgr)
                
                %Pass it on
                close(mgr(i), varargin{:});
        
            end
            
        end
        
        function bCanClose = canClose(obj, bForce, varargin)
            
            %Assume the worst
            bCanClose = false;
           
            %Allow for hg callback
            if nargin > 2 && ishghandle(bForce)
                
                %Pull nargin forwards
                if numel(varargin) > 1
                    bForce = varargin{2};
                else
                    bForce = [];
                end
                varargin(1) = [];
                
            end
            
            %If nothing has changed
            if ~obj.ChangedSinceLastSave
                
                %No worries - just crack on
                
            elseif nargin > 1 && ~isempty(bForce) && bForce
                
                %Caller is specifically saying force a close regardless of unsaved changes
                
            elseif ~confirm(obj, 'Close WITHOUT saving changes ?')
                
                %Bail out
                return;
                
            end
            
            %If we get this far
            bCanClose = true;
            
        end
        
        function varargout = save(obj, fn, bPrompt, varargin)
            
            %Do we need to prompt for save-file ?
            if nargin > 2
                
                %Caller is telling us explicitly
                if isempty(fn)
                    
                    %Go with current filename
                    fn = obj.SerializeFullFile;
                    
                    %But nothing is not really acceptable
                    if isempty(fn)
                        
                        %So override whatever caller told us (?? ok ??)
                        bPrompt = true;
                        
                    end
                    
                end
                
            elseif nargin > 1 && ~isempty(fn)
                
                %Caller is supplying a filename, so prompt is not needed
                bPrompt = false;
                
            else
                
                %Go with current filename
                fn = obj.SerializeFullFile;
                
                %And prompt if this is not yet specified
                bPrompt = isempty(fn);
                
            end
            
            %Check for mru option
            [bNoMRU, varargin] = checkOption('-nomru', false, varargin{:});
            
            %Helper function shared across classes and methods
            [fn, idx] = filepicker(fn, obj.SerializeSpec(:,1:2), bPrompt, true, dlgtitle(obj));
            if isempty(fn)
                
                %Bail out
                [varargout{1:nargout}] = deal([]);
                return;
                
            end
            
            %Pull index-specific save function from serialzation spec
            fcn = obj.SerializeSpec{idx,3};
            
            %Anything ?
            if isempty(fcn)
                
                %No - just save this object and its collection to specified file
                save(fn, '-mat', 'obj');

                %Are we bothered about view arrangement too ?
                if obj.SerializeViewArrangement
                    
                    %If invoked from a view derived from ViewManager
                    mgr = viewManager(obj);
                    if ~isempty(mgr)
                        
                        %Get the arrangement of views
                        ViewArrangement = mgr.ViewArrangement; %#ok<NASGU>
                        
                        %Append to file
                        save(fn, '-mat', '-append', 'ViewArrangement');
                        
                    end
                
                end
                
            else
                
                %Pass it on, assuming save function wants object and filename as input args
                [varargout{1:nargout}] = fcn(obj, fn, varargin{:});
                
            end
            
            %If we got this far,make a note of file details
            if bNoMRU
                obj.SerializeFullFile = fn;
            else
                obj.SerializeFullFileMRU = fn;
            end
            obj.SerializeIndex = idx;
            
            %And this must now be true
            obj.ChangedSinceLastSave = false;
            
        end
        
        function varargout = saveAs(obj, fn, varargin)
            
            %If no file specified
            if nargin < 2
                
                %Go with empty
                fn = [];
                
            end
            
            %Same as 'save' but with the 'prompt' flag set true
            [varargout{1:nargout}] = save(obj, fn, true, varargin{:});
            
        end
        
    end
    
    methods (Access = public)
        
        function mgr = viewManager(obj, hf)
            %
            % Return ViewManager MGR associated with instance OBJ and figure HF.
            % If HF not specified, use current figure
            
            %Needed later
            bAll = false;
            
            %Assume the worst
            mgr = [];
            
            %Figure not specified ?
            if nargin < 2
                
                %NO - this approach doesn't work when using programmatic interface
                %                 %Go with current figure (but using the custom function allowing for Java widgets)
                %                 hf = mvc.mixin.UiTools.currentFigure(gcbo);
                %
                %                 %Anything ?
                %                 if isempty(hf)
                %
                %                     %No
                %                     return;
                %
                %                 end
            
                %Better to look at all figures
                temp = get(0, 'ShowHiddenHandles');
                set(0, 'ShowHiddenHandles', 'on');
                hf = findobj(0, 'Type', 'figure');
                set(0, 'ShowHiddenHandles', temp);
            
            elseif ischar(hf) && strcmpi('all', hf)
                
                %Look at ALL figures
                bAll = true;
                temp = get(0, 'ShowHiddenHandles');
                set(0, 'ShowHiddenHandles', 'on');
                hf = findobj(0, 'Type', 'figure');
                set(0, 'ShowHiddenHandles', temp);
            
            elseif ishghandle(hf) || isa(hf, 'mvc.view.Container')
                
                %Get parent figure
                hf = ancestor(hf, 'figure');
                
            else
                error('bad input');
            end
            
            %Placeholder
            mgr = cell(size(hf));
            
            %For each figure handle
            for i = 1:numel(hf)
                
                %Get appdata from figure
                appdata = getappdata(hf(i));
                
                %Anything ?
                if isempty(appdata)
                    
                    %No
                    continue;
                    
                end
                
                %All, or just the one ?
                if bAll
                    
                    %Look for a ViewManager in appdata associated with ANY model
                    b = structfun(@(x)isobject(x) && isvalid(x) && isa(x, 'mvc.view.ViewManager') && isa(x.Model, class(obj)), appdata);
                    
                else
                    
                    %Look for a ViewManager in appdata associated with THIS model
                    %   - TODO : We get an error to do with logical scalar
                    %   values here.
                    b = structfun(@(x)isobject(x) && isvalid(x) && isa(x, 'mvc.view.ViewManager') && x.Model == obj, appdata);
                    %   - This fails because 'appdata' is taken from all
                    %   figures that matlab has open. So if the user has
                    %   any figures open other than the AWI framework the
                    %   logical check above will fail. 
                    %
                    %   - Fix. Add a simpler logical check first.
                end
                
                %Anything ?
                if any(b)
                    
                    %Get it
                    fld = fieldnames(appdata);
                    mgr{i} = appdata.(fld{b});
                    
                end
            
            end
            
            %Expand
            mgr = [mgr{:}];
            
        end
        
    end
    
end
