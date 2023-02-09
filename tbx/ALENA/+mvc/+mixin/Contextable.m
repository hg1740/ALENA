classdef (ConstructOnLoad) Contextable < matlab.mixin.SetGet
    %
    % Contextable provides functionality allowing objects, and object arrays, to present context-sensitive options.
    %
    % Examples: (where fcn is a constructor of a class that implements Contextable)
    %
    %  obj = fcn('Name', 'Peter', 'Type', 'Employee', 'Description', {'Peter is an employee'});
    %  obj.addContext('Label', cb);
    %  obj.context
        
    properties (Access = protected, Transient)
        
        %Storage for context details
        Context = struct('Label', {}, 'Callback', {}, 'Children', {}, 'Priority', {});
        
    end
    
    methods (Sealed) % Helper functions acting on a heterogrenous array of Contextables
        
        function varargout = context(obj, hc, ctx, bClean, varargin)

            %Maybe needed later
            bFigure = false;
            bJava = false;

            %What have we got ?
            if ~isvalid(obj)
                
                %Nothing
                return;
                
            elseif nargin == 1
                
                %Nothing at all
                hcm = [];
                
            elseif isempty(hc)
                
                %Nothing to do
                return;
                
            elseif ~all(ishandle(hc))
                
                %Not allowed
                error('arg1 must be a handle');
             
            elseif numel(hc) == 1 && isa(hc, 'matlab.ui.Figure')
                
                %Figure handle passed in directly
                bFigure = true;
                bClean = false;
                
            elseif numel(hc) == 1 && isa(hc, 'com.mathworks.hg.peer.utils.UIMJTree')
                
                %Java tree control passed in - nothing more to do
                hcm = hc;
                bJava = true;
                
            elseif numel(hc) == 1 && isa(hc, 'javahandle_withcallbacks.com.mathworks.hg.peer.UITreePeer')
                
                %MATLAB tree control passed in - get the underlying Java tree
                hcm = hc.Tree;
                bJava = true;
                
                %AND make a note of the originating figure (maybe needed later)
                obj.currentObject(hc.Parent);
            
            elseif numel(hc) == 1 && isa(hc, 'javax.swing.JMenu')
                
                %Java menu item passed in
                hcm = hc;
                bJava = true;
                
                %AND make a note of the originating figure (maybe needed later)
                %obj.currentObject(hc.getParent);
             
            elseif numel(hc) == 1 && isa(hc, 'javax.swing.JMenuItem')
                
                %Java menu item passed in
                hcm = hc;
                bJava = true;
           
            elseif numel(hc) == 1 && isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
                
                %Java menu item passed in
                hcm = hc;
                bJava = true;
                
            elseif numel(hc) == 1 && any(strcmpi(get(hc, 'Type'), {'UiContextMenu', 'UiMenu'}))
                
                %Passed in - nothing more to do
                hcm = hc;
                                
            else
                
                %If multiple handles being passed in
                if numel(hc) > 1
                    
                    %Get all ancestor figures
                    hff = arrayfun(@(x)ancestor(x, 'Figure'), hc, 'UniformOutput', false);
                    hff = unique([hff{:}]);
                    
                    %TODO: Could handle this more gracefully
                    assert(numel(hff) == 1, 'all handles must be parented by same figure');
                    
                else
                    
                    %Just get the figure
                    hff = ancestor(hc(1), 'Figure');
                    
                end
                
                %Create menu accordingly
                hcm = uicontextmenu('Parent', hff);
                set(hc, 'UiContextMenu', hcm);

            end
            
            %Clean sheet - optional
            if nargin < 4 || isempty(bClean)
                bClean = true;
            end
            
            %Clean, if necessary (not applicable to Figures or Java objects)
            if bFigure || bJava
                %Clean not applicable
            elseif isempty(hcm) || ~ishandle(hcm)
                %Nothing to do anyway
            elseif bClean
                delete(allchild(hcm));
            end

            %Caller may pass in context detail
            if nargin > 2
                
                %What have we got ?
                if isstruct(ctx)
                
                    %Just go with whatever passed in
    
                elseif ischar(ctx)
                    
                    %Use to down-select from built-in context
                    ctx = obj.Context(strcmp(ctx, {obj.Context.Label})).Children;
                    
                else
                    error('bad input');
                end
                    
            elseif isempty(obj)
                
                %Nothing else to do
                return;
                
            elseif numel(obj) == 1
                
                %Just the one - easy
                ctx = obj.Context;
                
            else
                
                %COULD just get all the options
                % ctx = [obj.Context];
                %
                % But ctx would then be a structure array, with many repeating entries.  Better to go with the common options, perhaps ?
                ctx = i_common(obj.Context);
                
            end
            
            %If originally called with nothing
            if nargin == 1
                
                %Force display of context options now
                [varargout{1:nargout}] = i_display(ctx);
             
            elseif bFigure
                
                %Context menu applied directly to a figure window
                i_figure(hc, ctx);
                
            else
            
                %If using Java menu
                if bJava
                    
                    %Menu already initialised ?
                    if isa(hc, 'javax.swing.JMenu') || isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
                        
                        %Yes
                        jmenu = hcm;
                        
                        %Clean required ?
                        if bClean
                            if jmenu.getItemCount > 0
                                jmenu.removeAll;
                            end
                        end
                        
                    else
                        
                        %No - initialise
                        jmenu = i_initJava(hcm);
                        
                    end
                    
                else
                    
                    %Nothing to do
                    jmenu = [];
                    
                end
                
                %Create context menu
                [varargout{1:nargout}] = i_create(hcm, ctx, jmenu);
                
            end
            
            function ctx = i_common(varargin)
            
                %Which labels are common to all inputs ?
                lab = {varargin{1}.Label};
                for i = 2:numel(varargin)
                    lab = intersect(lab, {varargin{i}.Label});
                end
                   
                %Down-select accordingly
                ctx = varargin{1}(ismember(lab, {varargin{1}.Label}));
                
                %There's more to do here, because of hierarchy and because not all methods will be applicable to multiple instances
                
            end
            
            function [str, fcn] = i_display(ctx, str, fcn, prefix)
                
                %Starting point supplied ?
                if nargin < 2
                    
                    %No - initialise
                    str = {};
                    fcn = {};
                    prefix = '';
                    
                end
                
                %For each element of structure array
                for i = 1:numel(ctx)
                
                    %Extend string
                    str{end+1} = [prefix, ctx(i).Label];
                    fcn{end+1} = ctx(i).Callback;
                    
                    %Any submenus required ?  Handle by recursion
                    [str, fcn] = i_display(ctx(i).Children, str, fcn, [' ', prefix]);
                    
                end
                
                %Are we done ?
                if nargout == 0
                
                    %If nothing available
                    if isempty(str)
                        
                        %Must be ok to display a modal dialog, because we were about to anyway with listdlg
                        uiwait(warndlg(obj, 'No context available'));
                        
                    else
                        
                        %Offer to user
                        sel = listdlg(obj, str);
                        
                        %Anything ?
                        if isempty(sel)
                            
                            %No
                            
                        else
                            
                            %Do it
                            i_eval([], [], fcn{sel});
                            
                        end
                    
                    end
                    
                end
                
            end
            
            function jmenu = i_initJava(hcm)
                
                %Start here
                jmenu = javax.swing.JPopupMenu;
                
                %Get handle to callback properties
                hcb = handle(hcm, 'CallbackProperties');
                
                % Set the tree mouse-click callback
                % Note: MousePressedCallback is better than MouseClickedCallback
                %       since it fires immediately when mouse button is pressed,
                %       without waiting for its release, as MouseClickedCallback does
                set(hcb, 'MousePressedCallback', {@mousePressedCallback,jmenu});
                
                % Set the mouse-press callback
                function mousePressedCallback(hTree, eventData, jmenu)
                    if eventData.isMetaDown  % right-click is like a Meta-button
                        % Get the clicked node
                        clickX = eventData.getX;
                        clickY = eventData.getY;
                        jtree = eventData.getSource;
                        treePath = jtree.getPathForLocation(clickX, clickY);
                        try
                            % Modify the context menu or some other element
                            % based on the clicked node. Here is an example:
                            node = treePath.getLastPathComponent;
                            nodeName = ['Current node: ' char(node.getName)];
                            item = jmenu.add(nodeName);
                            
                            % remember to call jmenu.remove(item) in item callback
                            % or use the timer hack shown here to remove the item:
                            timerFcn = {@removeItem,jmenu,item};
                            start(timer('TimerFcn',timerFcn,'StartDelay',0.2));
                        catch
                            % clicked location is NOT on top of any node
                            % Note: can also be tested by isempty(treePath)
                        end
                        
                        % Display the (possibly-modified) context menu
                        jmenu.show(jtree, clickX, clickY);
                        jmenu.repaint;
                        
                    end
                end
                
                % Remove the extra context menu item after display
                function removeItem(hObj,eventData,jmenu,item)
                    jmenu.remove(item);
                end
                
            end
            
            function i_figure(hf, ctx)
                %
                % Applying content to figure window directly is slightly special
                
                %What order to present context items ?  Based on Priority
                [~, odx] = sort([ctx.Priority]);
                
                %For each in turn
                for i = odx
                    
                    %If this item has no children
                    if isempty(ctx(i).Children)
                        
                        %Ignore it
                        continue;
                        
                    end
                    
                    %Label needs slight tweaking for main uimenu
                    lab = ctx(i).Label;
                    if lab(1) == '|'
                        lab(1) = [];
                    end
                    lab = ['&', lab];
                    
                    %Look for existing uimenu with same name
                    hm = findobj(allchild(hf), 'Type', 'uimenu', 'Label', lab);
                    
                    %Anything ?
                    if isempty(hm)
                        
                        %No - create
                        hm = uimenu('Parent', hf, ...
                            'Label', lab, ...
                            'Callback', {@i_eval, ctx(i).Callback});
                        
                    else
                        
                        %Clean it
                        delete(allchild(hm));
                        
                        %But ensure callback is set
                        set(hm, 'Callback', {@i_eval, ctx(i).Callback});
                        
                    end
                    
                    %Now pass it on
                    i_create(hm, ctx(i).Children, []);
                    
                    %If we are dealing with the 'File' menu, and the figure
                    % has a good old fashioned standard menu bar
                    if strcmpi(lab, '&File') && strcmpi('figure', get(hf, 'MenuBar'))
                    
                        %Link menubar callbacks to children of the File context menu
                        i_menubar(hf, ctx(i).Children);
                        
                    end
                    
                end
                
            end
            
            function i_menubar(hf, ctx)
            
                %Get existing menubar
                ha = allchild(hf);
                b = arrayfun(@(x)isa(x, 'matlab.ui.container.Toolbar'), ha);
                if ~any(b)
                    return;
                end
                hb = ha(b);
                
                %Which toolbar buttons are we interested in ?
                tag = {'Standard.NewFigure', 'Standard.FileOpen', 'Standard.SaveFigure'};
                hp = cellfun(@(x)findobj(allchild(hb), 'Tag', x), tag, 'UniformOutput', false);
               
                %Which menu functions do these map to ?
                map = {'New...', 'Open...', 'Save'};
                idx = cellfun(@(x)find(~cellfun(@isempty, strfind({ctx.Label}, x)), 1, 'first'), map);
                
                %Corresponding toolip strings (TODO: find a better way)
                tts = {'Create a new session', 'Open a previously saved session', 'Save current session to file'};
                
                %Link the pushtools to callbacks stored in context detail
                for i = find(idx > 0)
                    set(hp{i}, 'ClickedCallback', {@i_eval, ctx(idx(i)).Callback}, ...
                        'TooltipString', tts{i});
                end
                
            end
            
            function i_create(hm, ctx, jmenu)
                
                %What order to present context items ?  Based on Priority
                [~, odx] = sort([ctx.Priority]);
                
                %For each in turn
                for i = odx
                
                    %Separator ?
                    if ctx(i).Label(1) == '|'
                        
                        %Yes
                        sep = 'on';
                        lab = ctx(i).Label(2:end);
                        
                    else
                        
                        %No
                        sep = 'off';
                        lab = ctx(i).Label;
                        
                    end
                    
                    %Accelerator ?
                    idx = find(lab == '!', 1, 'last');
                    if isempty(idx)
                        
                        %Nothing
                        acc = '';
                        
                    else
                        
                        %Get it                    
                        acc = lab(idx+1:end);
                        lab(idx:end) = [];                        
                        
                    end
                    
                    %Enabled ?
                    bEn = ~isempty(ctx(i).Children) || ~isempty(ctx(i).Callback);
                    
                    %Create menu item
                    if isempty(jmenu)
                        
                        %Using regular handle-graphics
                        hmm = uimenu('Parent', hm, ...
                            'Label', lab, ...
                            'Separator', sep, ...
                            'Accelerator', acc, ...
                            'Enable', mvc.mixin.UiTools.bool2offon(bEn), ...
                            'Callback', {@i_eval, ctx(i).Callback});
                        
                        %Ensure parent menu item enabled
                        if isa(hm, 'matlab.ui.container.Menu')
                            set(hm, 'Enable', 'on');
                        end
                        
                        %In which case
                        jitem = [];
                    
                    else
                        
                        %Need separator ?
                        if strcmpi(sep, 'on')
                            
                            %Yes
                            jmenu.addSeparator;
                            
                        end
                        
                        %TODO PROPERLY: SPECIAL CASE of File->Recent
                        bSpecial = strcmp(lab, 'Recent') && strcmp(jmenu.getText, 'File');
                        
                        %TODO PROPERLY: SPECIAL CASE of File->Reimport
                        bSpecial(end+1) = strcmp(lab, 'Reimport') && strcmp(jmenu.getText, 'File');
                        
                        %TODO PROPERLY: SPECIAL CASE of Edit->Add
                        bSpecial(end+1) = strcmp(lab, 'Add') && strcmp(jmenu.getText, 'Edit');
                        
                        %Using Java menu, if this item has NO children and is NOT on of the special cases
                        if isempty(ctx(i).Children) && ~any(bSpecial)
                            
                            %Create a menu item with action
                            jitem = javax.swing.JMenuItem(lab);
                            
                            %If callback action specified
                            if isempty(ctx(i).Callback)
                                
                                %Disable it
                                set(jitem, 'Enable', false);
                                
                            else
                                
                                %Assign action
                                set(jitem, 'ActionPerformedCallback', {@i_eval, ctx(i).Callback});
                                
                            end
    
                        else
                            
                            %Create a menu (to which actual items are added later)
                            jitem = javax.swing.JMenu(lab);
                            
                            %Menu might also have an action (e.g. to control enable status of children)
                            if ~isempty(ctx(i).Callback)
                            
                                %Assign action
                                set(jitem, 'MenuSelectedCallback', {@i_eval, ctx(i).Callback});

                            end
                            
                        end
                        
                        %And add to menu
                        jmenu.add(jitem);

                        %In which case
                        hmm = [];
                        
                    end
                    
                    %Any submenus required ?
                    if ~isempty(ctx(i).Children)
                        
                        %Handle by recursion
                        i_create(hmm, ctx(i).Children, jitem);
    
                    end
                    
                end
                
            end

            function i_eval(hc, ~, fnam, varargin)
                
                %There might be nothing to do
                if nargin < 3 || isempty(fnam)
                    
                    %Bail
                    return;
                    
                end
                
                %Might take some time
                clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                
                %Allow for fnam and varargin to be specified by caller as a single cell array
                if iscell(fnam) && isempty(varargin)
                    varargin = fnam(2:end);
                    fnam = fnam{1};
                end
                
                %Careful
                try
                    
                    %Do what ?
                    if ischar(fnam)
                        
                        %Convert to function handle
                        fcn = str2func(fnam);
                        
                    elseif isa(fnam, 'function_handle')
                        
                        %Just do whatever
                        fcn = fnam;
                        
                    else
                        error('argument must either be a method name or function handle');
                    end
                    
                    %Do it, treating it as a method of the object
                    fcn(obj, varargin{:});
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj, err, 'modal'));
                   
                end
                
            end
            
        end
        
        function addContext(obj, varargin)
            
            %Caller MAY specify priority as first arg
            if isnumeric(varargin{1})
                
                %Pull it out
                pr = varargin{1};
                varargin(1) = [];
                bForce = true;
                
            else
                
                %Go with a nominal value for priority
                % Can we make a good guess based on super-class hierarchy ?
                % Would like to use evalin('caller', 'mfilename('class')')
                % but doesn't seem to work.
                %
                % This is extremely tedious
                dbs = dbstack('-completenames');
                dbs(3:end) = [];
                [p,n] = fileparts(dbs(end).file);
                tok = strsplit(p, filesep);
                tok(cellfun(@isempty, tok)) = [];
                tok = tok(cellfun(@(x)x(1) == '+', tok));
                tok = cellfun(@(x)x(2:end), tok, 'UniformOutput', false);
                tok{end+1} = n;
                cls = strjoin(tok, '.');

                %Now look into metaclass list
                mc = metaclass(obj);

                %Pass it on to helper function
                pr = i_findInSuperclassList(cls, mc);
                bForce = false;
                
            end
            
            %Where the recursive helper function looks like this
            function pr = i_findInSuperclassList(cls, mc, pr)

                %Starting point
                if nargin < 3
                    pr = 1;
                end
                
                %Anything ?
                if isempty(mc)
                    
                    %No
                    return;
                    
                end
                
                %Present at this level ?
                if ismember(cls, {mc.Name})
                    
                    %We're done
                    return;
                    
                else
                    
                    %Recurse
                    pr = i_findInSuperclassList(cls, vertcat(mc.SuperclassList), pr + 1);
                    
                end
                
            end
            
            %For each pair of arguments
            for i = 1:2:numel(varargin)
                
                %Extend context store for this object
                i_addContext(obj, varargin{i}, varargin{i + 1}, pr, bForce);
                
            end
            
            function i_addContext(obj, str, fcn, pr, bForce)
        
                %Pass it on
                obj.Context = obj.addToContextStruct(obj.Context, str, fcn, pr, bForce);
                                
            end
            
        end
        
        function rmContext(obj, varargin)
        
            %For each supplied argument
            for i = 1:numel(varargin)
                
                %%Pass it on
                obj.Context = obj.rmFromContextStruct(obj.Context, varargin{i});
                
            end
            
        end
        
        function C = getContext(obj, varargin)

            %Start here
            C = obj.Context;

            %No down-selection ?
            if nargin == 1
                
                %Send back everything
                return;
                
            end

            %Down-select by label
            if iscellstr(varargin)
                
                %Look for match
                [b, cdx] = ismember(varargin, {obj.Context.Label});
            
                %Apply down-selection
                C = C(cdx(b));
            
            end
            
        end
        
    end
    
    methods (Static)

        function S = rmFromContextStruct(S, tok)
            
            %If string passed in
            if ischar(tok)
                
                %Tokenise, allowing the character '>' to be used to define hierarchy in the context structure array
                tok = strsplit(tok, '>');
                
            end
            
            %Anything to do ?
            if isempty(S)
                
                %No
                return;
                
            end
            
            %Look for match with label
            [b, idx] = ismember(tok{1}, {S.Label});
            
            %Nothing ?
            if ~any(b)
                
                %Bail
                return;
                
            elseif numel(tok) == 1
                
                %Strip this value
                S(idx) = [];
                
            else
                
                %Recurse
                S(idx).Children = mvc.mixin.Contextable.rmFromContextStruct(S(idx).Children, tok(2:end));
                
            end
            
        end
        
        function S = addToContextStruct(S, tok, fcn, pr, bForce)
            
            %Defaults
            if nargin < 4
                pr = 1;
                bForce = false;
            end
            
            %Anything yet ?
            if isempty(S)
                
                %No - initialise
                S = struct('Label', {}, 'Callback', {}, 'Children', {}, 'Priority', {});
                
            end
            
            %If string passed in
            if ischar(tok)
                
                %Tokenise, allowing the character '>' to be used to define hierarchy in the context structure array
                tok = strsplit(tok, '>');
                
            end
            
            %Look for existing match
            idx = find(strcmp(tok{1}, {S.Label}));
            
            %Nothing ?
            if isempty(idx)
                
                %Extend structure
                idx = numel(S) + 1;
                
                %Assign label
                S(idx).Label = tok{1};
                
                %Assign priority
                S(idx).Priority = pr;
                
            elseif bForce
                
                %Caller is saying assign specified priority,
                % possibly overwriting whatever is already set
                S(idx).Priority = pr;
                
            end
            
            %Allow for hierarchy
            if numel(tok) == 1
                
                %This point is bottom of hierarchy
                S(idx).Callback = fcn;
                
            else
                
                %Recurse
                S(idx).Children = mvc.mixin.Contextable.addToContextStruct(S(idx).Children, tok(2:end), fcn, pr, bForce);
                
            end
            
        end
        
    end
    
end
