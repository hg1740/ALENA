classdef (ConstructOnLoad) UiTools < matlab.mixin.SetGet
    %
    % Implements a number of basic helper dialogs, mainly thin wrappers
    %  around standard MATLAB functions, such as msgbox, warndlg, listdlg, inputdlg, etc.
    %  The wrapper does little more than set the title of the helper dialog.
    %
    % Examples: where fcn is a class that implements UiTools
    %
    %  obj = fcn('Name', 'Peter', 'Type', 'Employee', 'Description', {'Peter is an employee'});
    %   val = inputdlg(obj)
    %   val = inputdlg(obj, 1)
    %   val = inputdlg(obj, 1, 'Enter a value:')
    %   val = inputdlg(obj, {1,2}, {'Enter a value:', 'Enter another value:'})
    %
    % Also implements a slightly thicker wrapper around standard MATLAB waitbar:
    %
    %   b = waitbar(obj) returns B as true/false if waitbar associated with this object exists or not.
    
    % ASSUMES existence of method DLGTITLE, which returns a meaningful string used as name assigned to helper dialog windows.
    
    properties (Access = private, Transient) % for internal use only
        
        %Handle to waitbar
        hWaitbar;
        
    end
    
    methods % construction / destruction
        
        function delete(obj)
            
            %Suicide pact with waitbar associated with this object
            waitbar(obj, 'close');
            
        end
        
    end
    
    methods (Sealed) % Thin wrappers around standard MATLAB functions
        
        function varargout = errordlg(obj, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified error';
                
            elseif isa(str, 'MException')
                
                %Error supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            [varargout{1:nargout}] = errordlg(str, tit, varargin{:});
            
        end
        
        function varargout = warndlg(obj, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified warning';
                
            elseif isa(str, 'MException')
                
                %Warning supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            [varargout{1:nargout}] = warndlg(str, tit, varargin{:});
            
        end
        
        function varargout = helpdlg(obj, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified help';
                
            elseif isa(str, 'MException')
                
                %Help supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            [varargout{1:nargout}] = helpdlg(str, tit, varargin{:});
            
        end
        
        function varargout = questdlg(obj, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified question';
                
            elseif isa(str, 'MException')
                
                %Question supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            [varargout{1:nargout}] = questdlg(str, tit, varargin{:});
            
        end
        
        function [val, sel] = choose(obj, str, varargin)
            %
            % Thin wrapper around questdlg, or listdlg, depending on number of choices
            
            %Need a message
            if nargin < 2 || isempty(str)
                
                %Reasonable default
                str = 'Choose from the following';
                
            end
            
            %How many choices ?
            if numel(varargin) == 4 && any(strcmp(varargin(1:3), varargin{4}))
                
                %Go with questdlg (user has passed in 3 options, and specified the the default answer)
                val = questdlg(obj, str, varargin{:});
                
            elseif numel(varargin) == 3 && any(strcmp(varargin(1:2), varargin{3}))
                
                %Go with questdlg (user has passed in 2 options, and specified the the default answer)
                val = questdlg(obj, str, varargin{:});
             
            elseif numel(varargin) <= 3
                
                %Go with questdlg (doubling up on last of varargin, as the default answer)
                val = questdlg(obj, str, varargin{:}, varargin{end});
                
            else
                
                %Go with listdlg
                [sel, bOK] = listdlg(obj, varargin, str, 'SelectionMode', 'single');
                
                %Cancelled ?
                if ~bOK
                    
                    %Send back nothing
                    val = [];
                    
                else
                    
                    %So the answer is
                    val = varargin{sel};
                    
                end
                
            end
            
            %Caller may want choice by index back
            if nargout > 1
                
                %Look it up
                sel = find(strcmp(val, varargin));
                
            end
            
        end
        
        function [b, val] = confirm(obj, str, varargin)
            %
            % Thin wrapper around choose
            
            %Need a message
            if nargin < 2 || isempty(str)
                
                %Reasonable default
                str = 'Are you sure';
                
            end
            
            %Need some options
            if nargin < 3 || isempty(varargin)
                
                %Reasonable defaults
                varargin = {'Yes', 'No'};
                
            end
               
            %Pass it on
            val = choose(obj, str, varargin{:});
            
            %And the answer is
            b = isequal(val, varargin{1});
            
        end
            
        function msgdlg(obj, str, varargin)
            %
            % Thin wrapper around msgbox, but assumes user wants the message box to be modal
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified message';
                
            elseif isa(str, 'MException')
                
                %Message supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            uiwait(msgbox(str, tit, varargin{:}, 'modal'));
                        
        end
        
        function varargout = msgbox(obj, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 2
                
                %Nothing much
                str = 'unspecified message';
                
            elseif isa(str, 'MException')
                
                %Message supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Now pass it on
            [varargout{1:nargout}] = msgbox(str, tit, varargin{:});
            
        end
        
        function varargout = listdlg(obj, val, str, varargin)
            %
            % Thin wrapper around builtin
            
            %What have we got ?
            if nargin < 3
                
                %Nothing much
                str = 'select from the following:';
                
            elseif isa(str, 'MException')
                
                %Message supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Pass it on
            [varargout{1:nargout}] = listdlg('ListString', val, ...
                'PromptString', str, ...
                'Name', tit, ...
                'ListSize', [300 300], ... % better default than MATLAB's
                varargin{:});
            
        end
        
        function varargout = sortdlg(obj, val, str, varargin)
            %
            % Thin wrapper around private sortdlg, which is itself
            %  a variation on the theme of the original OTS listdlg
            
            %What have we got ?
            if nargin < 3 || isempty(str)
                
                %Nothing much
                str = 'sort the following:';
                
            elseif isa(str, 'MException')
                
                %Message supplied as an exception object - convert accordingly
                str = err2str(obj, str);
                
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Pass it on
            [varargout{1:nargout}] = sortdlg('ListString', val, ...
                'PromptString', str, ...
                'Name', tit, ...
                varargin{:});
            
        end
        
        function val = inputdlg(obj, val, str, varargin)
            %
            % Thin wrapper around builtin
            
            %Used later
            bReturnSingleValue = false;
            
            %If no value supplied
            if nargin < 2
                
                %Go with empty
                val = [];
                
            end
            
            %If no prompt string supplied
            if nargin < 3
                
                %Go with this nominal prompt
                str = 'Enter a value:';
                
            end
            
            %If caller wants only a single input ?
            if ~iscell(val) && ischar(str)
                
                %Cast to cells
                val = {val};
                str = {str};
                
                %Make a note
                bReturnSingleValue = true;
                
            elseif iscell(val) && iscell(str)
                
                %Assert equal length
                assert(numel(val) == numel(str), 'require same number of inputs as prompts');
                
            else
                error('bad inputs');
            end
            
            %Good title for dialog ?
            tit = dlgtitle(obj);
            
            %Default numlines
            n = ones(size(val));
            n = n(:);
            
            %What is editable in this way ?
            bEmpty = cellfun(@isempty, val);
            bNumeric = cellfun(@isnumeric, val);
            bLogical = cellfun(@islogical, val);
            bChar = cellfun(@ischar, val);
            bCellstr = cellfun(@iscellstr, val);
            n(bCellstr) = cellfun(@numel, val(bCellstr)) + 1;
            bEditable = bNumeric | bLogical | bChar | bCellstr;
            assert(any(bEditable), 'no editable values');
            
            %Cast values accordingly - taking care to handle empty numerics gracefully
            val(bNumeric & ~bEmpty) = cellfun(@mat2str, val(bNumeric & ~bEmpty), 'UniformOutput', false);
            [val{bNumeric & bEmpty}] = deal('');
            val(bLogical) = cellfun(@mat2str, val(bLogical), 'UniformOutput', false);
            val(bCellstr) = cellfun(@(x)strjoin(x, newline), val(bCellstr), 'UniformOutput', false);
            
            %If no modifiers specified
            if isempty(varargin)
                
                %Resize is good
                varargin = {'on'};
                
            end
            
            %Pass it on
            newval = inputdlg(str(bEditable), tit, n(bEditable), val(bEditable), varargin{:});
            
            %Cancelled ?
            if isempty(newval)
                
                %Bail out
                val = {};
                return;
                
            end
            
            %Unpick
            val(bEditable) = newval;
            
            %Cast back - numerics that were NOT originally empty are easy
            val(bNumeric & ~bEmpty) = cellfun(@str2num, val(bNumeric & ~bEmpty), 'UniformOutput', false);
            
            %Numerics that WERE originally empty MIGHT be being set to numeric or char
            bWasEmptyIsNumeric = false(size(val));
            bWasEmptyIsNumeric(bNumeric & bEmpty) = cellfun(@(x)~isempty(x) && (x(1) == '.' || isstrprop(x(1), 'digit')), val(bNumeric & bEmpty));
            
            %Assign numerics accordingly, leaving chars alone
            val(bNumeric & bWasEmptyIsNumeric) = cellfun(@str2num, val(bNumeric & bWasEmptyIsNumeric), 'UniformOutput', false);
            
            %Logical and cellstr are preserved
            val(bLogical) = cellfun(@(x)logical(str2num(x)), val(bLogical), 'UniformOutput', false);
            val(bCellstr) = cellfun(@(x)cellstr(strtrim(x)), val(bCellstr), 'UniformOutput', false);
            
            %Return single value ?
            if bReturnSingleValue
                
                %Yes
                val = val{1};
                
            end
            
        end
        
    end
    
    methods (Sealed) % Slightly thicker wrappers around standard MATLAB functions
        
        function varargout = waitbar(obj, varargin)
            %
            %Wrapper around builtin function of same name
            
            %Does waitbar associated with this object (or any of these objects) currently exist ?
            bExist = arrayfun(@(x)~isempty(x.hWaitbar) && ishandle(x.hWaitbar), obj);
            
            %If called with no further arguments
            if isempty(varargin)
                
                %Caller is asking whether waitbar exists, i.e. if user has cancelled
                varargout{1} = any(bExist);
                
            elseif isequal(varargin, {'close'}) || isequal(varargin, {'delete'})
                
                %Caller is issuing explicit instruction to close existing waitbar(s)
                if any(bExist)
                    delete([obj(bExist).hWaitbar]);
                end
                
                %Might as well keep it clean too
                [obj.hWaitbar] = deal([]);
                
            else
                
                %If numeric value supplied as first argument
                if isnumeric(varargin{1})
                    
                    %Use it as the value of the progress bar
                    val = varargin{1};
                    
                    %Remove val from any remaining inputs
                    varargin(1) = [];
                    
                else
                    
                    %No value given, assume zero (?? could cache internally ??)
                    val = 0;
                    
                end
                
                %The message to display ?
                if isempty(varargin)
                    
                    %Make something up
                    msg = 'Please wait...';
                    
                elseif numel(varargin) == 1 && ischar(varargin{1})
                    
                    %Go with whatever
                    msg = varargin{1};

                elseif ischar(varargin{1})
                    
                    %Hope it is sprintf-able
                    msg = sprintf(varargin{:});
                
                else
                    error('bad input');
                end
                
                %If waitbar already exists
                if any(bExist)
                    
                    %Adjust content
                    if isempty(msg)
                        
                        %Value only
                        arrayfun(@(x)waitbar(val, x), [obj(bExist).hWaitbar]);
                        
                    else
                        
                        %Value and message
                        arrayfun(@(x)waitbar(val, x, msg), [obj(bExist).hWaitbar]);
                        
                    end
                    
                    %Allow for change of title
                    arrayfun(@(x)set(x, 'Name', dlgtitle(obj)), [obj(bExist).hWaitbar]);
                    
                else
                    
                    %Create
                    hwb = waitbar(val, msg, ...
                        'CreateCancelBtn', 'delete(gcbf)', ...
                        'Name', dlgtitle(obj));
                    
                    %This is almost certainly appropriate
                    ha = findobj(hwb, 'type', 'axes');
                    set(get(ha, 'Title'), 'Interpreter', 'none');
                    
                    %Assign to all instances (??)
                    [obj.hWaitbar] = deal(hwb);
                    
                end
                
                %If caller wants something back
                if nargout == 1
                
                    %Send back a clean-up object, so caller need not worry about closing waitbar
                    varargout{1} = onCleanup(@()waitbar(obj, 'close'));
                    
                end
                
            end
            
        end
        
        function hf = dialog(obj, varargin)
            
            %In principle, it's as simple as this
            hf = dialog('Visible', 'off', varargin{:});

            %But can we initialise position of dialog ? If we are a "proper" application object
            if isa(obj, 'mvc.mixin.Serializable')
                
                %AND we can find a view manager
                mgr = viewManager(obj);
                if ~isempty(mgr)
                    
                    %Then we can get a position
                    mpos = mgr.hFigure.Position;
                
                    %Hence centroid
                    xy = mpos(1:2) + mpos(3:4) ./ 2;
                    
                    %Get position of dialog
                    pos = hf.Position;
                    
                    %Adjust accordingly
                    pos(1:2) = [xy - pos(3:4) ./ 2];
                    
                    %Putback
                    hf.Position = pos;
                    
                end
                
            end
            
        end
        
        function [fcn, clu] = progressdlg(obj, tit)
            %
            % Creates a simple dialog comprising a listbox and a Cancel button.
            % Returns a function handle to an "appender" function, which takes a string
            % as input which is appended to the listbox.
            
            %Can we display progress in a (docked) window within View Manager ?
            % Only applicable if we are a "proper" application object, linked to View Manager
            if isa(obj, 'mvc.mixin.Serializable')
                mgr = obj.viewManager;
            else
                mgr = [];
            end
            if ~isempty(mgr)
                
                %Create the progress view TODO: handle multiple return gracefully
                vw = mgr(1).addView('mvc.view.Progress', [], 'allow');
                
                %And the update function is
                fcn = @vw.update;
                
            else
            
                %Start from simple dialog
                hf = dialog(obj, 'WindowStyle', 'modal', 'Name', tit);

                %Caller may want a cleanup object back
                if nargout > 1
                    clu = onCleanup(@()delete(hf(ishghandle(hf))));
                end
                
                %Containing a vbox
                hv = uiextras.VBox('Parent', hf, 'Padding', 6, 'Spacing', 6);
                
                %With a listbox to display progress messages
                hProgress = uicontrol('Parent', hv, 'Style', 'listbox', 'String', {});
                
                %And a "Cancel" button
                bCancel = false;
                hb = uiextras.HButtonBox('Parent', hv);
                hv.Heights(end) = 25;
                uicontrol('Parent', hb, ...
                    'Style', 'pushbutton', ...
                    'String', 'Cancel', ...
                    'Callback', @i_cancel);
                
                %Passed back to caller
                fcn = @i_update;
                
                %Show us
                figure(hf);
                
            end
            
            %Where the cancel function is
            function i_cancel(varargin)
                
                %Simply
                bCancel = true;
                
            end
            
            %Define a simple update function
            function [bOK, strs] = i_update(str)
                
                %Return flag indicating existence of listbox, and state of Cancel button
                bOK = ishghandle(hProgress) && ~bCancel;
                
                %Not ideal, but need a way to allow caller to explicitly close
                bClose = false;
                
                %Anything to add ?
                if ishghandle(hProgress) && nargin > 0 && ~isempty(str)

                    %Allow explicit close
                    if strcmpi(str, 'close')
                        
                        %Make a note for later
                        bClose = true;
                        
                    elseif isstrprop(str(1), 'punct')
                        
                        %A little presumptious, perhaps, but if the string starts with a punctuation character
                        % append existing last-line of listbox
                        %   - EDIT, C.Szczyglowski 27/11/2019
                        %       Need to call 'sprintf' to actually display
                        %       the punctuation!!
                        hProgress.String{end} = [hProgress.String{end}, sprintf(str)];
                        
                    else

                        %Add new line item to list
                        hProgress.String{end+1} = str;
                        hProgress.Value = numel(hProgress.String);
                        
                    end
                    
                    %Force an update
                    drawnow;
                    
                end
                
                %Caller might want the entire message back too
                if nargout > 1
                    
                    %Safety net
                    if ishghandle(hProgress)
                        strs = hProgress.String;
                    else
                        strs = {};
                    end
                
                    %And fair to assume that once retrieved, we're done with the dialog (??)
                    % delete(hf(ishghandle(hf))); NO
                    
                end
                
                %But need some way of allowing caller to close - this is not ideal
                if bClose
                    delete(hf(ishghandle(hf)));
                end
                
            end   

        end
        
    end
    
    methods (Sealed) % Helper functions
        
        function str = err2str(obj, err, varargin)
            
            %Build a string based on exception object
            str = [err.identifier, ': ', err.message];
            
            %IF the object implements debuggable
            if isa(obj, 'mvc.mixin.Debugable')
                
                %Add a space
                str(end+1) = 10;
                
                %For each element of the stack
                for i = 1:numel(err.stack)
            
                    %Extend string
                    str = [str, 10, ' line ', num2str(err.stack(i).line), ' of function ', err.stack(i).name, ' in file ', err.stack(i).file];
                
                end

            end
            
        end
        
    end
    
    methods (Static) %, Access = protected)
        
        function clu = pointer(hc, ptr)
            
            %If handle not specified
            if nargin < 1 || isempty(hc)
                
                %Go with this
                hc = mvc.mixin.UiTools.currentObject;
                
            end
            
            %If pointer not specified
            if nargin < 2
                
                %Go with watch
                ptr = 'watch';
                
            end
            
            %Get containing figure
            if isempty(hc)
                
                %Nothing to do
                hf = [];
                
            elseif ishghandle(hc)
                
                %Easy
                hf = ancestor(hc, 'Figure');
                
            elseif strncmpi('javahandle_withcallbacks', class(hc), numel('javahandle_withcallbacks')) % isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu') || isa(hc, 'javahandle_withcallbacks.javax.swing.JMenuItem')
                
                %Not quite so easy
                hf = ancestor(mvc.mixin.UiTools.currentObject, 'Figure');
                
            elseif isa(hc, 'mvc.view.internal.dndcontrol')
                
                %Not quite so easy
                hf = ancestor(mvc.mixin.UiTools.currentObject, 'Figure');
            
            elseif isa(hc, 'mvc.model.Application')
                
                %No action required
                hf = [];
                
            else
                warning(['unhandled class ''', class(hc), '''']);
                hf = [];
            end
            
            %Anything ?
            if isempty(hf)
                
                %No
                clu = [];
                return;
                
            elseif nargout > 0
            
                %Return clean-up object restoring old value
                oldptr = get(hf, 'Pointer');
                clu = onCleanup(@()set(hf(isvalid(hf)), 'Pointer', oldptr));
                
            end
            
            %Apply value and force a refresh
            set(hf, 'Pointer', ptr);
            drawnow;
            
        end

        function co = currentObject(co)
        
            %There is one and only one
            persistent CO;
            
            %Ensure remains valid
            if ~isempty(CO) && ~isvalid(CO)
                CO = [];
            end
            
            %May be set
            if nargin > 0
                CO = co;
            end
            
            %May be got
            if nargout > 0
                co = CO;
            end
            
        end
        
        function hf = currentFigure(hc)
            
            %Anything ?
            if isempty(hc)
                
                %No
                hf = [];
                
            elseif ishghandle(hc)
                
                %Get directly from object
                hf = ancestor(hc, 'Figure');
                
            elseif isa(hc, 'javahandle_withcallbacks.javax.swing.JMenuItem')
                
                %Get from helper
                hf = ancestor(mvc.mixin.UiTools.currentObject, 'Figure');
                
            elseif isa(hc, 'javahandle_withcallbacks.MLDropTarget')
                
                %Get from helper
                hf = ancestor(mvc.mixin.UiTools.currentObject, 'Figure');
                
            else
                warning(['class ''', class(hc), ''' not (yet) handled']);
                hf = [];
            end
            
        end
        
        function val = bool2offon(b, varargin)
            %
            % Yes I know you could use 'matlab.lang.OnOffSwitchState', but
            % that is not backward compat with 2015b, and it returns a
            % thing that still needs to be cast to char before being
            % assigned to a uicontrol
            
            %Choice of strings to return
            if nargin < 2
                str = {'off', 'on'};
            elseif nargin == 2 && iscell(varargin)
                str = varargin{1};
            else
                str = varargin;
            end
            
            %Allow caller to force return as a cell array, even if scalar passed in
            bCell = iscell(b);
            if bCell
                b = [b{:}];
            end
            
            %Allow for vectorisation
            if isscalar(b) && ~bCell
                
                %Extract value
                val = str{double(b) + 1};
                
            else
                
                %Extract values
                val = str(double(b) + 1);
                
            end
            
        end
        
        function fn = findfile(dn, mask)
            %
            %Starting from directory DN, looks for file matching MASK.
            % If not found in DN, recursively look in subfolders.
            
            %Placeholder for return value
            fn = [];
            
            %Look for match in this folder
            d = dir(fullfile(dn,mask));
            
            %Anything ?
            if ~isempty(d)
                
                %Yes - reconstruct full filename
                fn = fullfile(d(1).folder, d(1).name);
                
                %We're done
                return;
                
            else
                
                %Any subfolders availabile ?
                d = dir(dn);
                d(~[d.isdir]) = [];
                
                %Exclude any folders associated with OS, or SVN
                d(ismember({d.name},{'.', '..', '.svn'})) = [];
                
                %FOr each in turn
                for j = 1:numel(d)
                    
                    %Recurse
                    fn = mvc.mixin.UiTools.findfile(fullfile(dn, d(j).name), mask);
                    
                    %Anything ?
                    if ~isempty(fn)
                        
                        %Yes - we're done
                        return
                        
                    end
                    
                end
                
            end
            
        end

    end
    
end
