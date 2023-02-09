classdef (ConstructOnLoad) Tabulatable < matlab.mixin.SetGet
    %
    % The Tabulatable class provides functionality that allows an object, or an object array,
    % to draw itself in a number of different ways:
    %   i. In aggregate tabular form, suitable for presenting a summary of object array properties.
    %  ii. In detailed tabular form, suitable for presenting (and editing where appropriate) individual object contents.
    % iii. In detailed graphical form, suitable for presenting individual object content graphically.
    %  iv. In aggregate graphical form, suitable for presenting object array content graphically.
    %
    % The generic Tabulatable class implements a reasonable guess at suitable default behaviour for the above.
    %  Views can be overridden by specific subclasses as required.
    
    properties (Access = protected, Transient, AbortSet)
        
        %How is an array of objects of this class viewed ?
        MultiselectView = 'individual';
        
    end
    
    methods (Sealed)
        
        function [hf, ht, hv, hcc, refresh_content_fcn] = viewer(obj, prp, refresh_name_fcn, varargin)
            %
            % Generic viewer for a (heterogeneous) array of objects derived from Base.
            %
            % Presents object array metadata in tabular form, on LHS of GUI.
            % Presents object-specific view of selected item(s), on RHS of GUI.
             
            %What properties are we interested in ?
            if nargin < 2 || isempty(prp)
                
                %Go with these, by default
                Properties = allProperties(obj);
                
            else
                
                %Just go with whatever caller says
                Properties = prp;
                
            end
           
            %Caller can provide a refresh-name function
            if nargin < 3 || isempty(refresh_name_fcn)
                
                %But did not do so, so use this one by default
                refresh_name_fcn = @()dlgtitle(obj);
                
            end
            
            %Create a suitably sized figure
            [hf, ht, hv, hcc] = i_createFigure(obj);
            
            %Placeholder for current selection
            RowSelection = [];
            ColSelection = [];
            
            %Placeholders for listeners
            ln = {};
            
            %Use helper function to populate content (so easier to trigger a redraw when required)
            i_populatePropertyTable;
                        
            %Send back the refresh function to caller
            refresh_content_fcn = @i_populatePropertyTable;
            
            %--------------------------------------------------------------------------------
            %Where various helper functions are defined thus
            function [hf, ht, hv, hcc] = i_createFigure(obj)
                
                %Start from screensize
                sz = get(0, 'ScreenSize');
                
                %Create figure that fills most of the screen
                hf = figure('NumberTitle', 'off', ...
                    'HandleVisibility', 'callback', ...
                    'Name', refresh_name_fcn(), ...
                    'Position', [100, 100, sz(3:4) - 200]);
             
                %Split figure horizontally
                hh = uiextras.HBoxFlex('Parent', hf, 'Padding', 3, 'Spacing', 3);
                
                %Into two regions, starting with a left-hand panel
                hl = uiextras.Panel('Parent', hh, 'Title', 'Object array');
                
                %And a right-hand pane split vertically, and (initially) left blank
                hv = uiextras.VBoxFlex('Parent', hh, 'Padding', 3, 'Spacing', 3);
                
                %Set good initial relative sizes
                hh.Sizes = [-1, -3];
                
                %Left-hand panel contains a table into which property details are tabulated
                ht = uitable('Parent', hl, 'Units', 'normalized', 'Position', [0, 0, 1, 1], ...
                    'CellSelectionCallback', @i_propertyTableSelectionChange, ...
                    'CellEditCallback', @i_propertyTableCellEdit);

                %Context-menu provides right-click functionality
                hcm = uicontextmenu('Parent', hf);
                set(ht, 'UiContextMenu', hcm);
                
                %Add context menu starting point for selection
                hcc = uimenu(hcm, 'Label', 'Selection');
                
                %Add context menu starting point for the tabular property view
                hcv = uimenu(hcm, 'Label', 'Table');
                
                %For the left-hand table control, user can customise the properties to display
                uimenu(hcv, 'Label', 'Columns...', 'Callback', @i_selectPropertiesToTabulate);
           
                %And permit cell content copy / paste actions
                hm = uimenu(hcv, 'Label', 'Copy to');
                uimenu(hm, 'Label', 'Clipboard', 'Callback', {@i_propertyTableCellCopy, 'clipboard'});
                uimenu(hm, 'Label', 'Above', 'Callback', {@i_propertyTableCellCopy, 'above'});
                uimenu(hm, 'Label', 'All', 'Callback', {@i_propertyTableCellCopy, 'all'});
                uimenu(hm, 'Label', 'Below', 'Callback', {@i_propertyTableCellCopy, 'below'});
                hm = uimenu(hcv, 'Label', 'Paste from');
                uimenu(hm, 'Label', 'Clipboard', 'Callback', {@i_propertyTableCellPaste, 'clipboard'});
                
                %Listen for changes that might require the viewer's name to be updated
                ln{1} = addlisteners(obj, Properties, @i_refreshName);
                
            end
            
            function i_refreshName(varargin)
                
                %Simple
                set(hf, 'Name', refresh_name_fcn());
                
            end
            
            function i_populatePropertyTable(newobj, varargin)
                %
                %Populate the property table with content of object array.
                %
                % This calls 'set' on the uitable, therefore the entire contents are redrawn.
                % This in turn means that any detail of visible rows, selected rows, etc, is lost.
                
                %May be called with a new object array
                % (to support arrays that change size during lifetime of the viewer)
                if nargin > 0
                    
                    %Assign new object array
                    obj = newobj;
                    
                    %If properties being displayed not (yet) set
                    if isempty(Properties)
                        
                        %Go with all available
                        Properties = allProperties(obj);
                        
                    end
                    
                    %Listen for changes that might require the viewer's name to be updated
                    ln{1} = addlisteners(obj, Properties, @i_refreshName);
                    
                end
                
                %Get values of required properties from object array
                val = obj2cell(obj, Properties, 'n/a');
                
                %Convert to tabulatable content
                [val, bInPlaceTabulatable] = obj.tabulatable(val);
                
                %Ensure we only offer in-place edit when it is valid to do so
                bInPlaceEditable = all(isEditable(obj, Properties{:}),1);
                
                %Assign in table
                set(ht, 'Data', val, ...
                    'ColumnName', Properties, ...
                    'ColumnEditable', bInPlaceTabulatable & bInPlaceEditable);
                
                %Listen for changes made elsewhere that impact on content of property table
                ln{2} = addlisteners(obj, Properties, @i_refreshPropertyTable);
                
            end
            
            function i_refreshPropertyTable(varargin)
                %
                %Refresh the property table with content of object array.
                %
                % This uses low-level (java) helper functions, to avoid repopulating entire table.
                % This in turn means that any detail of visible rows, selected rows, etc, is retained.
                % But this is only suitable for use when individual values in existing table are modified,
                % not appropriate if row / columns change.
                
                %Refresh table with these properties from object array
                newval = obj2cell(obj, Properties, 'n/a');
                
                %Convert to tabulatable content
                newval = obj.tabulatable(newval);
                
                %What has NOT changed ?
                val = get(ht, 'Data');
                bSame = cellfun(@isequaln, val, newval);
                
                %Maybe nothing at all
                if ~any(bSame(:))
                    return;
                end
                
                %Enlist some undocumented help
                jScroll = findjobj(ht);
                jTable = jScroll.getViewport.getComponent(0);
                
                %For each change in turn
                [r,c] = find(~bSame);
                for i = 1:numel(r)
                    
                    %Refresh table, without rewriting 'Data' attribute, or changing selection, or anything
                    jTable.setValueAt(newval{r(i),c(i)}, r(i)-1, c(i)-1);
                    
                end
                
            end

            function i_selectPropertiesToTabulate(varargin)
                %
                %Present user with a listbox in which available properties are displayed,
                % the selection initialised with those properties that are actually being
                % tabulated in the left hand panel of the view.
                %
                % User can change the selection, the left-hand properties table redraws itself accordingly.
                
                %Options ?
                allprp = allProperties(obj);
                
                %Current selection
                pdx = find(ismember(allprp, Properties));
                
                %Offer to user
                pdx = listdlg(obj, allprp, 'Properties to view...', ...
                    'InitialValue', pdx);
                
                %Cancelled ?
                if isempty(pdx)
                    
                    %Bail out
                    return;
                    
                end

                %Update store accordingly
                Properties = allprp(pdx);
                
                %Redraw content of property table
                i_populatePropertyTable;
                
            end
            
            function i_propertyTableCellEdit(ht, evtdata)
                
                %Careful
                try
                    
                    %Can be called on a property written to table in abbreviated form (such as cellstr)
                    if ~ht.ColumnEditable(evtdata.Indices(2))
                        
                        %Do nothing
                        return;
                        
                    end
                    
                    %Can be called on a non-settable property (e.g. when table content is updated by other analyses)
                    if ~isEditable(obj, Properties{evtdata.Indices(2)})
                        
                        %Do nothing
                        return;
                        
                    end
                    
                    %Preserve data type going from object to table, and then back again
                    val = obj(evtdata.Indices(1)).(Properties{evtdata.Indices(2)});
                    if ischar(val)
                        
                        %Nothing to worry about
                        newval = evtdata.NewData;
                        
                    elseif isnumeric(val)
                        
                        %Convert from string
                        newval = str2num(evtdata.NewData);
                        
                        %Preserve data type
                        fcn = str2func(class(val));
                        newval = fcn(newval);
                        
                    end
                    
                    %Update underlying object
                    obj(evtdata.Indices(1)).(Properties{evtdata.Indices(2)}) = newval;
        
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj(evtdata.Indices(1)), err, 'modal'));
                    
                    %Revert content
                    ht.Data{evtdata.Indices(1), evtdata.Indices(2)} = evtdata.PreviousData;
                    
                end
                
            end
            
            function i_propertyTableSelectionChange(~, evtdata)

                %Make a note of current column selection
                ColSelection = unique(evtdata.Indices(:,2));
                
                %Make a note of previous and current row selection
                PrevRowSelection = RowSelection;
                RowSelection = unique(evtdata.Indices(:,1));
                
                %No change ?
                if isequal(PrevRowSelection, RowSelection)
                    
                    %Nothing to do
                    return;
                    
                end
                
                %Set context-menu according to current selection
                context(obj(RowSelection), hcc);
                
                %Might take a while
                clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                
                %Careful
                try
                    
                    %Selection view - clean sheet (TODO: may be expensive !)
                    delete(hv.Children);
                    
                    %How many different views required by selection ?
                    [nv, vdx] = views(obj(RowSelection));
                    
                    %For each view required by selection
                    for i = 1:nv
                        
                        %View in separate panel
                        hp = uiextras.Panel('Parent', hv, 'Padding', 3);
                        
                        %Pass it on to viewer
                        view(obj(RowSelection(vdx == i)), hp);
                        
                    end
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj(RowSelection), err, 'modal'));
                    
                end

            end
            
            function i_propertyTableCellCopy(ht, evtdata, dst)
                
                %Careful
                try
                    
                    %Get value(s) associated with current selection
                    val = get(obj(RowSelection), Properties(ColSelection));
                    
                    %Hence property / value pairs
                    pav = [Properties(ColSelection); val];
                    
                    %Copy to where ?
                    switch lower(dst)
                        
                        case 'clipboard'
                            
                        case 'below'
                            
                            %Apply to all after this one
                            set(obj(RowSelection+1:end), pav{:});
                            
                        case 'above'
                            
                            %Apply to all up to this one
                            set(obj(1:RowSelection-1), pav{:});
                            
                        case 'all'
                            
                            %Apply to all
                            set(obj, pav{:});
                            
                        otherwise
                            error('bad copy destination');
                    end
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj(evtdata.Indices(1)), err, 'modal'));
                    
                    %Revert content
                    ht.Data{evtdata.Indices(1), evtdata.Indices(2)} = evtdata.PreviousData;
                    
                end
                
            end
            
            function i_propertyTableCellPaste(ht, evtdata, src)
                
                %Careful
                try
                    
                    %Paste from where ?
                    switch src
                        
                        case 'clipboard'
                            
                        otherwise
                            error('bad paste source');
                    end
                    
                catch err
                    
                    %What went wrong ?
                    uiwait(errordlg(obj(evtdata.Indices(1)), err, 'modal'));
                    
                    %Revert content
                    ht.Data{evtdata.Indices(1), evtdata.Indices(2)} = evtdata.PreviousData;
                    
                end
                
            end
            
        end
        
        function [nv, vdx] = views(obj)
            %
            % Given one or more instances in OBJ, return:
            %  NV - the number of separate views required
            %  VDX - mapping each member of obj to one of the NV views
        
            %Anything to do ?
            if numel(obj) == 1
                
                %No
                nv = 1;
                vdx = 1;
                return;
                
            end
            
            %What are the options ?  TODO: Make this an enumeration
            opt = {obj.MultiselectView};
            bInd = strcmpi('individual', opt);
            bGrp = strcmpi('grouped', opt);
            
            %If everything is viewed individually
            if all(bInd)
                
                %It is easy
                nv = numel(obj);
                vdx = 1:numel(obj);
                
            elseif all(bGrp)
                
                %It is easy-ish, but need to worry about different classes (because the 'view' method is sealed)
                [uv, ~, vdx] = unique(arrayfun(@class, obj, 'UniformOutput', false), 'stable');
                nv = numel(uv);
            
            else
                
                %The individuals are easy
                nv = sum(bInd);
                vdx = find(bInd);
                
                %Add the group(s)
                [uv, ~, vdy] = unique(arrayfun(@class, obj(bGrp), 'UniformOutput', false), 'stable');
                
                %So the totals are
                vdx(end+1:end+numel(vdy)) = vdy + nv;
                nv = nv + numel(uv);
                
            end
            
        end       

    end
    
    methods % viewers
    
        function varargout = view(obj, hp, varargin)
            %
            % Default view method displays content of OBJ in panel HP.
            %
            % Distinguishes between different types of view (modelled closely on matlab.mixin.CustomDisplay):
            %
            %  When called on an empty instance, passes control on to viewEmptyInstance
            %  When called on a scalar instance, passes control on to viewScalarInstance
            %  When called on a non-scalar instance, passes control on to viewNonScalarInstance
            
            %Requires a valid panel
            assert(~(isempty(hp) || ~isa(hp, 'uiextras.Panel') || ~isvalid(hp)), 'invalid panel');
            
            %Start from a clean sheet (Expensive ??)
            delete(hp.Children);

            %What have we got ?
            if numel(obj) == 0
                
                %Pass it on
                [varargout{1:nargout}] = viewEmptyInstance(obj, hp, varargin{:});
                
            elseif numel(obj) == 1
                
                %Pass it on
                [varargout{1:nargout}] = viewScalarInstance(obj, hp, varargin{:});
                
            else
                
                %Pass it on
                [varargout{1:nargout}] = viewNonScalarInstance(obj, hp, varargin{:});
                
            end
            
        end
        
        function ht = viewEmptyInstance(obj, hp, varargin)
        
            %Set panel title
            set(hp, 'Title', dlgtitle(obj));
            
            %Display placeholder text
            ha = axes('Parent', hp, 'Visible', 'off');
            ht = text(0.5, 0.5, 'Make a selection to view content...', ...
                'Parent', ha, 'HorizontalAlignment', 'center');
            
        end
        
        function ht = viewScalarInstance(obj, hp, varargin)
            
            %Pass it on to propedit
            ht = propedit(obj, hp);
            
        end            
        
        function ht = viewNonScalarInstance(obj, hp, varargin)
            
            %Display table in first tab
            ht = uiextras.TabPanel('Parent', hp, 'Padding', 3);
            
            %Get all properties, in tabulatable form
            [val, prp] = obj2cell(obj);
            
            %Down-select to numeric-scalar content only
            bNumericScalar = all(cellfun(@(x)isnumeric(x) && isscalar(x), val),1);
            
            %We can cater for datetime-scalar too
            bDatetimeScalar = all(cellfun(@(x)isdatetime(x) && isscalar(x), val),1);
            
            %Pull out the potentially useful content
            val = val(:,bNumericScalar | bDatetimeScalar);
            prp = prp(:,bNumericScalar | bDatetimeScalar);
            
            %Tabulate data here
            hd = uitable('Parent', ht, ...
                    'ColumnName', prp, ...
                    'RowName', {obj.Name});
            ht.TabTitles{1} = 'Table';

            %But use a helper function to put data into the table itself
            i_tabulate;
            
            %Defined thus
            function i_tabulate(varargin)
                
                %Cater for possibility that listener has triggered this update on an axes that has since gone away
                if isempty(hd) || ~ishandle(hd)
                    
                    %Do nothing TODO: delete the listener that called this in the first place
                    return;
                    
                end
                
                %Refresh content to be tabulated (to allow for this function to be called by a change listener)
                val = obj2cell(obj, prp);
                
                %Need to maintain distinction between the data we hold, and what we show in table
                dat = val;
                dat(cellfun(@isdatetime, dat)) = cellfun(@char, dat(cellfun(@isdatetime, dat)), 'UniformOutput', false);
                
                %Display in table
                set(hd, 'Data', dat);                
                
            end

            %So we can easily update content in event of change
            addlisteners(obj, prp, @i_tabulate);
            
            %Add second tab as container for chart(s)
            hc = uiextras.Panel('Parent', ht);
            ht.TabTitles{2} = 'Chart';
            
            %Display in chart form - need to interject an additional container, else problem with legend later
            hcc = uicontainer('Parent', hc);
            ha = axes('Parent', hcc);
            
            %Placeholder for line handle
            hl = [];
            
            %Showing what against what ?
            idx = 1;
            idy = 2;
            
            %Do it
            i_chart;
            
            %Listen for changes - to anything tabulated
            addlisteners(obj, prp, @i_chart);

            %Allow user interaction
            hcm = uicontextmenu('Parent', ancestor(ha, 'Figure'));
            set(ha, 'UIContextMenu', hcm);
            uimenu(hcm, 'Label', 'X-axes...', 'callback', {@i_axis, 'x'});
            uimenu(hcm, 'Label', 'Y-axes...', 'callback', {@i_axis, 'y'});
            
            function i_axis(~, ~, opt)
                
                %What are we offering ?
                if strcmpi(opt, 'X')
                    
                    %Initialise selection
                    sel = idx;
                    
                    %Selection must be single
                    selmode = 'Single';
                    
                elseif strcmpi(opt, 'Y')
                    
                    %Initialise selection
                    sel = idy;
                    
                    %We can cope with multiple selection
                    selmode = 'Multiple';
                    
                else
                    error('bad option');
                end
                    
                %Offer options to user
                [sel, bOK] = listdlg(obj, prp, [opt, '-axis variable...'], ...
                    'SelectionMode', selmode, 'InitialValue', sel);
                
                %Cancelled ?
                if isempty(sel) || ~bOK
                    
                    %Bail out
                    return;
                    
                end
                
                %Apply selection
                if strcmpi(opt, 'X')
                    
                    %Initialise selection
                    idx = sel;
                    
                elseif strcmpi(opt, 'Y')
                    
                    %Initialise selection
                    idy = sel;
                    
                else
                    error('bad option');
                end
                
                %Refresh chart
                i_chart;
                
            end
            
            function i_chart(varargin)
                
                %Cater for possibility that listener has triggered this update on an axes that has since gone away
                if isempty(ha) || ~ishandle(ha)
                    
                    %Do nothing TODO: delete the listener that called this in the first place
                    return;
                    
                end
                
                %Anything to do ?
                if any(size(val,2) < [idx, idy])
                    
                    %No
                    cla(ha);
                    text(ha, 0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
                    
                else
                    
                    %If we simply use plot, then we'd lose any line formatting
                    if isempty(hl) || ~all(ishandle(hl)) || numel(hl) ~= numel(idy)
                        
                        %Deal arguments out to allow mix of numeric and datetimes
                        args = repmat({vertcat(val{:,idx})}, 1, numel(idy) .* 2);
                        for i = 1:numel(idy)
                            args{i .* 2} = vertcat(val{:,idy(i)});
                        end
                        
                        %Plot all in a one-er
                        hl = plot(ha, args{:});
                        
                    else
                        
                        %Modify existing lines, so we retain format etc
                        for i = 1:numel(idy)
                            set(hl(i), 'XData', vertcat(val{:,idx}), 'YData', vertcat(val{:,idy(i)}));
                        end
                        
                    end
                    
                    %What have we got ?
                    xlabel(ha, prp{idx});
                    ylabel(ha, strjoin(prp(idy), ', '));
                    if numel(idy) > 1
                        legend(ha, prp(idy));
                    else
                        legend(ha, 'off');
                    end
                
                end
                
            end
            
        end
        
    end
    
end
