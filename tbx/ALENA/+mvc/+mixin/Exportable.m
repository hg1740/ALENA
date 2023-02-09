classdef (ConstructOnLoad) Exportable < matlab.mixin.SetGet
    
    properties (Access = protected)
        
        %Which properties of the object are exportable ?
        ExportProp = {};
        
        %File and path from which data are exported
        ExportFile;
        ExportPath;
        
        %For use when exporting to Excel
        ExportSheet = '';
        ExportSheetTag = 'Name';
        
        %For use when exporting to XML
        ExportXmlTag = 'Type';
        
        %To fine tune the export process
        ExportSpec = { ...
            '*.xml',                'XML files (*.xml)',                        @export_xml; ...
            '*.xls;*.xlsx;*.xlsm',	'Excel workbooks (*.xls,*.xlsx,*.xlsm)',    @export_xls; ...
            '*.mat',                'MATLAB data files (*.mat)',                @export_mat};
        ExportIndex = [];
        
        %To track whether the export process is interrupted
        LastExportCompleted = false;
        
    end
    
    properties (Dependent)
        
        %For convenience
        ExportFullFile;
        ExportFileExtensions;
        
    end
    
    methods % get/set 
        
        function val = get.ExportFileExtensions(obj)
            
            %Start here
            val = obj.ExportSpec(:,1);
            
            %But break up any semicolon-separated lists
            val = cellfun(@(x)strsplit(x, ';'), val, 'UniformOutput', false);
            val = [val{:}];
            
            %And strip any wildcard characters
            val = strrep(val, '*.', '.');
            
        end
        
       function val = get.ExportFile(obj)
            
            %If nothing
            if isempty(obj.ExportFile)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.ExportFile;
                
            end
            
        end
        
        function val = get.ExportPath(obj)
            
            %If nothing
            if isempty(obj.ExportPath)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.ExportPath;
                
            end
            
        end
        
        function val = get.ExportFullFile(obj)
            
            %Combine and force to char
            val = fullfile(obj.ExportPath, obj.ExportFile);
            
        end
        
    end
    
    methods % construction

        function obj = Exportable
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, 'File>Export...', 'export')
            end
            
        end
        
    end
    
    methods

        function varargout = export(obj, fn, varargin)
            
            %Assume the worst
            obj.LastExportCompleted = false;
            
            %Do we need to prompt for file ?
            if nargin < 2
                
                %Yes
                bPrompt = true;
                fn = obj.ExportFullFile;
                
            else
                
                %No
                bPrompt = false;
                
            end
            
            %Helper function shared across classes and methods
            [fn, idx] = filepicker(fn, obj.ExportSpec(:,1:2), bPrompt, true, dlgtitle(obj));
            if isempty(fn)
                
                %Bail out
                [varargout{1:nargout}] = deal([]);
                return;
                
            end
            
            %Sheetname inferred from object ?
            if isempty(obj.ExportSheetTag)
                
                %Pass it on to helper
                [sn, varargin] = sheetpicker(fn, varargin, bPrompt, true, dlgtitle(obj));
    
            else
                
                %Tag sheetname with a dollar sign (?? good idea ??) to indicate it is
                % the field from which sheetname is taken, rather than sheetname itself
                sn = ['$', obj.ExportSheetTag];
                
            end
            
            %Make a note (before actually dong the export - good idea ??)
            [obj.ExportPath, n, e] = fileparts(fn);
            obj.ExportFile = [n, e];
            obj.ExportIndex = idx;
            obj.ExportSheet = sn;
            
            %What sort of export ?
            fcn = obj.ExportSpec{idx,3};
            
            %Anything ?
            if isempty(fcn)
                
                %No
                warning('no default export function defined for selected file');
                [varargout{1:nargout}] = deal([]);
                
            else
                
                %Might take a while
                clu = obj.pointer; %#ok<NASGU>
                
                %Do it
                [varargout{1:nargout}] = fcn(obj, varargin{:});
                
            end
            
            %We made it
            obj.LastExportCompleted = true;
            
        end
        
        function export_mat(obj, varargin)
            
            %Make a note of filename before overwriting obj
            fn = obj.ExportFullFile;
            
            %Starting where ?
            obj = export_startingPoint(obj, varargin{:});
            if isempty(obj)
                return;
            end
            
            %If we are serializable
            if isa(obj, 'mvc.mixin.Serializable')
                
                %Then for generic behaviour, export to mat might as well be considered same as "save",
                % but with option to suppress mru
                save(obj, fn, [], '-nomru', varargin{:});
                
            else
                
                %Better than nothing
                save(fn, 'obj');
                
            end
            
        end
        
        function export_xml(obj, varargin)
            
            %Starting where ?
            obj2export = export_startingPoint(obj, varargin{:});
            if isempty(obj2export)
                return;
            end
            
            %Pass it on to static helper function
            obj.obj2xml(obj2export, obj.ExportXmlTag, obj.ExportFullFile);
            
        end

        function export_xls(obj, varargin)
            
            %Starting where ?
            obj2export = export_startingPoint(obj, varargin{:});
            if isempty(obj2export)
                return;
            end
            
            %Pass it on to static helper function
            obj.obj2xls(obj2export, obj.ExportSheet, obj.ExportFullFile);
            
        end
        
        function obj2export = export_startingPoint(obj, varargin)

            %If we are a collection, and caller has specified some selection criteria
            if ~isempty(varargin) && isa(obj, 'mvc.mixin.Collectable')
                    
                %Pass it on
                obj2export = obj.selectFromAll(varargin{:});
                
            else
                
                %There's no choice
                obj2export = obj;
                
            end
            
        end
        
        function idx = sortExportSpec(obj, val, bRetain)
            %
            % Utility function allows "most likely" export file format to be brought to
            % top of list in export spec
            %
            % e.g. obj.sortExportSpec({'xls', 'xml', 'mat'}) will ensure that 'xls'
            %  is the default offering, when export file dialog is presented.
            %
            %  If the "bRetain" flag is specified TRUE, then any content in existing
            %  export spec that is not specified in VAL will be retained, but moved
            %  to end of export spec
        
            %Ensure we have a cell-array of values
            if ischar(val)
                val = {val};
            elseif ~iscellstr(val)
                error('value must be a char or a cellstr');
            end
            
            %What have we got at the moment ?
            ext = obj.ExportSpec(:,1);
            
            %Look for match in each of varargin with export spec file extensions
            idx = cellfun(@(x)min([Inf, find(~cellfun(@isempty, regexp(ext, x)), 1, 'first')]), val);
            
            %Optionally retain anything not found
            if nargin > 2 && bRetain
                idx = [idx, setdiff(1:numel(ext), idx)];
            end
            
            %Adjust order accordingly
            obj.ExportSpec = obj.ExportSpec(idx,:);
            
        end
        
    end
    
    methods (Static)
        
        function obj2xml(obj, fld, doc, par)

            %If document name rather than handle passed in
            if ischar(doc)
                
                %Make a note of the filename
                fn = doc;
                
                %Use the DOM interface to create starting node
                doc = com.mathworks.xml.XMLUtils.createDocument(obj.(fld));
                
                %Get root node
                par = doc.getDocumentElement;
                
                %Ensure we close handles and create document on exit
                clu = onCleanup(@()xmlwrite(fn, doc));

            end
            
            %Properties to export ?
            prp = allProperties(obj, 'Export', true);
            
            %Corresponding values
            val = get(obj, prp);
            
            %For each in turn
            for i = 1:numel(prp)
                
                %Convert to exportable
                [node, attr] = i_node(doc, val{i});
                
                %Nothing ?
                if isempty(node)
                    
                    %Move on
                    continue;
                    
                end
                
                %Create corresponding new element
                elem = doc.createElement(prp{i});
                
                %Set attributes if necessary
                for j = 1:2:numel(attr)
                    elem.setAttribute(attr{j:j+1});
                end
                
                %Append node (or nodes) to new element
                if iscell(node)
                    
                    %Append each in turn
                    cellfun(@(x)elem.appendChild(x), node, 'UniformOutput', false);
                    
                else
                    
                    %Just the one
                    elem.appendChild(node);
    
                end
                
                %Append new element to parent
                par.appendChild(elem);
                
            end
            
            %If we have any children
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %For each in turn
                for i = 1:numel(obj.Children)
                    
                    %Name of entry point determined by passed-in argument
                    nam = obj.Children(i).(fld);
                    
                    %UNLESS we are of type 'Collector' (TODO: properly)
                    if isa(obj.Children(i), 'mvc.model.Collector')
                        
                        %In which case, go with Name instead
                        nam = matlab.lang.makeValidName(obj.Children(i).Name);
                        
                    end
                    
                    %Create entry point
                    elem = doc.createElement(nam);
                    
                    %Recurse
                    mvc.mixin.Exportable.obj2xml(obj.Children(i), fld, doc, elem);
                    
                    %Append
                    par.appendChild(elem);
                    
                end
                
            end
            
            function [node, attr] = i_node(docNode, val)
                
                %What have we got ?
                switch class(val)
                                            
                    case 'table'
                        
                        %Start here
                        node = {};
                        attr = {'type', class(val), 'size', mat2str(size(val))};
                       
                        %For tables, take each column in turn
                        col = val.Properties.VariableNames;
                        
                        %For each column
                        for k = 1:numel(col)
                            
                            %Recurse
                            [node_k, attr_k] = i_node(docNode, val.(col{k}));
                            
                            %Nothing ?
                            if isempty(node_k)
                                
                                %Move on
                                continue;
                                
                            end
                            
                            %Create corresponding element
                            elem_k = docNode.createElement(col{k});
                            
                            %Set attributes if necessary
                            for kj = 1:2:numel(attr_k)
                                elem_k.setAttribute(attr_k{kj:kj+1});
                            end
                            
                            %Append node (or nodes) as appropriate
                            if iscell(node_k)
                                
                                %Append each in turn
                                cellfun(@(x)elem_k.appendChild(x), node_k, 'UniformOutput', false);
                                
                            else
                                
                                %Just the one
                                elem_k.appendChild(node_k);
                                
                            end
                            
                            %Add to list
                            node{end+1} = elem_k;
                            
                        end
                        
                    case 'struct'
                        
                        %Start here
                        node = {};
                        attr = {'type', class(val), 'size', mat2str(size(val))};
                       
                        %For structs, take each field in turn
                        col = fieldnames(val);
                        
                        %If it is a scalar structure
                        if numel(val) == 1
                            
                            %For each field as appropriate
                            for k = 1:numel(col)
                                
                                %Recurse
                                [node_k, attr_k] = i_node(docNode, val.(col{k}));
                                
                                %Nothing ?
                                if isempty(node_k)
                                    
                                    %Move on
                                    continue;
                                    
                                end
                                
                                %Create corresponding element
                                elem_k = docNode.createElement(col{k});
                                
                                %Set attributes if necessary
                                for kj = 1:2:numel(attr_k)
                                    elem_k.setAttribute(attr_k{kj:kj+1});
                                end
                                
                                %Append node (or nodes) as appropriate
                                if iscell(node_k)
                                    
                                    %Append each in turn
                                    cellfun(@(x)elem_k.appendChild(x), node_k, 'UniformOutput', false);
                                    
                                else
                                    
                                    %Just the one
                                    elem_k.appendChild(node_k);
                                    
                                end
                                
                                %Add to list
                                node{end+1} = elem_k;
                                
                            end
                        
                        else
                            
                            %For structure arrays, need an extra level of hierarchy,
                            % based on element number (rather like for cell array export)
                            for k = 1:numel(val)
                                
                                %Recurse
                                [node_k, attr_k] = i_node(docNode, val(k));
                                
                                %Nothing ?
                                if isempty(node_k)
                                    
                                    %Move on
                                    continue;
                                    
                                end
                                
                                %Create corresponding element
                                elem_k = docNode.createElement(['_', num2str(k)]);
                                
                                %Set attributes if necessary
                                for kj = 1:2:numel(attr_k)
                                    elem_k.setAttribute(attr_k{kj:kj+1});
                                end
                                
                                %Append node (or nodes) as appropriate
                                if iscell(node_k)
                                    
                                    %Append each in turn
                                    cellfun(@(x)elem_k.appendChild(x), node_k, 'UniformOutput', false);
                                    
                                else
                                    
                                    %Just the one
                                    elem_k.appendChild(node_k);
                                    
                                end
                                
                                %Add to list
                                node{end+1} = elem_k;
                                
                            end
                            
                        end
                        
                    case 'char'
                        
                        %Easy
                        node = docNode.createTextNode(val);
                        attr = {};
                    
                    case 'cell'
                        
                        %If a cellstr
                        if iscellstr(val)
                            
                            %Choose a "safe" delimiter
                            delim = ',;/_';
                            b = cellfun(@(x)ismember(delim, x), val, 'UniformOutput', false);
                            b = any(vertcat(b{:}),1);
                            delim(b) = [];
                            assert(~isempty(delim), 'unable to identify a "safe" delimiter');
                            delim = delim(1);
                            
                            %Easy
                            node = docNode.createTextNode(strjoin(val, delim));
                            attr = {'type', 'cellstr', 'size', mat2str(size(val)), 'delim', delim};
                            
                        else
                            
                            %Start here
                            node = {};
                            attr = {'type', 'cell', 'size', mat2str(size(val))};
                            
                            %Take each element in turn
                            for k = 1:numel(val)
                                
                                %Recurse
                                [node_k, attr_k] = i_node(docNode, val{k});
                                
                                %Nothing ?
                                if isempty(node_k)
                                    
                                    %Move on
                                    continue;
                                    
                                end
                                
                                %Create corresponding element
                                elem_k = docNode.createElement(['_', num2str(k)]);
                                
                                %Set attributes if necessary
                                for kj = 1:2:numel(attr_k)
                                    elem_k.setAttribute(attr_k{kj:kj+1});
                                end
                                
                                %Append
                                elem_k.appendChild(node_k);
                                
                                %Add to list
                                node{end+1} = elem_k;
                                
                            end
                            
                        end
                        
                    case {'double', 'logical'}
                        
                        
                        %Make a note of attributes as required
                        node = {};
                        sz = size(val);
                        attr = {'type', class(val)};
                        if ~isscalar(val)
                            attr(end+1:end+2) = {'size', mat2str(sz)};
                        end
                        
                        %Cater for n-d arrays
                        if numel(sz) > 2
                            
                            %Write value page by page
                            for k = 1:prod(sz(3:end))
                                
                                %Recurse
                                [node_k, attr_k] = i_node(docNode, val(:,:,k));
                                
                                %Nothing ?
                                if isempty(node_k)
                                    
                                    %Move on
                                    continue;
                                    
                                end
                                
                                %Create corresponding element
                                elem_k = docNode.createElement(['_', num2str(k)]);
                                
                                %Set attributes if necessary
                                for kj = 1:2:numel(attr_k)
                                    elem_k.setAttribute(attr_k{kj:kj+1});
                                end
                                
                                %Append node (or nodes) as appropriate
                                if iscell(node_k)
                                    
                                    %Append each in turn
                                    cellfun(@(x)elem_k.appendChild(x), node_k, 'UniformOutput', false);
                                    
                                else
                                    
                                    %Just the one
                                    elem_k.appendChild(node_k);
                                    
                                end
                                
                                %Add to list
                                node{end+1} = elem_k;  
                                
                            end
                            
                        else
                            
                            %Create node
                            node = docNode.createTextNode(mat2str(val));
                            
                        end                       
                        
                    case 'datetime'
                        
                        %Create node
                        node = docNode.createTextNode(strjoin(cellstr(val), ','));
                        
                        %Make a note of attributes as required
                        attr = {'type', 'datetime'};
                        if ~isscalar(val)
                            attr(end+1:end+2) = {'size', mat2str(size(val))};
                        end
                        
                    otherwise
                        
                        %IF the item is a member of a collection
                        if isa(val, 'mvc.mixin.Collectable')
                            
                            %Then export its full name, rather than the item itself
                            % (which leaves the problem of rebuilding the network properly
                            %  to the import routine, as and when this file is re-imported)
                            node = docNode.createTextNode(val.FullName);
                            attr = {'type', 'node'};
                            
                        else
                            
                            %Tell us about it (useful for debug/devel)
                            disp(['exporting class ''', class(val), ''' to xml not yet implemented']);
                            
                            %Not (yet) handled
                            node = [];
                            attr = {};
                        
                        end
                        
                end
                
            end
            
        end
        
        function sns = obj2xls(obj, sn, fn, sns)
            
            %Ensure we don't overwrite any extant sheets
            if nargin < 4
                if exist(fn, 'file') == 2
                    [~, sns] = xlsfinfo(fn);
                else
                    sns = {};
                end
            end
            
            %Properties to export ?
            prp = allProperties(obj);
            
            %Ensure we have a row
            prp = prp(:).';
            
            %Need to ensure that Name and Type are ALWAYS included in export detail (to help with subsequent re-import)
            fld = {'Name', 'Type'};
            b = ~ismember(fld, prp);
            prp(end+1:end+sum(b)) = fld(b);
            
            %Corresponding values
            val = get(obj, prp);
            
            %If we are part of a collection, and parent is defined
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasParent
                
                %Look for property detail
                [b, pdx] = ismember('Parent', prp);
                if ~b
                    
                    %Add it explicitly
                    prp{end+1} = 'Parent';
                    pdx = numel(prp);
                    
                end
                
                %Include the name of the parent, not the object itself (obviously)
                val{pdx} = obj.Parent.Name;

            end
            
            %Which of these are easy exports straight to Excel ?
            bChar = cellfun(@ischar, val);
            bNumericScalarOrEmpty = cellfun(@(x)isnumeric(x) && (numel(x) <= 1), val);
            bEasy = bChar | bNumericScalarOrEmpty;
            
            %Cellstr are straight forward too
            bCellstr = cellfun(@iscellstr, val);
            val(bCellstr) = cellfun(@(x)strjoin(x, newline), val(bCellstr), 'UniformOutput', false);
            bEasy = bEasy | bCellstr;
            
            %Rows and columns are reasonably straight forward
            bRowCol = cellfun(@(x)isnumeric(x) && (sum(size(x) == 1) == numel(size(x)) - 1), val);
            
            %But possible loss of precision here
            val(bRowCol) = cellfun(@mat2str, val(bRowCol), 'UniformOutput', false);
            bEasy = bEasy | bRowCol;
            
            %Pull the easy ones out
            dat = [prp(bEasy); val(bEasy)].';
            
            %There might be nothing to write
            if ~isempty(dat)
                
                %Sheetname might be specified by reference to a field
                if sn(1) == '$' && isProperty(obj, sn(2:end)) && ischar(obj.(sn(2:end)))
                    
                    %Use it
                    nam = obj.(sn(2:end));
                    
                else
                    
                    %Just go with whatever
                    nam = sn;
                    
                end
                    
                %Make sure this sheetname is new
                nam = matlab.lang.makeUniqueStrings(nam, sns);
                
                %Send it to Excel
                mvc.mixin.Exportable.xlswrite(fn, dat, nam);
        
                %Add sheetname to list
                sns{end+1} = nam;
                
            end
            
            %If we have any children
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %For each in turn
                for i = 1:numel(obj.Children)
                                        
                    %Recurse (recycling list of sneetnames to ensure all child sheets are uniquely named)
                    sns = mvc.mixin.Exportable.obj2xls(obj.Children(i), sn, fn, sns);
                    
                end
                
            end
                        
        end
        
        function [b, msg] = csvwrite(fn, dat, delim)
        
            %Hope for the best
            b = 1;
            msg = [];
            
            %Using delimeter
            if nargin < 3 || isempty(delim)
                delim = ',';
            end

            %Careful
            try
                
                %Open file for writing
                fid = fopen(fn, 'w');
                assert(fid ~= -1, ['unable to open file ''', fn, ''' for writing']);
                
                %Ensure we can't get out of here without closing file
                clu = onCleanup(@()fclose(fid));
                
                %Write content
                sz = size(dat);
                for r = 1:sz(1)
                    for c = 1:sz(2)
                        
                        %What have we got ?
                        if ischar(dat{r,c})
                            
                            %Write as string
                            fprintf(fid, '%s', dat{r,c});
                            
                        elseif isnumeric(dat{r,c})
                            
                            %Write as numeric
                            fprintf(fid, '%d', dat{r,c});
                            
                        else
                            
                            %Let us know we couldn't handle it
                            fprintf(fid, '%s', ['-', class(dat{r,c}), '-']);
                            
                        end
                        
                        %Delimeter
                        if c < sz(2)
                            fprintf(fid, delim);
                        else
                            fprintf(fid, newline);
                        end
                        
                    end
                    
                end
        
            catch err
                
                %Lower ok flag
                b = 0;
                
                %What went wrong ?
                msg = err.message;
                
            end
            
        end
        
        function varargout = xlswrite(varargin)
            %
            %Simple wrapper around built-in function, but suppresses warnings of type "AddSheet",
            % which are otherwise rather annoying
            
            %Make a note of current warning state
            st = warning('query', 'MATLAB:xlswrite:AddSheet');
            
            %Turn warnings off
            warning('off', 'MATLAB:xlswrite:AddSheet');
            
            %Now do whatever
            [varargout{1:nargout}] = xlswrite(varargin{:});
            
            %Reinstate warnings
            warning(st.state, 'MATLAB:xlswrite:AddSheet');

        end
        
    end
    
end
