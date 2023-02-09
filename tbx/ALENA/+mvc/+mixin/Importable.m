classdef (ConstructOnLoad) Importable < matlab.mixin.SetGet
    
    properties (Access = protected)
        
        %File and path from which data are imported
        ImportFile;
        ImportPath;
        
        %For use when importing from Excel
        ImportSheet;
        
    end
    
    properties (Access = protected, Transient)
        
        %To fine tune the import process
        ImportSpec = { ...
            '*.xml',               'XML files (*.xml)',                     @import_xml; ...
            '*.xls;*.xlsx;*.xlsm', 'Excel workbooks (*.xls,*.xlsx,*.xlsm)', @import_xls; ...
            '*.mat',               'MATLAB data files (*.mat)',             @import_mat};
        ImportIndex = [];
        
        %Whether property names in imported files are treated as case-sensitive
        ImportCaseSensitive = false;        
        
        %For use when importing from Excel
        ImportXLSParser = @mvc.mixin.Importable.xls2obj;
        
        %For use when importing from XML
        ImportXMLParser = @mvc.mixin.Importable.xml2obj;
        ImportXMLDefaultClass = @mvc.model.Thing;
        ImportXMLMapping;
        
        %To track whether the import process is interrupted
        LastImportCompleted = false;
        
        %For managing MRU
        ImportRecentMax = 5;
        
        %How much feedback does user want ?
        ImportShowProgress = true;
        
    end
    
    properties (Dependent)
        
        %For convenience
        ImportFullFile;
        ImportFileExtensions;
        
        %For managing MRU
        ImportRecent;
    
    end
    
    methods % get/set 
        
        function val = get.ImportRecent(obj)
            
            %Only applicable to "proper" application objects (i.e. serializable ones)
            if isa(obj, 'mvc.mixin.Serializable')
                
                %Pass it on
                val = getpref(obj, 'Recent Imports', {});
            
                %Ensure cellstr
                if isempty(val)
                    val = {};
                end
                
            else
                
                %Send back nothing
                val = {};
                
            end
            
        end
        
        function set.ImportRecent(obj, val)
            
            %Only applicable to "proper" application objects (i.e. serializable ones)
            if isa(obj, 'mvc.mixin.Serializable')
                
                %Ensure cellstr
                assert(iscellstr(val), 'must be a list of strings');
                
                %No point have anything appear twice on the list
                val = unique(val, 'stable');
                
                %Truncate, if necessary
                val(obj.ImportRecentMax+1:end) = [];
                
                %Pass it on
                setpref(obj, 'Recent Imports', val);
            
            end
            
        end
        
        function val = get.ImportFileExtensions(obj)
            
            %Start here
            val = obj.ImportSpec(:,1);
            
            %But break up any semicolon-separated lists
            val = cellfun(@(x)strsplit(x, ';'), val, 'UniformOutput', false);
            val = [val{:}];
            
            %And strip any wildcard characters
            val = strrep(val, '*.', '.');
            
        end
        
        function val = get.ImportFile(obj)
            
            %If nothing
            if isempty(obj.ImportFile)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.ImportFile;
                
            end
            
        end

        function val = get.ImportPath(obj)
            
            %If nothing
            if isempty(obj.ImportPath)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.ImportPath;
                
            end
            
        end
        
        function val = get.ImportFullFile(obj)
            
            %If file not set
            if isempty(obj.ImportFile)
                
                %Send back nothing
                val = [];
                
            else
                
                %Combine and force to char
                val = fullfile(obj.ImportPath, obj.ImportFile);
        
            end
            
        end
        
        function set.ImportFullFile(obj, val)
            
            %Split
            [p, n, e] = fileparts(char(val));
            
            %Make a note
            obj.ImportPath = p;
            obj.ImportFile = [n, e];
            
            %If we have anything
            if ~isempty(val)
                
                %Update MRU list with most recent file FIRST
                mru = obj.ImportRecent;
                mru = [{val}, mru];
                obj.ImportRecent = mru;
            
            end
            
        end
        
    end
    
    methods % construction

        function obj = Importable
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                
                %Add option to Import
                addContext(obj, 'File>|Import...', 'import');
                
                %If we are a "proper" application objects (i.e. serializable ones)
                if isa(obj, 'mvc.mixin.Serializable')
                
                    %Add option to Reimport from MRU
                    addContext(obj, 'File>Reimport', []);
                    
                else
                    
                    %Just add option to reimport last file only
                    addContext(obj, 'File>Reimport...', 'reimport');
                    
                end
                
            end

        end
        
    end
    
    methods

        function varargout = reimport(obj)
            
            %Pass it using existing details
            if isempty(obj.ImportSheet)
                [varargout{1:nargout}] = import(obj, obj.ImportFullFile);
            else
                [varargout{1:nargout}] = import(obj, obj.ImportFullFile, obj.ImportSheet);
            end
            
        end
        
        function [bImportable, strWhyNot] = isImportable(obj, fn)
            %
            % Helper function returns true of file FN exists,
            %  and has an extension that matches one defined in the import spec
            
            %Hope for the best
            bImportable = true;
            strWhyNot = [];
            
            %File exists ?
            if exist(fn, 'file') ~= 2
                
                %Simple
                bImportable = false;
                strWhyNot = 'file does not exist';
                
            else
            
                %Compare file extension with file spec
                [~, ~, e] = fileparts(fn);
                idxx = find(cellfun(@(x)~isempty(strfind(x,e)), obj.ImportSpec(:,1)));
                
                %Anything ?
                if isempty(idxx)
                    
                    %No
                    bImportable = false;
                    strWhyNot = 'file extension does not match any listed in import spec';                    
                    
                end
                
            end
            
        end
        
        function [newobj, logfcn] = import(obj, fn, varargin)
            
            %Placeholder, in case of premature exit before we create the dialog box accessed by logfcn
            logfcn = [];
            
            %Assume the worst
            obj.LastImportCompleted = false;
            
            %Do we need to prompt for file ?
            if nargin < 2
                
                %Yes
                bPrompt = true;
                fn = obj.ImportFullFile;
            
            else
                
                %Maybe
                bPrompt = isempty(fn);
                
            end
            
            %If caller specifies file as numeric
            if ~isempty(fn) && isnumeric(fn)
                
                %Assume it to be an index into the Recent Import list
                fn = obj.ImportRecent{fn};
                
            end
            
            %Helper function shared across classes and methods
            [fn, idx] = filepicker(fn, obj.ImportSpec(:,1:2), bPrompt, false, dlgtitle(obj));
            if isempty(fn)
                
                %Bail out
                newobj = [];
                return;
                
            end
            
            %Check that file exists
            if exist(fn, 'file') ~= 2
                
                %Ensure this file is not on recent list
                mru = obj.ImportRecent;
                mru(strcmpi(fn, mru)) = [];
                obj.ImportRecent = mru;
                
                %Throw
                error(['file ''', fn, ''' does not exist']);
                
            end
                        
            %If Excel format selected, cater for sheetname
            [sn, varargin] = sheetpicker(fn, varargin, bPrompt, false, dlgtitle(obj));
            if islogical(sn) && ~sn
                
                %Bail out
                newobj = [];
                return;
                
            end
            
            %Make a note (before actually dong the import - good idea ??)
            obj.ImportFullFile = fn;
            obj.ImportSheet = sn;
            obj.ImportIndex = idx;
              
            %Pull index-specific import function from import spec
            fcn = obj.ImportSpec{idx,3};
            
            %Anything ?
            assert(~isempty(fcn), 'no import function defined for selected file');
            
            %Use progress window during import ?
            if obj.ImportShowProgress
                
                %Something to look at during import
                logfcn = progressdlg(obj, ['Importing ''', obj.ImportFile, ''' from ''', obj.ImportPath, '''...']);
                
                %Do it
                newobj = fcn(obj, logfcn, varargin{:});
                
                %Escape route
                if isempty(newobj)
                    return
                end
                
                %Retain log for audit trail, if possible
                if isa(newobj, 'mvc.mixin.Auditable')
                    
                    %What happened ?
                    [~, str] = logfcn();
                    
                    %Yes
                    newobj.addAuditTrailEntry(str);
                    
                end
                                
                %Make a note
                obj.LastImportCompleted = true;
                newobj.LastImportCompleted = true;
                
                %If caller wants new object back
                if nargout
                    
                    %Take the name of the imported file
                    [~,nam] = fileparts(fn);
                    
                    %Append sheet, if appropriate
                    if ~isempty(sn)
                        nam = [nam, '_', sn];
                    end
                    
                    %Assign the DefaultName, so if user has not explicitly set it otherwise,
                    % name will be displayed based on import file,
                    % but if name has been explicitly set it will take precedence)
                    newobj.DefaultName = nam;
                    
                    %Does caller want log function back ?
                    if nargout < 2
                        
                        %No
                        logfcn('close');
                        
                    end
                    
                    %Caller taking responsibility, so do nothing more here
                    return;
                    
                end
                
                %If THIS has changed since last save
                if obj.ChangedSinceLastSave
                    
                    %Caller may force overwrite
                    if nargin > 1 && bForce
                        
                        %Crack on
                        
                    elseif ~obj.confirm('Continue WITHOUT saving changes to the current model ?')
                        
                        %Bail out
                        return;
                        
                    end
                    
                end
                
                %Assign imported children as children of this object
                obj.Children = newobj.Children;
                
                %Update audit trail, if appropriate
                if isa(obj, 'mvc.mixin.Auditable')
                    addAuditTrailEntry(obj, ['Imported content from ''', obj.ImportFullFile, '''']);
                end
                
            else
            
                %It's this simple
                fcn(obj);
                
            end
            
            %We made it
            obj.LastImportCompleted = true;

        end
        
        function newobj = import_mat(obj, logfcn, varargin)
            
            %If we are serializable
            if isa(obj, 'mvc.mixin.Serializable')
                
                %What's happening ?
                logfcn(['Reading file ''', obj.ImportFullFile, '''']);
                    
                %Careful
                try
                    
                    %Then for generic behaviour, import from mat might as well be considered same as "open"
                    % but with option to suppress mry
                    newobj = open(obj, obj.ImportFullFile, true, '-nomru', varargin{:});
        
                catch err
                    
                    %Make a note, for the record
                    logfcn(['ERROR: ', err.message]);
                    
                end
                
            else
                error('generic import from mat-file not available');
            end
            
        end
        
        function newobj = import_xls(obj, logfcn)
            
            %Pass it on to recursive parser
            newobj = obj.ImportXLSParser(obj.ImportFullFile, obj.ImportSheet, ...
                obj.ImportXMLDefaultClass, obj.ImportCaseSensitive, logfcn);
            
        end
        
        function newobj = import_xml(obj, logfcn)
            
            %Pass it on to recursive parser
            newobj = obj.ImportXMLParser(obj.ImportFullFile, 1, [], ...
                obj.ImportXMLDefaultClass, obj.ImportCaseSensitive, logfcn);
            
        end
                
        function import_xml_using_mapping_table(obj)
            
            %Update audit trail, if appropriate
            if isa(obj, 'mvc.mixin.Auditable')
                addAuditTrailEntry(obj, ['Importing from ''', obj.ImportFullFile, '''']);
            end
            
            %What's happening ?
            waitbar(obj, 0, ['Importing from ''', obj.ImportFullFile, '''...']);
            cuo = onCleanup(@()waitbar(obj, 'close'));
            
            %Open handle to file
            xDoc = xmlread(obj.ImportFullFile);
            
            %For each entry defined in XML mapping table
            for i = 1:size(obj.ImportXMLMapping,1)
                
                %Property name ?
                nam = obj.ImportXMLMapping{i,1};

                %What's happening ?
                waitbar(obj, i ./ size(obj.ImportXMLMapping,1), ...
                    ['Mapping tag ''', obj.ImportXMLMapping{i,2}, ' to property ', nam]);
                
                %What are we looking for ?
                tok = strsplit(obj.ImportXMLMapping{i,2}, '|');
                
                %Look for element specified by token{1}
                el = xDoc.getElementsByTagName(tok{1});
                
                %Anything ?
                if el.getLength == 0
                    
                    %No - could treat as error, or just crack on ?
                    disp([' tag ''', tok{1}, ''' not found']);
                    continue;
                    
                end
                
                %Get corresponding value, corresponding to
                if numel(tok) == 1
                    
                    %Just go with the element
                    val = el.item(0);
                    
                elseif numel(tok) == 2
                    
                    %Specified attribute
                    val = getNamedItem(get(el.item(0), 'Attributes'), tok{2});
                    
                else
                    
                    %Not (yet) supported
                    error('bad specification of XML mapping');
                    
                end
                
                %Get the value as text
                txt = get(val, 'TextContent');
                
                %Assign in object
                if ischar(obj.(nam))
                    obj.(nam) = txt;
                elseif isnumeric(obj.(nam))
                    obj.(nam) = str2num(txt);
                elseif islogical(obj.(nam))
                    obj.(nam) = logical(str2num(txt));
                elseif isdatetime(obj.(nam))
                    obj.(nam) = datetime(txt);
                else
                    error('unhandled mapping');
                end
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function cls = type2class(cls, typ)
            
            %Cater for cls supplied as function handle or string
            if isa(cls, 'function_handle')
                cls = func2str(cls);
            end
            
            %Look for a file in same location as cls, but named 'typ'
            if exist(fullfile(fileparts(which(cls)), typ), 'file') == 2
                
                %Modify the class name accordingly
                cls = strsplit(cls, '.');
                cls{end} = typ;
                cls = strjoin(cls, '.');
             
            elseif ~strcmp(cls, 'mvc.model.Thing')
                
                %Could not find a matching class in application package,
                % look in mvc package instead   
                src2 = 'mvc.model.Thing';
                cls2 = mvc.mixin.Importable.type2class(src2, typ);
                if ~strcmp(src2, func2str(cls2))
                    cls = cls2;
                end
                
            end
            
            %Always return a function handle
            if ~isa(cls, 'function_handle')
                cls = str2func(cls);
            end
            
        end
        
        function [obj, par] = xls2obj(fn, sn, cls, bCaseSensitive, logfcn)
            %xls2obj Parses workbook/worksheet to construct equivalent component OBJ.
            %
            % Also returns PAR, the name of the parent of this object (if found in sheet)
            
            %Allow multiple sheets
            if iscell(sn)
                
                %Pass it on
                [obj, par] = i_sheets(fn, sn, cls, logfcn);
                                
            else
                
                %Pass it on to helper
                [obj, par] = i_sheet(fn, sn, cls, logfcn);
                
            end
            
            function [obj, par] = i_sheets(fn, sn, cls, logfcn)
                
                %Import each in turn
                [obj, par] = cellfun(@(x)i_sheet(fn, x, cls, logfcn), ...
                    sn, 'UniformOutput', false);
                
                %Discard any empties
                b = cellfun(@isempty, obj);
                obj(b) = [];
                par(b) = [];
                
                %Rebuild parent/child relationships (if possible)
                nam = cellfun(@(x)x.Name, obj, 'UniformOutput', false);
                for i = 1:numel(obj)
                    
                    %Anything to do ?
                    if obj{i}.HasParent
                        
                        %No
                        continue;
                        
                    end
                    
                    %No parent named ?
                    if isempty(par{i})
                        
                        %Move on
                        continue;
                        
                    end
                    
                    %Look for match to left
                    pdx = find(strcmp(par{i}, nam(1:i-1)));
                    
                    %No match ?
                    if isempty(pdx)
                        
                        %Move on
                        continue;
                        
                    end
                    
                    %Go with LAST
                    pdx = pdx(end);
                    
                    %What's happening ?
                    logfcn([' adding ', obj{i}.Name, ' as child of ', obj{pdx}.FullName]);
                    
                    %Do it
                    obj{pdx}.add(obj{i});
                    
                end
                
                %Pull out anything that has not yet been parented
                obj(cellfun(@(x)x.HasParent, obj)) = [];
                
                %Expand and return
                obj = vertcat(obj{:});
                
            end
            
            function [obj, par] = i_sheet(fn, sn, cls, logfcn)
                
                %Placeholders
                obj = [];
                par = [];
                
                %Careful
                try
                    
                    %User cancelled ?
                    if ~logfcn()
                        
                        %Make a note
                        logfcn('Import cancelled by user');
                        
                        %And bail out
                        return;
                        
                    end
                    
                    %What's happening ?
                    logfcn(['Reading sheet ''', sn, ''' of file ''', fn, '''']);
                    
                    %Load the sheet
                    [~,~, raw] = xlsread(fn, sn);
                    
                    %Nothing ?
                    if isempty(raw) || all(all(cellfun(@(x)isscalar(x) && isnan(x), raw)))
                        logfcn(' no useful content');
                        return;
                    end
                    
                    %Clean it up
                    [fav, par, raw] = i_clean(raw);
                    
                    %Anything left ?
                    if ~isempty(raw)
                        
                        %TODO
                        logfcn('...includes additional content from Excel not (yet) implemented');
                        
                    end
                    
                    %Content may override default class
                    idx = find(strcmpi('Type', fav(1,:)));
                    if isempty(idx) || isempty(fav{2,idx})
                        
                        %Nothing
                        typ = [];
                        
                    else
                        
                        %Get detail
                        typ = fav{2,idx};
                        
                    end
                    
                    %If we get this far, we can create a component, which would by default be of the class
                    % specified by caller, but use Type to tell us more precisely what class we're looking for,
                    % (assumed to be in the same package as the default class)
                    fcn = mvc.mixin.Importable.type2class(cls, typ);
                    
                    %What's happening ?
                    logfcn([' ', typ]);
                    
                    %Create the object WITHOUT applying properties,  but setting name as per import sheet
                    % (so we've got something, although expecting it to be overwritten when fav applied)
                    obj = fcn('Name', sn);
                    
                    %Now let's see what is settable
                    b = isEditable(obj, fav{1,:});
                    
                    %Anything to apply ?
                    if any(b)
                        
                        %What's happening ?
                        logfcn([': ', strjoin(fav(1,b), ', ')]);
                        
                        %Do it
                        set(obj, fav{:,b});
                        
                    end
                
                catch err
                    
                    %Make a note for the record
                    logfcn(['FAILED with error: ', err.message]);
                    
                end
                    
            end
            
            function [fav, par, raw] = i_clean(raw)
                
                %Columns 1 and 2 comprise field / value pairs
                fav = raw(:,1:2);
                raw(:,1:2) = [];
                
                %Eliminate any NaNs, which correspond to no data read from Excel
                % (how do we know it was not supposed to be a Nan ?)
                fav(cellfun(@(x)isscalar(x) && isnan(x), fav(:,2)),:) = [];
                
                %Return transposed - easier for caller to bulk-apply with 'set'
                fav = fav.';
                
                %Treat 'Parent' as a special case
                b = strcmpi('Parent', fav(1,:));
                if any(b) && ischar(fav{2,b})
                    
                    %Pull it out
                    par = fav{2,b};
                    fav(:,b) = [];
                    
                else
                    
                    %Nothing
                    par = [];
                    
                end
                
                %This is very presumptious, but consistent with the corresponding export function
                b = cellfun(@(x)ischar(x) && x(1) == '[' && x(end) == ']', fav(2,:));
                
                %Try re-casting anything that looks like it was numeric but is now char
                fav(2,b) = cellfun(@str2num, fav(2,b), 'UniformOutput', false);
                
                %Similarly, identify any multi-line text
                b = cellfun(@(x)ischar(x) && size(x,1) > 1, fav(2,:));
                
                %Re-cast as cellstr
                fav(2,b) = cellfun(@(x)strtrim(cellstr(x)), fav(2,b), 'UniformOutput', false);
                
                %Or this type of multi-line text
                b = cellfun(@(x)ischar(x) && any(x == newline), fav(2,:));
                fav(2,b) = cellfun(@(x)strsplit(x, newline), fav(2,b), 'UniformOutput', false);
                
            end
            
        end
        
        function [obj, par] = xml2obj(item, nDepth, nDepthMin, cls, bCaseSensitive, logfcn, parpth)
            %xml2obj Parses XML node ITEM to construct equivalent component OBJ.
            %
            % Also returns PAR, the name of the parent of this object.
            
            %Parent pathname optional
            if nargin < 7
                parpth = [];
            end
            
            %What minimum depth are we interested in ?
            if nargin < 3 || isempty(nDepthMin)
                nDepthMin = 1;
            end
            
            %Placeholders for return arguments
            obj = [];
            par = [];
            
            %We can start from a file
            if ischar(item)
                
                %If we have an absolute file path
                pth = fileparts(item);
                if ~isempty(pth) && pth(1) ~= '.'
                    
                    %Make a note
                    parpth = pth;
                    
                elseif nargin > 6 && ~isempty(parpth) && isdir(parpth)
                    
                    %Does the file exist in the parent-path ?
                    if exist(fullfile(parpth, item), 'file') == 2
                        
                        %Go with this full file name
                        item = fullfile(parpth, item);
                        
                    end
                    
                end
                
                %If we STILL do NOT have a file path
                if isempty(fileparts(item))
                    
                    %Enlist help from 'which' (without losing name of original input, if no match found)
                    witem = which(item);
                    if ~isempty(witem)
                        item = witem;
                    end
                    
                end
                
                %File must exist to make progress - but is assertion too harsh ?
                %assert(exist(item, 'file') == 2, ['file ''', item, ''' not found']);
                if exist(item, 'file') ~= 2
                    
                    %What's happening ?
                    logfcn([repmat(' ', 1, nDepth), 'ERROR: file ''', item, ''' not found']);
                    
                    %And bail out
                    return;
                    
                end
                    
                %What's happening ?
                logfcn([repmat(' ', 1, nDepth), 'Reading file ''', item, '''']);
                
                %Open it
                item = xmlread(item);
                
                %Treat nDepth as min depth of interest
                if nargin > 1
                    nDepthMin = nDepth;
                end
                
                %And by definition we are now at depth 1
                nDepth = 1;
                
            elseif nargin < 2
                nDepth = 1;
            end
            
            %Anything to do ?
            if i_canIgnore(item)
                
                %No
                return;
                
            end
            
            %To help format messages
            prefix = ': ';

            %Get node name, and initialise object type from node name - may be overridden by XML content
            nam = strtrim(char(item.getNodeName));
            typ = nam;
            
            %Placeholders for other fields and values
            fld = {};
            val = {};
            
            %Placeholders for recursive content
            subobj = {};
            subpar = {};
            
            %For each piece of content
            for i = 1:item.getLength
                
                %User cancelled ?
                if ~logfcn()
                    
                    %Make a note
                    logfcn('Import cancelled by user');
                    
                    %And bail out
                    break;
                    
                end
                
                %Get the sub-item
                subitem = item.item(i - 1);
                
                %Anything to do ?
                if i_canIgnore(subitem)
                    
                    %No
                    continue;
                    
                end
                
                %Get name, text content and length, plus type and size
                subnam = strtrim(char(subitem.getNodeName));
                subtxt = char(subitem.getTextContent);
                sublen = subitem.getLength;
                subtyp = char(subitem.getAttribute('type'));
                subsz = str2num(char(subitem.getAttribute('size')));
                
                %Look for "special" tokens
                if strcmp(subnam, 'Type')
                    
                    %Make a note
                    typ = subtxt;
                    
                elseif strcmp(subnam, 'Parent')
                    
                    %Make a note
                    par = subtxt;
                    
                elseif sublen == 1 || ismember(subtyp, {'table', 'cell', 'struct'}) || ...
                        (ismember(subtyp, {'double', 'logical'}) && numel(subsz) > 2)
                    
                    %What's happening ?
                    logfcn([prefix, subnam]);
                    prefix = ', ';
                    
                    %This node has a single value associated with it - make a note
                    fld{end+1} = subnam;
                    val{end+1} = i_item2val(subtxt, subitem, nDepth, nDepthMin, cls, bCaseSensitive, logfcn, parpth);
                    
                else
                    
                    %What's happening ?
                    logfcn([repmat(' ', 1, nDepth), subnam]);
                    prefix = ': ';
                    
                    %Recurse
                    [subobj{end+1}, subpar{end+1}] = mvc.mixin.Importable.xml2obj(subitem, nDepth + 1, nDepthMin, ...
                        cls, bCaseSensitive, logfcn, parpth);
                    
                    %But keep it clean
                    if isempty(subobj{end})
                        subobj(end) = [];
                        subpar(end) = [];
                    end
                    
                end
                
            end
            
            %If we get this far, and we are at the top of the tree
            if nDepth <= nDepthMin
                
                %Just expand whatever we've found so far
                obj = [subobj{:}];
                
                %And we're done
                return;
                
            end
            
            %Careful
            try
                
                %If we get this far, we can create a component, which would by default be of the class
                % specified by caller, but use Type to tell us more precisely what class we're looking for,
                % (assumed to be in the same package as the default class)
                fcn = mvc.mixin.Importable.type2class(cls, typ);
                
                %Create the component
                obj = fcn();
                
                %Properties for objects of this class
                allPrp = allProperties(obj, true);
                
                %Only interested in settable properties
                prp = allPrp;
                prp(~isEditable(obj, prp{:})) = [];
                
                %Look for match in content read from XML
                if bCaseSensitive
                    [b, pdx] = ismember(fld, prp);
                else
                    [b, pdx] = ismember(lower(fld), lower(prp));
                end
                
                %Pull out the valid fields (and values)
                fav = [prp(pdx(b)); val(b)];
                
                %Assign in object - initialising DefaultName with value read from XML tag
                % (but XML content could itself specify a new name, which would take properity over DefaultName)
                set(obj, 'Name', nam, fav{:}); % TODO - SHOULD BE DefaultName, IDEALLY
                
                %Look for fields read from XML that are NOT properties
                % of the object (including non-editable properties)
                b = ~ismember(fld, allPrp);
                if any(b)
                    
                    %If we support addition of dynamic properties
                    if isa(obj, 'mvc.mixin.Dynamicable')
                        
                        %Pull out the remaining fields (and values)
                        fav = [fld(b); val(b)];
                        
                        %Add them as dynamic properties of object
                        obj.addDynamic('-additional', fav{:});
                        
                    end
                    
                end
                
                %Any recursive content to deal with ?
                if ~isempty(subobj)
                    
                    %Flatten list
                    subobj = [subobj{:}];
                    
                    %Ensure empties are empty strings
                    [subpar{cellfun(@isempty, subpar)}] = deal('');
                    
                    %Match parents across components
                    [b, pdx] = ismember(subpar, {subobj.Name});
                    
                    %For everything else
                    for j = find(b)
                        
                        %Add to parent
                        subobj(j) = subobj(pdx(j)).add(subobj(j), '-ignoreLockedFlags');
                        
                    end
                    
                    %TODO - Avoid this shocking and despicable for loop and
                    %allow a collectable object to add a selection of
                    %groupable and not-groupable objects
                    for kk = find(~b)
                        obj.add(subobj(kk), '-ignoreLockedFlags');
                    end
                        
                    %Add everything else to the parent object
%                     obj.add(subobj(~b), '-ignoreLockedFlags');
                    
                end
            
            catch err
            
                %Make a note for the record
                logfcn(['FAILED with error: ', err.message]);
                
                %And rethrow
                rethrow(err);
                
            end
            
            function val = i_item2val(txt, item, nDepth, nDepthMin, cls, bCaseSensitive, fcn, pth)
                %i_item2val Converts the literal text from the xml node
                %into the correct Matlab data type and validates the
                %attributes.
                
                %Allowable XML attributes
                validXMLAttr = {'type', 'size', 'attr', 'delim'};
                
                %Attributes that must be followed by a numeric as part of
                %the 'validateattributes' implementation.
                tokens = {'numel', 'ncols', 'nromws', 'ndims', '>', '>=', '<', '<='};
                
                %Start with name
                val = txt;
                
                %Get attributes
                attrs = item.getAttributes;
                
                nA = attrs.getLength;
                attr = cell(nA, 2);
                
                %Grab the attributes name-value pairs
                for k = 1:nA
                    a = attrs.item(k - 1);
                    attr{k, 1} = char(a.getName);
                    attr{k, 2} = char(a.getTextContent);
                end
                
                % check attribute list
                idx = ismember(attr(:, 1), validXMLAttr);
                if any(~idx)% no valid attributes have been defined 
                    
                    %Where are we ?
                    stack = i_getStack(item, char(item.getTagName));
                    
                    %Tell us what went wrong
                    fcn(sprintf(['ERROR: Unknown attributes specified in the xml ', ...
                        'file at the node ''%s''. The following '    , ...
                        'attributes are not recognised:\n\t%s']    , ...
                        stack, strjoin(attr(~idx, 1), ', ')));
                    
                    %Bail out
                    return;
                    
                end 
                if nA == 1 && strcmp(attr(idx, 1), 'attr')% only attributes defined but not data type
                    
                    %Where are we ?
                    stack = i_getStack(item, char(item.getTagName));
                    
                    %Tell us what's going on
                    fcn(['WARNING: Attribute "type" is not defined for XML ', ...
                        'node ''', stack, '''. Assuming it is of '  , ...
                        'type ''char''.']);
                    
                    %Bail out
                    return;
                    
                end 
                if nA == 0 % no attributes or data type has been provided 
                    % assume the data is of type char - pass the data on as
                    % literal text
                    attr = {'type', 'char'};
                end 
                type = attr{ismember(attr(:, 1), 'type'), 2};
                
                %Handle data types
                switch type
                    
                    case 'char'
                        
                        %Nothing to do
                    
                    case 'cellstr'
                        
                        %Get delimiter
                        delim = char(item.getAttribute('delim'));
                        
                        %Hence get cellstr
                        val = strsplit(val, delim);
                        
                        %Ensure size is correct
                        sz = str2num(char(item.getAttribute('size')));
                        val = reshape(val, sz);
                        
                    case 'string'
                        
                        %Treat as a comma-separated list of strings
                        val = strsplit(val, ','); % comma is the only valid delimiter
                        
                        %COULD treat this as a proper string, comma-separated
                        % val = arrayfun(@(x) string(x), val, 'Unif', false); 
                        % val = [val{:}];
                        %
                        %But doing so causes problems with 2015b, for which string type
                        % does not exist.  So could switch on verLessThan('matlab', '9.3'),
                        % but surely simpler just to cast input flagged as string into cellstr
                        type = 'cell';
                        
                    case {'double', 'logical'}
                        
                        %Check size, so we can allow for n-d arrays
                        sz = str2num(char(item.getAttribute('size')));

                        %What have we got ?
                        if numel(sz) < 3
                            
                            %Easy - convert as required
                            %   - Need to use str2num as it handles
                            %     comma-seperated text and
                            %     space-delimited text as standard.
                            val = str2num(val); %#ok<ST2NM>
                            if strcmp(type, 'logical')
                                val = logical(val);
                            end
                            
                        else
                            
                            %Initialise
                            val = nan(sz);
                            
                            %Retrieve content page by page
                            for k = 1:1:item.getLength
                            
                                %Unpick
                                item_k = item.item(k - 1);
                                
                                %Anything to do ?
                                if i_canIgnore(item_k)
                                    
                                    %No
                                    continue;
                                    
                                end
                                
                                %Get name and text
                                nam_k = strtrim(char(item_k.getNodeName));
                                txt_k = char(item_k.getTextContent);
                                
                                %Recurse
                                val_k = i_item2val(txt_k, item_k, nDepth + 1, nDepthMin, cls, bCaseSensitive, fcn, pth);
                                
                                %Add page to variable
                                val(:,:,str2num(nam_k(2:end))) = val_k;
                                
                            end
                            
                        end
                        
                    case 'datetime'
                        
                        %Convert as required - allowing for a comma-separated list of datetimes
                        val = datetime(strsplit(val, ','));
                        
                    case 'table'
                        
                        %Start here
                        val = table;
                        
                        %Get content variable by variable
                        for k = 1:item.getLength
                            
                            %Unpick
                            item_k = item.item(k - 1);
                            
                            %Anything to do ?
                            if i_canIgnore(item_k)
                                
                                %No
                                continue;
                                
                            end
                            
                            %Get name and text
                            nam_k = strtrim(char(item_k.getNodeName));
                            txt_k = char(item_k.getTextContent);
                            
                            %Recurse
                            val_k = i_item2val(txt_k, item_k, nDepth + 1, nDepthMin, cls, bCaseSensitive, fcn, pth);
                            
                            %Add to table
                            val.(nam_k) = val_k;
                            
                        end
                    
                    case 'struct'
                        
                        %Start here
                        sz = str2num(char(item.getAttribute('size')));
                        if isempty(sz) || isequal(sz, [1, 1])
                            
                            %Scalar struct
                            val = struct;
                            
                            %Get content by field
                            for k = 1:item.getLength
                                
                                %Unpick
                                item_k = item.item(k - 1);
                                
                                %Anything to do ?
                                if i_canIgnore(item_k)
                                    
                                    %No
                                    continue;
                                    
                                end
                                
                                %Get name and text
                                nam_k = strtrim(char(item_k.getNodeName));
                                txt_k = char(item_k.getTextContent);
                                
                                %Recurse
                                val_k = i_item2val(txt_k, item_k, nDepth + 1, nDepthMin, cls, bCaseSensitive, fcn, pth);
                                
                                %Add to struct
                                val.(nam_k) = val_k;
                            
                            end
                            
                        else
                            
                            %Non-scalar struct - store as cell, temporarily
                            val = {};                            

                            %Retrieve content element by element
                            for k = 1:item.getLength
                            
                                %Unpick
                                item_k = item.item(k - 1);
                                
                                %Anything to do ?
                                if i_canIgnore(item_k)
                                    
                                    %No
                                    continue;
                                    
                                end
                                
                                %Get name and text
                                nam_k = strtrim(char(item_k.getNodeName));
                                txt_k = char(item_k.getTextContent);
                                
                                %Recurse
                                val_k = i_item2val(txt_k, item_k, nDepth + 1, nDepthMin, cls, bCaseSensitive, fcn, pth);
                                
                                %Add to store
                                val{str2num(nam_k(2:end))} = val_k;
                                
                            end
                            
                            %And now expand
                            val = reshape([val{:}], sz);
                            
                        end
                        
                    case 'cell'
                        
                        %Start here
                        sz = str2num(char(item.getAttribute('size')));
                        val = cell(sz);
                        
                        %Get content by element
                        for k = 1:item.getLength
                            
                            %Unpick
                            item_k = item.item(k - 1);
                            
                            %Anything to do ?
                            if i_canIgnore(item_k)
                                
                                %No
                                continue;
                                
                            end                        
                            
                            %Get name and text
                            nam_k = strtrim(char(item_k.getNodeName));
                            txt_k = char(item_k.getTextContent);
                            
                            %Recurse
                            val_k = i_item2val(txt_k, item_k, nDepth + 1, nDepthMin, cls, bCaseSensitive, fcn, pth);
                            
                            %Add to cell
                            val{str2num(nam_k(2:end))} = val_k;
                            
                        end
                        
                    case 'import'
                        
                        %Imports data from the file
                        %   - Assumes 'val' is a valid matlab file
                        val = mvc.mixin.Importable.xml2obj(val, nDepth, nDepthMin, cls, bCaseSensitive, fcn, pth);
                        %It is assumed that 'val' will be a valkid object
                        %belonging to the awi class package so there is no
                        %need to do any more checks on the value --> Quit
                        %out to the invoking function
                        return
                        
                    case 'handle'
                        
                        %XML file contains UUIDs to denote handle
                        %references/pointers between objects.
                        %Handle this in the  same way as character data
                        %   -> i.e. do nothing
                    
                    otherwise
                        
                        %Tell us we don't know how to handle it
                        stack = i_getStack(item, char(item.getTagName));
                        fcn(sprintf('WARNING: Unknown/unhandled data type ''%s'' at ''%s''.\n', type, stack));
                        
                end
                
                attr = attr(ismember(attr(:, 1), 'attr'), 2);
                %If data attributes have not been defined then return value
                %assuming that the data is correct
                if isempty(attr)
                    return
                end
                
                %Construct attribute list for use with 'validateattributes'
                attrList = strsplit(attr{:}, ',');
                idx = ismember(attrList, tokens);   % find tokens
                idx = [false, idx(1:end-1)];
                attrList(idx) = cellfun(@(x) str2double(x), attrList(idx), 'Unif', false);
                
                %Remove any blanks spaces
                attrList(~idx) = strtrim(attrList(~idx));
                
                %For backward compatibility with R2015b
                if verLessThan('matlab', '9.3')
                    
                    %Attribute 'scalartext' not supported
                    attrList(cellfun(@(x)ischar(x) && strcmpi(x, 'scalartext'), attrList)) = [];
                    
                end
                
                %If we get this far then both "type" and "attr" has been
                %defined as attributes of the node and we can now parse
                %the node.
                stack = i_getStack(item, char(item.getTagName));
                validateattributes(val, {type}, attrList, 'xml2obj', stack);
            
            end
            
            function stack = i_getStack(item, stack)
                %i_getStack Gets the complete stack for the xml node 'item'
                
                if nargin < 2
                    stack = '';
                end
                
                % get current parent
                p = item.getParentNode;
                
                % call function recursively until we find top level node
                if ~isempty(p.getParentNode)
                    stack = [char(p.getTagName), '.', stack];
                    stack = i_getStack(p, stack);
                end
                
            end
                       
            function b = i_canIgnore(item)
                
                %Ignore if length is zero, or if nodename matches any on ignore list
                b = item.getLength == 0 || any(strcmp(strtrim(char(item.getNodeName)), {'#text', '#comment'}));
                
            end
            
        end
        
    end            

end
