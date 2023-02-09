classdef (ConstructOnLoad) Reportable < matlab.mixin.SetGet
    
    properties (Access = protected)
        
        %File and path to template from which report is created
        TemplateFile;
        TemplatePath;
        
        %To fine tune the template selection process
        TemplateSpec = { ...
            '*.dot', 'Word document templates (*.dot)'; ...
            '*.*', 'All files (*.*)'};
        TemplateIndex = [];
        
    end
    
    properties (Dependent)
        
        %For convenience
        TemplateFullFile;
        
    end
    
    methods       

        function val = get.TemplateFile(obj)
            
            %If nothing
            if isempty(obj.TemplateFile)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.TemplateFile;
                
            end
            
        end

        function val = get.TemplatePath(obj)
            
            %If nothing
            if isempty(obj.TemplatePath)
                
                %Send back empty char
                val = '';
                
            else
                
                %Just send back whatever
                val = obj.TemplatePath;
                
            end
            
        end
        
        function val = get.TemplateFullFile(obj)
            
            %Combine and force to char
            val = fullfile(obj.TemplatePath, obj.TemplateFile);
            
        end
        
        function obj = Reportable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, 'Report...', 'report')
            end

        end
        
        function fn = report(obj, fn, varargin)
            
            %Report to where ?
            if nargin < 2 || isempty(fn)
                
                %Ask the user
                [f, p, idx] = uiputfile(obj.ImportSpec, dlgtitle(obj), obj.ImportFullFile);
                
                %Cancelled ?
                if isempty(f) || isequal(f, 0)
                    
                    %Bail out
                    fn = [];
                    sn = [];
                    return;
                    
                end
                
                %Make a note of full filename
                fn = fullfile(p,f);
            
            else
                
                %Needed later
                idx = [];
                
            end
            
            %If Excel format selected
            [~,~,e] = fileparts(fn);
            if any(strcmpi(e, {'.xls', '.xlsx', '.xlsm'}))
                    
                %What have we got ?
                [~, shts] = xlsfinfo(fn);
                
                %Sheet specified by user ?
                if nargin < 2 || strcmp(sn, 'ask')
                    
                    %Ask the user
                    [sel, bOK] = listdlg(obj, shts, 'Select sheet(s) to import...');
                    
                    %Cancelled ?
                    if ~bOK
                        
                        %Bail out
                        fn = [];
                        sn = [];
                        return;
                        
                    end
                    
                    %Make a note
                    sn = shts(sel);
                    
                elseif isequal(sn, 'all')
                    
                    %Make a note
                    sn = shts;
                    
                else
                    
                    %Assert valid
                    assert(all(ismember(sn, shts)), 'sheet(s) not found in selected file');
                    
                end
                
            else
                
                %Sheetname field not relevant
                sn = [];
                
            end
            
            %Make a note
            [obj.ImportPath, n, e] = fileparts(fn);
            obj.ImportFile = [n, e];
            obj.ImportSheet = sn;
            
            %If we have no filemask index yet
            if isempty(idx)
                
                %Can we work it out ?
                idxx = find(cellfun(@(x)~isempty(strfind(x,e)), obj.ImportSpec(:,1)));
                if numel(idxx) == 1
                    
                    %Yes
                    idx = idxx;
                    
                end
                
            end

            %Make a note
            obj.ImportIndex = idx;
                
        end
        
    end
    
end
