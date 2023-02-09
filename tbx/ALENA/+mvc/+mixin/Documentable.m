classdef (ConstructOnLoad) Documentable < matlab.mixin.SetGet
    %
    % Documentable provides functionality associated with the Help menu, including
    %  Help->About to display an 'About' box
    %  Help->User Guide to open a user guide on screen
    
    properties (Transient)
    
        %The About message
        AboutString;
        
        %Location of relevant document(s)
        DocumentFolder;
        
    end
    
    methods % get/set
        
        function val = get.DocumentFolder(obj)
            
            %Get the underlying value
            val = obj.DocumentFolder;
            
            %If not set explicitly
            if isempty(val)
                
                %Better than nothing
                val = fullfile(fileparts(which(class(obj))), 'doc');
                
            end
            
        end
        
        function val = get.AboutString(obj)
            
            %Get the underlying value
            val = obj.AboutString;
            
            %If not set explicitly
            if isempty(val)
                
                %If we are Deployable, and the EntryPointFcn is set
                if isa(obj, 'mvc.mixin.Deployable') && ~isempty(obj.EntryPointFcn)
                    
                    %Get help from entry point
                    val = help(obj.EntryPointFcn);
                    
                else
                    
                    %Better than nothing
                    val = help(class(obj));
        
                end
                
                %Strip closing references to doc centre
                idx = strfind(val, 'Reference page in Doc Center');
                if ~isempty(idx)
                    val(idx:end) = [];
                end
                
            end
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Documentable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, Inf, ... % Specify a super-low priority, so this appears at end
                    'Help>About...', 'helpAbout', ...
                    'Help>Documents...', 'helpDocuments');
                
                %If we are running in MATLAB
                if ~isdeployed
                    
                    %Provide access to the MATLAB Doc Center
                    addContext(obj, Inf, 'Help>Documentation Centre...', 'helpDocCentre');
                
                    %AND if we are Debuggable
                    if isa(obj, 'mvc.mixin.Debugable')
                        addContext(obj, Inf, 'Help>|Classes...', 'codehelp');
                    end
                    
                end
                
            end

            %Hmmm TODO: Really need to add AboutString and DocumentFolder to
            % 'CopyOnNew' list, but ONLY if they are not empty, but can't test
            % for that by retrieving the property, because the getter will
            % intervene and provide the default.
                        
        end
        
    end
    
    methods (Sealed)
        
        function helpAbout(obj)
            
            %Dislay about string in simple dialog
            uiwait(msgbox(obj, obj.AboutString, 'modal'));
            
        end
        
        function helpDocuments(obj, varargin)
            
            %What have we got ?
            fn = dir(obj.DocumentFolder);
            fn([fn.isdir]) = [];
            if isempty(fn)
                
                %Nothing
                error(['no documentation found in folder ''', obj.DocumentFolder, '''']);
                
            elseif numel(fn) > 1
                
                %Ask user which one(s) are of interest
                sel = listdlg(obj, {fn.name}, 'Choose document(s) of interest...');
                
                %Cancelled ?
                if isempty(sel)
                    return;
                end
                
                %Apply selection
                fn = fn(sel);
                
            end
            
            %Open each selected document
            arrayfun(@(x)winopen(fullfile(x.folder, x.name)), fn);
            
        end
    
        function helpDocCentre(obj, varargin)
            
            %Do it
            doc(class(obj), varargin{:});
            
        end
        
    end
    
end
