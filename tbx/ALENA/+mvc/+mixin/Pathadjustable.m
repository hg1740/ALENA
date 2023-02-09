classdef (ConstructOnLoad) Pathadjustable < matlab.mixin.SetGet
    %
    % The Pathadjustable class provides functionality that allows an object
    % to know about, and adjust dynamically, the location of code on MATLAB path
    %
    % (probably only useful in development / debug environment, so perhaps could
    %  be combined with Debuggable ??)
    
    properties (AbortSet, SetObservable)
        CustomCodeFolder;
    end
    
    properties (Dependent)
        CustomCodeFolders = {};
        CustomCodeFolderIndex;
    end
    
    methods % get/set
        
        function set.CustomCodeFolders(~, val)
            
            %Rather tedious but need to avoid losing original path settings
            % even after "clear all"
            customCodeFolders(val);
            
        end
        
        function val = get.CustomCodeFolders(obj) %#ok<MANU>
            
            %Rather tedious but need to avoid losing original path settings
            % even after "clear all"
            val = customCodeFolders;
            
        end
        
        function val = get.CustomCodeFolderIndex(obj)
            
            %Search on list
            [~, val] = ismember(obj.CustomCodeFolder, obj.CustomCodeFolders);

        end
        
        function set.CustomCodeFolder(obj, val)
            
            %Not yet on our list ?
            if ~ismember(val, obj.CustomCodeFolders) %#ok<MCSUP>
                
                %Add it
                obj.CustomCodeFolders{end+1} = val; %#ok<MCSUP>
                
            end
            
            %Remove current folder from path (silently)
            st = warning('query', 'MATLAB:rmpath:DirNotFound');
            warning('off', 'MATLAB:rmpath:DirNotFound');
            rmpath(genpath(obj.CustomCodeFolder));
            warning(st.state, 'MATLAB:rmpath:DirNotFound');
            
            %Make a note
            obj.CustomCodeFolder = val;
           
            %Add new folder to path
            addpath(genpath(obj.CustomCodeFolder));
            
        end
        
        function val = get.CustomCodeFolder(obj)
        
            %Start here
            val = obj.CustomCodeFolder;
            
            %Nothing yet ?
            if isempty(val) && ~isempty(obj.CustomCodeFolders)
                
                %Return first in list of available folders
                val = obj.CustomCodeFolders{1};
                
            end
            
        end
        
    end
    
    methods % construction
    
        function obj = Pathadjustable(varargin)
             
            %Extend property groups
            if isa(obj, 'mvc.mixin.Nameable')
                obj.addPropertyGroup('General', ...
                    'CustomCodeFolder', 'Custom Code Folder', 'Location of Custom code', @obj.CustomCodeFolders);
            end
            
            %Consistency between multiple instances
            if isa(obj, 'mvc.model.Application')
                obj.CopyOnNew(end+1:end+2) = {'CustomCodeFolder'};
            end
            
        end
        
    end
    
end
