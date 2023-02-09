classdef DefaultBuild < matlab.mixin.SetGet
    %DefaultBuilds Handles a collection of static methods which allow
    %various default objects to be created.
    %
    % Detailed Description:
    %   - Each static method must return an instance of the object
    
    properties (SetAccess = protected)
        %Collection of default build methods.
        DefaultObjects = struct('Name', {}, 'Method', {}, 'Description', {});
    end
    
    properties (Dependent)
        %List of available defaults
        AvailableDefaults
    end
    
    methods % set / get
        function val = get.AvailableDefaults(obj)
            val = {obj.DefaultObjects.Name};
        end
    end
    
    methods % handling defaults
        function NewObj = runDefaultBuild(obj, name)
            %runDefaultBuild Executes a default build method and returns
            %the new object.
                        
            assert(ischar(name), ['Expected the build name to be a ', ...
                'row vector of characters']);
            idx = ismember(obj.AvailableDefaults, name);
            assert(nnz(idx) == 1, ['Ambiguous match for default build ', ...
                'name ''%s''. Make sure the name is one of the following ', ...
                'options:\n\t%s\n'], name, obj.AvailableDefaults);
            
            func = obj.DefaultObjects(idx).Method;
            
            NewObj = func();
            
        end
        function addDefaultBuild(obj, name, func_handle, varargin)
            %addDefaultBuild Adds a new default build option to the
            %collection.
            
            p = inputParser;
            addRequired(p, 'name'       , @isvarname);
            addRequired(p, 'func_handle', @(x) isa(x, 'function_handle'));
            addOptional(p, 'descr', ''  , @(x)validateattributes(x, {'char'}, {'row'}));
            parse(p, name, func_handle, varargin{:});
           
            default_methods = cellfun(@char, {obj.DefaultObjects.Method}, 'Unif', false);
                        
            assert(~any(ismember(obj.AvailableDefaults, name)), ['There ', ...
                'is already a default build with name ''%s''. Choose ', ...
                'another name and try again.'], name);
            assert(~any(ismember(default_methods, char(func_handle))), ...
                ['There is already a default build that uses the method ', ...
                '''%s''. Choose another method and try again.'], char(func_handle));
            
            S = struct( ...
                'Name'       , name, ...
                'Method'     , func_handle, ...
                'Description', p.Results.descr);
            
            obj.DefaultObjects = [obj.DefaultObjects, S];
            
        end
    end
    
end

