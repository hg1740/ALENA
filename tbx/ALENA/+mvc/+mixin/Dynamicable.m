classdef (ConstructOnLoad) Dynamicable < dynamicprops & matlab.mixin.Copyable
    %
    % Implements dynamicprops, providing a generic means by which, for example,
    % an XML file import routine can have somewhere to put "extra" properties
    % read from file, that should not be thrown away, yet we don't know
    % exactly where to put them in the object being imported.
    
    properties
        
        %To allow distinction between dynamic properties that we know in advance
        % are required by an object, and those added for convenience (e.g. during import)
        AdditionalPropertyNames = {};
        
    end
    
    properties (Dependent)
        
        DynamicProperties;
        DynamicPropertyNames
        ObservableDynamicPropertyNames
        
        AdditionalProperties;
        
    end
    
    properties (Access = public, Transient)
        
        hDynamicProperties;
        
    end
    
    methods % get/set
        
        function val = get.ObservableDynamicPropertyNames(obj)
            
            %Anything ?
            if isempty(obj.hDynamicProperties)
                
                %No
                val = {};
            else
                
                idx = [obj.hDynamicProperties.SetObservable];
                val = {obj.hDynamicProperties(idx).Name};
                
            end
            
            
        end
        
        function set.DynamicProperties(obj, varargin)
            
            %Caller can pass in a single table of properties
            if numel(varargin) == 1 && istable(varargin{1})
                
                %Expand
                varargin = table2cell(varargin{1}).';
               
            elseif numel(varargin) == 1 && iscell(varargin{1})
                
                %Expand
                varargin = varargin{1};
                
            end
            
            %And pass it on
            obj.addDynamic(varargin{:});
            
        end
        
        function val = get.AdditionalProperties(obj)
            %
            % Returns a subset of all dynamic properties,
            %  listing only those dynamic properties that have not been
            %  'claimed' by another property group
        
            %Start here
            val = obj.DynamicProperties;
            
            %Anything else to do ?
            if isempty(val)
                
                %No
                return;
                
            end
            
            %Down-select to only those declared as 'Additiona;'
            val(~ismember(val.Name, obj.AdditionalPropertyNames),:) = [];
            
        end
        
        function val = get.DynamicProperties(obj)
            
            %Anything ?
            if isempty(obj.hDynamicProperties)
                
                %No
                val = {};
                return;
                
            end
            
            %Return a cell-array of [name, value] for all dynamic properties in this instance
            nam = {obj.hDynamicProperties.Name};
            val = [nam; get(obj, nam)].';

            %Or return a table ?
            val = cell2table(val, 'VariableNames', {'Name', 'Value'});
            
        end
        
        function val = get.DynamicPropertyNames(obj)
           
            %Anything ?
            if isempty(obj.hDynamicProperties)
                
                %No
                val = {};
            else
                
                val = {obj.hDynamicProperties.Name};
                
            end
            
        end
        
    end
    
    methods % construction
        
        function obj = Dynamicable
            
            if isa(obj, 'mvc.mixin.Nameable')
                
                %Extend property groups
                obj.addPropertyGroup('Additional', ...
                    'AdditionalProperties', 'Additional Properties');
                
            end
            
        end
        
    end
    
    methods % property management

        function rmDynamic(obj, varargin)
            
            %For each input
            for i = 1:numel(varargin)
                
                %What have we got ?
                if isnumeric(varargin{i})
                    
                    %Treat as index into list of dynamic properties
                    idx = varargin{i};
                    
                elseif ischar(varargin{i})
                    
                    %Look for match with name
                    [b, idx] = ismember(varargin{i}, {obj.hDynamicProperties.Name});
                    assert(all(b), 'no match');
                    
                else
                    error('bad input');
                end
                
                %Remove the dynamic property(ies)
                p = obj.hDynamicProperties(idx);
                obj.hDynamicProperties(idx) = [];
                delete(p);
                
            end
            
        end
        
        function p = addDynamic(obj, varargin)
            
            %Check for the "additional" option
            [bAdd, varargin] = checkOption('-additional', false, varargin{:});
            
            %For pairs of inputs
            for i = 1:2:numel(varargin)
                
                try
                    
                    %Add dynamic property
                    p = addprop(obj, varargin{i});
                    
                    %Make a note
                    if isempty(obj.hDynamicProperties)
                        obj.hDynamicProperties = p;
                    else
                        obj.hDynamicProperties(end+1) = p;
                    end
                    
                catch ME
                    
                    %Check that the error is the one we are expecting
                    switch ME.identifier
                        
                        case 'MATLAB:class:PropertyInUse'
                            
                            %OK - let it go
                            
                        otherwise
                            rethrow(ME);
                            
                    end   
                    
                end
                
                %Assign the supplied value
                obj.(varargin{i}) = varargin{i + 1};
                
            end
                    
            %If flagged as "additional"
            if bAdd
                
                %Make a note
                obj.AdditionalPropertyNames(end+1:end+numel(varargin)/2) = varargin(1:2:end);
                
            end
            
        end
        
    end
    
    methods (Access = protected) % to make dynamic properties copyable
        
        function cpy = copyElement(obj)
        
            %Start with base-class
            cpy = copyElement@matlab.mixin.Copyable(obj);
            
            %Nuke any reference to extant dyanamics
            cpy.hDynamicProperties = [];
            
            %Explicitly copy any dynamic properties in object
            cpy.DynamicProperties = obj.DynamicProperties;
            
        end
        
    end
        
end
