classdef Rbe0 < hgsetget
    %RBE0 Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id   = [];
        grid = [];
        data = [];
        part = [];
    end
    
    properties(Hidden = true)
        outline = [];
        le      = [];
        te      = [];
    end
    
    methods
        
        function obj = Rbe0(~)
            
        end
        
        %% Set field values of object
        % Define a set method for this class to make sure relevant
        % fields are set and whole object not redefined in a similar
        % way to structure behaviour
        function obj = setFields(obj,Value)
            
            if isa(Value,'Rbe0')
                obj = Value;
            elseif isobject(Value) || isstruct(Value)
                fieldsStruct = fieldnames(Value);
                fieldsObject = fieldnames(obj);
                fieldsObjectReduced = regexprep(fieldsObject,'Obj','');
                
                [~,iVal,iObj] = intersect(fieldsStruct,fieldsObjectReduced);
                
                for i = 1:numel(iObj)
                    obj.(fieldsObject{iObj(i)}) = Value.(fieldsStruct{iVal(i)});
                end
            end
            
        end
                
        function obj = set.data(obj,Value)
            obj.data    = Value;
        end
        
        
        function obj = set.grid(obj,Value)
            obj.grid    = Value;
        end
        
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        function f = genNeocassDeck(obj)
            id   = sprintf('%-8d', obj.id);
            grid = sprintf('%-8d', obj.grid);
            n1   = sprintf('%-8d', obj.data(1));
            n2   = sprintf('%-8d', obj.data(2));
            f{1} = sprintf('RBE0    %s%s%s%s', id, grid, n1, n2);
        end
        
    end
end
