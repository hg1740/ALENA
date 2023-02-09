classdef Spline < hgsetget
    %SPLINE Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id   = [];
        aid  = [];
        p1   = [];
        p2   = [];
        w    = [];
        rmx  = [];
        cond = [];
        set  = [];
        part  = [];
    end
    
    methods
        
        function obj = Spline(~)
            
        end
        
        function obj = setFields(obj,Value)
            
            if isa(Value,'Spline')
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
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.aid(obj,Value)
            obj.aid     = Value;
        end
        
        
        function obj = set.p1(obj,Value)
            obj.p1     = Value;
        end
        
        
        function obj = set.p2(obj,Value)
            obj.p2     = Value;
        end
        
        
        function obj = set.w(obj,Value)
            obj.w      = Value;
        end
        
        
        function obj = set.rmx(obj,Value)
            obj.rmx     = Value;
        end
        
        
        function obj = set.cond(obj,Value)
            obj.cond    = Value;
        end
        
        
        function obj = set.set(obj,Value)
            obj.set     = Value;
        end
        
        function obj = set.part(obj,Value)
            obj.part     = Value;
        end
        
        function f = genNeocassDeck(obj)
            
            id   = sprintf('%-8d', obj.id);
            aid  = sprintf('%-8d', obj.aid);
            p1   = sprintf('%-8d', obj.p1);
            p2   = sprintf('%-8d', obj.p2);
            set  = sprintf('%-8d', obj.set);
            w    = sprintf('%-8d', obj.w);
            rmx  = sprintf('%-8d', obj.rmx);
            cond = sprintf('%-.2e', obj.cond);
            f{1} = sprintf('SPLINE2 %s%s%s%s%s%s%s%s',id,aid,p1,p2,set,w,rmx,cond);
            
            
        end
        
    end
end
