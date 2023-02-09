classdef PartId < hgsetget
    %PARTID Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id   = [];
        part = [];
        type = [];
        data = [];
    end
    
    methods
        
        function obj = PartId(~)
            
        end
        
        function obj = PartId.id(obj,Value)
            obj.id     = Value;
        end
         
        function obj = PartId.data(obj,Value)
            obj.data    = Value;
        end
        
        function obj = PartId.part(obj,Value)
            obj.data    = Value;
        end
        
        function obj = PartId.type(obj,Value)
            obj.data    = Value;
        end
        
        function f = genNeocassDeck(obj)
            setVec = cat(1,obj.data);
            setVec = sprintf('%-8d',setVec);
            setVec = cat(2,setVec,blanks(9 * 8));
            setID   = sprintf('%-8d',obj.id);
            setPart = sprintf('%-8s',obj.part);
            setType = sprintf('%-8s',obj.type);
            
            setVec = [setID,setPart,setType,setVec];
            
            f = {};
            while ~isempty(deblank(setVec))
                f = cat(1,f,sprintf('        %s',setVec(1:8 * 8)));
                setVec(1:8 * 8) = [];
            end
            
            f{1}(1:6) = 'PARTID';
        end
        
    end
end
