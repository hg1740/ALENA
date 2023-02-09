classdef Sets < hgsetget
    %SETS Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id   = [];
        data = [];
        part = [];
    end
    
    methods
        
        function obj = Sets(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.data(obj,Value)
            obj.data    = Value;
        end
        
        function f = genNeocassDeck(obj)
            setVec = cat(1,obj.id,obj.data);
            setVec = sprintf('%-8d',setVec);
            setVec = cat(2,setVec,blanks(9 * 8));
            
            f = {};
            while ~isempty(deblank(setVec))
                f = cat(1,f,sprintf('        %s',setVec(1:8 * 8)));
                setVec(1:8 * 8) = [];
            end
            
            f{1}(1:4) = 'SET1';
        end
        
    end
end
