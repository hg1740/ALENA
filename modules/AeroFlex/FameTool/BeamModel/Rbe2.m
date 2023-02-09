classdef Rbe2 < hgsetget
    %RBE2 Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id  = [];
        gn  = [];
        dof = [];
        gm  = [];
        parent  = [];
        child  = [];
    end
    
    methods
        
        function obj = Rbe2(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.gn(obj,Value)
            obj.gn     = Value;
        end
        
        
        function obj = set.dof(obj,Value)
            obj.dof     = Value;
        end
        
        
        function obj = set.gm(obj,Value)
            obj.gm     = Value;
        end
        
        function f = genNeocassDeck(obj)
            id   = sprintf('%-8d', obj.id);
            node = sprintf('%-8d', obj.gn);
            dof  = sprintf('%-8d', obj.dof);
            gm   = sprintf('%-8d', obj.gm);
            f{1} = sprintf('RBE2    %s%s%s%s',id,node,dof,gm);
        end
    end
end
