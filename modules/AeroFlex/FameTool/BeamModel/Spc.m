classdef Spc < hgsetget
    %SPC Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id    = [];
        grids = [];
        dof   = [];
        d     = [];
        part  = [];
    end
    
    methods
        
        function obj = Spc(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        function obj = set.grids(obj,Value)
            obj.grids    = Value;
        end
        
        function obj = set.dof(obj,Value)
            obj.dof     = Value;
        end
        
        function obj = set.d(obj,Value)
            obj.d     = Value;
        end
        
        % TODO: Check that it is SPC or SPC1
        function f = genNastranDeck(obj)
            f = {};
            f = cat(1,f,sprintf('SPC     %8.0f%8.0f%8.0f%8.1f',obj.id,obj.grids,obj.dof,obj.d)); 
        end
        
        function f = genNeocassDeck(obj)
            f = {};
            f = cat(1,f,sprintf('SPC1    %8.0f%8.0f%8.0f',obj.id,obj.dof,obj.grids));
        end
    end
end
