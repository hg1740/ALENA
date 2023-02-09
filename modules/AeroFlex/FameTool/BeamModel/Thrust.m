classdef Thrust < hgsetget
    %THRUST Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        lid = [];
        g   = [];
        cid = [];
        cx  = [];
        cy  = [];
        cz  = [];
        ox  = [];
        oy  = [];
        oz  = [];
    end
    
    methods
        
        function obj = Thrust(~)
            
        end
        
        function obj = set.lid(obj,Value)
            obj.lid     = Value;
        end
        
        
        function obj = set.g(obj,Value)
            obj.g      = Value;
        end
        
        
        function obj = set.cid(obj,Value)
            obj.cid     = Value;
        end
        
        
        function obj = set.cx(obj,Value)
            obj.cx     = Value;
        end
        
        
        function obj = set.cy(obj,Value)
            obj.cy     = Value;
        end
        
        
        function obj = set.cz(obj,Value)
            obj.cz     = Value;
        end
        
        
        function obj = set.ox(obj,Value)
            obj.ox     = Value;
        end
        
        
        function obj = set.oy(obj,Value)
            obj.oy     = Value;
        end
        
        
        function obj = set.oz(obj,Value)
            obj.oz     = Value;
        end
        
    end
end
