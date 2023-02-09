%% FORCE 
%  Detailed explanation goes here
classdef Force < hgsetget
    
    % Properties
    properties
        id       = [];
        gid      = [];
        cid      = [];
        forceMag = [];
        x        = [];
        y        = [];
        z        = [];
    end
    
    methods
        
        function obj = Force(~)
            
        end
        
        %% Set the id property
        % Set the id property here
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        %% Set the gid property
        function obj = set.gid(obj,Value)
            obj.gid     = Value;
        end
        
        
        function obj = set.cid(obj,Value)
            obj.cid     = Value;
        end
        
        
        function obj = set.forceMag(obj,Value)
            obj.forceMag = Value;
        end
        
        
        function obj = set.x(obj,Value)
            obj.x      = Value;
        end
        
        
        function obj = set.y(obj,Value)
            obj.y      = Value;
        end
        
        
        function obj = set.z(obj,Value)
            obj.z      = Value;
        end
        
        function f = genNastranDeck(obj)
            f = {};
            f = cat(1,f,sprintf('FORCE*  %16.0f%16.0f%16.0f%16.9f',obj.Force.id,obj.Force.gid,obj.Force.cid,obj.Force.Force));
            f = cat(1,f,sprintf('*       %16.9E%16.9f%16.9f',obj.Force.x,obj.Force.y,obj.Force.z));
        end
    end
end
