classdef Moment < hgsetget
    %FORCE Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id        = [];
        gid       = [];
        cid       = [];
        momentMag = [];
        x         = [];
        y         = [];
        z         = [];
    end
    
    methods
        
        function obj = Moment(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.gid(obj,Value)
            obj.gid     = Value;
        end
        
        
        function obj = set.cid(obj,Value)
            obj.cid     = Value;
        end
        
        
        function obj = set.momentMag(obj,Value)
            obj.momentMag = Value;
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
            f = cat(1,f,sprintf('MOMENT* %16.0f%16.0f%16.0f%16.9f',obj.Moment.id,obj.Moment.gid,obj.Moment.cid,obj.Moment.Force));
            f = cat(1,f,sprintf('*       %16.9E%16.9f%16.9f',obj.Moment.x,obj.Moment.y,obj.Moment.z));
        end
        
    end
end
