classdef Mat < hgsetget
    %MAT Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        id   = [];
        e    = [];
        g    = [];
        nu   = [];
        rho  = [];
        a    = [];
        tRef = 0;
        ge   = 0.02;
        st   = 0;
        sc   = 0;
        ss   = 0;
    end
    
    methods
        
        function obj = Mat(~)
            
        end
        
        function obj = set.id(obj,Value)
            obj.id     = Value;
        end
        
        
        function obj = set.e(obj,Value)
            obj.e      = Value;
        end
        
        
        function obj = set.g(obj,Value)
            obj.g      = Value;
        end
        
        
        function obj = set.nu(obj,Value)
            obj.nu     = Value;
            obj.g      = obj.e / (2 * (1 + obj.nu));
        end
        
        
        function obj = set.rho(obj,Value)
            obj.rho     = Value;
        end
        
        function obj = set.a(obj,Value)
            obj.a     = Value;
        end
        
        function obj = set.tRef(obj,Value)
            obj.tRef     = Value;
        end
        
        function obj = set.ge(obj,Value)
            obj.ge     = Value;
        end
        
        function obj = set.st(obj,Value)
            obj.st     = Value;
        end
        
        
        function obj = set.sc(obj,Value)
            obj.sc     = Value;
        end
        
        
        function obj = set.ss(obj,Value)
            obj.ss     = Value;
        end
        
        function f = genNeocassDeck(obj)
            
            id  = sprintf('%-8d', obj.id);
            e   = sprintf('%-.2e', obj.e);
            nu  = sprintf('%-8g', obj.nu);
            rho = sprintf('%-8g', obj.rho);
            st  = sprintf('%-8g', obj.st);
            sc  = sprintf('%-8g', obj.sc);
            ss  = sprintf('%-8g', obj.ss);
            
            f = {};
            f{1} = sprintf('MAT1    %s%s        %s%s%s%s%s', id, e, nu, rho);
            f{2} = sprintf('        %s%s%s', st,sc,sc);         
        end
    end
end
