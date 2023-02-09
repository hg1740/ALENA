classdef Info < hgsetget
    %INFO Summary of this class goes here
    %  Detailed explanation goes here
    
    % Properties
    properties
        nGrid    = [];
        nMgrid   = [];
        nPbar    = [];
        nBar     = [];
        nMat     = [];
        nForce   = [];
        nMoment  = [];
        nConm2   = [];
        nCaero   = [];
        nRbe0    = [];
        nSgrid   = [];
        nSets    = [];
        nEngines = [];
        nRbe2    = [];
        nThrust  = [];
        nSpc     = [];
        nLoads   = [];
    end
    
    methods
        
        function obj = Info(~)
            
        end
        
        function obj = setFields(obj,varargin)
            
            if numel(varargin) == 1
                INFO = varargin{1};
                f = fieldnames(INFO);
                
                for i = 1:numel(f)
                    obj.(f{i}) = INFO.(f{i});
                end
            else
                INFO1 = varargin{1};
                INFO2 = varargin{2};
                
                f = fieldnames(INFO1);
                
                for i = 1:numel(f)
                    if INFO1.(f{i}) >= INFO2.(f{i})
                        val = INFO1.(f{i});
                    else
                        val = INFO2.(f{i});
                    end
                    
                    if ~isempty(val)
                        obj.(f{i}) = val;
                    end
                end
                
            end
        end
        
        function obj = set.nGrid(obj,Value)
            obj.nGrid    = Value;
        end
        
        
        function obj = set.nMgrid(obj,Value)
            obj.nMgrid   = Value;
        end
        
        
        function obj = set.nPbar(obj,Value)
            obj.nPbar    = Value;
        end
        
        
        function obj = set.nBar(obj,Value)
            obj.nBar    = Value;
        end
        
        
        function obj = set.nMat(obj,Value)
            obj.nMat    = Value;
        end
        
        function obj = set.nForce(obj,Value)
            obj.nForce   = Value;
        end
        
        function obj = set.nMoment(obj,Value)
            obj.nMoment   = Value;
        end
        
        function obj = set.nConm2(obj,Value)
            obj.nConm2    = Value;
        end
        
        
        function obj = set.nCaero(obj,Value)
            obj.nCaero   = Value;
        end
        
        
        function obj = set.nRbe0(obj,Value)
            obj.nRbe0    = Value;
        end
        
        
        function obj = set.nSgrid(obj,Value)
            obj.nSgrid   = Value;
        end
        
        
        function obj = set.nSets(obj,Value)
            obj.nSets    = Value;
        end
        
        
        function obj = set.nEngines(obj,Value)
            obj.nEngines  = Value;
        end
        
        
        function obj = set.nRbe2(obj,Value)
            obj.nRbe2    = Value;
        end
        
        
        function obj = set.nThrust(obj,Value)
            obj.nThrust   = Value;
        end
        
        
        function obj = set.nSpc(obj,Value)
            obj.nSpc    = Value;
        end
        
        
        function obj = set.nLoads(obj,Value)
            obj.nLoads   = Value;
        end
        
    end
end
