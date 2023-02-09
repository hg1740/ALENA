classdef TransientResult < awi.model.Collector & awi.mixin.Beamable
    %TransientResult Collector class for all transient results sets.
    
    properties
    end
    
    methods % constructor
        function obj = TransientResult(varargin)
            %TransientResult Constructor for the 'TransientResult' class.
            
            %Start with base class
            obj@awi.model.Collector(varargin{:});
            
            %Add collectables
            obj.addCollectionSpec(@awi.model.ResultSet, 'Results Sets');
            
            %Could potentially have a lot of transient results so don't
            %draw children in the tree as standard
            obj.DrawChildrenInTree = false;
            
            %Configure dynamic properties?
            
            
        end
    end
    
end

