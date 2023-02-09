classdef (ConstructOnLoad) Collector < awi.model.Entity
    %
    % Collector represents a node in a hierarchical model that serves no purpose other than to collect child nodes.
    %  It implements everything supported by the Thing, which means it is:
    %
    %  Collectable, so that objects of this type can be arranged in a hierarchical fashion (with parent/child relationships)
    %  Searchable, so that the children in collects are searchable from higher in the tree
    %  Contextable, so that the content of the object array can present context-sensitive options to the user
    %
    % But on construction the Collector suppresses any detail of property groups and removes propedit from context.
    
    methods % construction / destruction
        
        function obj = Collector(varargin)

            %Start with base class
            obj@awi.model.Entity(varargin{:});
            
            %Suppress any detail of property groups
            obj.PropertyGroups(1 : numel(obj.PropertyGroups)) = [];
            
            %Remove propedit from context
            
        end
        
    end
    
end
