classdef (ConstructOnLoad) Collector < mvc.model.DrawableThing 
    %
    % Collector represents a node in a hierarchical model that serves no purpose other than to collect child nodes.
    %  It implements everything supported by the Thing, which means it is:
    %
    %  Collectable, so that objects of this type can be arranged in a hierarchical fashion (with parent/child relationships)
    %  Searchable, so that the children in collects are searchable from higher in the tree
    %  Contextable, so that the content of the object array can present context-sensitive options to the user
    %
    % But on construction the Collector suppresses any detail of property groups and removes propedit from context.
    
    properties (AbortSet, SetObservable, Dependent)
        
        %Does this collector draw its children in tree view ?
        % (this dependent attribute mirrors that on Collectable, hence the underscore in name)
        DrawChildrenInTree_ = true;

    end

    methods % set/get
        
        function val = get.DrawChildrenInTree_(obj)
            
            %Pass it on
            val = obj.DrawChildrenInTree;
            
        end
        
        function set.DrawChildrenInTree_(obj, val)
            
            %Pass it on
            obj.DrawChildrenInTree = val;
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Collector(varargin)

            %Start with base class
            obj@mvc.model.DrawableThing('Visible', false, varargin{:});
            
            %Suppress any detail of property groups
            obj.PropertyGroups(1:numel(obj.PropertyGroups)) = [];
            
            %But treat "Name" and "DrawInTree" as displayable property
            % (TODO: consider whether Name should be uneditable, for a collector ?)
            obj.addPropertyGroup('General', ...
                'Name', 'Name', ...
                'DrawChildrenInTree_', 'Draw Children in Tree');
            
        end
        
    end
    
end
