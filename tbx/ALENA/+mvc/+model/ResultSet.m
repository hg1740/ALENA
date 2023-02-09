classdef (ConstructOnLoad) ResultSet < mvc.model.Thing
    %
    % ResultSet represents the result of performing an analysis on a Thing

    properties (SetAccess = protected) % immutable)
        
        %Useful metadata
        CreatedAt;
        CreatedBy;
        
    end
    
    methods % construction / destruction
        
        function obj = ResultSet(varargin)
        
            %Default args
            args = {'DefaultName', 'ResultSet', ... just a name
                'IsLeafNode', true, ... Result Sets inherrently have no children (? always true ?)
                'StructureLocked', true, ... follows implicitly from IsLeafNode, not really necessary
                'PropertiesLocked', true, ... do not allow user to edit content
                'PropertiesLockable', false ... user is not granted direct access to change editability
                };
            
            %Pass it on
            obj@mvc.model.Thing(args{:}, varargin{:});
            
            %When, and by whom ?
            obj.CreatedBy = getenv('UserName');
            obj.CreatedAt = datetime();
            
            %Extend property groups
            obj.addPropertyGroup('General', ...
                'CreatedBy', 'CreatedBy', ...
                'CreatedAt', 'CreatedAt');
            
        end
        
    end
    
end
