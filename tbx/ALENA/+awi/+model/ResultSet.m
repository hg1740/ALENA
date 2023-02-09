classdef (ConstructOnLoad) ResultSet < awi.model.Entity
    %ResultSet Represents the result of performing an analysis on a model.
    
    %Meta Data   
    properties (SetAccess = protected)  
        %Time of creation
        CreatedAt;
        %Author
        CreatedBy;        
    end
    
    %The actual results objects...
    properties %BeamResults
        %Aggregate an array of 'awi.model.BeamResult' for viewing later
        BeamResults         
    end
    
    methods % set / get
        function set.BeamResults(obj, val) %set.BeamResults 
            if iscolumn(val)
                val = val';
            end
            validateattributes(val, {'awi.model.BeamResult'}, {'row'}, ...
                class(obj), 'BeamResults');
            obj.BeamResults = val;
        end
    end
    
    methods % construction / destruction
        
        function obj = ResultSet(varargin)
            %ResultSet Constructor for the class 'ResultSet'. 
            %
            % An instance of 'awi.model.ResultSet' is a leaf-node and its 
            % structure and properties are locked.
            
            %Default args
            args = { ... 'Name', 'ResultSet', ... just a name
                'IsLeafNode'        , true, ... Result Sets inherrently have no children (? always true ?)
                'StructureLocked'   , true, ... follows implicitly from IsLeafNode, not really necessary
                'PropertiesLocked'  , true, ... do not allow user to edit content
                'PropertiesLockable', false ... user is not granted direct access to change editability
                };
            
            %Pass it on
            obj@awi.model.Entity(args{:}, varargin{:});
            
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
