classdef PanelCollection < awi.fe.FECollection
    %PanelCollection Describes a collection of 'awi.fe.Panel' objects.
    %
    % See also awi.fe.Panel
    
    %Primary properties
    properties
        %Panel element ID number
        EID
        %Panel property ID number 
        PID
        %ID of the Nodes associated with this panel
        GridID
        %Orientation of the material in the panel
        THETA = 0;
        %Material coordinate system ID number
        MCID
        %Offset from the surface of the grid points to the element reference plan 
        ZOFFS
        %Logical flag specifying the meaning of the 'Thickness' data
        ThicknessFlag = true;
        %Thickness of the panel at the grid points
        Thickness
    end
    
    %Handle to 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.NodeCollection' object
        NodeCollection
        %Index number of the data in 'NodeCollection' that relates to the
        %grid points in 'GridID'.
        NodeIndex
        %Handle to the 'awi.fe.PanelPropCollection' object
        PanelPropCollection
        %Index number of the data in 'PanelPropCollection' that relates to
        %the elements
        PanelPropIndex
    end
    
    methods % set / get
        function set.EID(obj, val)            %set.EID
           validateIDvector(obj, val, 'EID');
           obj.IDNumbers = val;
        end
        function set.PID(obj, val)            %set.PID
            validateIDvector(obj, val, 'PID');
            obj.PID = val;
        end
        function set.GridID(obj, val)         %set.GridID
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'numeric'}, {'2d', 'nrows', 4, ...
                    'real', 'nonnan', 'finite', 'integer'}, class(obj), 'X');
            end
            obj.GridID = val;
        end
        function set.THETA(obj, val)          %set.THETA
            validateattributes(val, {'numeric'}, {'row', 'real', ...
                'nonnan', 'finite'}, class(obj), 'THETA');
            obj.THETA = val;
        end
        function set.MCID(obj, val)           %set.MCID
            validateIDvector(obj, val, 'MCID');
            obj.MCID = val;
        end
        function set.ZOFFS(obj, val)          %set.ZOFFS
            validateattributes(val, {'numeric'}, {'row', 'real', ...
                'nonnan', 'finite'}, class(obj), 'THETA');
            obj.ZOFFS = val;
        end
        function set.ThicknessFlag(obj, val)  %set.ThicknessFlag
            validateattributes(val, {'logical'}, {'scalar'}, ...
                class(obj), 'ThicknessFlag');
        end
        function set.Thickness(obj, val)      %set.Thickness
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'numeric'}, {'2d', 'nrows', 4, ...
                    'real', 'nonnan', 'finite', 'nonnegative'}, ...
                    class(obj), 'Thickness');
            end
            obj.Thickness = val;
        end
        function set.NodeCollection(obj, val) %set.NodeCollection
            validateattributes(val, {'awi.fe.NodeCollection'}, {'scalar'}, ...
                class(obj), 'NodeCollection');
            obj.NodeCollection = val;
        end
        function set.NodeIndex(obj, val)      %set.NodeIndex
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'numeric'}, {'2d', 'nrows', 4, ...
                    'real', 'nonnan', 'finite', 'integer', 'nonnegative'}, class(obj), 'X');
            end
            obj.NodeIndex = val;
        end
        function set.PanelPropCollection(obj, val) %set.PanelPropCollection
            validateattributes(val, {'awi.fe.PanelPropCollection'}, {'scalar'}, ...
                class(obj), 'PanelPropCollection');
            obj.PanelPropCollection = val;
        end
        function set.PanelPropIndex(obj, val)      %set.PanelPropIndex 
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'numeric'}, {'row', 'real', ...
                    'nonnan', 'finite', 'integer', 'nonnegative'}, class(obj), 'PanelPropIndex');
            end
            obj.PanelPropIndex = val;
        end
        function val = get.EID(obj)                %get.EID
            val = obj.IDnumbers;
        end
        function val = get.PID(obj)                %get.PID
            if isempty(obj.PanelPropCollection)
                val = obj.PID;
            else
                val = obj.PanelPropCollection.IDnumbers;
                if ~isempty(val) && ~isempty(obj.PanelPropIndex)
                   val = val(obj.PanelPropIndex); 
                end
            end
        end
        function val = get.GridID(obj)             %get.GridID
            if isempty(obj.NodeCollection)
                val = obj.GridID;
            else
                val = obj.NodeCollection.IDNumbers;
                if ~isempty(val) && ~isempty(obj.NodeIndex)
                    val = val(:, obj.NodeIndex);
                end
            end
        end
    end
    
    methods % construction
        function obj = PanelCollection 
        
            %Make a note of the property names
            addFEProp(obj, 'EID', 'PID', 'GridID', 'THETA', 'MCID', 'ZOFFS', 'ThicknessFlag', 'Thickness');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draws the PanelCollection object as a set of
            %patches in 3D space. A patch is used so that this function
            %returns only a single graphics handle for all the aerodynamic
            %panels in the collection.
            %
            % Accepts an array of objects.
            
            hg = [];
            
            %Parse
            NC  = get(obj, {'NodeCollection'});
            ind = get(obj, {'NodeIndex'});
            %   - Check for empties
            bEmpty = or(cellfun(@isempty, NC), cellfun(@isempty, ind));
            if bEmpty
                warning(['Unable to to draw the ''%s'' objects as the ', ...
                    '''NodeCollection'' and/or ''NodeIndex'' properties ', ...
                    'have not been defined.'], class(obj));
            end    
            %   - Check for lack of Node coordinate data
            NC  = horzcat(NC{:});            
            bEmpty = arrayfun(@(nc) isempty(nc.X), NC);
            if bEmpty
                warning(['Unable to to draw the ''%s'' objects as the ', ...
                    '''NodeCollection'' does not have coordinate data.'], class(obj));
            end
            
            %Set up patch vertices and faces
            vertices = [NC.X]';
            %   - Update the index number for the node collection
            nPanel = cellfun(@(x) size(x, 2), ind);
            bound  = cumsum([0, nPanel(1 : end - 1)]);
            ind    = arrayfun(@(i) ind{i} + bound(i), 1 : numel(bound), 'Unif', false);
            faces  = horzcat(ind{:})';
            
            %Draw as a single patch
            hg = patch(ha, ...
                'Faces', faces, 'Vertices', vertices, ...
                'EdgeColor', 'k', 'FaceColor', [52, 235, 185] ./ 255 , ... %gray
                'Tag'     , 'Structural Panels', ...                
                'SelectionHighlight', 'off');
            
            %Draw subgroups?
            
            %Draw material coordinate systems?            
            
        end
    end
    
end

