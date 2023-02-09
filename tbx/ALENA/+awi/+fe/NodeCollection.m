classdef NodeCollection < awi.fe.FECollection
    %NodeCollection Describes a collection of 'awi.fe.Node' objects.
    %
    % See also awi.fe.Node
    
    %Primary properties
    properties
        %Identification number
        GID 
        %Definition coordinates system identification number
        CP = 0;
        %Coordinates of the node in the coordinate system defined by CP
        X  = [0 ; 0 ; 0];
        %Output coordinate system identification number
        CD = 0;
    end
    
    methods % set / get
        function set.GID(obj, val) %set.GID 
           obj.IDNumbers = val;
        end
        function set.CP(obj, val)  %set.CP  
            validateIDvector(obj, val, 'CP');
            obj.CP = val;
        end
        function set.X(obj, val)   %set.X   
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                    'real', 'nonnan', 'finite'}, class(obj), 'X');
            end
            obj.X = val;
        end
        function set.CD(obj, val)  %set.CD  
            validateIDvector(obj, val, 'CD');
            obj.CD = val;
        end
    end
    
    methods % construction 
        function obj = NodeCollection
           %Make a note of the property names
            addFEProp(obj, 'GID', 'CP', 'X', 'CD');
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha, mode)
            %drawElement Draws a collection of nodes as discrete markers
            %and returns a single graphics handle for all the nodes in the
            %collection. 
            %
            % Accepts a vector of objects.
            
            if nargin < 3
                mode = 'undeformed';
            end

            %coords = getDrawCoords(obj, mode); 
            coords = [obj.X];
            hg = drawNodes(coords, ha);
            
        end
    end
    
end

