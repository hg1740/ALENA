classdef PointMass < awi.fe.FEBaseClass
    %PointMass Defines a mass located at a point in 3D space
    %
    % The definition of the 'PointMass' object matches that of the CONM2
    % bulk data type from MSC.Nastran
    
    %Primary Properties     
    properties
        %Identification number
        EID
        %Identification number of the node that the point mass is
        %associated with
        G
        %Identification number of the coordinate system defining the 
        %direction of the offsets
        CID = 0;
        %Mass of the point mass
        M = 0;
        %Offset distance in direction [1,2,3]
        X = [0 ; 0 ; 0];
        %Inertia value in direction 11
        I11 = 0;
        %Inertia value in direction 21
        I21 = 0;
        %Inertia value in direction 22
        I22 = 0;
        %Inertia value in direction 31
        I31 = 0;
        %Inertia value in direction 32
        I32 = 0;
        %Inertia value in direction 33
        I33 = 0;
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node. object that this mass relates to.
        Node    = [];
%         %Handle to the 'awi.fe.CoordSys' object
%         CoordSys = [];
    end
    
    methods % set / get
        function set.EID(obj, val)  %set.EID  
            %set.EID Set method for the property 'EID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.G(obj, val)    %set.G    
            validateID(obj, val, 'G');
            obj.G = val;
        end
        function set.CID(obj, val)  %set.CID  
            validateID(obj, val, 'CID')
            obj.CID = val;
        end
        function set.M(obj, val)    %set.M    
            %set.M Set method for the property 'M'.
            %
            % 'M' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'nonnan', 'finite', 'real'}, class(obj), 'M');
            obj.M = val;
        end
        function set.X(obj, val)    %set.X    
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'X');
            obj.X = val;
        end
        function set.I11(obj, val)  %set.I11  
            %set.I11 Set method for the property 'I11'.
            %
            % 'I11' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'nonnan', 'finite', 'real'}, class(obj), 'I11');
            obj.I11 = val;
        end
        function set.I21(obj, val)  %set.I21  
            %set.I21 Set method for the property 'M'.
            %
            % 'I21' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, class(obj), 'I21');
            obj.I21 = val;
        end
        function set.I22(obj, val)  %set.I22  
            %set.I22 Set method for the property 'I22'.
            %
            % 'I22' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'nonnan', 'finite', 'real'}, class(obj), 'I22');
            obj.I22 = val;
        end
        function set.I31(obj, val)  %set.I31  
            %set.I31 Set method for the property 'I31'.
            %
            % 'I31' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, class(obj), 'I31');
            obj.I31 = val;
        end
        function set.I32(obj, val)  %set.M    
            %set.M Set method for the property 'I32'.
            %
            % 'I32' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, class(obj), 'I32');
            obj.I32 = val;
        end
        function set.I33(obj, val)  %set.I33  
            %set.I33 Set method for the property 'I33'.
            %
            % 'I33' must be a positive, scalar, numeric.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'nonnan', 'finite', 'real'}, class(obj), 'I33');
            obj.I33 = val;
        end
        function set.Node(obj, val) %set.Node 
           %set.Node Set method for the property 'Node'.
           %
           % 'Node' must be a scalar 'awi.fe.PointMass' object
           
           if isempty(val)
               obj.Node = [];
               return
           end
           validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
               class(obj), 'Node');
           obj.Node = val;
        end
        function val = get.EID(obj) %get.EID  
            val = obj.ID;
        end
        function val = get.G(obj)   %get.G    
            %get.G Get method for the property 'G'.
            %
            % If the object has been assigned a handle to its 'awi.fe.Node'
            % object then always use their ID numbers, else use G.
            
            if isempty(obj.Node)
                val = obj.G;
            else
                val = obj.Node(1).GID;
            end
        end
    end
    
    methods % constructor
        function obj = PointMass
            
            %Make a note of the property names
            addFEProp(obj, 'EID', 'G', 'CID', 'M', 'X', 'i11', 'I21', ...
                'I22', 'I31', 'I32', 'I33');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draws the point mass objects as a set of discrete
            %markers. Returns a single graphics handle for all point masses
            %in the collection.
            %
            % Accepts an array of objects.
            
            hg = [];
            
            %Gather the nodes.
            nodes  = {obj.Node};
            
            %Only retain the masses that have an associated node
            idx    = cellfun(@(x) ~isempty(x), nodes);
            nodes  = horzcat(nodes{idx});
            
            if isempty(nodes)
                return
            end
            
            %Generate coordinates - include offset terms
            massCoords = [nodes.X] + [obj(idx).X];
            xm = massCoords(1, :);
            ym = massCoords(2, :);
            zm = massCoords(3, :);
            
            %Plot
            hg = line(ha, xm, ym, zm, ...
                'LineStyle'      , 'none'  , ...
                'Marker'         , '^'     , ...
                'MarkerFaceColor', 'b'     , ...
                'MarkerEdgeColor', 'k'     , ...
                'Tag'            , 'Masses', ...
                'SelectionHighlight', 'off');           
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Node' object into a
            %text file using the format of the MSC.Nastran 'GRID' bulk data
            %entry. 
            %
            % The following assumptions are made:
            %   * The PS and SEID properties are omitted. i.e. There are no
            %   permanenent single point constraints and all nodes are
            %   assumed to belong to the same super element.            
            
            %By default, do not close the file
            bClose = false;
                        
            if nargin < 2 %Ask the user for the file
                fName = awi.fe.FEBaseClass.getBulkDataFile;
                bClose = true;
                fid = fopen(fName, 'w');                
            end
            
            if nargin < 3 %Comments by standard
                bComment = true;
            end
                        
            if bComment %Helpful comments?
                comment = ['CONM2 : Defines a concentrated mass at ', ...
                    'a grid point'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
                  
            nObj = numel(obj);
            
            %Split up the coordinates            
            coords = [obj.X];
            X1 = num2cell(coords(1, :));
            X2 = num2cell(coords(2, :));
            X3 = num2cell(coords(3, :));
            
            %Card name
            nam    = repmat({'CONM2*'}  , [1, nObj]);
            blnks  = repmat({['*', blanks(7)]}, [1, nObj]);
            blnks_ = repmat({blanks(16)}, [1, nObj]);
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID}  ; {obj.G}   ; {obj.CID} ; {obj.M}   ; ...
                blnks ; X1        ; X2        ; X3        ; blnks_    ; ...
                blnks ; {obj.I11} ; {obj.I21} ; {obj.I22} ; {obj.I31} ; ...
                blnks ; {obj.I32} ; {obj.I33} ];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%-16i%-16i%#-16.8g\r\n'         , ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%-16s\r\n'   , ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

