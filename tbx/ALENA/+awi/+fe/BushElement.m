classdef BushElement < awi.fe.FEBaseClass
    %BushElement Describes a generic bush element connecting two nodes.
    %
    % The definition of the 'BushElement' object matches that of the CBUSH
    % bulk data type from MSC.Nastran
    
    %Primary Properties
    properties
        %Bush element identification number
        EID
        %Bush property idenfitication number
        PID
        %Identification number of the node at end A
        GA
        %Identification number of the node at end B
        GB
        %Orientation vector
        X
        %Coordinate System
        CID=0;
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node. object that this 'Bush' connects to
        Nodes    = [];
        %Handle to the 'awi.fe.BushProp' object
        BushProperty = [];
        %coordinate system
        CoordSys = [];
    end
    
    methods % set / get
        function set.EID(obj, val)          %set.EID          
            %set.EID Set method for the property 'EID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.PID(obj, val)          %set.PID          
            validateID(obj, val, 'PID')
            obj.PID = val;
        end
        function set.GA(obj, val)           %set.GA           
            validateID(obj, val, 'GA')
            obj.GA = val;
        end
        function set.GB(obj, val)           %set.GB           
            validateID(obj, val, 'GB')
            obj.GB = val;
        end
        function set.X(obj, val)            %set.X            
            validateattributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                'real', 'finite', 'nonnan'}, class(obj), 'X');
            obj.X = val;
        end
        function set.CoordSys(obj, val)    %set.OutputCoordSys
            %set.OutputCoordSys Set method for the property 'OutputCoordSys'.
            %
            % 'OutputCoordSys' must be a scalar instance of 'awi.fe.CoordSys'.
            validateattributes(val, {'awi.fe.CoordSys'}, {'scalar', ...
                'nonempty'}, class(obj), 'CoordSys');
            obj.CoordSys = val;
        end
        function set.Nodes(obj, val)        %set.Nodes        
            %set.Nodes Set method for the 'Nodes' property.
            %
            % 'Nodes' must be a [2,1] vector of 'awi.fe.Node' objects.
            
            assert(numel(val) <= 2, sprintf('An object of class ''%s'' ', ...
                'can only connect to 2 nodes.', class(obj))); 
            validateattributes(val, {'awi.fe.Node'}, {'column', 'nonempty'}, ...
                class(obj), 'Nodes')
            obj.Nodes = val;            
        end
        function set.BushProperty(obj, val) %set.BeamProperty 
            %set.BushProperty Set method for the property 'BushProperty'
            %
            % 'BushProperty' must be a scalar 'awi.fe.BushProp' object.
            
            validateattributes(val, {'awi.fe.BushProp'}, {'scalar', ...
                'nonempty'}, class(obj), 'BushProperty');
            obj.BushProperty = val;
        end
        function val = get.EID(obj)         %get.EID          
            %get.EID Get method for the property 'EID'.
            val = obj.ID;
        end
        function val = get.GA(obj)          %get.GA           
            %get.GA Get method for the property 'GA'.
            %
            % If the object has been assigned a handle to its 'awi.fe.Node'
            % object then always use their ID numbers, else use GA/GB.
            
            if isempty(obj.Nodes)
                val = obj.GA;
            else
                val = obj.Nodes(1).GID;                
            end
        end
        function val = get.GB(obj)          %get.GB           
            %get.GB Get method for the property 'GB'.
            %
            % If the object has been assigned a handle to its 'awi.fe.Node'
            % object then always use their ID numbers, else use GA/GB.
            
            if numel(obj.Nodes) == 2 
                val = obj.Nodes(2).GID;
            else
                val = obj.GB;
            end
        end
        function val = get.PID(obj)         %get.PID          
            %get.PID Get method for the property 'PID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.BushProp' object then always use their ID number, 
            % else use PID.
            
            if isempty(obj.BushProperty)
                val = obj.PID;
            else
                val = obj.BushProperty.ID;                
            end
        end
        function val = get.CID(obj)               %get.CID
            %get.CID Get method for the property 'CID'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.CoordSys' object then always use their ID numbers,
            % else use CID.
            
            if isempty(obj.CoordSys)
                val = obj.CID;
            else
                val = obj.CoordSys.CID;
            end
        end
    end
    
    methods % construction
        function obj = BushElement 
        
            %Make a note of the property names
            addFEProp(obj, 'EID', 'PID', 'GA', 'GB', 'X','CID');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draws the Bush objects as a set of lines in 3D
            %space. Pads the coordinates with nan terms in order to return
            %a single graphics handle for all beams in the collection.
            %markers.
            %
            % Accepts an array of objects.
            
            hg = [];
            
            %Grab the nodes - Do not plot any objects that do not have a
            %reference to a 'awi.fe.Node' object
            bNodes = {obj.Nodes};
            bNodes(cellfun(@isempty, bNodes)) = [];
            bNodes = horzcat(bNodes{:});
            if isempty(bNodes)
                return
            end
            
            %Grab coordinates - TODO Check for empty cells!
            bCoords = get(bNodes, 'X');
            bCoords = reshape(bCoords, size(bNodes));
            endA = horzcat(bCoords{1, :});
            endB = horzcat(bCoords{2, :});
            
            %Pad with NaN terms to allow vectorised plotting
            xBu  = obj.padCoordsWithNaN([endA(1, :) ; endB(1, :)]);
            yBu  = obj.padCoordsWithNaN([endA(2, :) ; endB(2, :)]);
            zBu  = obj.padCoordsWithNaN([endA(3, :) ; endB(3, :)]);
            
            %Plot the line!
            hg  = line(ha, xBu, yBu, zBu, ...
                'Color'    , 'b', ...
                'LineStyle', '--', ...
                'LineWidth', 2  , ...
                'Tag'      , 'Bush Elements');            
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.BushElement' object
            %into a text file using the format of the MSC.Nastran 'CBUSH'
            %bulk data entry.
            %
            % The following assumptions are made:
            %   * The element coordinate system is assumed to be aligned
            %   with the basic coordinate system.
            %   * The spring/damper is assumed to be located at the
            %   mid-point of the two nodes. (The default setting)
            %   * There is no offset of the spring/damper element.
            
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
                comment = ['CBUSH : Defines a generalized '         , ...
                    'spring-and-damper structural element that may ', ...
                    'be nonlinear or frequency dependent.'];
                awi.fe.FEBaseClass.writeComment(comment, fid); 
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            nObj = numel(obj);    
            
            %Card name
            nam   = repmat({'CBUSH'}   , 1, nObj);
            blnks = repmat({ blanks(8)}, 1, nObj);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.PID} ; {obj.GA} ; {obj.GB} ; blnks ; blnks ; blnks ; {obj.CID}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-8i%-8i%-8i%-8i%-8s%-8s%-8s%-8i\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

