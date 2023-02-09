classdef Beam < awi.fe.FEBaseClass
    %Beam Describes a flexible beam connecting two nodes for use in a
    %finite element model.
    %
    % The definition of the 'Beam' object matches that of the CBEAM bulk 
    % data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Beam element identification number
        EID
        %Beam property idenfitication number
        PID
        %Identification number of the node at end A
        GA
        %Identification number of the node at end B
        GB
        %Orientation vector
        X
        %Shear centre offset in the y-direction of the local beam
        %coordinate system
        SCy
        %Shear centre offset in the z-direction of the local beam
        %coordinate system
        SCz
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node' object that this 'Beam' connects to.
        Nodes    = [];
        %Handle to the 'awi.fe.BeamProp' object
        BeamProperty = [];
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
        function set.SCy(obj, val)          %set.SCy            
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'SCy');
            obj.SCy = val;
        end
        function set.SCz(obj, val)          %set.SCz            
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'SCz');
            obj.SCz = val;
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
        function set.BeamProperty(obj, val) %set.BeamProperty   
            %set.BeamProperty Set method for the property 'BeamProperty'
            %
            % 'BeamProperty' must be a scalar 'awi.fe.BeamProp' object.
            
            validateattributes(val, {'awi.fe.BeamProp', ...
                'awi.fe.BeamCrossSection'}, {'scalar', 'nonempty'}, ...
                class(obj), 'BeamProperty');
            obj.BeamProperty = val;
        end
        function val = get.EID(obj)         %get.EID            
            %get.EID Get method for the property 'EID'.
            val = obj.ID;
        end
        function val = get.PID(obj)         %get.PID            
            %get.PID Get method for the property 'PID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.BeamProp' object then always use their ID number, 
            % else use PID.
            
            if isempty(obj.BeamProperty)
                val = obj.PID;
            else
                val = obj.BeamProperty.ID;                
            end
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
    end
    
    methods % construction
        function obj = Beam 
        
            %Make a note of the property names
            addFEProp(obj, 'EID', 'PID', 'GA', 'GB', 'X', 'SCy', 'SCz');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha, mode)
            %drawElement Draws the Beam objects as a set of lines in 3D
            %space. Pads the coordinates with nan terms in order to return
            %a single graphics handle for all beams in the collection.
            %markers.
            %
            % Accepts an array of objects.
            
            if nargin < 3
               mode = []; 
            end
            
            %Gather node handles
            nodes = {obj.Nodes};
            
            %Remove any that are empty
            idx   = cellfun(@isempty, nodes(1, :));
            nodes = horzcat(nodes{:, ~idx});
            
            %Gather coordinates
            X_A = getDrawCoords(nodes(1, :), mode);
            X_B = getDrawCoords(nodes(2, :), mode);
            
            %Pad with nan terms to allow only a single graphics handle to
            %be used.
            x = obj.padCoordsWithNaN([X_A(1, :) ; X_B(1, :)]);
            y = obj.padCoordsWithNaN([X_A(2, :) ; X_B(2, :)]);
            z = obj.padCoordsWithNaN([X_A(3, :) ; X_B(3, :)]);
            
            %Plot
            hg = line(ha, x, y, z, ...
                'Color'    , 'k', ...
                'LineStyle', '-', ...
                'LineWidth', 2  , ...
                'Tag'      , 'Beam Elements', ...
                'SelectionHighlight', 'off');
            
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Beam' object into a
            %text file using the format of the MSC.Nastran 'CBEAM' bulk
            %data entry.
            %
            % The following assumptions are made:
            %   * The orientation vector is defined in the local
            %   displacement coordinate system of the node.
            %   * There are no pin flags.
            %   * The shear centre is only offset in the local (y,z)
            %   directions.
            %   * There is no warping of the beam.
            
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
                comment = ['CBEAM : Defines a beam element.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                 
            nObj = numel(obj);
            
            %Split up the coordinates            
            offset = [obj.X];
            X1 = num2cell(offset(1, :));
            X2 = num2cell(offset(2, :));
            X3 = num2cell(offset(3, :));
            
            %Grab the shear centre offsets
            SCy_ = [obj.SCy];
            SCz_ = [obj.SCz];
            SCyA = num2cell(SCy_(1, :));
            SCyB = num2cell(SCy_(2, :));
            SCzA = num2cell(SCz_(1, :)); 
            SCzB = num2cell(SCz_(2, :)); 
            
            %Card name
            nam   = repmat({'CBEAM'}, [1, nObj]);
            blnks = repmat({blanks(8)}, [1, nObj]);
            offt  = repmat({'GGG'}, [1, nObj]);
            zrs   = num2cell(zeros(1, nObj));
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.PID} ; {obj.GA} ; {obj.GB}; X1   ; X2  ; X3   ; blnks ; ... 
                blnks ; blnks    ; blnks     ; zrs      ; SCyA    ; SCzA ; zrs ; SCyB ; SCzB];
            
            %Write in 8-character column width 
            format = [ ...
                '%-8s%-8i%-8i%-8i%-8i%#-8.3g%#-8.3g%#-8.3g%-8s\r\n', ...
                '%-8s%-8s%-8s%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

