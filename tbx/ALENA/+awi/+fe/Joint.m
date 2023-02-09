classdef Joint < awi.fe.FEBaseClass
    %Joint Describes a joint connecting two coincident nodes for use in a
    %finite element model.
    %
    % The definition of the 'Joint' object matches that of the RJOINT bulk 
    % data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Element identification number
        EID
        %ID number of the GA node
        GA
        %ID number of the GB node
        GB
        %Constrained DOFs in the global coordinate system at GB
        CB
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node' objects that the joint connect to
        Nodes
    end
    
    methods % set / get
        function set.EID(obj, val)   %set.EID   
            %set.EID Set method for the property 'EID'.
            % 
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.GA(obj, val)    %set.GA    
            validateID(obj, val, 'GA')
            obj.GA = val;
        end
        function set.GB(obj, val)    %set.GB    
            validateID(obj, val, 'GB')
            obj.GB = val;
        end
        function set.CB(obj, val)    %set.CB    
           validateDOF(obj, val, 'CB')
           if ischar(val)
               val = str2double(val);
           end
           obj.CB = val;
        end
        function set.Nodes(obj, val) %set.Nodes 
            %set.Nodes Set method for the 'Nodes' property.
            %
            % 'Nodes' must be a [2,1] vector of 'awi.fe.Node' objects.
            
            assert(numel(val) <= 2, sprintf('An object of class ''%s'' ', ...
                'can only connect to 2 nodes.', class(obj))); 
            validateattributes(val, {'awi.fe.Node'}, {'column', 'nonempty'}, ...
                class(obj), 'Nodes')
            obj.Nodes = val;            
        end
        function val = get.GA(obj)   %get.GA    
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
        function val = get.GB(obj)   %get.GB    
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
        function obj = Joint 
        
            %Make a note of the property names
            addFEProp(obj, 'EID', 'GA', 'GB', 'CB');
            
        end
    end
    
    methods % visualisation 
        function hg = drawElement(obj, ha)
            %drawElement Draws the Joint objects as a set of discrete
            %markers.
            %
            % Accepts an array of objects.
            
            hg = [];
            
            %Grab the node objects
            jNodes  = {obj.Nodes};
            jNodes(cellfun(@isempty, jNodes)) = [];
            jNodes = horzcat(jNodes{:});
            
            if isempty(jNodes)
                return
            end
            
            jCoords = [jNodes(1, :).X];
            
            hg = line(ha, jCoords(1, :), jCoords(2, :), jCoords(3, :), ...
                'LineStyle', 'none', ...
                'Marker'   , 's'   , ...
                'MarkerFaceColor', 'm', ...
                'MarkerEdgeColor', 'k', ...
                'Tag', 'Joints');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Joint' object into
            %a text file using the format of the MSC.Nastran RJOINT bulk
            %data entry.           
            
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
                comment = ['RJOINT : Defines a rigid joint element ', ...
                    'connecting two coinciding grid points.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'RJOINT'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.GA} ; {obj.GB} ; {obj.CB}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-8i%-8i%-8i%-8i\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

