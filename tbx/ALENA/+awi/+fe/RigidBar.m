classdef RigidBar < awi.fe.FEBaseClass
    %RigidBar Describes a rigid connection a single independent grid point
    %and one or more dependent grid points for use in a finite element
    %model.
    %
    % The definition of the 'RigidBar' object matches that of the RBE2 bulk
    % data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Element identification number
        EID
        %Identification number of the independent node
        GN
        %Constrained degrees of freedom of the independent nodes
        CN
        %Identification number of the dependent nodes
        GM
    end
    
    %Store a handle to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.model.Node' object that is the INDEPENDENT node
        %(GA)
        NodesI
        %Handle(s) to the 'awi.model.Node' object(s) that are the DEPENDENT
        %nodes (GB)
        NodesD
    end
    
    methods % set / get
        function set.EID(obj, val)    %set.EID    
            %set.EID Set method for the property 'EID'.
            % 
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.GN(obj, val)     %set.GN     
            validateID(obj, val, 'GN')
            obj.GN = val;
        end
        function set.CN(obj, val)     %set.CN     
            validateDOF(obj, val, 'CN')
            obj.CN = val;
        end
        function set.GM(obj, val)     %set.GM     
            validateattributes(val, {'numeric'}, {'column', 'integer', ...
                'nonnegative', 'finite', 'real'}, class(obj), 'GM');
            obj.GM = val;
        end
        function set.NodesI(obj, val) %set.NodesI 
            %set.NodesI Set method for the property 'NodesI'.
            validateattributes(val, {'awi.fe.Node'}, {'scalar', ...
                'nonempty'}, class(obj), 'NodesI');
            obj.NodesI = val;
        end
        function set.NodesD(obj, val) %set.NodesD 
            %set.NodesD Set method for the property 'NodesD'.
            validateattributes(val, {'awi.fe.Node'}, {'column', ...
                'nonempty'}, class(obj), 'NodesD');
            obj.NodesD = val;
        end
        function val = get.EID(obj)   %get.EID    
            val = obj.ID;
        end
        function val = get.GN(obj)    %get.GN     
            %get.GN Get method for the property 'GN'.
            
            if isempty(obj.NodesI)
                val = obj.GN;
            else
                val = obj.NodesI.ID;
            end
            
        end
        function val = get.GM(obj)    %get.GM     
            %get.GM Get method for the property 'GM'.
            
            if isempty(obj.NodesD)
                val = obj.GM;
            else
                val = vertcat(obj.NodesD.ID);
            end
        end
    end
    
    methods % construction
        function obj = RigidBar
           
            %Make a note of the property names
            addFEProp(obj, 'EID', 'GN', 'CN', 'GM');
            
        end
    end

    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draws the Rigid Bar objects as a set of lines in 3D
            %space. Pads the coordinates with nan terms in order to return
            %a single graphics handle for all beams in the collection.
            %markers.
            %
            % Accepts an array of objects.
            
            hg = [];
            
            %Grab the nodes
            iNodes = {obj.NodesI};
            dNodes = {obj.NodesD}; 
            
            %Remove empties
            idx = or(cellfun(@isempty, iNodes), cellfun(@isempty, dNodes));
            iNodes(idx) = [];
            dNodes(idx) = [];
            iNodes = horzcat(iNodes{:});
            
            if isempty(iNodes) || isempty(dNodes) %Escape route
                return
            end
            
            %Grab coordinates
            iCoords = num2cell(getDrawCoords(iNodes), 1); 
            dCoords = cellfun(@(n) getDrawCoords(n), dNodes, 'Unif', false);
            iCoords = arrayfun(@(i) repmat(iCoords{i}, ...
                [1, size(dCoords{i}, 2)]), 1 : numel(iCoords), 'Unif', false);
            
            %Convert to matrix format
            iCoords = horzcat(iCoords{:});
            dCoords = horzcat(dCoords{:});
            
            %Pad data with nan terms
            x = obj.padCoordsWithNaN([iCoords(1, :) ; dCoords(1, :)]);
            y = obj.padCoordsWithNaN([iCoords(2, :) ; dCoords(2, :)]);
            z = obj.padCoordsWithNaN([iCoords(3, :) ; dCoords(3, :)]);
            
            %Plot
            hg = line(ha, x, y, z, ...
                'Color'    , 'r', ...
                'LineStyle', '-', ...
                'LineWidth', 1  , ...
                'Marker'   , 'o', ...
                'Tag'      , 'Rigid Elements', ...
                'SelectionHighlight', 'off');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.RigidBar' object
            %into a text file using the format of the MSC.Nastran 'RBE2'
            %bulk data entry.           
            
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
                comment = ['RBE2 : Defines a rigid body with independent ', ... 
                    'degrees-of-freedom that are specified at a single '  , ...
                    'grid point and with dependent degrees-of-freedom '   , ...
                    'that are specified at an arbitrary number of grid '  , ...
                    'points.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
               
            %How many objects?
            nObj = numel(obj);
                        
            %Plot in a loop for now...            
            for ii = 1 : nObj               
                %TODO - Vectorise this!
                awi.fe.FEBaseClass.writeListDataCard(fid, 'RBE2', obj(ii), {'GM'});                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

