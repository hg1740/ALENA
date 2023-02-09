classdef StructuralSet < awi.fe.FEBaseClass
    %StructuralSet Defines a collection of 'awi.fe.Node' objects.
    %
    % The definition of the 'StructuralSet' object matches that of the SET1
    % bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %ID number of this structural set
        SID 
        %ID numbers of the nodes belonging to this set
        NodeID
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node' objects
        Nodes
    end
    
    methods % set / get
        function set.SID(obj, val)     %set.SID    
            %set.SID Set method for the property 'SID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;            
        end
        function set.NodeID(obj, val)  %set.NodeID 
            %set.NodeID Set method for the propery 'NodeID'.
            validateattributes(val, {'numeric'}, {'column', 'integer', ...
                'nonnegative', 'finite', 'real'}, class(obj), 'NodeID');
            obj.NodeID = val;            
        end
        function set.Nodes(obj, val)   %set.Nodes  
            %set.Nodes Set method for the property 'Nodes'.
            validateattributes(val, {'awi.fe.Node'}, {'column', ...
                'nonempty'}, class(obj), 'Nodes');
            obj.Nodes = val;
        end
        function val = get.SID(obj)    %get.SID    
            %get.SID Get method for property 'SID'.
            val = obj.ID;
        end
        function val = get.NodeID(obj) %get.NodeID 
            %get.NodeID Get method for the property 'NodeID'.
            %
            % If the object has been assigned a handle to its 'awi.fe.Node' 
            % objects then always use their ID number, else use 'NodeID'.
            if isempty(obj.Nodes)
                val = obj.NodeID;
            else
                val = vertcat(obj.Nodes.ID);
            end
        end
    end
    
    methods % construction
        function obj = StructuralSet
            
            %Make a note of the property names
            addFEProp(obj, 'SID', 'NodeID');

        end

    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.StructuralSet'
            %object into a text file using the format of the MSC.Nastran
            %'SET1' bulk data entry.
            
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
                comment = ['SET1 : Defines a list of structural grid ', ...
                    'points or element identification numbers..'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            %How many objects?
            nObj = numel(obj);
            
            %Plot in a loop for now...
            for ii = 1 : nObj
                %TODO - Vectorise this!
                awi.fe.FEBaseClass.writeListDataCard(fid, 'SET1', obj(ii), {'NodeID'});
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

