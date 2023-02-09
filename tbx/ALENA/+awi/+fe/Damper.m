classdef Damper < awi.fe.FEBaseClass
    %Damper Describes a generic spring element connecting two nodes.
    %
    % The definition of the 'Damper' object matches that of the
    % CDAMP1 bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Damper element identification number
        EID
        %Damper property idenfitication number
        PID
        %Identification number of the node at end A
        GA
        %Connecting degree of freedom at end A
        CA
        %Identification number of the node at end B
        GB
        %Connecting degree of freedom at end A
        CB
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Node. object that this 'Damper' connects to
        Nodes          = [];
        %Handle to the 'awi.fe.DamperProp' object
        DamperProperty = [];
    end
    
    methods % set / get
        function set.EID(obj, val)            %set.EID            
            %set.EID Set method for the property 'EID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.PID(obj, val)            %set.PID            
            validateID(obj, val, 'PID')
            obj.PID = val;
        end
        function set.GA(obj, val)             %set.GA             
            validateID(obj, val, 'GA')
            obj.GA = val;
        end
        function set.CA(obj, val)             %set.CA             
           %set.CA Set method for the property 'CA'.
           %
           % 'CA' must be a valid degree of freedom.
           validateDOF(obj, val, 'CA');
           obj.CA = val;
        end
        function set.GB(obj, val)             %set.GB             
            validateID(obj, val, 'GB')
            obj.GB = val;
        end 
        function set.CB(obj, val)             %set.CB             
           %set.CB Set method for the property 'CB'.
           %
           % 'CB' must be a valid degree of freedom.
           validateDOF(obj, val, 'CB');
           obj.CB = val;
        end
        function set.Nodes(obj, val)          %set.Nodes          
            %set.Nodes Set method for the 'Nodes' property.
            %
            % 'Nodes' must be a [2,1] vector of 'awi.fe.Node' objects.
            
            assert(numel(val) <= 2, sprintf('An object of class ''%s'' ', ...
                'can only connect to 2 nodes.', class(obj))); 
            validateattributes(val, {'awi.fe.Node', 'awi.fe.ScalarPoint'}, {'column', 'nonempty'}, ...
                class(obj), 'Nodes')
            obj.Nodes = val;            
        end
        function set.DamperProperty(obj, val) %set.DamperProperty 
            %set.DamperProperty Set method for the property 'DamperProperty'
            %
            % 'DamperProperty' must be a scalar 'awi.fe.DamperProp' object.
            
            validateattributes(val, {'awi.fe.DamperProp'}, { ...
                'scalar', 'nonempty'}, class(obj), 'DamperProperty');
            obj.DamperProperty = val;
        end
        function val = get.EID(obj)           %get.EID            
            %get.EID Get method for the property 'EID'.
            val = obj.ID;
        end
        function val = get.GA(obj)            %get.GA             
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
        function val = get.GB(obj)            %get.GB             
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
        function val = get.PID(obj)           %get.PID            
            %get.PID Get method for the property 'PID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.DamperProp' object then always use their ID number, 
            % else use PID.
            
            if isempty(obj.DamperProperty)
                val = obj.PID;
            else
                val = obj.DamperProperty.ID;                
            end
        end
    end
    
    methods % construction
        function obj = Damper
            
            %Make a note of the property names
            addFEProp(obj, 'EID', 'PID', 'GA', 'CA', 'GB', 'CB');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.Damper' object into
            %a text file using the format of the MSC.Nastran CDAMP1 bulk
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
                comment = 'CDAMP1 : Defines a scalar damper element';
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'CDAMP1'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.PID} ; {obj.GA} ; {obj.CA} ; {obj.GB} ; {obj.CB}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-8i%-8i%-8i%-8i%-8i%-8i\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

