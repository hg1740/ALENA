classdef LoadSet < awi.fe.FEBaseClass
    %RigidBar Describes a combination of applied loads.
    %
    % The definition of the 'LoadSet' object matches that of the LOAD bulk
    % data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Global scale factor
        S0 = 1.
        %Scale factors for individual loads
        Si = 1;
        %ID number of applied loads
        LoadID
    end
    
    %Store a handle to the 'awi.fe' objects
    properties (Hidden  = true)
        %Handle to the 'awi.fe.PointLoad' objects
        Loads
    end
    
    methods % set / get
        function set.S0(obj, val)      %set.S0
            %set.S0 Set method for the property 'S0'
            
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'nonnan', 'real', 'nonenegative'}, class(obj), 'S0');
            obj.S0 = val;
            
        end
        function set.Si(obj, val)      %set.Si
            %set.Si Set method for the property 'Si'.
            
            validateattributes(val, {'numeric'}, {'row', 'finite', ...
                'nonnan', 'real', 'nonnegative'}, class(obj), 'S0');
            obj.Si = val;
            
        end
        function set.LoadID(obj, val)  %set.LoadID
            %set.LoadID Set method for the property 'LoadID'.
            
            validateattributes(val, {'numeric'}, {'integer', 'column', ...
                'nonnan', 'finite', 'real'}, class(obj), 'LoadID');
            obj.LoadID = val;
            
        end
        function set.Loads(obj, val)   %set.Loads
            %set.Loads Set method for the property 'Loads'.
            
            validateattributes(val, {'awi.fe.PointLoad'}, {'nonempty', ...
                'vector'}, class(obj), 'Loads');
            if iscolumn(val)
                val = val';
            end
            obj.Loads = val;
        end
        function val = get.LoadID(obj) %get.LoadID
            %get.LoadID Get method for the property 'LoadID'.
            
            if isempty(obj.Loads)
                val = obj.LoadID;
            else
                val = [obj.Loads.ID];
            end
        end
    end
    
    methods % construction
        function obj = LoadSet
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'S0', 'Si', 'LoadID');
            
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
                comment = ['LOAD : Defines a static load as a linear ', ...
                    'combination of load sets.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            %Write each object seperately as MASSSET is list-formatted
            for i = 1 : numel(obj)
                
                %Grab data 
                id  = obj(i).ID; 
                s0  = obj(i).S0;
                si  = [obj(i).Si]';
                lid = [obj(i).LoadID]';
                
                awi.fe.FEBaseClass.writeSuperpositionList(fid, 'LOAD', id, s0, si, lid);
                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

