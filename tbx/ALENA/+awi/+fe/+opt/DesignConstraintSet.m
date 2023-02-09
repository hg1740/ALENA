classdef DesignConstraintSet < awi.fe.FEBaseClass
    %DesignConstraintSet Describes a collection of design constraints for
    %use in an optimsitation routine.
    %
    % The definition of the 'awi.fe.opt.DesignConstraintSet' class matches
    % that of the DCONADD bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Identification number of the Design Constraint Set
        DCID
        %Identification of the design constraints
        DC
    end
    
    %Handle to hidden 'awi.fe.opt' objects
    properties (Hidden = true)
        %Handle to the underlying 'awi.fe.opt.DesignConstraint' objects
        DesignConstraints
    end
    
    methods % set / get
        function set.DCID(obj, val)              %set.DCID
            %set.DCID Set method for the property 'DCID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.DC(obj, val)                %set.DC
            %set.DC Set method for the property 'DC'.
            
            validateattributes(val, {'numeric'}, {'row', 'integer'}, ...
                class(obj), 'DC');
            obj.DC = val;
        end
        function set.DesignConstraints(obj, val) %set.DesignConstraints
            %set.DesignConstraints Set method for the property
            %'DesignConstraints'.
            
            if isempty(val) %Always okay
                obj.DesignConstraints = [];
                return
            end
            
            validateattributes(val, {'awi.fe.opt.DesignConstraint'}, ...
                {'row', 'nonempty'}, class(obj), 'DesignConstraints');
            obj.DesignConstraints = val;
        end
        function val = get.DCID(obj)             %get.DCID
            %get.DCID Get method for the property 'DCID'.
            
            val = obj.ID;
        end
        function val = get.DC(obj)               %get.DC
            %get.DC Get method for the property 'DC'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.DesignConstraint' object then always use their ID
            % number, else use DC.
            
            if isempty(obj.DesignConstraints)
                val = obj.DC;
            else
                val = [obj.DesignConstraints.ID];
            end
        end
    end
    
    methods % construction
        function obj = DesignConstraintSet
            
            %Make a note of the property names
            addFEProp(obj, 'DIC', 'DC');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the
            %'awi.fe.opt.DesignConstraintSet' object into a text file using
            %the format of the MSC.Nastran 'DCONADD' bulk data entry.
            %
            % The following assumptions apply:
            %   * 
            
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
                comment = ['DCONADD: Defines the design constraints ', ...
                    'for a subcase as a union of DCONSTR entries.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj  = numel(obj);
            blnks = blanks(8);
            
            %Have to use a loop as we cannot guarantee the same number of
            %DC entries for each object
            for ii = 1 : nObj
                %Constraint IDs
                DCid = obj(ii).DC;
                %No. design contraints
                nDC = numel(DCid);
                if nDC < 7
                    %Gather data and print
                    data = [{'DCONADD'}, {obj.ID}, num2cell(obj.DCID)];
                    fmt  = ['%-8s%-8i', repmat('%-8i', [1, nDC]), '\r\n'];
                    fprintf(fid, fmt, data{:});
                else
                    %Seperate data into the first line and subsequent lines
                    %line
                    row1    = DCid(1 : 7);
                    rowN    = DCid(8 : end);
                    nData   = numel(rowN);
                    nRows   = floor(nData ./ 8);
                    rowData = rowN(1 : nRows * 8);
                    rowData = reshape(rowData, [8, nRows]);
                    remData = rowN(nRows * 8 + 1 : end);
                    %Print the first line
                    fmt  = ['%-8s%-8i', repmat('%-8i', [1, 7]), '\r\n'];
                    data = [{'DCONADD'}, {obj.ID}, num2cell(row1)];
                    fprintf(fid, fmt, data{:});
                    %Print full-lines
                    data = [repmat({blnks}, [1, nRows]) ; num2cell(rowData)];
                    fmt  = ['%-8s', repmat('%-8i', [1, 8]), '\r\n'];
                    fprintf(fid, fmt, data{:});
                    %Print final line
                    data = [{blanks(8)}, num2cell(remData)];
                    fmt  = ['%-8s', repmat('%-8i', [1, numel(remData)]), '\r\n'];
                    fprintf(fid, fmt, data{:});
                end
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
    
end

