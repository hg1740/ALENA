classdef DesignConstraint < awi.fe.FEBaseClass
    %DesignConstraint Defines a constraint to be used in an optimisation.
    %
    % The definition of the 'awi.fe.opt.DesignConstraint' class matches
    % that of the DCONSTR bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Identification number of the design constraint
        DCID
        %Identification of the design response
        RID
        %Lower bound of the response quantity
        LowerBound = -1e20;
        %Upper bound of the response quantity
        UpperBound = 1e20;
    end
    
    %Store a reference to the 'awi.fe.opt' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.DesignResponse' object that this 'Beam'
        %connects to.
        DesignResponse
    end
    
    methods % set / get
        function set.DCID(obj, val)           %set.DCID
            %set.DCID Set method for the property 'DCID'.
            
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.RID(obj, val)            %set.RID
            %set.RID Set method for the property 'RID'.
            
            validateID(obj, val, 'RID');
            obj.RID = val;
        end
        function set.LowerBound(obj, val)     %set.LowerBound
            %set.LowerBound Set method for the property 'LowerBound'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, class(obj), 'LowerBound');
            obj.LowerBound = val;
        end
        function set.UpperBound(obj, val)     %set.UpperBound
            %set.UpperBound Set method for the property 'UpperBound'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, class(obj), 'UpperBound');
            obj.UpperBound = val;
        end
        function set.DesignResponse(obj, val) %set.DesignResponse
            %set.DesignResponse Set method for the property
            %'DesignResponse'.
            
            validateattributes(val, {'awi.fe.opt.DesignResponse', ...
                'awi.fe.opt.DesignResponseEqn'}, {'scalar'}, class(obj), ...
                'DesignResponse');
            obj.DesignResponse = val;
        end
        function val = get.DCID(obj)          %get.DCID
            %get.DCID Get method for the property 'DCID'.
            
            val = obj.ID;
        end
        function val = get.RID(obj)           %get.RID
            %get.RID Get method for the property 'RID'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.opt.DesignResponse' object then always use their ID
            % number, else use RID.
            
            if isempty(obj.DesignResponse)
                val = obj.RID;
            else
                val = obj.DesignResponse.ID;
            end
        end
    end
    
    methods % construction
        function obj = DesignConstraint
            
            %Make a note of the property names
            addFEProp(obj, 'DCID', 'RID', 'LowerBound', 'UpperBound');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the
            %'awi.fe.opt.DesignConstraint' object into a text file using
            %the format of the MSC.Nastran 'DCONSTR' bulk data entry.
            %
            % The following assumptions are made:
            %   * The 'LOWFQ', 'HIGHFQ' values are ignored.
            
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
                comment = ['DCONSTR: Defines design constraints.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                 
            nObj = numel(obj);
            
            %Card name
            nam    = repmat({'DCONSTR'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.RID} ; {obj.LowerBound} ; {obj.UpperBound}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-8i%-8i%#-8.2g%#-8.2g\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

