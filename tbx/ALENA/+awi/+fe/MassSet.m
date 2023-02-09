classdef MassSet < awi.fe.FEBaseClass
    %MassSet Describes a collection of mass groups.
    %
    % The definition of the 'MassSet' object matches that of the 'MASSSET'
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Global scale factor
        S0 = 1.
        %Scale factors for individual mass groups
        Si = 1;
        %ID number of mass groups in this 'MassSet'
        MassGroupID
    end
    
    %Flag to include base mass
    properties
        %Logical flag to include base mass in this mass set
        IncludeBaseMass = true;
    end
    
    %Handle to 'awi.fe' objects
    properties (Hidden)
        %Handle to the 'awi.fe.MassModel' object(s)
        MassModel
    end
    
    methods % set / get
        function set.S0(obj, val)              %set.S0
            %set.S0 Set method for the property 'S0'
            
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'nonnan', 'real', 'nonnegative'}, class(obj), 'S0');
            obj.S0 = val;
            
        end
        function set.Si(obj, val)              %set.Si
            %set.Si Set method for the property 'Si'.
            
            validateattributes(val, {'numeric'}, {'column', 'finite', ...
                'nonnan', 'real', 'nonenegative'}, class(obj), 'S0');
            obj.Si = val;
            
        end
        function set.MassGroupID(obj, val)     %set.MassGroupID 
            %set.MassGroupID Set method for the property 'MassGroupID'.
            
            validateattributes(val, {'numeric'}, {'integer', 'column', ...
                'nonnan', 'finite', 'real'}, class(obj), 'MassGroupID');
            obj.MassGroupID = val;
            
        end
        function set.IncludeBaseMass(obj, val) %set.IncludeBaseMass 
           %set.IncludeBaseMass Set method for the property
           %'IncludeBaseMass'.
           
           validateattributes(val, {'logical'}, {'scalar'}, class(obj), ...
               'IncludeBaseMass');
           obj.IncludeBaseMass = val;
           
        end
        function set.MassModel(obj, val)       %set.MassModel   
            %set.MassModel Set method for the property 'MassModel'.
            
            validateattributes(val, {'awi.fe.MassModel'}, {'column'}, ...
                class(obj), 'MassModel');
            obj.MassModel = val;
            
        end
        function val = get.MassGroupID(obj)    %get.MassGroupID 
           %get.MassGroupID Get method for the property 'MassGroupID'.
           
           if isempty(obj.MassModel)
               val = obj.MassGroupID;
           else
               val = vertcat(obj.MassModel.ID);
           end
           
        end
    end
    
    methods % construction
        function obj = MassSet
            
            %Make a note of the property names
            addFEProp(obj, 'S0', 'Si', 'MassGroupID');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.MassSet' object
            %into a text file using the format of the MSC.Nastran 'MASSSET'
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
                comment = ['MASSSET : Defines a linear combination ', ...
                    'of mass cases to form the subcase dependent mass..'];
                awi.fe.FEBaseClass.writeComment(comment, fid);
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            %Write each object seperately as MASSSET is list-formatted
            for i = 1 : numel(obj)
                
                %Grab data 
                id  = obj(i).ID; 
                s0  = obj(i).S0;
                si  = [obj(i).Si]';
                mid = [obj(i).MassGroupID]';
                
                %Include base mass?
                if obj(i).IncludeBaseMass
                    si  = [si , 1]; %#ok<*AGROW>
                    mid = [mid , 0];
                end
                
                awi.fe.FEBaseClass.writeSuperpositionList(fid, 'MASSSET', id, s0, si, mid);
                
%                 if numel(si) ~= numel(mid) %Escape route
%                     continue
%                 end
%                 
%                 %Print the data
%                 i_printBulkData(fid, id, s0, si, mid);
                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
            function i_printBulkData(fid, id, s0, si, mid)
                %i_printBulkData Prints the bulk data for the 'MASSSET'
                %card.
                
                %Set up data for vectorised writing
                blnks = {blanks(8)};
                nData = numel(mid);
                if nData > 2
                    line1 = mid(1 : 3);
                    Si1   = si(1 : 3);
                else
                    line1 = mid(1 : end);
                    Si1   = si(1 : end);
                end
                rem    = nData - 3;
                nLines = ceil(rem / 4);
                if nLines > 1
                    ub       = 3 + (nLines - 1) * 4;
                    lineData = mid(4 : ub);
                    SiData   = si(4 : ub);
                    lineEnd  = mid(ub + 1 : end);
                    SiEnd    = si(ub + 1 : end);
                    ub    = cumsum(repmat(8, [1, nLines - 1]));
                    lb    = [1, ub(1 : end - 1) + 1];
                else
                    lineData = [];
                    SiData   = [];
                    lineEnd  = mid(4 : end);
                    SiEnd    = si(4 : end);
                    ub = numel(lineData) * 2;
                    lb = 1;
                end
                
                %First line
                data  = [id ; s0];
                data1 = [Si1 ; line1];
                data  = [{'MASSSET'} ; num2cell(data) ; num2cell(data1(:))];
                
                %Intermediate lines - pad with blanks
                dataN = [SiData ; lineData];
                dataN = dataN(:);
                dataN = arrayfun(@(i) num2cell(dataN(lb(i) : ub(i))), 1 : nLines - 1, 'Unif', false);
                dataN = [dataN ; repmat(blnks, [1, nLines - 1])];
                dataN = vertcat(dataN{:});
                
                %Final line
                dataEnd = [SiEnd ; lineEnd];
                dataEnd = num2cell(dataEnd(:));
                
                %Format for printing
                data = [data ; blnks ; dataN ; dataEnd ];
                fmt  = [ ...
                    '%-8s%-8i%#-8.4g', ...
                    repmat('%#-8.4g%-8i', [1, numel(line1)]), '\r\n', ...
                    repmat('%-8s%#-8.4g%-8i%#-8.4g%-8i%#-8.4g%-8i%#-8.4g%-8i\r\n', [1, nLines - 1]), ...
                    '%-8s', repmat('%#-8.4g%-8i', [1, numel(lineEnd)]), '\r\n'];
                
                %Print it
                fprintf(fid, fmt, data{:});
                
            end
            
        end
    end
    
end

