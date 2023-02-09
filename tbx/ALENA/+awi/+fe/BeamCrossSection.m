classdef BeamCrossSection < awi.fe.FEBaseClass
    %BeamCrossSection Describes the beam properties of a flexible beam for
    %use in a finite element model by defining the cross-section shape and
    %its dimensions.
    %
    % The definition of the 'BeamCrossSecion' object matches that of the
    % PBEAML bulk data type from MSC.Nastran.
    
    %Primary properites
    properties
        %Property identification number
        PID
        %Material identification number
        MID
        %Name of the cross-section
        TYPE
        %Value of the cross-section dimensions at ends A & B
        Dimensions
        %Non-structural mass at ends A & B
        NSM
    end
    
    %Store a handle to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Material' object
        Material     
    end
    
    properties (Constant, Hidden = true)
        %Valid cross-section names
        ValidCrossSectionNames = {'BOX'};
    end
    
    methods % set / get
        function set.PID(obj, val)        %set.PID        
            %set.PID Set method for the property 'PID'.
            % 
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.MID(obj, val)        %set.MID        
            validateID(obj, val, 'MID')
            obj.MID = val;
        end
        function set.TYPE(obj, val)       %set.Type       
           %set.Type Set method for the property 'Type'. 
           
           validatestring(val, obj.ValidCrossSectionNames, class(obj), 'TYPE');
           obj.TYPE = val;
        end
        function set.Dimensions(obj, val) %set.Dimensions 
           %set.Dimensions Set method for the property 'Dimensions'.
           
           validateattributes(val, {'numeric'}, {'2d', 'nrows', 2, ...
               'positive'}, class(obj), 'Dimensions');
           obj.Dimensions = val;
        end
        function set.NSM(obj, val)        %set.NSM        
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'NSM');
            obj.NSM = val;
        end
        function set.Material(obj, val)   %set.Material   
            %set.Material Set method for the property 'Material'.
            %
            % 'Material' must be a scalar 'awi.fe.Material' object.
            
            validateattributes(val, {'awi.fe.Material'}, {'scalar', ...
                'nonempty'}, class(obj), 'Material');
            obj.Material = val;
        end
        function val = get.PID(obj)       %get.PID        
           %get.PID Get method for the property 'PID'.
           val = obj.ID;
        end
        function val = get.MID(obj)       %get.MID        
            %get.MID Get method for the property 'MID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.Material' object then always use their ID number, 
            % else use MID.
            
            if isempty(obj.Material)
                val = obj.MID;
            else
                val = obj.Material.ID;
            end
                
        end
    end
    
    methods % construction
        function obj = BeamCrossSection
            
            %Make a note of the property names
            addFEProp(obj, 'PID', 'MID', 'TYPE', 'Dimensions', 'NSM');
            
        end
    end
    
    methods % writing to a .txt file
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.BeamCrossSection'
            %object into a text file using the format of the MSC.Nastran
            %'PBEAML' bulk data entry.
            %
            % The following assumptions are made:
            %   * Each beam only has 2 stations where the properties are
            %     calculated, ends A and ends B. Not intermediate sections
            %     are defined.
                        
            %Check that all object are of the same type - i.e. the same
            %number of cross-sectional dimensions
            type = {obj.TYPE};
            assert(numel(unique(type)) == 1, ['Expected the ', ...
                '''BeamCrossSection'' objects to all be of the same ', ...
                'type. Unable to continue.']);
            
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
                comment = ['PBEAML : Defines the properties of a beam ', ...
                    'element by cross-sectional dimensions.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                 
            nObj = numel(obj);
            nDim = size(obj(1).Dimensions, 2);
            
            %Grab the beam properties
            NSM_ = num2cell([obj.NSM]);
                        
            %Card name, blanks, zero/ones arrays & stress recovery flag  
            nam   = repmat({'PBEAML'} , [1, nObj]);
            blnks = repmat({blanks(8)}, [1, nObj]);
            zrs   = num2cell(zeros(1, nObj));
            onz   = num2cell(ones(1, nObj));
            SO    = repmat({'YES'}, [1, nObj]);
            
            %Set up the cross-section dimensions & NSM /SO data
            %   - Assume we don't know how many we have
            dimA = arrayfun(@(i) obj(i).Dimensions(1, :)', 1 : nObj, 'Unif', false);
            dimB = arrayfun(@(i) obj(i).Dimensions(2, :)', 1 : nObj, 'Unif', false);
            dimA = num2cell(horzcat(dimA{:}));
            dimB = num2cell(horzcat(dimB{:}));
            dim  = [dimA ; NSM_(1, :) ; SO ; onz ; dimB ; NSM_(2, :)];
            
            %Define the format-spec for 'fprintf'
            fStr = {'%#-8.5f'};
            fmt  = [repmat(fStr, [1, nDim]), fStr, {'%-8s'}, fStr, repmat(fStr, [1, nDim]), fStr];                
            
            %Can have 8 columns per row
            nCol  = 8;
            nData = size(dim, 1);
            nRows = floor(nData / nCol); %Number of complete rows
            rem   = mod(nData, nCol);    %Number of remaining data
                        
            rowData = dim(1 : nRows * nCol, :);
            remData = dim(end - (rem - 1) : end, :);     
            fmtData = fmt(1 : nRows * nCol);
            fmtRem  = fmt(end - (rem - 1) : end);
            
            %Pad with blanks or end-of-line characters
            ub = (nCol : nCol : nCol * nRows);
            lb = cumsum(ub) - nCol + 1;
            rowData = arrayfun(@(i) rowData(lb(i) : ub(i), :), 1 : nRows, 'Unif', false);
            fmtData = arrayfun(@(i) fmtData(lb(i) : ub(i))   , 1 : nRows, 'Unif', false);
            rowData = [rowData ; repmat({blnks}, [1, nRows])];
            rowData = vertcat(rowData{:});
            fmtData = [fmtData ; repmat({'\r\n%-8s'}, [1, nRows])];
            fmtData = horzcat(fmtData{:});
            
            %Define data and format for printing
            data   = [ ...
                nam     ; {obj.ID} ; {obj.MID} ; blnks ; {obj.TYPE} ; ...
                blnks   ; blnks    ; blnks     ; blnks ; blnks      ; ...
                rowData ; remData];
            format = ['%-8s%-8i%-8i%-8s%-8s%-8s%-8s%-8s%-8s\r\n%-8s', ...
                horzcat(fmtData{:}), horzcat(fmtRem{:}), '\r\n'];            
          
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

