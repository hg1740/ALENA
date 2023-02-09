classdef DesignResponse < awi.fe.FEBaseClass
    %DesignResponse Describes an output quantity of the analysis that can
    %be used to formulate a constraint or objective function.
    %
    % The definition of the 'awi.fe.opt.DesignResponse' class matches that
    % of the DRESP1 bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %String for identifying this design response
        LABEL
        %String describing the type of response
        ResponseType
        %String describing the propery name or the keywork "ELEM"
        PropertyType
        %Region identifer for constraint screening
        Region
        %Response attributes
        ResponseAttributes
    end
    
    methods % set / get
        function set.LABEL(obj, val)              %set.LABEL
            %set.LABEL Set method for the property 'LABEL'.
            
            validateLabel(obj, val, 'LABEL')
            obj.LABEL = val;
        end
        function set.ResponseType(obj, val)       %set.ResponseType
            %set.ResponseType Set method for the property 'ResponseType'.
            
            validateLabel(obj, val, 'ResponseType')
            obj.ResponseType = val;
        end
        function set.PropertyType(obj, val)       %set.PropertyType
            %set.PropertyType Set method for the property 'PropertyType'.
            
            validateLabel(obj, val, 'PropertyType')
            obj.PropertyType = val;
        end
        function set.Region(obj, val)             %set.Region
            %set.Region Set method for the property 'Region'
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
                class(obj), 'Region');
            obj.Region = val;
        end
        function set.ResponseAttributes(obj, val) %set.ResponseAttributes
            %set.ResponseAttributes Set method for the property
            %'ResponseAttributes'
            
            validateattributes(val, {'cell'}, {'row', 'nonempty'}, ... 
                class(obj), 'ResponseAttributes');
            obj.ResponseAttributes = val;
        end
    end
    
    methods % construction
        function obj = DesignResponse
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'LABEL', 'ResponseType', 'PropertyType', 'Region', 'ResponseAttributes');
            
            %Update the card name
            obj.CardName = 'DRESP1';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.opt.DesignResponse'
            %object into a text file using the format of the MSC.Nastran
            %'DRESP1' bulk data entry.
            %
            % The following assumptions are made:
            %   * The 'REGION' entry for the DRESP1 is left blank.
            
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
                comment = ['DRESP1: Defines a set of structural ' , ...
                    'responses that is used in the design either ', ...
                    'as constraints or as an objective.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                      
            %Can only vectorise if every object has a single ATTi 
            attr = vertcat(obj.ResponseAttributes);
            nAtr = cellfun(@numel, attr(:, 3));
            
            %Write the data to the file
            idx = nAtr > 1;
            if any(idx)
                
                %How many and where are the responses with > 1 attribute?
                n   = nnz(idx);
                ind = find(idx);
                
                %Loop through and write to file
                for i = 1 : n
                    
                    nAtr_ = nAtr(ind(i));
                    o     = obj(ind(i));
                    attr_ = attr(ind(i), :);
                    
                    nLines = ceil((nAtr_ - 1) / 8);
                    
                    if nLines > 1
                        error('Update code for arbitrary number of lines');
                    end
                    
                    %Line 1          
                    fmt = '%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s\r\n';
                    dat = [{'DRESP1', o.ID, o.LABEL, o.ResponseType, o.PropertyType, blanks(8)}, attr_(1 : 2), attr_{3}(1)];
                    dat = cellfun(@num2str, dat, 'Unif', false);
                    dat(cellfun(@isempty, dat)) = {blanks(8)};
                    fprintf(fid, fmt, dat{:});
                    
                    %Continuation Lines
                    dat = num2cell(attr_{3}(2 : end));
                    fmt = ['%-8s', repmat('%-8i', [1, numel(dat)]), '\r\n'];
                    fprintf(fid, fmt, blanks(8), dat{:});
                    
                end
                
                %Vectorise all the others
                attr       = attr(~idx, :);              
                i_vectorised(fid, obj(~idx), attr);
                
            else
                
                %Simple vectorised writing function for the case where only
                %one attribute is defined.
                i_vectorised(fid, obj, attr);
                
            end

            if bClose %Close the file?
                fclose(fid);
            end
            
            function i_vectorised(fid, obj_, attr)
                %i_vectorised Vecrorised writing of the data to a file
                
                nObj = numel(obj_);
                
                %Convert attribute data to string format
                attr(:, 1) = cellfun(@num2str, attr(:, 1), 'Unif', false);
                attr(:, 2) = cellfun(@num2str, attr(:, 2), 'Unif', false);
                attr(:, 3) = cellfun(@num2str, attr(:, 3), 'Unif', false);
                
                %Card name
                nam    = repmat({'DRESP1'}  , [1, nObj]);
                blnks  = repmat({ blanks(8)}, [1, nObj]);           
                         
                %Check 'attr' for empties
                attr(cellfun(@isempty, attr(:, 2)), 2) = {blanks(8)};
                
                %Set up the format for printing
                data = [ ...
                    nam   ; {obj_.ID} ; {obj_.LABEL} ; {obj_.ResponseType} ; {obj_.PropertyType} ; blnks ; attr'];          
                
                %Write in 8-character column width as standard
                format = '%-8s%-8i%-8s%-8s%-8s%-8s%-8s%-8s%-8s\r\n';
                
                %Write the data to the file
                fprintf(fid, format, data{:});
                
            end
            
        end
    end
    
end

