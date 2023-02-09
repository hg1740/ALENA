classdef DesignTable < awi.fe.FEBaseClass
    %DesignTable Describes a set of constants that are used in design
    %equations.
    %
    % The definition of the 'awi.fe.opt.DesigntABLE' class matches that of
    % the DTABLE bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Property names
        Prop = {};
        %Value of the properties
        Val  = {};
    end
    
    methods % set / get
        function set.Prop(obj, val) %set.Prop
            %set.Prop Set method for the property 'Prop'.
            
            assert(iscellstr(val), ['Expected ''Prop'' to be a cell ', ...
                'array of strings']);
            idx = cellfun(@(x) numel(x) < 9, val);
            assert(all(idx), ['Expected each property name in ', ...
                '''Prop'' to be less than 8 characters.']);
            if iscolumn(val)
                val = val';
            end
            obj.Prop = val;
            
        end
        function set.Val(obj, val)  %set.Val
            %set.Val Set method for the property 'Val'.
            
            if isnumeric(val)
                validateattributes(val, {'double'}, {'row', 'nonempty'}, ...
                    class(obj), 'Val');
                val = cellstr(num2str(val', '%-#.3g'))';
            elseif iscellstr(val)
                %Assume the user is defining exactly what they want in the
                %DTABLE.
            else
                error('Expected ''Val'' to be a cell-array of strings or a numeric row vector');
            end
            obj.Val = val;
            
        end
    end
    
    methods % construction
        function obj = DesignTable
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'Prop', 'Val');
            
            %Update the card name
            obj.CardName = 'DTABLE';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.opt.DesignTable'
            %object into a text file using the format of the MSC.Nastran
            %'DTABLE' bulk data entry.
            
            %Check an equal number of properties and vlaues have been
            %defined.
            idx = arrayfun(@(o) numel(o.Prop) == numel(o.Val), obj);
            assert(all(idx), ['Expected the number of values to equal ', ...
                'the number of properties on the DesignTable entries']);
            
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
                comment = ['DTABLE: Defines a table of real constants ', ...
                    'that are used in equations.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj  = numel(obj);
            
            %Have to use a loop as each DTABLE could have an arbitrary
            %number of entries
            for ii = 1 : nObj
                nPrp = numel(obj(ii).Prop);
                nCol = 4;
                nRow = floor(nPrp ./ nCol);
                if nRow == 0
                    %Convert values to strings;
                    data = [obj(ii).Prop ; obj(ii).Val];
%                     fprintf(fid, '%-8s%-8s\r\n', 'DTABLE', strjoin(data(:), ','));                 
                    fprintf(fid, ['%-8s', repmat('%-8s%-8s', [1, nPrp]), '\r\n'], ...
                        'DTABLE', data{:});
                else
                    %Print the complete rows
                    n      = nRow * 4;
                    prpStr = obj(ii).Prop(1 : n);
                    valStr = obj(ii).Val(1 : n);
                    data   = [prpStr ; valStr];
                    data   = reshape(data(:), [8, nRow]);
                    data   = [ ...
                        [{'DTABLE'}, repmat({blanks(8)}, [1, nRow - 1])] ; ...
                        data]; %#ok<AGROW>
                    fmt    = ['%-8s', repmat('%-8s%-8s', [1, 4]), '\r\n'];
                    fprintf(fid, fmt, data{:});
                    %Deal with the remainder
                    rem = nPrp - (nRow * 4);
                    if rem ~= 0
                        prpStr = obj(ii).Prop(end - rem + 1 : end);
                        valStr = obj(ii).Val(end - rem + 1 : end);
                        data   = [prpStr ; valStr];
                        data   = data(:);
                        fmt    = ['%-8s', repmat('%-8s%-8s', [1, numel(prpStr)]), '\r\n'];
                        fprintf(fid, fmt, blanks(8), data{:});
                    end
                end
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

