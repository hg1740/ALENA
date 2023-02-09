classdef DesignEquation < awi.fe.FEBaseClass
    %DesignEquation Describes an equation used as part of an optimisation.
    %
    % The definition of the 'awi.fe.opt.DesignEquation' class matches
    % that of the DEQATN bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Equation text 
        Eqn
    end
    
    methods % set / get
        function set.Eqn(obj, val) %set.Eqn
            %set.Eqn Set method for the property 'Eqn'.
            
            validateattributes(val, {'char'}, {'row', 'nonempty'}, ...
                class(obj), 'Eqn');
            obj.Eqn = val;
            
        end
    end
    
    methods % construction
        function obj = DesignEquation
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'Eqn');
            
            %Update the card name
            obj.CardName = 'DEQATN';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.opt.DesignEquation'
            %object into a text file using the format of the MSC.Nastran
            %'DEQATN' bulk data entry.
            
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
                comment = ['DEQATN: Defines a table of real constants ', ...
                    'that are used in equations.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj  = numel(obj);
            
            %Have to use a loop as each DEQATN could have arbitrary length
            for ii = 1 : nObj
                %How many characters in the equation?
                nChar = numel(obj(ii).Eqn);
                if nChar > 56
                    %56 characters on the first line, 72 afterwards
                    nLines = floor((nChar - 56)/64) + 1;
                    %Number of characters on subsequent lines
                    nChar_ = nChar - 56;
                    if nChar_ < 64
                        nChar_per_line = [56, nChar_];
                    else
                        rem_ = nChar - 56 - (64*(nLines - 1));
                        %How many characters on each line?
                        nChar_per_line = [56, repmat(64, [1, nLines - 1]), rem_];
                    end
                else
                    %Simple...
                    nLines = 1;
                    %How many characters on each line?
                    nChar_per_line = nChar;
                end

                %Calculate bounds
                ub = cumsum(nChar_per_line);
                lb = [1, ub(1 : end - 1) + 1];
                
                %Split into cell-arrays - One per line
                data = arrayfun(@(i) obj(ii).Eqn(lb(i) : ub(i)), 1 : numel(nChar_per_line), 'Unif', false);
                
                %Add preamble to the first line
                data{1} = sprintf('%-8i%s', obj(ii).ID, data{1});
                
                %Format data for printing
                data = [ ...
                    ['DEQATN', repmat({blanks(8)}, [1, numel(data) - 1])] ; ...
                    data]; %#ok<AGROW>
%                 data = [data ; repmat({blanks(8)}, size(data))];
                
                %Print
                fprintf(fid, '%-8s%-s\r\n', data{:});
                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

