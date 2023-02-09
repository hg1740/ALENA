classdef DesignVariable < awi.fe.FEBaseClass
    %DesignVariable Describes a design variable to be used in a generic
    %optimisation problem.
    %
    % The definition of the 'awi.fe.opt.DesignVariable' class matches that
    % of the DESVAR bulk data type from MSC.Nastran.

    %Primary Properties
    properties
        %String for indentifying this design variable 
        LABEL
        %Initial value
        XINIT
        %Lower bound
        XLB = -1e20;
        %Upper bound
        XUB = 1e20;
        %Fractional change allowed for the design variable during
        %approximate optimization
        DELXV
        %Identification number of a 'awi.
        %DDVAL
    end
    
    methods % set / get
        function set.LABEL(obj, val) %set.LABEL 
            %set.LABEL Set method for the property 'LABEL'.
            
            validateLabel(obj, val, 'LABEL')
            obj.LABEL = val;
        end
        function set.XINIT(obj, val) %set.XINIT 
            %set.XINIT Set method for the propery 'XINIT'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'XINIT');
            obj.XINIT = val;
        end
        function set.XLB(obj, val)   %set.XLB   
            %set.XLB Set method for the property 'XLB'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'XLB');
            obj.XLB = val;
        end
        function set.XUB(obj, val)   %set.XUB   
            %set.XUB Set method for the property 'XUB'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'XUB');
            obj.XUB = val;
        end
        function set.DELXV(obj, val) %set.DELXV 
           %set.DELXV Set method for the property 'DELXV'.
           
           validateattributes(val, {'numeric'}, {'positive', '<=', 1}, ...
               class(obj), 'DELXV');
           obj.DELXV = val;
        end
    end
    
    methods % construction
        function obj = DesignVariable
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'LABEL', 'XINIT', 'XLB', 'XUB', 'DELXV'); %, 'DDVAL');
            
            %Update the card name
            obj.CardName = 'DESVAR';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.opt.DesignVariable'
            %object into a text file using the format of the MSC.Nastran
            %'DESVAR' bulk data entry.
            %
            % The following assumptions are made:
            %   * The design variable is continuous i.e. DDVAL is blank
            %   * The entry is written in large-field format to allow the
            %     maximum amount of data to be written.
            
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
                comment = ['DESVAR: Defines a design variable for ', ...
                    'design optimization.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
                 
            nObj = numel(obj);
            
            %Card name
            nam    = repmat({'DESVAR*'}, [1, nObj]);
            blnks  = repmat({['*', blanks(7)]}, [1, nObj]);
            blnks_ = repmat({blanks(16)}, [1, nObj]);
            
            %Check the user has defined the initial conditions
            x0 = {obj.XINIT};
            if any(cellfun(@isempty, x0))
                error('Initial conditions for design variables not set');
            end
            
            %Check that the user has defined the move constraints (DELXV)
            dxv = {obj.DELXV};
            idx = cellfun(@isempty, dxv);
            dxv(idx) = {blanks(16)};
            dxv(~idx) = cellfun(@(x) num2str(x, '%#-16.8g'), dxv(~idx), 'Unif', false);
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.LABEL} ; x0 ; {obj.XLB} ; ...
                blnks ; {obj.XUB} ; dxv ; blnks_];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%-16s%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%-16s%-16s\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

