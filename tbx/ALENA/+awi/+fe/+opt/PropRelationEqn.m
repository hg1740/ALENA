classdef PropRelationEqn < awi.fe.FEBaseClass
    %PropRelationEqn Defines a relationship between a design variable and
    %an element property using a user-specified equation.
    %
    % The definition of the 'awi.fe.opt.PropRelation' class matches that of
    % the DVPREL2 bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Name of the property entry
        TYPE
        %ID number of the property entry
        PID
        %Name of the property that the is related to the design variable or
        %the ID number of the field 
        PropNameOrField
        %Minimum value of the property
        PMIN
        %Maximum value of the property
        PMAX
        %Constant term of relation
        EQID
        %ID number of the design variables
        DVID
        %Name of any design constants used by the DEQATN
        VarLabel
    end
    
    %Reference to 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.BeamProperty' or 'awi.fe.BeamCrossSection'
        %object
        Property
        %Handle to the 'awi.fe.opt.DesignVariable' objects
        DesignVariables
        %Handle to the 'awi.fe.opt.DesignEquation' object
        DesignEquation
    end
    
    methods % set / get
        function set.TYPE(obj, val)            %set.TYPE
            %set.TYPE Set method for the property 'TYPE'.
            
            validateLabel(obj, val, 'TYPE')
            obj.TYPE = val;
        end
        function set.PID(obj, val)             %set.PID
            %set.PID Set method for the property 'ID'.
            
            %Pass it on
            validateID(obj, val, 'PID');
            obj.PID = val;
            
        end
        function set.PropNameOrField(obj, val) %set.PropNameOrField
            %set.PropNameOrField Set method for the property
            %'PropNameOrField'.
            
            if ischar(val)
                validateLabel(obj, val, 'PropNameOrField')
            elseif isnumeric(val)
                validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
                    class(obj), 'PropNameOrField');
            else
                validateattributes(val, {'numeric', 'char'}, {'scalar'}, ...
                    class(obj), 'PropNameOrField');
            end
            obj.PropNameOrField = val;
        end
        function set.PMIN(obj, val)            %set.PMIN
            %set.PMIN Set method for the property 'PMIN'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'PMIN');
            obj.PMIN = val;
        end
        function set.PMAX(obj, val)            %set.PMAX
            %set.PMAX Set method for the property 'PMAX'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'PMAX');
            obj.PMAX = val;
        end
        function set.EQID(obj, val)            %set.EQID
            %set.EQID Set method for the property 'EQID'.
            
            validateID(obj, val, 'EQID');
            obj.EQID = val;
        end
        function set.DVID(obj, val)            %set.DVID
            %set.DVID Set method for the property 'DVID'.
            
            validateattributes(val, {'numeric'}, {'row', 'integer'}, ...
                class(obj), 'DVID');
            obj.DVID = val;
        end
        function set.VarLabel(obj, val)        %set.VarLabel
           %set.VarLabel Set method for the property 'VarLabel'.
           
           assert(iscellstr(val), ['Expected the ''VarLabel'' to be a', ...
               'cell-array of strings.']);
           obj.VarLabel = val;
        end
        function set.Property(obj, val)        %set.Property
            %set.Property Set method for the property 'Property'.
            %
            % 'Property' must be a handle to an object of class
            % 'awi.fe.BeamProp' or 'awi.fe.BeamCrossSection'
            %
            % - TODO : In the future this will need to chance to a generic
            % FE property class.
            
            validateattributes(val, {'awi.fe.BeamProp', ...
                'awi.fe.BeamCrossSection'}, {'scalar'}, class(obj), 'Property');
            obj.Property = val;
        end
        function set.DesignVariables(obj, val) %set.DesignVariables
            %set.DesignVariabels Set method for the property
            %'DesignVariables'. 
            %
            % 'DesignVariables' must be a row vector of
            % 'awi.fe.optp.DesignVariable' objects.
            
            validateattributes(val, {'awi.fe.opt.DesignVariable'}, ...
                {'vector'}, class(obj), 'DesignVariables');
            if iscolumn(val)
                val = val';
            end
            obj.DesignVariables = val;
        end
        function set.DesignEquation(obj, val)  %set.DesignEquation
            %set.DesignEquation Set method for the property
            %'DesignEquation'.
            
            validateattributes(val, {'awi.fe.opt.DesignEquation'}, { ... 
                'scalar', 'nonempty'}, class(obj), 'DesignEquation');
            obj.DesignEquation = val;
        end
        function val = get.PID(obj)            %get.PID 
           %get.PID Get method for the property 'PID'.
           %
           % If the object has been assigned a handle to its
           % underlying property object then always use their ID number,
           % else use PID.
           
           if isempty(obj.Property)
               val = obj.PID;
           else
               val = obj.Property.ID;
           end
        end
        function val = get.DVID(obj)           %get.DVID
           %get.DVID Get method for the property 'DVID'.
           %
           % If the object has been assigned a handle to its
           % 'awi.fe.opt.DesignVariable' objects then always use their ID number,
           % else use DVID.
           
           if isempty(obj.DesignVariables)
               val = obj.DVID;
           else
               val = [obj.DesignVariables.ID];
           end
           
        end
        function val = get.EQID(obj)           %get.EQID)
            %get.EQID Get method for the property 'DVID'.
           %
           % If the object has been assigned a handle to its
           % 'awi.fe.opt.DesignEquation' objects then always use their ID
           % number, else use EQID.
           
           if isempty(obj.DesignEquation)
               val = obj.DVID;
           else
               val = obj.DesignEquation.ID;
           end
            
        end
    end
    
    methods % construction
        function obj = PropRelationEqn
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'TYPE', 'PID', 'PropNameOrField', ...
                'EQID', 'DVID', 'VarLabel');
            
            %Update the card name
            obj.CardName = 'DVPREL2';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.opt.PropRelation'
            %object into a text file using the format of the MSC.Nastran
            %'DVPREL1' bulk data entry.
                        
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
                comment = ['DVPREL2: Defines the relation between an ', ...
                    'analysis model property and design variables '   , ...
                    'with a user-supplied equation.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj = numel(obj);
                 
            %Can only vectorise if every object has a single design
            %variable attached
            nDV = arrayfun(@(o) numel(o.DesignVariables), obj);
            nDT = arrayfun(@(o) numel(o.VarLabel)       , obj);
            
            %Write the data to the file 
            if and(range(nDV) == 0, range(nDT) == 0) && nDV(1) < 7
                %Vectorised...
                
                %Card name
                nam    = repmat({'DVPREL2'}, [1, nObj]);
                dv     = repmat({'DESVAR'} , [1, nObj]);
                tbl    = repmat({'DTABLE'} , [1, nObj]);
                blnks  = repmat({blanks(8)}, [1, nObj]); 
                
                %Grab data - Convert PMIN/PMAX to strings to vectorise
                pid  = get([obj.Property], {'ID'})';
                pnam = {obj.PropNameOrField};
                idx  = cellfun(@(x) isnumeric(x), pnam);
                pnam(idx) = cellfun(@(x) num2str(X), pnam(idx), 'Unif', false); 
                pmin = cellfun(@num2str, {obj.PMIN}, 'Unif',false);
                pmin(cellfun(@isempty, pmin)) = {blanks(8)};
                pmax = cellfun(@num2str, {obj.PMAX}, 'Unif',false);
                pmax(cellfun(@isempty, pmax)) = {blanks(8)};
                eqid = get([obj.DesignEquation], {'ID'})';
                
                %Get DESVAR ID & DTABLE labels
                dvid = arrayfun(@(o) sprintf('%-8i', [o.DesignVariables.ID]), obj, 'Unif', false);
                lab  = arrayfun(@(o) sprintf('%-8s', o.VarLabel{:})         , obj, 'Unif', false);
                
                %Set up the format for printing
                data = [ ...
                    nam   ; {obj.ID} ; {obj.TYPE} ; pid ; pnam ; pmin ; pmax ; eqid ; ...
                    blnks ; dv       ; dvid       ; ...
                    blnks ; tbl      ; lab];          
                
                %Write in 16-character column width as standard
                format = ['%-8s%-8i%-8s%-8i%-8s%-8s%-8s%#-8i\r\n', ...
                    '%-8s%-8s%-s\r\n', ...
                    '%-8s%-8s%-s\r\n'];
                
                %Write the data to the file
                fprintf(fid, format, data{:});

            else
                %Use loop...
                for iO = 1 : nObj
                    
                    nDV = numel(obj(iO).DesignVariables);
                    nDT = numel(obj(iO).VarLabel);
                    
                    if nDV > 7 || nDT > 7
                        error('Update code');
                    end
                    
                    pmin = num2str(obj(iO).PMIN);
                    if isempty(pmin), pmin = blanks(8); end
                    pmax = num2str(obj(iO).PMAX);
                    if isempty(pmax), pmax = blanks(8); end
                    vl = obj(iO).VarLabel;
                    
                    fmt_dv = ['%-8s%-8s', repmat('%-8i', [1, nDV]), '\r\n'];
                    fmt_dt = ['%-8s%-8s', repmat('%-8s', [1, nDT]), '\r\n'];
                    
                    fprintf(fid, '%-8s%-8i%-8s%-8i%-8s%-8s%-8s%#-8i\r\n', ...
                        'DVPREL2', obj(iO).ID, obj(iO).TYPE, obj(iO).Property.ID, ...
                        obj(iO).PropNameOrField, pmin, pmax, obj(iO).DesignEquation.ID);
                    fprintf(fid, fmt_dv, blanks(8), 'DESVAR', obj(iO).DesignVariables.ID);
                    fprintf(fid, fmt_dt, blanks(8), 'DTABLE', vl{:});
                end
            end

            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

