classdef PropRelation < awi.fe.FEBaseClass
    %PropRelation Defines a relationship between a design variable and an
    %element property.
    %
    % The definition of the 'awi.fe.opt.PropRelation' class matches that of
    % the DVPREL1 bulk data type from MSC.Nastran.
    
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
        C0 = 0;
        %ID number of the design variables
        DVID
        %Coefficients for each design variable
        COEFi
    end
    
    %Reference to 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.BeamProperty' or 'awi.fe.BeamCrossSection'
        %object
        Property
        %Handle to the 'awi.fe.opt.DesignVariable' objects
        DesignVariables
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
        function set.C0(obj, val)              %set.C0
            %set.C0 Set method for the property 'C0'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty'}, ...
                class(obj), 'C0');
            obj.C0 = val;
        end
        function set.DVID(obj, val)            %set.DVID
            %set.DVID Set method for the property 'DVID'.
            
            validateattributes(val, {'numeric'}, {'row', 'integer'}, ...
                class(obj), 'DVID');
            obj.DVID = val;
        end
        function set.COEFi(obj, val)           %set.COEFi
           %set.COEFi Set method for the property 'COEFi'.
           
           validateattributes(val, {'numeric'}, {'row', 'real', ...
               'nonnan', 'finite'}, class(obj), 'COEFi');
           obj.COEFi = val;
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
                'awi.fe.BeamCrossSection', 'awi.fe.DamperProp', ...
                'awi.fe.SpringProp', 'awi.fe.ScalarMassProp'} , ...
                {'scalar'}, class(obj), 'Property');
            obj.Property = val;
        end
        function set.DesignVariables(obj, val) %set.DesignVariables
            %set.DesignVariabels Set method for the property
            %'DesignVariables'. 
            %
            % 'DesignVariables' must be a row vector of
            % 'awi.fe.optp.DesignVariable' objects.
            
            validateattributes(val, {'awi.fe.opt.DesignVariable'}, ...
                {'row'}, class(obj), 'DesignVariables');
            obj.DesignVariables = val;
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
    end
    
    methods % construction
        function obj = PropRelation
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'LABEL', 'XINIT', 'XLB', 'XUB', 'DELXV'); %, 'DDVAL');
            
            %Update the card name
            obj.CardName = 'DVPREL1';
            
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
                comment = ['DVPREL1: Defines the relation between an ', ...
                    'analysis model property and design variables.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj = numel(obj);
                 
            %Can only vectorise if every object has a single design
            %variable attached
            nDV = arrayfun(@(o) numel(o.DesignVariables), obj);
            
            %Write the data to the file
            if any(nDV > 1)
                %Simple loop
                error('Update the code for case where the DVPREL1 references multiple DESVAR');
            else
                
                %Card name
                nam    = repmat({'DVPREL1'} , [1, nObj]);
                blnks  = repmat({ blanks(8)}, [1, nObj]); 
                
                %Grab data - Convert PMIN/PMAX/Ci to strings to vectorise
                pid  = get([obj.Property], 'ID')';
                pnam = {obj.PropNameOrField};
                idx  = cellfun(@(x) isnumeric(x), pnam);
                pnam(idx) = cellfun(@(x) num2str(X), pnam(idx), 'Unif', false); 
                pmin = cellfun(@(x) num2str(x, '%#-8.3g'), {obj.PMIN}, 'Unif',false);
                pmin(cellfun(@isempty, pmin)) = {blanks(8)};
                pmax = cellfun(@(x) num2str(x, '%#-8.3g'), {obj.PMAX}, 'Unif',false);
                pmax(cellfun(@isempty, pmax)) = {blanks(8)};
                coef = {obj.COEFi};
                idx  = cellfun(@isempty, coef);
                coef(idx) = {'1.0'};
                coef(~idx) = cellfun(@(x) num2str(x, '%#-8.3g'), coef(~idx), 'Unif', false);
                dvid = get([obj.DesignVariables], 'ID')';
                
                %Set up the format for printing
                data = [ ...
                    nam   ; {obj.ID} ; {obj.TYPE} ; pid ; pnam ; pmin ; pmax ; {obj.C0} ; ...
                    blnks ; dvid     ; coef ];          
                
                %Write in 8-character column width as standard
                format = ['%-8s%-8i%-8s%-8i%-8s%-8s%-8s%#-8.4g\r\n', ...
                    '%-8s%-8i%-8s\r\n'];
                
                %Write the data to the file
                fprintf(fid, format, data{:});
                
            end

            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

