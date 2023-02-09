classdef DesignResponseEqn < awi.fe.FEBaseClass
    %DesignResponseEqn Describes a design response as a combination of
    %various structural design responses and other output quantities.
    %
    % The definition of the 'awi.fe.opt.DesignResposeEqn' class matches
    % that of the DRESP2 bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %String for identifying this design response
        LABEL
        %ID number of the design equation
        EQID
        %Name of the built-in FORTAN function that MSC.Nastran will use
        FUNC
        %Region identifer for constraint screening
        Region
        %Method type for FUNC = BETA
        METHOD
        %Constants used when FUNC = BETA or FUNC = MATCH
        Constants = [1, 0.005, 10];        
        %Name of constants relating to a DTABLE entry
        DesignConstants = {};
    end
    
    %Handle to 'awi.fe' objects
    properties (Hidden = true)
        %Handle to a single 'awi.fe.opt.DesignEquation' object
        DesignEquation
        %Handle to the underlying 'awi.fe.opt.DesignVariable' objects
        DesignVariables
        %Handle to the underlying 'awi.fe.opt.PropRelationEqn' objects
        PropRelationEqns
        %Handle to the underlying 'awi.fe.opt.DesignResponse' objects
        DesignResponses
        %Handle to the underlying 'awi.fe.opt.DesignResponseEqn' object
        DesignResponseEqns
    end
    
    methods % set / get
        function set.LABEL(obj, val)              %set.LABEL
            %set.LABEL Set method for the property 'LABEL'.
            
            validateLabel(obj, val, 'LABEL')
            obj.LABEL = val;
        end
        function set.EQID(obj, val)               %set.EQID
            %set.EQID Set method for the property 'EQID'.
            
            validateID(obj, val, 'EQID');
            obj.EQID = val;
            
        end
        function set.FUNC(obj, val)               %set.FUNC
            %set.FUNC Set method for the property 'FUNC'.
            
            validateLabel(obj, val, 'FUNC');
            obj.FUNC = val;
        end
        function set.Region(obj, val)             %set.Region
            %set.Region Set method for the property 'Region'
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
                class(obj), 'Region');
            obj.Region = val;
        end
        function set.METHOD(obj, val)             %set.METHOD
            %set.METHOD Set method for the property 'METHOD'.
            
            validateLabel(obj, val, 'METHOD')
            obj.METHOD = val;
        end
        function set.Constants(obj, val)          %set.Constants
           %set.Constants Set method for the property 'Constants'.
           
           validateattributes(val, {'numeric'}, {'row', 'numel', 3}, ...
               class(obj), 'Constants');
           obj.Constants = val;
        end
        function set.DesignEquation(obj, val)     %set.DesignEquation
            %set.DesignEquation Set method for the property
            %'DesignEquation'.
            
            validateattributes(val, {'awi.fe.opt.DesignEquation'}, ...
                {'scalar', 'nonempty'}, class(obj), 'DesignEquation');
            obj.DesignEquation = val;
        end
        function set.DesignVariables(obj, val)    %set.DesignVariables
            %set.DesignVariables Set method for the property
            %'DesignVariables'.
            
            validateattributes(val, {'awi.fe.opt.DesignVariable'}, {'row', ...
                'nonempty'}, class(obj), 'DesignVariables');
            obj.DesignVariables = val;
        end
        function set.PropRelationEqns(obj, val)   %set.PropRelationEqn
            %set.PropRelationEqn Set method for the property
            %'PropRelationEqns'.
            
            validateattributes(val, {'awi.fe.opt.PropRelationEqn'}, {'row', ...
                'nonempty'}, class(obj), 'PropRelationEqns');
            obj.PropRelationEqns = val;
        end
        function set.DesignConstants(obj, val)    %set.DesignConstants
            %set.DesignConstants Set method for the property
            %'DesignConstants'.
            
            assert(iscellstr(val), ['Expected ''DesignConstants'' to ', ...
                'be a cell array of strings']);
            idx = cellfun(@(x) numel(x) < 9, val);
            assert(all(idx), ['Expected each property name in ', ...
                '''DesignConstants'' to be less than 8 characters.']);
            obj.DesignConstants = val;
        end
        function set.DesignResponses(obj, val)    %set.DesignResponses
            %set.DesignResponses Set method for the property
            %'DesignResponses'.
            
            validateattributes(val, {'awi.fe.opt.DesignResponse'}, ...
                {'row', 'nonempty'}, class(obj), 'DesignResponses');
            obj.DesignResponses = val;
        end
        function set.DesignResponseEqns(obj, val) %set.DesignResponseEqns
            %set.DesignResponseEqns Set method for the property
            %'DesignResponseEqns'.
            
            validateattributes(val, {'awi.fe.opt.DesignResponseEqn'}, ...
                {'row', 'nonempty'}, class(obj), 'DesignResponseEqns');
            obj.DesignResponseEqns = val;
        end
        function val = get.EQID(obj)              %get.EQID
           %get.EQID Get method for the property 'EQID'.
           %
           % If the object has been assigned a handle to its
           % 'awi.fe.opt.DesignEquation' object then always use their ID
           % number, else use EQID.
           
           if isempty(obj.DesignEquation)
               val = obj.EQID;
           else
               val = obj.DesignEquation.ID;
           end
           
        end
    end
    
    methods % construction
        function obj = DesignResponseEqn
            
            %Make a note of the property names
            addFEProp(obj, 'ID', 'LABEL', 'EQID', 'FUNC', 'Region', 'METHOD', 'Constants');
            
            %Update the card name
            obj.CardName = 'DRESP2';
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the
            %'awi.fe.opt.DesignResponseEqn' object into a text file using
            %the format of the MSC.Nastran 'DRESP2' bulk data entry.
            %
            % The following assumptions apply:
            %   * The 'REGION' property is ignored and set to blank.
            
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
                comment = ['DRESP2: Defines equation responses that '  , ...
                    'are used in the design, either as constraints or ', ...
                    'as an objective.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            nObj  = numel(obj);
            blnks = blanks(8);
            
            %Have to use a loop as there are too many permutations of the
            %DRESP2 card to possibly vectorise            
            for ii = 1 : nObj
                %Data for the header
                eqn = num2str(obj(ii).EQID);
                if isempty(eqn)
                    eqn = obj(ii).FUNC;
                    C   = cellstr(num2str(obj(ii).Constants'));
                else
                    C   = repmat({blanks(8)}, [1, 3]);
                end
                %Header
                fprintf(fid, '%-8s%-8i%-8s%-8s%-8s%-8s%-8s%-8s\r\n', ...
                    'DRESP2', obj(ii).ID, obj(ii).LABEL, eqn, blnks, C{1}, C{2}, C{3});
                %DESVAR (Design Variables)
                if ~isempty(obj(ii).DesignVariables)
                    data = get(obj(ii).DesignVariables, {'ID'});
                    i_writeList(fid, 'DESVAR', data, '%-8i');
                end
                %DTABLE (Design Constants)
                if ~isempty(obj(ii).DesignConstants)
                    data = obj(ii).DesignConstants;
                    i_writeList(fid, 'DTABLE', data, '%-8s');
                end
                %DRESP1 (Design Responses)
                if ~isempty(obj(ii).DesignResponses)
                    data = get(obj(ii).DesignResponses, {'ID'});
                    i_writeList(fid, 'DRESP1', data, '%-8i');
                end
                %DVPREL2 (Design Property Relation Equations)
                if ~isempty(obj(ii).PropRelationEqns)
                    data = get(obj(ii).PropRelationEqns, {'ID'});
                    i_writeList(fid, 'DVPREL2', data, '%-8i');
                end
                %DRESP2 (Design Response Equations
                if ~isempty(obj(ii).DesignResponseEqns)
                    data = get(obj(ii).DesignResponseEqns, {'ID'});
                    i_writeList(fid, 'DRESP2', data, '%-8i');
                end
            end                            

            if bClose %Close the file?
                fclose(fid);
            end            
            
            function i_writeList(fid, token, data, fmt)
                %i_writeList Writes the data in 'data' into small-field
                %format and splits the data over multiple lines to ensure
                %that is meets the 80-character width.
                
                %Assume we have loads of data to write
                nCol = 7;
                nData = numel(data);
                nRows = floor(nData / nCol);
                nRem  = mod(nData, nCol);
                
                %Force row vector
                if iscolumn(data)
                    data = data';
                end
                
                %Split data up
                if nRows == 0
                    fmt  = ['%-8s%-8s', repmat(fmt, [1, nData]), '\r\n'];
                    data = [{blanks(8)}, token, data];
                    fprintf(fid, fmt, data{:});
                else
                    error('Update code')
                    %TODO - Update the code with the method from
                    %'awi.fe.BeamCrossSection'
                    %Split into rows
                    rowData = data(1 : nRows * nCol, :);
                    remData = data(end - (rem - 1) : end, :);
                end
                
            end
            
        end
    end
    
end

