classdef AeroControlSurf < awi.fe.FEBaseClass
    %AeroControlSurf Defines an aerodynamic control surface.
    
    %Primary properties
    properties
        %Name of the control surface
        LABEL 
        %ID number of the coordinate system that defines the axis of
        %rotation of the control surface
        CID
        %ID number of the AELIST entry that defines the aerodynamic panel
        %numbers
        ALID
        %Control surface effectiveness
        EFF = 1.0;
        %Linear (LDW) or nonlinear downwash (NOLDW) key
        LDW = 'LDW';
    end    
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)  
        %Handle to the 'awi.fe.CoordSys' object that defines the axis of
        %rotation of the control surface
        CoordSys
        %Handle to the 'awi.fe.AeroPanelSet' object that defines the panel
        %ID numbers
        AeroPanelSet        
    end
    
    methods % set / get
        function set.LABEL(obj, val)        %set.LABEL       
            validateattributes(val, {'char'}, {'row', 'nonempty'}, ...
                class(obj), 'LABEL');
            obj.LABEL = val;
        end
        function set.CID(obj, val)          %set.CID         
            validateID(obj, val, 'CID');
            obj.CID = val;
        end
        function set.ALID(obj, val)         %set.ALID        
            validateID(obj, val, 'ALID');
            obj.ALID = val;
        end
        function set.EFF(obj, val)          %set.EFF         
           validateattributes(val, {'numeric'}, {'<=', 1, '>', 0}, ...
               class(obj), 'EFF');
           obj.EFF = val;
        end
        function set.LDW(obj, val)          %set.LDW         
            validatestring(val, {'LDW', 'NOLDW'}, class(obj), 'LDW');
            obj.LDW = val;
        end
        function set.CoordSys(obj, val)     %set.CoordSys    
           %set.CoordSys Set method for the property 'CoordSys'. 
            %
            % 'CoordSys' must be a valid instance of the 'awi.fe.CoordSys'
            % class.
            validateattributes(val, {'awi.fe.CoordSys'}, {'scalar', ...
                'nonempty'}, class(obj), 'CoordSys');
            obj.CoordSys = val;
        end
        function set.AeroPanelSet(obj, val) %set.AeroPanelSet
            %set.CoordSys Set method for the property 'AeroPanelSet'.
            %
            % 'AeroPaneSet' must be a valid instance of the
            % 'awi.fe.AeroPanelSet' class.
            validateattributes(val, {'awi.fe.AeroPanelSet'}, {'scalar', ...
                'nonempty'}, class(obj), 'AeroPanelSet');
            obj.AeroPanelSet = val;
        end
        function val = get.CID(obj)         %get.CID         
            %get.CID Get method for the property 'CID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.CoordSys' object then always use their ID number, 
            % else use CID.
            if isempty(obj.CoordSys)
                val = obj.CID;
            else
                val = obj.CoordSys.ID;
            end
        end
        function val = get.ALID(obj)        %get.ALID        
            %get.ALID Get method for the property 'ALID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.AeroPanelSet' object then always use their ID number, 
            % else use ALID.
            if isempty(obj.AeroPanelSet)
                val = obj.ALID;
            else
                val = obj.AeroPanelSet.ID;
            end
        end
    end
    
     methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.AeroControlSurf'
            %object into a text file using the format of the MSC.Nastran 
            %'AESURF' bulk data entry.
            %
            % As well as writing the AESURF card, the following cards are
            % also written:
            %   * CAERO1 - The CAERO1 
            %            % The following assumptions are made:
            %   * The continuation entry is neglected. Therefore the
            %     following properties have their default values:
            %       - CREFC  : 1.0
            %       - CREFS  : 1.0
            %       - PLLIM  : -pi/2
            %       - PULIM  : +pi/2
            %       - HMLLIM : inf
            %       - HMULIM : inf
            %       - TGLLIM : na (No table referenced)
            %       - TGULIM : na (No table referenced)
                        
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
                comment = ['AESURF : Specifies an aerodynamic control ', ...
                    'surface as a member of the set of aerodynamic '   , ...
                    'extra points.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            %How many objects?
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'AESURF'} , [1, nObj]);
            blnks = repmat({blanks(8)}, [1, nObj]);
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.LABEL} ; {obj.CID} ; {obj.ALID} ; blnks ; blnks ; {obj.EFF} ; {obj.LDW}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-8i%-8s%-8i%-8i%-8s%-8s%-8.3f%-8s\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
     end
    
end

