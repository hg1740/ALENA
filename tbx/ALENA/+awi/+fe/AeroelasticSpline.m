classdef AeroelasticSpline < awi.fe.FEBaseClass
    %AeroelasticSpline Defines a generic spline between some structural
    %and aerodynamic degrees-of-freedom.
    
    %Primary Properties
    properties
        %Identification number for this spline
        EID
        %List of ID numbers for 'awi.fe.Node' objects belonging to this
        %spline
        StrucSetID
        %List of aerodynamic panel ID numbers for 'awi.fe.Node' objects
        %belonging to this spline
        AeroPanelSetID
        %ID number of the coordinate system defining the y-axis of the
        %spline
        CID
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.StructuralSet' object
        StrucSet
        %Handle to the 'awi.fe.AeroPanel' object
        AeroPanelSet
        %Handle to the 'awi.fe.CoordSys' object
        CoordSys
    end
    
    methods % set / get
        function set.EID(obj, val)             %set.EID
            %set.EID Set method for the property 'EID'.
            
            %Pass it straight to the hidden 'ID' property
            obj.ID = val;
        end
        function set.StrucSetID(obj, val)      %set.StructSetID
            %set.StrucSetID Set method for the property 'StrucSetID'
            validateID(obj, val, 'StrucSetID');
            obj.StructSetID = val;
        end
        function set.AeroPanelSetID(obj, val)  %set.AeroPanelSetID
            %set.StrucSetID Set method for the property 'StrucSetID'
%             validateattributes(val, {'numeric'}, {'column', 'integer', ...
%                 'nonnegative', 'finite', 'real'}, class(obj), 'AeroPanelSetID');
            validateID(obj, val, 'AeroPanelSetID');
            obj.AeroPanelSetID = val;
        end
        function set.CID(obj, val)             %set.CID
            %set.CID Set method for the property 'CID'.
            validateID(obj, val, 'CID');
            obj.CID = val;
        end
        function set.StrucSet(obj, val)        %set.StrucSet
            %set.StrucSet Set method for the property 'StrucSet'
            validateattributes(val, {'awi.fe.StructuralSet'}, {'scalar', ...
                'nonempty'}, class(obj), 'StrucSet');
            obj.StrucSet = val;
        end
        function set.AeroPanelSet(obj, val)    %set.AeroPanelSet
            %set.AeroPanelSet Set method for the property 'AeroPanelSet'
            validateattributes(val, {'awi.fe.AeroPanelSet'}, {'scalar', ...
                'nonempty'}, class(obj), 'AeroPanelSet');
            obj.AeroPanelSet = val;
        end
        function val = get.EID(obj)            %get.EID
            %get.EID Get method for the property 'EID'.
            
            val = obj.ID;
            
        end
        function val = get.StrucSetID(obj)     %get.StrucSetID
            %get.StrucSetID Get method for the property 'StrucSetID'.
            %
            % If the object has been assigned a handle to its
            %'awi.fe.StructuralSet' object then always use its ID number,
            % else use 'StrucSet'.
            if isempty(obj.StrucSet)
                val = obj.StrucSetID;
            else
                val = obj.StrucSet.ID;
            end
        end
        function val = get.AeroPanelSetID(obj) %get.AeroPanelID
            %get.AeroPanelSetID Get method for the property 'AeroPanelID'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.AeroPanel' object then always use its ID number, else
            % use 'AeroPanelID'.
            if isempty(obj.AeroPanelSet)
                val = obj.AeroPanelSetID;
            else
                val = obj.AeroPanelSet.ID;
            end
        end
        function val = get.CID(obj)            %get.CID
            %get.CID Get method for the property 'CID'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.CoordSys' object then always use its ID number, else
            % use 'CID'.
            if isempty(obj.CoordSys)
                val = obj.CID;
            else
                val = obj.CoordSys.ID;
            end
        end
    end
    
    methods % construction
        function obj = AeroelasticSpline
            
            addFEProp(obj, 'StrucSetID', 'AeroPanelSetID', 'CID');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.AeroelasticSpline'
            %object into a text file using the format of the MSC.Nastran
            %'SPLINE7' bulk data entry.
            %
            % The following assumptions are made:
            %   * The linear attachment flexibility is assumed to be zero,
            %   i.e. DZ = 0.0
            %   * The ratio of beam bending stiffness to the beam
            %   torsionsal stiffness is assumed to be 1.0, e.g. DTOR = 1.0
            %   * The spline is assumed to map both forces and
            %   displacements, i.e. USAGE = BOTH
            %   * The spline method is assumed to be 'FBS6'.
            %   * The rotational attachment flexibility is assumed to be
            %   0.0, i.e. DZR = 0.0
            %   * The ratio of the beam bending stiffness to the beam
            %   extensional stiffness is assumed to be 1.0, i.e. IA2 = 1.0
            %   * The ratio of the minimum beam length to the total beam
            %   length is assumed to be 0.01, i.e. EPSM = 0.01
            
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
                comment = ['SPLINE7 : Defines a 6DOF finite beam ', ...
                    'spline for interpolating motion and/or '     , ...
                    'forces between two meshes.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            %How many objects?
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'SPLINE7'}, [1, nObj]);
            blnks = repmat({blanks(8)}, [1, nObj]);
            
            %Set up data
            DZ    = num2cell(zeros([1, nObj]));
            DTOR  = num2cell(ones([1, nObj]));
            usge  = repmat({'BOTH'}, [1, nObj]);
            meth  = repmat({'FBS6'}, [1, nObj]);
            DZR   = num2cell(zeros([1, nObj]));
            IA2   = num2cell(ones([1, nObj]));
            EPSBM = repmat({0.01}  , [1, nObj]);
            
            %Grab the set object
            Struc  = [obj.StrucSet];
            Panels = [obj.AeroPanelSet];            
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; blnks ; {Panels.ID} ; blnks ; {Struc.ID} ; DZ  ; DTOR ; {obj.CID} ; ...
                blnks ; blnks    ; blnks ; blnks       ; usge  ; meth       ; DZR ; IA2  ; EPSBM ];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-8i%-s%#-8i%-8s%-8i%#-8.3g%#-8.3g%-8i\r\n', ...
                '%-8s%-8s%-8s%-8s%-8s%-8s%#-8.3g%#-8.3g%#-8.3g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});

%             %Card name
%             nam   = repmat({'SPLINE4'}, [1, nObj]);
%             blnks = repmat({blanks(8)}, [1, nObj]);
%             
%             %Set up data
%             DZ    = num2cell(zeros([1, nObj]));
%             usge  = repmat({'BOTH'}, [1, nObj]);
%             meth  = repmat({'IPS'}, [1, nObj]);            
%             
%             %Grab the set object
%             Struc  = [obj.StrucSet];
%             Panels = [obj.AeroPanelSet];
%             
%             
%             %Set up the format for printing
%             data = [nam ; {obj.ID} ; blnks ; {Panels.ID} ; blnks ; {Struc.ID} ; DZ  ; meth ; usge];
%             
%             %Write in 16-character column width as standard
%             format = '%-8s%-8i%-s%#-8i%-8s%-8i%#-8.3g%#-8s%-8s\r\n';
%             
%             %Write the data to the file
%             fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

