classdef BeamProp < awi.fe.FEBaseClass
    %BeamProp Describes the beam properties of a flexible beam for use in 
    %a finite element model.
    %
    % The definition of the 'BeamProp' object matches that of the PBEAM 
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Property identification number
        PID
        %Material identification number
        MID
        %Cross-sectional area 
        A
        %Second moment of area about plane-1
        I11
        %Second moment of area about plane-2
        I22
        %Cross-second moment of area
        I12
        %Polar moment of area
        J
        %Non-structural mass
        NSM
        %Nonstructural inertia
        NSI
        %Centre of mass offset in the y-direction of the local beam
        %coordinate system
        CMy
        %Centre of mass offset in the z-direction of the local beam
        %coordinate system
        CMz
        %Neutral axis offset in the y-direction of the local beam
        %coordinate system
        NAy
        %Neutral axis offset in the z-direction of the local beam
        %coordinate system
        NAz
        
        %Stress recovery point, top right corner
        Cy
        Cz
        %Stress recovery point, bottom right corner
        Dy
        Dz
        %Stress recovery point, bottom left corner
        Ey
        Ez
        %Stress recovery point, top left corner
        Fy
        Fz 
    end
    
    %Store a handle to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.Material' object
        Material     
    end
    
    methods % set / get
        function set.PID(obj, val)      %set.PID
            %set.PID Set method for the property 'PID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.MID(obj, val)      %set.MID
            validateID(obj, val, 'MID')
            obj.MID = val;
        end
        function set.A(obj, val)        %set.A
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'A');
            obj.A = val;
        end
        function set.I11(obj, val)      %set.I11
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'I11');
            obj.I11 = val;
        end
        function set.I22(obj, val)      %set.I22
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'I22');
            obj.I22 = val;
        end
        function set.I12(obj, val)      %set.I12
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'I12');
            obj.I12 = val;
        end
        function set.J(obj, val)        %set.J
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'J');
            obj.J = val;
        end
        function set.NSM(obj, val)      %set.NSM
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'NSM');
            obj.NSM = val;
        end
        function set.NSI(obj, val)      %set.NSI
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'NSI');
            obj.NSI = val;
        end
        function set.CMy(obj, val)      %set.CMy
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'CMy');
            obj.CMy = val;
        end
        function set.CMz(obj, val)      %set.CMz
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'CMz');
            obj.CMz = val;
        end
        function set.NAy(obj, val)      %set.NAy
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'NAy');
            obj.NAy = val;
        end
        function set.NAz(obj, val)      %set.NAz
            val(isnan(val)) = 0;    %Replace any NaN terms with zeros
            validateBeamProp(obj, val, 'NAz');
            obj.NAz = val;
        end
        function set.Cy(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Cy');
            obj.Cy = val;
        end
        function set.Cz(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Cz');
            obj.Cz = val;
        end
        function set.Dy(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Dy');
            obj.Dy = val;
        end
        function set.Dz(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Dz');
            obj.Dz = val;
        end
        function set.Ey(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Ey');
            obj.Ey = val;
        end
        function set.Ez(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Ez');
            obj.Ez = val;
        end
        function set.Fy(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Fy');
            obj.Fy = val;
        end
        function set.Fz(obj, val)
            val(isnan(val)) = 0;
            validateBeamProp(obj, val, 'Fz');
            obj.Fz = val;
        end
        function set.Material(obj, val) %set.Material
            %set.Material Set method for the property 'Material'.
            %
            % 'Material' must be a scalar 'awi.fe.Material' object.
            
            validateattributes(val, {'awi.fe.Material'}, {'scalar', ...
                'nonempty'}, class(obj), 'Material');
            obj.Material = val;
        end
        function val = get.PID(obj)     %get.PID
            %get.PID Get method for the property 'PID'.
            
            val = obj.ID;
        end
        function val = get.MID(obj)     %get.MID
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
        function obj = BeamProp
            
            %Make a note of the property names
            addFEProp(obj, 'PID', 'MID', 'A', 'I11', 'I22', 'I12', 'J', ...
                'NSM', 'NSI', 'CMy', 'CMz', 'NAy', 'NAz', ...
                'Cy', 'Cz', 'Dy', 'Dz', 'Ey', 'Ez', 'Fy', 'Fz');
            
        end
    end
    
    methods % writing to a .txt file
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.BeamProp' object
            %into a text file using the format of the MSC.Nastran 'PBEAM'
            %bulk data entry.
            %
            % The following assumptions are made:
            %   * The beam is assumed to conform to Euler-Bernoulli beam 
            %   theory, i.e. "plane sections remain plane" - K1 = K2 = 0 &
            %   there is no warping of the cross-section - CW = 0.
            %   * The beam is assumed to ignore any shear relief due to
            %   taper of the beam. i.e. S1 = S2 = 0.
            %   * The direct stress due to bending/axial behaviour of the
            %   beam is not recovered.
            
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
                comment = ['PBEAM : Defines the properties of a ', ...
                    'tapered beam element'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
                 
            nObj = numel(obj);
            
            %Grab the beam properties
            A_   = num2cell([obj.A]);
            I11_ = num2cell([obj.I11]);
            I22_ = num2cell([obj.I22]);
            I12_ = num2cell([obj.I12]);
            J_   = num2cell([obj.J]);
            NSM_ = num2cell([obj.NSM]);
            NSI_ = num2cell([obj.NSI]);
            M1   = num2cell([obj.CMy]);
            M2   = num2cell([obj.CMz]);
            N1   = num2cell([obj.NAy]);
            N2   = num2cell([obj.NAz]);
            
            %Stress recovery points relative to the shear center
            Cy_   = num2cell([obj.Cy]); % Top right corner
            Cz_   = num2cell([obj.Cz]);
            Dy_   = num2cell([obj.Dy]); % Bottom right corner
            Dz_   = num2cell([obj.Dz]);
            Ey_   = num2cell([obj.Ey]); % Bottom left corner
            Ez_   = num2cell([obj.Ez]);
            Fy_   = num2cell([obj.Fy]); % Top left corner
            Fz_   = num2cell([obj.Fz]);          
                        
            %Card name, blanks, zero/ones arrays & stress recovery flag  
            nam   = repmat({'PBEAM*'}  , [1, nObj]);
            blnks = repmat({['*', blanks(7)]}, [1, nObj]);
            zrs   = num2cell(zeros(1, nObj));
            onz   = num2cell(ones(1, nObj));
            SO    = repmat({'YES'}, [1, nObj]);
            
            %Set up the format for printing
%             data = [ ...
%                 nam   ; {obj.ID}   ; {obj.MID}  ; A_(1, :) ; I11_(1, :) ; ... % PBEAM,  PID, MID, A(A), I1(A)
%                 blnks ; I22_(1, :) ; I12_(1, :) ; J_(1, :) ; NSM_(1, :) ; ... %         I2(A), I12(A), J(A), NSM(A)
%                 blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         C1(A), C2(A), D1(A), D2(A)      [1 = y, 2 = z]
%                 blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         E1(A), E2(A), F1(A), F2(A)      [1 = y, 2 = z]
%                 blnks ; SO         ; onz        ; A_(2, :) ; I11_(2, :) ; ... %         SO, X/XB, A, I1
%                 blnks ; I22_(2, :) ; I12_(2, :) ; J_(2, :) ; NSM_(2, :) ; ... %         I2, I12, J, NSM
%                 blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         C1, C2, D1, D2                  [1 = y, 2 = z]
%                 blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         E1, E2, F1, F2                  [1 = y, 2 = z]
%                 blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         K1, K2, S1, S2                  [1 = y, 2 = z]
%                 blnks ; NSI_(1, :) ; NSI_(2, :) ; zrs      ; zrs        ; ... %         NSI(A), NSI(B), CW(A), CW(B)
%                 blnks ; M1(1, :)   ; M2(1, :)   ; M1(2, :) ; M2(2, :)   ; ... %         M1(A), M2(A), M1(B), M2(B)      [1 = y, 2 = z]
%                 blnks ; N1(1, :)   ; N2(1, :)   ; N1(2, :) ; N2(2, :)  ];     %         N1(A), N2(A), N1(B), N2(B)      [1 = y, 2 = z]
            
            data = [ ...
                nam   ; {obj.ID}   ; {obj.MID}  ; A_(1, :) ; I11_(1, :) ; ... % PBEAM,  PID, MID, A(A), I1(A)
                blnks ; I22_(1, :) ; I12_(1, :) ; J_(1, :) ; NSM_(1, :) ; ... %         I2(A), I12(A), J(A), NSM(A)
                blnks ; Cy_(1, :)  ; Cz_(1, :)  ; Dy_(1, :); Dz_(1, :)  ; ... %         C1(A), C2(A), D1(A), D2(A)      [1 = y, 2 = z]
                blnks ; Ey_(1, :)  ; Ez_(1, :)  ; Fy_(1, :); Fz_(1, :)  ; ... %         E1(A), E2(A), F1(A), F2(A)      [1 = y, 2 = z]
                blnks ; SO         ; onz        ; A_(2, :) ; I11_(2, :) ; ... %         SO, X/XB, A, I1
                blnks ; I22_(2, :) ; I12_(2, :) ; J_(2, :) ; NSM_(2, :) ; ... %         I2, I12, J, NSM
                blnks ; Cy_(2, :)  ; Cz_(2, :)  ; Dy_(2, :); Dz_(2, :)  ; ... %         C1, C2, D1, D2                  [1 = y, 2 = z]
                blnks ; Ey_(2, :)  ; Ez_(2, :)  ; Fy_(2, :); Fz_(2, :)  ; ... %         E1, E2, F1, F2                  [1 = y, 2 = z]
                blnks ; zrs        ; zrs        ; zrs      ; zrs        ; ... %         K1, K2, S1, S2                  [1 = y, 2 = z]
                blnks ; NSI_(1, :) ; NSI_(2, :) ; zrs      ; zrs        ; ... %         NSI(A), NSI(B), CW(A), CW(B)
                blnks ; M1(1, :)   ; M2(1, :)   ; M1(2, :) ; M2(2, :)   ; ... %         M1(A), M2(A), M1(B), M2(B)      [1 = y, 2 = z]
                blnks ; N1(1, :)   ; N2(1, :)   ; N1(2, :) ; N2(2, :)  ];     %         N1(A), N2(A), N1(B), N2(B)      [1 = y, 2 = z]
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%-16i%#-16.8g%#-16.8g\r\n'      , ... 
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%-16s%#-16.8g%#-16.8g%#-16.8g\r\n'   , ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

