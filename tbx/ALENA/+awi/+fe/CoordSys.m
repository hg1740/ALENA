classdef CoordSys < awi.fe.FEBaseClass
    %CoordSys Describes a cartesian coordinate system.
    %
    % The definition of the 'CoordSys' object matches that of the 'CORD2R'
    % bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Identification number of this coordinate system
        CID        
        %Coordinates of the origin
        A
        %Coordinates of a point on the local z-axis
        B
        %Coordinates of a point in the XZ-plane
        C
    end
    
    properties (Constant)
        %Identification number of the coordinate system - Always 0 (basic)
        RID = 0;
    end
    
    methods % set / get
        function set.CID(obj, val)  %set.CID 
            %set.RID Set method for the property 'CID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.A(obj, val)    %set.A   
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'A');
            obj.A = val;
        end
        function set.B(obj, val)    %set.B   
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'B');
            obj.B = val;
        end
        function set.C(obj, val)    %set.C   
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'C');
            obj.C = val;
        end
        function val = get.CID(obj) %get.CID 
            val = obj.ID;
        end
    end
    
    methods % construction
        function obj = CoordSys
            
            %Make a note of the property names
            addFEProp(obj, 'CID', 'RID', 'A', 'B', 'C');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            %drawElement Draw method for the 'awi.fe.CoordSys' class.
            
            %Collect data (accounts for object arrays)
            A_ = [obj.A];
            
            %Calculate axis vectors
            [~, eX, eY, eZ] = getRotationMatrix(obj);
            
            %Create coordinates
            OX = A_ + eX;
            OY = A_ + eY;
            OZ = A_ + eZ;
            
            %Plot eX
            eXx = obj.padCoordsWithNaN([A_(1, :) ; OX(1, :)]);
            eXy = obj.padCoordsWithNaN([A_(2, :) ; OX(2, :)]);
            eXz = obj.padCoordsWithNaN([A_(3, :) ; OX(3, :)]);            
            hg(1) = line(ha, eXx, eXy, eXz, ...
                'Color'    , 'r', ...
                'LineWidth', 2  , ...
                'Tag', 'Coordinate Systems (X)');
            
            %Plot eY
            eYx = obj.padCoordsWithNaN([A_(1, :) ; OY(1, :)]);
            eYy = obj.padCoordsWithNaN([A_(2, :) ; OY(2, :)]);
            eYz = obj.padCoordsWithNaN([A_(3, :) ; OY(3, :)]);             
            hg(2) = line(ha, eYx, eYy, eYz, ...
                'Color'    , 'g', ...
                'LineWidth', 2  , ...
                'Tag', 'Coordinate Systems (Y)');
            
            %Plot eZ
            eZx = obj.padCoordsWithNaN([A_(1, :) ; OZ(1, :)]);
            eZy = obj.padCoordsWithNaN([A_(2, :) ; OZ(2, :)]);
            eZz = obj.padCoordsWithNaN([A_(3, :) ; OZ(3, :)]);             
            hg(3) = line(ha, eZx, eZy, eZz, ...
                'Color'    , 'b', ...
                'LineWidth', 2  , ...
                'Tag', 'Coordinate Systems (Z)');
            
            %Return a column vector of handles
            hg = hg(:);

        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.CoordSys' object
            %into a text file using the format of the MSC.Nastran 'CORD2R'
            %bulk data entry.
            %
            % The following assumptions are made:
            %   * The coordinate system defining the points A, B & C is
            %   assumed to the basic coordinate system. e.g. RID = 0.
            
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
                comment = ['CORD2R: Defines a rectangular coordinate ', ...
                    'system using the coordinates of three points.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'large');
                      
            nObj = numel(obj);
            
            %Split up the coordinates            
            A_ = [obj.A];
            B_ = [obj.B];
            C_ = [obj.C];
            A1 = num2cell(A_(1, :));
            A2 = num2cell(A_(2, :));
            A3 = num2cell(A_(3, :));
            B1 = num2cell(B_(1, :));
            B2 = num2cell(B_(2, :));
            B3 = num2cell(B_(3, :));
            C1 = num2cell(C_(1, :));
            C2 = num2cell(C_(2, :));
            C3 = num2cell(C_(3, :));
            
            %Card name
            nam   = repmat({'CORD2R*'}  , [1, nObj]);
            blnks = repmat({['*', blanks(7)]}, [1, nObj]);
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.RID} ; A1 ; A2 ; ...
                blnks ; A3       ; B1        ; B2 ; B3 ; ...
                blnks ; C1       ; C2        ; C3];
            
            %Write in 16-character column width as standard
            format = [ ...
                '%-8s%-16i%-16i%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g%#-16.8g\r\n', ...
                '%-8s%#-16.8g%#-16.8g%#-16.8g\r\n'];

%                     %Set up the format for printing
%             data = [ ...
%                 nam   ; {obj.ID} ; {obj.RID} ; A1 ; A2 ; A3 ; B1 ; B2 ; B3 ; ...
%                 blnks ; C1       ; C2        ; C3];
%             
%             %Write in 8-character column width
%             format = [ ...
%                 '%-8s%-8i%-8i%-8.3g%-8.3g%-8.3g%-8.3g%-8.3g%-8.3g\r\n', ...
%                 '%-8s%-8.3g%-8.3g%-8.3g\r\n'];

            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
    methods % helper functions
        function [R, eX, eY, eZ] = getRotationMatrix(obj)
            %getRotationMatrix Returns the rotation matrix [eX, eY, eZ]
            %described by this coordinate system.
            
            %Collect data (accounts for object arrays)
            A_ = [obj.A];
            B_ = [obj.B];
            C_ = [obj.C];
            
            %Calculate axis vectors
            eZ = B_ - A_;
            AC = C_ - A_;            
            eY = cross(eZ, AC);
            eX = cross(eY, eZ);
            
            %Force unit vectors 
            eX = eX ./ vecnorm(eX);
            eY = eY ./ vecnorm(eY);
            eZ = eZ ./ vecnorm(eZ);
                      
            %Collect
            R = cat(3, eX, eY, eZ);
            R = permute(R, [1, 3, 2]);
            
        end
    end
    
end

