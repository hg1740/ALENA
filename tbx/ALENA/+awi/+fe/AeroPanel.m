classdef AeroPanel < awi.fe.FEBaseClass
    %AeroPanel Defines a set of aerodynamic panels for use in a finite
    %element model.
    %
    % The definition of the 'AeroPanel' object matches that of the CAERO1
    % bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Element ID number of the first aeroydnamic panel in each set
        EID
        %Property ID number of the aerodynamic property group
        PID
        %Number of chordwise panels in each set
        NCHORD
        %Number of spanwise panels in each set
        NSPAN
        %(X,Y,Z) coordinates of the inboard LE corner in each set
        X1
        %(X,Y,Z) coordinates of the outboard LE corner in each set
        X4
        %Chord values in the global X-axis in each set
        CHORD
    end
    
    %Downwash distribution over the panel
    properties
        %Normalised distibution of Angle-of-Attack (AoA) along the panel
        %set (X1 -> X4)
        AoA_eta
        %AoA values along the panel set at the points given by 'AoA_eta'.
        AoA
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.AeroProp' object that describes the
        %aerodynamic properties.
        AeroProp
    end
    
    properties (Dependent) % NumPanels, PanelSetCornerCoords
        %Total number of panels in each section
        NumPanels
        %The corner coordinates of each panel set
        PanelSetCornerCoords
    end
    
    properties (SetAccess = private)
        PanelCornerCoords
    end
    
    methods % set / get
        function set.EID(obj, val)        %set.EID       
            %set.EID Set method for the property 'EID'.
            % 
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function set.PID(obj, val)        %set.PID       
            validateID(obj, val, 'PID')
            obj.PID = val;
        end
        function set.NCHORD(obj, val)     %set.NCHORD    
            validateattributes(val, {'numeric'}, {'integer', 'scalar', ...
                'positive', 'nonnan', 'finite', 'real'}, class(obj), 'NCHORD');
            obj.NCHORD = val;
        end
        function set.NSPAN(obj, val)      %set.NSPAN     
            validateattributes(val, {'numeric'}, {'integer', 'scalar', ...
                'positive', 'nonnan', 'finite', 'real'}, class(obj), 'NSPAN');
            obj.NSPAN = val;
        end
        function set.X1(obj, val)         %set.X1        
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'X1');
            obj.X1 = val;
        end
        function set.X4(obj, val)         %set.X4        
            validateattributes(val, {'numeric'}, {'column', 'numel', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'X4');
            obj.X4 = val;
        end
        function set.CHORD(obj, val)      %set.CHORD     
            validateattributes(val, {'numeric'}, {'column', 'numel', 2, ...
                'nonnegative', 'finite', 'real', 'nonnan'}, class(obj), 'CHORD');
            obj.CHORD = val;
        end
        function set.AeroProp(obj, val)   %set.AeroProp  
            %set.AeroProp Set method for the property 'AeroProp'. 
            %
            % 'AeroProp' must be a valid instance of the 'awi.fe.AeroProp'
            % class.
            validateattributes(val, {'awi.fe.AeroProp'}, {'scalar', ...
                'nonempty'}, class(obj), 'AeroProp');
            obj.AeroProp = val;
        end
        function val = get.EID(obj)       %get.EID       
           %get.EID Get method for the property 'EID'. 
           val = obj.ID;
        end
        function val = get.PID(obj)       %get.PID       
            %get.PID Get method for the property 'PID'.
            %
            % If the object has been assigned a handle to its 
            % 'awi.fe.AeroProp' object then always use their ID number, 
            % else use PID.
            if isempty(obj.AeroProp)
                val = obj.PID;
            else
                val = obj.AeroProp.ID;
            end
        end
        function val = get.NumPanels(obj) %get.NumPanels 
            
            if isempty(obj.NCHORD) || isempty(obj.NSPAN)
                val = [];
                return
            end
            
            val = obj.NCHORD .* obj.NSPAN;
            
        end
        function val = get.PanelSetCornerCoords(obj) %get.PanelSetCoords
            %get.PanelSetCornerCoords Get method for the dependent property
            %'PanelSetCornerCoords'.
            %
            % Returns the corner coordinates of each panel set. This data
            % is returned in a structure to reduce multiple calls to a
            % dependent variable.
            %
            % Each (x,y,z) data is returned in a matrix of size [2,2,nSets]
            
            %Pass it on...
            val = getPanelCornerCoords(obj);                     
            
        end
    end
    
    methods % construction
        function obj = AeroPanel
            
            %Make a note of the property names
            addFEProp(obj, 'EID', 'PID', 'NCHORD', 'NSPAN', 'X1', ...
                'X4', 'CHORD','AoA','AoA_eta');
            
        end
    end
    
    methods % defining the individual panels 
        function CornerCoords = getPanelCornerCoords(obj)
            %getPanelCornerCoords Returns the corner coordinates of the
            %various panel objects 'obj' in a MATLAB structure with 'X',
            %'Y', 'Z' fields. 
            %           
            % Each (x,y,z) data is returned in a matrix of size [2,2,nObj]
            % as this enables the plotting of the panels to be vectorised.
            
            %TODO - Check that the panel data has been defined.
            
            %Gather data
            x1 = [obj.X1];
            x4 = [obj.X4];
            c  = [obj.CHORD];
            
            %Construct the coordinates
            CornerCoords.X = [ ...
                permute([x1(1, :) ; x1(1, :) + c(1, :)], [1, 3, 2]), ...
                permute([x4(1, :) ; x4(1, :) + c(2, :)], [1, 3, 2])];
            CornerCoords.Y = [ ...
                permute([x1(2, :) ; x1(2, :)], [1, 3, 2]), ...
                permute([x4(2, :) ; x4(2, :)], [1, 3, 2])];
            CornerCoords.Z = [ ...
                permute([x1(3, :) ; x1(3, :)], [1, 3, 2]), ...
                permute([x4(3, :) ; x4(3, :)], [1, 3, 2])];
            
        end        
        function AeroPanel = definePanels(obj)
            %definePanels Returns a structure containing the following
            %information:
            %
            %  - The coordinates in the global coordinate system of the 
            % verticies of each individual aerodynamic panel in all of the
            % aerodynamic panel sets.
            %
            %  - The centre coordinate in the global coordinate system of 
            % each individual panel in all of the aerodynamic panel sets.
            %
            %  - The three-quarter chord position in the centre of each
            %  panel - collocation point for VLM.
            %
            %  - The ID number of each individual panel in all of the
            % aerodynamic panel sets.
            %
            % * * Having the data in this format allows the plotting of * *
            % * * the aerodynamic panels to be vectorised with a single * *
            % * * call to patch.                                        * *
            %
            % The aerodynamic panel coordinates are ordered starting at the
            % inboard LE section and moving counter clockwise, e.g.
            %
            %     1           LE          4
            %       *--------------------*
            %       |                    |
            %       |                    |
            %       |                    |
            % y = i |                    | y = i + 1
            %       |                    |
            %       |                    |
            %       |                    |
            %       *--------------------*
            %      2          TE          3
            
            %Number of aerodynamic panel sets and number of panels in each
            %set
            nPanelSet = numel(obj);
            nPanels   = [obj.NumPanels];  
            
            %Preallocate
            panel   = cell(nPanelSet, 1); % cell array containing the corner coordinates of the panels
            panelID = cell(nPanelSet, 1); % cell array containing the panel numbers
            panelAoA = cell(nPanelSet, 1); % cell array containing the panel numbers
            
            %Define corner points as a [2,2, nPanelSet] matrix
            panelCorners = getPanelCornerCoords(obj);
            x = panelCorners.X;
            y = panelCorners.Y;
            z = panelCorners.Z;
            
            for iC = 1 : nPanelSet
                %Chordwise panel coordiante
                [X, Y, Z] = i_chordwisePanelCoords(x(:, :, iC), y(:, :, iC), z(:, :, iC), obj(iC).NSPAN);
                %All panel corner coordinates
                [xDat, yDat, zDat] = i_allPanelCoords(X, Y, Z, obj(iC).NCHORD);
                %Define panels [5, nPanel, 3]
                panel{iC} = i_panelVerticies(xDat, yDat, zDat, obj(iC).NCHORD, obj(iC).NSPAN, zeros(5, nPanels(iC), 3));
                %Define panel numbers
                if isempty(obj(iC).EID)
                    continue
                end
                panelID{iC} = obj(iC).EID + (0 : nPanels(iC) - 1)';
                
                % get AoA distribution
                coord = permute(panel{iC}, [2, 1, 3]);
                centre = mean(coord(1 : 4, :, :), 1);
                ys_eta = (centre(1,:,2)-obj(iC).X1(2))./(obj(iC).X4(2)-obj(iC).X1(2));
                if isempty(obj(iC).AoA)
                    panelAoA{iC} = zeros(size(ys_eta(:)));
                else
                    panelAoA{iC} = interp1(obj(iC).AoA_eta,obj(iC).AoA,ys_eta(:));
                end
            end
            
            %Combine all panels into a single matrix
            AeroPanel.coord   = permute(cell2mat(panel), [2, 1, 3]);
            AeroPanel.centre  = mean(AeroPanel.coord(1 : 4, :, :), 1);
            AeroPanel.panelID = cell2mat(panelID);
            AeroPanel.AoA = cell2mat(panelAoA);               
            
            % ... local functions
            function [X, Y, Z] = i_chordwisePanelCoords(x, y, z, nSpanPanels)
                %i_chordwisePanelCoords Defines the (x,y,z) coordinates of 
                %each panel corner. Should return 3 matricies (X,Y,Z) of 
                %size [2, nSpanPanels + 1]
                
                % normalised incremental position of spanwise panels
                dS = linspace(0, 1, nSpanPanels + 1);
                
                % difference in x & y & z across panel
                dX = diff(x, [], 2);
                dY = diff(y, [], 2);
                dZ = diff(z, [], 2);
                
                % chordwise lines -- plotted at intermediate locations therefore new x
                % ... and z values are required!
                X = repmat(x(:, 1), [1, nSpanPanels + 1]) + dX * dS;
                Y = repmat(y(:, 1), [1, nSpanPanels + 1]) + dY * dS;
                Z = repmat(z(:, 1), [1, nSpanPanels + 1]) + dZ * dS;
                
            end
            
            function [xDat, yDat, zDat] = i_allPanelCoords(X, Y, Z, nChord)
                %i_allPanelCoords Defines the coordinates at each vertex of
                %the panel.
                
                xDat = i_linspaceMat(X(1, :), X(2, :), nChord + 1);
                yDat = repmat(Y(1, :), [nChord + 1, 1]);
                zDat = repmat(Z(1, :), [nChord + 1, 1]);
                
                
                function M = i_linspaceMat(d1, d2, n)
                    %linspaceMat Linearly spaced matrix
                    
                    dx      = (d2 - d1) / (n - 1);
                    M       = repmat(dx, n, 1);
                    M(1, :) = d1;
                    M       = cumsum(M, 1);
                    
                end
            end
            
            function panel = i_panelVerticies(xDat, yDat, zDat, nChord, nSpan, panel)
                
                k = 1;                  % counter for the panel ID
                for j = 1 : nSpan       % <-- loop through spanwise points
                    for i = 1 : nChord  % <-- loop through chordwise points
                        % panel x-coordinates
                        panel(1, k, 1) = xDat(i    , j);
                        panel(2, k, 1) = xDat(i + 1, j);
                        panel(3, k, 1) = xDat(i + 1, j + 1);
                        panel(4, k, 1) = xDat(i    , j + 1);
                        panel(5, k, 1) = xDat(i    , j);
                        % panel y-coordinates
                        panel(1, k, 2) = yDat(i    , j);
                        panel(2, k, 2) = yDat(i + 1, j);
                        panel(3, k, 2) = yDat(i + 1, j + 1);
                        panel(4, k, 2) = yDat(i    , j + 1);
                        panel(5, k, 2) = yDat(i    , j);
                        % panel z-coordinates
                        panel(1, k, 3) = zDat(i    , j);
                        panel(2, k, 3) = zDat(i + 1, j);
                        panel(3, k, 3) = zDat(i + 1, j + 1);
                        panel(4, k, 3) = zDat(i    , j + 1);
                        panel(5, k, 3) = zDat(i    , j);
                        % next counter
                        k = k + 1;
                    end
                end
                
                panel = permute(panel, [2, 1, 3]);
                
            end
            
        end
        function AeroMesh = defineAeroMeshGeometry(obj)
           %defineAeroMeshGeometry Calculate the mesh geometry for the VLM
           %aerodynamic method based on the panel geometry data. Returns
           %the following properties:
           %
           %    - Coordinates of the centre of each panel
           %    - Coordinates of the 1/4 and 3/4 position of each panel
           %    - Vortex corner coordinates (horseshoe or ring vortex)
           %
           % The aerodynamic panel coordinates are ordered starting at the
           % inboard LE section and moving counter clockwise, e.g.
           %
           %     1           LE          4
           %       *--------------------*
           %       |                    |
           %       |                    |
           %       |                    |
           % y = i |                    | y = i + 1
           %       |                    |
           %       |                    |
           %       |                    |
           %       *--------------------*
           %      2          TE          3    
           %
           % Detailed description
           %    - Any variable that begins with 'r' is a vector.
           %    - Any variable that begins with 'l' is the length/magnitude
           %      of the vector.
           
           AeroMesh = [];
           
           %Get panel geometry
           AeroPanel = definePanels(obj);
           
           %Get panel vertex coordinates
           coords = AeroPanel.coord(1 : 4, :, :); clear AeroPanel           
           nPanel = size(coords, 2);
           
           %Calculate unit vectors and length for each panel edge
           %    - Vector
           r12 = coords(2, :, :) - coords(1, :, :);           
           r23 = coords(3, :, :) - coords(2, :, :); 
           r34 = coords(4, :, :) - coords(3, :, :); 
           r41 = coords(1, :, :) - coords(4, :, :); 
           %    - Length
           l12 = sqrt(sum(r12.^2, 3));
           l23 = sqrt(sum(r23.^2, 3));
           l34 = sqrt(sum(r34.^2, 3));
           l41 = sqrt(sum(r41.^2, 3));
           %    - Unit vector
           r12_bar = r12 ./ l12;
           r23_bar = r23 ./ l23;
           r34_bar = r34 ./ l34;
           r41_bar = r41 ./ l41;
           
           %Calculate mid-points of each panel [0;1/4;1/2;3/4;1] x chord
           midCoords = zeros(5, nPanel, 3);
           %    - LE & TE
           midLE = coords(4, :, :) + (r41_bar .* l41 ./ 2);
           midTE = coords(2, :, :) + (r23_bar .* l23 ./ 2);
           midCoords(1, :, :)   = midLE;
           midCoords(end, :, :) = midTE;
           %    - Vector along middle of panel
           rMid2Mid  = midTE - midLE;
           lMid2Mid  = sqrt(sum(rMid2Mid.^2, 3));
           rMid2Mid_ = rMid2Mid ./ lMid2Mid;
           %    - 1/4, 1/2, 3/4 chord
           midCoords(2, :, :) = midLE + rMid2Mid_ .* (lMid2Mid ./ 4);
           midCoords(3, :, :) = midLE + rMid2Mid_ .* (lMid2Mid ./ 2);
           midCoords(4, :, :) = midLE + rMid2Mid_ .* (3 .* lMid2Mid ./ 4);
           
           %Plotting (if required...)
           hF     = figure;
           hAx    = axes('Parent', hF, 'NextPlot', 'add');
           hPanel = drawElement(obj, hAx);
           hM = gobjects(1, 5);
           str = {'LE', '1/4 chord', '1/2 chord', '3/4 chord', 'TE'};
           mkr = {'o'; 's'; 'd'; '^'; '*'};
           clr = {'r'; 'g'; 'b'; 'm'; 'k'};
           for i = 1 : 5
              hM(i) = plot3(hAx, ...
                  midCoords(i, :, 1), midCoords(i, :, 2), midCoords(i, :, 3), ...
                  'LineStyle', 'none');
           end
           set(hM, ...
               'MarkerEdgeColor'  , 'k', ...
               {'MarkerFaceColor'}, clr, ...
               {'Marker'}         , mkr);
           legend(hAx, [hPanel, hM], [{'AeroPanels'}, str]);           
           for i = 1 : size(coords, 1)
                text(coords(i, :, 1), coords(i, :, 2), coords(i, :, 3), num2str(i));
           end
           
        end
    end
    
    methods % visualisation
        function [hg, varargout] = drawElement(obj, ha, varargin)
            %drawElement Draws the AeroPanel objects as a set of patches
            %in 3D space. A patch is used so that this function returns
            %only a single graphics handle for all the aerodynamic panels
            %in the collection.
            %
            % Accepts an array of objects.
            
            %Sort the 'AeroPanel' objects in ascending order of ID numbers
            %   - Makes it easier to assign results data later on
            id = [obj.ID];
            if ~isempty(id)
                [~, index] = sort(id);
                obj = obj(index);
            end
            
            %Define the aerodynamic panel coordinates
            AeroPanelData = definePanels(obj);
            
            %Plot! - Attach the panel ID in the user data
            hg = patch(ha, ...
                AeroPanelData.coord(:, :, 1), ...
                AeroPanelData.coord(:, :, 2), ...
                AeroPanelData.coord(:, :, 3), ...
                'k'       , 'FaceColor', [0.88, 0.88, 0.88] , ... %gray
                'Tag'     , 'Aerodynamic Panels', ...
                'UserData', {'PanelID', AeroPanelData.panelID}, ...
                'SelectionHighlight', 'off');  
            
            %Any additional arguments?
            if ~isempty(varargin)
                tok = varargin(1 : 2 : end);
                val = varargin(2 : 2 : end);
                set(hg, tok, val);
            end
            
            %Could define the colour data at each vertex to match the
            %defined AoA distribution?
            
            if nargout > 1
                varargout{1} = AeroPanelData;
            end
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.AeroPanel' object
            %into a text file using the format of the MSC.Nastran 'CAERO1'
            %bulk data entry.
            %
            % The following assumptions are made:
            %   * The coordinate system for locating the corner points of
            %   the panel is assumed to be the basic coordinate system.
            %   i.e. CP = 0.
            %   * The panels are assumed to have a constant mesh
            %   distribution. i.e. NSPAN & NCHORD are used instead of LSPAN
            %   and LCHORD.
            %   * The interference group ID number (IGID) is assumed to be
            %   the same for all CAERO1 panels.
            
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
                comment = ['CAERO1 : Defines an aerodynamic macro ', ...
                    'element (panel) in terms of two leading edge ', ...
                    'locations and side chords.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
            
            %How many objects?
            nObj = numel(obj);
                        
            %Split up the coordinates & chord values        
            coords1 = [obj.X1];
            coords4 = [obj.X4];
            chords  = [obj.CHORD];
            X1_ = num2cell(coords1(1, :));
            Y1_ = num2cell(coords1(2, :));
            Z1_ = num2cell(coords1(3, :));
            X12 = num2cell(chords(1, :));
            X4_ = num2cell(coords4(1, :));
            Y4_ = num2cell(coords4(2, :));
            Z4_ = num2cell(coords4(3, :));
            X43 = num2cell(chords(2, :));
            
            %Card name
            nam   = repmat({'CAERO1'}  , [1, nObj]);
            blnks = repmat({blanks(8)}, [1, nObj]);
            zrs   = num2cell(zeros(1, nObj));
            onez  = num2cell(ones(1, nObj));
            
            %Set up the format for printing
            data = [ ...
                nam   ; {obj.ID} ; {obj.PID} ; zrs ; {obj.NSPAN} ; {obj.NCHORD} ; blnks ; blnks ; onez ; ...
                blnks ; X1_      ; Y1_       ; Z1_ ; X12         ; X4_          ; Y4_   ; Z4_   ; X43];
            
            %Write in 16-character column width as standard
%             format = [ ...
%                 '%-8s%-8i%-8i%-8i%-8i%-8i%-8s%-8s%-8i\r\n', ...
%                 '%-8s%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g%#-8.3g\r\n'];

            format = [ ...
                '%-8s%-8i%-8i%-8i%-8i%-8i%-8s%-8s%-8i\r\n', ...
                '%-8s%#-8.5g%#-8.4g%#-8.5g%#-8.5g%#-8.5g%#-8.4g%#-8.5g%#-8.5g\r\n'];
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            %Write the W2GJ file...
            aoa = {obj.AoA};
            if ~any(cellfun(@isempty, aoa))
                
                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end

    methods (Static)
        function n = calculateFaceNormals(hP)
            %calculateFaceNormals Calculates the normal vector for each
            %patch in the patch object.            
            %
            %https://uk.mathworks.com/matlabcentral/fileexchange/24330-patch-normals
            
            validateattributes(hP, {'matlab.graphics.primitive.Patch'}, ...
                {'scalar'}, 'calcualteFaceNormals', 'hP');
            
            %Face corners index
            ia = hP.Faces(:, 1);
            ib = hP.Faces(:, 2);
            ic = hP.Faces(:, 3);
            
            %Grab the vertices
            v = hP.Vertices;
            
            %Face normals (area weighted)
            v1 = v(ia, :) - v(ib,:);
            v2 = v(ic,:) - v(ia,:);
            n  = cross(v1, v2);
            
            %Unit vector
            n = n ./ sqrt(n(:, 1).^2 + n(:, 2).^2 + n(:, 3).^2);
            
            %             %vertice normals
            %             N = zeros(size(FV.vertices)); %init vertix normals
            %             for i = 1:size(FV.faces,1) %step through faces (a vertex can be reference any number of times)
            %                 N(A(i),:) = N(A(i),:)+n(i,:); %sum face normals
            %                 N(B(i),:) = N(B(i),:)+n(i,:);
            %                 N(C(i),:) = N(C(i),:)+n(i,:);
            %             end
            
        end
    end
    
end

