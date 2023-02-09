classdef CrossSection < awi.model.Component & awi.mixin.BeamPropertyObject
    %CrossSection Defines a generic planar cross-section that can be
    %rotated and plotted in 3-dimensions using the 'Orientation' property.
    %
    % TODO - Add cosine-spacing for the circle and ellipse
    
    %Basic geometry parameters
    properties
        %X-Coordinates of the cross-section in a 2D-plane
        X
        %Z-Coordinates of the cross-section in a 2D-plane
        Z
        %Rotation matrix describing the orientation of the plane in the
        %global coordinate system. The x- and z-axes lie in-plane and the
        %y-axes is normal to the plane.
        Orientation
    end
    
    %Using cross-section library
    properties
        %Name of present property
        CrossSectionName = '';
        %Number of points in the x-direction for the default cross-section
        NumPointsX = 50;
    end
    
    %Circle properties
    properties
        Radius = nan;
    end
    
    %Ellipse properties
    properties
        %Length of the major axis (a)
        MajorAxis = nan;
        %Length of the minor axis (b)
        MinorAxis = nan;
    end
    
    %Helper properties
    properties (Constant)
        %Name of valid cross-sections
        CrossSectionLibrary = {'Circle', 'Ellipse'};
        %Function handles of methods for defining coordinates
        CrossSectionMethods = {@getCircleCoords, @getEllipseCoords};
        %Name of properties that relate to specific pre-set cross-sections
        CrossSectionProps = {{'Radius'}, {'MajorAxis', 'MinorAxis'}};
    end
    
    methods % set / get
        function set.X(obj, val)                %set.X           
            %set.X Set method for the property 'X'.
            %
            % 'X' must be a row vector
            
            validateattributes(val, {'numeric'}, {'row', 'real'}, ...
                class(obj), 'X');
            obj.X = val;
        end
        function set.Z(obj, val)                %set.Z           
            %set.Z Set method for the property 'Z'.
            %
            % 'Z' must be a row vector
            
            validateattributes(val, {'numeric'}, {'row', 'real'},  ...
                class(obj), 'Z');
            obj.Z = val;
        end
        function set.Orientation(obj, val)      %set.Orientation 
            %set.Orientation Set method for the property 'Orientation'.
            %
            % 'Orientation' must be a valid instance of the
            % 'awi.model.CoordinateSystem' class.
            
            validateattributes(val, {'awi.model.CoordSys'}, ...
                {'scalar', 'nonempty'}, class(obj), 'Orientation');
            obj.Orientation = val;
            
        end    
        function set.CrossSectionName(obj, val) %set.CrossSectionName
            %set.CrossSectionNames Set method for the property
            %'CrossSectionName'.
            
            validatestring(val, obj.CrossSectionLibrary, class(obj), ...
                'CrossSectionName');
            obj.CrossSectionName = val;
            
        end
        function set.NumPointsX(obj, val)       %set.NumPointsX
            %set.NumPointsX Set method for the property 'NumPointsX'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnan', 'finite', 'real', 'nonempty'}      , ...
                class(obj), 'NumPointsX');
            obj.NumPointsX = val;
        end
        function set.Radius(obj, val)           %set.Radius 
            %set.Radius Set method for the property 'Radius'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'Radius');
            obj.Radius = val;
            
        end
        function set.MajorAxis(obj, val)        %set.MajorAxis 
            %set.MajorAxis Set method for the property 'MajorAxis'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'MajorAxis');
            obj.MajorAxis = val;
            
        end
        function set.MinorAxis(obj, val)        %set.MinorAxis 
            %set.MinorAxis Set method for the property 'MinorAxis'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'MinorAxis');
            obj.MinorAxis = val;
            
        end
    end
    
    methods % construction
        function obj = CrossSection
            
            %Need to populate the 'Orientation' property
            obj.Orientation = awi.model.CoordSys;
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht, varargin)
            %drawElement Draw method for the 'awi.model.CrossSection'
            %class.
            
            hg = [];
            
            p = inputParser;
            addParameter(p, 'FillCrossSection', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %Get the 3D coordinates of the cross-section
            [x, y, z] = calculateGlobalCoords(obj);
            
            if isempty(x) %Escape route
                return
            end
            
            nObj = numel(obj);
            
            %Pad profile coordinates with NaN so we can vectorise the plot
            xP = [x' ; nan(1, nObj)];
            yP = [y' ; nan(1, nObj)];
            zP = [z' ; nan(1, nObj)];
            
            %Plot it - TODO: Replace with vectorised call to 'patch'
            hg = surf(x, y, z, ...
                'Parent'   , ht    , ...
                'EdgeColor', 'none', ...
                'FaceColor', [8, 181, 48] ./ 255, ...
                'FaceAlpha', 0.25  , ...
                'Tag'      , 'Cross-Section');
 
            %Profiles
            hg(end + 1) = plot3(xP(:), yP(:), zP(:), 'Parent', ht, ...
                'LineStyle', '-', ...
                'Color'    , 'k', ...
                'Tag', 'Cross-Section Profiles');
            
            if p.Results.FillCrossSection
                %Vectorised patch
                [nCS, nPoints] = size(x);
                ub = cumsum(repmat(nPoints, [1, nCS]));
                lb = [1, ub(1 : end - 1) + 1];
                faces = arrayfun(@(i) lb(i) : ub(i), 1 : nCS, 'Unif', false);
                faces = vertcat(faces{:});                
                x      = num2cell(x', 1);
                y      = num2cell(y', 1);
                z      = num2cell(z', 1);
                vertex = [vertcat(x{:}), vertcat(y{:}), vertcat(z{:})];                
                %Draw as a single patch object
                hg(end + 1) = patch(ht, 'Faces', faces, 'Vertices', vertex, 'FaceColor', 'b');
            end
            
            %Coordinate systems - Likely to plot in the incorrect place as
            %the Orientation CoordSys objects do not have knowledge of the
            %beam they are attached to.
%             hg_ = drawElement([obj.Orientation], ht, 'Cross-Section Profile Orienation');
        end
    end    
   
    methods % public facing methods
        %Generating 2D cross-section using library
        function generateCoordsFromLibrary(obj)
            %generateCoordsFromLibrary Generates the 2D cross-section
            %coordinates using the pre-defined library of cross-sections.
                        
            %Filter the cross-section names
            name  = {obj.CrossSectionName};
            idx = cellfun(@isempty, name);
            if all(idx) %Escape route
                return
            end
            name = name(~idx);
            obj  = obj(~idx);
            
            %Loop through unique cross-section names and calculate coords
            uName = unique(name);            
            for iCS = 1 : numel(uName)
                
                %Down-select
                CS = obj(ismember(name, uName{iCS}));
                
                %Get prop names and method handle for this cross-section
                idx  = ismember(CS(1).CrossSectionLibrary, uName);
                prps = CS(1).CrossSectionProps{idx};
                fn   = CS(1).CrossSectionMethods{idx};
                
                %Check correct properties have been defined
                val  = get(CS, prps);
                
                if any(any(cellfun(@isempty, val)))
                    continue
                end
                
                %Update the properties
                fn(CS);                
                
            end
            
        end
        function generateCoordsFromNACA(obj, naca)
            %generateCoordsFromNACA Creates the (x,z) coordinates for the
            %CrossSection using the NACA 4-series or reading the coordinate
            %data from a file.
            %
            % Inputs:
            %   - 'naca': A cell-string of NACA aerofoil names or filepaths
            %             pointing to a valid coordinate file.
            % Example:
            % >> obj(1) = awi.model.CrossSection;
            % >> obj(2) = awi.model.CrossSection;
            % >> obj(3) = awi.model.CrossSection;
            % >> generateCoordsFromNACA(obj, {'NACA0012', 'NACA0015',
            % 'mycoordinatefile.txt})
            
            nObj = numel(obj);
            validateattributes(naca, {'cell'}, {'nonempty'}, class(obj), 'naca');
            assert(iscellstr(naca), ['Expected the profiles to be ', ...
                'defined in a cell array of characters (cellstr).']); %#ok<ISCLSTR>
            assert(nObj == numel(naca), ['Expected the number ', ...
                'of objects to match the number of NACA aerofoils.']);
            
            %What have we got?
            idxNACA = contains(naca, 'NACA');
            idxFile = contains(naca, '.');  %File is denoted by '.'
            idxNACA(idxFile) = false;       %Override NACA if '.' is present
            if ~any(or(idxFile, idxNACA))
                error(['Unable to parse NACA aerofoil codes or filenames ', ...
                    'in the aerofoil list\n\n%s\n\nConsult the AWI ', ...
                    'documentation for help.'], strjoin(naca, '\n'));
            end
            indNACA = find(idxNACA);
            indFile = find(~idxNACA);
            
            %How many points?
            nP = [obj.NumPointsX];
            if range(nP) > 0
                warning(['Ambiguous number of points for cross-section ' , ...
                    '(x,z) coordinates -> assuming the finest available ', ...
                    'distribution of points.']);
                nP = max(nP);
            else
                nP = nP(1);
            end
            
            %Preallocate - 2 x nPoints accounts for upper and lower surface
            x = zeros(nObj, nP * 2); 
            z = zeros(nObj, nP * 2);
            x_ = cell(nObj, 1);
            z_ = cell(nObj, 1);
            
            %Assign the (x,z) coordinate data
            for iCS = 1 : numel(indFile)
                ii = indFile(iCS);
                [x_{ii}, z_{ii}] = loadCrossSectionFromFile(naca{ii});      
            end
            %Check for compatibility with other files
            if ~isempty(indFile)
                nPoints = cellfun(@numel, x_);
                assert(range(nPoints) == 0, ['Expected the cross-sections ', ...
                    'to be defined using the same number of points. Check ', ...
                    'the input files and try again.']);
                nPoints = nPoints(1) / 2;
                if nPoints ~= nP
                    nP = nPoints(1);
                    x  = zeros(nObj, nP * 2);
                    z  = zeros(nObj, nP * 2);
                    set(obj, 'NumPointsX', nP);
                end
                x(indFile, :) = vertcat(x_{:});
                z(indFile, :) = vertcat(z_{:});
            end
            for iCS = 1 : numel(indNACA)
                ii = indNACA(iCS);
                [x(ii, :), z(ii, :)] = calculateNACACoords(naca{ii}, nP);
            end
            
            %Check that the x-vectors are the same, if not then interpolate
            xUpper = x(:, 1 : nP);
            xLower = x(:, nP + 1 : end);
            zUpper = z(:, 1 : nP);
            zLower = z(:, nP + 1 : end);   
            if or( ...
                    any(any(abs(diff(xUpper, [], 1)) > 1e-10)), ...
                    any(any(abs(diff(xLower, [], 1)) > 1e-10)))
                
                %Use cosine spacing
                x_c  = cosspace(0, 1, nP, 'bRow', true);
                x_c_ = fliplr(x_c);
                
                indU = find(~all((xUpper - x_c), 2));
                indL = find(~all((xLower - x_c_), 2));    
                for ii = 1 : numel(indU)
                    zUpper(indU(ii), :) = interp1(xUpper(indU(ii), :), zUpper(indU(ii), :), x_c, 'spline');
                end
                for ii = 1 : numel(indL)
                    zLower(indL(ii), :) = interp1(xLower(indL(ii), :), zLower(indL(ii), :), x_c_, 'spline');
                end                
                xUpper(indU, :) = repmat(x_c , [numel(indU), 1]);
                xLower(indL, :) = repmat(x_c_, [numel(indU), 1]);
                x = [xUpper, xLower];
                z = [zUpper, zLower];
                warning(['Each cross-section must be defined over the '  , ...
                    'same normalised set of points, otherwise it is not ', ...
                    'possible to interpolate the profiles. The '         , ...
                    'x-vectors for some profiles have been adjusted to ' , ...
                    'use a cossine spacing and the profiles have been '  , ... 
                    'interpolated to this new set of points using a '    , ...
                    '''spline'' interpolation.']);                
            end
            
            set(obj, {'X'}, num2cell(x, 2), {'Z'}, num2cell(z, 2), {'Name'}, naca(:));

        end
        %Calculating 3D coordinates of cross-section 
        function [x, y, z] = calculateGlobalCoords(obj)
            %calculateGlobalCoords Calculates the (x,y,z) coordinates of
            %the cross-sections, accounting for the position of the
            %cross-section along beam and the orientation of the
            %cross-section.
            %
            % The resulting coordinates are expressed in the beam
            % coordinate system, i.e. with respect to the root of the beam.
            
            x = [];
            y = [];
            z = [];
            
            %Check we only have one parent beam
            pb = [obj.BeamHandle];
            if isempty(pb)
                return
            end
%             assert(numel(unique(pb)) == 1, ['Cannot handle the case ', ...
%                 'where we have cross sections from multiple beams. ' , ...
%                 'Update code.']); %TODO - Update the 'sort' method for
%                 all MVC objects as it calls the custom sort method during
%                 a call to unique!!!
            
            nObj = numel(obj);
            
            %For objects of class 'BeamPropertyObject' the 'eta_flag' is
            %always set to 'R'. i.e. Along the axis of the line -> Get the
            %length of the line then calculate the (x,y,z) coordinates of
            %the object along the beam
            [xd, yd, zd] = xyzdata(pb(1));
            if isempty(xd)
                return
            end
            r = awi.model.Stick.getLineLength(xd, yd, zd);
            eta  = r ./ r(end);
            etaC = [obj.BeamEta];
            
            %Can do it in one call to 'interp1'
            xyzC = interp1(eta, [xd ; yd ; zd]', etaC);
            
            %Get the coordinates of the cross-sections
            %   - Note: The beam coordinate system assumes that the x-axis
            %     lies along the beam, so we need to swap the cross-section
            %     x-axis with the beam y-axis.
            %   - TODO: Define cross-section coordinate to match (x,y) of
            %     beam
            y = vertcat(obj.X);
            z = vertcat(obj.Z);
            x = zeros(size(y));
            
            if isempty(x) || isempty(z) %Escape route
                return
            end
            
            %Rotate through the 'Orientation' rotation matrix
            RMat   = get([obj.Orientation], {'RMatrix'});
            coords = arrayfun(@(i) ...
                RMat{i} * [x(i, :) ; y(i, :) ; z(i, :)], 1 : nObj, 'Unif', false);
            coords = cat(3, coords{:});
            coords = permute(coords, [3, 2, 1]);
            x = coords(:, :, 1);
            y = coords(:, :, 2);
            z = coords(:, :, 3);            
            
            %Add the offsets along the beam
            x = x + xyzC(:, 1);
            y = y + xyzC(:, 2);
            z = z + xyzC(:, 3);
            
        end
        %Interpolating cross-section objects
        function NewObj = interpolateCrossSection(obj, etaQ)
            %interpolateCrossSection Interpolates the cross-section data at
            %new eta position 'etaQ' and returns a new set of
            %'awi.model.CrossSection' objects.
            
            NewObj = [];
            
            [xU, xL, zU, zL, ~] = interpolateXZData(obj, etaQ);
            rMatrix             = interpolateOrientationData(obj, etaQ);
            
            x = [xU, fliplr(xL)];
            z = [zU, fliplr(zL)];
            
            NewObj    = arrayfun(@(~) awi.model.CrossSection, 1 : numel(etaQ));
            NewOrient = arrayfun(@(~) awi.model.CoordSys    , 1 : numel(etaQ));
            
            set(NewObj, {'X'}          , num2cell(x, 2));
            set(NewObj, {'Z'}          , num2cell(z, 2));
            set(NewObj, 'BeamHandle'   , obj(1).BeamHandle);
            set(NewObj, {'BeamEta'}    , num2cell(etaQ)');
            set(NewObj, {'Orientation'}, num2cell(NewOrient)');
            set(NewOrient, {'RMatrix'} , rMatrix');
            
        end
        %Interpolating cross-section coordinates
        function [xU, xL, zU, zL, etaQ2D] = interpolateXZData(obj, etaQ)
            %interpolateXZData Interpolates the (x,z) data of the 
            %CrossSection objects at positions 'etaQ' and returns the 
            %FULL-SCALE coordinates of the profiles.
            %
            % % See also: interpolateProfileData normaliseProfileData scaleProfileData
            
            %Interpolate the coordinates to the new eta positions
            [xU, xL, zU, zL, etaQ2D, xzData] = interpolateNormXZData(obj, etaQ);
            
            %Scale the profiles back to full scale
            [xU, zU] = scaleProfileData(xU, zU, xzData);
            [xL, zL] = scaleProfileData(xL, zL, xzData);
            
        end
        function [xU, xL, zU, zL, etaQ2D, xzData] = interpolateNormXZData(obj, etaQ)
            %interpolateNormXZData Interpolates the (x,z) data of the 
            %CrossSection objects at positions 'etaQ' and returns the 
            %NORMALISED coordinates of the profiles.
            %
            % See also: interpolateProfileData, normaliseProfileData
                        
            parse(obj);
            validateattributes(etaQ, {'numeric'}, {'vector', 'nonnegative', ...
                '<=', 1}, class(obj), 'etaQ');
            
            %Grab data
            x    = vertcat(obj.X);
            z    = vertcat(obj.Z);
            eta  = horzcat(obj.BeamEta);
            n    = size(x, 2) / 2;
            
            %Normalise
            [xNorm, zNorm, xzDataCS] = normaliseProfileData(x, z);
            
            %Interpolate
            [xU, xL, zU, zL, etaQ2D] = interpolateProfileData(xNorm, zNorm, eta, etaQ, n);
            xzData                   = interp1(eta, xzDataCS, etaQ);
            
        end
        %Normalising the cross-section profile
        function [xNorm, zNorm, xzData] = normaliseXZData(obj)
            %normaliseXZData Returns the normalised and centred (x,z)
            %coordinates for the CrossSection objects as well as a vector
            %containing the scaling valus and centre-offset values. 
            %
            % See also: normaliseProfileData
            
            parse(obj);
            
            x = vertcat(obj.X);
            z = vertcat(obj.Z);
            
            [xNorm, zNorm, xzData] = normaliseProfileData(x, z);            
                           
        end
        %Interpolating the cross-section orientations
        function rMatrix = interpolateOrientationData(obj, etaQ)
            %interpolateXZData Interpolates the rotation matrix data of the 
            %CrossSection objects at positions 'etaQ' and returns a set of 
            %N rotation matrices where N = numel(etaQ).
            
            rMatrix = [];            
            Orient  = get(obj, {'Orientation'});            
            if any(cellfun(@isempty, Orient))
                return
            end
            
            parse(obj);
            
            Orient = horzcat(Orient{:});
            rData  = [Orient.RMatrix_];
            eta    = [obj.BeamEta];
            
            rMatrix = interp1(eta, rData', etaQ, 'linear', 'extrap')';
            rMatrix = arrayfun(@(i) reshape(rMatrix(:, i), [3, 3]), 1 : numel(etaQ), 'Unif', false);
            
        end
    end
            
    methods (Access = protected) %  calculating cross-section coordinates
        %Methods for generating cross-sections from the library
        function getCircleCoords(obj)
            %getCircleCoords Calculates the coordinates of a circle using
            %the radius.
            %
            % Eqn of a circle: x^2 + z^2 = r^2
            
            r = vertcat(obj.Radius);
            
            %Must have consistent number of points in X direction
            nx = obj(1).NumPointsX * 2;
            th = linspace(0, 2 * pi, nx);
            
            %Use polar coordinates
            x = r * cos(th);
            z = r * sin(th);
            
            %Assign to the object
            set(obj, {'X'}, num2cell(x, 2));
            set(obj, {'Z'}, num2cell(z, 2));
            
        end
        function getEllipseCoords(obj)
            %getEllipseCoords Calculates the coordinates of an ellipse
            %using the lengths of the major and minor axis.
            %
            % Eqn of an ellipse: x^2 / a^2 + z^2 / b^2 = 1
           
            %Gather data
            a = vertcat(obj.MajorAxis) ./ 2;
            b = vertcat(obj.MinorAxis) ./ 2;
            
            %Must have consistent number of points in X direction
            nx = obj.NumPointsX;
            x  = arrayfun(@(i) linspace(-a(i), a(i), nx), 1 : numel(a), 'Unif', false);
            x  = vertcat(x{:});
            
            %Calculate y-coorindates
            z  = sqrt((b.^2 .* (a.^2 - x.^2)) ./ a.^2);
            
            %Apply symmetry
            x = [x, fliplr(x)];
            z = [z, -fliplr(z)];
            
            %Assign to the object
            set(obj, {'X'}, num2cell(x, 2));
            set(obj, {'Z'}, num2cell(z, 2));
            
        end
        %Parsing the object
        function parse(obj)
            %parse Asserts that each CrossSection object has the required
            %(x,z) data defined.
            
            xData = get(obj, {'X'});
            zData = get(obj, {'Z'});
            eta   = get(obj, {'BeamEta'});
            assert( ...
                all(and(~cellfun(@isempty, xData), ~cellfun(@isempty, zData))), ...
                ['The CrossSection objects must have the ''X'' and '   , ...
                '''Z'' properties defined in order to interpolate the ', ...
                'profile data.']);
            assert(all(~cellfun(@isempty, eta)), ['Each CrossSection ', ...
                'object must have the ''BeamEta'' property defined. ' , ...
                'Ensure that all CrossSection objects are associated ', ...
                'with a Beam object and try again.']);
            nData = [cellfun(@numel, xData) ; cellfun(@numel, zData)];
            assert(range(nData(:)) == 0, ['The number of points '        , ...
                'defining the ''X'' and ''Z'' data for the CrossSection ', ...
                'object must be the same in order to interpolate the '   , ...
                'profile data.']);
            xData = vertcat(xData{:});
            zData = vertcat(zData{:});
            if mod(nData(1), 2) ~= 0
                bEqual = and( ... 
                    isequal(xData(:, 1), xData(:, end)), ...
                    isequal(zData(:, 1), zData(:, end)));
                assert(bEqual, ['Expected the (x,z) data points in the ' , ...
                    'CrossSection objects to describe a continuous loop ', ...
                    'around the perimeter. When an odd number of data '  , ...
                    'points are provided the first and last data points ', ...
                    'must be equal.']);
            end
            
        end
    end
    
end

