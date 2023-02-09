classdef Compartment < mvc.mixin.Nameable
    %Compartment Defines a compartment within a generic volume defined by a
    %beam and some cross-sections.
    %
    % TODO - The payload is assume to fill in the local Z-coordinate
    % system. This means there is no consideration of how gravity would
    % effect the distribution of the payload, say for example, if the
    % payload was a liquid - such as for fuel loading.        
    %
    % TODO - Do not inherit from Nameable! The 'addPropertyGroup' method
    % causes massive overheads when doing the golden section.
    
    %Defining the compartment geometry
    properties
        %Normalised position along the beam where the compartment is
        %defined. N.B. This 'eta' refers to the straightline distance along
        %the parent beam (i.e. eta_flag = 'R')
        EtaLocations = [0 ; 1];
        %Normalised positon along the local cross-section x-axes where the
        %compartment starts 
        LocalXStart = 0;
        %Normalised positon along the local cross-section x-axes where the
        %compartment ends
        LocalXEnd = 1;
        %Normalised positon along the local cross-section z-axes where the
        %compartment starts
        LocalZStart = 0;
        %Normalised positon along the local cross-section z-axes where the
        %compartment ends
        LocalZEnd = 1;
        %Number of points on the upper/lower surfaces
        NumPoints = 50;
    end
    
    %Defining the compartment mass
    properties
        %Fraction of compartment which is filled with the payload
        PayloadFraction = 0;
        %Density of the contents of the compartment [kg/m^3]
        MassDensity = 0;
    end
    
    %Associating with other objects
    properties
        %Handle to the 'awi.model.Beam' object that this compartment is
        %attached to 
        %TODO - Access to this should be protected somehow
        BeamHandle
    end
    
    %Visualising the compartment
    properties
        %Sets the 'FaceColor' property of the compartment
        CompartmentColour = [39, 126, 219] ./ 255;
        %Sets the 'FaceColor' property of the payload
        PayloadColour = [123, 184, 11] ./ 255;
    end
    
    %3D geometry
    properties (SetAccess = private)
        %Vector of 'awi.model.CrossSection' objects that define the
        %compartment
        CompartmentCrossSections = awi.model.Compartment.empty;
        %Vector of 'awi.model.CrossSection' objects that define the payload
        PayloadCrossSections = awi.model.Compartment.empty;
        %Coordinates of the compartment (parent coordinate system)
%         CompartmentCoords
        %Coordinates of the payload within the compartment
%         PayloadCoords
        %Vector of 'awi.model.CrossSection' objects that define this
        %compartment
        %Vector of 'awi.model.CrossSection' objects that define the payload
        %Total volume of the compartment
        Volume
        %Total volume of the payload
        PayloadVolume
        %Triangular mesh of the compartment geometry (parent coordinate
        %system)
        CompartmentMesh
        %Triangular mesh of the payload geometry (parent coordinate system)
        PayloadMesh
    end
    
    %Mass properties
    properties (SetAccess = private)
        %Total mass of the payload within the compartment
        PayloadMass
        %(x,y,z) position of the centre of gravity of the payload
        PayloadCoG
        %3 x 3 inertia tensor of the payload
        PayloadInertia
        %Compartment rigid body geometric properties as calculated by
        %'RigidBodyParams' for each section in the compartment
        CompartmentRigidBodyPropsBySection
        %Payload rigid body geometric properties as calculated by
        %'RigidBodyParams' for each section in the comparment
        PayloadRigidBodyPropsBySection
        %Payload mass properties split by section
        PayloadRigidBodyMassPropsBySection        
    end        
        
    methods % set / get
        function set.EtaLocations(obj, val)    %set.EtaLocations    
            %set.EtaLocations Set method for the property 'EtaLocations'.
            
            if isrow(val)
                val = val';
            end
            validateattributes(val, {'numeric'}, {'column', 'nonempty', ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1    , ...
                'increasing'}, class(obj), 'EtaLocations');
            obj.EtaLocations = val;
        end
        function set.LocalXStart(obj, val)     %set.LocalXStart     
            %set.LocalXStart Set method for the property 'LocalXStart'.
            
            validateattributes(val, {'numeric'}, {'column', 'nonempty', ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'LocalXStart');
            obj.LocalXStart = val;
        end
        function set.LocalXEnd(obj, val)       %set.LocalXEnd       
            %set.LocalXEnd Set method for the property 'LocalXEnd'.
            
            validateattributes(val, {'numeric'}, {'column', 'nonempty', ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'LocalXEnd');
            obj.LocalXEnd = val;
        end
        function set.LocalZStart(obj, val)     %set.LocalZStart     
            %set.LocalZStart Set method for the property 'LocalZStart'.
            
            validateattributes(val, {'numeric'}, {'column', 'nonempty', ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'LocalZStart');
            obj.LocalZStart = val;
        end
        function set.LocalZEnd(obj, val)       %set.LocalZEnd       
            %set.LocalZEnd Set method for the property 'LocalXEnd'.
            
            validateattributes(val, {'numeric'}, {'column', 'nonempty', ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'LocalZEnd');
            obj.LocalZEnd = val;
        end
        function set.NumPoints(obj, val)       %set.NumPoints       
            %set.NumPoints Set method for the property 'NumPoints'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
                class(obj), 'NumPoints');
            obj.NumPoints = val;
        end
        function set.PayloadFraction(obj, val) %set.PayloadFraction 
            %set.PayloadFraction Set method for the property
            %'PayloadFraction'.
            
            validateattributes(val, {'numeric'}, {'scalar','finite'  , ...
                'real', 'nonnan', 'nonempty', 'nonnegative', '<=', 1}, ...
                class(obj), 'PayloadFraction');
            obj.PayloadFraction = val;
            
        end
        function set.MassDensity(obj, val)     %set.MassDensity     
            %set.MassDensity Set method for the property 'MassDensity'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real', 'nonnan'}, class(obj), 'MassDensity');
            obj.MassDensity = val;
        end
        function set.BeamHandle(obj, val)      %set.BeamHandle      
            %set.BeamHandle Set method for the property 'BeamHandle'.
            
            validateattributes(val,{'awi.model.Beam'}, {'scalar', ...
                'nonempty'}, class(obj), 'BeamHandle');
            obj.BeamHandle = val;            
        end
    end    
   
    methods (Sealed) % visualisation
        function hg = drawElement(obj, ht, varargin)
            %drawElement Draw method for the 'awi.model.Compartment'
            %class.
            
            hg = [];
            
            if numel(obj) > 1
               hg = arrayfun(@(o) drawElement(o, ht, varargin{:}), obj, 'Unif', false); 
               hg = horzcat(hg{:});
               return
            end
            
            %Parse input
            prp = varargin(1 : 2 : end);
            idx = cellfun(@iscell, prp);
            prp(idx) = cellfun(@(x) x{1}, prp(idx), 'Unif', false);
            idx = ismember(prp, 'DrawPayload');
            if any(idx)
               ind = find(idx) * 2;
               val = varargin{ind};
               if islogical(val)
                   drawPL = val;
               elseif ischar
                   switch val
                       case 'yes'
                           drawPL = true;
                       case 'no'
                           drawPL = false;
                       otherwise
                           validatestring(val, {'yes', 'no'}, 'drawElement', 'DrawPayload');
                   end
               end
               varargin([ind - 1, ind]) = [];
            else
                drawPL = true;
            end
            
            %Draw the compartment
            if isempty(obj.CompartmentCrossSections)
                hg{1} = [];
            else
                hg{1} = drawElement(obj.CompartmentCrossSections, ht);
            end
            if ~isempty(hg{1})
                set(hg{1}, 'Tag', 'Compartment');
                set(hg{1}(1), ...
                    'FaceColor', obj.CompartmentColour, ...
                    'FaceAlpha', 0.8                  , ...
                    'EdgeColor', 'k');
            end
            
            %Draw the payload
            if isempty(obj.PayloadCrossSections) || ~drawPL
                hg{2} = [];
            else 
                hg{2} = drawElement(obj.PayloadCrossSections, ht);
            end
            if ~isempty(hg{2})
                set(hg{2}, 'Tag', 'Payload');
                set(hg{2}(1), ...
                    'FaceColor', obj.PayloadColour, ...
                    'EdgeColor', 'k');
            end
            
            %Collapse hgobjects
            hg = horzcat(hg{:});
            
            if ~isempty(hg)
                hgS = hg(arrayfun(@(o) isa(o, ...
                    'matlab.graphics.chart.primitive.Surface'), hg));
                %Apply properties
                arrayfun(@(i) set(hgS, varargin{i}, varargin{i + 1}), 1 : 2 : numel(varargin));
            end
            
            %Draw the cabin - TODO : Update the drawing methods of a cabin
            if isa(obj, 'awi.model.Cabin')
                hg(end + 1 : end + 2) = drawCabinLayout(obj, ht);
            end
            
        end        
        function hg = drawMesh(obj, ht, varargin)
            %drawMesh Plots the triangular mesh using 'trisurf'.
            
            if numel(obj) > 1
                hg = arrayfun(@(o) drawMesh(o, ht, varargin{:}), obj, 'Unif', false);
                hg = horzcat(hg{:});
                return
            end
            
            %Get the mesh
            updateGeometry(obj);            
            TR = obj.CompartmentMesh;
            
            %Plot it
            hg = arrayfun(@(tr) trisurf(tr.faces, ...
                tr.vertices(:, 1), tr.vertices(:, 2), tr.vertices(:, 3), ...
                'Parent', ht), TR);
            
            %Update properties
            set(hg, varargin(1 : 2 : end), varargin(2 : 2 : end));
            
        end     
        function hg = drawMassLocations(obj, ht, varargin)
            %drawMassLocations Annotates the locations of the mass CoG for
            %each section in the compartment and for the total compartment.
            
            if numel(obj) > 1
               hg = arrayfun(@(o) drawMassLocations(o, ht, varargin{:}), obj, 'Unif', false); 
               hg = horzcat(hg{:});
               return
            end
            
            if isa(ht, 'matlab.graphics.axis.Axes')
                ht = hgtransform('Parent', ht);
                M  = makehgtform('translate', obj.BeamHandle.AbsPosition);
                set(ht, 'Matrix', M);
            end
            
            MassProps = obj.PayloadRigidBodyMassPropsBySection;
            
            %Draw the compartments first
            hg = drawElement(obj, ht, 'DrawPayload', true, ...
                {'FaceColor'}, {'red' ; 'blue'}, ...
                'FaceAlpha'  , 0.2, ...
                'EdgeColor'  , 'none');
            
            if isempty(MassProps) %Escape route
                return
            end
            
            %Draw markers at the CoG locations
            cg = vertcat(MassProps.CoG);
            hg(end + 1) = plot3(cg(:, 1), cg(:, 2), cg(:, 3), 'Parent', ht, ...
                'Marker'         , 'o', ...
                'MarkerFaceColor', 'm', ...
                'MarkerEdgeColor', 'k', ...
                'LineStyle'      , 'none', ...
                'MarkerSize'     , 10);
            hg(end + 1) = plot3( ...
                obj.PayloadCoG(1), obj.PayloadCoG(2), obj.PayloadCoG(3), ...
                'Parent', ht, ...
                'Marker'         , 's', ...
                'MarkerFaceColor', 'g', ...
                'MarkerEdgeColor', 'k', ...
                'LineStyle'      , 'none', ...
                'MarkerSize'     , 12);
        end
    end
    
    methods (Sealed) % public methods for updating the geometry/inertia properties
        function updateGeometry(obj)
            %updateGeometry Updates the geomertry properties of the
            %compartment. 
            %
            % Properties updated include:
            %   - Coordinates (compartment & payload)
            %   - Triangular mesh (compartment & payload)
            %   - Volume (compartment & payload)
                        
            if isempty(obj)
                return
            end
            if numel(obj) > 1 %Loop through object arrays... 
                arrayfun(@(o) updateGeometry(o), obj);
                return
            end
            
            %Update the cross-sections of the compartment
            generateCrossSections(obj)
            
            %Grab the global coordinates so we can create mesh
            [xC , yC , zC ] = calculateGlobalCoords(obj.CompartmentCrossSections);
            [xPL, yPL, zPL] = calculateGlobalCoords(obj.PayloadCrossSections);
            
            %Translate into the actual global frame - not the beam frame
            gPos = obj.BeamHandle.AbsPosition;
            xC   = xC + gPos(1);
            yC   = yC + gPos(2);
            zC   = zC + gPos(3);
            xPL  = xPL + gPos(1);
            yPL  = yPL + gPos(2);
            zPL  = zPL + gPos(3);         

            %Shrink factor for the 'boundary' function (0 = tight fit, 1 = convex hull)
            sf = 0;  
            
            %Split the compartment into bays and generate a triangular mesh             
            TR    = i_generateTriangularMesh(xC , yC , zC , sf);
            TR_PL = i_generateTriangularMesh(xPL, yPL, zPL, sf);
            obj.CompartmentMesh = TR;
            obj.PayloadMesh     = TR_PL;
            
            %Update volume
            obj.Volume        = sum([TR.volume]);
            obj.PayloadVolume = sum([TR_PL.volume]);
                       
            function TR = i_generateTriangularMesh(x, y, z, sf)
                %i_generateTriangularMesh Creates a triangular mesh for
                %each bay in the compartment coordinates.
                %                
                % TODO - Vectorise this operation! There is a problem with
                % generating a triangulation of the entire compartment
                % which is why it has been split into bays. 
  
                %Split into bays
                nBay = size(x, 1) - 1;
                
                %Preallocate
                TR  = repmat(struct('faces', [], 'vertices', [], 'volume', []), [1, nBay]);
                
                if any(any(isnan(x))) || any(any(isnan(y))) || any(any(isnan(z)))
                    return
                end
                
                for iB = 1 : nBay
                    
                    %Grab data
                    x_ = x([iB, iB + 1], :);
                    y_ = y([iB, iB + 1], :);
                    z_ = z([iB, iB + 1], :);
                    
                    %Remove duplicate points
                    data = [x_(:), y_(:), z_(:)];
                    data = unique(data, 'rows');
                    
                    %Calculate boundary triangulation
                    [k, v] = boundary(data, sf);
                    
                    %Assign to structure
                    TR(iB).faces    = k;
                    TR(iB).vertices = data;
                    TR(iB).volume   = v;                    
                    
                end
                
            end
            
        end        
        function updateMassAndInertia(obj)
            %updateMassAndInertia Updates the mass properties of the
            %compartment. 
            %
            % Properties updated include:
            %   - Mass
            %   - CoG
            %   - InertiaMatrix
            %   - BayInertialProperties
            
            if numel(obj) > 1 %Loop through object arrays... 
                arrayfun(@(o) updateMassAndInertia(o), obj);
                return
            end
            
            %Update the geometry first
            updateGeometry(obj);
            
            %Sometimes the 'boundary; triangulation does not use all the
            %points - This is not a problem so suppress the warning.
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId')
            
            %Calculate the rigid body properties
            calculateRigidBodyProperties(obj);
            
            %Calculate payload mass using rigid body parameters
            calculatePayloadMassProperties(obj)
            
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId')
            
        end
        function pf = calculatePayloadFraction(obj, mass)
            %calculatePayloadFraction Calculates the 'PayloadFraction'
            %required to fit a payload with mass 'mass' into the
            %compartment.
                       
            %Parse
            validateattributes(mass, {'numeric'}, {'scalar', 'finite', ...
                'nonnegative', 'real'}, class(obj), 'mass');
            if obj.MassDensity == 0 || mass == 0
                pf = 0;
                return
            end
            
            %What fraction of the total possible mass is being requested?
            compMass     = obj.Volume * obj.MassDensity;
            compMassFrac = mass / compMass;            
            
            %Cannot allocate more mass than the compartment has volume for
            if compMassFrac > 1
                pf = nan;
                return
            end
            
            %Retain the current payload fraction
            pf0 = obj.PayloadFraction;
            
            %Use a line-search method to calculate the payload fraction
            %which yields the desired mass
            pf  = compMassFrac;
            res = i_Jx(obj, pf, mass);
            function res = i_Jx(obj, pf, mass)
                %i_Jx Calculates the mass residual for a given payload
                %fraction.
                
                obj.PayloadFraction = pf;
                updateMassAndInertia(obj);
                res = mass - obj.PayloadMass;
                
            end
            %Bounds
            xlb = 0;
            xub = 1;
            %Golden Ratio
            phi = 0.5 * (1 + sqrt(5));  
            %Search criteria
            tol  = 10;
            kMax = 20;
            %Do the line-search (Golden Section)
            k = 1;
            while abs(res) > tol
                if k > kMax
                    break
                end
                %Search direction?
                if sign(res) < 1
                    xub  = pf;
                    delx = (xub - xlb) / phi;
                    pf   = xub - delx;
                else
                    xlb  = pf;
                    delx = (xub - xlb) / phi;
                    pf   = xlb + delx;
                end
                %Evalauate mass
                res = i_Jx(obj, pf, mass);
                %Counter
                k = k + 1;
            end
            
            %Return the final payload fraction & reset compartment
            pf = obj.PayloadFraction;
            obj.PayloadFraction = pf0;
            updateMassAndInertia(obj);
            
        end
        function generateCrossSections(obj)
            %generateCrossSections Creates the 'awi.model.CrossSection'
            %objects that represent this compartment. The cross-sections
            %are calculated by interpolating the existing cross-sections
            %that have been assigned to the parent 'awi.model.Beam' object.
            %
            % The 'Compartment' object simply defines the (x,z) extent of
            % the compartment, therefore, we still need to check for
            % existing cross-sections that occupy the 'eta' range of
            % interest and then retain the profiles at these points.
            
            %Can't define a compartment without a cross-section
            if isempty(obj.BeamHandle)
                warning(['An ''awi.model.Compartment'' object must '       , ...
                    'be assigned to a Beam object in order to define a '   , ...
                    'compartment. Use the method ''addCompartment'' '      , ...
                    'of the Beam class to assign a Compartment to a Beam.'], obj.BeamHandle.Name);
                return
            end
            if isempty(obj.BeamHandle.CrossSection)
                warning(['The ''awi.model.Beam'' object "%s" must ', ...
                    'have a cross-section in order to define a ' , ...
                    'compartment.'], obj.BeamHandle.Name);
                return
            end
            
            %Grab local(x,z) coordinates of the cross-section
            [CS, etaCS] = i_getCrossSectionCoords(obj);
            etaC = obj.EtaLocations;
            
            function [CS, etaCS] = i_getCrossSectionCoords(obj)
                %i_getCrossSectionCoords Retrieves the coordinates of the
                %cross-section objects which lie within, or immediately
                %either side of, the comparment 'eta' range.
                
                CS_    = obj.BeamHandle.CrossSection;
                etaCS_ = obj.BeamHandle.CrossSection_eta;
                
                %Only retain data that is in the 'eta' range of interest
                idx_ = and(etaCS_ >= min(obj.EtaLocations), etaCS_ <= max(obj.EtaLocations));
                
                %Grab additional cross-section so no need to extrapolate
                if ~all(idx_)
                    ind_ = find(idx_, 1, 'first') - 1;
                    if ind_ > 0
                        idx_(ind_) = true;
                    end
                    ind_ = find(idx_, 1, 'last') + 1;
                    if ind_ <= numel(idx_)
                        idx_(ind_) = true;
                    end
                end
                
                %Check if we actually have enough data to proceed!
                if nnz(idx_) < 2
                    %Get the cross-sections immediately bordering the eta
                    %range of the compartment
                    ind(1) = find(etaCS_ <= min(obj.EtaLocations), 1, 'last');
                    ind(2) = find(etaCS_ >= max(obj.EtaLocations), 1, 'first');
                    idx_ = false(size(idx_));
                    idx_(ind) = true;
                end
                
                %Get coordinates of the cross-section
                etaCS = etaCS_(idx_);
                CS    = CS_(idx_);
                
            end
            
            %Normalise each cross-section
            [xNormCS, zNormCS, xzDataCS] = normaliseXZData(CS);
            
            %Interpolate the CrossSection (x,z) data to new eta positions
            %of the Compartment
            [xUpper, xLower, zUpper, zLower] = interpolateNormXZData(CS, etaC);
            xNormC = [xUpper, fliplr(xLower)]; %Compartment data
            zNormC = [zUpper, fliplr(zLower)];
            
            %Get data for the start/end of the compartment
            [x1, xN, z1, zN] = getOffsetData(obj);
            compDataC = [x1, xN, z1, zN];
            
            %Return data at the eta positions of interest
            [xNorm, zNorm, eta, xzData, compData] = i_indexCrossSection( ...
                xNormCS, zNormCS, xNormC, zNormC, etaCS, etaC, xzDataCS, compDataC);
            
            %Crop the cross-section coordinates to get the coordinates of
            %the compartment
            nX = obj.NumPoints;
            [xCoords, zCoords] = i_cropCoordinates(xNorm, zNorm, compData, nX);
            
            %Calculate coordinates of payload using payload fraction
            pf = obj.PayloadFraction;
            if pf == 0
                %No payload? Just use NaN so we don't plot anything
                xCoordsPL = nan(size(xCoords));
                zCoordsPL = xCoordsPL;
            else
                %Caclulate new upper/lower limits for the z-coordinates
                compData(:, 4) = ((compData(:, 4) - compData(:, 3)) .* pf) + compData(:, 3);
                [xCoordsPL, zCoordsPL] = i_cropCoordinates(xNorm, zNorm, compData, nX);
            end
            
            %Scale the compartment coordinates & add offset
            [xCoords, zCoords]     = scaleProfileData(xCoords  , zCoords  , xzData);
            [xCoordsPL, zCoordsPL] = scaleProfileData(xCoordsPL, zCoordsPL, xzData);
            
            %Generate 'awi.model.CrossSection' objects
            nCS       = numel(eta);
            CompCS    = arrayfun(@(~) awi.model.CrossSection, 1 : nCS);
            PayloadCS = arrayfun(@(~) awi.model.CrossSection, 1 : nCS);
            
            %Interpolate the cross-section orientation (rotation matrix)
            CoordSystems = [obj.BeamHandle.CrossSection.Orientation];
            etaCSys      = [obj.BeamHandle.CrossSection.BeamEta]';
            RMat  = get(CoordSystems, 'RMatrix_');
            RMat  = horzcat(RMat{:})';
            RMat_ = interp1(etaCSys, RMat, eta, 'linear');
            RMat  = arrayfun(@(i) reshape(RMat_(i, :)', [3, 3]), 1 : nCS, 'Unif', false);
            
            %Assign data
            set([CompCS  , PayloadCS], 'BeamHandle', obj.BeamHandle);
            set(CompCS   , {'BeamEta'} , num2cell(eta)');
            set(PayloadCS, {'BeamEta'} , num2cell(eta)');
            set(CompCS   , {'X'}       , num2cell(xCoords', 2));
            set(CompCS   , {'Z'}       , num2cell(zCoords', 2));
            set(PayloadCS, {'X'}       , num2cell(xCoordsPL', 2));
            set(PayloadCS, {'Z'}       , num2cell(zCoordsPL', 2));
            set([PayloadCS.Orientation], {'RMatrix'}, RMat');
            set([CompCS.Orientation]   , {'RMatrix'}, RMat');
            
            %Assign to object
            obj.CompartmentCrossSections = CompCS;
            obj.PayloadCrossSections     = PayloadCS;
            
            function [xNorm, zNorm, eta, xzData, offsetData] = i_indexCrossSection(xNormCS, zNormCS, xNormC, zNormC, etaCS, etaC, xzDataCS, compDataC)
                %i_indexCrossSections Combines the coordinates & data that
                %is defined at the cross-section 'eta' positions & the
                %compartment 'eta' positions and returns the data at only
                %the eta positions that sit within the bounds defined by
                %the 'EtaLocations' property of the Compartment.
                
                %Interpolate the scale/translation data to compartment eta
                xzDataC = interp1(etaCS, xzDataCS, etaC);
                
                %Interpolate the comaprtment start/end data to the
                %cross-section eta
                offsetDataCS = interp1(etaC, compDataC, etaCS', 'linear', 'extrap');
                
                %Gather coordinates
                eta        = [etaCS, etaC'];
                offsetData = [offsetDataCS ; compDataC];
                xzData     = [xzDataCS ; xzDataC];
                xNorm      = [xNormCS ; xNormC];
                zNorm      = [zNormCS ; zNormC];
                
                %Sort by 'eta' positions and retain unique values only
                [eta, ind] = unique(eta);
                xzData     = xzData(ind, :);
                offsetData = offsetData(ind, :);
                xNorm      = xNorm(ind, :);
                zNorm      = zNorm(ind, :);
                
                %Remove data outside compartment 'eta' range now that
                %interpolation is complete
                idx = and(eta >= min(etaC), eta <= max(etaC));
                eta = eta(idx);
                xzData     = xzData(idx, :);
                offsetData = offsetData(idx, :);
                xNorm   = xNorm(idx, :);    %Note - 'xNorm' & 'zNorm' are in
                zNorm   = zNorm(idx, :);    %plotting format (i.e. one continuous loop)
                
            end
            
            function [xCoords, zCoords] = i_cropCoordinates(xNorm, zNorm, compData, nX)
                %i_cropCoordinates Removes any coordinate data outside the
                %(x,z) range specified by x1 -> xN and z1 -> zN
                
                x1 = compData(:, 1);
                xN = compData(:, 2);
                z1 = compData(:, 3);
                zN = compData(:, 4);
                
                %How many points in the current distribution?
                nX_ = size(xNorm, 2) / 2;
                
                %Preallocate
                nProfile = size(xNorm, 1);
                xCoords  = zeros((2 * nX) + 1, nProfile);
                zCoords  = xCoords;
                
                %Crop the cross-section at the desired locations
                for iC = 1 : nProfile
                    %Split into upper and lower
                    xU = xNorm(iC, 1 : nX_);
                    xL = xNorm(iC, nX_ + 1 : end);
                    zU = zNorm(iC, 1 : nX_);
                    zL = zNorm(iC, nX_ + 1 : end);
                    %Interpolate the z values at x1/xN
                    zU1N = interp1(xU, zU, [x1(iC), xN(iC)]);
                    zL1N = interp1(xL, zL, [x1(iC), xN(iC)]);
                    %Crop the x-coordinates
                    idxU = and(xU >= x1(iC), xU<= xN(iC));
                    idxL = and(xL >= x1(iC), xL<= xN(iC));
                    xU = [x1(iC), xU(idxU), xN(iC)];
                    xL = [xN(iC), xL(idxL), x1(iC)];
                    zU = [zU1N(1), zU(idxU), zU1N(2)];
                    zL = [zL1N(2), zL(idxL), zL1N(1)];
                    %Crop the z-coordinates
                    idxU = and(zU >= z1(iC), zU <= zN(iC));
                    idxL = and(zL >= z1(iC), zL <= zN(iC));
                    zU   = [zN(iC), zU(idxU), zN(iC)];
                    zL   = [z1(iC), zL(idxL), z1(iC)];
                    xU   = [x1(iC), xU(idxU), xN(iC)];
                    xL   = [xN(iC), xL(idxL), x1(iC)];
                    %Remove points outside cross-section
                    [inU, onU] = inpolygon(xU, zU, xNorm(iC, :), zNorm(iC, :));
                    [inL, onL] = inpolygon(xL, zL, xNorm(iC, :), zNorm(iC, :));
                    idxU = or(inU, onU);
                    idxL = or(inL, onL);
                    xU = xU(idxU);
                    zU = zU(idxU);
                    xL = xL(idxL);
                    zL = zL(idxL);
                    %Remove duplicate data
                    %                     [xU, ind] = unique(xU, 'stable');
                    [xU, ind] = unique(xU);
                    zU = zU(ind);
                    %                     [xL, ind] = unique(xL, 'stable');
                    [xL, ind] = unique(xL);
                    zL = zL(ind);
                    %Resample the x-locations and interpolate the z-coordinates
                    xEnd = [];
                    zEnd = [];
                    if isempty(xU)
                        xU_ = linspace(xL(1), xL(end), nX);
                        zU_ = repmat(zL(1), [1, nX]);
                    else
                        xU_  = linspace(min(xU), max(xU), nX);
                        zU_  = interp1(xU, zU, xU_);
                        xEnd = xU_(1);
                        zEnd = zU_(1);
                    end
                    if isempty(xL)
                        xL_ = linspace(xU(1), xU(end), nX);
                        zL_ = repmat(zU(1), [1, nX]);
                    else
                        xL_ = linspace(min(xL), max(xL), nX);
                        zL_ = interp1(xL, zL, xL_);
                    end
                    if isempty(xEnd)
                        xEnd = xL_(1);
                        zEnd = zL_(1);
                    end
                    %Assign to matrix
                    xCoords(:, iC) = [xU_ , fliplr(xL_), xEnd]';
                    zCoords(:, iC) = [zU_ , fliplr(zL_), zEnd]';
                end
                
            end
            
        end
    end
        
    methods (Sealed, Access = private) %Calculating the compartment inertia properties
        function calculateRigidBodyProperties(obj)
            %calculateRigidBodyProperties Caluclates the rigid body
            %properties such as volume, inertia tensor, 3D moments, etc.
            %using the File Exchange Function 'RigidBodyParams'.
            %
            % Note, these properties are irrespective of mass.
            
            if numel(obj) > 1 %Loop through object arrays... 
                arrayfun(@(o) calculateRigidBodyProperties(o), obj);
                return
            end
            
            TR   = obj.CompartmentMesh;
            TRPL = obj.PayloadMesh;
            nBay = numel(TR);
            
            if isempty(TR) || isempty(TRPL) %Escape route
                return
            end
            
            %Preallocate
            fNames = {'volume', 'centroid', 'inertia_tensor', 'PAI', 'eigs', 'moments'};
            RBP    = repmat(cell2struct(cell(size(fNames)), fNames, 2), [1, nBay]);
            RBP_PL = RBP;
            
            %Calculate rigid-body properties 
            for iB = 1 : nBay
                [RBP(iB), ~] = RigidBodyParams(TR(iB));
            end
            if obj.PayloadFraction ~= 0
                for iB = 1 : nBay
                    [RBP_PL(iB), ~] = RigidBodyParams(TRPL(iB));
                end
            end
            
            %Assign to the object
            obj.CompartmentRigidBodyPropsBySection = RBP;
            obj.PayloadRigidBodyPropsBySection     = RBP_PL;
            
        end
        function calculatePayloadMassProperties(obj)
            %calculatePayloadMassProperties Calculates the mass properties
            %of the payload, split by section and for the compartment. If
            %the compartment payload fraction ('PayloadFraction') is equal
            %to zero then then no further mass data is generated after the
            %mass has been set to zero.
            %
            % Mass properties:
            %   - Mass (section & compartment)
            %   - Centre of gravity (section & compartment)
            %   - Inertia Tensor (section & compartment)
            
            if numel(obj) > 1 %Loop through object arrays... 
                arrayfun(@(o) calculatePayloadMassProperties(o), obj);
                return
            end
            
            RBP_PL   = obj.PayloadRigidBodyPropsBySection;            
            massdata = [obj.PayloadFraction, obj.MassDensity];
            
            if any(massdata == 0) || isempty(RBP_PL) %Escape route
                obj.PayloadMass = 0;
                return
            end
            
            %% Calculate mass properties of individual sections
                
            %Preallocate
            nBay      = numel(RBP_PL);
            MassProps = repmat(struct('Mass', [], 'CoG', [], ...
                'InertiaTensor', []), [1, nBay]);
            
            %Multiply by the mass density to get mass & inertia properties 
            mass = cellfun(@(x) obj.MassDensity .* x, {RBP_PL.volume}, 'Unif', false);
            cg   = {RBP_PL.centroid};
            Iij  = cellfun(@(x) obj.MassDensity .* x, {RBP_PL.inertia_tensor}, 'Unif', false);
           
            %Assign to structure array
            [MassProps.Mass]          = deal(mass{:});
            [MassProps.CoG]           = deal(cg{:});
            [MassProps.InertiaTensor] = deal(Iij{:});
            
            %Assign to the object
            obj.PayloadRigidBodyMassPropsBySection = MassProps;
            
            %% Calculate the mass properties of the entire compartment
            
            %Collect
            mass = vertcat(mass{:});
            cg   = vertcat(cg{:});
            
            %Calculate the CoG of the payload masses in each section
            M   = sum(mass);
            cgm = cg .* mass;
            CoG  = sum(cgm, 1) ./ M;
            
            %Calculate contributions to the inertia tensor from the
            %parallel axis theorem
            dxyz  = CoG - cg;
            dxyz2 = dxyz.^2;
            
            dIxx = sum(mass .* (dxyz2(:, 2) + dxyz(:, 3)));
            dIyy = sum(mass .* (dxyz2(:, 1) + dxyz(:, 3)));
            dIzz = sum(mass .* (dxyz2(:, 1) + dxyz(:, 2)));
            dIxy = sum(mass .* dxyz(:, 1) .* dxyz(:, 2));
            dIxz = sum(mass .* dxyz(:, 1) .* dxyz(:, 3));
            dIyz = sum(mass .* dxyz(:, 2) .* dxyz(:, 3));
            
            %Total inertia tensor (no parallel axis)
            Iij  = sum(cat(3, Iij{:}), 3);
            
            %Contribution form parallel axis theorem
            % - TODO - Check if we need negative signs here...
            dIij = [ ...
                dIxx, dIxy, dIxz ; ...
                dIxy, dIyy, dIyz ; ...
                dIxz, dIyz, dIzz ];
            
            %Total inertia tensor
            Iij = Iij + dIij;
            
            %Assign to the object
            obj.PayloadMass    = M;
            obj.PayloadCoG     = CoG;
            obj.PayloadInertia = Iij;            
            
        end
    end
    
    methods %Helper functions
        function [x1, xN, z1, zN] = getOffsetData(obj)
            %getOffsetData Returns the values of the start/end points in
            %the local (x,z) coordinate system with the correct number of
            %elements in each variable, according to the number of points
            %defined in 'obj.EtaLocations'.
            %
            % Each variable is padded using the final value.

            n = numel(obj.EtaLocations);
            
            x1 = obj.LocalXStart;
            xN = obj.LocalXEnd;
            z1 = obj.LocalZStart;
            zN = obj.LocalZEnd;
            
            x1 = [x1 ; repmat(x1(end), [n - numel(x1), 1])];
            xN = [xN ; repmat(xN(end), [n - numel(xN), 1])];
            z1 = [z1 ; repmat(z1(end), [n - numel(z1), 1])];
            zN = [zN ; repmat(zN(end), [n - numel(zN), 1])];
        
        end
    end
    
end

%% old code

         
%         function [Coords, PayloadCoords] = getCoords(obj)
%             %getCoords Calculates the coordinates of the compartment in the
%             %global coordinate system.
%             %
%             % TODO - Interpolate the cross-section "Orientation" property
%             % as well.
%             %
%             % The 'Compartment' object simply defines the (x,z) extent of
%             % the compartment, therefore, we still need to check for
%             % existing cross-sections that occupy the 'eta' range of
%             % interest and then retain the profiles at these points.
%             
%             Coords        = [];
%             PayloadCoords = [];
%             
%             %Can't define a compartment without a cross-section
%             if isempty(obj.BeamHandle.CrossSection)
%                 warning(['The ''awi.model.Beam'' object "%s" must ', ...
%                     'have a cross-section in order to define a ' , ...
%                     'compartment.'], obj.BeamHandle.Name);
%                 return
%             end
%                                     
%             %Grab local(x,z) coordinates of the cross-section
%             [x, z, etaCS] = i_getCrossSectionCoords(obj);
%             etaC = obj.EtaLocations;
%             
%             %Normalise each cross-section
%             [xNormCS, zNormCS, xzDataCS] = i_normaliseCrossSection(x, z);
%             
%             %Assume the cross-section have equal number of coordinates on
%             %upper and lower curves
%             nX = size(xNormCS, 2) / 2;
%             
%             %Interpolate the aerofoil coordinates, scale factors & (x,z)
%             %origin
%             [xUpper, xLower, zUpper, zLower, ~] = ...
%                 awi.model.LiftingSurface.interpolateAerofoilProfile( ...
%                 xNormCS', zNormCS', etaCS, nX, etaC);   
%             xNormC = [xUpper, fliplr(xLower)]; %Compartment data
%             zNormC = [zUpper, fliplr(zLower)];
%         
%             %Get data for the start/end of the compartment
%             [x1, xN, z1, zN] = getOffsetData(obj);
%             compDataC = [x1, xN, z1, zN];
%             
%             %Return data at the eta positions of interest
%            [xNorm, zNorm, eta, xzData, compData] = indexCrossSection( ...
%                xNormCS, zNormCS, xNormC, zNormC, etaCS, etaC, xzDataCS, compDataC);
%                         
%             %Crop the cross-section coordinates to get the coordinates of
%             %the compartment
%             nX = obj.NumPoints;
%             [xCoords, zCoords] = i_cropCoordinates(xNorm, zNorm, compData, nX);
%             
%             %Calculate coordinates of payload using payload fraction
%             pf = obj.PayloadFraction;
%             if pf == 0
%                 %No payload? Just use NaN so we don't plot anything
%                 xCoordsPL = nan(size(xCoords));
%                 zCoordsPL = xCoordsPL;
%             else
%                 %Caclulate new upper/lower limits for the z-coordinates
%                 compData(:, 4) = ((compData(:, 4) - compData(:, 3)) .* pf) + compData(:, 3);
%                 [xCoordsPL, zCoordsPL] = i_cropCoordinates(xNorm, zNorm, compData, nX);
%             end
%                         
%             %Scale the compartment coordinates & add offset
%             [xCoords, zCoords]     = i_scaleCoords(xCoords  , zCoords  , xzData);
%             [xCoordsPL, zCoordsPL] = i_scaleCoords(xCoordsPL, zCoordsPL, xzData);
%             
%             %Generate 'awi.model.CrossSection' objects
%             nCS       = numel(eta);
%             CompCS    = arrayfun(@(~) awi.model.CrossSection, 1 : nCS);
%             PayloadCS = arrayfun(@(~) awi.model.CrossSection, 1 : nCS);
%             
%             %Interpolate the cross-section orientation (rotation matrix)
%             CoordSystems = [obj.BeamHandle.CrossSection.Orientation];
%             etaCSys      =[obj.BeamHandle.CrossSection.BeamEta]';
%             RMat = get(CoordSystems, 'RMatrix_');
%             RMat = horzcat(RMat{:})';
%             RMat_ = interp1(etaCSys, RMat, eta, 'linear');
%             RMat  = arrayfun(@(i) reshape(RMat_(i, :)', [3, 3]), 1 : nCS, 'Unif', false);
%             
%             %Assign data
%             set([CompCS  , PayloadCS], 'BeamHandle', obj.BeamHandle);
%             set(CompCS   , {'BeamEta'} , num2cell(eta)');
%             set(PayloadCS, {'BeamEta'} , num2cell(eta)');
%             set(CompCS   , {'X'}       , num2cell(xCoords', 2));
%             set(CompCS   , {'Z'}       , num2cell(zCoords', 2));
%             set(PayloadCS, {'X'}       , num2cell(xCoordsPL', 2));
%             set(PayloadCS, {'Z'}       , num2cell(zCoordsPL', 2));
%             set([PayloadCS.Orientation], {'RMatrix'}, RMat');
%             set([CompCS.Orientation]   , {'RMatrix'}, RMat');
%             
%             %Assign to object
%             obj.CompartmentCrossSections = CompCS;
%             obj.PayloadCrossSections     = PayloadCS;
%      
%             function [x, z, etaCS] = i_getCrossSectionCoords(obj)
%                 %i_getCrossSectionCoords Retrieves the coordinates of the
%                 %cross-section objects which lie within, or immediately
%                 %either side of, the comparment 'eta' range.
%                 
%                 CS    = obj.BeamHandle.CrossSection;
%                 etaCS = obj.BeamHandle.CrossSection_eta;
% 
%                 %Only retain data that is in the 'eta' range of interest
%                 idx_ = and(etaCS >= min(obj.EtaLocations), etaCS <= max(obj.EtaLocations));
%                 
%                 %Grab additional cross-section so no need to extrapolate
%                 if ~all(idx_)
%                     ind_ = find(idx_, 1, 'first') - 1;
%                     if ind_ > 0
%                         idx_(ind_) = true;
%                     end
%                     ind_ = find(idx_, 1, 'last') + 1;
%                     if ind_ <= numel(idx_)
%                         idx_(ind_) = true;
%                     end
%                 end
%                 etaCS = etaCS(idx_);
%                 
%                 %Get coordinates of the cross-section
%                 x = vertcat(CS(idx_).X);
%                 z = vertcat(CS(idx_).Z);
%                 
%             end
%             
%             function [xNorm, zNorm, xzData] = i_normaliseCrossSection(x, z)
%                 %i_normaliseCrossSection Nornalises the cross-section
%                 %coordinates between [0 : 1] and centres the cross-section
%                 %of the (0.5, 0.5) coordinate.
%                 %
%                 % Returns the normalised coordinates and a matrix of data,
%                 % each column represents the offset and scale data for the
%                 % (x,z) axes.
%                 
%                 xMin   = min(x, [], 2);
%                 xMax   = max(x, [], 2);
%                 zMin   = min(z, [], 2);
%                 zMax   = max(z, [], 2);
%                 
%                 xScale = (xMax - xMin);
%                 zScale = (zMax - zMin);
%                 
%                 xNorm  = (x - xMin) ./ xScale;
%                 zNorm  = (z - zMin) ./ zScale;
%                 
%                 xzData = [xMin, zMin, xScale, zScale];
%                 
%             end
%             
%             function [xNorm, zNorm, eta, xzData, offsetData] = indexCrossSection(xNormCS, zNormCS, xNormC, zNormC, etaCS, etaC, xzDataCS, compDataC)
%                 %indexCrossSections Combines the coordinates & data that is
%                 %defined at the cross-section 'eta' positions & the
%                 %compartment 'eta' positions and returns the data at only
%                 %the eta positions that sit within the bounds defined by
%                 %the 'EtaLocations' property of the Compartment.
%                 
%                 %Interpolate the scale/translation data to compartment eta
%                 xzDataC = interp1(etaCS, xzDataCS, etaC);
%                 
%                 %Interpolate the comaprtment start/end data to the
%                 %cross-section eta
%                 offsetDataCS = interp1(etaC, compDataC, etaCS', 'linear', 'extrap');
%                 
%                 %Gather coordinates
%                 eta        = [etaCS, etaC'];
%                 offsetData = [offsetDataCS ; compDataC];
%                 xzData     = [xzDataCS ; xzDataC];
%                 xNorm      = [xNormCS ; xNormC];
%                 zNorm      = [zNormCS ; zNormC];
%                 
%                 %Sort by 'eta' positions and retain unique values only
%                 [eta, ind] = unique(eta);
%                 xzData     = xzData(ind, :);
%                 offsetData = offsetData(ind, :);
%                 xNorm      = xNorm(ind, :);
%                 zNorm      = zNorm(ind, :);
%                 
%                 %Remove data outside compartment 'eta' range now that
%                 %interpolation is complete
%                 idx = and(eta >= min(etaC), eta <= max(etaC));
%                 eta = eta(idx);
%                 xzData     = xzData(idx, :);
%                 offsetData = offsetData(idx, :);
%                 xNorm   = xNorm(idx, :);    %Note - 'xNorm' & 'zNorm' are in
%                 zNorm   = zNorm(idx, :);    %plotting format (i.e. one continuous loop)
%                 
%             end
%             
%             function [xCoords, zCoords] = i_cropCoordinates(xNorm, zNorm, compData, nX)
%                 %i_cropCoordinates Removes any coordinate data outside the
%                 %(x,z) range specified by x1 -> xN and z1 -> zN
% 
%                 x1 = compData(:, 1);
%                 xN = compData(:, 2);
%                 z1 = compData(:, 3);
%                 zN = compData(:, 4);
%                 
%                 %Preallocate
%                 nProfile = size(xNorm, 1);
%                 xCoords  = zeros((2 * nX) + 1, nProfile);
%                 zCoords  = xCoords;
%                 
%                 %Crop the cross-section at the desired locations
%                 for iC =1 : nProfile
%                     %Split into upper and lower
%                     xU = xNorm(iC, 1 : nX);
%                     xL = xNorm(iC, nX + 1 : end);
%                     zU = zNorm(iC, 1 : nX);
%                     zL = zNorm(iC, nX + 1 : end);
%                     %Interpolate the z values at x1/xN
%                     zU1N = interp1(xU, zU, [x1(iC), xN(iC)]);
%                     zL1N = interp1(xL, zL, [x1(iC), xN(iC)]);
%                     %Crop the x-coordinates
%                     idxU = and(xU >= x1(iC), xU<= xN(iC));
%                     idxL = and(xL >= x1(iC), xL<= xN(iC));
%                     xU = [x1(iC), xU(idxU), xN(iC)];
%                     xL = [xN(iC), xL(idxL), x1(iC)];
%                     zU = [zU1N(1), zU(idxU), zU1N(2)];
%                     zL = [zL1N(2), zL(idxL), zL1N(1)];
%                     %Crop the z-coordinates
%                     idxU = and(zU >= z1(iC), zU <= zN(iC));
%                     idxL = and(zL >= z1(iC), zL <= zN(iC));
%                     zU   = [zN(iC), zU(idxU), zN(iC)];
%                     zL   = [z1(iC), zL(idxL), z1(iC)];
%                     xU   = [x1(iC), xU(idxU), xN(iC)];
%                     xL   = [xN(iC), xL(idxL), x1(iC)];
%                     %Remove points outside cross-section
%                     [inU, onU] = inpolygon(xU, zU, xNorm(iC, :), zNorm(iC, :));
%                     [inL, onL] = inpolygon(xL, zL, xNorm(iC, :), zNorm(iC, :));
%                     idxU = or(inU, onU);
%                     idxL = or(inL, onL);
%                     xU = xU(idxU);
%                     zU = zU(idxU);
%                     xL = xL(idxL);
%                     zL = zL(idxL);
%                     %Remove duplicate data
%                     [xU, ind] = unique(xU, 'stable');
%                     zU = zU(ind);
%                     [xL, ind] = unique(xL, 'stable');
%                     zL = zL(ind);
%                     %Resample the x-locations and interpolate the z-coordinates
%                     xU_ = linspace(xU(1), xU(end), nX);
%                     xL_ = linspace(xL(1), xL(end), nX);
%                     zU_ = interp1(xU, zU, xU_);
%                     zL_ = interp1(xL, zL, xL_);
%                     %Assign to matrix
%                     xCoords(:, iC) = [xU_ , xL_, xU_(1)]';
%                     zCoords(:, iC) = [zU_ , zL_, zU_(1)]';
%                 end
%                 
%             end
%             
%             function [xCoords, zCoords] = i_scaleCoords(xCoords, zCoords, xzData)
%                 %i_scaleCoords Returns the coordinates to their original
%                 %dimensions using the scale/translate data in 'xzData'.
%                 
%                 n   = size(xCoords, 1);
%                 xSc = repmat(xzData(:, 3)', [n, 1]);
%                 zSc = repmat(xzData(:, 4)', [n, 1]);
%                 x0  = xzData(:, 1)';
%                 z0  = xzData(:, 2)';
%                 xCoords = xCoords .* xSc + x0;
%                 zCoords = zCoords .* zSc + z0;
%                 
%             end
% 
%         end
%         function [GlobalCoords_Comp, GlobalCoords_PL] = getGlobalCoords(obj)
%             %getGlobalCoords Returns the coordinates in the global
%             %coordinate system.
%             %
%             % TODO - This needs updating so we can account for beams that
%             % do not have their span-vector in the y-axes.
%                         
%             %Get coordinates of the compartments
%             [CompCoords, PLCoords] = getCoords(obj);            
%             
%             %TODO - Fix the issue with unique in mvc class
%             %Check we only have one parent beam
%             pb = [obj.BeamHandle];
% %             assert(numel(unique(pb)) == 1, ['Cannot handle the case ', ...
% %                 'where we have compartments from multiple beams. ' , ...
% %                 'Update code.']);
% 
%             %For objects of class 'BeamPropertyObject' the 'eta_flag' is
%             %always set to 'R'. i.e. Along the axis of the line -> Get the
%             %length of the line then calculate the (x,y,z) coordinates of
%             %the object along the beam
%             [xd, yd, zd] = xyzdata(pb(1));
%             r = awi.model.Stick.getLineLength(xd, yd, zd);
%             eta  = r ./ r(end);
%             etaC = CompCoords.Eta;
%             
%             %Can do it in one call to 'interp1'
%             xyzC = interp1(eta, [xd ; yd ; zd]', etaC);
%             
%             %Gather coordinates and rotate through the 'Orientation'
%             %rotation matrix
%             
%             y = zeros(size(CompCoords.X));
%            
%             %Add the offsets along the beam
%             GlobalCoords_Comp.X = CompCoords.X + xyzC(:, 1)';
%             GlobalCoords_Comp.Y = y + xyzC(:, 2)';
%             GlobalCoords_Comp.Z = CompCoords.Z + xyzC(:, 3)';   
%             %
%             GlobalCoords_PL.X = PLCoords.X + xyzC(:, 1)';
%             GlobalCoords_PL.Y = y + xyzC(:, 2)';
%             GlobalCoords_PL.Z = PLCoords.Z + xyzC(:, 3)';
%             
%         end     


% methods (Static) %Analytic functions
%         function volume = calculateVolume(x, y, z, method)
%             %calculateVolume Calculates the volume of the shape given by
%             %coordinates (x, y, z) using the method specified in 'method'.
%             
%             validmethods = {'deluanay'};
%             
%             if nargin < 4 || isempty(method)
%                 method = 'delaunay';
%             end
%             
%             switch method 
%                 case 'delaunay'
%                     %Split the shape into a set of triangles
%                     
%                     hF = figure; 
%                     hAx = axes('Parent', hF, 'NextPlot', 'add');
%                     hg = drawElement(Comp, hAx);
%                     set(hg, 'EdgeColor', 'none');
%                     axis equal 
%                     view([0, 90]);
%                     
%                     %Simple
%                     figure; grid on, hold on, box on
%                     tri1 = delaunay(x, y, z);                    
%                     trisurf(tri1, x, y, z);
%                     title('delaunay -> trisurf');
%                     axis equal
%                     view([0, 90]);
%                     
%                     %Remove duplicate points first
%                     figure; grid on, hold on, box on
%                     data = unique([x, y, z], 'rows');
%                     tri2 = delaunay(data);
%                     trisurf(tri2, data(:, 1), data(:, 2), data(:, 3));
%                     title('delaunay -> trisurf (unique points)');
%                     axis equal
%                     view([0, 90]);
%                     
% %                     %Trianguluation
% %                     figure
% %                     tri3 = delaunay(x, y, z);
% %                     TR = triangulation(tri3, x, y, z);
% % %                     [K, vol1] = convexHull(TR);
% %                     triplot(TR);
% %                     title('Convex Hull - TRIANGULATION');
%                     
%                     %Convex Hull
%                     figure; grid on, hold on, box on
%                     DT = delaunayTriangulation(x, y, z);
%                     [K, vol] = convexHull(DT);
%                     trisurf(K, DT.Points(:, 1), DT.Points(:, 2), DT.Points(:, 3));
%                     title('delaunayTriangulation -> convexHull -> trisurf');
%                     axis equal
%                     view([0, 90]);
%                     
%                     %Convex Hull
%                     figure; grid on, hold on, box on
%                     DT = delaunayTriangulation(data(:, 1), data(:, 2), data(:, 3));
%                     [K, vol] = convexHull(DT);
%                     trisurf(K, DT.Points(:, 1), DT.Points(:, 2), DT.Points(:, 3));
%                     title('delaunayTriangulation -> convexHull -> trisurf (unique points)');
%                     axis equal
%                     view([0, 90]);
%                     
%                     %Boundary
%                     figure; grid on, hold on, box on
%                     [k, v] = boundary(x, y, z, 0);
%                     trisurf(k, x, y, z);
%                     title('boundary (shrink factor = 0');
%                     axis equal
%                     view([0 90]);
%                     
%                     %Calculate rigid body parameters
%                     data = unique([x, y, z], 'rows');
%                     tri2 = delaunay(data);
%                     TR = triangulation(tri2, data(:, 1), data(:, 2), data(:, 3));                    
%                     [RBP, DT]=RigidBodyParams(TR);
%                     
%                     disp('Wait');
%                 otherwise
%                 validatestring(method, validmethods, 'calculateVolume', 'method');
%             end
%             
%         end
%     end