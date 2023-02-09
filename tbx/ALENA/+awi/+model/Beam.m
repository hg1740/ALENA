classdef (ConstructOnLoad) Beam < awi.model.Connector & awi.mixin.Beamable
    %Beam Defines a line object that has additional properties describing
    %its orientation, stiffness, mass and inertia, and material properties.
    %A Beam object can also aggregate other object belonging to the AWI
    %package, including Compartment and LoadDistribution objects.
    %
    % *Orientation*
    %
    % The orientation of the beam is handled by the 'CoordSys' class. This
    % class defines a generic rectangular coordinate system using a number
    % of parameterisation schemes, for more information see the
    % documentation for the 'CoordSys' class. (TODO - Link to go here)
    %
    % By default, each 'CoordSys' class defines a coordinate system which
    % is aligned with the global coordinate system. (e.g. a 3x3 identity
    % matrix). As new coordinate systems are added to the beam it is
    % necessary to translate each coordinate system so that they are
    % positioned correctly along the beam axis.
    %
    % *Cross-Section*
    %
    % A beam can have a generic cross-section...
    %
    % *Materials*
    % A beam can have a number of materials defined along its length
    %
    % *Beam Properties*
    %
    % 'Beam' is a subclass of the AWI mixin class 'Beamable'. This class
    % allows generic properties to be defined along the length of the beam.
    %
    % During the construction of the 'Beam' class various dynamic
    % properties are added to the object which describe these properties
    % and their distribution along the length of the beam.
    %
    % The following quantities are used to describe the 'Beam' class:
    %
    %   - 'I11'
    %   - 'I22'
    %   - 'I12'
    %   - 'A'
    %   - 'J'
    %   - 'NSM'
    %   - 'NSI'
    %   - 'SCy'
    %   - 'SCz'
    %   - 'NAy'
    %   - 'NAz'
    %   - 'CMy'
    %   - 'CMz'
    %
    % These properties are all of type 'FE BeamProp'.
    
    %Compartments
    properties (SetAccess = private)
        %Handle to any 'awi.model.Compartment' objects
        Compartments = awi.model.Compartment.empty;
        %Names of classes which can be added as a Compartment.
        ValidCompartments = {'awi.model.Compartment'};
    end
    
    %Load distribution
    properties (SetAccess = private)
        %Loads/motions applied to this beam
        AppliedLoads = awi.loads.LoadDistribution.empty;
    end
    
    %Mass from the NSM beam property
    properties (Dependent)
        %Mass associated with the non-structural mass of the beam
        NonStructuralMass;
    end
    
    methods % set / get
        function val = get.NonStructuralMass(obj)
            %get.NonStructuralMass Get method for the dependent property
            %'NonStructuralMass'.
            
            val = [];
            
            %Check if special beam property 'NSM' has been assigned to the
            %hidden beam properties.
            bpNames = {obj.BeamProperties.Name};
            idxNSM  = ismember(bpNames, 'NSM');
            if ~any(idxNSM)
                return
            end
            
            %Grab beam property quantities
            BP  = obj.BeamProperties(idxNSM);
            nsm = getBPV(obj, 'NSM');
            eta = getEta(obj, BP.Type, 'AxisFlag', BP.AxisFlag);
            
            %If the user has not set the value then it will return nan
            if any(isnan(nsm))
                val = 0;
                return
            end
            
            %Length along beam geometry
            xyzLength = obj.RData;
            
            %Get actual length of beam segments
            len = eta .* xyzLength(end);
            dL  = diff(len);
            
            %Calculate mass assuming a linear variation
            %   - i.e. Using the trapezium rule
            if ~strcmp(BP.Variation, 'linear')
                warning(['Unable to calculate ''NonStructuralMass'' for ', ...
                    'mass distributions that are not piecewise-linear.']);
                return
            end
            val = sum(dL .* (nsm(1 : end - 1) + nsm(2 : end)) ./ 2);
            
        end
    end
    
    methods % constructor
        function obj = Beam(varargin)
            %Beam Class constructor for the 'awi.model.Beam' class.
            %
            % Actions performed:
            %   - Call superclass constructor
            %   - Update property groups
            %   - Add (dynamic) beam properties
            %   - Add (dynamic) beam property objects
            
            %Pass it on
            obj@awi.model.Connector(varargin{:});
            
            if isa(obj, 'awi.mixin.FEable') %Update Property Groups
                obj.addPropertyGroup('FE Mesh Props', ...
                    'BeamLength' , 'Beam Element Length'    , 'Defines the length of a FE beam element.', [], ...
                    'NumBeamElem', 'Number of Beam Elements', 'Defines the number of beam elements along the beam length. Overrides ''BeamLength''.', []);
            end
            
            %Add beam properties
            addBeamProperty(obj, 'A'  , 'Type', 'FE BeamProp'); %Cross-sectional area
            addBeamProperty(obj, 'I11',  'Type', 'FE BeamProp'); %Second moment of area in plane 1
            addBeamProperty(obj, 'I22', 'Type', 'FE BeamProp'); %Second moment of area in plane 2
            addBeamProperty(obj, 'I12', 'Type', 'FE BeamProp'); %Second moment of area in plane 1-2
            addBeamProperty(obj, 'J'  , 'Type', 'FE BeamProp'); %Polar second moment of area
            addBeamProperty(obj, 'NSM', 'Type', 'FE BeamProp'); %Non-structural mass (per unit length)
            addBeamProperty(obj, 'NSI', 'Type', 'FE BeamProp'); %Non-structural inertia (per unit length)
            addBeamProperty(obj, 'RI1', 'Type', 'FE BeamProp'); %Rotational inertia about the 1st axis
            addBeamProperty(obj, 'RI2', 'Type', 'FE BeamProp'); %Rotational inertia about the 2nd axis
            addBeamProperty(obj, 'R12', 'Type', 'FE BeamProp'); %Rotational inertia about the 3rd axis
            addBeamProperty(obj, 'SCy', 'Type', 'FE BeamProp'); %Shear centre offset in the y-direction of beam coordinate system
            addBeamProperty(obj, 'SCz', 'Type', 'FE BeamProp'); %Shear centre offset in the z-direction of beam coordinate system
            addBeamProperty(obj, 'NAy', 'Type', 'FE BeamProp'); %Neutral axis offset in the y-direction of beam coordinate system
            addBeamProperty(obj, 'NAz', 'Type', 'FE BeamProp'); %Neutral axis offset in the z-direction of beam coordinate system
            addBeamProperty(obj, 'CMy', 'Type', 'FE BeamProp'); %Centre-of-mass offset in the y-direction of beam coordinate system
            addBeamProperty(obj, 'CMz', 'Type', 'FE BeamProp'); %Centre-of-mass offset in the z-direction of beam coordinate system
            
            %Add stress recovery points relative to the shear center
            addBeamProperty(obj, 'Cy'  , 'Type', 'FE BeamProp'); % Top right corner
            addBeamProperty(obj, 'Cz'  , 'Type', 'FE BeamProp');
            addBeamProperty(obj, 'Dy'  , 'Type', 'FE BeamProp'); % Bottom right corner
            addBeamProperty(obj, 'Dz'  , 'Type', 'FE BeamProp');
            addBeamProperty(obj, 'Ey'  , 'Type', 'FE BeamProp'); % Bottom left corner
            addBeamProperty(obj, 'Ez'  , 'Type', 'FE BeamProp');
            addBeamProperty(obj, 'Fy'  , 'Type', 'FE BeamProp'); % Top left corner
            addBeamProperty(obj, 'Fz'  , 'Type', 'FE BeamProp');
            
              
            
            %Add beam property objects...
            
            %Coordinate systems describing the orientation along the beam
            %    - Must use 'previous' interpolation method. The beam is
            %      assumed to be piecewise linear so each section has the
            %      same oriention until a kink point which is driven by the
            %      geometry of the 'awi.model.Stick' object, not the beam.
            addBeamPropObj(obj, 'awi.model.CoordSys', 'previous', ...
                'Orientation', 'RMatrix_');
            %Material properties along the beam
            %   - Must use 'previous' interpolation method. Using anything
            %     else doesn't make sense as material properties are
            %     assumed to be constant along the material!
            addBeamPropObj(obj, 'awi.model.Material', 'previous', ...
                'Material', 'E', 'G', 'Nu', 'Type', 'FE BeamProp');
            %Generic cross-section along the beam
            %   - Assume a linear variation of cross-sections
            addBeamPropObj(obj, 'awi.model.CrossSection', 'linear', ...
                'CrossSection', 'X', 'Z');
            %Box-Beam cross-sections
            %   - Assume a linear variation in box properties
            addBeamPropObj(obj, 'awi.model.BoxBeam', 'linear', ...
                'BoxBeam', 'Ixx', 'Izz', 'Ixz', 'Jbb', 'Abb');
            
            %A Beam is a primary component in a FE model
            obj.IsFEComponent = true;
            
        end
    end
    
    methods % converting to FE model
        function FEModel = convertThisToFE(obj, FEModel, varargin)
            %convertThisToFE Converts the 'awi.model.Beam' object to a
            %collection of Finite Element (FE) entities.
            %
            % A beam has the following FE entities:
            %   - Coordinate Systems:
            %       * A local coordinate system is defined at each kink
            %       point along the beam.
            %    - Nodes :
            %       * A node is defined at the start and end of each
            %       beam element.
            %       * The number of beam elements is controlled by the
            %       mesh parameters 'NumBeamElem' and 'BeamElemLength'.
            %       * The position of the nodes along the coordinates of
            %       the'awi.model.Beam' object is controlled by the
            %       non-dimensionsal 'eta' positions of the various beam
            %       properties that belong to this 'awi.model.Beam' object.
            %       * Additional arguments can be provided which specify
            %       coordinates along the beam which MUST be included in
            %       the FE model. The code will then make sure that these
            %       nodes are connected by beam elements and the mesh will
            %       adjust accordingly.
            %    - Beams :
            %       * A beam element is created between each node.
            %       * The beam orientation is aligned with the local node
            %       output coordinate system.
            %    - Beam Properties :
            %       * Properties are automatically assigned to each
            %       beam element based on the values that have been set for
            %       the various beam properties that belong to this
            %       'awi.model.Beam' object.
            %       * The beam properties are interpolated using the
            %       'awi.mixin.Beamable' & 'awi.mixin.BeamProperty'
            %       methods. See these class definitions for more details.
            %    - Materials :
            %       * Material objects are defined based on the
            %       distribution of 'awi.model.Material' objects along the
            %       beam.
            %       * The material properties are taken directly from the
            %       underlying material objects.
            
            %Start with the superclass
            FEModel = convertThisToFE@awi.model.Component(obj, FEModel, varargin{:});
            
            %FEM depends on type of discretisation
            switch obj.ModelType
                case '1D'
                    FEModel = generate1DBeamModel(obj, FEModel, varargin{:});
                case '2D'
                    FEModel = generate2DFEModel(obj, FEModel, varargin{:});
            end
            
        end
    end
    
    methods % converting to FE model
        function tf = canConvertToFE(obj)
            %canConvertToFE Check is the Beam object can be converted to a
            %FE representation.
            
            tf = true;
            
            %Must have material defined in order to make FEM
            if isempty(obj.Material) || isempty(obj.Orientation)
               tf = false;
               return
            end
            
        end
    end
    methods (Access = protected) % converting to FE model
        %1D beam model 
        function FEModel = generate1DBeamModel(obj, FEModel, varargin)
            %generate1DBeamModel Generates the 'awi.fe.' data necessary to
            %describe a 1D representation of a Beam object.
            
            %Cannot create the FE equivalent model if we have not defined
            %the orientation or material properties of the beam, or if the
            %beam does not have geometry data.
            [xd, yd, zd] = obj.xyzdata;
            if isempty(xd) || isempty(yd) || isempty(zd)
                i_printWarning(obj.Name, 'beam (x,y,z) geometry');
                return
            end
            if isempty(obj.Orientation)
                i_printWarning(obj.Name, 'Orientation');
                return
            end
            if isempty(obj.Material)
                i_printWarning(obj.Name, 'Material');
                return
            end
            
            function i_printWarning(nam, beamObjName)
                warning(['Unable to create a FE representation of the'   , ...
                    'Beam object ''%s'' as the %s has not been defined.'], ...
                    nam, beamObjName);
            end
                        
            %Parse inputs and return the required coordinates to be included
            %in the beam mesh
            [etaG, xG, yG, zG] = i_parseInputs(obj, varargin{:});
            
            function [etaR, xR, yR, zR] = i_parseInputs(obj, varargin)
                
                p = inputParser;
                addParameter(p, 'etaR', [0, 1], @(x)validateattributes(x, ...
                    {'numeric'}, {'row', 'finite', 'real', 'nonnegative', '<=', 1}));
                addParameter(p, 'xR', obj.XData([1, end]), @(x)validateattributes(x, {'numeric'}, {'row', 'finite', 'real', 'nonnan'}));
                addParameter(p, 'yR', obj.YData([1, end]), @(x)validateattributes(x, {'numeric'}, {'row', 'finite', 'real', 'nonnan'}));
                addParameter(p, 'zR', obj.ZData([1, end]), @(x)validateattributes(x, {'numeric'}, {'row', 'finite', 'real', 'nonnan'}));
                parse(p, varargin{:});
                
                %Select results and send them up
                etaR = p.Results.etaR;
                xR   = p.Results.xR;
                yR   = p.Results.yR;
                zR   = p.Results.zR;
                
            end
            
            %Return the set of coordinates 'R' that will be used to define
            %the FE beam coordinates.
            R = i_combineBeamCoords(obj, etaG, xG, yG, zG);
            
            function R = i_combineBeamCoords(obj, etaR, xR, yR, zR)
                %i_combineBeamCoords Returns the coordinates that will be
                %used to define the FE beam geometry.
                %
                % Combination of:
                %   - Required coordinates (xR, yR, zR) & etaR
                %   - FE property coordinates (xB, yB, zB) & etaB
                %   - Connection point with any 'awi.model.Connector'
                %     objects
                
                bUseFEProp = false;
                
                %Points from beam properties
                if bUseFEProp
                    etaB = getEta(obj, 'FE BeamProp');
                    xB   = interp1(etaR, xR, etaB);
                    yB   = interp1(etaR, yR, etaB);
                    zB   = interp1(etaR, zR, etaB);
                else
                    etaB = [];
                    xB   = [];
                    yB   = [];
                    zB   = [];
                end
                
                %Points from beam property objects
                etaBP = [obj.Orientation_eta, obj.Material_eta];
                if ~isempty(obj.AppliedLoads)
                    etaBP = [etaBP, horzcat(obj.AppliedLoads.EtaDistribution)];
                end
                if isempty(etaBP)
                    xBP = [];
                    yBP = [];
                    zBP = [];
                else
                    etaBP = unique(etaBP);
                    r = obj.RData;
                    r_bar = r ./ r(end);
                    xBP = interp1(r_bar, obj.XData, etaBP);
                    yBP = interp1(r_bar, obj.YData, etaBP);
                    zBP = interp1(r_bar, obj.ZData, etaBP);
                end
                
                %Points from 'awi.fe.Connector' objects
                Con = obj.ChildConnectors;
                if isempty(Con)
                    xC   = [];
                    yC   = [];
                    zC   = [];
                    etaC = [];
                else
                    %Only include data from Connectors that connect
                    %directly to the beam axis (i.e. they have no tip
                    %(x,y,z) offsets.
                    off  = cell2mat(get(Con, ...
                        {'TipXOffset', 'TipYOffset', 'TipZOffset'}));
                    Con  = Con(all(off == 0, 2));
                    if isempty(Con)
                        xC   = [];
                        yC   = [];
                        zC   = [];
                        etaC = [];
                    else
                        etaC = cell2mat(get(Con, {'TipNormPos'}));
                        etaC = etaC(:)';
                        xC   = interp1(etaR, xR, etaC);
                        yC   = interp1(etaR, yR, etaC);
                        zC   = interp1(etaR, zR, etaC);
                    end
                end
                
                %Combine and find unique values based on 'eta'
                eta  = [etaR, etaB, etaC, etaBP];
                x_   = [xR, xB, xC, xBP];
                y_   = [yR, yB, yC, yBP];
                z_   = [zR, zB, zC, zBP];
                [~, ind] = unique(eta);
                
                %Reduced set that will be used for defining beam coords.
                R = [x_(ind) ; y_(ind) ; z_(ind)];
                
                %Final check for uniqueness (some properties have different
                %'eta' axis directions but describe the same coordinates.
                [~, ind] = unique(R', 'rows', 'stable');
                R = R(:, ind);
                
            end
                        
            %Define in the global coordinate system
            R = R + repmat(obj.AbsPosition', [1, size(R, 2)]);
            
            %Calculate the mesh for the beam elements
            [x, y, z, etaNodes] = i_defineBeamMesh(obj, R(1, :), R(2, :), R(3, :));
            
            function [xB, yB, zB, etaB] = i_defineBeamMesh(obj, x, y, z)
                %defineBeamMesh Defines a series of points (xB, yB, zB)
                %along the line through the coordinates (x,y,z)
                
                %No beam less than 0.1% of the segment length
                tol = 1e-3;
                
                %Straight line distance along the beam
                s = awi.model.Stick.getLineLength(x, y, z);
                segLengths = diff(s);           %Length of each segment
                nSeg       = numel(segLengths); %Number of segments
                
                %Calculate length of beam element
                if isempty(obj.NumBeamElem)
                    beamLength = obj.BeamLength;
                else
                    beamLength = s(end) ./ obj.NumBeamElem;
                end
                
                %Determine number of elements along the beam
                %   - Don't bother generating any beams that are going to
                %   less than a certain fraction of the segment.
                temp        = segLengths ./ beamLength;
                idx_        = (temp - floor(temp)) < (tol .* segLengths);
                nBeam       = ceil(temp);
                nBeam(idx_) = floor(temp(idx_));
                
                %Use 'linspace' to define new points
                beamNodes = arrayfun(@(i) linspace(s(i), s(i + 1), nBeam(i) + 1), 1 : nSeg, 'Unif', false);
                
                %Collect points and only keep unique values
                allBeamNodes    = horzcat(beamNodes{:});
                uniqueBeamNodes = uniquetol(allBeamNodes, tol);
                
                %Return the normalised distribution so that beam properties
                %can be interpolated
                etaB = uniqueBeamNodes ./ uniqueBeamNodes(end);
                
                %Go back to (x,y,z) coordinates
                xB = interp1(s, x, uniqueBeamNodes);
                yB = interp1(s, y, uniqueBeamNodes);
                zB = interp1(s, z, uniqueBeamNodes);
                
            end
            
            %Coordinate Systems
            CoordSys = makeFECoordSys(obj);
            
            %Material Properties
            Mat = makeFEMaterials(obj);
            
            %Nodes
            Nodes = makeFENodes(obj, x, y, z, etaNodes, CoordSys);
            
            %Beams
            [Beams, BeamProps] = makeFEBeams(obj, etaNodes, Nodes, Mat);
            
            %Loads
            Loads = makeFELoads(obj, etaNodes, Nodes);
            
            %Define cross-section using GRID and RBE?
            if obj.GenerateCrossSectionNodes
                if isempty(obj.CrossSection)
                    warning(['Unable to generate cross-section nodes ', ...
                        'as the beam object ''%s'' does not have any ', ...
                        'cross-section data.'], obj.Name)
                else
                    [CSNodes, CSRigidBar] = i_generateCrossSectionBulk(obj, etaNodes, FEModel);
                    set(CSRigidBar, {'NodesI'}, num2cell(Nodes)');
                end
            else
                CSNodes    = [];
                CSRigidBar = [];
            end
            
            function [CSNodes, CSRigidBar] = i_generateCrossSectionBulk(obj, etaNodes, FEModel)
                
                %Interpolate the profiles to the beam positions
                
                %Could just do this...
                %
                %NewBeamCS = interpolateCrossSection(obj.CrossSection, etaNodes);
                %[x, y, z] = calculateGlobalCoords(BeamCrossSection);
                %
                %But in order to avoid overheads from generating the
                %objects we will duplicate some code from the
                %'awi.model.CrossSection' object.
                
                [xU, xL, zU, zL, ~] = interpolateXZData(obj.CrossSection, etaNodes);
                rMatrix             = interpolateOrientationData(obj.CrossSection, etaNodes);
                
                %* * * Duplicate code * * * TODO - Tidy this up
                [xd_, yd_, zd_] = xyzdata(obj);
                r_    = awi.model.Stick.getLineLength(xd_, yd_, zd_);
                eta_  = r_ ./ r_(end);
                xyzC  = interp1(eta_, [xd_ ; yd_ ; zd_]', etaNodes);
                
                %Swap coordinate system (see also
                %awi.model.CrossSection.getGlobalCoords)
                y_ = [xU, fliplr(xL)];
                z_ = [zU, fliplr(zL)];
                x_ = zeros(size(y_));
                
                coords = arrayfun(@(i) ...
                    rMatrix{i} * [x_(i, :) ; y_(i, :) ; z_(i, :)], 1 : numel(etaNodes), 'Unif', false);
                coords = cat(3, coords{:});
                coords = permute(coords, [3, 2, 1]);
                x_ = coords(:, :, 1);
                y_ = coords(:, :, 2);
                z_ = coords(:, :, 3);
                               
                x_ = x_ + xyzC(:, 1);
                y_ = y_ + xyzC(:, 2);
                z_ = z_ + xyzC(:, 3);
                
                pos_ = obj.AbsPosition;
                x_   = x_ + pos_(1);
                y_   = y_ + pos_(2);
                z_   = z_ + pos_(3);
                
                %Check for duplicate data and if so trim the final node
                if isequal([x_(:, 1), y_(:, 1), z_(:, 1)], [x_(:, end), y_(:, end), z_(:, end)])
                    x_ = x_(:, 1 : end - 1);
                    y_ = y_(:, 1 : end - 1);
                    z_ = z_(:, 1 : end - 1);
                end
                
                [nCS, nNodePerCS] = size(x_);
                
                CSNodes    = arrayfun(@(~) awi.fe.Node    , 1 : numel(x_));
                CSRigidBar = arrayfun(@(~) awi.fe.RigidBar, 1 : nCS);
                set(CSRigidBar, 'CN', 123456); %Transfer all translation/rotation

                ub = cumsum(repmat(nNodePerCS, [1, nCS]));
                lb = [1, ub(1 : end - 1) + 1];
                for iCS = 1 : nCS
                    ind = lb(iCS) : ub(iCS);
                    set(CSNodes(ind)   , {'X'}, num2cell([x_(iCS, :) ; y_(iCS, :) ; z_(iCS, :)], 1)');
                    set(CSRigidBar(iCS), 'NodesD', CSNodes(ind)');
                    addPart(FEModel, ['CrossSectionNodes_', num2str(iCS)], ...
                        CSNodes(ind));
                end
                
            end
            
            %Set up parts
            addPart(FEModel, 'BeamNodes', Nodes);
            addPart(FEModel, 'Beams'    , Beams);
            addPart(FEModel, 'BeamProps', BeamProps);
            addPart(FEModel, 'Materials', Mat);
            addPart(FEModel, 'Loads'    , Loads);
            
            %Add the objects to the FE model
            addFEData(FEModel, CoordSys, Nodes, Beams, BeamProps, Mat, Loads, CSNodes, CSRigidBar);
            
            %PointMasses
            PointMass = FEModel.PointMasses; %Any point masses been generated at Component level?
            if ~isempty(PointMass)
                
                %Which point masses have not been allocated a node?
                idx = cellfun(@(x) isempty(x), {PointMass.Node});
                PointMass = PointMass(idx);
                
                %Grab corresponding 'awi.model.PointMass' objects
                %   - Mass objects are collected in 'awi.model.PointMasses'
                %     objects.
                b = arrayfun(@(x)isa(x, 'awi.model.PointMasses'), obj.Children);
                if ~any(b)
                    return
                end
                mass = vertcat(obj.Children(b).Children);
                mass = mass(idx);
                
                %What are their positions?
                pos = vertcat(mass.AbsPosition);
                
                %Assign the 'awi.fe.Node' objects
                assignMassNodes(FEModel, PointMass, pos);
                
            end
                  
        end
        function CoordSys = makeFECoordSys(obj, cs)
            %makeFECoordSys Creates the 'awi.fe.CoordSys' objects using the
            %'awi.model.CoordSys' objects that are contained in the
            %'Orientation' property of this object or are provided in the
            %variable 'cs'.
            
            if isempty(obj.Orientation) && nargin < 2
                return
            end
            
            if nargin < 2
                %Always disregard the last orientation entry as the orientation
                %properties use the 'previous' interpolation method.
                cs = obj.Orientation(1 : end - 1);
            end
            
            nCS = numel(cs);
            
            %Make the coordinate system
            CoordSys = arrayfun(@(~) awi.fe.CoordSys, 1 : nCS, 'Unif', false);
            CoordSys = horzcat(CoordSys{:});
            
            %TODO - Really we want to define B and C such that it will
            %yield the original eX, eY & eZ vectors of the CoordSys
            %object. Right now the axis is flipped because of pesky matrix
            %algebra.
            
            %Grab data
            A  = (vertcat(cs.AbsPosition) + obj.AbsPosition)'; %Origin - Account for beam position
            r  = horzcat(cs.RMatrix_);     %Column format of 3x3 matrix --> [eX ; eY ; eZ]
            B  = A + r([7, 8, 9], :);      %Point along the z-axis
            C  = A + r([1, 2, 3], :);      %Point along the x-axis
            
            %Assign data
            set(CoordSys, {'A'}, num2cell(A, 1)');
            set(CoordSys, {'B'}, num2cell(B, 1)');
            set(CoordSys, {'C'}, num2cell(C, 1)');
            
        end
        function Nodes = makeFENodes(obj, x, y, z, etaNodes, CoordSys)
            %makeFENodes Creates the 'awi.fe.Node' objects using the
            %(x,y,z) coordinates defined by 'x', 'y' & 'z'. Optionally
            %stores a reference to the output coordinate system defined by
            %the 'Orientation' property.
            
            %Create the objects
            Nodes  = arrayfun(@(~) awi.fe.Node, 1 : numel(x), 'Unif', false);
            Nodes  = horzcat(Nodes{:});
            
            %Assign data - Always assume the coordinates are defined in the
            %basic coordinate system.
            set(Nodes, 'CP' , 0);
            set(Nodes, {'X'}, num2cell([x ; y ; z], 1)');
            
            %Any orientation defined?
            if nargin <  6
                return
            end
            
            %Assign a reference to the output coordinate system
            idx = arrayfun(@(i) and(etaNodes >= obj.Orientation_eta(i), ...
                etaNodes < obj.Orientation_eta(i + 1)), ...
                1 : (numel(obj.Orientation) - 1) , 'Unif', false);
            idx{end}(end) = true; %The last node always has the final output coordinate system
            
            %Assign the material to the beam properties
            for i = 1 : (numel(obj.Orientation) - 1)
                set(Nodes(idx{i}), 'OutputCoordSys', CoordSys(i));
            end
            
        end
        function [Beams, BeamProps] = makeFEBeams(obj, etaNodes, Nodes, Mat)
            %makeFEBeams Creates the 'awi.fe.Beam' objects using the
            %'awi.fe.Node' objects and the beam properties of the
            %'awi.model.Beam'.
            
            %How many objects are required?
            nBeams = numel(etaNodes) - 1;
            
            %Determine how the beam properties are to be defined
            if obj.UseBeamCrossSection && ~isempty(obj.BoxBeam)
                bpNames = {'NSM', 'SCy', 'SCz'};
                BeamProps = arrayfun(@(~) awi.fe.BeamCrossSection, 1 : nBeams, 'Unif', false);
            else
                bpNames = {'A', 'I11', 'I22', 'I12', 'J', 'NSM', 'NSI', ...
                    'SCy', 'SCz', 'CMy', 'CMz', 'NAy', 'NAz',...
                    'Cy', 'Cz', 'Dy', 'Dz', 'Ey', 'Ez', 'Fy', 'Fz'};
                BeamProps = arrayfun(@(~) awi.fe.BeamProp, 1 : nBeams, 'Unif', false);
            end
            
            %Interpolate the beam properties to the beam element positions
            bpData  = cellfun(@(x) getBPV(obj, x, etaNodes, [], ...
                'AxisFlag', 'R'), bpNames, 'Unif', false);
            
            %Remove NaN data
            function matrix = replace_nan(matrix)
                matrix(isnan(matrix)) = 0;
            end
            bpData = cellfun(@replace_nan, bpData, 'Unif', false);
            
            %Check the beam data to see if we can remove any duplicate data
            % TODO
            
            %Reshape the data into the format required by 'awi.fe.BeamProp'
            %    i.e. [END-A ; END-B]
            bpData   = cellfun(@(x) [x(1 : end - 1) ; x(2 : end)], bpData, 'Unif', false);
            etaBeams = [etaNodes(1 : end - 1) ; etaNodes(2 : end)];
            
            %Make the initial objects
            Beams     = arrayfun(@(~) awi.fe.Beam    , 1 : nBeams, 'Unif', false);
            Beams     = horzcat(Beams{:});
            BeamProps = horzcat(BeamProps{:});
            
            %Assign data
            set(Beams, {'Nodes'}       , num2cell([Nodes(1 : end - 1) ; Nodes(2 : end)], 1)');
            set(Beams, {'BeamProperty'}, num2cell(BeamProps)');
            set(Beams, 'X'             , [0 ; 1 ; 0]);
            set(Beams, {'SCy'}         , num2cell(bpData{ismember(bpNames, 'SCy')}, 1)');
            set(Beams, {'SCz'}         , num2cell(bpData{ismember(bpNames, 'SCz')}, 1)');
            
            %Strip 'SCy' & 'SCz' from the list of beam property names
            idx = ismember(bpNames, {'SCy', 'SCz'});
            bpNames(idx) = [];
            bpData(idx)  = [];
            
            %Beam Properties
            for i = 1 : numel(bpNames)
                set(BeamProps, bpNames(i), num2cell(bpData{i}, 1)');
            end
            
            %Assign the dimensions and type of the cross-section
            if obj.UseBeamCrossSection && ~isempty(obj.BoxBeam)
                i_assignCrossSectionData(obj, BeamProps, etaNodes);
            end
            
            %Determine which material belongs to which beam property
            idx = arrayfun(@(i) and(etaBeams(1, :) >= obj.Material_eta(i), ...
                etaBeams(1, :) < obj.Material_eta(i + 1)), 1 : numel(Mat) , ...
                'Unif', false);
            
            %Assign the material to the beam properties
            for i = 1 : numel(Mat)
                set(BeamProps(idx{i}), 'Material', Mat(i));
            end
            
            function i_assignCrossSectionData(obj, BeamCrossSection, etaBeam)
                %i_assignCrossSectionData Assigns the type and dimensions
                %of the beam cross-section based on the 'BoxBeam' objects
                %defined along the span.
                
                %Grab the 'BoxBeam' objects
                BoxBeam = [obj.BoxBeam];
                
                %Determine the type of the 'BoxBeam'
                type = unique({BoxBeam.BoxType});
                assert(numel(type) == 1, sprintf(['Expected all of the ' , ...
                    'Box-Beam objects attached to the Beam ''%s'' to be ', ...
                    'of the same type.' ], obj.Name));
                
                %Define the MSC.Nastran beam library type & no. dimensions
                idx_ = ismember(BoxBeam(1).ValidBoxTypes, type);
                nPrp = numel(BoxBeam(1).BoxProperties{idx_});
                code = BoxBeam(1).NastranBoxCode{idx_};
                set(BeamCrossSection, 'TYPE', code);
                
                d = zeros(nPrp, numel(etaBeam));
                
                %Interpolate the 'Height'  & 'Width' values
                d(1, :) = interp1(obj.BoxBeam_eta, [BoxBeam.Width] , etaBeam);
                d(2, :) = interp1(obj.BoxBeam_eta, [BoxBeam.Height], etaBeam);
                
                %Interpolate remaining dimensions
                switch type{:}
                    case 'SymmetricBox'
                        d(3, :) = interp1(obj.BoxBeam_eta, [BoxBeam.CoverThickness], etaBeam);
                        d(4, :) = interp1(obj.BoxBeam_eta, [BoxBeam.SparThickness] , etaBeam);
                    otherwise
                        error(['Unknown BoxType for the Box-Beam ' , ...
                            'objects attached to the Beam ''%s''.'], obj.Name);
                end
                
                %Rearrange into the correct format (END-A ; END-B)
                dim = arrayfun(@(i) cat(3, d(i, 1 : end - 1), d(i, 2 : end)), 1 : nPrp, 'Unif', false);
                dim = permute(vertcat(dim{:}), [3, 1, 2]);
                dim = arrayfun(@(i) dim(:, :, i), ...
                    1 : numel(BeamCrossSection), 'Unif', false)';
                
                %Assign to the object
                set(BeamCrossSection, {'Dimensions'}, dim);
                
            end
            
        end
        function Loads = makeFELoads(obj, etaNodes, Nodes)
            %makeFELoads Creates the 'awi.fe.Load' objects using load
            %distribution defined by 'obj.AppliedLoads'.
            
            Loads = [];
            
            %Tolerance for matching loads with nodes and accounting for
            %floating point error
            tol = 1e-5;                       
            
            %Anything?
            ActiveLoads = obj.AppliedLoads;
            if isempty(ActiveLoads)
                return
            end
            
            %Where are loads applied and which nodes are the loads
            %associated with?
            etaLoads = horzcat(ActiveLoads.EtaDistribution);  
            ind = find(any((abs(etaNodes - etaLoads') < tol), 1));
            
            %Orientation coordinate system
            coordSys = arrayfun(@(L) repmat({L.PointLoadCoordSys}, ...
                [numel(L.EtaDistribution), 1]), ActiveLoads, 'Unif', false);
            coordSys = horzcat(coordSys{:});
            
            %Global or follower?
            type = arrayfun(@(L) repmat({L.PointLoadBehaviour}, ...
                [numel(L.EtaDistribution), 1]), ActiveLoads, 'Unif', false);
            type = horzcat(type{:});
            
            %Grab force/moment data
            pl = horzcat(ActiveLoads.PointLoads);
            f  = pl(1 : 3, :, 1);
            m  = pl(4 : 6, :, 1); 
            
            %Check if load data has been populated
            idxF     = any(f);
            idxM     = any(m);
            f        = f(:, idxF);
            m        = m(:, idxM);
            
            %Make the force and moment objects
            Forces  = arrayfun(@(~) awi.fe.PointLoad, 1 : nnz(idxF));
            Moments = arrayfun(@(~) awi.fe.PointLoad, 1 : nnz(idxM));     
            
            %Assign data
            i_assignLoadData(Forces , f, Nodes(ind(idxF)), coordSys(idxF), type(idxF));
            i_assignLoadData(Moments, m, Nodes(ind(idxM)), coordSys(idxM), type(idxM));            
            set(Forces , 'LoadType', 'force');
            set(Moments, 'LoadType', 'moment');
            
            %Return
            Loads = [Forces, Moments];
            
            function i_assignLoadData(Load, loadData, Nodes, coordSys, type)
                %i_assignLoadData Assigns the node, load magnitude and load
                %orientation data to the 'Load' objects.

                %Magnitude and orientation of the load
                mag    = vecnorm(loadData);
                orient = loadData ./ mag; 
                set(Load,  ...
                    {'Node'}         , num2cell(Nodes)'    , ...
                    {'Magnitude'}    , num2cell(mag)'      , ...
                    {'Orientation'}  , num2cell(orient, 1)', ...
                    {'LoadCoordSys'} , coordSys            , ...
                    {'LoadBehaviour'}, type);
                
            end
            
        end
        function Mat = makeFEMaterials(obj)
            %makeFEMaterials Creates the 'awi.fe.Material' objects using
            %the 'awi.model.Material' objects that are contained in the
            %'Material' property of this object.
            
            %Check for materials
            if isempty(obj.Material) || numel(obj.Material) == 1
                Mat = [];
                warning(['No material properties assigned for the ', ...
                    'object ''%s''. Unable to generate FEM.'], obj.Name);
                return
            end
            %             assert(~isempty(obj.Material), sprintf(['No material '   , ...
            %                 'properties assigned for the object ''%s''. Unable ', ...
            %                 'to generate FEM.'], obj.Name));
            
            %Always disregard the last material entry as the material
            %properties use the 'previous' interpolation method.
            nMat = numel(obj.Material) - 1;
            
            %Make the 'awi.fe.Material' objects
            Mat = arrayfun(@(i) awi.fe.Material, 1 : nMat, 'Unif', false);
            Mat = horzcat(Mat{:});
            
            %Populate the data
            matProps = {'E', 'G', 'Nu', 'Rho'};
            val = get(obj.Material(1 : end - 1), matProps);
            set(Mat, matProps, val);
            
        end
        %2D panel element model
        function FEModel  = generate2DFEModel(obj, FEModel, varargin)
            
            %Something like
            %   - Get cross-sections
            %   - Inteprolate to get mesh 
            %   - Defines nodes
            %   - Define panel elements
            %   - Add materials
            %   - Go from there...
            error('Update code');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht, tag)
            %drawElement Visualisation method for the 'awi.model.Beam'
            %class.
            
            %Caller supplied tag explicitly ?
            if nargin < 3
                
                %No - apply a tag indicating that we are just drawing some lines
                tag = 'Beam Offsets';
                
            end
            
            %Pass it on
            hg = drawElement@awi.model.Stick(obj, ht);
            
            %Draw the coordinate systems (if we have any)
            if ~isempty(obj.Orientation)
                for iCS = 1 : numel(obj.Orientation)
                    hg{end + 1} = drawElement(obj.Orientation(iCS), ht); %#ok<AGROW>
                end
            end
            
            %Draw the cross-section (if we have any)
            if ~isempty(obj.CrossSection)
                hg{end + 1} = drawElement(obj.CrossSection, ht);
            end
            
            %Draw the box-beam objects (if we have any)
            if ~isempty(obj.BoxBeam)
                %                 hg{end + 1} = drawElement(obj.BoxBeam, ht, 'global');
            end
            
            %Draw compartments (if we have any)
            if ~isempty(obj.Compartments)
                drawElement(obj.Compartments, ht);
            end
            
            %Draw loads
            if ~isempty(obj.AppliedLoads)
                drawElement(obj.AppliedLoads, ht);
            end
            
            %Transform offsets from the local to the global coord-system
            
            %Add offsets to beam position
            
            %Plot
            
        end        
        function hg = drawCrossSections(obj, ht, varargin)
            %drawCrossSections 
            
            hg = [];
            
            Crs = [obj.CrossSection];
            if isempty(Crs)
                return
            end            
            
            x   = vertcat(Crs.X)';
            z   = vertcat(Crs.Z)';
            nam = strrep(matlab.lang.makeUniqueStrings({Crs.Name}), '_', ' ');
            
            %Plot the full-scale profiles
            hg = plot(ht, x, z);
            set(hg, {'DisplayName'}, nam(:));            
            
            %TODO - Add option to plot 3D coordinates and plot normalised
            %coordinates
           
            
        end
    end
    
    methods % class building
        function updateModel(obj)
            %updateModel Runs any final update that must take place before
            %exiting the build of this object.
            %
            % Actions performed:
            %   1. Makes the coordinate system objects along the beam.
            
            %Make the coordinate system objects - Can be overridden at
            %subclass level.
            [cs, eta] = createCoordSysObjects(obj);
            
            %Assign the coordinate systems to the beam as
            %'BeamPropertyObjects'.
            obj.replaceBeamObject(cs, eta);
            
            %Update the compartments - TODO
            updateGeometry(obj.Compartments);
            
            %Update cross-sections
            createCrossSectionObjects(obj);   
            
        end
    end
    
    methods % analytic methods
        function [cs, eta] = createCoordSysObjects(obj, dG2, varargin)
            %createCoordSysObjects Creates the 'awi.model.CoordSys' objects
            %along the span of the 'Beam'.
            %
            % Builds the coordinate systems along the beam using the 'gSet'
            % parameterisation methods - Assume that the calling function
            % is providing a matrix of offset coordinates for calculating
            % the position of G2.
            
            if nargin < 2 %Escape route
                return
            end
            
            %Allow specification of 'GridVector' & 'GridPlane' values
            p = inputParser;
            addParameter(p, 'GridVector', 'x' , @(x)validatestring(x, {'x', 'y', 'z'}));
            addParameter(p, 'GridPlane' , 'xy', @(x)any(validatestring(x, {'xy', 'xz', 'yz'})));
            parse(p, varargin{:});
            
            assert(all(size(dG2) == [3, numel(obj.XData)]), ['The offset ', ...
                'coordinates matrix must match the size of the '     , ...
                'beam coordinates.']);
            
            %Where do we actually need to place coordinate systems...
            index = findChangeInOrientation(obj);
            
            %What are the eta positions where the coordinate systems will
            %be placed?
            eta = obj.Eta_(index);
            
            %How many kinks?
            nKink = numel(eta);
            
            %Coordinates of the kink points will be the G1 coordinates for
            %the 'gSet' parameterisation.
            kinkCoords = s2pos(obj, eta, 'R');
            
            %Calculate unit orientation vector of the beam for each
            %piecewise linear segment
            v     = diff(kinkCoords, [], 1);
            v_mod = sqrt(sum(v.^2, 2));
            v     = v ./ repmat(v_mod, [1, 3]);
            
            %Make sure the coordinate system axis defined by 'GridVector'
            %points along the beam
            G1 = [kinkCoords(2 : end , :) ; kinkCoords(end, :) + v(end, :)];
            
            %The coordinate system given 'GridPlane' will be in the same
            %plane as a vector from the kink positions to the node defined
            %by the offset vector.
            G2 = kinkCoords +  dG2(:, index)';
            
            %Make the coordinate system objects
            cs = arrayfun(@(i) awi.model.CoordSys, 1 : nKink, 'Unif', false);
            cs = horzcat(cs{:});
            
            %Assign data
            set(cs, 'GridVector', p.Results.GridVector);
            set(cs, 'GridPlane' , p.Results.GridPlane);
            set(cs, 'ActiveSet' , 'gSet');
            set(cs, {'XOffset'} , num2cell(kinkCoords(:, 1)));
            set(cs, {'YOffset'} , num2cell(kinkCoords(:, 2)));
            set(cs, {'ZOffset'} , num2cell(kinkCoords(:, 3)));
            set(cs, {'G1'}      , num2cell(G1, 2));
            set(cs, {'G2'}      , num2cell(G2, 2));
            
            %Build the object
            arrayfun(@(i) build(cs(i)), 1 : nKink);
            
        end
        function index = findChangeInOrientation(obj, tol)
            %findChangeInOrientation Returns the logical index of positions
            %along the beam where there is a change in the orientation of
            %the beam, subjec to a tolerance.
            
            if nargin < 2
                tol = 1e-10;
            end
            
            %Gather coordinates
            r = [obj.XData ; obj.YData ; obj.ZData];
            
            %Calculate unit-vector for each piecewise linear segment
            v  = diff(r, [], 2);
            v_ = sqrt(sum(v.^2));
            v  = v ./ v_;
            v  = [v, v(:, end)];
            
            %Change in orientation vector - The (x,y,z) data is provided
            %along the beam so we need to search left to right.
            dv = abs(fliplr(diff(fliplr(v), [], 2)));
            
            %Logical indexing to find changes in orientation > 'tol'
            idx = any(dv > tol);
            
            %Always include the final eta position!
            index = [true, idx];
            index(end) = true;
            
        end
        function [crs, eta] = createCrossSectionObjects(obj, varargin)
            %createCrossSectionObjects Rebuilds the cross-section objects
            %assigned to this beam.
            %
            % Specific updating of CrossSection objects to be done at the
            % subclass level.
            
            %Sensible defaults
            crs = obj.CrossSection;
            eta = [];
            
            assert(numel(obj) == 1, ['Method ''createCrossSectionObjects'' ', ...
                'is not valid for object arrays.']);
                        
        end
    end
    
    methods % handling 'awi.model.BoxBeam' objects
        function createBoxBeamObjects(obj, boxType, boxDims, boxHeightMethod)
            %generateBoxBeamObjects Creates a set of 'awi.model.BoxBeam' objects along the beam
            
            assert(numel(obj) == 1, ['Method ''createBoxBeamObjects'' ', ...
                'is not valid for handle arrays.']);
            
            p = inputParser;
            addRequired(p, 'boxType');
            addRequired(p, 'boxDims');
            addRequired(p, 'boxHeightMethod');
            parse(p, boxType, boxDims, boxHeightMethod);
            
            %Get the length of the beam (segment-wise linear distance long axis)
            [xd, yd, zd] = xyzdata(obj);
            
            %As 'awi.model.BoxBeam' objects are a beam property object we
            %need to define the normalised 'eta' positions along the beam.
            %   - Define them based on the rib pitch
            r = awi.model.Stick.getLineLength(xd, yd, zd);
            eta = [r(1) : obj.RibPitch : r(end), r(end)] ./ r(end);
            
            %Make the objects and assign them to the 'BoxBeam' property of
            %part
            BoxBeam = arrayfun(@(~) awi.model.BoxBeam, 1 : numel(eta));
            assignBeamObject(obj, BoxBeam, eta, 'replace');
            
            %Set thicknesses & box type
            set(BoxBeam, 'BoxType', boxType);
            prpName = BoxBeam(1).CurrentPropNames;
            prpName(ismember(prpName, {'Width', 'Height'})) = [];
            set(BoxBeam, prpName, num2cell(boxDims));
            
            %Update the width and height of the box beam based on the spar
            %locations and the spar height method
            updateBoxBeamDimensions(obj, 'SparHeightMethod', boxHeightMethod);
            
            %Calculate new sectional properties
            getGeometricProps(obj.BoxBeam);
            
        end
        function updateBoxBeamDimensions(obj, varargin)
            %updateBoxBeamDimensions Updates the 'Width' and 'Height'
            %parameters of the 'BoxBeam' objects based on the cross-section
            %of the beam.
            
            error('Update code to use new interpolation methods in the CrossSection object');
            
            assert(numel(obj) == 1, ['Method ''updateBoxBeamDimensions'' ', ...
                'is not valid for handle arrays.']);
            
            if isempty(obj.BoxBeam) %Escape route
                return
            end
            
            p = inputParser;
            addParameter(p, 'SparHeightMethod', 'average', @(x)any(validatestring(x, ...
                {'frontspar', 'rearspar', 'max', 'average', 'torenbreek'})));
            parse(p, varargin{:});
            
            %If the object is of type 'awi.model.LiftingSurface' then look
            %for 'Spar' objects. If none have been defined then assume the
            %box-beam runs from 15% to 65% of the chord
            if isa(obj, 'awi.model.LiftingSurface')
                sp = findall(obj, 'Type', 'Spar');
                if isempty(sp) || numel(sp) < 2
                    %
                    xLoc = [0.15, 0.65];
                else
                    error('Update code');
                end
                
                %Calculate the chord length at the 'eta' positions and
                %calculate the box heights at the front spar, rear spar and
                %the maximum height
                etaBB  = obj.BoxBeam_eta;
                chrdBB = interp1(obj.Chord_eta, obj.Chord, etaBB);
                swpBB  = interp1(obj.Sweep_eta, obj.Sweep, etaBB, 'previous');
                width  = chrdBB .* abs(diff(xLoc)) .* cosd(swpBB);
                
                nBB = numel(etaBB);
                nSp = numel(xLoc);
                
                %Interpolate the normalised profile data to the desired
                %'eta' positions.
                [xUpper, xLower, zUpper, zLower, etaInterp] = ...
                    interpolateProfileData(NormCoords.X, NormCoords.Z, ...
                    obj.Aerofoil_eta, etaBB, obj.NumAerofoilPointsX);
                
                %Scale by chord
                chrdBB_ = chrdBB';
                xSpar   = chrdBB_ * xLoc;
                xUpper  = xUpper .* chrdBB_;
                xLower  = xLower .* chrdBB_;
                zUpper  = zUpper .* chrdBB_;
                zLower  = zLower .* chrdBB_;
                
                %Preallocate
                zSparU = zeros(nBB, nSp);
                zSparL = zeros(nBB, nSp);
                
                %Calculate the height of the wing at the spar locations and
                %calculate max height
                % - TODO : Vectorise this!!
                for i = 1 : nBB
                    zSparU(i, :) = interp1(xUpper(i, :), zUpper(i, :), xSpar(i, :));
                    zSparL(i, :) = interp1(xLower(i, :), zLower(i, :), xSpar(i, :));
                end
                hSpar = zSparU - zSparL;
                hMax  = max(zUpper - zLower, [], 2);
                
                %Define the box height based on the 'SparHeightMethod'
                switch  p.Results.SparHeightMethod
                    case 'frontspar'
                        height = hSpar(:, 1);
                    case 'rearspar'
                        height = hSpar(:, end);
                    case 'max'
                        height = hMax;
                    case 'average'
                        height = mean(hSpar, 2);
                    case 'torenbreek'
                        height = ((2 .* hSpar(:, 1) + 2 .* hMax + 2 .* hSpar(:, end)) ./ 4) ./ 2;
                end
                
                %Set the position of the local box-beam origin so that it
                %coincides with the front-spar when plotted in global CS.
                xB  = interp1(obj.Eta_, obj.BeamLoc_i, etaBB)';
                dX  = xSpar(:, 1) - (xB .* chrdBB_);
                dZ  = -height ./ 2;
                org = [dX, dZ]';
                
                %Assign data to the 'awi.model.BoxBeam' objects
                set(obj.BoxBeam, {'Height'}, num2cell(height));
                set(obj.BoxBeam, {'Width'} , num2cell(width'));
                set(obj.BoxBeam, {'SectionOrigin'}, num2cell(org, 1)');
                
            else
                error(['Unable to update the box-beam geometry for ', ...
                    'objects that are not of class ', ...
                    '''LiftingSurface''. Update code.']);
            end
            
            
        end
    end
    
    methods % handling 'awi.model.Compartment' objects
        function varargout = addCompartment(obj, Comp, varargin)
            %addCompartment Adds a compartment object to this beam and
            %assigns any additional parameters in 'varargin' to the object.
            
            %Parse
            validateattributes(Comp, {'awi.model.Compartment'}, {'vector'}, ...
                class(obj), 'Compartment');
            if iscolumn(Comp)
                Comp = Comp';
            end
            
            %Assign values to object
            %prp    = varargin(1 : 2 : end);
            %pNames = properties(Comp);
            %idx    = ismember(pNames, prp);
            %prp_   = pNames(idx);
            %index  = find(ismember(prp, prp_)) * 2;
            %val    = varargin(index);
            %set(Comp, prp_, val);
            set(Comp, varargin{:});
            
            %Assign a reference to this 'awi.model.Beam' object
            set(Comp, 'BeamHandle', obj);
            
            %Add to the collection
            obj.Compartments = [obj.Compartments, Comp];
            
            if nargout == 1
                varargout{1} = Comp;
            end
            
        end
        function removeCompartment(obj, Comp)
            %removeCompartment Removes a particular compartment object from
            %the current collection.
            
            idx = ismember(obj.Compartments, Comp);
            obj.Compartments = obj.Compartments(~idx);
            
        end
        function removeAllCompartments(obj)
            %removeAllCompartments Removes all compartments from the beam
            %object but does not delete the handle.
            
            %Simple
            obj.Compartments = [];
            
        end
        function deleteCompartments(obj)
            %deleteCompartments Removes the compartments from the beam
            %object and deletes the current 'awi.model.Compartment'
            %objects.
            
            c = obj.Compartments;
            removeCompartments(obj);
            delete(c);
            
        end
    end
    
    methods % handling 'awi.loads.LoadDistribution' objects
        function varargout = assignLoadToBeam(obj, Load)
            %assignLoadToBeam Assigns the load distribution to a beam
            %object and vice-versa.
            
            %Parse
            validateattributes(Load, {'awi.loads.LoadDistribution'}, {'vector'}, ...
                class(obj), 'Load');
            if iscolumn(Load)
                Load = Load';
            end
            
            if nargout == 1
                varargout{1} = Load;
            end
            
            %Had the load already been assigned to another beam? If so
            %remove the current association and replace with this beam
            hBeam = Load.BeamHandle;
            if ~isequal(obj, hBeam)
                %This will recursively call 'assignLoadToBeam' on the Beam
                %object so we can quit out after this
                assignLoadToBeam(Load, obj);
                return
            end
            
            %Add the load to the collection but check it first
            parse(Load);
            obj.AppliedLoads = [obj.AppliedLoads, Load];
            
        end
        function removeLoadFromBeam(obj, Load)
            %removeLoadFromBeam Removes a particular LoadDistribution from
            %the 'Beam' object.
            
            idx = ismember(obj.AppliedLoads, Load);
            obj.AppliedLoads = obj.AppliedLoads(idx);
            
        end
        function removeAllLoads(obj)
            %removeAllLoads Removes all loads from the Beam object.
            
            if isempty(obj.AppliedLoads)
                return
            end
            
            %Break association between load and beam
            assignLoadToBeam(obj.AppliedLoads, []);
            
            %Remove from beam
            obj.AppliedLoads = [];
            
            
        end
    end
    
    methods (Access = protected) % handling 'awi.model.Compartment' objects
        function updateCompartmentType(obj, compType)
            %updateCompartmentType Adds the classnames specified in
            %'compType' to the list of allowable compartments specified in
            %'ValidCompartments'.
            
            %Parse
            if ~iscell(compType)
                compType = {compType};
            end
            assert(iscellstr(compType), ['Expected the list of ', ...
                'allowable compartments to be a cell-array of ' , ...
                'character vectors']);
            if isrow(compType)
                compType = compType';
            end
            
            %Check we are defining valid 'awi.model.Compartment' types
            func = cellfun(@str2func, compType, 'Unif', false);
            o    = cellfun(@(x) x(), func, 'Unif', false);
            idx  = cellfun(@(x) isa(x, 'awi.model.Compartment'), o);
            compType = compType(idx);
            
            %Update the list
            obj.ValidCompartments = [obj.ValidCompartments ; compType];
            
        end
    end
    
end

