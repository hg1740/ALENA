classdef SlenderBody < matlab.mixin.SetGet
    %SlenderBody Describes the basic geometry and system matrices of a
    %slender body in the Intrinsic, Strain-Based, Finite Element (ISBFE)
    %beam formulation.
    %
    % Detailed Description:
    %   * This object describes only a single slender body (which itself is
    %     a collection of beam elements). To represent an aircraft multiple
    %     objects must be defined and their system matrices combined.
    %   * The properties are split into two categories:
    %       1. Public facing - These properties are visible to the user and
    %       are validated each time they are changed. All public facing
    %       properties have an equivalent internal property which can be
    %       set by adjusting the public facing value.
    %       2. Internal - These are the actual property values that are
    %       used by the solver. These properties are hidden and are not
    %       validated each time they are set, therefore they are intended
    %       for use by an advanced user only.
    
    %System states (public)
    properties
        %Forces and moments in the local (intrinsic) beam coordinate sytem
        LocalBeamLoads
        LocalVelocities
        GlobalVelocities
    end
    
    %System states (internal)
    properties (Hidden, SetAccess = private)
        %Local beam loads [6 * nElem]
        x_f
        %Velocities in the global frame [6, 1]
        x_vg
        %Velocities in the aircraft frame [6, 1]
        x_va
    end
    
    %System matrices (public)
    properties
        %Beam element mass matrix (excluding FE shape functions)
        MassMatrix
        %Beam element stiffness matrix (excluding FE shape functions)
        StiffnessMatrix
        %Beam element compliance matrix (excluding FE shape functions)
        ComplianceMatrix
        %Direction vector (e1)
        ElementDirectionVector
        %Loads applied in the local coordinate system at the FE nodes
        LocalNodalLoads
        %Loads applied in the aircraft coordinate system at the FE nodes
        AircraftNodalLoads
        %Loads applied in the local coordinate system at the midpoint of
        %the FE beam elements
        LocalElementLoads
    end
    
    %System matrices (internal)
    properties (Hidden, SetAccess = private)
        %Beam element mass matrix (excluding FE shape functions)
        M
        %Beam element stiffness matrix (excluding FE shape functions)
        K
        %Beam element compliance matrix (excluding FE shape functions)
        C
        %Full compliance matrix (block diagonal of element compliance)
        C_full
        %Full matrix of strain offset terms for evaluating loads about
        %shear axis
        CFFbar_full
        %Direction vector (e1)
        E1
        %Loads applied in the local coordinate system at the FE nodes
        f_n
        %Loads applied in the aircraft coordinate system at the FE nodes
        f_a
        %Loads applied in the local coordinate system at the midpoint of
        %the FE beam elements
        f_e
    end
    
    %Body geometry (public)
    properties
        %Initial strain in the undeformed state
        PreStrain
        %Initial curvature in the undeformed state
        PreCurvature
        %Offset distance from root of body to reference node
        BodyXYZOffset
        %Straight line distance along the beam nodes
        BeamAxisSegments
        %Length of each beam elements
        BeamElemLength
        %Positions along the body where the coordinates/orientations are
        %recovered
        BeamEvalPoints
        %Number of beam elements in this body
        NumBeamElem
        %Number of beam nodes in this body
        NumBeamNode
        %Number of evaluation points in this body
        NumEvalPoints
    end
    
    %Body geometry (internal)
    properties (Hidden, SetAccess = private)
        %Initial strain in the undeformed state
        Eps0
        %Initial curvature in the undeformed state
        Kappa0
        %Offset distance from root of body to reference node
        P0
        %Straight line distance along the beam nodes
        s_node
        %Length of each beam element
        ds_node
        %Positions along the body where the coordinates/orientations are
        %recovered
        s
        %Number of beam elements in this body
        NumElem
        %Number of beam nodes in this body
        NumNode
        %Number of evaluation points in this body
        NumPoint
    end
    
    %Transformation matrices
    properties
        %Transformation from the aircraft to the body frame
        BodyTransform
        %Transformation from the aircraft to the element frame
        ElementTransform
        %Transformation from one element to its next neighbour along the
        %beam
        Element2ElementTransform
    end
    
    %Transformation matrices (internal)
    properties (Hidden, SetAccess = private)
        %Transformation from the aircraft to the body frame
        CaB0
        %Transformation from the aircraft to the element frame
        CaBi
        %Transformation from one element to its next neighbour along the
        %beam
        CBB
    end
    
    %Indexing, shape functions, boundary conditions (internal)
    properties (SetAccess = private)
        %Structure describing the shape function terms
        ShapeFunctions
        %Structure describing the boundary conditions
        BoundaryConditions
        %Structure containing additional mass data (e.g. from point masses)
        MassData
        %Structure describing the element indices
        MatrixIndexing
    end
    
    %Helper properties
    properties (Constant, Hidden)
        %Direction of axial extension
        AxialStrainDirection = [1 ; 0 ; 0];
    end
    
    methods % set / get
        %   - System states
        function set.LocalBeamLoads(obj, val)    %set.LocalBeamLoads
            validateattributes(val, {'numeric'}, {'column', 'real', ...
                'finite', 'nonnan'}, class(obj), 'LocalBeamLoads');
            obj.x_f = val; %#ok<*MCSUP>
        end
        function val = get.LocalBeamLoads(obj)   %get.LocalBeamLoads
            val = obj.x_f;
        end
        %   - System matrices
        function set.MassMatrix(obj, val)              %set.MassMatrix
            validateattributes(val, {'numeric'}, {'3d', 'ncols', 6, ...
                'nrows', 6, 'real', 'finite', 'nonnan'}, class(obj), 'MassMatrix');
            obj.M = val;
        end
        function val = get.MassMatrix(obj)             %get.MassMatrix
            val = obj.M;
        end
        function set.StiffnessMatrix(obj, val)         %set.StiffnessMatrix
            validateattributes(val, {'numeric'}, {'3d', 'ncols', 6, ...
                'nrows', 6, 'real', 'finite', 'nonnan'}, class(obj), 'StiffnessMatrix');
            obj.K = val;
        end
        function val = get.StiffnessMatrix(obj)        %get.StiffnessMatrix
            val = obj.K;
        end
        function set.ComplianceMatrix(obj, val)        %set.ComplianceMatrix
            validateattributes(val, {'numeric'}, {'3d', 'ncols', 6, ...
                'nrows', 6, 'real', 'finite', 'nonnan'}, class(obj), 'ComplianceMatrix');
            obj.C = val;
        end
        function val = get.ComplianceMatrix(obj)       %get.ComplianceMatrix
            val = obj.C;
        end
        function set.ElementDirectionVector(obj, val)  %set.ElementDirectionVector
            validateattributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'ElementDirectionVector');
            obj.E1 = val;
        end
        function val = get.ElementDirectionVector(obj) %get.ElementDirectionVector
            val = obj.E1;
        end
        function set.LocalNodalLoads(obj, val)         %set.LocalNodalLoads
            validateattributes(val, {'numeric'}, {'column', 'nonnan', ...
                'finite', 'real'}, class(obj), 'LocalNodalLoads');
            obj.f_n = val;
        end
        function val = get.LocalNodalLoads(obj)        %get.LocalNodalLoads
            val = obj.f_n;
        end
        function set.AircraftNodalLoads(obj, val)      %set.AircraftNodalLoads
            validateattributes(val, {'numeric'}, {'column', 'nonnan', ...
                'finite', 'real'}, class(obj), 'AircraftNodalLoads');
            obj.f_a = val;
        end
        function val = get.AircraftNodalLoads(obj)     %get.AircraftNodalLoads
            val = obj.f_a;
        end
        function set.LocalElementLoads(obj, val)       %set.LocalElementLoads
            validateattributes(val, {'numeric'}, {'column', 'nonnan', ...
                'finite', 'real'}, class(obj), 'LocalElementLoads');
            obj.f_n = val;
        end
        function val = get.LocalElementLoads(obj)      %get.LocalElementLoads
            val = obj.f_e;
        end
        %   - Body geometry
        function set.PreStrain(obj, val)         %set.PreStrain
            validateatttributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'PreStrain');
            obj.Eps0 = val;
        end
        function val = get.PreStrain(obj)        %get.PreStrain
            val = obj.Eps0;
        end
        function set.PreCurvature(obj, val)      %set.PreCurvature
            validateatttributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'Curvature');
            obj.Eps0 = val;
        end
        function val = get.PreCurvature(obj)     %get.PreCurvature
            val = obj.Kappa0;
        end
        function set.BodyXYZOffset(obj, val)     %set.BodyXYZOffset
            validateattributes(val, {'numeric'}, {'column', 'nrows', 3, ...
                'real', 'nonnan', 'finite'}, class(obj), 'BodyXYZOffset');
            obj.P0 = val;
        end
        function val = get.BodyXYZOffset(obj)    %get.BodyXYZOffset
            val = obj.P0;
        end
        function set.BeamAxisSegments(obj, val)  %set.BeamAxisSegments
            validateattributes(val, {'numeric'}, {'row', 'increasing', ...
                'nonnan', 'finite', 'real'}, class(obj), 'BeamAxisSegments');
            obj.s_node = val;
        end
        function val = get.BeamAxisSegments(obj) %get.BeamAxisSegments
            val = obj.s_node;
        end
        function set.BeamElemLength(obj, val)    %set.BeamElemLength
            validateattributes(val, {'numeric'}, {'row', 'positive', ...
                'nonnan', 'finite', 'real'}, class(obj), 'BeamElemLength');
            obj.ds_node = val;
        end
        function val = get.BeamElemLength(obj)   %get.BeamElemLength
            val = obj.ds_node;
        end
        function set.BeamEvalPoints(obj, val)    %set.BeamEvalPoints
            validateattributes(val, {'numeric'}, {'row', 'increasing', ...
                'nonnan', 'finite', 'real'}, class(obj), 'BeamEvalPoints');
            obj.s = val;
        end
        function val = get.BeamEvalPoints(obj)   %get.BeamEvalPoints
            val = obj.s;
        end
        function set.NumBeamElem(obj, val)       %set.NumBeamElem
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'real', 'finite', 'nonnan'}, class(obj), 'NumBeamElem');
            obj.NumElem = val;
        end
        function val = get.NumBeamElem(obj)      %get.NumBeamElem
            val = obj.NumElem;
        end
        function set.NumBeamNode(obj, val)       %set.NumBeamNode
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'real', 'finite', 'nonnan'}, class(obj), 'NumBeamNode');
            obj.NumNode= val;
        end
        function val = get.NumBeamNode(obj)      %get.NumBeamNode
            val = obj.NumNode;
        end
        function set.NumEvalPoints(obj, val)     %set.NumEvalPoints
            validateattributes(val, {'numeric'}, {'scalar', 'positive', ...
                'real', 'finite', 'nonnan'}, class(obj), 'NumEvalPoints');
            obj.NumPoint = val;
        end
        function val = get.NumEvalPoints(obj)    %get.NumEvalPoints
            val = obj.NumPoint;
        end
        %   - Transformation matrices
        function set.BodyTransform(obj, val)             %set.BodyTransform
            valdiateattributes(val, {'numeric'}, {'2d', 'nrows', 3, ...
                'ncols', 3, 'nonnan', 'finite', 'real'}, class(obj), 'BodyTransform');
            obj.CaB0 = val;
        end
        function val = get.BodyTransform(obj)            %get.BodyTransform
            val = obj.CaB0;
        end
        function set.ElementTransform(obj, val)          %set.ElementTransform
            valdiateattributes(val, {'numeric'}, {'3d', 'nrows', 3, ...
                'ncols', 3, 'nonnan', 'finite', 'real'}, class(obj), 'ElementTransform');
            obj.CaBi = val;
        end
        function val = get.ElementTransform(obj)         %get.ElementTransform
            val = obj.CaBi;
        end
        function set.Element2ElementTransform(obj, val)  %set.Element2ElementTransform
            valdiateattributes(val, {'numeric'}, {'3d', 'nrows', 3, ...
                'ncols', 3, 'nonnan', 'finite', 'real'}, class(obj), 'Element2ElementTransform');
            obj.CaBB = val;
        end
        function val = get.Element2ElementTransform(obj) %get.Element2ElementTransform
            val = obj.CaBB;
        end
    end
    
    methods % constructor
        function obj = Model(FEM, varargin)
            %Model Constructor for the 'awi.methods.isfe.Model' class.
            %
            % Actions performed:
            %   1. If an instance of ''awi.fe.FEModel'' has been passed as
            %    the first argument then it will be passed to the ''setup''
            %    method.
            
            if nargin < 1
                FEM = [];
            end
            if isa(FEM, 'awi.fe.FEModel')
                setup(obj, FEM);
            end
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha)
            
            if nargin < 2 || isempty(ha)
                hF = figure('Name', 'Instrinsic, Strain-Based Model Geoemtry');
                ha = axes('Parent', hF, 'NextPlot', 'add', ...
                    'XLim'    , [-20, 30], ...
                    'YLim'    , [0, 85]  , ...
                    'ZLim'    , [-20, 60]  , ...
                    'View'    , [150, 35], ...
                    'Box'     , 'on'     , ...
                    'XGrid'   , 'on'     , ...
                    'YGrid'   , 'on'     , ...
                    'ZGrid'   , 'on');
                set([ha.XLabel], 'String', 'X [m]');
                set([ha.YLabel], 'String', 'Y [m]');
                set([ha.ZLabel], 'String', 'Z [m]');
            end
            
            %Get the current coordinates
            [coords, CaB] = strains2coords(obj);
            
            %Draw the nodes
            nodeCoords = coords(:, obj.MatrixIndexing.NodeInd);                        
            hg{1} = plot3(ha, nodeCoords(1, :), nodeCoords(2, :), nodeCoords(3, :), ....
                'Marker'         , 'o', ...                
                'MarkerFaceColor', 'g', ...
                'MarkerEdgeColor', 'k', ...
                'LineStyle'      , '-', ...
                'Color'          , 'k', ...
                'LineWidth'      , 2);
            
            %Draw the coordinate systems
            midElemCoords = coords(:, obj.MatrixIndexing.NodeMidInd);   
            nodeOrient    = CaB(:, :, obj.MatrixIndexing.NodeMidInd);
            hg{2} = i_drawCoordSys(ha, midElemCoords, nodeOrient);
            
            %Return column vector of hg-objects
            hg = vertcat(hg{:});
            
            function hg = i_drawCoordSys(ha, origin, rotMatrix)
                
                %Grab axis vectors
                eX = squeeze(rotMatrix(:, 1, :));
                eY = squeeze(rotMatrix(:, 2, :));
                eZ = squeeze(rotMatrix(:, 3, :));
                
                %Create coordinates
                OX = origin + eX;
                OY = origin + eY;
                OZ = origin + eZ;
                
                %Plot eX
                eXx = i_padCoordsWithNaN([origin(1, :) ; OX(1, :)]);
                eXy = i_padCoordsWithNaN([origin(2, :) ; OX(2, :)]);
                eXz = i_padCoordsWithNaN([origin(3, :) ; OX(3, :)]);
                hg(1) = line(ha, eXx, eXy, eXz, ...
                    'Color'    , 'r', ...
                    'LineWidth', 2  , ...
                    'Tag', 'Coordinate Systems');
                
                %Plot eY
                eYx = i_padCoordsWithNaN([origin(1, :) ; OY(1, :)]);
                eYy = i_padCoordsWithNaN([origin(2, :) ; OY(2, :)]);
                eYz = i_padCoordsWithNaN([origin(3, :) ; OY(3, :)]);
                hg(2) = line(ha, eYx, eYy, eYz, ...
                    'Color'    , 'g', ...
                    'LineWidth', 2  , ...
                    'Tag', 'Coordinate Systems');
                
                %Plot eZ
                eZx = i_padCoordsWithNaN([origin(1, :) ; OZ(1, :)]);
                eZy = i_padCoordsWithNaN([origin(2, :) ; OZ(2, :)]);
                eZz = i_padCoordsWithNaN([origin(3, :) ; OZ(3, :)]);
                hg(3) = line(ha, eZx, eZy, eZz, ...
                    'Color'    , 'b', ...
                    'LineWidth', 2  , ...
                    'Tag', 'Coordinate Systems');
                
                %Return a column vector of handles
                hg = hg(:);
                
            end
            
            function x = i_padCoordsWithNaN(x)
                %padCoordsWithNaN Accepts a matrix of [2, N] sets of
                %coordinates which represent the coordinate of a series of
                %lines from end-A to end-B and returns a single vector with all
                %of the coordinates padded by NaN terms.
                %
                % This function enables the plotting of line objects to be
                % vectorised.
                
                %Convert to cell so we retain the pairs of coordinates in the
                %correct order
                x  = num2cell(x, 1);
                
                %Preallocate
                x_ = cell(1, 2 * numel(x));
                
                %Assign the data and NaN terms
                x_(1 : 2 : end - 1) = x;
                x_(2 : 2 : end)     = {nan};
                
                %Return a column vector
                x = vertcat(x_{:});
                
            end
            
        end
    end
    
    methods % defining the Model
        function setup(obj, FEM, type, logfcn)
            %setup Populates the object properties using the data in the FE
            %model ''FEM''.
            %
            % Detailed Description:
            %   - The properties are dependent on the type of model that is
            %   requested, i.e. type = {'static', 'dynamic'}.
            %   - The status of the function will be output using the the
            %   function handle ''logfcn''.
            
            %Parse
            if nargin < 3
                type = 'static';
            end
            if nargin < 4
                logfcn = @(str) fprintf('%s\n', str);
            end
            p = inputParser;
            addRequired(p, 'FEM'   , @(x)validateattributes(x, ...
                {'awi.fe.FEModel'}, {'nonempty'}, 'FEM', 'setup'));
            addRequired(p, 'type'  , @(x)any(validatestring(x, ...
                {'static', 'dynamic'}, 'type', 'setup')));
            addRequired(p, 'logfcn', @(x) isa(x, 'function_handle'));
            parse(p, FEM, type, logfcn);
            if ~isempty(FEM.Children)
                error('Update code to allow converting AWI FEM into a mult-body system.');
            end
            assert(all(arrayfun(@(o) isa(o, 'awi.model.Beam'), ...
                [FEM.GeometryObject])), ['The Intrinsic, Strain-Based, ', ...
                'Finite Element formulation is only valid for models '  , ...
                'which are a collection of Beam objects.']);
            assert(numel(obj) == numel(FEM), sprintf(['Expected the ', ...
                'number of ''%s'' objects and ''%s'' objects to be ', ...
                'the same.'], class(obj), class(FEM)));
            
            pad = blanks(2);
            
            %Static terms
            setupStatic(obj, FEM, logfcn, pad);
            
            %Dynamic terms
            if strcmp(type, 'dynamic')
                setupDynamic(obj, FEM, logfcn, pad);
            end
            
        end
    end
    
    methods (Access = private) % defining the Model
        function setupStatic(obj, FEM, logfcn, pad)
            %setupStatic Populates the object properties which are required
            %for a static analysis.
            
            logfcn([pad, 'Parsing FE model...']);
            
            %   - Check the beams are defined in order along the beam
            BeamNodes = i_parseBeamNodes(FEM);
            
            function BeamNodes = i_parseBeamNodes(FEM)
                %i_parseBeamNodes Ensures that the 'awi.fe.Beam' objects
                %are defined in order along the span of the parent beam.
                %That is, end-B of the i^th beam is end-A of the (i+1)^th
                %beam.
                %
                % Returns only the unique nodes connecting the beams.
                
                BeamNodes = FEM.BeamNodes; %[endA ; endB]
                
                temp = [BeamNodes(1, 2 : end) ; BeamNodes(2, 1 : end - 1)];
                idx  = arrayfun(@(i) isequal(temp(1, i), temp(2, i)), 1 : size(temp, 2));
                assert(all(idx), ['Expected the ''awi.fe.Beam'' objects to be ', ...
                    'defined in one continuous line along the beam.']);
                
                BeamNodes = [BeamNodes(1, :), BeamNodes(2, end)];
            end
            
            logfcn(sprintf('%sDefining properties for body ''%s''...', pad, FEM.Name));
            
            %% Define reference node
            
            logfcn([pad, 'Defining reference node...']);
            
            if isempty(FEM.Children)
                RefNode = BeamNodes(1);
            else
                error('Update case for when the model has multiple bodies.');
            end
            
            %% Structural and Aero Nodes
            %   - A beam element will defined for each Beam object in the FEM.
            %   - This block defines position of the nodes in the 's-domain', which is
            %     the straight line distance along the FE beams. It also defines the
            %     mid-point positions which is where the local strains and curvatures
            %     are defined.
            %   - If aerodynamic data has been provided the application points for the
            %     aerodynamic loads will also be defined. TODO - Actually do this bit!
            
            logfcn([pad, 'Setting up structural nodes and elements...']);
            
            %Grab the node coordinates and set up segment-wise geometry
            nodeCoords = [BeamNodes.X];
            %   - Straight line distance along the nodes (hypotenuse)
            s_node    = awi.model.Stick.getLineLength(nodeCoords(1, :), nodeCoords(2, :), nodeCoords(3, :)); %#ok<*PROPLC>
            %   - Mid point of the elements
            s_node_mp = (s_node(1 : end - 1) + s_node(2 : end)) ./ 2;
            %   - Length of each element
            ds_node   = diff(s_node);
            
            %Aero?
            if isempty(FEM.AeroPanels)
                s_aero    = nan;
                s_aero_mp = nan;
            else
                %TODO - Implement generic aero codes
                s_aero    = nan;
                s_aero_mp = nan;
                %error('Update method to account for aerodynamic panels');
            end
            
            %Unique points
            s           = unique([s_node, s_node_mp, s_aero, s_aero_mp]);
            s(isnan(s)) = [];
            nNode       = numel(s_node);
            nElem       = numel(BeamNodes) - 1;
            nPoint      = numel(s);
            
            %Indexing so we can easily recover these points later
            indNode         = find(any((s_node'    - s) == 0, 1));
            indNodeMidPoint = find(any((s_node_mp' - s) == 0, 1));
            indAero         = find(any((s_aero'    - s) == 0, 1));
            indAeroMidPoint = find(any((s_aero_mp' - s) == 0, 1));
            
            obj.MatrixIndexing = struct('NodeInd', indNode, ...
                'NodeMidInd', indNodeMidPoint, 'AeroInd', indAero, ...
                'AeroMidInd', indAeroMidPoint);
            
            %Offset of this beam from the reference node?
            p0 = nodeCoords(:, 1) - RefNode.X;
            
            %% Beam rotation matrices
            %   - Defines the following properties for each body/beam element:
            %       * 'CaB0' - The body orientation defined as a 3x3 rotation matrix.
            %       * 'CBBi' - The elem orientation defined as a 3x3 rotation matrix.
            %       * 'CBB'  - The rotation matrix from one beam element to another.
            %   - N.B. For a rotation matrix inv(R) = transpose(R)
            
            logfcn([pad, 'Defining element and body rotation matrices...']);
            
            %Define element orientation
            %   - Assume local x-axis is along the beam axis
            eX = nodeCoords(:, 2 : end) - nodeCoords(:, 1 : end - 1);
            %   - Orientation of local xy-plane is taken from the 'Beam' object
            %     orienation vector, which itself is defined in the output coordinate
            %     system of the node at end-A.
            eY_tmp_local  = [FEM.Beams.X];
            rotMatrix     = getRotationMatrix([BeamNodes(1 : end - 1).OutputCoordSys]);
            eY_tmp_global = arrayfun(@(ii) rotMatrix(:, :, ii) * eY_tmp_local(:, ii), 1 : nElem, 'Unif', false);
            eY_tmp_global = horzcat(eY_tmp_global{:});
            %   - Construct rotation matrix
            CaBi = generateRotationMatrix(eX, eY_tmp_global);
            
            %Orienation at the root of the body
            CaB0 = CaBi(:, :, 1);
            
            %Rotation from one beam element to the next
            %   - CBB = inv(R_ii-1) * inv(R_ii)
            %   - N.B. First beam is always aligned with body orienation, hence eye(3).
            CBB = arrayfun(@(ii) (CaBi(:, :, ii - 1)' * CaBi(:, :, ii))', 2 : nElem, 'Unif', false);
            CBB = cat(3, eye(3), CBB{:});
            
            %% Beam stiffness/compliance matrix
            %   - This formulation does not allow tapered beam sections, therefore the
            %     A/I/J values for each beam element are the average of the values at
            %     end-A and end-B.
            
            logfcn([pad, 'Grabbing beam sectional properties...']);
            
            %Cross-section & mass properties of the beam
            SectionalProps = getBeamSectionalProps(FEM);
            
            function SectionalProps = getBeamSectionalProps(FEM)
                %getBeamSectionalProps Returns a MATLAB structure describing the
                %A/I11/I22/I12/J values along the beam, the radius of gyration
                %terms and the E/G/rho distribution from the material.
                %
                % N.B. This formulation does not allow tapered beam sections,
                % therefore the A/I/J values for each beam element are the average
                % of the values at end-A and end-B.
                
                %Grab the 'awi.fe' objects describing the beam props and materials
                BeamProps    = FEM.BeamProps;
                BeamMaterial = [BeamProps.Material];
                
                %Material
                E   = [BeamMaterial.E];
                G   = [BeamMaterial.G];
                rho = [BeamMaterial.Rho];
                
                %Sectional terms
                A   = mean([BeamProps.A]  , 1);
                I11 = mean([BeamProps.I11], 1);
                I22 = mean([BeamProps.I22], 1);
                I12 = mean([BeamProps.I12], 1);
                J   = mean([BeamProps.J]  , 1);
                
                %Radius of gyration
                h = sqrt(sqrt(I11 ./ I22) ./ A);
                b = A ./ h;
                
                SectionalProps = struct('A', A, 'I11', I11, 'I22', I22, ...
                    'I12', I12, 'J', J, 'h', h, 'b', b, 'E', E, 'G', G, 'Rho', rho);
                
            end
            
            logfcn([pad, 'Setting up beam stiffness matrix...']);
            
            %Stiffness matrix without including shear offsets
            stiffMatrix = defineBeamStiffnessMatrix(SectionalProps, nElem);
            
            function K = defineBeamStiffnessMatrix(SectionalProps, nBeams)
                %defineBeamStiffnessMatrix Returns a [6, 6, nBeams] matrix
                %describing the stiffness properties of the beams.
                %
                % N.B. This formulation does not allow tapered beam sections,
                % therefore the A/I/J values for each beam element are the average
                % of the values at end-A and end-B.
                %
                % N.B. Vector of strains and curvatures {sc} is ordered as follows:
                %   [eps_x ; eps_y ; eps_z ; kappa_x ; kappa_y ; kappa_z]
                %   therefore the local forces and moments are given by [K] * {sc}
                
                K = zeros(6, 6, nBeams);
                
                K(1, 1, :) = SectionalProps.E .* SectionalProps.A;    %axial stiffness
                K(2, 2, :) = SectionalProps.G .* SectionalProps.A;    %timoshenko terms (???)
                K(3, 3, :) = K(2, 2, :);
                K(4, 4, :) = SectionalProps.G .* SectionalProps.J;    %torsional stiffness
                K(5, 5, :) = SectionalProps.E .* SectionalProps.I22;  %bending stiffness (plane 2)
                K(6, 6, :) = SectionalProps.E .* SectionalProps.I11;  %bending stiffness (plane 1)                
                K(5, 6, :) = SectionalProps.E .* SectionalProps.I12;  %bending stiffness (cross-plane)
                K(6, 5, :) = K(5, 6, :);
                
            end
            
            logfcn([pad, 'Accounting for shear offsets...']);
            
            %Stiffness matrix with effect of shear offsets
            [stiffMatrix, CFFbar] = adjustStiffnessForShearOffsets(FEM, stiffMatrix, nElem, nodeCoords, ds_node, CaBi, rotMatrix, eY_tmp_global);
            
            function [K, CFFbar] = adjustStiffnessForShearOffsets(FEM, K, nBeams, nodeCoords, ds, CaBi, rotMatrix, eY_tmp_global)
                %adjustStiffnessForShearOffsets Accounts for the effects of offset
                %shear centres in the stiffness matrices.
                %
                % Example of shear centre offsets when looking down on the beam...
                %
                %        x <-- Shear offset at end-A
                %        |
                %        |                           x <---- Shear offset at end-B
                %        |                           |
                % (Node) O-+-+-+-+-+-+-+-+-+-+-+-+-+-O (Node)
                %        A             ^             B
                %                      |
                %                      |
                %                   Beam axis
                
                %Shear offsets
                %   - Local frame ('awi.fe.BeamProp' has no offset in local x, hence 0.)
                SCx = zeros(2, nBeams);
                SCy = [FEM.Beams.SCy];
                SCz = [FEM.Beams.SCz];
                drSC_A = [SCx(1, :) ; SCy(1, :) ; SCz(1, :)];
                drSC_B = [SCx(2, :) ; SCy(2, :) ; SCz(2, :)];
                %   - Global frame
                drSC_A = arrayfun(@(ii) rotMatrix(:, :, ii) * drSC_A(:, ii), 1 : nBeams, 'Unif', false);
                drSC_A = horzcat(drSC_A{:});
                drSC_B = arrayfun(@(ii) rotMatrix(:, :, ii) * drSC_B(:, ii), 1 : nBeams, 'Unif', false);
                drSC_B = horzcat(drSC_B{:});
                %   - Average of ends A and B
                drSC_bar = (drSC_A + drSC_B) ./ 2;
                
                %Coordinates of shear centres at end-A & end-B in the global frame
                rSC_A = nodeCoords(:, 1 : end - 1) + drSC_A;
                rSC_B = nodeCoords(:, 2 : end)     + drSC_B;
                
                %Construct rotation matrix for the shear axis
                eX_    = rSC_B - rSC_A;
                CB_SC  = generateRotationMatrix(eX_, eY_tmp_global);
                
                %Construct rotation matrix from the element to the shear coordinate system
                CBBbar = arrayfun(@(ii) CaBi(:, :, ii)' * CB_SC(:, :, ii), 1 : nBeams, 'Unif', false);
                %CBBbar = cat(3, CBBbar{:});
                
                %Construct 'R_tilde'
                %   - 'R_tilde' = [1, skew(CaB * rSC_bar) ; 0, 1] (i.e. 6 x 6 x nBeam)
                %   - 'R_tilde' accounts for axial load induced by shear centre offset when
                %     beam element bends (cross product of shear centre moment arm).
                temp    = arrayfun(@(ii) CaBi(:, :, ii)' * drSC_bar(:, ii), 1 : nBeams, 'Unif', false);
                temp    = horzcat(temp{:});
                temp    = skewmat(temp);
                eyz     = repmat(eye(3), [1, 1, nBeams]);
                zeroz   = zeros(3, 3, nBeams);
                R_tilde       = [eyz, temp ; zeroz, eyz];
                R_tilde_trans = [eyz, zeroz ; temp, eyz];
                
                dsprime     = vecnorm(eX_);
                length_fact = ds ./ dsprime;
                
                %Construct block diagonal matrices
                CBBbar_blk       = cellfun(@(x) blkdiag(x, x), CBBbar, 'Unif', false);
                CBBbar_trans_blk = cellfun(@(x) blkdiag(x', x'), CBBbar, 'Unif', false);
                
                %Pre and post multiply stiffness matrix by rotation matrices in
                %order to evaluate loads about the beam axis and not the shear axis
                %(rotate and shift loads)
                %   K = R_t' * CBBbar * K * CBBbar' * R_t * ds
                K = arrayfun(@(ii) ...
                    R_tilde_trans(:, :, ii) * CBBbar_blk{ii} * K(:, :, ii) *    ...
                    CBBbar_trans_blk{ii} * R_tilde(:, :, ii) * length_fact(ii), ...
                    1 : nBeams, 'Unif', false);
                
                %Calculate rotation matrix to transform loads from beam axis to shear axis
                CFFbar = arrayfun(@(ii) inv(R_tilde_trans(:, :, ii) * CBBbar_blk{ii}), ...
                    1 : nBeams, 'Unif', false);
                
                %Stack as 3D arrays
                K      = cat(3, K{:});
                CFFbar = cat(3, CFFbar{:});
                
            end
            
            logfcn([pad, 'Inverting stiffness matrix to obtain compliance matrix...']);
            
            %Get full-sized compliance matrix and shear-axis-to-beam-axis
            %matrix
            [compMatrix, compMatrix_full, CFFbar_full] = getComplianceMatrix(stiffMatrix, CFFbar, nElem);
            
            function [compMatrix, compMatrix_full, CFFbar_full] = getComplianceMatrix(K, CFFbar, nElem)
                
                %Compliance matrix
                compMatrix = arrayfun(@(ii) inv(K(:, :, ii)), 1 : nElem, 'Unif', false);
                
                %Block diagonal of compliance matrices (no shape functions)
                %   - Multiply the force & moment vector by 'Cfull' to obtain the
                %     curvatures & strains                
                compMatrix_full = blkdiag(compMatrix{:});                
                compMatrix      = cat(3, compMatrix{:});
                
                %Block diagonal of shear stiffness matrices (no shape functions)
                CFFdiag     = mat2cell(CFFbar, 6, 6, ones(1, nElem));
                CFFbar_full = blkdiag(CFFdiag{:});
                
            end
            
            %Pre-strain and pre-curvature
            %   - TODO: Add this in to the beam property definition somehow...
            eps0   = zeros(3, nElem);
            kappa0 = zeros(3, nElem);
            
            %% Beam mass matrix
            %   - TODO: Add point mass terms
            %   - TODO: Get reference for radius of gyration calculation
            
            logfcn([pad, 'Calculating contributions from point masses...UPDATE THIS!!!!']);
            
            %Calculate contributions from point mass data - UPDATE THIS!!!
            [PointMass, PointInertia, PointCoG] = getPointMassContributions(FEM);
            
            M_pt_store = zeros(6, 6, nElem);
            
            function [PointMass, PointInertia, PointCoG] = getPointMassContributions(FEM)
                %getPointMassContributions Returns the mass, inertia and CoG terms
                %related to point masses in the model.
                
                PointMass    = [];
                PointInertia = [];
                PointCoG     = [];
                
            end
            
            logfcn([pad, 'Setting up beam mass matrix...']);
            
            %Define the mass matrix (including point mass terms)
            massMatrix = defineBeamMassMatrix(SectionalProps, PointMass, PointInertia, PointCoG, ds_node, nElem);
            
            function massMatrix = defineBeamMassMatrix(SectionalProps, PointMass, PointInertia, PointCoG, ds, nBeams)
                %defineMassMatrix Defines the mass matrix for the beam using
                %contributions from the beam sectional properties as well as the
                %point mass data.
                
                %Point mass data
                m_pt0  = 0;
                J_pt0  = zeros(3, 3, nBeams);
                cg_pt0 = zeros(3, nBeams);
                
                %Sectional properties
                A   = SectionalProps.A;
                b   = SectionalProps.b;
                h   = SectionalProps.h;
                rho = SectionalProps.Rho;
                
                J = zeros(3, 3, nBeams);
                M = zeros(3, 3, nBeams);
                
                %Mass terms
                m   = rho .* A + m_pt0 ./ ds; %Smear point mass over the length of the element
                i1  = rho .* A ./ 12 .* (b.^2 + h.^2) ;
                i2  = rho .* A .* (b.^2 + ds.^2) ./ 12;%1e-5;%
                i3  = rho .* A .* (h.^2 + ds.^2) ./ 12;%1e-5;%
                %J   = J + (J_pt0 ./ ds);
                cg  = zeros(3, nBeams) + cg_pt0;
                
                %Matrix form
                M(1, 1, :) = m;
                M(2, 2, :) = m;
                M(3, 3, :) = m;
                J(1, 1, :) = i1;
                J(2, 2, :) = i2;
                J(3, 3, :) = i3;
                
                %Get skew of cg
                skewCG = skewmat(cg);
                
                %Form the mass matrix
                %   [M] = [m, -m x cg ; m x cg, J]
                massMatrix = arrayfun(@(ii) [ ...
                    M(:, :, ii)             , -m(ii) * skewCG(:, :, ii) ; ...
                    m(ii) * skewCG(:, :, ii), J(:, :, ii)], 1 : nBeams, 'Unif', false);
                massMatrix = cat(3, massMatrix{:});
                
            end
            
            %Stash additional mass information for use downstream
            obj.MassData = struct('PointMassTerms', M_pt_store);
            
            %% Applied loads
            %   - Defines the load vector for loads applied at the nodes and midpoints
            %     of the beam elements.
            %   - There are 6 loads for each node, defined in the following order:
            %       [F_x ; F_y ; F_z ; M_x ; M_y ; M_z]
            
            %Initialise as zeros
            f        = zeros(nNode * 6, 1); %follower loads applied at nodes
            f_a      = zeros(nNode * 6, 1); %global loads applied at nodes
            f_mp_pnt = zeros(nElem * 6, 1); %follower loads applied at element mid-point
            
            %Assign load data
            AppliedLoads = FEM.PointLoads;
            
            if ~isempty(AppliedLoads)
                
                %Account for follower and non-follower loads
                idxFollower  = ismember(get(AppliedLoads, {'LoadBehaviour'}), 'follower');
                
                %Follower loads
                f = populateLoadVector(AppliedLoads(idxFollower), BeamNodes, f);
                
                %Global loads
                f_a = populateLoadVector(AppliedLoads(~idxFollower), BeamNodes, f_a);                
                
            end
            
            function loadVec = populateLoadVector(PointLoads, BeamNodes, loadVec)
                %popualteLoadVector Assigns the applied loads to the correct
                %position in the load vector.
                
                if isempty(PointLoads)
                    return
                end
                
                %Get load properties
                loadType   = {PointLoads.LoadType};
                loadOrient = [PointLoads.Orientation];
                loadMag    = [PointLoads.Magnitude];
                
                %Get index number of the nodes
                %   - TODO : Possibly overload 'isequal' for awi.fe object in order to
                %   allow [Node1 ; Node2] == [nan, Node1, nan, nan, nan, Node2] type
                %   operations
                ind      = arrayfun(@(n) find(BeamNodes == n), [PointLoads.Node]);
                loadInd  = cell(1, numel(PointLoads));
                loadInd(ismember(loadType, 'force'))  = {[1 ; 2 ; 3]};
                loadInd(ismember(loadType, 'moment')) = {[4 ; 5 ; 6]};
                loadInd = horzcat(loadInd{:});
                loadInd = loadInd + ((ind - 1) * 6);
                
                %Force moment in (x, y, z) directions
                load =  loadOrient .* loadMag;
                
                %Assign
                loadVec(loadInd(:)) = load(:);
                
            end
            
            %% Define matrix terms
            
            logfcn([pad, 'Defining finite element terms...']);
            
            %Constant terms
            e1_   = [ones(1, nElem) ; zeros(5, nElem)] + [eps0 ; kappa0];
            %   - Sets boundary conditions by eliminating certain rows and columns from
            %   the system matrices. This formulation will define a cantilever beam
            %   with a fixed root.
            R1_bc   = [zeros(nElem * 6, 6), eye(nElem * 6)];
            
            %Preallocate
            e1      = zeros(6, nElem);
            Ab      = zeros(12, 6, nElem);
            Btot_mp = zeros(nNode * 6, nElem * 6);
            Dtot    = zeros(nNode * 6, nElem * 6);
            
            %Easier to do this in a loop otherwise we need numerous calls to arrayfun
            for ii = 1 : nElem
                
                indRow = (1 : 12) + (ii - 1) * 6;
                indCol = (1 : 6)  + (ii - 1) * 6;
                
                %Define vector 'e1'
                %   - Defines a direction vector for integrating along the beam element
                %   - See Eqn 2. of [1] or Eqn. 8 of [2].
                e1(:, ii) = stiffMatrix(:, :, ii) * e1_(:, ii);
                
                %Account for change in orientation from one beam element to the next
                CBB_tot = blkdiag(CBB(:, :, ii), CBB(:, :, ii), eye(6));
                
                %Shape functions for finite elements
                %   - See Eqn. 18 of [2].
                A1 = [-eye(6) * ds_node(ii), zeros(6)] / (-ds_node(ii) * CBB_tot);
                A2 = [ eye(6)              , -eye(6) ] / (-ds_node(ii) * CBB_tot);
                
                %Integral of shape function terms
                Ab(:, :, ii)  = A1' * ds_node(ii) + A2' * ds_node(ii)^2/2;
                
                %Differentiation matrix (for calculating tangent matrix)
                D = -A2' * eye(6) * ds_node(ii);
                
                %Collate terms
                Btot_mp(indRow, indCol) = Btot_mp(indRow, indCol) + Ab(:, :, ii);
                Dtot(indRow, indCol)    = Dtot(indRow, indCol) + D;
                
            end
            
            %Collect
            obj.ShapeFunctions = struct('Spatial', Ab, ...
                'Derivative', Dtot, 'Integral', Btot_mp);
            obj.BoundaryConditions = struct('R1_bc', R1_bc);
            
            logfcn([pad, 'Finite element matrices populated!']);
            
            %% Assign to object
            
            %System matrices
            obj.M      = massMatrix;
            obj.K      = stiffMatrix;
            obj.C      = compMatrix;
            obj.C_full = compMatrix_full;
            obj.E1     = e1;
            
            %Applied loads
            obj.f_n = f;
            obj.f_a = f_a;
            obj.f_e = f_mp_pnt;
            
            %Body geometry
            obj.P0       = p0;
            obj.s_node   = s_node;
            obj.ds_node  = ds_node;
            obj.s        = s;
            obj.Eps0     = eps0;
            obj.Kappa0   = kappa0;
            obj.NumNode  = nNode;
            obj.NumElem  = nElem;
            obj.NumPoint = nPoint;
            
            %Orientation
            obj.CaB0        = CaB0;
            obj.CaBi        = CaBi;
            obj.CBB         = CBB;
            obj.CFFbar_full = CFFbar_full;
            
            logfcn([pad, 'Body definition complete!']);
                     
            drawElement(obj);
            
        end
        function setupDynamic(obj, FEM, type, logfcn, pad)
            
            error('Update code');
            
        end
    end
    
    methods % setting up the system states
        function initialiseStates(obj, Options)
            %initialiseStates Updates the state values with sensible
            %defaults depending on the desired analysis type.
            
            %Loop through object arrays
            if numel(obj) > 1
                arrayfun(@(o) initialiseStates(o, Options), obj);
                return
            end
            
            %Define the states
            switch Options.AnalysisType
                
                case 'static'
                    
                    %Set initial beam forces & moments to zero
                    if isempty(obj.x_f)
                        obj.x_f = zeros(6 * obj.NumElem, 1);
                    end
                    
                otherwise
                    
                    error('Update code');
                    
            end
            
        end
    end
    
    methods % analytical methods
        function [coords, CaB] = strains2coords(obj, x_f)
            %strains2coords Converts beam states into coordinates and
            %orientations in the body frame.
            
            if nargin < 2
                x_f = obj.x_f;
            end
            if isempty(x_f)
                x_f = zeros(obj.NumElem * 6, 1);
            end
            
            %Loop through object arrays
            if numel(obj) > 1
                temp   = arrayfun(@(o) strains2coords(o), obj, 'Unif', false);
                coords = cellfun(@(x) x{1}, temp, 'Unif', false);
                CaB    = cellfun(@(x) x{2}, temp, 'Unif', false);
                return
            end
            
            %Preallocate
            coords = zeros(3, obj.NumPoint);
            CaB    = zeros(3, 3, obj.NumPoint);
                        
            %Initial orientation of body w.r.t aircraft
            CaB(:, :, 1) = obj.CaB0;
            
            %Indexing
            counter  = 1;
            counter2 = 1;
            ub = cumsum(repmat(6, [1, obj.NumElem]));
            lb = [1, ub(1 : end - 1) + 1];
            
            %Multiply forces & moments by compliance matrix to obtain strains and
            %curvatures
            % -> straincurv = [eps_x ; eps_y ; eps_z ; kappa_x ; kappa_y ; kappa_z]
            straincurv = obj.C_full * x_f;
            e1         = obj.AxialStrainDirection;
            
            for ii = 1 : obj.NumElem
                
                %Grab strains and curvatures for this element with pre-curve
                eps_i   = straincurv(lb(ii) : ub(ii) - 3)      + obj.Eps0(:,ii)   + e1;
                kap_i   = skew(straincurv(ub(ii) - 2 : ub(ii)) + obj.Kappa0(:,ii))    ;
                
                %Which nodes lie in the element
                idx = and(obj.s > obj.s_node(ii), obj.s <= obj.s_node(ii + 1));
                
                %Calculate distance along each element for evaluation points in this element
                ds_eval = (obj.s(idx) - obj.s_node(ii));
                
                %Calculate new coordinates and orientations for each node
                %   - N.B. Strains and curvatures are constant across the element
                for kk = 1 : nnz(idx)                    
                    %Exponential term
                    expon   = [zeros(1, 4); eps_i, kap_i] * ds_eval(kk);
                    
                    %         expon_tot = expon_tot + expon;
                    coord_CaB = expm(expon') * [ ...
                        coords(:,counter2)'; ...
                        obj.CBB(:, :, ii) * CaB(:, :, counter2)'];
                    CaB(:, :, counter+1) = coord_CaB(2 : end, :)';
                    coords(:, counter+1) = coord_CaB(1, :)';
                    
                    counter = counter + 1;
                end
                counter2 = counter;
            end
            
        end        
    end
    
end

