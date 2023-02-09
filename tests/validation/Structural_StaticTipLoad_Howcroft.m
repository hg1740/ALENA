function Structural_StaticTipLoad_Howcroft
%Structural_StaticTipLoad_Howcroft Large deformation static analysis of a
%pre-curved beam for validating the structural solvers in the AWI software.
%
% Detailed description:
%   * The test case is documented in the paper "On the Geometrically Exact
%   Low Order Modelling of a Flexible Beam: Formulation and Numerical
%   Tests", Howcroft et. al, Jan 2018. The section related to this test
%   case is titled "Large Static Deformation"
%       - Geometry: The beam is a pre-curved beam which forms a shape of
%         1/8 of a circle with radius 100m.
%       - Stiffness: The beam has a constant stiffness distribution of
%         EI_xx = EI_zz = 25/3 Nm^2 and GJ = 7.02885 x 10^5 Nm/rad.
%       - Loading: A 600N load is applied to the tip node - both follower
%         and global loads are considered.

close all

%% Parameters

%Beam parameters
E         = 1e7;        %[N/m^2] , Young's modulus - required because the ALENA beam object accepts I/A/J and E/G values, not EI/EA/GJ.
A         = 1;          %[m]     , cross-sectional area
EIxx      = 25/3 * 1e7;       %[Nm^2]  , in-plane bending stiffness
EIzz      = EIxx;       %[Nm^2]  , out-of-plane bending stiffness
GJ        = 7.20885e5;  %[Nm/rad], torsional stiffness
R         = 100;        %[m]     , radius of circle
theta_bar = 1/8;        %[-]     , fraction of the circle represented by the beam

%Applied loads
loadMag   = 600;        %[N]     , magnitude of applied load

%Structural discretization
nStructuralNodes = 25;  %[-]     , number of structural nodes representing the curve

%Dependent parameters
kappa0 = 1 / R;         %[1/m]   , curvature
G      = E;             %[N/m^2] , shear modulus
I11    = EIxx / E;      %[m^4]   , 2nd moment of area in plane 1
I22    = EIzz / E;      %[m^4]   , 2nd moment of area in plane 2
J      = GJ   / G;      %[m^4]   , polar moment of area

%% Set up graphics objects

%Geometry drawing
hF = figure('Name', 'Geometry Model');
hAx = axes('Parent', hF);

%Analysis Model drawing
hF(2)  = figure('Name', 'Analysis Model');
hAx(2) = axes('Parent', hF(2));

%Set axes appearance
set(hAx, ...
    'NextPlot', 'add'    , ...
    'XLim'    , [-20, 30], ...
    'YLim'    , [0, 85]  , ...
    'ZLim'    , [-20, 60]  , ...
    'View'    , [150, 35], ...
    'Box'     , 'on'     , ...
    'XGrid'   , 'on'     , ...
    'YGrid'   , 'on'     , ...
    'ZGrid'   , 'on');
set([hAx.XLabel], 'String', 'X [m]');
set([hAx.YLabel], 'String', 'Y [m]');
set([hAx.ZLabel], 'String', 'Z [m]');

%% Calculate beam geometry
%   - Beam forms an arc 1/8 of a circle in the lower-right quadrant of the
%     (x,y) cartesian plane.

%x = r x cos(theta), y = r x sin(theta)
thetaVec = linspace(1.5 * pi, 1.5 * pi + theta_bar * 2 * pi, nStructuralNodes);
yd       = R .* cos(thetaVec);
xd       = R .* sin(thetaVec);
zd       = zeros(1, nStructuralNodes);

%Set origin to (0, 0, 0)
xd = xd - xd(1);
yd = yd - yd(1);
zd = zd - zd(1);

%Plot to check
plot3(hAx(1), xd, yd, zd, '-');

%% Generate the ALENA objects
%   - 1 x awi.model.Beam             : For describing the geometry and stiffness properties
%   - 1 x awi.model.Material         : For describing the material properties
%   - nNode x awi.model.CoordSys     : For describing the beam orientation
%   - 1 x awi.model.LoadDistribution : For describing the loading conditions
% Beam   = awi.model.Beam;
Beam = awi.model.LiftingSurface;
Mat    = awi.model.Material;
Orient = arrayfun(@(~) awi.model.CoordSys, 1 : nStructuralNodes);
Load   = awi.loads.LoadDistribution;

%Beam geometry
set(Beam, 'XData', xd, 'YData', yd, 'ZData', zd, 'Origin', [0, 0, 0]);

%LiftingSurface data
Beam.SpanVector = 'Y';
% Beam.Chord = [2, 2];
Beam.LE = [xd - 1 ; yd ; zd];
Beam.TE = [xd + 1 ; yd ; zd];
Beam.ActiveSet = 'cSet';
build(Beam);
drawElement(Beam, hAx(1));

%Set up coordinate systems
%   - Position along the beam
r     = Beam.RData;
r_bar = r ./ r(end);
set(Orient, {'SOffset'}, num2cell(r_bar)');
%   - Orientation
%       * Local x along beam
%       * Local y in same plane as global X axis
v1 = [diff(xd) ; diff(yd) ; diff(zd)];
v1 = num2cell([v1 , v1(:, end)]', 2);
v2 = [-1, 0, 0];
set(Orient, {'V1'}, v1, 'V2', v2, 'GridVector', 'x', 'GridPlane', 'xy', 'ActiveSet', 'vSet');
%   - Build the object and assign to the beam
arrayfun(@build, Orient);
assignBeamObject(Beam, Orient, r_bar);

%Define stiffness properties (constant)
set(Beam, ...
    'A'  , [A  , A  ], 'A_eta'  , [0, 1], ...
    'I11', [I11, I11], 'I11_eta', [0, 1], ...
    'I22', [I22, I22], 'I22_eta', [0, 1], ...
    'J'  , [J  , J  ], 'J_eta'  , [0, 1]);

%Set material properties and add to the object
set(Mat, 'E', E, 'G', G);
assignBeamObject(Beam, [Mat, Mat], [0, 1]);

%Applied loads
Load.EtaDistribution     = [0.5, 0.7, 0.9, 1];                               %tip load
Load.PointLoads          = [0 ; 0 ; loadMag ; 0 ; 0 ; 0];   %positive Z direction
Load.PointLoads          = repmat(Load.PointLoads, [1, numel(Load.EtaDistribution)]);
Load.PointLoadCoordSys   = 'local';                         %orientation is in local system
Load.PointLoadBehaviour  = 'follower';                      %follower force
Load.AxisFlag            = 'R';                             %defined along beam piecewise length
Load.MaximumVectorLength = 5;                               %scale the plot so we can clearly see it
%   - Attach the load to the beam object
assignLoadToBeam(Load, Beam);

%Draw to check
drawElement(Beam, hAx(1));

%% Convert to analysis model

%Ensure consistent discretization of the geometry
Beam.NumBeamElem = (nStructuralNodes - 1);

%Make the model
FEM = convertToFE(Beam);

%Change visualisation of load
set([FEM.PointLoads], 'MaxVectorLength', Load.MaximumVectorLength);

%Draw to check
%   - Overlay the nodes to make sure FE geometry matches the ALENA geometry
draw(FEM, hAx(1), 'PartList', {'Nodes'});
draw(FEM, hAx(2));

%% Run analysis

%Set up analysis options
StaticOptions = awi.methods.Options;
StaticOptions.StructuralResponse = 'nonlinear';

resDir = awi.methods.Nastran.makeDefaultAnalysisDirectory('nonlinear_static', pwd);

%Nastran analysis
Nastran = awi.methods.Nastran;
Nastran.AnalysisModel = FEM;
NasResults = static(Nastran, StaticOptions, resDir);

% runOriginalRobbie(FEM, awi.methods.Options);

%Intrinsic, strain-based, FE
StrainFE = awi.methods.IntrinsicStrainFE;
StrainFE.AnalysisModel = FEM;
StrainResults = static(StrainFE, StaticOptions);


end
