function Structural_TipLoad_Weihau_Su_Thesis
%Structural_TipLoad_Weihau_Su_Thesis Large deformation static analysis of a
%straight beam for validating the structural solvers in the AWI software.
%
% Detailed description:
%   * The test case is documented in the thesis "Coupled Nonlinear 
%     Aeroelasticity and Flight Dynamics of Fully Flexible Aircraft", Su, 
%     2008. The section related to this test case is "Section 4.1 - 
%     Cantilever Beam Configuration" 
%       - Geometry: The beam is straight with length 1m and cantilever
%         boundary conditions.
%       - Stiffness: The beam has a constant stiffness distribution of
%           * EA    = 1    x 10^6 N
%           * EI_xx = 5    x 10^1 Nm^2
%           * EI_zz = 1.25 x 10^3 Nm^2
%           * GJ    = 8    x 10^1 Nm^2
%       - Mass: The beam has a constant mass distribution of 
%           * NSM = 0.1 kg/m
%           * Ixx = 1.3  x 10^-4 kgm
%           * Iyy = 5    x 10^-6 kgm
%           * Izz = 1.25 x 10^-4 kgm
%       - Loading: A series of global tip loads in the z-direction are
%         applied up to 150N.

close all

%% Parameters

%Beam parameters
length    = 1;          %[m]     , length of the beam
E         = 1e7;        %[N/m^2] , Young's modulus - required because the ALENA beam object accepts I/A/J and E/G values, not EI/EA/GJ.
EA        = 1e6;        %[N]     , axial stiffness
EIxx      = 5e1;        %[Nm^2]  , flat-plate bending stiffness
EIzz      = 1.25e3;     %[Nm^2]  , out-of-plane bending stiffness
GJ        = 8e1;        %[Nm/rad], torsional stiffness
NSM       = 0.1;        %[kg/m]  , mass-per-unit-length
rotIxx    = 1.3e-4;     %[kgm]   , rotational inertia about xx-axis
rotIyy    = 5e-6;       %[kgm]   , rotational inertia about yy-axis
rotIzz    = 1.25e-4;    %[kgm]   , rotational inertia about zz-axis

%Applied loads
loadMag = 500;        %[N]     , magnitude of max applied load
nLoads  = 1;         %[-]     , number of loads to be applied

%Structural discretization
nStructuralNodes = 25;  %[-]     , number of structural nodes representing the curve

%Dependent parameters
G       = E;                            %[N/m^2] , shear modulus
A       = EA   / E;                     %[m^2]   , cross-sectional area
I11     = EIxx / E;                     %[m^4]   , 2nd moment of area in plane 1
I22     = EIzz / E;                     %[m^4]   , 2nd moment of area in plane 2
J       = GJ   / G;                     %[m^4]   , polar moment of area
dLoadMag = loadMag / nLoads;            %[N]     , increment for defining multiple loads
loadVec = dLoadMag : dLoadMag : loadMag;%[N]     , vector of applied loads (different cases, not all at once)

%% Set up graphics objects

%Geometry drawing
hF = figure('Name', 'Geometry Model');
hAx = axes('Parent', hF);

%Analysis Model drawing
% hF(2)  = figure('Name', 'Analysis Model');
% hAx(2) = axes('Parent', hF(2));

%Set axes appearance
set(hAx, ...
    'NextPlot', 'add'    , ...
    'XLim'    , [-3, 3], ...
    'YLim'    , [0, 3]  , ...
    'ZLim'    , [-5, 5]  , ...
    'View'    , [150, 35], ...
    'Box'     , 'on'     , ...
    'XGrid'   , 'on'     , ...
    'YGrid'   , 'on'     , ...
    'ZGrid'   , 'on');
set([hAx.XLabel], 'String', 'X [m]');
set([hAx.YLabel], 'String', 'Y [m]');
set([hAx.ZLabel], 'String', 'Z [m]');

%% Calculate beam geometry

%Straight beam
xd = [0, 0];
yd = [0, length];
zd = [0, 0];

%Plot to check
plot3(hAx(1), xd, yd, zd, '-');

%% Generate the ALENA objects
%   - 1 x awi.model.Beam             : For describing the geometry and stiffness properties
%   - 1 x awi.model.Material         : For describing the material properties
%   - nNode x awi.model.CoordSys     : For describing the beam orientation
%   - 1 x awi.model.LoadDistribution : For describing the loading conditions
% Beam   = awi.model.Beam;
Beam  = awi.model.LiftingSurface;
Mat   = awi.model.Material;
Load  = arrayfun(@(~) awi.loads.LoadDistribution, 1 : numel(loadVec));

%Ensure consistent discretization of the geometry
Beam.NumBeamElem = (nStructuralNodes - 1);

%Beam geometry
set(Beam, 'XData', xd, 'YData', yd, 'ZData', zd, 'Origin', [0, 0, 0]);

%LiftingSurface data
Beam.SpanVector = 'Y';
% Beam.Chord = [2, 2];
Beam.LE = [xd - 1 ; yd ; zd];
Beam.TE = [xd + 1 ; yd ; zd];
Beam.ActiveSet = 'cSet';
build(Beam);
[cs, eta] = createCoordSysObjects(Beam);
Beam.replaceBeamObject(cs, eta);
drawElement(Beam, hAx(1));

%Define stiffness properties (constant)
set(Beam, ...
    'A'  , [A  , A  ], 'A_eta'  , [0, 1], ...
    'I11', [I11, I11], 'I11_eta', [0, 1], ...
    'I22', [I22, I22], 'I22_eta', [0, 1], ...
    'J'  , [J  , J  ], 'J_eta'  , [0, 1], ...
    'NSM', [NSM, NSM], 'NSM_eta', [0, 1]);

%Set material properties and add to the object
set(Mat, 'E', E, 'G', G);
assignBeamObject(Beam, [Mat, Mat], [0, 1]);

%Applied loads
loadVal       = zeros(6, numel(loadVec));
loadVal(3, :) = loadVec;
set(Load, {'PointLoads'}      , num2cell(loadVal, 1)');
set(Load, 'EtaDistribution'   , 1);
set(Load, 'PointLoadCoordSys' , 'local');
set(Load, 'PointLoadBehaviour', 'non-follower');
set(Load, 'MaximumVectorLength', 2);

%Draw to check
drawElement(Beam, hAx(1));

%% Convert to analysis model

%FEM = convertToFE(Beam);

%Change visualisation of load
%set([FEM.PointLoads], 'MaxVectorLength', Load(1).MaximumVectorLength);

%Draw to check
%   - Overlay the nodes to make sure FE geometry matches the ALENA geometry
%draw(FEM, hAx(1), 'PartList', {'Nodes'});
%draw(FEM, hAx(2));

%% Run analysis

%Global applied load
i_runAnalysis(Beam, Load, loadVec, 'non-follower');

%Follower load
i_runAnalysis(Beam, Load, loadVec, 'follower');

end

function i_runAnalysis(Beam, Load, loadVec, loadType)

set(Load, 'PointLoadBehaviour', loadType);

%Set up analysis options
StaticOptions = awi.methods.Options;
StaticOptions.StructuralResponse = 'nonlinear';

%Don't model gravity
StaticOptions.bModelGravity = false;

%For each load...
nLoads  = numel(Load);
tipDisp = zeros(2, nLoads, 3);
for iL = 1 : nLoads
    
    %Remove any loads
    removeAllLoads(Beam);
    
    %Attach the load to the beam object
    assignLoadToBeam(Load(iL), Beam);
    
    %Make the model (including loads this time)
    FEM = convertToFE(Beam);

    %Make a directory for Nastran results
    resDir = awi.methods.Nastran.makeDefaultAnalysisDirectory('nonlinear_static', pwd);
    
    %Nastran analysis
    Nastran = awi.methods.Nastran;
    Nastran.AnalysisModel = FEM;
    NasResults = static(Nastran, StaticOptions, resDir);
    
    %Intrinsic, strain-based, FE
    StrainFE = awi.methods.IntrinsicStrainFE;
    StrainFE.AnalysisModel = FEM;
    StrainResults = static(StrainFE, StaticOptions, loadVec(iL));
    
    %Extract the tip displacement of the beam    
    nStructuralNodes = numel(FEM.BeamNodes) / 2 + 1;
    %   - Nastran
    tipDisp(1, iL, 1) = NasResults.Translation(1, nStructuralNodes, end);
    tipDisp(1, iL, 2) = NasResults.Translation(2, nStructuralNodes, end);
    tipDisp(1, iL, 3) = NasResults.Translation(3, nStructuralNodes, end);
    %   - Intrinsic Strain FE
    tipDisp(2, iL, 1) = StrainResults.Translation(1, end);
    tipDisp(2, iL, 2) = StrainResults.Translation(2, end);
    tipDisp(2, iL, 3) = StrainResults.Translation(3, end);
end

tipDisp = tipDisp .* 100;

%Plot the difference 
hF  = figure('Name', sprintf('Comparison of structural solvers for load type %s', upper(loadType)));
hg  = gobjects(size(tipDisp, 1), 3);
hAx = arrayfun(@(i) subplot(1, 3, i, 'Parent', hF), 1 : 3);
for i = 1 : 3    
    hg(:, i) = plot(hAx(i), loadVec, tipDisp(:, :, i));
    set(hg(:, i), {'MarkerFaceColor'}, get(hg(:, i), {'Color'}), ...
        'MarkerEdgeColor', 'k', {'Marker'}, {'s' ; '^'});
end
set(hg(:, end), {'DisplayName'}, {'Nastran' ; 'Instrinsic Strain FE'});
legend(hAx(end));
set([hAx.XLabel], 'String', 'Tip Load Magnitude [N]');
ylabel(hAx(1), 'x-displacement [%]');
ylabel(hAx(2), 'y-displacement [%]');
ylabel(hAx(3), 'z-displacement [%]');

end
