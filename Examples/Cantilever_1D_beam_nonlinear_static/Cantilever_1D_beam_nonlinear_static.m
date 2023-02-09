% Example for a static nonlinear (SOL 400) 1D tapered cantilever beam with
% a distribution of point loads and masses. This will show how to define
% stress recovery points and extract stress results from the ALENA output.

% Gobal coordinate system: +x forward, +y spanwise, +z down 

clear,clc,close all


%% Parameters
root_section_h = 0.015; % [m]
root_section_w = 0.020; % [m]
tip_section_h  = 0.003; % [m]
tip_section_w  = 0.010; % [m]

L1 = 0.05; % [m] start of section transition
L2 = 0.15; % [m] end of section transition


%% Properties

% Aluminium (6061)
E         	= 68.9E9; 	% [N/m^2] Young's modulus - required because the ALENA beam object accepts I/A/J and E/G values, not EI/EA/GJ.
G        	= 25.9E9;  	% [N/m^2] Shear modulus
rho     	= 2700;   	% [kg/m^3] Density
nu       	= 0.3;    	% Poisson's ratio

% Rectangular beam section properties a*b^3/12
root_beam_Ixx = root_section_w*root_section_h^3/12;     % [m^4] flapping second moment of area
root_beam_Izz = root_section_h*root_section_w^3/12;     % [m^4] lead-lag second moment of area
root_beam_J = 0.281*root_section_w*root_section_h^3;    % [m^4]
root_crossSectionArea = root_section_w*root_section_h; 	% [m^2] beam cross section area
tip_beam_Ixx = tip_section_w*tip_section_h^3/12;        % [m^4] flapping second moment of area
tip_beam_Izz = tip_section_h*tip_section_w^3/12;        % [m^4] lead-lag second moment of area
tip_beam_J = 0.281*tip_section_w*tip_section_h^3;       % [m^4]
tip_crossSectionArea = tip_section_w*tip_section_h;   	% [m^2] beam cross section area


% Point load distribution
loadPositions = [0.2, 0.4, 0.6, 0.7, 0.8]; % [m] row array of locations of loads application on beam shear centre
loads = [   0   0   0   0   0
            0   0   0   0   0
            -1  -1  -1  -1  -1
            0   0   0   0   0
            0   0   0   0   0
            0   0   0   0   0   ]; % [forces: N, moments: Nm] row array of local loads [Fx;Fy;Fz;Mx;My;Mz] acting at loadPositions

% Point masses
masses = [0.1, 0.1, 0.1, 0.1, 0.1]; % [kg] row array of point masses
massPositions = [   -0.010	-0.010  -0.010 -0.010   -0.010
                     0.100   0.300   0.400  0.500    0.700
                     0.000   0.000   0.000  0.000    0.000  ]; % [m] row array of [x;y;z] locations for point masses
     
% Tapered beam properties distribution
L                   = 0.8;                      % [m] Beam length
nStructuralNodes    = 50;                       % [] Number of structural nodes / grid points
taperPositions      = [0, L1, L2, L];     	% [m] row array of locations of taper properties
A                   = [root_crossSectionArea, root_crossSectionArea, tip_crossSectionArea, tip_crossSectionArea];     % [m^2] Cross-sectional area at taperPositions
I11                 = [root_beam_Ixx, root_beam_Ixx, tip_beam_Ixx, tip_beam_Ixx];  	% [Nm^2] x second moment of area at taperPositions (element Izz) about the neutral axis
I22                 = [root_beam_Izz, root_beam_Izz, tip_beam_Izz, tip_beam_Izz]; 	% [Nm^2] z second moment of area at taperPositions (element Iyy) about the neutral axis
J                   = [root_beam_J, root_beam_J, tip_beam_J, tip_beam_J];        	% [Nm^2] Torsional constant at taperPositions
SC                  = [ 0.000   0.000   0.000   0.000
                        0.000   0.000   0.000   0.000   ];  % [m] row array of shear centre vectors (y(vertial);z(horizontal)) relative to centroid/neutral axis at taperPositions
NSM                 = [ 0.0     0.0     0.0     0.0 	]; 	% [kg/m] 'NSM' %Non-structural mass (per unit length)
NSI                 = [ 0.0     0.0     0.0     0.0 	]; 	% [kg.m] 'NSI' %Non-structural inertia (per unit length)
% Calculates the non-structural mass offset for Nastran 
CMy                 = [ 0.0     0.0     0.0     0.0 	]; 	% [m] 'CMy' % Non-structural centre-of-mass offset in the y-direction of beam coordinate system (positive upward)
CMz                 = [ 0.0     0.0     0.0     0.0 	];  % [m] 'CMz' % Non-structural centre-of-mass offset in the z-direction of beam coordinate system (positive forward)
bar_width        	= [root_section_w root_section_w tip_section_w tip_section_w];           % [m] rectangular section dimensions for stress recovery points at taperPositions
bar_height        	= [root_section_h root_section_h tip_section_h tip_section_h];           % [m] rectangular section dimensions for stress recovery points at taperPositions
% Stress recovery points (relative to neutral axis): z is horizontal axis (+ve towards trailing edge), y is vertical axis (+ve up)
Cz = -bar_width./2;  Cy = bar_height./2;       % Top right corner
Dz = -bar_width./2;  Dy = -bar_height./2;      % Bottom right corner
Ez = bar_width./2; Ey = -bar_height./2;        % Bottom left corner
Fz = bar_width./2; Fy = bar_height./2;         % Top left corner


%% Set up graphics objects

% %Geometry drawing
% hF = figure('Name', 'Geometry Model');
% hAx = axes('Parent', hF);
% 
% %Analysis Model drawing
% hF(2)  = figure('Name', 'Analysis Model');
% hAx(2) = axes('Parent', hF(2));
% 
% %Set axes appearance
% set(hAx, ...
%     'NextPlot', 'add'    , ...
%     'XLim'    , [-1, 1], ...
%     'YLim'    , [0, 1]  , ...
%     'ZLim'    , [-1, 1]  , ...
%     'View'    , [150, 35], ...
%     'Box'     , 'on'     , ...
%     'XGrid'   , 'on'     , ...
%     'YGrid'   , 'on'     , ...
%     'ZGrid'   , 'on');
% set([hAx.XLabel], 'String', 'X [m]');
% set([hAx.YLabel], 'String', 'Y [m]');
% set([hAx.ZLabel], 'String', 'Z [m]');

%% Calculate beam geometry

yd       = linspace(0,L,nStructuralNodes);
xd       = linspace(0,0,nStructuralNodes);
zd       = linspace(0,0,nStructuralNodes);

%Plot to check
% plot3(hAx(1), xd, yd, zd, '-');

%% Generate the ALENA objects
%   - 1 x awi.model.Beam             : For describing the geometry and stiffness properties
%   - 1 x awi.model.Material         : For describing the material properties
%   - nNode x awi.model.CoordSys     : For describing the beam orientation
%   - 1 x awi.model.LoadDistribution : For describing the loading conditions
Beam   = awi.model.Beam;
% Beam = awi.model.LiftingSurface;
Mat    = awi.model.Material;
Orient = arrayfun(@(~) awi.model.CoordSys, 1 : nStructuralNodes);
Load   = awi.loads.LoadDistribution;

%Beam geometry
set(Beam, 'XData', xd, 'YData', yd, 'ZData', zd, 'Origin', [0, 0, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %LiftingSurface data
% Beam.SpanVector = 'Y';
% % Beam.Chord = [2, 2];
% Beam.LE = [xd - 0.2 ; yd ; zd];
% Beam.TE = [xd + 0.2 ; yd ; zd];
% Beam.ActiveSet = 'cSet';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build(Beam);
% drawElement(Beam, hAx(1));

%Set up coordinate systems
%   - Position along the beam
r     = Beam.RData;
r_bar = r ./ r(end); % Normalising
set(Orient, {'SOffset'}, num2cell(r_bar)');

% %   - Orientation
% %       * Local x along beam
% %       * Local y in same plane as global X axis
% v1 = [diff(xd) ; diff(yd) ; diff(zd)]; % v1 = 
% v1 = num2cell([v1 , v1(:, end)]', 2);
% v2 = [-1, 0, 0];
% set(Orient, {'V1'}, v1, 'V2', v2, 'GridVector', 'x', 'GridPlane', 'xy', 'ActiveSet', 'vSet');

%   - Orientation
%       * Local z same as global axis (G1)
%       * Local x same as global axis (G2), used to define xz plane
G1 = [0, 0, 1];
G2 = [1, 0, 0];
set(Orient, 'G1', G1, 'G2', G2, 'ActiveSet', 'gSet');


%   - Build the object and assign to the beam
arrayfun(@build, Orient);
assignBeamObject(Beam, Orient, r_bar);

taperPositionsNorm = taperPositions./L;

% Shear centres along beam (relative to grid point, set by CBEAM)
SC_y = SC(1,:);
SC_z = SC(2,:);

% Neutral axis (relative to section shear centre, set by PBEAM)
NA = -SC; % Shift NA in the opposite to SC so that neutral axis is at grid point
NA_y = NA(1,:);
NA_z = NA(2,:);

% Non-structural mass/inertia properties (alternative is to use individual point loads)
NSM = NSM;          % [kg/m] 'NSM' %Non-structural mass (per unit length)
NSI = NSI;          % [kg.m] 'NSI' %Non-structural inertia (per unit length)
CMy = CMy-SC_y;     % [m] 'CMy' % Non-structural centre-of-mass offset in the y-direction of beam coordinate system (relative to shear centre, but I substract the shear centre location so points are relative to neutral axis)
CMz = CMz-SC_z;     % [m] 'CMz' % Non-structural centre-of-mass offset in the z-direction of beam coordinate system (relative to shear centre, but I substract the shear centre location so points are relative to neutral axis)

% Stress recovery points relative to shear centre (I substract the shear centre location so points are relative to neutral axis)
Cy = Cy-SC_y;
Cz = Cz-SC_z;
Dy = Dy-SC_y;
Dz = Dz-SC_z;
Ey = Ey-SC_y;
Ez = Ez-SC_z;
Fy = Fy-SC_y;
Fz = Fz-SC_z;

%Define stiffness properties (constant)
set(Beam, ...
    'A'  , A    , 'A_eta'  , taperPositionsNorm, ...
    'I11', I11  , 'I11_eta', taperPositionsNorm, ...
    'I22', I22  , 'I22_eta', taperPositionsNorm, ...
    'J'  , J    , 'J_eta'  , taperPositionsNorm, ...
    'SCy', SC_y , 'SCy_eta', taperPositionsNorm, ...
    'SCz', SC_z , 'SCz_eta', taperPositionsNorm, ...
    'NAy', NA_y , 'NAy_eta', taperPositionsNorm, ...
    'NAz', NA_z , 'NAz_eta', taperPositionsNorm, ...
    'NSM', NSM  , 'NSM_eta', taperPositionsNorm, ...
    'NSI', NSI  , 'NSI_eta', taperPositionsNorm, ...
    'CMy', CMy  , 'CMy_eta', taperPositionsNorm, ...
    'CMz', CMz  , 'CMz_eta', taperPositionsNorm, ...
    'Cy'  , Cy  , 'Cy_eta'  , taperPositionsNorm, ...
    'Cz'  , Cz  , 'Cz_eta'  , taperPositionsNorm, ...
    'Dy'  , Dy  , 'Dy_eta'  , taperPositionsNorm, ...
    'Dz'  , Dz  , 'Dz_eta'  , taperPositionsNorm, ...
    'Ey'  , Ey  , 'Ey_eta'  , taperPositionsNorm, ...
    'Ez'  , Ez  , 'Ez_eta'  , taperPositionsNorm, ...
    'Fy'  , Fy  , 'Fy_eta'  , taperPositionsNorm, ...
    'Fz'  , Fz  , 'Fz_eta'  , taperPositionsNorm);

% %Full list of beam properties:
% 'A'   %Cross-sectional area
% 'I11' %Second moment of area in plane 1
% 'I22' %Second moment of area in plane 2
% 'I12' %Second moment of area in plane 1-2
% 'J'   %Polar second moment of area
% 'NSM' %Non-structural mass (per unit length)
% 'NSI' %Non-structural inertia (per unit length)
% 'RI1' %Rotational inertia about the 1st axis
% 'RI2' %Rotational inertia about the 2nd axis
% 'R12' %Rotational inertia about the 3rd axis
% 'SCy' %Shear centre offset in the y-direction of beam coordinate system
% 'SCz' %Shear centre offset in the z-direction of beam coordinate system
% 'NAy' %Neutral axis offset in the y-direction of beam coordinate system
% 'NAz' %Neutral axis offset in the z-direction of beam coordinate system
% 'CMy' %Centre-of-mass offset in the y-direction of beam coordinate system
% 'CMz' %Centre-of-mass offset in the z-direction of beam coordinate system
% %Add stress recovery points relative to the shear center
% 'Cy' % Top right corner
% 'Cz'
% 'Dy' % Bottom right corner
% 'Dz'
% 'Ey' % Bottom left corner
% 'Ez'
% 'Fy' % Top left corner
% 'Fz'
            
%Set material properties and add to the object
set(Mat, 'E', E, 'G', G, 'Rho', rho, 'Nu', nu);
assignBeamObject(Beam, [Mat, Mat], [0, 1]);

%Applied loads
Load.EtaDistribution     = loadPositions./L;                % Normalised load distribution
Load.PointLoads          = loads;                           % row array of load vectors
Load.PointLoadCoordSys   = 'local';                         %orientation is in local system
Load.PointLoadBehaviour  = 'follower';                      %follower force
Load.AxisFlag            = 'R';                             %defined along beam piecewise length
Load.MaximumVectorLength = 5;                               %scale the plot so we can clearly see it
%   - Attach the load to the beam object
assignLoadToBeam(Load, Beam);

%Point masses
for i = 1:length(masses)
    pointMass = awi.model.PointMass;
    set(pointMass,'Mass', masses(i));
    set(pointMass,'Origin', massPositions(:,i)');
    Beam.add(pointMass);
end

%Draw to check
% drawElement(Beam, hAx(1));

%% Convert to analysis model

%Ensure consistent discretization of the geometry
Beam.NumBeamElem = (nStructuralNodes - 1);

%Make the model
FEM = convertToFE(Beam);


%Setting my own element orientation vection (v) which defines the element xy plane in Nastran...
%The x axis is already known by the CBEAM element as the A-to-B grid points vector...
%So I only need v to be the element y axis which should be [0; 0; -1] so it
%points up in awi.fe.Beam (FEM variable here) which takes X from awi.model.Beam (Beam variable here)
%... so positive y is upwards, and positive z is backwards (negative x direction)
elementOrientationVector = [0 ; 0 ; -1];
for i = 1:length(FEM.Beams)
    FEM.Beams(i).X = elementOrientationVector;
end

%Change visualisation of load
set([FEM.PointLoads], 'MaxVectorLength', Load.MaximumVectorLength);

%Draw to check
%   - Overlay the nodes to make sure FE geometry matches the ALENA geometry
% draw(FEM, hAx(1), 'PartList', {'Nodes'});
% draw(FEM, hAx(2));

%% Run analysis

%Set up analysis options
StaticOptions = awi.methods.Options;
StaticOptions.StructuralResponse = 'nonlinear';
StaticOptions.bModelGravity = true;
StaticOptions.GravityVector = [0; 0; 9.81; 0; 0; 0];
% StaticOptions.CTRLDEF = 'MILDLY'; % Option for controlling load iteration in SOL 400, {'', 'QLINEAR', 'MILDLY', 'SEVERELY'}


%Nastran analysis
Nastran = awi.methods.Nastran;
Nastran.AnalysisModel = FEM;
fileDir = './Examples/Cantilever_1D_beam_nonlinear_static';
NasResults = static(Nastran, StaticOptions, fileDir); % Creates bdf files and runs Nastran


loadStepEnd = NasResults.displacement.TimeVector(end);
if loadStepEnd ~= 1
    disp('Nastran loading did not converge! (check .f06 for fatal error or divergence)');
    return
end


%% Displacement
tipDisplacementDelta = NasResults.displacement.Translation(:,end,end);
tipInitialPosition = NasResults.displacement.Nodes(end).X;
tipDisplacement = tipInitialPosition + tipDisplacementDelta;
fprintf('Tip deflection vector [x y z] = [%0.2f %0.2f %0.2f] mm\n',-tipDisplacement(1)*1000,tipDisplacement(2)*1000,-tipDisplacement(3)*1000);                         


%% Maximum equivalent stress (NSE)
stressResults = NasResults.stress(end);

% Initial variables
maxStressGRID = [];
maxStressPosition = [];
maxStress = 0;

% Looping through A grid points
GridA = stressResults.GRIDA;
for i = 1:length(GridA)
    Grid = GridA(i);
    stressC = stressResults.NSECA(i); if stressC > maxStress; maxStress = stressC; maxStressGRID = Grid; maxStressPosition = 'C'; end % Position C
    stressD = stressResults.NSEDA(i); if stressD > maxStress; maxStress = stressD; maxStressGRID = Grid; maxStressPosition = 'D'; end % Position D
    stressE = stressResults.NSEEA(i); if stressE > maxStress; maxStress = stressE; maxStressGRID = Grid; maxStressPosition = 'E'; end % Position E
    stressF = stressResults.NSEFA(i); if stressF > maxStress; maxStress = stressF; maxStressGRID = Grid; maxStressPosition = 'F'; end % Position F
end

% Looping through B grid points
GridB = stressResults.GRIDB;
for i = 1:length(GridB)
    Grid = GridB(i);
    stressC = stressResults.NSECB(i); if stressC > maxStress; maxStress = stressC; maxStressGRID = Grid; maxStressPosition = 'C'; end % Position C
    stressD = stressResults.NSEDB(i); if stressD > maxStress; maxStress = stressD; maxStressGRID = Grid; maxStressPosition = 'D'; end % Position D
    stressE = stressResults.NSEEB(i); if stressE > maxStress; maxStress = stressE; maxStressGRID = Grid; maxStressPosition = 'E'; end % Position E
    stressF = stressResults.NSEFB(i); if stressF > maxStress; maxStress = stressF; maxStressGRID = Grid; maxStressPosition = 'F'; end % Position F
end

% Position along beam
Grids = cell2mat({NasResults.displacement.Nodes.GID})';
idxMaxStressGRID = find(Grids == maxStressGRID);
GridPositions = cell2mat({NasResults.displacement.Nodes.X})';
maxStressBeamPosition = GridPositions(idxMaxStressGRID,:);

% Printing out results
fprintf('\nMax stress = %0.3G N/m^2 @ beam location [x,y,z] = [%0.3f %0.3f %0.3f] m (stress recovery point "%s") \n\n',maxStress,maxStressBeamPosition(1),maxStressBeamPosition(2),maxStressBeamPosition(3),maxStressPosition);




