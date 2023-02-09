function [Aircraft, FuelDistribution] = FameModel2AWIFramework(fame, fameFEModel, logfcn)
%FameModel2AWIFramework The following function is in charge of converting
%a FAME structure, derived directly from the *.fm4 file, into one that is
%useable by the AWI framework.
%
% TODO - Add the HTP to the VTP if it is present
% TODO - Sort out fuselage cross-section generation

if nargin < 3
    logfcn = @(s) fprintf('%s\n', s);
end

%Check fame model matches expected format
%   - At a minimum we require information about the wing
i_parseFameGeometry(fame.Geometry)

    function i_parseFameGeometry(Geometry)
        expFields = {'Wing'};
        fNames    = fieldnames(Geometry);
        assert(all(contains(expFields, fNames)), sprintf(['FAME model ', ...
            'does not match the expected format. Expected the FAME'    , ...
            'geometry to contain data about the %s components.'], strjoin(expFields)));
    end

%% Aircraft
logfcn([blanks(8) 'Initial Aircraft created ...']);
Aircraft = awi.model.Aircraft('Name','Aircraft');

logfcn([blanks(8) 'Fuselage added ...']);
F = Aircraft.add(@awi.model.BluffBody,'Name','Fuselage');

map = { ...
    'mtow', 'MTOM';...
    'mzfw', 'MZFM';...
    'mlw' , 'MLM';...
    'owe' , 'OEM'};

fnames = fieldnames(fame.Geometry.Weights);
[b,~]  = ismember(map(:,1),fnames);
iCS = find(b)';
for i = iCS
    Aircraft.(map{i,2}) = fame.Geometry.Weights.(fnames{i});
end

%% Lifting Surfaces

%Assign a Starboard Wing to the Fuselage
W = F.add(@awi.model.LiftingSurface, 'Name', 'StbdWing');
assignLiftingSurfaceData(W, fame.Geometry.Wing, 'Y', W.Name);
logfcn([blanks(8) 'Starboard wing successfully initialised.']);

%Assign a PortWing to the Fuselage
Wp = F.add(@awi.model.LiftingSurface,'Name','PortWing');
assignLiftingSurfaceData(Wp, fame.Geometry.Wing, 'Y', Wp.Name);
setupPortLiftingSurface(Wp);
logfcn([blanks(8) 'Port wing successfully initialised.']);

%VTP
if isfield(fame.Geometry, 'VTP')
    V = F.add(@awi.model.LiftingSurface, 'Name', 'VTP');
    assignLiftingSurfaceData(V, fame.Geometry.VTP, 'Z', V.Name);
    %Make the assumption that the VTP is always vertical
    set(V, 'Dihedral', zeros(size(V.Dihedral_eta)));
    %Double the span as there is not port component for the VTP
    V.Span = 2 * V.Span;
    logfcn([blanks(8) 'VTP successfully initalised.']);
end

%HTP
if isfield(fame.Geometry, 'HTP')
    
    %Starboard
    H = F.add(@awi.model.LiftingSurface, 'Name', 'StbdHTP');
    assignLiftingSurfaceData(H, fame.Geometry.HTP, 'Y', H.Name);
    logfcn([blanks(8) 'Starboard HTP successfully initialised.']);
    
    %Port
    Hp = F.add(@awi.model.LiftingSurface, 'Name', 'PortHTP');
    assignLiftingSurfaceData(Hp, fame.Geometry.HTP, 'Y', Hp.Name);
    setupPortLiftingSurface(Hp);
    logfcn([blanks(8) 'Port HTP successfully initialised.']);
    
end

%% FUSELAGE PROPERTIES
%   - FAME input provides limited information about the fuselage. Need to
%     infer some basic properties based on where the current lifting
%     surface objects are placed.
%
% TODO - Overhaul this section when the fuselage geometrt parameterisation
% is established.

%Generate fuselage data - Use TUX xml data if available
FuselageParam = getFuselageData(fame, F, W, logfcn);

% %Data from Dario's fuselage model
% eta = [0;0.005268;0.010536;0.015805;0.021073;...
%     0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
%     0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
%     0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
%     0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';
% rad = [0,0.1818,0.2672,0.3348,0.3928,...
%     0.4447,0.4921,0.5361,0.5774,...
%     0.6165,0.6537,0.6893,0.7234,...
%     0.7563,0.7881,0.8190,0.8489,...
%     0.8780,0.9063,0.934,1,1,1,1,1,1,1,1,1,...
%     1,1,1,1,1,0.8985,0.8275,0.7566,...
%     0.6620,0.5674,0.4729,0];
% 
% %Want to determine scaling values for the fore, mid & aft fuselage sections
% idxMid  = rad == 1;
% indFore = 1 : find(idxMid, 1, 'first') - 1;
% indAft  = find(idxMid, 1, 'last') + 1 : numel(eta);
% etaFore = (eta(indFore) - eta(indFore(1))) ./ (eta(indFore(end)) - eta(indFore(1)));
% etaAft  = (eta(indAft)  - eta(indAft(1)))  ./ (eta(indAft(end))  - eta(indAft(1)));
% radFore_bar = rad(indFore);
% radAft_bar  = rad(indAft);
% rad_bar     = [radFore_bar, radAft_bar];
% etaMid  = eta(idxMid);
% etaMid  = [etaMid(1), etaMid(end)];

%Set up fuselage parameters
% totalLength = FuselageParam.ForeLength + FuselageParam.MidLength + FuselageParam.AftLength;
% xVec = [ ...
%     (etaFore .* FuselageParam.ForeLength), ...
%     (etaAft .* FuselageParam.AftLength) + FuselageParam.ForeLength + FuselageParam.MidLength];
% etaVec = xVec ./ totalLength;

%Set up fuselage parameters
xVec   = cumsum([0, FuselageParam.ForeLength, FuselageParam.MidLength, FuselageParam.AftLength]);
etaVec = xVec ./ xVec(end);

% %Assign to the object
% F.Origin = FuselageParam.Origin;
% F.Length = xVec(end);
% F.Eta    = etaVec;
% 
% %Make cross-section objects for the mid-section only
% etaMid = etaVec(ismember(xVec, [FuselageParam.ForeLength, FuselageParam.ForeLength + FuselageParam.MidLength]));
% Ellipse = arrayfun(@(~) awi.model.CrossSection, 1 : numel(etaMid));
% set(Ellipse, 'CrossSectionName', 'Ellipse');
% set(Ellipse, 'MajorAxis', FuselageParam.Height);
% set(Ellipse, 'MinorAxis', FuselageParam.Width);
% generateCoordsFromLibrary(Ellipse);
% assignBeamObject(F, Ellipse, etaMid, 'replace');

GuessFuseLength = xVec(end);
GuessFuseRadius = FuselageParam.Width / 2;
%Set up properties
F.Length = GuessFuseLength;
F.Eta = [0;0.005268;0.010536;0.015805;0.021073;...
    0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
    0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
    0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
    0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1]';

F.Radius = [0;0.3844030;0.565081;0.707928;0.830682;0.940375;...
    1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
    1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
    2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
    1.75;1.6;1.4;1.2;1.0;0]';

F.Radius = GuessFuseRadius*[0,0.1818,0.2672,0.3348,0.3928,...
    0.4447,0.4921,0.5361,0.5774,...
    0.6165,0.6537,0.6893,0.7234,...
    0.7563,0.7881,0.8190,0.8489,...
    0.8780,0.9063,0.934,1,1,1,1,1,1,1,1,1,...
    1,1,1,1,1,0.8985,0.8275,0.7566,...
    0.6620,0.5674,0.4729,0];

F.Origin = [0,0,0];

logfcn([blanks(8) 'Fuselage successfully initialised.']);

%% ENGINE PROPERTIES
if ~isempty(fame.Geometry.Engine)
    
    E = W.add(@awi.model.BluffBody,'Name','StbdEngine');
    E.Length  = 5;
    E.SOffsetFlag = W.SpanVector;
    E.SOffset = fame.Geometry.Engine.eta;
    E.XOffset = fame.Geometry.Engine.pylon_x_loc - E.Length/2;
    E.ZOffset = fame.Geometry.Engine.engine_z_loc;
    E.Eta     = [0,1];
    E.Radius  = [1.5,1.5];
    E = Wp.add(@awi.model.BluffBody,'Name','PortEngine');
    E.Length  = 5;
    E.SOffsetFlag = Wp.SpanVector;
    E.SOffset = fame.Geometry.Engine.eta;
    E.XOffset = fame.Geometry.Engine.pylon_x_loc - E.Length/2;
    E.ZOffset = fame.Geometry.Engine.engine_z_loc;
    E.Eta     = [0,1];
    E.Radius  = [1.5,1.5];
end

logfcn([blanks(8) 'Aircraft definition complete.']);

%% Stiffness/mass properties

%Assign beam stiffness
%   - TODO : Add non-structural mass
assignBeamProperties(Aircraft, fameFEModel);

%Assign the lumped mass objects to the AWI model
% TODO - Make this work for the case where a FAME FE Model is not generated
if ~isempty(fameFEModel.Grid)
    [Aircraft, FuelDistribution] = assignLumpedMass(Aircraft, fameFEModel, logfcn);
end

end

%% Local functions

%   - Lifting surfaces
function assignLiftingSurfaceData(LiftSurf, FameGeometry, spanVector, compName)
%assignLiftingSurfaceData Maps the FAME geometry data to the ALENA
%geometry and generates control surface objects which are added to
%the wing.

%Control Surface token-object mapping
tokObj = { ...
    'slat'   , 'awi.model.Slat' ; ...
    'flap'   , 'awi.model.Flap' ; ...
    'aileron', 'awi.model.Aileron' ; ...
    'spoiler', 'awi.model.Spoiler'};

%Get spar data
eta = unique([FameGeometry.fspar(:, 1) ; FameGeometry.aspar(:, 1)]);
frontSpar = interp1(FameGeometry.fspar(:, 1), FameGeometry.fspar(:, 2), eta);
rearSpar  = interp1(FameGeometry.aspar(:, 1), FameGeometry.aspar(:, 2), eta);

%Place the beam at the midpoint between the spars
beamLoc = (frontSpar + rearSpar) ./ 2;

%Calculate offset from the leading edge to the beam location
dxLE = [(beamLoc(1) * FameGeometry.chord(1, 2)), 0, 0];

%Assign data to the lifting surface
%   - Set the 'SpanVector' first to propogate changes to all
%     '_Flag' properties belonging to the 'Planform Property' type.
%     This is to do with the event/listener behaviour which is
%     inherent to ALENA LiftingSurface objects
LiftSurf.SpanVector   = spanVector;
LiftSurf.Span         = FameGeometry.span / 2;
LiftSurf.Origin       = FameGeometry.ref_point + dxLE;  %FAME origin is at the LE whereas ALENA origin is at the beam location
LiftSurf.Chord_eta    = FameGeometry.chord(:, 1);
LiftSurf.Chord        = FameGeometry.chord(:, 2);
LiftSurf.Sweep_eta    = FameGeometry.sweep(:, 1);
LiftSurf.Sweep        = FameGeometry.sweep(:, 2);
LiftSurf.SweepLoc     = FameGeometry.sweep(1, 3);
LiftSurf.Dihedral_eta = FameGeometry.dihedral(:, 1);
LiftSurf.Dihedral     = FameGeometry.dihedral(:, 2);
LiftSurf.BeamLoc_eta  = eta;                            %Define the beam location at every eta position
LiftSurf.BeamLoc      = beamLoc;                        %Assue beam is at middle of both spar locations

%Make 'awi.model.Spar' objects for the front and rear spar
%   - Front spar
Fs      = LiftSurf.add(@awi.model.Spar, 'Name', [compName, 'FrontSpar']);
Fs.Eta  = FameGeometry.fspar(:, 1);
Fs.XLoc = FameGeometry.fspar(:, 2);
%   - Rear spar
Rs      = LiftSurf.add(@awi.model.Spar, 'Name', [compName, 'RearSpar']);
Rs.Eta  = FameGeometry.aspar(:, 1);
Rs.XLoc = FameGeometry.aspar(:, 2);

%Any control surfaces?
if ~isfield(FameGeometry, 'ControlSurfaces')
    return
end

%If we get this far then we have to define the control surfaces
csTypes = fieldnames(FameGeometry.ControlSurfaces);

%Loop through each control surface type and generate objects & data
for iCS = 1 : numel(csTypes)
    
    ControlSurfData = FameGeometry.ControlSurfaces.(csTypes{iCS});
    
    %Determine the type of control surface object that should be made
    idx = ismember(tokObj(:, 1), csTypes{iCS});
    if any(idx)
        %Make the specific control surface type according to the
        %mapping defined in 'tokObj'.
        fn = str2func(tokObj{idx, 2});
    else
        %Make a generic control surface is there is no token match
        fn = @awi.model.ControlSurface;
    end
    
    %How many of this type are required?
    nCS = structfun(@numel, FameGeometry.ControlSurfaces.(csTypes{iCS}));
    assert(range(nCS) == 0, sprintf(['Incorrect amount of ', ...
        'data in the FAME control surface (''%s'') geometry for ', ...
        'the %s component.'], csTypes{iCS}, compName));
    
    %Make the objects and name them in number order
    ControlSurf = arrayfun(fn, 1 : nCS(1));
    csName      = arrayfun(@(i) [csTypes{iCS}(1 : 3), num2str(i)], 1 : nCS(1), 'Unif', false);
    set(ControlSurf, {'Name'}, csName');
    
    %Setup eta and chord data for each control surface
    eta = [ControlSurfData.eta_beg_le, ControlSurfData.eta_end_le];
    if isfield(ControlSurfData, 'chrd_beg_le')
        xLE = [ControlSurfData.chrd_beg_le, ControlSurfData.chrd_end_le];
    else
        xLE = zeros(nCS(1), 2);
    end
    if isfield(ControlSurfData, 'chrd_beg_te')
        xTE = [ControlSurfData.chrd_beg_te, ControlSurfData.chrd_end_te];
    else
        xTE = ones(nCS(1), 2);
    end
    
    %Assign to object
    set(ControlSurf, {'Eta'}, num2cell(eta, 2));
    set(ControlSurf, {'xLE'}, num2cell(xLE, 2));
    set(ControlSurf, {'xTE'}, num2cell(xTE, 2));
    
    %Add to the parent lifting surface
    LiftSurf.add(ControlSurf);
    
end

end

function setupPortLiftingSurface(LiftSurf)
%setupPortLiftingSurface Update the LiftingSurface properties to
%specify a port-side object, instead of a starboard object.
%
%   - Sets the span to be a negative value
%   - Prefaces the control surfaces with the letter 'P' to ensure
%     unique control surface names.

LiftSurf.Span         = -LiftSurf.Span;
PortControlSurf = LiftSurf.ControlSurfaces;
set(PortControlSurf, {'Name'}, strcat('P', get(PortControlSurf, {'Name'})));
end

%   - Fuselage
function FuselageParam = getFuselageData(fame, Fuselage, StbdWing, logfcn)
%getFuselageData Returns a MATLAB structure describing the basic
%parameters of the fuselage. This data is either extracted from the
%TUX xml of estimated based on typical values/current geometry.

%Traverse the TUX hierarchy until we get to the fuselage parameters
while true
    if ~isfield(fame, 'TUXInput')
        break
    end
    if ~isfield(fame.TUXInput, 'airplane')
        break
    end
    if ~isfield(fame.TUXInput.airplane, 'fuselage')
        break
    end
    if ~all(isfield(fame.TUXInput.airplane.fuselage, {'width', 'height', 'length', 'fore', 'aft'}))
        break
    end
    %Safe to proceed
    FuselageParam = i_extractFuselageFromTUX(fame.TUXInput.airplane.fuselage);
    logfcn([blanks(8) 'Fuselage parameters extracted from the TUX .xml file.']);
    return
end

% error('Update code');

%If we get to this point the TUX xml is not present or does not
%match expected formatting. Guess fuselage dimensions instead...

%Make an initial guess about the length of the fuselage
%   - Assume it is at least as long as the rearmost component
pos = vertcat(Fuselage.Children.Position);
if size(pos, 1) == 1
    GuessFuseLength = 2 * max(pos(:,1));
else
    GuessFuseLength = max(pos(:,1));
end

%Make an initial guess about the radius of the fuselage
if isfield(fame.Geometry.UserVariables, 'root_wing')
    GuessFuseRadius = StbdWing.Span * (fame.Geometry.UserVariables.root_wing);
elseif isfield(fame.Geometry.UserVariables, 'root_wing_ib')
    GuessFuseRadius = StbdWing.Span * (fame.Geometry.UserVariables.root_wing_ib);
else
    GuessFuseRadius = 3;
end

%Assign to structure
FuselageParam.ForeLength = GuessFuseLength * 0.2;
FuselageParam.MidLength  = GuessFuseLength * 0.5;
FuselageParam.AftLength  = GuessFuseLength * 0.3;
FuselageParam.Origin     = [0 0 0];
FuselageParam.Height     = GuessFuseRadius;
FuselageParam.Width      = GuessFuseRadius;

logfcn([blanks(8) 'TUX .xml file did not match expected format ', ...
    'guessing the fuselage geometry instead.'])

    function FuselageParam = i_extractFuselageFromTUX(TuxXmlData)
        %i_extractFuselageFromTUX Maps the data in the TUX .xml to the
        %ALENA Fuselage object.
        %
        % A fuselage is parameterised by the following properties:
        %   - origin
        %   - width
        %   - height
        %   - length
        %   - fore length
        %   - mid length
        %   - aft length
        
        L = str2double(TuxXmlData.length.Attributes.value);
        
        FuselageParam.Origin     = str2num(TuxXmlData.AxisOffset.Attributes.value); %#ok<ST2NM>
        FuselageParam.Width      = str2double(TuxXmlData.width.Attributes.value);
        FuselageParam.Height     = str2double(TuxXmlData.height.Attributes.value);
        FuselageParam.ForeLength = str2double(TuxXmlData.fore.Attributes.value);
        FuselageParam.AftLength  = str2double(TuxXmlData.aft.Attributes.value);
        FuselageParam.MidLength  = L - FuselageParam.ForeLength - FuselageParam.AftLength;
    end

end

%   - Generic beam properties
function assignBeamProperties(Aircraft, fameFEModel)
%assignBeamProperties Defines the beam stiffness properties (I1/I2/A/J) of
%the various beam objects using the FEM data from FAME and also adds
%material properties along the beam.

%Get all the 'awi.model.BeamObjects'
Parts = flatlist(Aircraft);
Beams = Parts(arrayfun(@(o) isa(o, 'awi.model.Beam'), Parts));
clear Parts
beamNames = get(Beams, {'Name'});

%What parts are in the FAME FEM?
partNames = {fameFEModel.PartId.part};
idx = ~ismember(beamNames, partNames);
beamNames(idx) = [];
Beams(idx)     = [];

%Only concerned with grid & element data
partType = {fameFEModel.PartId.type};
idx      = ismember(partType, {'GRID', 'CBAR'});
PartId   = fameFEModel.PartId(idx);

%Construct beam geometry and property distribution for each beam in the FEM
for iB = 1 : numel(beamNames)
    %Find the grids and beam elements
    idx = ismember({PartId.part}, beamNames{iB});
    BeamPart = PartId(idx);
    BeamGrid = BeamPart(ismember({BeamPart.type}, 'GRID'));
    BeamBar  = BeamPart(ismember({BeamPart.type}, 'CBAR'));
    if isempty(BeamGrid) || isempty(BeamBar)
        warning(['Unable to find FE data for the %s component in the ', ...
            'FAME FEM. Skipping the generation of beam properties for ', ...
            'this component.'], Beams(iB).Name);
        continue
    end
    %Get FE data -> associated via ID number (just like Nastran bulk data)
    Grid = fameFEModel.Grid(ismember([fameFEModel.Grid.id], BeamGrid.data));
    Bar  = fameFEModel.Cbar(ismember([fameFEModel.Cbar.id], BeamBar.data));
    PBar = fameFEModel.Pbar(ismember([fameFEModel.Pbar.id], [Bar.pid]));
    Mat  = fameFEModel.Mat(ismember([fameFEModel.Mat.id] , unique([PBar.mat])));
    %Assume bar objects are defined in order along the beam
    coords = vertcat(Grid.coord);
    %Normalise the distribution w.r.t the straight line distance along each
    %beam element
    r    = awi.model.Stick.getLineLength(coords(:, 1), coords(:, 2), coords(:, 3));
    etaR = r ./ r(end);
    %Get beam properties
    a = vertcat(PBar.a);
    i = vertcat(PBar.i);
    j = vertcat(PBar.j);
    %Assign to beam
    %   - Swap 11 & 22 planes as the Fame FEM defines the first plane
    %   parallel with the global z-direction. Typical convention is to
    %   define it  parallel with global x. (i.e. in-plane bending)
    set(Beams(iB), ...
        'I11', [i(:, 2) ; i(end, 2)], 'I11_eta', etaR, ...
        'I22', [i(:, 1) ; i(end, 1)], 'I22_eta', etaR, ...
        'I12', [i(:, 3) ; i(end, 3)], 'I12_eta', etaR, ...
        'J'  , [j ; j(end)]         , 'J_eta'  , etaR, ...
        'A'  , [a ; a(end)]         , 'A_eta'  , etaR);
    %Materials?
    if numel(Mat) > 1
        error('Update code for multiple material properties along the beam');
    end
    MatObj    = awi.model.Material;
    MatObj.E  = Mat.e;
    MatObj.G  = Mat.g;
    MatObj.Nu = Mat.nu;
    assignBeamObject(Beams(iB), [MatObj, MatObj], [0, 1]);
end

end

%	- Fuel/Mass distribution
function [AWIModel, FuelDistr] = assignLumpedMass(AWIModel, fameFEModel, logfcn)
%assignLumpedMass Assigns the lumped mass objects to the AWI objects based
%on the lumped mass distribution in the FAME FE model.
%
% Mass is assigned to the model based on the 'type' property of the 'Conm2'
% objects.
%
%   'Fuel' - Mass of type 'Fuel' is added as 'awi.model.FuelMass' objects
%            within an 'awi.model.FuelDistribution' object.
%
%   Otherwise - Mass is added as 'awi.model.PointMass' object to the
%               'awi.model.Beam' object defined by the 'part' property of
%               the 'Conm2' objects.
%
% This function returns the handle to the AWI model (typically an instance
% of 'awi.model.Aircraft') and also returns a cell array of
% 'awi.model.FuelDistribution' objects for pairing with any
% 'awi.model.LoadCase' objects.

%Grab all fuel masses - These are treated seperately to other masses
allTypes  = {fameFEModel.Conm2.type};
idxFuel   = ismember(allTypes, 'Fuel');
FuelMass  = fameFEModel.Conm2(idxFuel);
OtherMass = fameFEModel.Conm2(~idxFuel);

%Generate non-fuel objects
AWIModel = i_convertConm2_2_AWI(AWIModel, fameFEModel, OtherMass, logfcn);

%Generate 'awi.model.FuelDistribution' objects
FuelDistr = i_generateFuelDistribution(AWIModel, fameFEModel, FuelMass, logfcn);

    function FuelDistr = i_generateFuelDistribution(AWIModel, fameFEModel, FuelMass, logfcn)
        %i_generateFuelDistribution Generate the
        %'awi.model.FuelDistribution' objects for the different fuel cases
        %defined in the 'fameFEModel'.
        
        %Fuel cases are defined based on the 'file' property of the Conm2 object
        allFuelCases = {FuelMass.file};
        uFuelCases   = unique(allFuelCases);
        
        FuelDistr = cell(1, numel(uFuelCases));
        
        %Generate fuel objects
        for iF = 1 : numel(uFuelCases) %Loop through fuel cases
            
            %Grab all fuel masses for this fuel case
            idx = ismember(allFuelCases, uFuelCases{iF});
            fm  = FuelMass(idx);
            
            %Generate the AWI object but don't assign them to their parents!
            [~, MassData_] = i_convertConm2_2_AWI(AWIModel, fameFEModel, fm, logfcn, false);
            
            %Make the 'awi.model.FuelDistribution' classes
            fd = arrayfun(@(i) awi.model.FuelDistribution( ...
                'ParentBeam', MassData_(i).AWIParent, ...
                'FuelMasses', MassData_(i).AWIMass  , ...
                'FameFuelFile' , uFuelCases{iF}     , ...
                'PartName'     , MassData_(i).AWIParent.Name , ...
                'AeroFlexConm2', MassData_(i).AeroFlexMass)  , ...
                1 : numel(MassData_), 'Unif', false);
            FuelDistr{iF} = horzcat(fd{:});
            
        end
    end

    function [AWIModel, varargout] = i_convertConm2_2_AWI(AWIModel, fameFEModel, Conm2, logfcn, addFlag)
        %i_convertConm2_2_AWI Converts the AeroFlex Conm2 objects into
        %'awi.model.PointMass' objects and assigns these objects to the
        %relative parts of the model 'AWIModel'.
        %
        % If the user has set addFlag to false then this function will not
        % add the mass objects to the AWI model but will instead return a
        % structure 'MassData' containing the vector of handles to the AWI
        % objects parenting the mass objects and a cell array containing
        % the respective mass objects.
        
        if nargin < 5
            addFlag = true;
        end
        
        %Grab all names of parent structures - excluding fuel masses
        allMassParents = {Conm2.part};
        massParents    = unique(allMassParents);
        nMassPar       = numel(massParents);
        
        %Preallocate
        MassData   = repmat(struct('AWIParent', [], 'AWIMass', [], ...
            'AeroFlexMass', []), [1, nMassPar]);
        
        %Generate non-fuel objects
        for ii = 1 : numel(massParents)
            
            %Find the AWI Framework object related to this part...
            awiObj = findall(AWIModel.Children, 'Name', massParents{ii});
            
            %Escape route
            if isempty(awiObj)%Tell the user how much mass they are skipping!
                
                %Which masses?
                idx = ismember(allMassParents, massParents{ii});
                
                %How much mass?
                M_skip = sum([Conm2(idx).m]);
                
                %Tell them
                logfcn([blanks(8) sprintf(['No part with name ''%s'' present ', ...
                    'in the FAME FE model. * * Skipping lumped mass ', ...
                    'generation for %.10g Kg * *'], massParents{ii}, M_skip)]);
                
                continue
                
            end
            
            %Grab the masses belonging to this part
            idx  = ismember(allMassParents, massParents{ii});
            type = {Conm2(idx).type};
            M    = Conm2(idx);
            
            %Down-select the masses based on type
            uType = unique(type);
            m     = cell(1, numel(uType));
            
            %Assign a mass object for each type
            for iT = 1 : numel(uType)
                %Logical indexing
                idx = ismember(type, uType{iT});
                %Assign the lumped masses
                [~, m{iT}] = i_assignLumpedMass(awiObj, fameFEModel, M(idx), uType{iT}, logfcn, addFlag);
            end
            
            %Store the handles
            MassData(ii).AWIParent    = awiObj;         %AWI component object
            MassData(ii).AWIMass      = horzcat(m{:});  %AWI point mass objects
            MassData(ii).AeroFlexMass = M;              %Original AeroFlex objects
            
        end
        
        %Return the 'MassData' if the user asks for it
        if nargout == 2
            varargout{1} = MassData;
        end
        
    end

    function [awiObj, varargout] = i_assignLumpedMass(awiObj, fameFEModel, M, type, logfcn, addFlag)
        %i_assignLumpedMass Creates the required number of
        %'awi.model.PointMass' or 'awi.model.FuelMass' objects and assigns
        %them to a parent AWI object. Also, tags the mass objects based on
        %their 'type'.
        
        if nargin < 6
            addFlag = true;
        end
        
        nMass  = numel(M);
        
        %User message and mass object type is dependent on type of mass
        switch type
            case 'Fuel' %Fuel masss is special!
                %Use the 'fuel' object so we can collect these in the tree
                fn = @awi.model.FuelMass;
                %Tell the user which fuel case
                descr = sprintf('(Fuel case ''%s'')', M(1).file);
            otherwise %For all other types (e.g. Wing, HTP, VTP, etc.)
                %Just go with the normal point mass object
                fn = @awi.model.PointMass;
                %No additional description
                descr = '';
        end
        
        %Tell the user what is happening...
        logfcn([blanks(8) sprintf(['Generating ''%s'' lumped masses ', ...
            'for the ''%s'' object. %s'], type, awiObj.Name), descr]);
        
        %Grab the coordinates of the masses in the global coordinate system
        massCoords = i_getMassCoords(M, fameFEModel);
        
        %Position the masses in the coordinate frame of this part
        massLocalCoords = massCoords - awiObj.AbsPosition;
        
        %Collect the (3 x 3) inertia matrix
        I = cat(3, M.i);
        
        %Extract inertia data - Assume a symmetric inertia matrix
        I11 = num2cell(squeeze(I(1, 1, :)));
        I12 = num2cell(squeeze(I(1, 2, :)));
        I13 = num2cell(squeeze(I(1, 3, :)));
        I22 = num2cell(squeeze(I(2, 2, :)));
        I23 = num2cell(squeeze(I(2, 3, :)));
        I33 = num2cell(squeeze(I(3, 3, :)));
        
        %Create the mass objects
        MassObj = arrayfun(@(i) fn(), 1 : nMass, 'Unif', false);
        MassObj = [MassObj{:}];
        
        %Construct the 'MassGroup' string
        typeStr = [type, ' Masses'];
        
        %Assign the data to the masses
        set(MassObj, {'XOffset'}  , num2cell(massLocalCoords(:, 1)));
        set(MassObj, {'YOffset'}  , num2cell(massLocalCoords(:, 2)));
        set(MassObj, {'ZOffset'}  , num2cell(massLocalCoords(:, 3)));
        set(MassObj, {'Mass'}     , {M.m}');
        set(MassObj, {'Inertia11'}, I11);
        set(MassObj, {'Inertia12'}, I12);
        set(MassObj, {'Inertia13'}, I13);
        set(MassObj, {'Inertia22'}, I22);
        set(MassObj, {'Inertia23'}, I23);
        set(MassObj, {'Inertia33'}, I33);
        set(MassObj, 'MassGroup'  , typeStr);
        
        %Add to the parent AWI Framework object
        if addFlag
            awiObj.add(MassObj);
        end
        
        %Maybe the user wants the masses back?
        if nargout > 1
            varargout{1} = MassObj;
        end
        
    end

    function massCoords = i_getMassCoords(M, fameFEModel)
        %i_getMassCoords Grabs the coordinates of the FAME masses in the
        %global coordinate system using the GRID ID numbers.
        
        %Relate the masses to the GRID objects
        [~, gind] = ismember([M.grid], [fameFEModel.Grid.id]);
        
        %Grab the coordinates of these GRID objects
        gridCoords = vertcat(fameFEModel.Grid(gind).coord);
        
        %Add offsets to find the absolute position of the masses
        massCoords = gridCoords + vertcat(M.offset);
        
    end

end

