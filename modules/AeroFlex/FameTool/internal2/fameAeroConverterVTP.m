%% FAMEAEROCONVERTER   Extract FAME aerodynamic definitions and create aero
%                      panels
%
% Extract the FAME aerodynamic definition and create aerodynamic panels for
% NEOCASS and NASTRAN. A FAME model object is passed to this function from
% the FAME2MAT object.
%
% See NEOCASS and NASTRAN manuals for further explanations of the
% aerodynamic decks.
%
% Make sure the following inputs are defined in the FAME2MAT object
%
%     obj.Dirs.fameFolder       : FAME folder
%     obj.Inp.aeroPlanformFiles : auto generated
%     obj.Opts.Aero.nSpan       : Number of spanwise panels
%     obj.Opts.Aero.nChord      : Number of chordwise panels
%
%   Numbering conventions:
%     Right hand wing : 710001
%     Left hand wing  : 720001
%
%   See also |getFameAero| method within Fame2mat object.

%   Copyright 2016 University of Bristol
%   Private function.
function obj = fameAeroConverterVTP(obj)

Part = 'VTP';

if strcmp(Part,'HTP')
    startNumber = obj.Opts.Struct.HTPCAEROID;
    nSpan       = 1;
    nChord      = 12;
    TWIST       = 0;
elseif strcmp(Part,'Wing')
    startNumber = obj.Opts.Struct.WingCAEROID;
    nSpan       = obj.Opts.Aero.wing_nSpan;
    nChord      = obj.Opts.Aero.wing_nChord;
    TWIST       = obj.Opts.Aero.wing_addTwist;
elseif strcmp(Part,'VTP')
    startNumber = obj.Opts.Struct.VTPCAEROID;
    nSpan       = 1;
    nChord      = 12;
    TWIST       = 0;
end

Outline     = obj.Fame.Geometry.(Part).Outline;
RootNode    = obj.Fame.Geometry.(Part).ref_point;

Spanin      = min(Outline(:,3));
Spanout     = max(Outline(:,3));
Span = Spanout-Spanin;
npoints     = size(Outline,1)/2;

% % Find the control surfaces in order to interpolate the mesh correctly
%

CSEdges = uniquetol([obj.Fame.Geometry.(Part).ControlSurfaces.rudder.eta_beg_le;...
    obj.Fame.Geometry.(Part).ControlSurfaces.rudder.eta_end_le;(Outline(:,3)-Outline(1,3))/Span],0.005);

npointsNew = length(CSEdges);

ControlOutline = [];
ControlOutline(1:npointsNew,1) = interp1((Outline(1:npoints,3)-Outline(1,3))/Span,Outline(1:npoints,1),CSEdges);
ControlOutline(1:npointsNew,2) = interp1((Outline(1:npoints,3)-Outline(1,3))/Span,Outline(1:npoints,2),CSEdges);
ControlOutline(1:npointsNew,3) = Span*CSEdges;

ControlOutline(npointsNew+1: 2*npointsNew,1) = interp1((Outline(npoints+1:2*npoints,3)-Outline(1,3))/Span,Outline(npoints+1:2*npoints,1),flipud(CSEdges));
ControlOutline(npointsNew+1: 2*npointsNew,2) = flipud(ControlOutline(1:npointsNew,2));
ControlOutline(npointsNew+1: 2*npointsNew,3) = flipud(ControlOutline(1:npointsNew,3));

% Separate out the kinks
[a,b,CaeroObj]   = unique(ControlOutline(:,3),'rows');

% Extract the number of sections
nSections = length(a) - 1;
Span      = max(a);

% CSEdges
% Outline
% nSections
% Intdata = unique(round(10000*[CSEdges;Outline(:,2)])/10000)
% interp1(Outline(1:nSections-1,1),Outline(1:nSections-1,2),Intdata)
% pause;
% Work out the number of panels to go in each section
secFrac   = diff(ControlOutline(1:nSections + 1,3))./ Span;
nyPanels  = secFrac * nSpan;
nyPanels  = round(nyPanels);
nyPanels(nyPanels == 0) = 1;

accPanels = cumsum(nyPanels);
ip        = [0;accPanels(1:end - 1)];

% The definition of sweep in the fame file is at the 1/4 chord. This has
% nothing to do with where the beam sits as the outline is used to directly
% place the mesh in 3D space
sweepLoc  = 0.25;

% The following for loop creates a single card for each spanwise panel in
% each section - to accomodate for finer twist distribution. Twist is
% defined on the panel card but not for plotting purposes of the structure.

for nCaero = 1:nSections
    chordIndexIn  = find(CaeroObj == nCaero);
    chordIndexOut = find(CaeroObj == nCaero + 1);
    
    % Linearly interpolate and evenly distribute points along the kinks
    tempx = linspace(ControlOutline(b(nCaero),1),ControlOutline(b(nCaero + 1),1),nyPanels(nCaero) + 1);
    tempy = linspace(ControlOutline(b(nCaero),2),ControlOutline(b(nCaero + 1),2),nyPanels(nCaero) + 1);
    tempz = linspace(ControlOutline(b(nCaero),3),ControlOutline(b(nCaero + 1),3),nyPanels(nCaero) + 1);
    
    % Populate the CAERO cards with the inboard points
    CAERO.geo.startX(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = tempx(1:nyPanels(nCaero));
    CAERO.geo.startY(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = tempy(1:nyPanels(nCaero));
    CAERO.geo.startZ(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = tempz(1:nyPanels(nCaero));
    
    % Interpolate the twist along the span of the wing
    CAERO.geo.twist1(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = 0;% - rootangle;
    CAERO.geo.twist2(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = 0;% - rootangle;
    
    % Linearly interpolate and evenly distribute the chord along the
    % section
    temp_chord =  linspace(...
        ControlOutline(chordIndexIn(2),1) - ControlOutline(chordIndexIn(1),1), ...
        ControlOutline(chordIndexOut(2),1) - ControlOutline(chordIndexOut(1),1),nyPanels(nCaero) + 1);
    
    % Assign the chord to the structure
    CAERO.geo.c(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = temp_chord(1:nyPanels(nCaero));
    
    % Assign and calculate the taper ratio
    CAERO.geo.t(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = temp_chord(2:nyPanels(nCaero) + 1)./ temp_chord(1:nyPanels(nCaero));
    
    % Assign and calculate the section span
    temp_span = (ControlOutline(nCaero + 1,3) - ControlOutline(nCaero,3)) / nyPanels(nCaero);
    CAERO.geo.b(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = temp_span * ones(nyPanels(nCaero),1);
    
    % Determine the leading edge sweep
    temp = 180 * atan((ControlOutline(nCaero + 1,1) + sweepLoc * (ControlOutline(chordIndexOut(2)) - ControlOutline(chordIndexOut(1)))...
        - ControlOutline(nCaero,1) - sweepLoc * (ControlOutline(chordIndexIn(2)) - ControlOutline(chordIndexIn(1)))) / (ControlOutline(nCaero + 1,3) - ControlOutline(nCaero,3))) / pi;
    CAERO.geo.sw(ip(nCaero) + 1:ip(nCaero) + nyPanels(nCaero),1) = temp * ones(nyPanels(nCaero),1);
    
    for i = 1:nyPanels(nCaero)
        CAERO.geo.rootAirfoil{ip(nCaero) + i,1} = blanks(8);
        CAERO.geo.tipAirfoil{ip(nCaero) + i,1}  = blanks(8);
        %CAERO.geo.rootAirfoil{ip(nCaero) + i,1} = ['fame' num2str(ip(nCaero) + i)];
        %CAERO.geo.tipAirfoil{ip(nCaero) + i,1}  = ['fame' num2str(ip(nCaero) + i + 1)];
        CAERO.geo.dih(ip(nCaero) + i,1)         = 90; % TODO added dihedral entry
        CAERO.geo.meshType(ip(nCaero) + i,1)    = 1;
        CAERO.geo.cid(ip(nCaero) + i,1)         = 0;
    end
    
end

%% Allocate number of panels to the mesh
CAERO.id(1) = 1 + startNumber;

CAERO.geo.ny = ones(accPanels(end),1);
CAERO.geo.nx = nChord * ones(accPanels(end),1);

for i = 2:accPanels(end)
    CAERO.id(i,1) = CAERO.id(i - 1,1) + CAERO.geo.ny(i - 1,1) * CAERO.geo.nx(i - 1,1);
end

% offset the aero mesh to line up with the root
CAERO.geo.startZ = CAERO.geo.startZ + RootNode(3);
Outline(:,3)     = Outline(:,3);
ControlOutline(:,3)  = ControlOutline(:,3) + RootNode(3);
% Update the Outline entry:
obj.Fame.Geometry.(Part).Outline        = Outline;
obj.Fame.Geometry.(Part).ControlOutline = ControlOutline;

obj.Fame.Geometry.(Part).chord    = [(CAERO.geo.startZ - CAERO.geo.startZ(1))/(CAERO.geo.startZ(end)+CAERO.geo.b(end) - CAERO.geo.startZ(1));1];
obj.Fame.Geometry.(Part).chord(:,2) = [CAERO.geo.c;CAERO.geo.c(end)];

obj.Fame.Geometry.(Part).dihedral = [(CAERO.geo.startZ - CAERO.geo.startZ(1))/(CAERO.geo.startZ(end)+CAERO.geo.b(end) - CAERO.geo.startZ(1));1];
obj.Fame.Geometry.(Part).dihedral(:,2) = [CAERO.geo.dih;CAERO.geo.dih(end)];

obj.Fame.Geometry.(Part).sweep = [(CAERO.geo.startZ - CAERO.geo.startZ(1))/(CAERO.geo.startZ(end)+CAERO.geo.b(end) - CAERO.geo.startZ(1));1];
obj.Fame.Geometry.(Part).sweep(:,2) = [CAERO.geo.sw;CAERO.geo.sw(end)];
obj.Fame.Geometry.(Part).sweep(:,3) = 0.25;

nkinks   = size(ControlOutline,1) / 2;
leKinks  = ControlOutline(1:nkinks,:);
teKinks  = ControlOutline(nkinks + 1:end,:);

[~,idx] = sort(leKinks(:,2),'ascend');
leKinks = leKinks(idx,:);
[~,idx] = sort(teKinks(:,2),'ascend');
teKinks = teKinks(idx,:);

for i = 1:numel(CAERO.id)
    CaeroObj             = Caero;
    CaeroObj.id          = CAERO.id(i);
    CaeroObj.startX      = CAERO.geo.startX(i);
    CaeroObj.startY      = CAERO.geo.startY(i);
    CaeroObj.twist1      = CAERO.geo.twist1(i);
    CaeroObj.twist2      = CAERO.geo.twist2(i);
    CaeroObj.startZ      = CAERO.geo.startZ(i);
    CaeroObj.c           = CAERO.geo.c(i);
    CaeroObj.t           = CAERO.geo.t(i);
    CaeroObj.b           = CAERO.geo.b(i);
    CaeroObj.sw          = CAERO.geo.sw(i);
    CaeroObj.rootAirfoil = CAERO.geo.rootAirfoil{i};
    CaeroObj.tipAirfoil  = CAERO.geo.tipAirfoil{i};
    CaeroObj.dih         = CAERO.geo.dih(i);
    CaeroObj.meshType    = CAERO.geo.meshType(i);
    CaeroObj.cid         = CAERO.geo.cid(i);
    CaeroObj.ny          = CAERO.geo.ny(i);
    CaeroObj.nx          = CAERO.geo.nx(i);
    CaeroObj.csType      = [];
    CaeroObj.csId        = [];
    CaeroObj.csData      = [];
    %CaeroObj.outline     = Outline;
    CaeroObj.le          = leKinks;
    CaeroObj.te          = teKinks;
    CaeroObj.part        = 'VTP';
    CaeroObj.csLink      = 1;
    obj.Mdl.Caero = cat(2,obj.Mdl.Caero,CaeroObj);
end

% Add the ids to the part object
PartObj      = PartId;
PartObj.id   = 1;
PartObj.part = Part;
PartObj.type = 'CAERO';
PartObj.data = CAERO.id';

obj.Mdl.PartId  = cat(2,obj.Mdl.PartId,PartObj);
end


