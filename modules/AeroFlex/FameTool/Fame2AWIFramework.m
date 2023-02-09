% The following function is in charge of converting a FAME structure,
% derived directly from the *.fm4 file, into one that is useable by
% the AWI framework.

function C = Fame2AWIFramework(FAME, logfcn)

if nargin < 2
   logfcn = @(s) fprintf('%s\n', s); 
end

logfcn([blanks(8) 'Initial Aircraft created ...']);
C = awi.model.Aircraft('Name','Aircraft');

logfcn([blanks(8) 'Fuselage added ...']);
F = C.add(@awi.model.BluffBody,'Name','Fuselage');

% EVALUATE USER VARIABLES
User_var = fieldnames(FAME.USER_VARIABLES);

for i = 1:length(User_var)
    eval([User_var{i,1} '= str2double(FAME.USER_VARIABLES.' User_var{i,1} '{1,1});']);
end

% Assign some of the global mass variables to the Aircraft Object
map = {'mtow','MTOM';...
    'mzfw','MZFM';...
    'mlw','MLM';...
    'owe','OEM'};

[b,~] = ismember(map(:,1),User_var);
for i = find(b)'
    C.(map{i,2}) = eval(map{i,1});
end

% Extract the Wing Reference point
ref_point = [str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.X_REF{1})/1000,...
    0,str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.Z_REF{1})/1000];

%% WING PLANFORM PROPERTIES

% -----  QCSweep ------ %
if isfield(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_LE')
    SweepType = 'SWEEP_LE';
    SweepLoc  = 0.0;
elseif isfield(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_25')
    SweepType = 'SWEEP_25';
    SweepLoc  = 0.25;
else
    error('Neither SWEEP_LE nor SWEEP_25 have been found for the wing'); 
end

QCsweep   = zeros(length(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType)),2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType))
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    QCsweep(i,:) = [Eta,str2double(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType){i})];
end
QCsweep(end+1,:) = [1,QCsweep(end,2)];

% -----  Thickness to Chord Ratio ------- %
thickness = zeros(length(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA),2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    thickness(i,:)  = [Eta,str2double(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i})];
end

% ------ Dihedral ----- %
dihedral = zeros(length(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG),2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    dihedral(i,:) = [Eta,str2double(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i})];
end
dihedral(end+1,:) = [1,dihedral(end,2)];

% ------- Chord -------- %
chord = zeros(length(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD),2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    chord(i,:) = [Eta,str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i})/1000];
end

% -------- Front Spar ------- %
fspar = zeros(length(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA),2);
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA)
    check_var = strtrim(regexp(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    fspar(i,:) = [Eta,str2double(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.CHORD_REL{i})];
end

% ---------Rear Spar ------- %
aspar = zeros(length(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA),2);
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA)
    check_var = strtrim(regexp(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    aspar(i,:) = [Eta,str2double(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.CHORD_REL{i})];
end

logfcn([blanks(8) 'Starboard wing added ...']);

% Assign a StbdWing to the Fuselage
W = F.add(@awi.model.LiftingSurface,'Name','StbdWing');

%Calculate the normalised beam location as the average of the front and aft
%spar positions
beamLoc = 0.5*(fspar(1,2)+aspar(1,2));

% Assign values to the StbdWing
W.Span          = str2double(FAME.WING.GLOBAL_GEOMETRY.SPAN{1})/2000;
W.SpanVector    = 'Y';
W.Origin        = ref_point + [(beamLoc * chord(1, 2)), 0, 0];
W.Chord_eta     = chord(:,1)';
W.Chord         = chord(:,2)';
W.Chord_flag    = W.SpanVector;
W.Sweep_eta     = QCsweep(:,1)';
W.Sweep         = QCsweep(:,2)';
W.SweepLoc      = SweepLoc;
W.Dihedral_eta  = dihedral(:,1)';
W.Dihedral      = dihedral(:,2)';
W.BeamLoc_eta   = getEta(W, 'Planform Property');   %Define the beam location at every eta position
W.BeamLoc       = repmat(beamLoc, size(W.BeamLoc_eta)); %Assue beam is at mid-chord

% Assign a front spar to the wing
Fs      = W.add(@awi.model.Spar,'Name','StbdWingFrontSpar');
Fs.Eta  = fspar(:,1)';
Fs.XLoc = fspar(:,2)';

% Assign a rear spar to the wing
Rs      = W.add(@awi.model.Spar,'Name','StbdWingRearSpar');
Rs.Eta  = aspar(:,1)';
Rs.XLoc = aspar(:,2)';

% ------ Control Surface  ------- %
count = 0;
splcount = 0;

% --------- SPOILERS ------- %
spoiler.eta_beg_le = [];
spoiler.eta_beg_te = [];
spoiler.eta_end_le = [];
spoiler.eta_end_te = [];
spoiler.chrd_beg_le = [];
spoiler.chrd_beg_te = [];
spoiler.chrd_end_le = [];
spoiler.chrd_end_te = [];

for i = 1:length(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.NO)
    
    count = count + 1;
    splcount = splcount + 1;
    Eta_beg_le  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_BEG_LE{i});
    Eta_beg_te  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_BEG_TE{i});
    Eta_end_le  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_END_LE{i});
    Eta_end_te  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_END_TE{i});
    Chrd_beg_le = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_BEG_LE{i});
    Chrd_beg_te = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_BEG_TE{i});
    Chrd_end_le = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_END_LE{i});
    Chrd_end_te = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_END_TE{i});
    
    Cs = W.add(@awi.model.ControlSurface,'Name',['spl' num2str(count)]);
    Cs.Eta = [Eta_beg_le,Eta_end_le];
    Cs.xLE = [Chrd_beg_le,Chrd_end_le];
    Cs.xTE = [Chrd_beg_te,Chrd_end_te];
    
end

% --------- AILERONS ------- %
aileron.eta_beg_le = [];
aileron.eta_beg_te = [];
aileron.eta_end_le = [];
aileron.eta_end_te = [];
aileron.chrd_beg_le = [];
aileron.chrd_end_le = [];

ailcount = 0;

for i = 1:length(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.NO)
    
    count       = count + 1;
    ailcount    = ailcount + 1;
    
    Eta_beg_le  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_BEG_LE{i});
    Eta_beg_te  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_BEG_TE{i});
    Eta_end_le  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_END_LE{i});
    Eta_end_te  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_END_TE{i});
    Chrd_beg_le = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.CHRD_BEG_LE{i});
    Chrd_end_le = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.CHRD_END_LE{i});
    
    aileron.eta_beg_le = [aileron.eta_beg_le;Eta_beg_le];
    aileron.eta_beg_te = [aileron.eta_beg_te;Eta_beg_te];
    aileron.eta_end_le = [aileron.eta_end_le;Eta_end_le];
    aileron.eta_end_te = [aileron.eta_end_te;Eta_end_te];
    aileron.chrd_beg_le = [aileron.chrd_beg_le;Chrd_beg_le];
    aileron.chrd_end_le = [aileron.chrd_end_le;Chrd_end_le];
    
    Cs = W.add(@awi.model.ControlSurface,'Name',['ail' num2str(ailcount)]);
    Cs.Eta = [Eta_beg_le,Eta_end_le];
    Cs.xLE = [Chrd_beg_le,Chrd_end_le];
    Cs.xTE = [1,1];
    
end

logfcn([blanks(8) 'Starboard wing successfully initalised.']);

%% VTP PROPERTIES

% Front Spar Position
fspar = [0,0.1;1,0.1];

% Rear Spar Position
aspar = [0,0.65;1,0.65];


%Calculate the normalised beam location as the average of the front and aft
%spar positions
beamLoc = 0.5*(fspar(1,2)+aspar(1,2));

logfcn([blanks(8) 'VTP added ...']);
% Define a new wing
V = F.add(@awi.model.LiftingSurface,'Name','VTP');

V.Span          = 5.875;
V.SpanVector    = 'Z';
Origin          = [38.76 + 0.5*(0.1+0.65)*5.0354, 0, 2.1145];
V.Origin        = Origin;
V.Chord_eta     = [0,1];
V.Chord         = 5.0334*[1,0.9985];
V.Sweep_eta     = [0,1];
V.Sweep         = [37.4879,37.4879];
V.SweepLoc      = 0.25;
V.BeamLoc_eta   = getEta(V, 'Planform Property');   %Define the beam location at every eta position
V.BeamLoc       = repmat(beamLoc, size(V.BeamLoc_eta)); %Assue beam is at mid-chord

V_Fs      = V.add(@awi.model.Spar,'Name','VTPFrontSpar');
V_Fs.Eta  = fspar(:,1)';
V_Fs.XLoc = fspar(:,2)';

V_Rs      = V.add(@awi.model.Spar,'Name','VTPRearSpar');
V_Rs.Eta  = aspar(:,1)';
V_Rs.XLoc = aspar(:,2)';

% CONTROL SURFACE PROPERTIES

%% Control Surface Information
%ELEVATOR
V_Cs = V.add(@awi.model.ControlSurface,'Name','rudr1');
V_Cs.Eta = [0.05,1];
V_Cs.xLE = [0.65,0.65];
V_Cs.xTE = [1,1];

logfcn([blanks(8) 'VTP successfully initialised.']);

%% HTP PROPERTIES

% -------- QC Sweep ------- %
QCsweep  = zeros(length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25),2);
for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25)
    check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    QCsweep(i,:) = [Eta,str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25{i})];
end
QCsweep(end+1,:) = [1,QCsweep(end,2)];

% -------- Thickness to Chord ------- %
thickness  = zeros(length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL),2);
for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
    check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    thickness(i,:) = [Eta,str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i})];
end

% -------- Dihedral ------- %
dihedral  = zeros(length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL),2);
for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
    check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    dihedral(i,:) = [Eta,str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i})];
end
dihedral(end+1,:) = [1,dihedral(end,2)];

% -------- Chord ------- %
chord  = zeros(length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD),2);
for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
    check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    chord(i,:) = [Eta,str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i})/1000];
end

% -------- Front Spar ------- %
fspar = [0,0.15;1,0.15];

% -------- Rear Spar ------- %
aspar = [0,0.65;1,0.65];

%Calculate the normalised beam location as the average of the front and aft
%spar positions
beamLoc = 0.5*(fspar(1,2)+aspar(1,2));

logfcn([blanks(8) 'HTP added ...']);
% Define a new wing
H = F.add(@awi.model.LiftingSurface,'Name','HTP');

H.Span = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SPAN{1})/2000;
H.SpanVector    = 'Y';
H.Origin        = [42.884 + 0.5*(0.15+0.65)*chord(1,2), 0, 7.061];
H.Chord_eta     = chord(:,1)';
H.Chord         = chord(:,2)';
H.Sweep_eta     = QCsweep(:,1)';
H.Sweep         = QCsweep(:,2)';
H.SweepLoc      = 0.25;
H.BeamLoc_eta   = getEta(H, 'Planform Property');   %Define the beam location at every eta position
H.BeamLoc       = repmat(beamLoc, size(H.BeamLoc_eta)); %Assue beam is at mid-chord
H.Dihedral_eta  = dihedral(:,1)';
H.Dihedral      = dihedral(:,2)';


H_Fs      = H.add(@awi.model.Spar,'Name','HTPFrontSpar');
H_Fs.Eta  = fspar(:,1)';
H_Fs.XLoc = fspar(:,2)';

H_Rs      = H.add(@awi.model.Spar,'Name','HTPRearSpar');
H_Rs.Eta  = aspar(:,1)';
H_Rs.XLoc = aspar(:,2)';

% CONTROL SURFACE PROPERTIES

%% Control Surface Information
%ELEVATOR
H_Cs = H.add(@awi.model.ControlSurface,'Name','elev1');
H_Cs.Eta = [0.05,1];
H_Cs.xLE = [0.65,0.65];
H_Cs.xTE = [1,1];

logfcn([blanks(8) 'HTP successfully initialised.']);

%% FUSELAGE PROPERTIES
F.Length = 43.913;
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

F.Origin = [0,0,0];

logfcn([blanks(8) 'Fuselage successfully initialised.']);
logfcn([blanks(8) 'Aircraft definition complete.']);

end